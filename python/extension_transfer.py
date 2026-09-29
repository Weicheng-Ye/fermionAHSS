"""Lazy resolution coordinates for the fixed normalized-bar stacking model.

The comparison supplied by GAP must be an integral strong deformation
retraction. Universal formulas are imported through the package worker;
their calibration is unchanged. Branch predicates test whether a cochain
vanishes on every basis element of the resolution in its degree, through
the comparison chains g(e_j); no bar simplices are enumerated.

This implements the conditional cochain transfer, not a proof of completeness
of the transferred gauge relation or a ko classification.
Copyright (c) 2026 koAHSS contributors; MIT license.
"""
from functools import lru_cache
import gc
import itertools
import json
import sys
import weakref

from extension_worker import api, BoundedCache
import extension_acceleration as acceleration
acceleration.install()
import off_shell_beta as off
import stacking_lower as lower
import all_cochain_upper as upper
import extension_native_upper as native
from extension_three_local import ThreeLocalModel
from coherent_low_commutative import (DegreeOneCommutativeStacking,
                                     DegreeTwoCommutativeStacking)
from a0_gamma import beta2
from cochain_tools import LinearCochain

p = api.p
structural_zero = acceleration.is_zero
FIELDS = "ABCD"


class TransferResourceLimit(RuntimeError):
    pass


def degrees(k):
    return k - 3, k - 2, k - 1, k + 1


def integer(value):
    if hasattr(value, "denominator") and value.denominator != 1:
        raise ArithmeticError("nonintegral value in the transferred stacking model")
    if not isinstance(value, int) and not hasattr(value, "denominator"):
        raise TypeError("an exact integral cochain value is required")
    return int(value)


class LazyState:
    """A triangular state whose upper formulas are built only when needed.

    build(state, layer) and, for products, cross(state, layer) receive the
    state they fill, so no closure refers back to it. A state and its
    cochain memo tables are then released by reference counting once
    unused; the worker runs without automatic garbage collection.
    """
    def __init__(self, k, build, cross=None):
        self.k = k
        self.build = build
        self.layers = {}
        self.cross_terms = cross
        self.crosses = {}

    def __getitem__(self, layer):
        if layer not in self.layers:
            self.layers[layer] = self.build(self, layer)
        return self.layers[layer]

    def cross(self, layer):
        """The formula term of a product layer beyond its two summands."""
        if layer not in self.crosses:
            self.crosses[layer] = self.cross_terms(self, layer)
        return self.crosses[layer]


class TransferredModel:
    def __init__(self, setup, transport):
        if setup.get("schema", 1) != 1:
            raise ValueError("unsupported extension transfer schema")
        # Vertices are integer labels. With a multiplication table they are
        # group elements and are normalized here; otherwise GAP normalizes
        # every simplex it receives, and equal labels mean equal vertices.
        self.vertex_mode = setup.get("vertexMode", "table")
        if self.vertex_mode == "table":
            self.mul = setup["multiplication"]
            self.size = len(self.mul)
            if not self.size or any(len(row) != self.size for row in self.mul):
                raise ValueError("invalid group multiplication table")
            if any(type(i) is not int or not 0 <= i < self.size
                   for row in self.mul for i in row):
                raise ValueError("invalid group element index")
            self.inverse = [next(h for h in range(self.size)
                                 if self.mul[g][h] == self.mul[h][g] == 0)
                            for g in range(self.size)]
        elif self.vertex_mode == "labels":
            self.mul = self.size = self.inverse = None
        else:
            raise ValueError("unknown vertex mode")
        self.ranks = setup["ranks"]
        self.matrices = {False: setup["ordinary"], True: setup["signed"]}
        self.max_degree = len(self.ranks) - 1
        self.max_debug_simplices = setup.get("maxDebugSimplices", 8192)
        if type(self.max_debug_simplices) is not int or self.max_debug_simplices < 1:
            raise ValueError("maxDebugSimplices must be positive")
        if any(type(n) is not int or n < 0 for n in self.ranks):
            raise ValueError("invalid resolution ranks")
        # The comparison chains g(e_j): supplied in full, or requested from
        # GAP one basis element at a time when first needed (gMode "lazy").
        self._chains = {}
        if "g" in setup:
            g = setup["g"]
            if len(g) != len(self.ranks):
                raise ValueError("g must cover every supplied resolution degree")
            for n, rank in enumerate(self.ranks):
                if len(g[n]) != rank:
                    raise ValueError("g has the wrong resolution rank")
                for j, chain in enumerate(g[n]):
                    self._chains[n, j] = chain
        elif setup.get("gMode") != "lazy":
            raise ValueError("setup needs the comparison chains g or gMode 'lazy'")
        # Pairings <c, g(e_j)> per cochain (weight 0: modulo two, weight 1:
        # exact), kept without keeping the cochains alive, and the retraction
        # identities checked per degree (see _pairing).
        self._pairings = (weakref.WeakKeyDictionary(), weakref.WeakKeyDictionary())
        self._identities = {}
        for signed, matrices in self.matrices.items():
            if len(matrices) < self.max_degree:
                raise ValueError("missing resolution coboundary matrices")
            for n in range(self.max_degree):
                if len(matrices[n]) != self.ranks[n] or any(
                    len(row) != self.ranks[n + 1]
                    or any(type(v) is not int for v in row)
                    for row in matrices[n]):
                    raise ValueError("invalid resolution coboundary dimensions")
        self.transport = transport
        self._transport = BoundedCache(32768)
        self._answers = BoundedCache(256)
        self._phis = BoundedCache(128)
        self._kappas = BoundedCache(256)
        self._cores = BoundedCache(256)
        self.low2 = DegreeTwoCommutativeStacking()
        self.low1 = DegreeOneCommutativeStacking(successor=self.low2)
        # Requests carrying ``prime: 3`` use the two-layer three-local model.
        self.three_local = ThreeLocalModel(self)
        self.s_vector = self.check_vector(1, setup["s"], False)
        self.w_vector = self.check_vector(2, setup["omega"], False)
        if any(v % 2 for v in self.coboundary(1, self.s_vector, False)):
            raise ValueError("s is not a binary cocycle on the resolution")
        if any(v % 2 for v in self.coboundary(2, self.w_vector, False)):
            raise ValueError("omega is not a binary cocycle on the resolution")
        self.s = self.lift(1, self.s_vector, False)
        self.omega = self.lift(2, self.w_vector, False)

    def dimension(self, n):
        if n < 0:
            return 0
        if n > self.max_degree:
            raise ValueError("requested degree exceeds supplied resolution")
        return self.ranks[n]

    def check_vector(self, n, vector, signed):
        if not isinstance(vector, (list, tuple)) or len(vector) != self.dimension(n):
            raise ValueError("incorrect resolution cochain dimension")
        if any(type(v) is not int for v in vector):
            raise ValueError("resolution cochains must have integral coordinates")
        if not signed and any(v not in (0, 1) for v in vector):
            raise ValueError("binary resolution cochains must contain zero or one")
        return tuple(vector)

    def check_state(self, k, data):
        return {f: self.check_vector(n, data[f], f in "AD")
                for f, n in zip(FIELDS, degrees(k))}

    def zero(self, k):
        return {f: [0] * self.dimension(n) for f, n in zip(FIELDS, degrees(k))}

    def coboundary(self, n, vector, signed):
        if n < 0:
            return [0] * self.dimension(n + 1)
        if n >= self.max_degree:
            raise ValueError("missing successor resolution degree")
        return [sum(v * self.matrices[signed][n][i][j]
                    for i, v in enumerate(vector))
                for j in range(self.dimension(n + 1))]

    def normalize(self, vertices):
        if not vertices:
            return ()
        if any(type(g) is not int or g < 0 or (self.size is not None and g >= self.size)
               for g in vertices):
            raise ValueError("transport accepts only vertex labels")
        if any(a == b for a, b in zip(vertices, vertices[1:])):
            return None
        if self.mul is None:
            return tuple(vertices)
        first = self.inverse[vertices[0]]
        return tuple(self.mul[first][g] for g in vertices)

    def terms(self, kind, vertices):
        sigma = self.normalize(vertices)
        if sigma is None:
            return ()
        key = kind, sigma
        if key not in self._transport:
            answer = self.transport(kind, len(sigma) - 1, sigma)
            if isinstance(answer, dict):
                if answer.get("status") != "computed":
                    raise TransferResourceLimit(answer.get("reason", "transport refused"))
                answer = answer["terms"]
            self._transport[key] = tuple(tuple(t) for t in answer)
        return self._transport[key]

    def chain(self, n, j):
        """g(e_j) for the basis element j (0-based) of R_n, as [c, c*chi, labels] terms."""
        found = self._chains.get((n, j))
        if found is None:
            if not 0 <= n <= self.max_degree or not 0 <= j < self.ranks[n]:
                raise ValueError("comparison chain outside the supplied resolution basis")
            answer = self.transport("g", n, j)
            if isinstance(answer, dict):
                if answer.get("status") != "computed":
                    raise TransferResourceLimit(answer.get("reason", "transport refused"))
                answer = answer["terms"]
            found = tuple((t[0], t[1], tuple(t[2])) for t in answer)
            self._chains[n, j] = found
        return found

    def lift(self, n, vector, signed):
        if n < 0 or not any(vector):
            return p.zero(n)
        weight = 2 if signed else 1
        def value(vertices):
            result = sum(t[weight] * vector[t[0]] for t in self.terms("f", vertices))
            return result if signed else result % 2
        answer = p.Cochain(n, value)
        answer.native_lift = tuple(vector), signed
        return answer

    def project(self, cochain, signed):
        n = cochain.degree
        if n < 0:
            return []
        if structural_zero(cochain):
            return [0] * self.ranks[n]
        values = self._pairing(cochain, 1 if signed else 0)
        return [integer(v) if signed else integer(v) % 2 for v in values]

    def _pairing(self, cochain, weight):
        """<cochain, g(e_j)> for every basis element j of R in its degree:
        exact for weight 1, correct modulo two for weight 0.

        Pairing is linear, so a sum is paired term by term, and a pairing is
        computed once per cochain. A lift of v pairs to v and a homotopy image
        pairs to zero (f g = 1 and H g = 0), each identity checked once per
        degree on the chains themselves; otherwise the cochain is evaluated.
        """
        table = self._pairings[weight]
        found = table.get(cochain)
        if found is not None:
            return found
        n = cochain.degree
        rank = self.ranks[n]
        if structural_zero(cochain):
            values = [0] * rank
        elif isinstance(cochain, LinearCochain) and (
                cochain.modulus is None or (cochain.modulus == 2 and weight == 0)):
            values = [0] * rank
            for term, k in cochain.terms:
                part = self._pairing(term, weight)
                values = [a + k * b for a, b in zip(values, part)]
            if weight == 0:
                values = [v % 2 for v in values]
        else:
            values = self._identity_pairing(cochain, weight)
            if values is None:
                values = [sum(t[weight] * cochain(tuple(t[2])) for t in self.chain(n, j))
                          for j in range(rank)]
        table[cochain] = values
        return values

    def _identity_pairing(self, cochain, weight):
        n = cochain.degree
        native = getattr(cochain, "native_lift", None)
        if native is not None:
            vector, signed = native
            if (signed or weight == 0) and self._identity_holds("lift", n, weight, signed):
                return [v % 2 for v in vector] if weight == 0 else list(vector)
            return None
        signed = getattr(cochain, "homotopy_image", None)
        if signed is not None and (signed or weight == 0) and \
                self._identity_holds("homotopy", n, weight, signed):
            return [0] * self.ranks[n]
        return None

    def _identity_holds(self, kind, n, weight, signed):
        """Whether <lift(e_i), g(e_j)> = delta_ij (kind "lift") or the chain
        h(g(e_j)) vanishes (kind "homotopy") in degree n, as the pairing of
        that weight sees them; computed once from the f or h terms of every
        simplex of the chains of degree n."""
        key = kind, n, weight, signed
        answer = self._identities.get(key)
        if answer is None:
            answer = True
            inner = (2 if signed else 1) if kind == "lift" else (1 if signed else 0)
            for j in range(self.ranks[n]):
                total = {}
                for t in self.chain(n, j):
                    for u in self.terms("f" if kind == "lift" else "h", tuple(t[2])):
                        target = u[0] if kind == "lift" else tuple(u[2])
                        total[target] = total.get(target, 0) + t[weight] * u[inner]
                if weight == 0 or not signed:
                    total = {x: v % 2 for x, v in total.items()}
                total = {x: v for x, v in total.items() if v}
                expected = {j: 1} if kind == "lift" else {}
                if total != expected:
                    answer = False
                    break
            self._identities[key] = answer
        return answer

    def homotopy(self, cochain, signed):
        if cochain.degree <= 0 or structural_zero(cochain):
            return p.zero(cochain.degree - 1)
        weight = 1 if signed else 0
        def value(vertices):
            result = sum(t[weight] * cochain(tuple(t[2]))
                         for t in self.terms("h", vertices))
            return result if signed else result % 2
        answer = p.Cochain(cochain.degree - 1, value)
        answer.homotopy_image = signed
        return answer

    @lru_cache(None)
    def simplices(self, n):
        """All normalized bar simplices; only for developer value exports."""
        if n < 0:
            return ()
        if self.mul is None:
            raise ValueError("bar simplices can be enumerated only with a finite group table")
        count = (self.size - 1) ** n
        if count > self.max_debug_simplices:
            raise TransferResourceLimit(
                f"bar value export in degree {n} needs {count} simplices "
                f"(limit {self.max_debug_simplices})")
        result = []
        for increments in itertools.product(range(1, self.size), repeat=n):
            vertices = [0]
            for g in increments:
                vertices.append(self.mul[vertices[-1]][g])
            result.append(tuple(vertices))
        return tuple(result)

    def is_zero(self, cochain, binary):
        """The cochain vanishes on g(e_j) for every basis element e_j of R.

        Binary cochains are paired modulo two, integral ones with the local
        sign coefficients of the comparison chains.
        """
        n = cochain.degree
        if n < 0 or structural_zero(cochain):
            return True
        if n >= len(self.ranks):
            raise ValueError("zero test exceeds the supplied resolution degrees")
        weight = 0 if binary else 1
        known = self._pairings[weight].get(cochain)
        if known is not None:
            return all((v % 2 if binary else v) == 0 for v in known)
        # Chains after the first nonzero pairing are never requested.
        for j in range(self.ranks[n]):
            value = sum(t[weight] * cochain(tuple(t[2])) for t in self.chain(n, j))
            if (value % 2 if binary else value) != 0:
                return False
        return True

    def differential(self, cochain, signed):
        if cochain.degree < 0:
            return p.zero(cochain.degree + 1)
        return p.ds(cochain, self.s) if signed else p.binary(p.differential(cochain))

    def legal_pair(self, state):
        if state.k == 2:
            return self.is_zero(self.differential(state[1], False), True)
        a, b = off.curvature(state[0], state[1], self.s, self.omega)
        return self.is_zero(a, False) and self.is_zero(b, True)

    def triple(self, state, full=True):
        legal = self.legal_pair(state)
        if not full:
            return upper.Triple(state[0], state[1], state[2], legal)
        pure = (self.is_zero(state[0], False) and self.is_zero(state[1], True)
                and self.is_zero(self.differential(state[2], False), True))
        complete = legal and self.is_zero(p.binary(self.differential(state[2], False)
            + self.nonlinear(state, 2)), True)
        return upper.Triple(state[0], state[1], state[2], legal, pure, complete)

    def nonlinear(self, state, layer):
        key = ("nonlinear", layer)
        if key in state.layers:
            return state.layers[key]
        k = state.k
        degree = degrees(k + 1)[layer]
        if k <= 0 or layer == 0 or (k <= 2 and layer == 1) or (k == 1 and layer == 2):
            value = p.zero(degree)
        elif k == 1:
            value = api.cochain(self.low1.g(*map(api.local, (state[2], self.s, self.omega))))
        elif k == 2:
            if layer == 2:
                value = p.QD(state[1], self.s, self.omega)
            else:
                value = api.cochain(self.low2.g(
                    *map(api.local, (state[1], state[2], self.s, self.omega)),
                    b_closed=self.is_zero(self.differential(state[1], False), True)))
        elif k in (3, 4, 5, 6):
            if layer == 1:
                value = off.primary(state[0], self.s, self.omega)
            elif layer == 2:
                value = off.current_f(off.LowerPair(state[0], state[1],
                    self.legal_pair(state)), self.s, self.omega)
            else:
                # Native curvature reaches this layer only after the lower
                # layers vanish, so the degree-six section branch of the note
                # is never needed.
                value = native.g(k, self.triple(state, full=False), self.s, self.omega,
                                 a_zero=self.is_zero(state[0], False),
                                 b_zero=self.is_zero(state[1], True))
        else:
            raise ValueError("transferred nonlinear differential covers degrees 0 through 6")
        state.layers[key] = value
        return value

    def phi(self, k, data):
        checked = self.check_state(k, data)
        d_coordinates = checked["D"]
        checked = dict(checked, D=(0,) * self.dimension(k + 1))
        key = k, tuple(checked[f] for f in "ABC")
        if key not in self._phis:
            def build(result, layer):
                signed = layer in (0, 3)
                value = self.lift(degrees(k)[layer], checked[FIELDS[layer]], signed)
                if layer:
                    value = value - self.homotopy(self.nonlinear(result, layer), signed)
                return value if signed else p.binary(value)
            self._phis[key] = LazyState(k, build)
        base = self._phis[key]
        if not any(d_coordinates):
            return base
        # Every nonlinear term reads ABC only. Reuse its entire cochain DAG
        # while preserving the exact signed integral D lift separately.
        d_lift = self.lift(k + 1, d_coordinates, True)
        return LazyState(k, lambda _, layer: base[layer] if layer < 3 else base[3] + d_lift)

    def bar_d(self, state):
        def build(_, layer):
            value = self.differential(state[layer], layer in (0, 3))
            value = value + self.nonlinear(state, layer)
            return value if layer in (0, 3) else p.binary(value)
        return LazyState(state.k + 1, build)

    def bar_product(self, left, right):
        if left.k != right.k or left.k not in (1, 2, 3, 4, 5, 6):
            raise ValueError("transferred products cover equal degrees 1 through 6")
        k = left.k
        def low_cross(result, layer):
            # The degree-one and degree-two products of the complete-bar
            # reference: (C+C', D+D'+gamma1) and (B+B', C+C'+beta2, D+D'+gamma2).
            if layer < 2 or (k == 1 and layer == 2):
                return p.zero(degrees(k)[layer])
            if layer == 2:
                return api.cochain(beta2(*map(api.local, (left[1], right[1], self.s))))
            if k == 1:
                return api.cochain(self.low1.gamma(*map(api.local,
                    (left[2], right[2], self.s, self.omega))))
            return api.cochain(self.low2.gamma(*map(api.local,
                (left[1], left[2], right[1], right[2], self.s, self.omega)),
                b_closed=self.legal_pair(left), bp_closed=self.legal_pair(right),
                sum_closed=self.legal_pair(result)))
        def cross(result, layer):
            if k <= 2:
                return low_cross(result, layer)
            if layer == 0:
                return p.zero(degrees(left.k)[0])
            if layer == 1:
                return lower.alpha(left[0], right[0], self.s)
            legal = self.legal_pair(result)
            if layer == 2:
                return lower.all_cochain_beta(
                    off.LowerPair(left[0], left[1], self.legal_pair(left)),
                    off.LowerPair(right[0], right[1], self.legal_pair(right)),
                    legal, self.s, self.omega)
            pure = (self.is_zero(result[0], False) and self.is_zero(result[1], True)
                    and self.is_zero(self.differential(result[2], False), True))
            a_zero = (self.is_zero(left[0], False), self.is_zero(right[0], False))
            b_zero = self.is_zero(left[1], True) and self.is_zero(right[1], True)
            return native.gamma(k, self.triple(left), self.triple(right),
                                legal, pure, a_zero, self.s, self.omega, b_zero)
        def build(result, layer):
            value = left[layer] + right[layer] + result.cross(layer)
            return value if layer in (0, 3) else p.binary(value)
        return LazyState(left.k, build, cross)

    def kappa(self, k, data, upto=3):
        checked = self.check_state(k, data)
        key = k, upto, tuple(checked[f] for f in "ABC")
        if key not in self._kappas:
            core = dict(checked, D=(0,) * self.dimension(k + 1))
            image = self.phi(k, core)
            answer = self.zero(k + 1)
            for layer in range(upto + 1):
                f = FIELDS[layer]
                signed = layer in (0, 3)
                linear = self.coboundary(degrees(k)[layer], core[f], signed)
                nonlinear = self.project(self.nonlinear(image, layer), signed)
                answer[f] = [a + b if signed else (a + b) % 2
                             for a, b in zip(linear, nonlinear)]
                if any(answer[f]):
                    break
            self._kappas[key] = tuple(tuple(answer[f]) for f in FIELDS)
        answer = {f: list(v) for f, v in zip(FIELDS, self._kappas[key])}
        if upto == 3 and not any(any(answer[f]) for f in "ABC"):
            carry = self.coboundary(k + 1, checked["D"], True)
            answer["D"] = [a + b for a, b in zip(answer["D"], carry)]
        return answer

    def require_flat(self, k, data, upto=3):
        if upto >= 0 and any(any(v) for v in self.kappa(k, data, upto).values()):
            raise ValueError("transferred operation requires flat inputs through its requested layers")

    def reflect(self, state, upto=3):
        k = state.k
        native = self.zero(k)
        gauge = LazyState(k - 1, lambda _, layer: p.zero(degrees(k - 1)[layer]))
        # This embedding reads native coordinates only after they have been
        # determined. Each cached layer depends only on earlier coordinates.
        def build(embedded, layer):
            signed = layer in (0, 3)
            value = self.lift(degrees(k)[layer], tuple(native[FIELDS[layer]]), signed)
            if layer:
                value = value - self.homotopy(self.nonlinear(embedded, layer), signed)
            return value if signed else p.binary(value)
        embedded = LazyState(k, build)
        boundary = self.bar_d(gauge)
        reconstructed = self.bar_product(boundary, embedded)
        for layer in range(upto + 1):
            signed = layer in (0, 3)
            image, image_product = embedded, reconstructed
            if layer == 3:
                # The D terms read embedded layers A, B, C only. These equal
                # the cached embedding of the same coordinates, so its
                # formulas are shared with every other use of that state.
                image = self.phi(k, dict(native, D=[0] * self.dimension(k + 1)))
                image_product = self.bar_product(boundary, image)
            correction = self.nonlinear(gauge, layer)
            if layer:
                correction = correction - self.homotopy(
                    self.nonlinear(image, layer), signed)
            correction = correction + image_product.cross(layer)
            residual = state[layer] - correction
            if not signed:
                residual = p.binary(residual)
            native[FIELDS[layer]] = self.project(residual, signed)
            gauge.layers[layer] = self.homotopy(residual, signed)
        return native, gauge

    def product(self, k, left, right, upto=3):
        left, right = self.check_state(k, left), self.check_state(k, right)
        if not any(any(v) for v in left.values()) or not any(any(v) for v in right.values()):
            source = right if not any(any(v) for v in left.values()) else left
            answer = self.zero(k)
            for layer in range(upto + 1):
                answer[FIELDS[layer]] = list(source[FIELDS[layer]])
            return answer, LazyState(k - 1, lambda _, layer: p.zero(degrees(k - 1)[layer]))
        return self.reflect(self.bar_product(self.phi(k, left), self.phi(k, right)), upto)

    def act(self, k, gauge, canonical, upto=3):
        gauge = self.check_state(k - 1, gauge)
        canonical = self.check_state(k, canonical)
        if not any(any(v) for v in gauge.values()):
            answer = self.zero(k)
            for layer in range(upto + 1):
                answer[FIELDS[layer]] = list(canonical[FIELDS[layer]])
            return answer, LazyState(k - 1, lambda _, layer: p.zero(degrees(k - 1)[layer]))
        boundary = self.bar_d(self.phi(k - 1, gauge))
        return self.reflect(self.bar_product(boundary, self.phi(k, canonical)), upto)

    def divide_left(self, k, left, total, verify=True, upto=3):
        """The right factor with left * right = total through the layers up to ``upto``."""
        left = self.check_state(k, left)
        total = self.check_state(k, total)
        required = 3 if upto == 3 else upto - 1
        if verify:
            self.require_flat(k, left, required)
            self.require_flat(k, total, required)
        right = self.zero(k)
        for layer, f in enumerate(FIELDS[:upto + 1]):
            product, _ = self.product(k, left, right, upto=layer)
            if layer == 3 and not any(left["D"]):
                # Before its D coordinate is set, right has the final A, B, C.
                # This is the exact core of the verifying product left*right.
                key = ("xtimes", k, 3, tuple(left[g] for g in "ABC"),
                       tuple(tuple(right[g]) for g in "ABC"))
                if key not in self._cores:
                    self._cores[key] = tuple(tuple(product[g]) for g in FIELDS)
            right[f] = [a + b - c for a, b, c in zip(right[f], total[f], product[f])]
            if layer in (1, 2):
                right[f] = [v % 2 for v in right[f]]
        if verify:
            self.require_flat(k, right, required)
            product, _ = self.product(k, left, right, upto)
            if any(tuple(product[f]) != total[f] for f in FIELDS[:upto + 1]):
                raise ArithmeticError("transferred left division failed its exact product equality")
        return right

    def affine_core(self, operation, k, left, right, upto=3):
        """Cache ABC-dependent formulas, retaining D carries outside the key."""
        left_degree = k - 1 if operation == "act" else k
        x = self.check_state(left_degree, left)
        y = self.check_state(k, right)
        key = operation, k, upto, tuple(x[f] for f in "ABC"), tuple(y[f] for f in "ABC")
        if key not in self._cores:
            x0 = dict(x, D=(0,) * self.dimension(left_degree + 1))
            y0 = dict(y, D=(0,) * self.dimension(k + 1))
            if operation == "xtimes":
                result, _ = self.product(k, x0, y0, upto)
            elif operation == "act":
                result, _ = self.act(k, x0, y0, upto)
            else:
                result = self.divide_left(k, x0, y0, verify=False, upto=upto)
            self._cores[key] = tuple(tuple(result[f]) for f in FIELDS)
        answer = {f: list(v) for f, v in zip(FIELDS, self._cores[key])}
        if upto == 3:
            if operation == "act":
                carry = self.coboundary(k, x["D"], True)
                answer["D"] = [a + b + c for a, b, c in zip(answer["D"], carry, y["D"])]
            elif operation == "xtimes":
                answer["D"] = [a + b + c for a, b, c in zip(answer["D"], x["D"], y["D"])]
            else:
                answer["D"] = [a + b - c for a, b, c in zip(answer["D"], y["D"], x["D"])]
        return answer

    def values(self, state, simplices=None):
        return {f: [integer(state[layer](tuple(sigma))) for sigma in
                    (simplices[f] if simplices is not None else self.simplices(n))]
                for layer, (f, n) in enumerate(zip(FIELDS, degrees(state.k)))}

    def calculate(self, request):
        operation, k = request["operation"], request["degree"]
        if k not in (0, 1, 2, 3, 4, 5, 6) or (k == 0 and operation not in ("d", "phi_values")):
            raise ValueError("transferred states cover degrees 1 through 6; gauges include degree 0")
        key = json.dumps(request, sort_keys=True, separators=(",", ":"))
        if key in self._answers:
            return json.loads(self._answers[key])
        if "prime" in request:
            answer = self.three_local.calculate(request)
            answer["status"] = "computed"
            self._answers[key] = json.dumps(answer, separators=(",", ":"))
            return answer
        # Layer-limited requests: ``upto`` is the last layer index computed
        # (A=0 to D=3); the answer has zeros above it, and the inputs need to
        # be flat only through the layer below it.
        upto = request.get("upto", 3)
        if type(upto) is not int or not 0 <= upto <= 3:
            raise ValueError("upto must be a layer index from zero through three")
        required = 3 if upto == 3 else upto - 1
        if operation == "d":
            answer = {"state": self.kappa(k, request["state"], upto)}
        elif operation == "xtimes":
            self.require_flat(k, request["state"], required)
            self.require_flat(k, request["other"], required)
            result = self.affine_core("xtimes", k, request["state"], request["other"], upto)
            answer = {"state": result}
        elif operation == "act":
            self.require_flat(k, request["state"], required)
            result = self.affine_core("act", k, request["gauge"], request["state"], upto)
            answer = {"state": result}
        elif operation == "divide_left":
            self.require_flat(k, request["state"], required)
            self.require_flat(k, request["other"], required)
            result = self.affine_core("divide_left", k, request["state"], request["other"], upto)
            self.require_flat(k, result, required)
            product = self.affine_core("xtimes", k, request["state"], result, upto)
            if any(tuple(product[f]) != tuple(request["other"][f]) for f in FIELDS[:upto + 1]):
                raise ArithmeticError("transferred left division failed its exact product equality")
            answer = {"state": result}
        elif operation == "phi_values":
            answer = {"state": self.values(self.phi(k, request["state"]),
                                            request.get("simplices"))}
        elif operation == "gauge_values":
            if "gauge" in request:
                result, gauge = self.act(k, request["gauge"], request["state"])
            else:
                result, gauge = self.product(k, request["state"], request["other"])
            answer = {"state": result, "gauge": self.values(gauge, request.get("simplices"))}
        else:
            raise ValueError("unknown extension transfer operation")
        answer["status"] = "computed"
        self._answers[key] = json.dumps(answer, separators=(",", ":"))
        return answer


def serve():
    model = None
    acceleration.persist()
    def transport(kind, degree, vertices):
        # A comparison chain is addressed by its 0-based basis index.
        message = (dict(operation="transport", kind="g", degree=degree, basis=vertices)
                   if kind == "g" else
                   dict(operation="transport", kind=kind, degree=degree, vertices=vertices))
        print(json.dumps(message), flush=True)
        line = sys.stdin.readline()
        if not line:
            raise RuntimeError("GAP transport stopped during a callback")
        return json.loads(line)
    for line in sys.stdin:
        try:
            request = json.loads(line)
            if request["operation"] == "setup":
                model = TransferredModel(request, transport)
                answer = {"status": "computed"}
            elif model is None:
                raise ValueError("transferred model is not initialized")
            else:
                answer = model.calculate(request)
        except Exception as exc:
            unresolved = isinstance(exc, (TransferResourceLimit, native.NativeDegreeSixLimit))
            answer = {"status": "unresolved" if unresolved else "error",
                      "exception": type(exc).__name__, "reason": str(exc)}
        # GAP may stop the worker right after the answer; store values first.
        acceleration.flush()
        print(json.dumps(answer, separators=(",", ":")), flush=True)
        acceleration.recycle_sources()


if __name__ == "__main__":
    # Worker heaps are large, long-lived memo tables with almost no reference
    # cycles, so automatic collection only rescans them.
    gc.disable()
    serve()
