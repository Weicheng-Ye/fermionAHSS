"""Lazy resolution coordinates for the fixed normalized-bar stacking model.

The comparison supplied by GAP must be an integral strong deformation
retraction. Universal formulas are imported through the checksum-verified
worker; their source and calibration are unchanged. Global branch predicates
are evaluated on the complete bar only in their required lower degrees.

This implements the conditional cochain transfer, not a proof of completeness
of the transferred gauge relation or a ko classification.
Copyright (c) 2026 koAHSS contributors; MIT license.
"""
from functools import lru_cache
import itertools
import json
import sys

from extension_worker import api, BoundedCache
import off_shell_beta as off
import stacking_lower as lower
import all_cochain_upper as upper
from coherent_low_commutative import DegreeTwoCommutativeStacking

p = api.p
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
    """A triangular state whose upper formulas are built only when needed."""
    def __init__(self, k, build):
        self.k = k
        self.build = build
        self.layers = {}

    def __getitem__(self, layer):
        if layer not in self.layers:
            self.layers[layer] = self.build(layer)
        return self.layers[layer]


class TransferredModel:
    def __init__(self, setup, transport):
        if setup.get("schema", 1) != 1:
            raise ValueError("unsupported extension transfer schema")
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
        self.ranks = setup["ranks"]
        self.g = setup["g"]
        self.matrices = {False: setup["ordinary"], True: setup["signed"]}
        self.max_degree = len(self.ranks) - 1
        self.max_flags = setup.get("maxFlagSimplices", 8192)
        if type(self.max_flags) is not int or self.max_flags < 1:
            raise ValueError("maxFlagSimplices must be positive")
        if any(type(n) is not int or n < 0 for n in self.ranks):
            raise ValueError("invalid resolution ranks")
        if len(self.g) != len(self.ranks):
            raise ValueError("g must cover every supplied resolution degree")
        for n, rank in enumerate(self.ranks):
            if len(self.g[n]) != rank:
                raise ValueError("g has the wrong resolution rank")
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
        if any(type(g) is not int or not 0 <= g < self.size for g in vertices):
            raise ValueError("transport accepts only group vertex labels")
        if any(a == b for a, b in zip(vertices, vertices[1:])):
            return None
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

    def lift(self, n, vector, signed):
        if n < 0 or not any(vector):
            return p.zero(n)
        weight = 2 if signed else 1
        def value(vertices):
            result = sum(t[weight] * vector[t[0]] for t in self.terms("f", vertices))
            return result if signed else result % 2
        return p.Cochain(n, value)

    def project(self, cochain, signed):
        n = cochain.degree
        if n < 0:
            return []
        weight = 1 if signed else 0
        values = [sum(t[weight] * cochain(tuple(t[2])) for t in chain)
                  for chain in self.g[n]]
        return [integer(v) if signed else integer(v) % 2 for v in values]

    def homotopy(self, cochain, signed):
        if cochain.degree <= 0:
            return p.zero(cochain.degree - 1)
        weight = 1 if signed else 0
        def value(vertices):
            result = sum(t[weight] * cochain(tuple(t[2]))
                         for t in self.terms("h", vertices))
            return result if signed else result % 2
        return p.Cochain(cochain.degree - 1, value)

    @lru_cache(None)
    def simplices(self, n):
        if n < 0:
            return ()
        count = (self.size - 1) ** n
        if count > self.max_flags:
            raise TransferResourceLimit(
                f"exact bar zero test in degree {n} needs {count} simplices "
                f"(limit {self.max_flags})")
        result = []
        for increments in itertools.product(range(1, self.size), repeat=n):
            vertices = [0]
            for g in increments:
                vertices.append(self.mul[vertices[-1]][g])
            result.append(tuple(vertices))
        return tuple(result)

    def is_zero(self, cochain):
        return all(cochain(sigma) == 0 for sigma in self.simplices(cochain.degree))

    def differential(self, cochain, signed):
        if cochain.degree < 0:
            return p.zero(cochain.degree + 1)
        return p.ds(cochain, self.s) if signed else p.binary(p.differential(cochain))

    def legal_pair(self, state):
        if state.k == 2:
            return self.is_zero(self.differential(state[1], False))
        return all(self.is_zero(c) for c in
                   off.curvature(state[0], state[1], self.s, self.omega))

    def triple(self, state, full=True):
        legal = self.legal_pair(state)
        if not full:
            return upper.Triple(state[0], state[1], state[2], legal)
        pure = (self.is_zero(state[0]) and self.is_zero(state[1])
                and self.is_zero(self.differential(state[2], False)))
        complete = legal and self.is_zero(p.binary(self.differential(state[2], False)
            + self.nonlinear(state, 2)))
        return upper.Triple(state[0], state[1], state[2], legal, pure, complete)

    def nonlinear(self, state, layer):
        key = ("nonlinear", layer)
        if key in state.layers:
            return state.layers[key]
        k = state.k
        degree = degrees(k + 1)[layer]
        if layer == 0 or (k == 2 and layer == 1):
            value = p.zero(degree)
        elif k == 2:
            if layer == 2:
                value = p.QD(state[1], self.s, self.omega)
            else:
                value = api.cochain(self.low2.g(
                    *map(api.local, (state[1], state[2], self.s, self.omega)),
                    b_closed=self.is_zero(self.differential(state[1], False))))
        elif k in (3, 4, 5):
            if layer == 1:
                value = off.primary(state[0], self.s, self.omega)
            elif layer == 2:
                value = off.current_f(off.LowerPair(state[0], state[1],
                    self.legal_pair(state)), self.s, self.omega)
            else:
                value = upper.g(self.triple(state, full=False), self.s, self.omega)
        else:
            raise ValueError("transferred nonlinear differential covers degrees 2 through 5")
        state.layers[key] = value
        return value

    def phi(self, k, data):
        checked = self.check_state(k, data)
        d_coordinates = checked["D"]
        checked = dict(checked, D=(0,) * self.dimension(k + 1))
        key = k, tuple(checked[f] for f in "ABC")
        if key not in self._phis:
            def build(layer):
                signed = layer in (0, 3)
                value = self.lift(degrees(k)[layer], checked[FIELDS[layer]], signed)
                if layer:
                    value = value - self.homotopy(self.nonlinear(result, layer), signed)
                return value if signed else p.binary(value)
            result = LazyState(k, build)
            self._phis[key] = result
        base = self._phis[key]
        if not any(d_coordinates):
            return base
        # Every nonlinear term reads ABC only. Reuse its entire cochain DAG
        # while preserving the exact signed integral D lift separately.
        d_lift = self.lift(k + 1, d_coordinates, True)
        return LazyState(k, lambda layer: base[layer] if layer < 3 else base[3] + d_lift)

    def bar_d(self, state):
        def build(layer):
            value = self.differential(state[layer], layer in (0, 3))
            value = value + self.nonlinear(state, layer)
            return value if layer in (0, 3) else p.binary(value)
        return LazyState(state.k + 1, build)

    def bar_product(self, left, right):
        if left.k != right.k or left.k not in (3, 4, 5):
            raise ValueError("transferred products cover equal degrees 3 through 5")
        def cross(layer):
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
            pure = (self.is_zero(result[0]) and self.is_zero(result[1])
                    and self.is_zero(self.differential(result[2], False)))
            return upper.gamma(self.triple(left), self.triple(right),
                               legal, pure, self.s, self.omega)
        def build(layer):
            value = left[layer] + right[layer] + result.cross(layer)
            return value if layer in (0, 3) else p.binary(value)
        result = LazyState(left.k, build)
        result.cross = lru_cache(None)(cross)
        return result

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
        gauge = LazyState(k - 1, lambda layer: p.zero(degrees(k - 1)[layer]))
        # This embedding reads native coordinates only after they have been
        # determined. Each cached layer depends only on earlier coordinates.
        def build(layer):
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
            correction = self.nonlinear(gauge, layer)
            if layer:
                correction = correction - self.homotopy(
                    self.nonlinear(embedded, layer), signed)
            correction = correction + reconstructed.cross(layer)
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
            return answer, LazyState(k - 1, lambda layer: p.zero(degrees(k - 1)[layer]))
        return self.reflect(self.bar_product(self.phi(k, left), self.phi(k, right)), upto)

    def act(self, k, gauge, canonical):
        gauge = self.check_state(k - 1, gauge)
        canonical = self.check_state(k, canonical)
        if not any(any(v) for v in gauge.values()):
            return ({f: list(canonical[f]) for f in FIELDS},
                    LazyState(k - 1, lambda layer: p.zero(degrees(k - 1)[layer])))
        boundary = self.bar_d(self.phi(k - 1, gauge))
        return self.reflect(self.bar_product(boundary, self.phi(k, canonical)))

    def divide_left(self, k, left, total, verify=True):
        left = self.check_state(k, left)
        total = self.check_state(k, total)
        if verify:
            self.require_flat(k, left)
            self.require_flat(k, total)
        right = self.zero(k)
        for layer, f in enumerate(FIELDS):
            product, _ = self.product(k, left, right, upto=layer)
            right[f] = [a + b - c for a, b, c in zip(right[f], total[f], product[f])]
            if layer in (1, 2):
                right[f] = [v % 2 for v in right[f]]
        if verify:
            self.require_flat(k, right)
            product, _ = self.product(k, left, right)
            if any(tuple(product[f]) != total[f] for f in FIELDS):
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
                result, _ = self.act(k, x0, y0)
            else:
                result = self.divide_left(k, x0, y0, verify=False)
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
        if k not in (2, 3, 4, 5) or (k == 2 and operation not in ("d", "phi_values")):
            raise ValueError("transferred states cover degrees 3 through 5; gauges include degree 2")
        key = json.dumps(request, sort_keys=True, separators=(",", ":"))
        if key in self._answers:
            return json.loads(self._answers[key])
        if operation == "d":
            answer = {"state": self.kappa(k, request["state"])}
        elif operation == "xtimes":
            upto = request.get("upto", 3)
            if type(upto) is not int or not 0 <= upto <= 3:
                raise ValueError("upto must be a layer index from zero through three")
            required = 3 if upto == 3 else upto - 1
            self.require_flat(k, request["state"], required)
            self.require_flat(k, request["other"], required)
            result = self.affine_core("xtimes", k, request["state"], request["other"], upto)
            answer = {"state": result}
        elif operation == "act":
            self.require_flat(k, request["state"])
            result = self.affine_core("act", k, request["gauge"], request["state"])
            answer = {"state": result}
        elif operation == "divide_left":
            self.require_flat(k, request["state"])
            self.require_flat(k, request["other"])
            result = self.affine_core("divide_left", k, request["state"], request["other"])
            self.require_flat(k, result)
            product = self.affine_core("xtimes", k, request["state"], result)
            if any(tuple(product[f]) != tuple(request["other"][f]) for f in FIELDS):
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
    def transport(kind, degree, vertices):
        print(json.dumps(dict(operation="transport", kind=kind,
                              degree=degree, vertices=vertices)), flush=True)
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
            answer = {"status": "unresolved" if isinstance(exc, TransferResourceLimit) else "error",
                      "exception": type(exc).__name__, "reason": str(exc)}
        print(json.dumps(answer, separators=(",", ":")), flush=True)


if __name__ == "__main__":
    serve()
