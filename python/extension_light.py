"""Light extension rows of koFull: residues of transported defining data.

The light rows replace the measurement of a relation in the transferred
stacking model by one fixed residue cochain, evaluated on bar defining data
built from cochains on the resolution R and paired back to R along the
comparison chains g(e_j) (doc/extensions.md, "Light rows"). A generic
non-closed defining cochain is the primitive

    P(z; r) = Lambda r + H z,   delta_R r = Pi z  (solved in GAP),

of a closed source z, so delta P(z; r) = z holds literally on the bar
(1 - Lambda Pi = delta H + H delta). Branch flags of the stacking formulas
are those of the construction; nothing here tests a cochain for zero on R.

A request names a task and carries every R vector it needs; the worker keeps
no state between requests besides the transport caches. Rational cochains
whose coboundary is wanted are paired on the chains of their own degree and
the twisted coboundary of R is applied to the pairing (Pi delta_s = delta_s Pi
on normalized invariant cochains), as the page engine does for T.

Primary defining systems instead use Lambda r + K_D, with delta_R r=D_R
and delta K_D=QD Lambda+Lambda D_R. All later tasks use that same marking;
the normalized bar homotopy is absent from these defining cochains.

Answers are always computed in band: a resource limit is returned as
"refused", any other failure as "error"; the stream stays open.

Copyright (c) 2026 koAHSS contributors; MIT license.
"""
from fractions import Fraction
import json
import os
import traceback
import weakref

from extension_worker import api, BoundedCache
import extension_native_upper as native
import all_cochain_upper as upper
import closed_ab_upper as closed
import closed_a_upper as sharp
import off_shell_beta as off
import stacking_lower as lower
import low_phases as low
import h_tau_primitive as hp
from upper_phase_diagnostic import production_phase

p = api.p
TAIL = {4: 5, 5: 4, 6: 4}


def _local(c):
    return native._local(c)


def _times(c, n):
    """An integral multiple with integer values (chi needs integers)."""
    if n == 1:
        return c
    return p.Cochain(c.degree, lambda f, c=c, n=n: n * c(f))


def _product(x, y):
    """The pointwise product of two binary cochains of one degree."""
    return p.Cochain(x.degree, lambda f, x=x, y=y: (x(f) * y(f)) % 2)


def _binom2(n):
    return (n * (n - 1) // 2) % 2


class LightEvaluator:
    def __init__(self, model):
        self.model = model
        self.s, self.w = model.s, model.omega
        # Defining data and their expensive formula cochains, kept with their
        # memo tables between the requests of one relation; keyed by the R
        # vectors that determine them.
        self._memo = BoundedCache(48)
        self._normalized = weakref.WeakKeyDictionary()
        self._audit_normalization = os.environ.get('FERMIONAHSS_LIGHT_NORMALIZATION_CHECKS') == '1'
        self._reduce_transport = os.environ.get('FERMIONAHSS_LIGHT_TRANSPORT_REDUCTION', '1') != '0'

    def memo(self, key, build, objects=()):
        """A cached value; a key naming cochains by identity keeps them alive
        with the value, so the identities stay unique while it is cached."""
        found = self._memo.get(key)
        if found is None or any(x is not y for x, y in zip(found[0], objects)):
            found = (tuple(objects), build())
            self._memo[key] = found
        return found[1]

    # ------------------------------------------------------------------
    # transport and pairings
    def lift(self, n, vector, signed=False):
        # Unsigned lifts already have binary values. Retain native_lift so
        # Pi Lambda can be evaluated without consuming bar cochain values.
        return self.model.lift(n, tuple(vector), signed)

    def prim(self, source, vector, signed=False):
        """P(z; r) for a closed source, checking normalization before using H."""
        return self._primitive(source, vector, signed, check_normalization=True)

    def _prim_normalized(self, source, vector, signed=False):
        """Internal sources normalized by construction; see doc/extensions.md.

        Only the fixed light-task constructors below may use this path.
        Arbitrary Cochain callables must go through prim instead.
        """
        return self._primitive(source, vector, signed,
                               check_normalization=self._audit_normalization)

    def _primitive(self, source, vector, signed, *, check_normalization):
        """P(z; r) = Lambda r + H z, with delta P = z for a normalized closed z."""
        # Preserve the transfer engine's structural-zero rule: H(0)=0.
        # Wrapping it in an opaque cochain would also hide the zero from
        # all the formula builders that consume this primitive.
        if source.degree <= 0 or getattr(source, 'structural_zero', False):
            return self.lift(source.degree - 1, vector, signed)
        homotopy = self._homotopy(source, signed, check_normalization)
        value = self.model.lift(source.degree - 1, tuple(vector), signed) + homotopy
        return value if signed else p.binary(value)

    def _homotopy(self, source, signed=False, check_normalization=False):
        """H on a constructed normalized cochain; the source need not be closed."""
        if source.degree <= 0 or getattr(source, 'structural_zero', False):
            return p.zero(source.degree - 1)
        weight = 1 if signed else 0
        def h_value(vertices):
            # A normalized contraction drops degeneracies. Check the source
            # on the degeneracies of each face it actually consumes, as well
            # as those of the input, before using that contraction.
            if check_normalization:
                self.normalized_degeneracies(source, tuple(vertices))
            terms = self.model.terms(self.model.homotopy_kind(signed), vertices)
            if check_normalization:
                for term in terms:
                    simplex = tuple(term[2])
                    for i in range(len(simplex)):
                        self.normalized_degeneracies(source, simplex[:i] + simplex[i + 1:])
            value = sum(t[weight] * source(tuple(t[2])) for t in terms)
            return value if signed else value % 2
        homotopy = p.Cochain(source.degree - 1, h_value)
        # Retain the transfer engine's checked H g = 0 fast path. When that
        # chain identity holds, no source values are consumed by the pairing.
        homotopy.homotopy_image = signed
        return homotopy

    def b_primitive(self, k, bv, P, vector, C=None, primary_gauge=None):
        """The original P(hD(b,b)+C+P;vector), using a lower-degree carry.

        For the closed b= rho Lambda_s bv, kappa=(Lambda_s bv-b)/2,
        hD(b,b)=Lambda_2 rho beta_s^R bv + delta rho kappa. Thus
        H hD(b,b)=rho kappa-Lambda_2 Pi rho kappa-delta H rho kappa.
        Keeping both correction terms preserves the original marking exactly.
        C is a binary sum of native lifts, so H C=0.
        """
        b = self.lift(k - 2, bv)
        if not self._reduce_transport:
            source = p.hD(b, b, self.s) + P
            if C is not None:
                source = source + C
            return self._prim_normalized(p.binary(source), vector)
        carry, projected, correction, _ = self.b_carry(k, bv)
        base = self.lift(k - 2, [(a + v) % 2 for a, v in zip(vector, projected)])
        if primary_gauge is not None and self.model.primary_comparison:
            # P=QD Lambda y0. Keep its exact H image, using K_D one degree
            # lower: H P=K_D-Lambda Pi K_D-delta H K_D (binary signs).
            K, _ = self.primary_comparison(k - 3, primary_gauge)
            hp_ = p.binary(K + self.lift(k - 2, self.bin(K)) + p.differential(
                self._homotopy(K, check_normalization=self._audit_normalization)))
        else:
            hp_ = self._homotopy(P, check_normalization=self._audit_normalization)
        return p.binary(base + carry + correction + hp_)

    def b_carry(self, k, bv):
        """kappa, Pi kappa, delta H kappa and rho beta_s^R b, cached together."""
        def build():
            q = self.exact([Fraction(v, 2) for v in
                            self.model.coboundary(k - 2, bv, True)], 'B Bockstein')
            lifted = self.lift(k - 2, bv, True)
            b = self.lift(k - 2, bv)
            carry = p.binary(p.divide(lifted - b, 2, 'Bockstein comparison carry'))
            projected = self.bin(carry)
            hc = self._homotopy(carry, check_normalization=self._audit_normalization)
            return carry, projected, p.binary(p.differential(hc)), [v % 2 for v in q]
        return self.memo(('B carry', k, tuple(bv)), build)

    def cup_comparison(self, i, n, a, m, b):
        """K_i(a,b) for the fixed f and the native binary higher diagonals."""
        if i < 0 or not any(a) or not any(b):
            return p.zero(n + m - i - 1)
        def value(vertices):
            return sum(a[t[2]] * b[t[4]] for t in
                       self.model.terms('cupK' + str(i), vertices)
                       if t[1] == n and t[3] == m) % 2
        return p.Cochain(n + m - i - 1, value)

    def native_cup(self, i, n, a, m, b):
        if i < 0 or not any(a) or not any(b):
            return [0] * self.model.dimension(n + m - i)
        answer = self.model.transport('cupR', i, dict(degrees=[n, m], vectors=[a, b]))
        if answer.get('status') != 'computed':
            from extension_transfer import TransferResourceLimit
            raise TransferResourceLimit(answer.get('reason', 'native cup refused'))
        return answer['terms']

    def primary_comparison(self, n, vector):
        """K_D and D_R using cup Sq1, with delta K_D=QD Lambda+Lambda D_R.

        Both parts of the defining system use this cup representative;
        neither substitutes the Bockstein representative of nativePrimary.
        """
        v = [x % 2 for x in vector]
        def build():
            if any(x % 2 for x in self.model.coboundary(n, v, False)):
                raise ArithmeticError('primary comparison requires a native cocycle')
            sv, wv = self.model.s_vector, self.model.w_vector
            K = self.cup_comparison(n - 2, n, v, n, v) + self.cup_comparison(0, 2, wv, n, v)
            D = [a + b for a, b in zip(self.native_cup(n - 2, n, v, n, v),
                                       self.native_cup(0, 2, wv, n, v))]
            if any(sv):
                sq1 = self.native_cup(n - 1, n, v, n, v)
                K = K + self.cup_comparison(0, 1, sv, n + 1, sq1) + p.cup(
                    self.s, self.cup_comparison(n - 1, n, v, n, v))
                D = [a + b for a, b in zip(D, self.native_cup(0, 1, sv, n + 1, sq1))]
            K = p.binary(K)
            return K, [x % 2 for x in D]
        return self.memo(('primary comparison', n, tuple(v)), build)

    def qd_source(self, n, vector):
        """The defining right-hand side, native D_R or the generic Pi QD."""
        if not self._reduce_transport or not self.model.primary_comparison:
            return self.bin(p.QD(self.lift(n, vector), self.s, self.w))
        _, D = self.primary_comparison(n, vector)
        return D

    def qd_primitive(self, n, vector, primitive):
        """Construct a primary defining cochain with no normalized bar H.

        The native vector solves delta primitive=D_R vector. This is a
        coherent choice of marking, used by every atom, reference and gauge;
        it is not the old P with a homotopy correction simply discarded.
        """
        if not self._reduce_transport or not self.model.primary_comparison:
            return self._prim_normalized(p.QD(self.lift(n, vector), self.s, self.w), primitive)
        K, _ = self.primary_comparison(n, vector)
        return p.binary(self.lift(n + 1, primitive) + K)

    def check_normalized(self, c, vertices):
        checked = self._normalized.get(c)
        if checked is None:
            checked = BoundedCache(4096)
            self._normalized[c] = checked
        if vertices not in checked:
            if c(vertices) != 0:
                raise ArithmeticError(f'light cochain of degree {c.degree} is not normalized at {vertices}')
            checked[vertices] = True

    def normalized_degeneracies(self, c, vertices):
        """Check c on degeneracies of a simplex one degree below it."""
        for i in range(len(vertices)):
            self.check_normalized(c, vertices[:i] + (vertices[i],) + vertices[i:])

    def bin(self, c):
        # Projection already reduces modulo two; an opaque binary wrapper
        # would hide native lifts and homotopy images from the SDR identities.
        return [int(v) % 2 for v in self.model.project(c, False)]

    def integer(self, c):
        return [int(v) for v in self.model.project(c, True)]

    def rational(self, c):
        return [Fraction(v) for v in self.model._pairing(c, 1)]

    def rational_ds(self, c):
        """Pair a rational coboundary, checking omitted degenerate faces."""
        return self._rational_coboundary(c, check_normalization=True)

    def _rational_ds_normalized(self, c):
        """Internal potential built from normalized lifts and formula constructors."""
        return self._rational_coboundary(c, check_normalization=self._audit_normalization)

    def _rational_coboundary(self, c, *, check_normalization):
        """Pi(delta_s c) = delta_s^R Pi(c), from the chains of the degree of c."""
        # The normalized chain boundary omits these faces. A nonzero value
        # would invalidate the chain-map identity rather than change a residue
        # by an allowed integral lift.
        if check_normalization:
            for j in range(self.model.dimension(c.degree + 1)):
                for term in self.model.chain(c.degree + 1, j):
                    simplex = tuple(term[2])
                    for i in range(1, len(simplex) - 1):
                        if simplex[i - 1] == simplex[i + 1]:
                            self.check_normalized(c, simplex[:i] + simplex[i + 1:])
        return [Fraction(v) for v in self.model.coboundary(c.degree, self.rational(c), True)]

    @staticmethod
    def exact(values, label):
        answer = []
        for v in values:
            v = Fraction(v)
            if v.denominator != 1:
                raise ArithmeticError(f'{label} is not integral on R: {v}')
            answer.append(int(v))
        return answer

    @staticmethod
    def modular(values, m, prime, label):
        answer = []
        for v in values:
            v = Fraction(v)
            if v.denominator % prime == 0:
                raise ArithmeticError(f'{label} has a denominator divisible by {prime}: {v}')
            answer.append(v.numerator * pow(v.denominator, -1, m) % m)
        return answer

    # ------------------------------------------------------------------
    # defining data
    def atom(self, k, d):
        def build():
            b = self.lift(k - 2, d['b'])
            c = self.qd_primitive(k - 2, d['b'], d['c']) if 'c' in d else None
            return b, c
        return self.memo(('atom', k, tuple(d['b']), tuple(d.get('c', ()))), build)

    def combination(self, k, terms):
        """The lower part (B, C) of an integer combination of A=0 atoms (b, c),
        by the lower product law: B = sum e b, C = sum e c + sum binom(n,2) hD(b,b)
        + sum_{i<j} e_i e_j hD(b_i,b_j), e = n mod 2."""
        s = self.s
        B, C = p.zero(k - 2), p.zero(k - 1)
        odd = []
        for n, bv, cv in terms:
            b, c = self.atom(k, {'b': bv, 'c': cv})
            if n % 2:
                B, C = B + b, C + c
                odd.append(b)
            if _binom2(n):
                C = C + p.hD(b, b, s)
        for i in range(len(odd)):
            for j in range(i + 1, len(odd)):
                C = C + p.hD(odd[i], odd[j], s)
        return p.binary(B), p.binary(C)

    def pure_c(self, k, vectors):
        """C = sum C_j and N_C/2 of the pure-C reference (light-transport-proof (30))."""
        Cs = [self.lift(k - 1, v) for v in vectors]
        C = p.zero(k - 1)
        for x in Cs:
            C = C + x
        C = p.binary(C)
        if len(Cs) < 2:
            return C, p.zero(k + 1)
        P2 = p.zero(k)
        for i in range(len(Cs)):
            for j in range(i + 1, len(Cs)):
                P2 = P2 + p.polarization(Cs[i], Cs[j], 2)
        P2 = p.binary(P2)
        total = p.zero(k + 1)
        for x in Cs:
            total = total + p.E(x, self.w)
        total = total - p.E(C, self.w) + p.ds(P2, self.s)
        return C, p.divide(total, 2, 'pure-C reference carry')

    def reference(self, k, ref):
        """(B0, C0, -N_C/2, pure) of a marked lower reference."""
        return self.memo(('ref', k, json.dumps(ref, sort_keys=True)),
                         lambda: self._reference(k, ref))

    def _reference(self, k, ref):
        kind = ref['kind']
        if kind == 'zero':
            return p.zero(k - 2), p.zero(k - 1), p.zero(k + 1), True
        if kind == 'atom':
            b, c = self.atom(k, ref)
            return b, c, p.zero(k + 1), False
        if kind == 'pureC':
            C, half = self.pure_c(k, ref['C'])
            return p.zero(k - 2), C, p.scale(half, -1), True
        if kind == 'lower':
            B, C = self.combination(k, ref['terms'])
            return B, C, p.zero(k + 1), False
        raise ValueError('unknown light reference')

    def b_gauge(self, k, g):
        """(u, y, P) of a B relation gauge: delta_s u = 0, delta y = Q_D(rho u), P = f#(u, y)."""
        s, w = self.s, self.w
        if g.get('u'):
            u = self.lift(k - 4, g['u'], True)
            y = self.qd_primitive(k - 4, g['u'], g['yR'])
            if g.get('y1'):
                y = p.binary(y + self.lift(k - 3, g['y1']))
            return u, y, sharp.fsharp(u, y, s, w), 'tau'
        u = p.zero(k - 4)
        if g.get('y0'):
            y = self.lift(k - 3, g['y0'])
            return u, y, p.QD(y, s, w), 'D'
        return u, p.zero(k - 3), p.zero(k - 1), 'pure'

    def gauge_phase(self, k, u, y, pi, kind):
        """The actual gauge phase Omega_{k-1}(u,y,pi) of the page form (R0, R1, R2)."""
        if k == 6 and kind != 'tau':
            rule = native._a0(5)
            return p.Cochain(6, rule.omega(_local(y), _local(pi), _local(self.s), _local(self.w)))
        return production_phase(k - 4, u, y, pi, self.s, self.w)

    def gauge_potential(self, k1, u, y, pi, P, t, kind):
        """The potential of the gauge curvature J_{k1}(u,y,pi) = delta_s(pot) - h(E t), k1 = 4, 5."""
        s, w = self.s, self.w
        if k1 == 5 and kind != 'tau':
            if kind == 'pure':
                return p.half(p.E(pi, w))
            rule = native._a0(5)
            return (p.Cochain(6, rule.omega(_local(y), _local(pi), _local(s), _local(w)))
                    + p.half(p.binary(p.cup(t, P, 4))))
        return closed.potential(u, y, pi, s, w)

    def a0_potential(self, k, b, c):
        """The A=0 potential with J_k(0,b,c) = delta_s of it (C curvature zero), k = 4,5,6."""
        def build():
            if k == 4:
                return closed.potential(p.zero(1), b, c, self.s, self.w)
            rule = native._a0(k)
            return p.Cochain(k + 1, rule.omega(_local(b), _local(c), _local(self.s), _local(self.w)))
        return self.memo(('a0', k, id(b), id(c)), build, (b, c))

    # ------------------------------------------------------------------
    # A lower systems and relation gauges
    def a_key(self, k, d):
        return ('a', k, tuple(d['A']), tuple(d.get('BR', ())), tuple(d.get('CR', ())),
                json.dumps(d.get('star'), sort_keys=True))

    def a_data(self, k, d):
        """(A, a, B, C) with delta B = Q_D(a), delta C = f#(A, B), and the
        admissible correction by the A=0 combination star when given."""
        return self.memo(self.a_key(k, d), lambda: self._a_data(k, d))

    def _a_data(self, k, d):
        s, w = self.s, self.w
        A = self.lift(k - 3, d['A'], True)
        a = p.binary(A)
        B = self.qd_primitive(k - 3, d['A'], d['BR']) if 'BR' in d else None
        C = self._prim_normalized(sharp.fsharp(A, B, s, w), d['CR']) if 'CR' in d else None
        star = d.get('star')
        if star:
            bs, cs = self.combination(k, star.get('terms', []))
            for z in star.get('z', []):
                cs = p.binary(cs + self.lift(k - 1, z))
            beta = lower.legal_beta(A, B, p.zero(k - 3), bs, s, w)
            if C is not None:
                C = p.binary(C + cs + beta)
            B = p.binary(B + bs)
        return A, a, B, C

    def x_power(self, k, A, a, B, m, e):
        """(X_B, X_C) of the m-th lower power (light-transport-proof (17))."""
        s, w = self.s, self.w
        def build():
            alpha = hp.hD(a, a, s)
            if e == 1:
                return alpha, p.binary(lower.legal_beta(A, B, A, B, s, w))
            half = _times(A, m // 2)
            previous = alpha if e == 2 else p.zero(k - 2)
            return p.zero(k - 2), p.binary(lower.legal_beta(half, previous, half, previous, s, w))
        return self.memo(('x', k, id(A), id(B), m), build, (A, B))

    def a_gauge(self, k, d, A, a, ref):
        """U, the Y source and Y of an A relation gauge (relative note (4.2), (11.8))."""
        s, w = self.s, self.w
        m, e = d['m'], d['e']
        g = d['U']
        # delta_s u=m A_R and H Lambda_s=0, literally for the fixed SDR.
        U = (self.lift(k - 4, g['u'], True) if self._reduce_transport else
             self._prim_normalized(_times(A, m), g['u'], True))
        if g.get('v'):
            U = U + self.lift(k - 4, g['v'], True)
        XB = hp.hD(a, a, s) if e == 1 else p.zero(k - 2)
        B0 = ref[0]
        source = p.binary(p.QD(p.binary(U), s, w) + XB + B0)
        native_source = None
        if self._reduce_transport and self.model.primary_comparison:
            uv = list(g['u'])
            if g.get('v'):
                uv = [x + y for x, y in zip(uv, g['v'])]
            K, D = self.primary_comparison(k - 4, uv)
            carry, q = p.zero(k - 3), [0] * self.model.dimension(k - 2)
            if e == 1:
                av = [x % 2 for x in d['A']]
                carry = p.binary(p.divide(self.lift(k - 3, av, True) - a, 2,
                                          'A polarization comparison carry'))
                q = [x % 2 for x in self.exact([Fraction(x, 2) for x in
                     self.model.coboundary(k - 3, av, True)], 'A polarization Bockstein')]
            native_source = [(x + y + z) % 2 for x, y, z in zip(D, q, self.bin(B0))]
            Y = p.binary(self.lift(k - 3, g['YR']) + K + carry) if 'YR' in g else None
        else:
            Y = self._prim_normalized(source, g['YR']) if 'YR' in g else None
        if Y is not None and g.get('tau'):
            v = self.lift(k - 4, g['tau']['v'], True)
            yv = self.qd_primitive(k - 4, g['tau']['v'], g['tau']['yR'])
            Y = p.binary(Y + yv + hp.hD(p.binary(U), p.binary(v), s))
            U = U + v
        if Y is not None and g.get('yD'):
            Y = p.binary(Y + self.lift(k - 3, g['yD']))
        return U, Y, XB, source, native_source

    def r_c(self, k, A, a, B, U, Y, ref, m, e):
        """R_C = X_C + H_{k-1}(U,Y) + beta(mA, X_B+B0; 0, B0) + C0 (relative note (4.3))."""
        s, w = self.s, self.w
        XB, XC = self.x_power(k, A, a, B, m, e)
        B0, C0 = ref[0], ref[1]
        beta = self.memo(('betaL', k, id(A), id(XB), id(B0), m), lambda: p.binary(
            lower.legal_beta(_times(A, m), p.binary(XB + B0), p.zero(k - 3), B0, s, w)),
            (A, XB, B0))
        return p.binary(XC + off.raw_H(U, Y, s, w) + beta + C0)

    def diagonal_phase(self, k, A, B, C):
        if k == 6:
            import unary_gamma6
            return unary_gamma6.diagonal_phase(A, B, C, self.s, self.w)
        return closed.phase(A, B, C, A, B, C, self.s, self.w)

    def a_residue(self, k, d, A, a, B, C, ref):
        """rho_m Pi Psi_rel without the reference D component (light-transport-proof
        (14) with the bounded tail (25)); the D marking is added in GAP."""
        s, w = self.s, self.w
        m, e = d['m'], d['e']
        U, Y, _, _, _ = self.a_gauge(k, d, A, a, ref)
        R = self.r_c(k, A, a, B, U, Y, ref, m, e)
        W = self._prim_normalized(R, d['U']['WR'])
        B0, C0, carry, pure = ref
        z = [p.zero(j) for j in range(k + 2)]
        # degree k+1 terms
        omega_a = self.memo(('pot', id(A), id(B), id(C)), lambda: closed.potential(A, B, C, s, w),
                            (A, B, C))
        omega_0 = production_phase(k - 3, z[k - 3], z[k - 2], z[k - 1], s, w)
        if k == 4:
            omega_ref = production_phase(1, z[1], B0, C0, s, w)
        else:
            rule = native._a0(k)
            omega_ref = p.Cochain(k + 1, rule.omega(_local(B0), _local(C0), _local(s), _local(w)))
        direct = _times(omega_a, m) - omega_0 - omega_ref + carry
        # degree k terms, paired before the coboundary of R
        r = upper.r(upper.Triple(z[k - 3], B0, C0, True, pure), s)
        Ui, Yi, Wi = hp.interval(U, True), hp.interval(Y, True), hp.interval(W, True)
        si, wi = hp.interval(s), hp.interval(w)
        acyl, bcyl = off.curvature(Ui, Yi, si, wi)
        ccyl = p.binary(p.differential(Wi) + off.raw_H(Ui, Yi, si, wi))
        cylinder = hp.prism(production_phase(k - 3, acyl, bcyl, ccyl, si, wi))
        BF = p.binary(p.differential(Y) + p.QD(p.binary(U), s, w))
        CF = p.binary(p.differential(W) + off.raw_H(U, Y, s, w))
        mA = _times(A, m)
        cross = closed.phase(mA, BF, CF, z[k - 3], B0, C0, s, w)
        tail = p.zero(k)
        alpha = hp.hD(a, a, s)
        for i in range(max(0, e - TAIL[k]), e):
            if i == 0:
                x = (A, B, C)
            elif i == 1:
                x = (_times(A, 2), alpha, p.binary(lower.legal_beta(A, B, A, B, s, w)))
            else:
                h = _times(A, 2 ** (i - 1))
                previous = alpha if i == 2 else z[k - 2]
                x = (_times(A, 2 ** i), z[k - 2],
                     p.binary(lower.legal_beta(h, previous, h, previous, s, w)))
            phase = self.memo(('diag', k, id(A), id(B), id(C), i),
                              lambda x=x: self.diagonal_phase(k, *x), (A, B, C))
            tail = tail + _times(phase, 2 ** (e - 1 - i))
        delta = r - cylinder - cross + tail
        values = [x + y for x, y in zip(self.rational(direct), self._rational_ds_normalized(delta))]
        return self.modular(values, m, 2, 'the A-over-D residue')

    # ------------------------------------------------------------------
    # tasks
    def task_qd(self, k, d):
        return {'QD': self.qd_source(k - 2, d['b'])}

    def task_atom_curvature(self, k, d):
        """Pi J_k(0,b,c) in degree k+2 (integral)."""
        b, c = self.atom(k, d)
        s, w = self.s, self.w
        if k >= 4:
            return {'J': self.exact(self._rational_ds_normalized(self.a0_potential(k, b, c)), 'A=0 curvature')}
        z = [p.zero(j) for j in range(k + 2)]
        if k == 3:
            J = native.g(3, upper.Triple(z[0], b, c, True, False, True), s, w)
        else:
            J = api.cochain(self.model.low2.g(*map(api.local, (b, c, s, w)), b_closed=True))
        return {'J': self.integer(J)}

    def task_c_mark(self, k, d):
        """The unchanged pure-C curvature and the cohomology class of its square.

        gamma is consumed only in the closed row 2D+gamma. Replacing it by
        e+beta_s Sq1(c), e=Pi tilde E(Lambda c), changes that row by an
        integral coboundary and leaves the chosen D completion unchanged.
        """
        s, w = self.s, self.w
        C = self.lift(k - 1, d['c'])
        if self._reduce_transport:
            e = self.integer(p.E(C, w))
            curvature = self.model.coboundary(k + 1, e, True)
            sq1 = self.exact([Fraction(x, 2) for x in
                             self.model.coboundary(k - 1, d['c'], False)], 'Sq1 lift')
            beta = self.model.coboundary(k, [x % 2 for x in sq1], True)
            return {'J': self.exact([Fraction(x, 2) for x in curvature], 'pure-C curvature'),
                    'gamma': self.exact([x + Fraction(y, 2) for x, y in zip(e, beta)],
                                        'pure-C square class')}
        z = [p.zero(j) for j in range(k + 2)]
        if k >= 3:
            T = upper.Triple(z[k - 3], z[k - 2], C, True, True)
            J = native.g(k, T, s, w, a_zero=True, b_zero=True)
            gamma = native.gamma(k, T, T, True, True, (True, True), s, w, b_zero=True)
        elif k == 2:
            L = api.local
            J = api.cochain(self.model.low2.g(L(z[0]), L(C), L(s), L(w), b_closed=True))
            gamma = api.cochain(self.model.low2.gamma(L(z[0]), L(C), L(z[0]), L(C), L(s), L(w),
                                                       b_closed=True, bp_closed=True, sum_closed=True))
        else:
            L = api.local
            J = api.cochain(self.model.low1.g(L(C), L(s), L(w)))
            gamma = api.cochain(self.model.low1.gamma(L(C), L(C), L(s), L(w)))
        return {'J': self.integer(J), 'gamma': self.integer(gamma)}

    def task_gauge(self, k, d):
        """Projections that select a B relation gauge (u, y, pi)."""
        s, w = self.s, self.w
        if self._reduce_transport:
            _, projected, _, q = self.b_carry(k, d['b'])
            base = [x + y for x, y in zip(q, self.model.coboundary(k - 2, projected, False))]
            for v in d.get('Cref', []):
                base = [x + y for x, y in zip(base, v)]
            base = [x % 2 for x in base]
        else:
            b, _ = self.atom(k, d)
            C, _ = self.pure_c(k, d.get('Cref', []))
            base = self.bin(p.binary(p.hD(b, b, s) + C))
        g = d.get('gauge', {})
        out = {}
        want = d['want']
        if 'xC' in want:
            out['xC'] = base
        if 'QDu' in want:
            out['QDu'] = self.qd_source(k - 4, g['u'])
        if 'xCp' in want or 'target' in want:
            _, _, P, _ = self.b_gauge(k, g)
            out['xCp' if 'xCp' in want else 'target'] = [(x + y) % 2 for x, y in
                                                       zip(base, self.bin(P))]
        return out

    def task_bd_page(self, k, d):
        """The page form Phi_P(b,P;c,pi) of a kernel B atom (light-completion-tau (9))."""
        s, w = self.s, self.w
        b, c = self.atom(k, d)
        g = d['gauge']
        u, y, P, kind = self.b_gauge(k, g)
        pi = self.b_primitive(k, d['b'], P, g['pi'], primary_gauge=g.get('y0'))
        n = b.degree
        A = p.Cochain(n, lambda f, b=b, pi=pi: b(f) % 2 + 2 * (pi(f) % 2))
        if kind == 'pure':
            return {'Phi': self.bin(low.secondary(A, c, s, w))}
        phase = self.gauge_phase(k, u, y, pi, kind) - p.half(p.E(pi, w))
        S = p.source_splitting(A, s, w)
        F = p.binary(p.E(c, w) + S['FA'])
        G = p.binary(S['g'] + p.cup(s, c) + p.cup(S['e'], P, n))
        q = p.cup(w, S['B'], integral=True) + p.cup(S['B'], S['B'], n - 1, integral=True)
        L = p.integral(p.scale(q - p.ds(G, s), Fraction(1, 2)) - p.ds(phase, s), 'page-form carry')
        Phi = p.binary(F + p.binary(L) + p.cup(p.cup(p.cup(s, s), s), S['a'])
                       + p.cup(p.binary(p.cup(s, s) + w), P))
        return {'Phi': self.bin(Phi)}

    def task_bd_exact(self, k, d):
        """Pi of the exact B-over-D residue gamma(b^,b^) - J_{k-1}(u,y,pi)
        - gamma((0,0,t),(0,0,C)) - N_C/2 (light-transport-proof (11))."""
        s, w = self.s, self.w
        b, c = self.atom(k, d)
        g = d['gauge']
        u, y, P, kind = self.b_gauge(k, g) if k >= 3 else (None, None, p.zero(k - 1), 'pure')
        x = p.hD(b, b, s)
        C, half = self.pure_c(k, d.get('Cref', []))
        pi = self.b_primitive(k, d['b'], P, g['pi'], C, primary_gauge=g.get('y0'))
        t = p.binary(p.differential(pi) + P)
        z = [p.zero(j) for j in range(k + 2)]
        L = api.local
        if k <= 3:
            if k == 3:
                Tb = upper.Triple(z[0], b, c, True, False, True)
                gamma = native.gamma(3, Tb, Tb, True, True, (True, True), s, w)
                J = api.cochain(self.model.low2.g(L(y), L(pi), L(s), L(w), b_closed=True))
                unit = native.gamma(3, upper.Triple(z[0], z[1], t, True, True),
                                    upper.Triple(z[0], z[1], C, True, True),
                                    True, True, (True, True), s, w)
            else:
                gamma = api.cochain(self.model.low2.gamma(L(b), L(c), L(b), L(c), L(s), L(w),
                                                           b_closed=True, bp_closed=True, sum_closed=True))
                J = api.cochain(self.model.low1.g(L(pi), L(s), L(w)))
                unit = api.cochain(self.model.low2.gamma(L(z[0]), L(t), L(z[0]), L(C), L(s), L(w),
                                                          b_closed=True, bp_closed=True, sum_closed=True))
            return {'Z': self.integer(gamma - J - unit - half)}
        unit = native.gamma(k, upper.Triple(z[k - 3], z[k - 2], t, True, True),
                            upper.Triple(z[k - 3], z[k - 2], C, True, True),
                            True, True, (True, True), s, w, b_zero=True)
        omega_b = self.a0_potential(k, b, c)
        omega_x = self.a0_potential(k, z[k - 2], x)
        phase = self.memo(('eb', k, id(b), id(c)),
                          lambda: closed.phase(z[k - 3], b, c, z[k - 3], b, c, s, w), (b, c))
        rx = p.half(p.binary(p.cup(s, x)))
        direct = _times(omega_b, 2) - omega_x - unit - half
        delta = phase + rx
        if k == 4:
            direct = direct - native.g(3, upper.Triple(u, y, pi, True), s, w)
        else:
            direct = direct + p.half(p.E(t, w))
            delta = delta - self.gauge_potential(k - 1, u, y, pi, P, t, kind)
        values = [x_ + y_ for x_, y_ in zip(self.rational(direct), self._rational_ds_normalized(delta))]
        return {'Z': self.exact(values, 'the exact B-over-D residue')}

    def task_a_step(self, k, d):
        """Projections of an A relation: its lower system, gauge, C cocycle and residue."""
        s, w = self.s, self.w
        want = d['want']
        A, a, B, C = self.a_data(k, d)
        out = {}
        if 'QDa' in want:
            out['QDa'] = self.qd_source(k - 3, d['A'])
        if 'fsharp' in want:
            out['fsharp'] = self.bin(sharp.fsharp(A, B, s, w))
        if 'pot' in want:
            potential = self.memo(('pot', id(A), id(B), id(C)),
                                  lambda: closed.potential(A, B, C, s, w), (A, B, C))
            out['pot'] = self.exact(self._rational_ds_normalized(potential), 'the A curvature')
        ref = self.reference(k, d['ref']) if 'ref' in d else None
        if 'Ysrc' in want or 'RC' in want or 'QDv' in want:
            U, Y, XB, source, native_source = self.a_gauge(k, d, A, a, ref)
            if 'Ysrc' in want:
                out['Ysrc'] = self.bin(source) if native_source is None else native_source
            if 'QDv' in want:
                out['QDv'] = self.qd_source(k - 4, d['U']['tau']['v'])
            if 'RC' in want:
                out['RC'] = self.bin(self.r_c(k, A, a, B, U, Y, ref, d['m'], d['e']))
        if 'ypp' in want or 'page' in want:
            U = (self.lift(k - 4, d['U']['u'], True) if self._reduce_transport else
                 self._prim_normalized(_times(A, d['m']), d['U']['u'], True))
            if d['U'].get('v'):
                U = U + self.lift(k - 4, d['U']['v'], True)
            uv = list(d['U']['u'])
            if d['U'].get('v'):
                uv = [x + y for x, y in zip(uv, d['U']['v'])]
            if 'ypp' in want:
                out['ypp'] = self.qd_source(k - 4, uv)
            if 'page' in want:
                ypp = self.qd_primitive(k - 4, uv, d['U']['ypp'])
                t = p.binary(p.divide(A - a, 2, 'A carry'))
                ea = p.binary(p.divide(p.differential(a), 2, 'binary Bockstein'))
                st = p.cup(s, t)
                s2 = p.cup(s, s)
                page = (p.hD(B, B, s) + _product(p.E(a, w), p.cup(s, ea)) + p.hD(st, st, s)
                        + p.cup(s2, t) + _product(p.binary(p.differential(st)), p.cup(s2, a)))
                out['page'] = self.bin(p.binary(off.raw_H(U, ypp, s, w) + page + p.cup(w, a)))
        if 'thmC' in want:
            out['thmC'] = self.bin(p.hD(B, B, s))
        if 'residue' in want:
            out['residue'] = self.a_residue(k, d, A, a, B, C, ref)
        return out

    def task_p3_power(self, k, d):
        """Y with rho_3 Y = P^1_s rho_3 A: the cube in degree five, P^1 mod 3 in degree six."""
        s = self.s
        A = self.lift(k - 3, d['A'], True)
        if k == 5:
            return {'Y': self.integer(low.transported(low.transported(A, A, s), A, s))}
        import mod3_power
        return {'P1': [v % 3 for v in self.integer(mod3_power.reduced_power_1(A, s))]}

    TASKS = {'qd': task_qd, 'atom_curvature': task_atom_curvature, 'c_mark': task_c_mark,
             'gauge': task_gauge, 'bd_page': task_bd_page, 'bd_exact': task_bd_exact,
             'a_step': task_a_step, 'p3_power': task_p3_power}

    def calculate(self, request):
        """One light task; never raises (resource limits and failures in band)."""
        from extension_transfer import TransferResourceLimit
        task = request.get('task')
        try:
            if task not in self.TASKS:
                raise ValueError(f'unknown light task {task!r}')
            result = self.TASKS[task](self, request['degree'], request.get('data', {}))
            return {'status': 'computed', 'task': task, 'result': result}
        except (MemoryError, TransferResourceLimit, native.NativeDegreeSixLimit) as exc:
            return {'status': 'computed', 'task': task,
                    'refused': {'code': 'resource-limit', 'exception': type(exc).__name__,
                                'reason': str(exc)}}
        except Exception as exc:
            return {'status': 'computed', 'task': task,
                    'error': {'exception': type(exc).__name__, 'reason': str(exc),
                              'traceback': traceback.format_exc()[-1500:]}}
