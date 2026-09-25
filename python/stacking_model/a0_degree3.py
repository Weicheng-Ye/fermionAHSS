"""Strict A=0 degree-three stacking on arbitrary B,C,D cochains.

The differential is unchanged. Fixed universal contractors supply the
closed-B phase gauge, and an explicit successor phase rephasing accounts
for the global branch. Use successor_gamma for its compatible k4 product;
do not silently replace that successor by a differently calibrated gamma4.
"""
from fractions import Fraction as F
from functools import lru_cache

from cochains import Cochain, binary_sum, cup, differential, signed_differential, zero
from compatible_sector import E, QD, alpha, integral, polarization
from lower_stacking import pullback_interval, prism, raw_H_zero_A
from a0_gamma import omega3
from a0_degree4 import DegreeFourA0Stacking


TREE = ("b", ("s", "omega"))
SOURCE_BASES = ((1, (0, (2,))), (1, (1, (1,))), (1, (3, ())),
                (2, (0, (1,))), (2, (2, ())), (3, (1, ())), (4, (0, ())))


def lower_d(b, c, s, omega):
    return differential(b).mod2(), binary_sum(differential(c).mod2(), raw_H_zero_A(b, s, omega))


def beta(b, bp, s, omega):
    """Natural all-cochain beta3; equals h_D on closed B,B'."""
    bi, bpi = pullback_interval(b, True), pullback_interval(bp, True)
    si, wi = pullback_interval(s), pullback_interval(omega)
    residual = binary_sum(raw_H_zero_A(binary_sum(bi, bpi), si, wi),
                          raw_H_zero_A(bi, si, wi), raw_H_zero_A(bpi, si, wi),
                          alpha(differential(bi).mod2(), differential(bpi).mod2(), si))
    return prism(residual).mod2()


class DegreeThreeA0Stacking:
    def __init__(self, production_root):
        self.successor = DegreeFourA0Stacking(production_root)
        self.chains = self.successor.chains
        self.calibration = tuple(self.chain_value(self.inclusion(basis)) for basis in SOURCE_BASES)
        if any(value % 1 for value in self.calibration):
            raise ArithmeticError("the adjusted degree-three source is no longer zero modulo1")

    def raw_potential(self, b, c, s, omega):
        bi, ci = pullback_interval(b, True), pullback_interval(c, True)
        si, wi = pullback_interval(s), pullback_interval(omega)
        db, dc = lower_d(bi, ci, si, wi)
        return prism(self.successor.omega(db, dc, si, wi))

    def legal_pair_potential(self, b, c, s, omega):
        tau = QD(b, s, omega)
        t = binary_sum(differential(c).mod2(), tau)
        return omega3(b, c, s, omega)+F(1, 2)*cup(t, tau, 2)

    def potential(self, b, c, s, omega, *, b_closed):
        return (self.legal_pair_potential(b, c, s, omega) if b_closed
                else self.raw_potential(b, c, s, omega))

    def adjusted_difference(self, b, c, s, omega):
        """Closed phase Psi on closed B, with completely unrestricted C."""
        t = binary_sum(differential(c).mod2(), QD(b, s, omega))
        return (self.legal_pair_potential(b, c, s, omega)-self.raw_potential(b, c, s, omega)
                -F(1, 2)*cup(s, t))

    def _tensor_maps(self, pair, left, right):
        out = {}
        for a, coefficient in left(pair[0]).items():
            for b, weight in right(pair[1]).items():
                self.chains.add(out, {(a, b): coefficient*weight})
        return out

    @lru_cache(None)
    def inclusion(self, basis, tree=TREE):
        cm = self.chains
        if isinstance(tree, str):
            if tree == "omega":
                return cm.G("c2", basis)
            if tree == "b":
                return {cm.Atom("c2", 0, (1,)*basis): 1}
            return {cm.Diag("zsign", tuple((1, (0,)*basis) for _ in range(basis))): 1}
        return cm.linear(self._tensor_maps(basis, lambda x: self.inclusion(x, tree[0]),
                                          lambda x: self.inclusion(x, tree[1])), cm.product_shuffle)

    @lru_cache(None)
    def forward(self, simplex, tree=TREE):
        cm = self.chains
        if isinstance(tree, str):
            if tree == "omega":
                return cm.F(simplex)
            return {} if cm._degen_indices(simplex) else {cm._degree(simplex): 1}
        return cm.linear(cm.product_AW(simplex), lambda pair: self._tensor_maps(
            pair, lambda x: self.forward(x, tree[0]), lambda x: self.forward(x, tree[1])))

    @lru_cache(None)
    def homotopy(self, simplex, tree=TREE):
        cm = self.chains
        if isinstance(tree, str):
            return cm.H(simplex) if tree == "omega" else {}
        out = dict(cm.product_homotopy(simplex))
        for (a, b), coefficient in cm.product_AW(simplex).items():
            tensor = {(aa, b): c for aa, c in self.homotopy(a, tree[0]).items()}
            projected = cm.linear(self.forward(a, tree[0]), lambda key: self.inclusion(key, tree[0]))
            for aa, c in projected.items():
                for bb, d in self.homotopy(b, tree[1]).items():
                    cm.add(tensor, {(aa, bb): (-1)**cm._degree(a)*c*d})
            cm.add(out, cm.linear(tensor, cm.product_shuffle), coefficient)
        return out

    @lru_cache(None)
    def source_value(self, simplex):
        bb, (ss, ww) = simplex
        b = Cochain(1, lambda face: sum(bb.edges[face[0]:face[1]]) % 2)
        s = Cochain(1, lambda face: sum(row[0] for row in ss.rows[face[0]:face[1]]) % 2)
        omega = Cochain(2, lambda face: sum(ww.rows[i][1][j]
                        for i in range(face[0], face[1]) for j in range(face[1], face[2])) % 2)
        return self.adjusted_difference(b, zero(2), s, omega)(tuple(range(5)))

    def chain_value(self, chain):
        return sum(coefficient*self.source_value(simplex) for simplex, coefficient in chain.items())

    def section(self, b, s, omega, vertices):
        cm, n = self.chains, len(vertices)-1
        matrix = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                value = omega((vertices[i], vertices[i+1], vertices[j+1]))
                if j > i+1:
                    value ^= omega((vertices[i], vertices[i+1], vertices[j]))
                matrix[i][j] = matrix[j][i] = value
        bd = cm.Atom("c2", 0, tuple(b((vertices[i], vertices[i+1])) for i in range(n)))
        sd = cm.Diag("zsign", tuple((s((vertices[i], vertices[i+1])), (0,)*n) for i in range(n)))
        wd = cm.Diag("c2", tuple((0, tuple(row)) for row in matrix))
        return bd, (sd, wd)

    def closed_pair_gauge(self, b, c, s, omega):
        """q3 with Phi_on−Phi_raw=delta_s q3+h(s t), modulo1."""
        bc, cc = pullback_interval(b), pullback_interval(c, True)
        sc, wc = pullback_interval(s), pullback_interval(omega)
        c_primitive = prism(self.adjusted_difference(bc, cc, sc, wc))
        source_primitive = Cochain(3, lambda vertices: self.chain_value(
            self.homotopy(self.section(b, s, omega, vertices))))
        return c_primitive+source_primitive

    def raw_phase(self, b, c, bp, cp, s, omega):
        bi, ci = pullback_interval(b, True), pullback_interval(c, True)
        bpi, cpi = pullback_interval(bp, True), pullback_interval(cp, True)
        si, wi = pullback_interval(s), pullback_interval(omega)
        db, dc = lower_d(bi, ci, si, wi)
        dbp, dcp = lower_d(bpi, cpi, si, wi)
        eta = binary_sum(beta(bi, bpi, si, wi), pullback_interval(beta(b, bp, s, omega), True))
        total_c = binary_sum(dc, dcp, alpha(db, dbp, si))
        transport = F(1, 2)*binary_sum(E(eta, wi), polarization(total_c, differential(eta).mod2()))
        return -prism(self.successor.phase(db, dc, dbp, dcp, si, wi)+transport)

    @staticmethod
    def _successor_rephasing(c, s, b_zero):
        return F(1, 2)*cup(s, c) if b_zero else zero(4)

    def successor_phase(self, b, c, bp, cp, s, omega, *, b_zero, bp_zero, sum_zero):
        """Legal k4 phase in the calibration compatible with this k3 rule.

        Each flag denotes zero of the entire degree-two B cochain, not
        its value on a simplex. Since the input is legal, selected C is closed.
        """
        total_c = binary_sum(c, cp, alpha(b, bp, s))
        r = self._successor_rephasing
        return (self.successor.phase(b, c, bp, cp, s, omega)+r(total_c, s, sum_zero)
                -r(c, s, b_zero)-r(cp, s, bp_zero))

    def phase(self, b, c, bp, cp, s, omega, *, b_closed, bp_closed, sum_closed):
        """Degree-three phase coordinate for arbitrary lower inputs."""
        total_b, total_c = binary_sum(b, bp), binary_sum(c, cp, beta(b, bp, s, omega))
        def q(bb, cc, flag):
            return self.closed_pair_gauge(bb, cc, s, omega) if flag else zero(3)
        return (self.raw_phase(b, c, bp, cp, s, omega)-q(b, c, b_closed)
                -q(bp, cp, bp_closed)+q(total_b, total_c, sum_closed))

    def gamma(self, b, c, bp, cp, s, omega, *, b_closed, bp_closed, sum_closed):
        total_b, total_c = binary_sum(b, bp), binary_sum(c, cp, beta(b, bp, s, omega))
        db, dc = lower_d(b, c, s, omega)
        dbp, dcp = lower_d(bp, cp, s, omega)
        phase = self.phase(b, c, bp, cp, s, omega, b_closed=b_closed,
                           bp_closed=bp_closed, sum_closed=sum_closed)
        successor = self.successor_phase(db, dc, dbp, dcp, s, omega,
                        b_zero=b_closed, bp_zero=bp_closed, sum_zero=sum_closed)
        return integral(self.potential(b, c, s, omega, b_closed=b_closed)
                        +self.potential(bp, cp, s, omega, b_closed=bp_closed)
                        -self.potential(total_b, total_c, s, omega, b_closed=sum_closed)
                        +signed_differential(phase, s)+successor, "all-cochain A=0 gamma3")

    def successor_gamma(self, b, c, bp, cp, s, omega, *, b_zero, bp_zero, sum_zero):
        total_b, total_c = binary_sum(b, bp), binary_sum(c, cp, alpha(b, bp, s))
        phase = self.successor_phase(b, c, bp, cp, s, omega,
                        b_zero=b_zero, bp_zero=bp_zero, sum_zero=sum_zero)
        return integral(self.successor.omega(b, c, s, omega)+self.successor.omega(bp, cp, s, omega)
                        -self.successor.omega(total_b, total_c, s, omega)
                        +signed_differential(phase, s), "compatible legal A=0 gamma4")

    def g(self, b, c, s, omega, *, b_closed):
        db, dc = lower_d(b, c, s, omega)
        return integral(signed_differential(self.potential(b, c, s, omega, b_closed=b_closed), s)
                        -self.successor.omega(db, dc, s, omega), "existing A=0 g3")
