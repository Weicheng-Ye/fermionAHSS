"""Legal A=0 degree-four stacking in the repository's exact R1 gauge.

The phase primitive uses a fixed four-factor universal chain contractor,
not a cochain solver on X. The production modules are loaded explicitly
from the supplied repository so their calibrated helpers remain intact.
This module does not claim an all-cochain degree-four stacking law.
"""
from fractions import Fraction as F
from functools import lru_cache
import importlib
import importlib.util
from pathlib import Path
import sys

from cochains import Cochain, binary_sum, cup, differential, signed_differential
from compatible_sector import E, Q, QD, alpha, integral, polarization
from lower_stacking import pullback_interval, prism


TREE = (("b", "bp"), ("s", "omega"))
SOURCE_BASES = ((((1,), (2,)), (0, ())),
                (((2,), (1,)), (0, ())),
                (((1,), (1,)), (1, ())))


class DegreeFourA0Stacking:
    """Fixed production-gauge formula for closed B,B' and legal C,C'."""

    def __init__(self, production_root):
        root = Path(production_root).resolve()
        directory = root / "python"
        required = ("__init__.py", "phase_eval.py", "chain_models.py")
        if any(not (directory / name).is_file() for name in required):
            raise ValueError("production_root must contain the fermionAHSS Python package")
        # A dedicated package avoids collisions with unrelated modules named
        # python/phase_eval. Reusing the same path also reuses its fixed caches.
        name = "_stacking_fermionahss_" + str(abs(hash(str(root))))
        if name not in sys.modules:
            spec = importlib.util.spec_from_file_location(
                name, directory / "__init__.py",
                submodule_search_locations=[str(directory)])
            module = importlib.util.module_from_spec(spec)
            sys.modules[name] = module
            spec.loader.exec_module(module)
        self.production = importlib.import_module(name + ".phase_eval")
        self.chains = importlib.import_module(name + ".chain_models")
        self.package_name = name
        # These three generators exhaust the normalized source in degree5.
        # An integral value is zero in Q/Z, so the small primitive is zero.
        self.calibration = tuple(self.chain_value(self.inclusion(basis))
                                 for basis in SOURCE_BASES)
        if any(value % 1 for value in self.calibration):
            raise ArithmeticError("the fixed degree-four source calibration is no longer zero modulo1")

    def boundary_phase(self, b, s, omega):
        """Exact A=0 R1 production phase, including the suspension terms."""
        if (b.degree, s.degree, omega.degree) != (2, 1, 2):
            raise ValueError("expected B,s,omega degrees (2,1,2)")
        bi = pullback_interval(b, True)
        si, wi = pullback_interval(s), pullback_interval(omega)
        suspended = self.production.theta(differential(bi).mod2(), QD(bi, si, wi), si, wi)
        return -prism(suspended) + F(1, 2)*binary_sum(E(Q(b, 1), omega),
                                                    cup(s, QD(b, s, omega)))

    def omega(self, b, c, s, omega):
        return F(1, 2)*E(c, omega) + self.boundary_phase(b, s, omega)

    def source(self, b, bp, s, omega):
        """Closed Q/Z degree5 source depending only on B,B',s,omega."""
        tau, taup = QD(b, s, omega), QD(bp, s, omega)
        beta = alpha(b, bp, s)
        dbeta = differential(beta).mod2()
        change = (self.boundary_phase(binary_sum(b, bp), s, omega)
                  - self.boundary_phase(b, s, omega)
                  - self.boundary_phase(bp, s, omega))
        correction = binary_sum(polarization(tau, taup), E(beta, omega),
                                  polarization(binary_sum(tau, taup), dbeta))
        return change + F(1, 2)*correction

    def _tensor_maps(self, pair, left, right):
        out = {}
        for a, c in left(pair[0]).items():
            for b, d in right(pair[1]).items():
                self.chains.add(out, {(a, b): c*d})
        return out

    @lru_cache(None)
    def forward(self, simplex, tree=TREE):
        cm = self.chains
        if isinstance(tree, str):
            if tree != "s":
                return cm.F(simplex)
            return {} if cm._degen_indices(simplex) else {cm._degree(simplex): 1}
        return cm.linear(cm.product_AW(simplex), lambda pair: self._tensor_maps(
            pair, lambda x: self.forward(x, tree[0]), lambda x: self.forward(x, tree[1])))

    @lru_cache(None)
    def inclusion(self, basis, tree=TREE):
        cm = self.chains
        if isinstance(tree, str):
            if tree != "s":
                return cm.G("c2", basis)
            # Pure-sign zsign diagonal. Its pullback, including initial-frame
            # transport, is supplied by cm._pull; no untwisted Atom is used.
            return {cm.Diag("zsign", tuple((1, (0,)*basis) for _ in range(basis))): 1}
        tensor = self._tensor_maps(basis, lambda x: self.inclusion(x, tree[0]),
                                   lambda x: self.inclusion(x, tree[1]))
        return cm.linear(tensor, cm.product_shuffle)

    @lru_cache(None)
    def homotopy(self, simplex, tree=TREE):
        cm = self.chains
        if isinstance(tree, str):
            return {} if tree == "s" else cm.H(simplex)
        out = dict(cm.product_homotopy(simplex))
        for (a, b), coefficient in cm.product_AW(simplex).items():
            tensor = {(h, b): c for h, c in self.homotopy(a, tree[0]).items()}
            projected = cm.linear(self.forward(a, tree[0]),
                                  lambda key: self.inclusion(key, tree[0]))
            for aa, c in projected.items():
                for bb, d in self.homotopy(b, tree[1]).items():
                    cm.add(tensor, {(aa, bb): (-1)**cm._degree(a)*c*d})
            cm.add(out, cm.linear(tensor, cm.product_shuffle), coefficient)
        return out

    def _inputs(self, simplex):
        (bb, bpb), (sb, ww) = simplex
        s = Cochain(1, lambda face: sum(row[0] for row in sb.rows[face[0]:face[1]]) % 2)
        def rectangular_value(diag):
            return Cochain(2, lambda face: sum(
                diag.rows[i][1][j] for i in range(face[0], face[1])
                for j in range(face[1], face[2])) % 2)
        return rectangular_value(bb), rectangular_value(bpb), s, rectangular_value(ww)

    @lru_cache(None)
    def evaluate(self, simplex):
        return self.source(*self._inputs(simplex))(tuple(range(self.chains._degree(simplex)+1)))

    def chain_value(self, chain):
        return sum(coefficient*self.evaluate(simplex) for simplex, coefficient in chain.items())

    def section(self, b, bp, s, omega, vertices):
        """Natural symmetric binary section into the diagonal bar models."""
        cm, n = self.chains, len(vertices)-1
        def diagonal(q):
            matrix = [[0]*n for _ in range(n)]
            for i in range(n):
                for j in range(i+1, n):
                    value = q((vertices[i], vertices[i+1], vertices[j+1]))
                    if j > i+1:
                        value ^= q((vertices[i], vertices[i+1], vertices[j]))
                    matrix[i][j] = matrix[j][i] = value
            return cm.Diag("c2", tuple((0, tuple(row)) for row in matrix))
        sd = cm.Diag("zsign", tuple((s((vertices[i], vertices[i+1])), (0,)*n)
                                    for i in range(n)))
        return ((diagonal(b), diagonal(bp)), (sd, diagonal(omega)))

    def primitive(self, b, bp, s, omega):
        """Universal degree4 V with delta_s V=source modulo1."""
        return Cochain(4, lambda vertices: self.chain_value(
            self.homotopy(self.section(b, bp, s, omega, vertices))))

    def phase(self, b, c, bp, cp, s, omega):
        """Legal degree-four phase; requires closed B and delta C=QD B."""
        if (b.degree, c.degree, bp.degree, cp.degree) != (2, 3, 2, 3):
            raise ValueError("expected k4 input degrees B=2,C=3")
        beta = alpha(b, bp, s)
        direct = F(1, 2)*binary_sum(polarization(c, cp),
                                    polarization(binary_sum(c, cp), beta))
        return direct + self.primitive(b, bp, s, omega)

    def gamma(self, b, c, bp, cp, s, omega):
        """Integral degree5 correction on full legal lower data."""
        total_b = binary_sum(b, bp)
        total_c = binary_sum(c, cp, alpha(b, bp, s))
        return integral(self.omega(b, c, s, omega) + self.omega(bp, cp, s, omega)
                        - self.omega(total_b, total_c, s, omega)
                        + signed_differential(self.phase(b, c, bp, cp, s, omega), s),
                        "legal A=0 gamma4")
