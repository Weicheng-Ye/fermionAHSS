"""Degree-one product coherently calibrated to the selected degree-two rule.

Only C0,D2 are present at this degree. The phase is regenerated from the
actual successor instead of assuming the independent pure-C calibration.
The absent C layer in degree zero has ordinary integral addition.
"""
from fractions import Fraction as F

from cochains import Cochain, binary_sum, differential, signed_differential, zero
from compatible_sector import E, integral, pure_c_g
from lower_stacking import pullback_interval, prism
from a0_degree2 import DegreeTwoA0Stacking


def edge_phase(c, cp, s):
    """Exact nine-monomial rational lift of the selected degree-one phase.

    The variables are binary endpoint values, and products are ordinary
    rational products, not mod-two cochain products. Its equality to the
    successor-prism expression has a complete 32-configuration certificate.
    """
    if c.degree != 0 or cp.degree != 0 or s.degree != 1:
        raise ValueError('expected degree-zero C,Cprime and degree-one s')
    def value(face):
        x, y = c((face[0],)), c((face[1],))
        u, v = cp((face[0],)), cp((face[1],))
        z = s(face)
        return (F(3, 2)*x*u+F(1, 2)*y*u-x*y*u+x*v-x*y*v
                -2*x*u*v-y*u*v+2*x*y*u*v+z*y*v)
    return Cochain(1, value)


class DegreeOneA0Stacking:
    def __init__(self, production_root=None, *, successor=None):
        if successor is None:
            successor = DegreeTwoA0Stacking(production_root)
        self.successor = successor

    def successor_phase(self, c, cp, s, omega):
        return self.successor.phase(zero(0), c, zero(0), cp, s, omega)

    def phase_from_successor(self, c, cp, s, omega):
        """Finite prism definition retained to certify the chosen calibration."""
        if c.degree != 0 or cp.degree != 0:
            raise ValueError('degree-one C cochains have degree zero')
        ci, cpi = pullback_interval(c, True), pullback_interval(cp, True)
        si, wi = pullback_interval(s), pullback_interval(omega)
        return -prism(self.successor_phase(differential(ci).mod2(),
                                           differential(cpi).mod2(), si, wi))

    @staticmethod
    def phase(c, cp, s, omega):
        return edge_phase(c, cp, s)

    def gamma(self, c, cp, s, omega):
        total = binary_sum(c, cp)
        return integral(F(1, 2)*(E(c, omega)+E(cp, omega)-E(total, omega))
                        +signed_differential(self.phase(c, cp, s, omega), s)
                        +self.successor_phase(differential(c).mod2(),
                                              differential(cp).mod2(), s, omega),
                        'coherent all-cochain gamma1')

    def successor_gamma(self, c, cp, s, omega):
        return self.successor.gamma(zero(0), c, zero(0), cp, s, omega,
                                    b_closed=True, bp_closed=True, sum_closed=True)

    @staticmethod
    def g(c, s, omega):
        return pure_c_g(c, s, omega)
