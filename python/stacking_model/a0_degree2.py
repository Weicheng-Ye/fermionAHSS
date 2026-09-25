"""Degree-two integral product coherently using a0_degree3 as its successor.

All A=0 degree-two cochains are allowed. The previous standalone a0_gamma
construction uses another valid successor calibration; this module selects
the phase actually used by DegreeThreeA0Stacking, so the two degrees compose.
"""
from cochains import binary_sum, differential, signed_differential
from compatible_sector import integral
from lower_stacking import pullback_interval, prism
from a0_gamma import (beta2, beta3_closed, lower_d2, phi2, g2,
                      phase3_C_transport)
from a0_degree3 import DegreeThreeA0Stacking


class DegreeTwoA0Stacking:
    def __init__(self, production_root=None, *, successor=None):
        if successor is None:
            if production_root is None:
                raise ValueError('provide production_root or an existing DegreeThreeA0Stacking')
            successor = DegreeThreeA0Stacking(production_root)
        self.successor = successor

    def successor_phase(self, b, c, bp, cp, s, omega):
        return self.successor.phase(b, c, bp, cp, s, omega,
                                    b_closed=True, bp_closed=True, sum_closed=True)

    def phase(self, b, c, bp, cp, s, omega):
        bi, ci = pullback_interval(b, True), pullback_interval(c, True)
        bpi, cpi = pullback_interval(bp, True), pullback_interval(cp, True)
        si, wi = pullback_interval(s), pullback_interval(omega)
        db, dc = lower_d2(bi, ci, si, wi)
        dbp, dcp = lower_d2(bpi, cpi, si, wi)
        eta = binary_sum(beta2(bi, bpi, si), pullback_interval(beta2(b, bp, s), True))
        total_c = binary_sum(dc, dcp, beta3_closed(db, dbp, si))
        return -prism(self.successor_phase(db, dc, dbp, dcp, si, wi)
                      +phase3_C_transport(total_c, eta, wi))

    def gamma(self, b, c, bp, cp, s, omega, *, b_closed, bp_closed, sum_closed):
        total_b = binary_sum(b, bp)
        total_c = binary_sum(c, cp, beta2(b, bp, s))
        db, dc = lower_d2(b, c, s, omega)
        dbp, dcp = lower_d2(bp, cp, s, omega)
        change = (phi2(b, c, s, omega, b_closed=b_closed)
                  +phi2(bp, cp, s, omega, b_closed=bp_closed)
                  -phi2(total_b, total_c, s, omega, b_closed=sum_closed))
        return integral(change+signed_differential(self.phase(b, c, bp, cp, s, omega), s)
                        +self.successor_phase(db, dc, dbp, dcp, s, omega),
                        'coherent all-cochain A=0 gamma2')

    def successor_gamma(self, b, c, bp, cp, s, omega):
        return self.successor.gamma(b, c, bp, cp, s, omega,
                                    b_closed=True, bp_closed=True, sum_closed=True)

    @staticmethod
    def g(b, c, s, omega, *, b_closed):
        return g2(b, c, s, omega, b_closed=b_closed)
