"""Fixed legal pair phases for the prescribed Theta boundary helper.

These are auxiliary source formulas, not a replacement calibration of the
production upper differential. They cover q-degrees 3,4,5 with generic
twists, using the reduced binary pair contraction.
"""
from cochains import Cochain, zero
from a0_high_gamma import HigherA0Stacking, p


class ThetaPairPhase(HigherA0Stacking):
    def __init__(self, m):
        if m not in (3, 4, 5):
            raise ValueError('Theta pair source is implemented in degrees 3,4,5')
        self.m, self.k = m, m+2
        if m == 3 and self.minimal_period() % 1:
            raise ArithmeticError('the minimal Theta pair period is nonintegral')

    def boundary_phase(self, b, s, omega):
        return Cochain(self.m+3, p.theta(b, zero(self.m+1), s, omega))

    def omega(self, b, c, s, omega):
        # Keep Theta's actual binary-sum lift, including its integral carry
        # relative to the separated h(Ec)+Theta(b,0) representation.
        return Cochain(self.m+3, p.theta(b, c, s, omega))
