"""Degrees0--2 with the commutative legal A=0 degree-three successor.

The natural legal successor phase is the explicit repaired-paper transport.
The existing prism constructions regenerate all-cochain gamma2 and gamma1.
The earlier production calibration remains in coherent_low_production.py.
"""
from fractions import Fraction
from cochains import Cochain
from a0_degree2 import DegreeTwoA0Stacking
from a0_degree1 import DegreeOneA0Stacking
import paper_commutative_production as legal


class CommutativeLegalDegreeThree:
    @staticmethod
    def _flags(b_closed,bp_closed,sum_closed):
        if not(b_closed and bp_closed and sum_closed):
            raise ValueError('the selected degree-three successor requires legal data')

    def phase(self,b,c,bp,cp,s,omega,*,b_closed,bp_closed,sum_closed):
        self._flags(b_closed,bp_closed,sum_closed)
        return legal.phase(b,c,bp,cp,s,omega)

    def gamma(self,b,c,bp,cp,s,omega,*,b_closed,bp_closed,sum_closed):
        self._flags(b_closed,bp_closed,sum_closed)
        return legal.gamma(b,c,bp,cp,s,omega)


class DegreeTwoCommutativeStacking(DegreeTwoA0Stacking):
    def __init__(self):
        super().__init__(successor=CommutativeLegalDegreeThree())


class DegreeOneCommutativeStacking(DegreeOneA0Stacking):
    def __init__(self,*,successor=None):
        super().__init__(successor=successor or DegreeTwoCommutativeStacking())

    def phase(self,c,cp,s,omega):
        if c.degree!=0 or cp.degree!=0:
            raise ValueError('degree-one C cochains have degree zero')
        return Cochain(1,lambda face:Fraction(c(face[:1])
                        *((cp(face[:1])+cp(face[1:]))%2),2))
