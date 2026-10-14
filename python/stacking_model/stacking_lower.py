"""Explicit lower stacking API with the calibrated universal source primitives.

alpha covers all input degrees n=0..3. legal_beta covers k=3..6 and
the auxiliary k7 successor. all_cochain_beta covers k=3..6 using these
explicit successive formulas, with no restriction on the original A,B.
These are the first three coordinates, not a full integral stacking product.
The degree-five legal calibration includes rho A cup rho A'; this same
function is used in the degree-four off-shell successor prism.
"""
import cochains
import n0_beta
import off_shell_beta as off
from off_shell_beta import LowerPair, p


def alpha(A, Aprime, s):
    return off.alpha(A, Aprime, s)


def legal_beta(A, B, Aprime, Bprime, s, omega):
    """Explicit beta on signed-closed A and delta B=P(A), generic twists."""
    n = A.degree
    if (Aprime.degree, B.degree, Bprime.degree, s.degree, omega.degree) != (n,n+1,n+1,1,2):
        raise ValueError('inconsistent cochain degrees')
    if n == 0:
        def convert(c):
            return cochains.Cochain(c.degree, c)
        answer = n0_beta.beta_legal(*(convert(c) for c in (A,B,Aprime,Bprime,s,omega)))
        return p.Cochain(answer.degree, answer)
    if n == 1:
        from v1_pair_shared import legal_beta as formula
    elif n == 2:
        from v2_pair_shared import legal_beta as formula
    elif n == 3:
        from v3_pair_shared import legal_beta as formula
    elif n == 4:
        from v4_pair_shared import legal_beta as formula
    else:
        raise NotImplementedError('the legal beta source is implemented in degrees n=0..4')
    return formula(A,B,Aprime,Bprime,s,omega)


def all_cochain_beta(left, right, sum_legal, s, omega):
    """Strict lower compatibility for arbitrary inputs in k=3,4,5,6.

    left and right are LowerPair(A,B,global_legal). sum_legal is global
    legality of (A+A',B+B'+alpha). No user-supplied unknown operation occurs.
    """
    if left.A.degree not in (0,1,2,3):
        raise NotImplementedError('all-cochain beta is implemented through k6')
    return off.beta(left,right,sum_legal,s,omega,legal_beta,legal_beta)
