"""Legal A=0 k3 phase transported from the proved repaired paper exchange.

This is the selected both-zero-rank full-legal pair-phase calibration.
The global all-cochain branch is in degree3_legal_recalibration.py.
"""
from fractions import Fraction as F
from cochains import binary_sum,cup,differential,signed_differential
from compatible_sector import integral
import a0_gamma as production
import paper_stacking as paper


def reorder_gauge(b,s,omega):
    """u_B with r_B=(B+s)Tau+Rprime+delta u_B, modulo two."""
    return binary_sum(cup(cup(b,s,1),cup(b,b)),
                      cup(cup(b,omega,1),b))


def coordinate_gauge(b,c,s,omega):
    return (production._eighth_cube(b)-production._paper_c_gauge(c,s)
            +F(1,2)*binary_sum(cup(binary_sum(b,s),c),reorder_gauge(b,s,omega)))


def additive_remainder(b,s,omega):
    rw=integral(F(1,2)*differential(omega),'omega Bockstein').mod2()
    return F(1,2)*binary_sum(cup(omega,cup(b,b)),cup(cup(s,omega),b),cup(rw,b))


def phase(b,c,bp,cp,s,omega):
    corrected=paper.repaired_paper_corrections(3,b,c,bp,cp,s,omega)
    bs,cs=binary_sum(b,bp),binary_sum(c,cp,corrected.beta)
    lam=cup(b,bp,1)
    return (corrected.phase+coordinate_gauge(bs,cs,s,omega)
            -coordinate_gauge(b,c,s,omega)-coordinate_gauge(bp,cp,s,omega)
            +F(1,2)*cup(omega,lam))


def gamma(b,c,bp,cp,s,omega):
    beta=paper.repaired_paper_corrections(3,b,c,bp,cp,s,omega).beta
    bs,cs=binary_sum(b,bp),binary_sum(c,cp,beta)
    return integral(production.omega3(b,c,s,omega)+production.omega3(bp,cp,s,omega)
                    -production.omega3(bs,cs,s,omega)+signed_differential(phase(b,c,bp,cp,s,omega),s),
                    'commutative paper-calibrated production gamma3')


def transport(b,c,lam,s,omega):
    shifted=binary_sum(c,differential(lam).mod2())
    return (paper.phase_gauge_transport_k3(c,lam,s,omega)
            +coordinate_gauge(b,shifted,s,omega)-coordinate_gauge(b,c,s,omega))


def exchange(b,c,bp,cp,s,omega):
    """Explicit lambda,K,sigma and integral carries L,M for prime exchange."""
    beta=paper.repaired_paper_corrections(3,b,c,bp,cp,s,omega).beta
    bs,cs=binary_sum(b,bp),binary_sum(c,cp,beta)
    lam=cup(b,bp,1)
    sigma=paper.phase_swap_gauge_k3(b,c,bp,cp,s)
    K=transport(bs,cs,lam,s,omega)
    opposite=binary_sum(cs,differential(lam).mod2())
    L=integral(production.omega3(bs,opposite,s,omega)-production.omega3(bs,cs,s,omega)
               -signed_differential(K,s),'production C-transport carry')
    M=integral(phase(bp,cp,b,c,s,omega)-phase(b,c,bp,cp,s,omega)
               -K-signed_differential(sigma,s),'production phase-exchange carry')
    return lam,K,sigma,L,M
