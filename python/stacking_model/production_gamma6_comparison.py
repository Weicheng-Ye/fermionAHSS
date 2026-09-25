"""Explicit legal-data comparison of production Omega6 and auxiliary Phi.

The source/coboundary identity is separate from solving the resulting
A-only stacking source. Production phases and calibration data are read only.
"""
from fractions import Fraction as F
from functools import lru_cache
import closed_a_upper as upper
import production_upper_binary_comparison as binary
from cochains import Cochain
from theta_pair_phase import ThetaPairPhase
p=upper.p


@lru_cache(None)
def theta_pair(m=5):
    return ThetaPairPhase(m)


def adapt(c):
    return Cochain(c.degree,c)


def phase_theta(P,k,Pprime,kprime,s,omega):
    op=theta_pair(P.degree).phase(*(adapt(x) for x in (P,k,Pprime,kprime,s,omega)))
    return p.Cochain(op.degree,op)


def carry(A,B,s,omega):
    return upper.hp.splitting_carry(A,B,s,omega)


def C_transport(c,L,omega):
    return p.half(p.binary(p.E(L,omega)+upper.hE(c,p.binary(p.differential(L)))))


def A_phase(A,s,omega,include_three=True):
    """R_A in Omega_production-Phi=R_A+h J0(B)+delta_s Q."""
    if A.degree!=3:
        raise ValueError('this comparison covers input degree three')
    import high_phase,high_calibration
    data=upper.hp.uniform_source(A,s,omega);P,k=data['p'],data['k0']
    r=p.binary(p.Q(P,1)+p.cup(s,P))
    value=p.half(p.binary(p.cup(s,k)+p.Q(r,1)+p.binary(p.differential(k))
                        +p.E(P,omega)+upper.hE(p.Q(P,1),p.cup(s,P))))
    theta=theta_pair(P.degree).primitive(*(adapt(x) for x in(P,P,s,omega)))
    value=value+p.Cochain(theta.degree,theta)
    constants=high_calibration.constants();terms=high_phase.corrections(A,s,omega)
    for name,coefficient in [('L',constants['c4']),
            ('N',constants['cN']+constants['epsilon_c']),
            ('O',constants['cO']+constants['epsilon_c']),
            ('M',constants['cM'])]:
        value=value+p.scale(terms[name],coefficient)
    if include_three:
        value=value+high_phase.tertiary_three_primary_phase(A,s)
    return value


def comparison_gauge(A,B,s,omega):
    data=upper.hp.uniform_source(A,s,omega);P,k=data['p'],data['k0'];m=A.degree+2
    lam=carry(A,B,s,omega)
    Q=binary.comparison_terms(B,k,lam,s,omega)[3]
    I=upper.hp.interval
    Ai,Bi,si,wi=I(A),I(B,True),I(s),I(omega)
    Pi,ki=I(P),I(k);u=p.binary(p.differential(Bi));y=p.QD(Bi,si,wi)
    c0=p.binary(ki+y+upper.hp.hD(Pi,u,si))
    L=p.binary(carry(Ai,Bi,si,wi)+p.Q(Bi,1)+p.cup(si,Bi)+p.cup(Pi,u,m))
    trans=phase_theta(Pi,ki,u,y,si,wi)+C_transport(c0,L,wi)
    return p.half(Q)-upper.hp.prism(trans)


def reduced_production(A,B,C,s,omega):
    """Omega6 with the common -V3 removed, for a cheap comparison check."""
    data=upper.hp.uniform_source(A,s,omega);k=data['k0'];y=p.QD(B,s,omega);m=A.degree+2
    lam=carry(A,B,s,omega);dl=p.binary(p.differential(lam));tau=p.binary(y+k+dl)
    core=p.half(p.binary(p.cup(y,k,m)+p.cup(k,p.binary(p.differential(y)),m+1)
               +p.Q(k,1)+p.E(lam,omega)+p.cup(p.binary(y+k),dl,m)))
    I=upper.hp.interval;Bi=I(B,True);si,wi=I(s),I(omega)
    Pi=upper.hp.prism(p.theta(p.binary(p.differential(Bi)),p.QD(Bi,si,wi),si,wi))
    value=p.half(p.E(C,omega))+core-Pi+p.half(p.binary(p.E(p.Q(B,1),omega)+p.cup(s,tau)))
    return value


def reduced_auxiliary(A,B,C,s,omega):
    I=upper.hp.interval;Ai,Bi,si,wi=I(A),I(B,True),I(s),I(omega)
    bi=p.binary(p.differential(Bi)+upper.primary(Ai,si,wi))
    fi=upper.fsharp(Ai,Bi,si,wi);f=upper.fsharp(A,B,s,omega)
    return p.half(p.binary(p.E(C,omega)+upper.hE(p.binary(p.differential(C)),f)))-upper.hp.prism(p.theta(bi,fi,si,wi))


def reduced_A_phase(A,s,omega):
    data=upper.hp.uniform_source(A,s,omega);P,k=data['p'],data['k0']
    r=p.binary(p.Q(P,1)+p.cup(s,P))
    theta=theta_pair(P.degree).primitive(*(adapt(x) for x in(P,P,s,omega)))
    return p.Cochain(theta.degree,theta)+p.half(p.binary(p.cup(s,k)+p.Q(r,1)
        +p.binary(p.differential(k))+p.E(P,omega)+upper.hE(p.Q(P,1),p.cup(s,P))))
