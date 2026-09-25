"""Explicit prime exchange for the existing nonzero-even-rank k3 phase.

At least one rank must be nonzero on every component. Its sign supplies
an explicit primitive of the old phase's h(s B B') exchange class.
The both-zero-rank phase requires the separately calibrated paper rule.
"""
from fractions import Fraction as F
from cochains import Cochain,binary_sum,cup,differential,signed_differential
from compatible_sector import E,integral,polarization
import even_a_gamma3 as even

# Entry 8*g+h, g=(B01,B'01,s01), h=(B12,B'12,s12), bit order 0,1,2.
# This fixed table differentiates to W0-h(s B B') for every closed omega.
_SOURCE_EIGHTHS=(0,0,0,0,0,0,0,0,
 0,0,5,3,0,0,1,7,
 0,0,0,4,0,4,0,0,
 0,2,5,1,4,2,5,5,
 0,0,0,7,0,0,0,7,
 0,0,3,4,0,0,7,0,
 0,2,0,1,0,6,0,5,
 0,0,7,0,4,0,7,4)


def source_primitive(B,Bprime,s):
    def value(face):
        g=B(face[:2])+2*Bprime(face[:2])+4*s(face[:2])
        h=B(face[1:])+2*Bprime(face[1:])+4*s(face[1:])
        return F(_SOURCE_EIGHTHS[8*g+h],8)
    return Cochain(2,value,'fixed even-rank exchange source')


def rank_orientation(A,Aprime):
    def value(face):
        a=A(face)
        if not a:a=Aprime(face)
        if not a:raise ValueError('require a nonzero even rank on every component')
        if a%2:raise ValueError('require both ranks even')
        return int(a<0)
    return Cochain(0,value,'nonzero rank sign')


def lower_gauge(A,B,Aprime,Bprime,s):
    t=integral(F(1,2)*A,'A/2').mod2()
    tp=integral(F(1,2)*Aprime,'Aprime/2').mod2()
    return binary_sum(cup(B,Bprime,1),cup(cup(t,tp),s))


def exchange(A,B,C,Aprime,Bprime,Cprime,s,omega):
    """Return lambda,K,sigma,L,M with gamma_op-gamma+L=delta_s M."""
    d=even.canonical_data(A,B,C,s,omega)
    dp=even.canonical_data(Aprime,Bprime,Cprime,s,omega)
    m=even.alpha(B,Bprime,s)
    beta=even.beta_legal(A,B,Aprime,Bprime,s,omega)
    betap=even.beta_legal(Aprime,Bprime,A,B,s,omega)
    Bs=binary_sum(B,Bprime)
    Cs=binary_sum(C,Cprime,beta)
    Cop=binary_sum(C,Cprime,betap)
    total=even.canonical_data(A+Aprime,Bs,Cs,s,omega)
    totalop=even.canonical_data(A+Aprime,Bs,Cop,s,omega)
    q=cup(B,Bprime,1)
    tt=cup(d.t,dp.t)
    ell=cup(tt,s)
    lam=binary_sum(q,ell)
    dl=differential(lam).mod2()
    K=(F(1,2)*binary_sum(E(lam,omega),polarization(total.cbar,dl))
       +totalop.gauge-total.gauge)
    orientation=rank_orientation(A,Aprime)
    sigma=(F(1,2)*cup(d.cbar,dp.cbar,2)+source_primitive(B,Bprime,s)
           +F(1,2)*cup(cup(orientation,B),Bprime)
           +F(1,4)*cup(tt,omega,integral=True)
           +F(1,2)*binary_sum(cup(tt,cup(omega,s,1)),polarization(q,ell)))
    L=integral(even.omega3(A+Aprime,Bs,Cop,s,omega)
               -even.omega3(A+Aprime,Bs,Cs,s,omega)
               -signed_differential(K,s),'even-rank C transport')
    M=integral(even.phase3(Aprime,Bprime,Cprime,A,B,C,s,omega)
               -even.phase3(A,B,C,Aprime,Bprime,Cprime,s,omega)
               -K-signed_differential(sigma,s),'even-rank exchange carry')
    return lam,K,sigma,L,M
