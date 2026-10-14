"""Legal degree-three phase and integral stacking with nonzero even A.

Inputs obey delta_s A=0, A even, delta B=0 and delta C=Tau0(A;B).
Both twists are arbitrary closed binary cochains. Every expression is local
and explicit; no primitive on X is solved. Odd A and off-shell data are not
covered by this module.
"""
from dataclasses import dataclass
from fractions import Fraction as F

from cochains import Cochain,binary_sum,cup,differential,signed_differential,zero
from compatible_sector import E,Q,QD,alpha,integral,polarization
from a0_gamma import theta_degree_one,phase3_production
from n0_beta import beta_legal


@dataclass(frozen=True)
class CanonicalData:
    K: Cochain
    t: Cochain
    lam: Cochain
    H: Cochain
    z: Cochain
    cbar: Cochain
    gauge: Cochain


def background_carries(omega):
    R=integral(F(1,2)*differential(omega),'background integral Bockstein')
    r=R.mod2()
    carry=integral(F(1,2)*(R-r),'background Bockstein carry').mod2()
    return r,carry


def canonical_data(A,B,C,s,omega):
    if (A.degree,B.degree,C.degree,s.degree,omega.degree)!=(0,1,2,1,2):
        raise ValueError('expected degrees (0,1,2,1,2)')
    K=integral(F(1,2)*A,'even A/2')
    t=K.mod2()
    quarter=integral(F(1,2)*(K-t),'quarter A carry')
    Z=cup(t,omega)
    sb=cup(s,B)
    carry=integral(F(1,2)*(Z+sb-binary_sum(Z,sb)),'source sum carry')
    lam=binary_sum(cup(quarter,omega,integral=True),carry)
    H=QD(B,s,omega)
    r,_=background_carries(omega)
    z=cup(t,r)
    cbar=binary_sum(C,lam)
    gauge=F(1,2)*binary_sum(polarization(C,lam),Q(lam,1))
    return CanonicalData(K,t,lam,H,z,cbar,gauge)


def pontryagin(omega):
    return cup(omega,omega,integral=True)+cup(omega,differential(omega),1,integral=True)


def omega3(A,B,C,s,omega):
    """Exact even-branch production Omega3 representative."""
    d=canonical_data(A,B,C,s,omega)
    return (F(1,2)*E(C,omega)+theta_degree_one(B,s,omega)
            -F(1,8)*cup(d.K,pontryagin(omega),integral=True)
            +F(1,2)*binary_sum(cup(d.H,d.z,2),E(d.lam,omega),
                               cup(binary_sum(d.H,d.z),differential(d.lam).mod2(),2)))


def canonical_omega(d,B,s,omega):
    return (F(1,2)*E(d.cbar,omega)+theta_degree_one(B,s,omega)
            -F(1,8)*cup(d.K,pontryagin(omega),integral=True)
            +F(1,2)*cup(d.H,d.z,2))


def zero_a_source_phase(B,Bprime,s,omega):
    """Explicit A=0 source primitive independent of the C representatives."""
    H,Hp=QD(B,s,omega),QD(Bprime,s,omega)
    m=alpha(B,Bprime,s)
    return (phase3_production(B,zero(2),Bprime,zero(2),s,omega)
            +F(1,2)*binary_sum(cup(binary_sum(H,Hp),m,2),cup(H,Hp,3)))


def phase3(A,B,C,Aprime,Bprime,Cprime,s,omega):
    """Explicit e3 with delta_s e3=Delta Omega3 modulo one on legal data."""
    d=canonical_data(A,B,C,s,omega)
    dp=canonical_data(Aprime,Bprime,Cprime,s,omega)
    beta=beta_legal(A,B,Aprime,Bprime,s,omega)
    total=canonical_data(A+Aprime,binary_sum(B,Bprime),binary_sum(C,Cprime,beta),s,omega)
    m=alpha(B,Bprime,s)
    direct=F(1,2)*binary_sum(polarization(d.cbar,dp.cbar),
                             polarization(binary_sum(d.cbar,dp.cbar),m))
    _,carry=background_carries(omega)
    extra=F(1,2)*binary_sum(cup(d.z,dp.H,3),
                            cup(binary_sum(d.z,dp.z),differential(m).mod2(),3),
                            cup(cup(d.t,dp.t),carry))
    return (direct+zero_a_source_phase(B,Bprime,s,omega)+extra
            +total.gauge-d.gauge-dp.gauge)


def gamma3(A,B,C,Aprime,Bprime,Cprime,s,omega):
    """Integral degree-four correction for the stated even-A legal domain."""
    beta=beta_legal(A,B,Aprime,Bprime,s,omega)
    totalB=binary_sum(B,Bprime)
    totalC=binary_sum(C,Cprime,beta)
    return integral(omega3(A,B,C,s,omega)+omega3(Aprime,Bprime,Cprime,s,omega)
                    -omega3(A+Aprime,totalB,totalC,s,omega)
                    +signed_differential(phase3(A,B,C,Aprime,Bprime,Cprime,s,omega),s),
                    'even-A legal gamma3')
