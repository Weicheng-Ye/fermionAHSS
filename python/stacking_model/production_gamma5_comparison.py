"""Explicit n=2 production comparison, including the calibrated cubic phase."""
from fractions import Fraction as F
import production_gamma6_comparison as odd
p,upper=odd.p,odd.upper


def UA(P,k,s,omega):
    return p.binary(p.zeta1(s,P)+p.cup(s,k)+p.cup(p.cup(omega,s,1),P))


def reduced_A_phase(A,s,omega):
    if A.degree!=2:
        raise ValueError('this comparison covers input degree two')
    data=upper.hp.uniform_source(A,s,omega);P,k=data['p'],data['k0']
    return (odd.reduced_A_phase(A,s,omega)+p.half(p.binary(p.cup(s,k)+UA(P,k,s,omega)))
            -p.scale(p.QD(P,s,omega),F(1,4)))


def A_phase(A,s,omega):
    import low_phases,low_calibration
    alpha=low_calibration.constants()['alpha']
    linear=low_phases.transported(low_phases.pontryagin(omega),A,s)
    cube=low_phases.transported(low_phases.transported(A,A,s),A,s)
    return reduced_A_phase(A,s,omega)+p.scale(linear,F(alpha,4))-p.scale(cube,F(1,4))


def B_phase(B,s,omega):
    return p.half(p.cup(s,p.QD(B,s,omega)))


def comparison_gauge(A,B,s,omega):
    lam=odd.carry(A,B,s,omega);y=p.QD(B,s,omega)
    I=upper.hp.interval;Bi,wi=I(B,True),I(omega)
    return (odd.comparison_gauge(A,B,s,omega)+p.half(p.cup(s,lam))
            -upper.hp.prism(p.half(p.E(p.Q(Bi,1),wi)))+p.scale(y,F(1,4)))


def reduced_production(A,B,C,s,omega):
    data=upper.hp.uniform_source(A,s,omega);P,k=data['p'],data['k0'];y=p.QD(B,s,omega)
    lam=odd.carry(A,B,s,omega);dl=p.binary(p.differential(lam));m=4
    core=p.half(p.binary(p.cup(y,k,m)+p.cup(k,p.binary(p.differential(y)),m+1)
                  +p.Q(k,1)+p.E(lam,omega)+p.cup(p.binary(y+k),dl,m)))
    I=upper.hp.interval;Bi,si,wi=I(B,True),I(s),I(omega)
    Pi=upper.hp.prism(p.theta(p.binary(p.differential(Bi)),p.QD(Bi,si,wi),si,wi))
    return p.half(p.E(C,omega))+core+Pi-p.half(p.E(p.cup(s,B),omega))+p.half(UA(P,k,s,omega))


def affine_B_source(A,Aprime,s,omega):
    alpha=upper.hp.hD(p.binary(A),p.binary(Aprime),s)
    P,Pp=upper.primary(A,s,omega),upper.primary(Aprime,s,omega)
    q=p.binary(p.differential(alpha))
    return p.half(p.cup(s,p.binary(p.QD(alpha,s,omega)+upper.hp.hD(P,Pp,s)
                        +upper.hp.hD(p.binary(P+Pp),q,s))))


def affine_B_primitive(A,B,Aprime,Bprime,s,omega):
    alpha=upper.hp.hD(p.binary(A),p.binary(Aprime),s)
    return p.half(p.cup(s,p.binary(upper.hp.hD(B,Bprime,s)
                        +upper.hp.hD(p.binary(B+Bprime),alpha,s))))
