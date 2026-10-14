"""Production n=1 phase comparison through its explicit primary rephasing."""
from fractions import Fraction as F
import production_gamma6_comparison as odd
import production_upper_binary_comparison as binary
p,upper=odd.p,odd.upper


def source(A,s,omega):
    data=p.source(A,s,omega)
    return data['p'],data['q'],data['k'],p.cup(s,data['t'])


def shifted_lambda(A,B,s,omega):
    P,q,k,e=source(A,s,omega);de=p.binary(p.differential(e))
    N=odd.carry(A,B,s,omega)
    return p.binary(N+p.cup(B,e,1)+p.cup(P,e,2)+p.cup(P,de,3)
                    +p.cup(s,p.binary(p.cup(B,e,2)+p.cup(P,e,3)+p.cup(P,de,4))))


def fshift(A,B,s,omega):
    P,q,k,e=source(A,s,omega);bp=p.binary(B+e);u=p.binary(p.differential(bp))
    return p.binary(p.QD(bp,s,omega)+k+p.differential(shifted_lambda(A,B,s,omega))
                    +upper.hp.hD(u,q,s)+p.differential(p.binary(p.Q(bp,1)+p.cup(s,bp))))


def curvature_gauge(A,B,s,omega):
    P,q,k,e=source(A,s,omega);b=p.binary(p.differential(B)+P)
    return p.binary(p.cup(b,e,2)+p.cup(b,p.binary(p.differential(e)),3))


def source_prism(A,s,omega):
    P,q,k,e=source(A,s,omega);uniform=upper.hp.uniform_source(A,s,omega)
    I=upper.hp.interval;ei=I(e,True);de=p.binary(p.differential(ei));si,wi=I(s),I(omega)
    Pi,ki=I(P),I(uniform['k0'])
    kt=p.binary(ki+p.QD(ei,si,wi)+upper.hp.hD(de,Pi,si))
    return upper.hp.prism(p.theta(p.binary(Pi+de),kt,si,wi))


def reduced_phi_shift(A,B,C,s,omega):
    P,q,k,e=source(A,s,omega);bp=p.binary(B+e)
    I=upper.hp.interval;Ai,si,wi=I(A),I(s),I(omega)
    bpt=I(bp,True);Bt=p.binary(bpt+I(e));u=p.binary(p.differential(bpt))
    f=fshift(A,B,s,omega);ft=fshift(Ai,Bt,si,wi)
    initial=odd.C_transport(k,shifted_lambda(A,e,s,omega),omega)
    cross=p.half(p.binary(p.E(C,omega)+upper.hE(p.binary(p.differential(C)),f)))
    return cross-initial-upper.hp.prism(p.theta(p.binary(u+I(q)),ft,si,wi))


def closed_difference(A,B,C,s,omega):
    """Phi_shift+Ctransport-Phi_uniform, modulo one closed for all B,C."""
    b,c=upper.curvature(A,B,C,s,omega);L=curvature_gauge(A,B,s,omega)
    return (reduced_phi_shift(A,B,C,s,omega)+odd.C_transport(c,L,omega)
            -odd.reduced_auxiliary(A,B,C,s,omega)-source_prism(A,s,omega))


def shifted_A_phase(A,s,omega):
    P,q,k,e=source(A,s,omega);r=p.binary(p.Q(q,1)+p.cup(s,q))
    theta=odd.theta_pair(3).primitive(*(odd.adapt(x) for x in(q,q,s,omega)))
    return p.Cochain(theta.degree,theta)+p.half(p.binary(p.cup(s,k)+p.Q(r,1)
        +p.binary(p.differential(k))+p.E(q,omega)+upper.hE(p.Q(q,1),p.cup(s,q))))


def reduced_A_phase(A,s,omega):
    P,q,k,e=source(A,s,omega);de=p.binary(p.differential(e))
    return (closed_difference(A,p.zero(2),p.zero(3),s,omega)+shifted_A_phase(A,s,omega)
            +p.half(p.binary(binary.J0(e,s,omega)+binary.j0_homotopy(P,de,s,omega))))


def A_phase(A,s,omega):
    import low_phases,low_calibration
    epsilon=low_calibration.constants()['epsilon']
    linear=low_phases.transported(low_phases.pontryagin(omega),A,s)
    return reduced_A_phase(A,s,omega)-p.scale(linear,F(epsilon,4))


def shifted_comparison_gauge(A,B,s,omega):
    P,q,k,e=source(A,s,omega);bp=p.binary(B+e);lam=shifted_lambda(A,B,s,omega)
    Q=binary.comparison_terms(bp,k,lam,s,omega)[3]
    I=upper.hp.interval;Ai,si,wi=I(A),I(s),I(omega)
    bpt=I(bp,True);Bt=p.binary(bpt+I(e));u=p.binary(p.differential(bpt));y=p.QD(bpt,si,wi)
    qi,ki=I(q),I(k);c0=p.binary(ki+y+upper.hp.hD(qi,u,si))
    L=p.binary(shifted_lambda(Ai,Bt,si,wi)+p.Q(bpt,1)+p.cup(si,bpt)+p.cup(qi,u,3))
    trans=odd.phase_theta(qi,ki,u,y,si,wi)+odd.C_transport(c0,L,wi)
    return p.half(Q)-upper.hp.prism(trans)+p.half(p.Q(lam,1))


def comparison_gauge(A,B,C,s,omega):
    P,q,k,e=source(A,s,omega);I=upper.hp.interval
    ramp=closed_difference(I(A),I(B,True),I(C,True),I(s),I(omega))
    return (shifted_comparison_gauge(A,B,s,omega)+upper.hp.prism(ramp)
            +p.half(binary.j0_homotopy(B,e,s,omega)))


def reduced_production_difference(A,B,C,s,omega):
    """Omega4-Phi_uniform with only the linear A calibration removed."""
    P,q,k,e=source(A,s,omega);bp=p.binary(B+e);y=p.QD(bp,s,omega)
    lam=shifted_lambda(A,B,s,omega);dl=p.binary(p.differential(lam));tau=p.binary(y+k+dl)
    pol=p.binary(p.cup(y,k,3)+p.cup(k,p.binary(p.differential(y)),4)+p.Q(k,1))
    core=p.half(p.binary(pol+p.E(lam,omega)+p.cup(tau,dl,3)
                   +p.E(p.Q(bp,1),omega)+p.cup(s,tau)))
    I=upper.hp.interval;bpt=I(bp,True);si,wi=I(s),I(omega)
    Pi=upper.hp.prism(p.theta(p.binary(p.differential(bpt)),p.QD(bpt,si,wi),si,wi))
    return (p.half(p.E(C,omega))+core-Pi-odd.reduced_auxiliary(A,B,C,s,omega)
            -source_prism(A,s,omega))
