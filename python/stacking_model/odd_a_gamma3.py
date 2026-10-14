"""Legal k3 phase when at least one signed A0 is odd.

A shared odd anchor trivializes s and omega. Explicit legal cylinders
transport the production Omega to a pure-C endpoint. No cochain equation
on X is solved. This is a legal-domain construction, not an off-shell one.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from cochains import Cochain,binary_sum,cup,differential,signed_differential,zero
from compatible_sector import E,Q,alpha,integral,polarization
from lower_stacking import pullback_interval as pull,prism
from n0_beta import beta_legal
import even_a_gamma3 as even
import h_tau_primitive as hp


@dataclass(frozen=True)
class Triple:
    A: Cochain
    B: Cochain
    C: Cochain


def tau(A,B,s,omega):
    answer=hp.low.secondary(*(hp.p.Cochain(c.degree,c) for c in (A,B,s,omega)))
    return Cochain(answer.degree,answer)


def endpoint(c,t):
    return Cochain(c.degree,lambda vertices:c(tuple((v,t) for v in vertices)))


def interval_contract(c):
    """Contraction to t=0 on X×I: delta J+J delta=1-pull ev0."""
    def value(vertices):
        mapped=tuple((v[0][0],v[0][1]*v[1]) for v in vertices)
        if any(a==b for a,b in zip(mapped,mapped[1:])):
            return 0
        return c(mapped)
    return prism(Cochain(c.degree,value))


def stack(left,right,s,omega):
    return Triple(left.A+right.A,
                  binary_sum(left.B,right.B,alpha(left.A,right.A,s)),
                  binary_sum(left.C,right.C,beta_legal(left.A,left.B,right.A,right.B,s,omega)))


def omega3(A,B,C,s,omega):
    """Exact production Omega3, rewriting its internal sign prism explicitly."""
    even_part=even.omega3(A,B,C,s,omega)
    absolute=Cochain(0,lambda f:abs(A(f)))
    t=integral(F(1,2)*(absolute-absolute.mod2()),'odd absolute carry').mod2()
    lam=binary_sum(cup(B,B),cup(t,Q(B,1)))
    cube=cup(cup(B,B,integral=True),B,integral=True)
    coefficient=Cochain(0,lambda f:F(3+2*t(f),8))
    oriented=-F(1,2)*E(lam,omega)+cup(coefficient,differential(cube),integral=True)
    ell=Cochain(0,lambda f:int(A(f)<0))
    gauge=pull(ell,True)
    si=binary_sum(pull(s),differential(gauge).mod2())
    ai=Cochain(0,lambda f:(-1)**gauge(f)*pull(A)(f))
    psi=tau(ai,pull(B),si,pull(omega))
    correction=F(1,2)*prism(E(psi,pull(omega))).mod2()
    ordinary=F(1,2)*E(C,omega)
    def value(vertices):
        if A(vertices[:1])%2==0:
            return even_part(vertices)
        return (ordinary(vertices)+(-1)**ell(vertices[:1])*oriented(vertices)
                +correction(vertices))
    return Cochain(4,value)


def orient_path(x,ell,s,omega):
    gauge=pull(ell,True)
    si=binary_sum(pull(s),differential(gauge).mod2())
    ai=Cochain(0,lambda f:(-1)**gauge(f)*pull(x.A)(f))
    bi=pull(x.B)
    ci=binary_sum(pull(x.C),interval_contract(tau(ai,bi,si,pull(omega))).mod2())
    return Triple(ai,bi,ci),si,pull(omega)


def flatten_lambda(A,B,U,omega):
    """delta lambda=Tau0 at s=0, delta U=omega, delta B=(rho A)omega."""
    a=A.mod2()
    t=integral(F(1,2)*(A-a),'signed A carry').mod2()
    bbar=binary_sum(B,cup(a,U))
    carry=integral(F(1,2)*(differential(U)-omega),'exact omega carry').mod2()
    return binary_sum(cup(U,B),cup(cup(a,bbar),U),cup(t,carry))


def flatten_path(x,U,omega):
    a=x.A.mod2()
    bbar=binary_sum(x.B,cup(a,U))
    cbar=binary_sum(x.C,flatten_lambda(x.A,x.B,U,omega))
    ui=binary_sum(pull(U),pull(U,True))
    wi=differential(ui).mod2()
    ai=pull(x.A)
    bi=binary_sum(pull(bbar),cup(pull(a),ui))
    ci=binary_sum(pull(cbar),flatten_lambda(ai,bi,ui,wi))
    return Triple(ai,bi,ci),zero(1),wi


def endpoint_triple(x,t):
    return Triple(*(endpoint(c,t) for c in (x.A,x.B,x.C)))


def omega_difference(left,right,s,omega):
    total=stack(left,right,s,omega)
    return (omega3(left.A,left.B,left.C,s,omega)+omega3(right.A,right.B,right.C,s,omega)
            -omega3(total.A,total.B,total.C,s,omega))


def phase3(A,B,C,Aprime,Bprime,Cprime,s,omega):
    """Explicit legal e3 with delta_s e3=-Delta Omega modulo one."""
    left,right=Triple(A,B,C),Triple(Aprime,Bprime,Cprime)
    def anchor_value(f):
        if A(f)%2:
            return A(f)
        if Aprime(f)%2:
            return Aprime(f)
        raise ValueError('odd-A phase requires at least one odd A on every component')
    anchor=Cochain(0,anchor_value)
    ell=Cochain(0,lambda f:int(anchor(f)<0))
    u=Cochain(0,lambda f:(-1)**ell(f))
    U=Cochain(1,lambda f:B(f) if A(f[:1])%2 else Bprime(f))
    li,si,wi=orient_path(left,ell,s,omega)
    ri,_,_=orient_path(right,ell,s,omega)
    orient_phase=prism(omega_difference(li,ri,si,wi))
    lo,ro=endpoint_triple(li,1),endpoint_triple(ri,1)
    lf,sf,wf=flatten_path(lo,U,omega)
    rf,_,_=flatten_path(ro,U,omega)
    flat_phase=prism(omega_difference(lf,rf,sf,wf))
    le,re=endpoint_triple(lf,1),endpoint_triple(rf,1)
    final_phase=F(1,2)*polarization(le.C,re.C)
    return orient_phase+cup(u,flat_phase+final_phase,integral=True)


def gamma3(A,B,C,Aprime,Bprime,Cprime,s,omega):
    phase=phase3(A,B,C,Aprime,Bprime,Cprime,s,omega)
    return integral(omega_difference(Triple(A,B,C),Triple(Aprime,Bprime,Cprime),s,omega)
                    +signed_differential(phase,s),'odd-A legal gamma3')
