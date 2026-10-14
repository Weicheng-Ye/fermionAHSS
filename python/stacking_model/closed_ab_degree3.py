"""Actual k3 production stacking on legal A,B and arbitrary C,D.

The odd-anchor paths preserve the C curvature. Their transport uses the
closed potential polarization, retaining its pure-C successor phase.
On fully legal data this is literally the earlier even/odd legal rule.
"""
from fractions import Fraction as F
from cochains import Cochain,binary_sum,cup,differential,signed_differential,zero
from compatible_sector import integral,polarization
from lower_stacking import prism
import odd_a_gamma3 as odd
import even_a_gamma3 as even


def curvature(x,s,omega):
    return binary_sum(differential(x.C).mod2(),odd.tau(x.A,x.B,s,omega))


def potential(x,s,omega):
    tau=odd.tau(x.A,x.B,s,omega)
    t=binary_sum(differential(x.C).mod2(),tau)
    return odd.omega3(x.A,x.B,x.C,s,omega)+F(1,2)*cup(t,tau,2)


def source(left,right,s,omega):
    total=odd.stack(left,right,s,omega)
    return (potential(left,s,omega)+potential(right,s,omega)-potential(total,s,omega)
            +F(1,2)*polarization(curvature(left,s,omega),curvature(right,s,omega)))


def odd_phase(left,right,s,omega):
    def anchor_value(f):
        if left.A(f)%2:return left.A(f)
        if right.A(f)%2:return right.A(f)
        raise ValueError('an odd-anchor component needs at least one odd A')
    anchor=Cochain(0,anchor_value)
    ell=Cochain(0,lambda f:int(anchor(f)<0))
    unit=Cochain(0,lambda f:(-1)**ell(f))
    U=Cochain(1,lambda f:left.B(f)if left.A(f[:1])%2 else right.B(f))
    li,si,wi=odd.orient_path(left,ell,s,omega)
    ri,_,_=odd.orient_path(right,ell,s,omega)
    orient=prism(source(li,ri,si,wi))
    lo,ro=odd.endpoint_triple(li,1),odd.endpoint_triple(ri,1)
    lf,sf,wf=odd.flatten_path(lo,U,omega)
    rf,_,_=odd.flatten_path(ro,U,omega)
    flatten=prism(source(lf,rf,sf,wf))
    le,re=odd.endpoint_triple(lf,1),odd.endpoint_triple(rf,1)
    endpoint=even_phase(odd.Triple(zero(0),le.B,le.C),
                        odd.Triple(zero(0),re.B,re.C),zero(1),zero(2))
    return orient+cup(unit,flatten+endpoint,integral=True)


def even_phase(left,right,s,omega):
    canonical=even.canonical_data(left.A,left.B,left.C,s,omega)
    tau= binary_sum(canonical.H,canonical.z)
    return (even.phase3(left.A,left.B,left.C,right.A,right.B,right.C,s,omega)
            +F(1,2)*cup(tau,curvature(right,s,omega),3))


def phase(A,B,C,Ap,Bp,Cp,s,omega):
    """Componentwise even/odd phase, with no condition on either C."""
    args=tuple(Cochain(c.degree,c)for c in (A,B,C,Ap,Bp,Cp,s,omega))
    A,B,C,Ap,Bp,Cp,s,omega=args
    left,right=odd.Triple(A,B,C),odd.Triple(Ap,Bp,Cp)
    even_formula=even_phase(left,right,s,omega)
    odd_formula=odd_phase(left,right,s,omega)
    return Cochain(3,lambda f:(odd_formula if A(f[:1])%2 or Ap(f[:1])%2
                              else even_formula)(f))


def gamma(A,B,C,Ap,Bp,Cp,s,omega):
    args=tuple(Cochain(c.degree,c)for c in (A,B,C,Ap,Bp,Cp,s,omega))
    A,B,C,Ap,Bp,Cp,s,omega=args
    left,right=odd.Triple(A,B,C),odd.Triple(Ap,Bp,Cp)
    return integral(source(left,right,s,omega)+signed_differential(phase(*args),s),
                    'closed-AB k3 production gamma')
