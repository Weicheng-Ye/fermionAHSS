"""Legal odd-rank exchange from explicit orientation and flattening prisms."""
from fractions import Fraction as F
from cochains import Cochain,binary_sum,cup,differential,signed_differential
from compatible_sector import E,integral,polarization
from lower_stacking import prism
import odd_a_gamma3 as odd


def lower_gauge(A,B,Aprime,Bprime,s):
    a=A.mod2();ap=Aprime.mod2()
    t=integral(F(1,2)*(A-a),'A parity carry').mod2()
    tp=integral(F(1,2)*(Aprime-ap),'Aprime parity carry').mod2()
    return binary_sum(cup(B,Bprime,1),cup(a,Bprime),
                      cup(binary_sum(cup(t,tp),cup(a,tp)),s))


def transport(left,right,s,omega):
    total=odd.stack(left,right,s,omega)
    lam=lower_gauge(left.A,left.B,right.A,right.B,s)
    return F(1,2)*binary_sum(E(lam,omega),polarization(total.C,differential(lam).mod2()))


def phase_correction(A,B,Aprime,Bprime):
    coefficient=Cochain(0,lambda f:((A(f)*A(f)-1)//8)%2
                         if A(f)%2 and Aprime(f)%2 else 0)
    b=binary_sum(B,Bprime)
    return F(1,2)*cup(coefficient,cup(cup(b,b),b))


def phase(A,B,C,Aprime,Bprime,Cprime,s,omega):
    return (odd.phase3(A,B,C,Aprime,Bprime,Cprime,s,omega)
            +phase_correction(A,B,Aprime,Bprime))


def gamma(A,B,C,Aprime,Bprime,Cprime,s,omega):
    left,right=odd.Triple(A,B,C),odd.Triple(Aprime,Bprime,Cprime)
    return integral(odd.omega_difference(left,right,s,omega)
                    +signed_differential(phase(A,B,C,Aprime,Bprime,Cprime,s,omega),s),
                    'recalibrated odd-rank gamma3')


def residual(left,right,s,omega):
    return (phase(right.A,right.B,right.C,left.A,left.B,left.C,s,omega)
            -phase(left.A,left.B,left.C,right.A,right.B,right.C,s,omega)
            -transport(left,right,s,omega))


def flat_primitive(left,right):
    """Primitive when the chosen odd anchor has B=0 and s=omega=0."""
    def correction(face):
        anchor,other=(left,right)if left.A(face[:1])%2 else(right,left)
        if other.A(face[:1])%2==0:return 0
        av,apv=anchor.A(face[:1]),other.A(face[:1])
        coefficient=((apv+1)//2)%2
        ta=((av-1)//2)%2;tabs=((abs(apv)-1)//2)%2
        orient=ta*tabs+int(apv<0)*(1+tabs)
        g,h=other.B(face[:2]),other.B(face[1:])
        value=coefficient*g*h*anchor.C(face)+orient*(g+h)*(anchor.C(face)+other.C(face))
        return F(value%2,2)
    return F(1,2)*cup(left.C,right.C,2)+Cochain(2,correction)


def phase_primitive(left,right,s,omega):
    def anchor_value(face):
        if left.A(face)%2:return left.A(face)
        if right.A(face)%2:return right.A(face)
        raise ValueError('require an odd rank on every component')
    anchor=Cochain(0,anchor_value)
    ell=Cochain(0,lambda f:int(anchor(f)<0))
    u=Cochain(0,lambda f:(-1)**ell(f))
    U=Cochain(1,lambda f:left.B(f)if left.A(f[:1])%2 else right.B(f))
    li,si,wi=odd.orient_path(left,ell,s,omega)
    ri,_,_=odd.orient_path(right,ell,s,omega)
    orient=prism(residual(li,ri,si,wi))
    lo,ro=odd.endpoint_triple(li,1),odd.endpoint_triple(ri,1)
    lf,sf,wf=odd.flatten_path(lo,U,omega)
    rf,_,_=odd.flatten_path(ro,U,omega)
    flatten=prism(residual(lf,rf,sf,wf))
    le,re=odd.endpoint_triple(lf,1),odd.endpoint_triple(rf,1)
    return -orient+cup(u,flat_primitive(le,re)-flatten,integral=True)


def exchange(A,B,C,Aprime,Bprime,Cprime,s,omega):
    left,right=odd.Triple(A,B,C),odd.Triple(Aprime,Bprime,Cprime)
    total=odd.stack(left,right,s,omega)
    opposite=odd.stack(right,left,s,omega)
    lam=lower_gauge(A,B,Aprime,Bprime,s)
    K=transport(left,right,s,omega)
    sigma=phase_primitive(left,right,s,omega)
    L=integral(odd.omega3(opposite.A,opposite.B,opposite.C,s,omega)
               -odd.omega3(total.A,total.B,total.C,s,omega)
               -signed_differential(K,s),'odd-rank C transport')
    M=integral(residual(left,right,s,omega)-signed_differential(sigma,s),
               'odd-rank exchange carry')
    return lam,K,sigma,L,M
