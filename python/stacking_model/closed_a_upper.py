"""Explicit auxiliary differential on signed-closed A and arbitrary B,C.

This preserves the calibrated Tau on legal B, and supplies an integral top
output using a fixed production A-only source primitive. Matching this
auxiliary top representative to production on legal inputs is a separate
comparison; no complete upper stacking law is asserted here.
"""
from fractions import Fraction as F
from functools import lru_cache
import h_tau_primitive as hp
p,low=hp.p,hp.low


@lru_cache(None)
def _high_source_value(pair):
    import high_phase as high
    rc,rs=high.rc,high.rs
    coefficients=high.source_coefficients()
    return (rs.evaluate_r(rs.phi(),rc.Htot(pair))
            +sum(c*coefficients.get(b,0) for b,c in rc.Ftot(pair).items()))


def fsharp(A,B,s,omega):
    src=hp.uniform_source(A,s,omega)
    P=src['p']
    total=src['g']+p.cup(s,B)
    lam=p.binary(p.divide(total-p.binary(total),2,'splitting carry'))
    return p.binary(p.QD(B,s,omega)+src['k0']+p.differential(lam)
                    +hp.hD(p.binary(p.differential(B)),P,s)
                    +p.differential(p.binary(p.Q(B,1)+p.cup(s,B))))


def beta_sharp(A,B,Aprime,Bprime,s,omega):
    import stacking_lower
    alpha=hp.hD(p.binary(A),p.binary(Aprime),s)
    q=p.binary(p.differential(alpha))
    u,up=p.binary(p.differential(B)),p.binary(p.differential(Bprime))
    b=p.binary(u+p.QD(p.binary(A),s,omega))
    bp=p.binary(up+p.QD(p.binary(Aprime),s,omega))
    base=stacking_lower.legal_beta(A,B,Aprime,Bprime,s,omega)
    m=A.degree+2
    return p.binary(base+p.cup(b,up,m)+p.cup(q,p.binary(b+bp),m))



def hE(x,y):
    return p.binary(p.cup(x,y,x.degree-1)
                    +p.cup(p.binary(p.differential(x)),y,x.degree))


def primary(A,s,omega):
    return p.QD(p.binary(A),s,omega)


def curvature(A,B,C,s,omega):
    return p.binary(p.differential(B)+primary(A,s,omega)),p.binary(p.differential(C)+fsharp(A,B,s,omega))


def source_primitive(A,s,omega):
    """V_P with delta_s V_P=Theta(P,k0) modulo integers, n=1,2,3."""
    n=A.degree
    if n not in (1,2,3):
        raise ValueError('fixed A-only source primitive covers n=1,2,3')
    @lru_cache(None)
    def value(vertices):
        registered=p.make_source(n,len(vertices)-1,
             lambda f:s(tuple(vertices[i]for i in f)),
             lambda f:A(tuple(vertices[i]for i in f)),
             lambda f:omega(tuple(vertices[i]for i in f)))
        if n==1:
            return low.V1(registered)
        if n==2:
            return low.V2_pair(p.to_diags(registered))
        import high_phase as high
        return _high_source_value(high.rs.to_diags(registered))
    V=p.Cochain(n+4,value)
    if n!=1:
        return V
    data=hp.uniform_source(A,s,omega)
    a=p.binary(A);t=p.binary(p.divide(A-a,2,'A carry'))
    e0=p.cup(s,t)
    eI=hp.interval(e0,True)
    deI=p.binary(p.differential(eI))
    si,wi=hp.interval(s),hp.interval(omega)
    P,k=hp.interval(data['p']),hp.interval(data['k0'])
    shiftedP=p.binary(P+deI)
    shiftedk=p.binary(k+p.QD(eI,si,wi)+hp.hD(deI,P,si))
    return V-hp.prism(p.theta(shiftedP,shiftedk,si,wi))


def phi(A,B,C,s,omega):
    """Natural potential satisfying delta_s Phi=-Theta(curvature), mod1."""
    Ai,si,wi=hp.interval(A),hp.interval(s),hp.interval(omega)
    Bi=hp.interval(B,True)
    bi=p.binary(p.differential(Bi)+primary(Ai,si,wi))
    fi=fsharp(Ai,Bi,si,wi)
    f=fsharp(A,B,s,omega)
    cross=p.half(p.binary(p.E(C,omega)+hE(p.binary(p.differential(C)),f)))
    return cross-source_primitive(A,s,omega)-hp.prism(p.theta(bi,fi,si,wi))


def top(A,B,C,s,omega):
    b,c=curvature(A,B,C,s,omega)
    raw=p.ds(phi(A,B,C,s,omega),s)+p.theta(b,c,s,omega)
    def value(vertices):
        q=F(raw(vertices))
        if q.denominator!=1:
            raise ArithmeticError(f'auxiliary signed-closed-A top output is not integral: {q}')
        return q.numerator
    return p.Cochain(A.degree+5,value)
