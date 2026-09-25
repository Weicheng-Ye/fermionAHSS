"""Experimental A-only source for the signed-closed-A auxiliary top potential.

Fixed universal evaluations only. Production comparison is not asserted.
"""
from fractions import Fraction as F
from functools import lru_cache
import closed_a_upper as upper
import r1_pair_chain_signed as rc1
import v1_pair_shared as v1
from theta_pair_phase import ThetaPairPhase
from cochains import Cochain
p=upper.p


@lru_cache(None)
def theta_phase(n):
    return ThetaPairPhase(n+2)


def as_cochain(c):
    return Cochain(c.degree,c)


def source(A,Aprime,s,omega):
    n=A.degree
    B,C=p.zero(n+1),p.zero(n+2)
    alpha=upper.hp.hD(p.binary(A),p.binary(Aprime),s)
    beta=upper.beta_sharp(A,B,Aprime,B,s,omega)
    P,Pprime=upper.primary(A,s,omega),upper.primary(Aprime,s,omega)
    k,kprime=upper.fsharp(A,B,s,omega),upper.fsharp(Aprime,B,s,omega)
    theta=theta_phase(n).phase(*(as_cochain(c)for c in(P,k,Pprime,kprime,s,omega)))
    return (upper.phi(A,B,C,s,omega)+upper.phi(Aprime,B,C,s,omega)
            -upper.phi(A+Aprime,alpha,beta,s,omega)-p.Cochain(theta.degree,theta))


def from_model1(diag):
    ordinary=v1.rc.Borel(diag.signs,diag.left,diag.right)
    return v1.from_diag(ordinary)


@lru_cache(None)
def value1(pair):
    if not any(pair[0].left) or not any(pair[0].right):
        return F(0)
    A,Aprime,s=from_model1(pair[0])
    return source(A,Aprime,s,v1.from_omega(pair[1]))(tuple(range(6)))


def evaluate1(chain):
    return sum(c*value1(x) for x,c in chain.items())


if __name__=='__main__':
    import time
    for b in rc1.total_basis(5):
        started=time.time()
        value=evaluate1(rc1.Gtot(b))
        print(b,'value',value,'mod1',value%1,'seconds',time.time()-started,flush=True)


@lru_cache(None)
def certificate1():
    lower,upper_basis=rc1.total_basis(4),rc1.total_basis(5)
    coefficients=(F(1,8),F(1,8))
    indices={b:i for i,b in enumerate(lower)}
    values=tuple(evaluate1(rc1.Gtot(b)) for b in upper_basis)
    rows=tuple(tuple(rc1.total_boundary(b).get(a,0) for a in lower) for b in upper_basis)
    if any((sum(a*c for a,c in zip(row,coefficients))-value)%1
           for row,value in zip(rows,values)):
        raise ArithmeticError('the fixed auxiliary upper pair calibration failed')
    return dict(lower=lower,upper=upper_basis,rows=rows,values=values,coefficients=coefficients)


def primitive1(A,Aprime,s,omega):
    if (A.degree,Aprime.degree,s.degree,omega.degree)!=(1,1,1,2):
        raise ValueError('auxiliary upper pair primitive expects degrees (1,1,1,2)')
    cert=certificate1();small=dict(zip(cert['lower'],cert['coefficients']))
    def value(vertices):
        diag=v1.to_diag(A,Aprime,s,vertices)
        x=(rc1.Borel(diag.signs,diag.left,diag.right),v1.to_omega(omega,vertices))
        return evaluate1(rc1.Htot(x))+sum(c*small.get(b,0) for b,c in rc1.Ftot(x).items())
    return p.Cochain(4,value)


def full_source(A,B,C,Aprime,Bprime,Cprime,s,omega):
    n=A.degree
    alpha=upper.hp.hD(p.binary(A),p.binary(Aprime),s)
    beta=upper.beta_sharp(A,B,Aprime,Bprime,s,omega)
    b,c=upper.curvature(A,B,C,s,omega)
    bp,cp=upper.curvature(Aprime,Bprime,Cprime,s,omega)
    theta=theta_phase(n).phase(*(as_cochain(x)for x in(b,c,bp,cp,s,omega)))
    return (upper.phi(A,B,C,s,omega)+upper.phi(Aprime,Bprime,Cprime,s,omega)
            -upper.phi(A+Aprime,p.binary(B+Bprime+alpha),p.binary(C+Cprime+beta),s,omega)
            -p.Cochain(theta.degree,theta))


def phase4(A,B,C,Aprime,Bprime,Cprime,s,omega):
    """Auxiliary phase for signed-closed A,A' and arbitrary B,C."""
    I=upper.hp.interval
    along=full_source(I(A),I(B,True),I(C,True),I(Aprime),I(Bprime,True),I(Cprime,True),I(s),I(omega))
    return p.scale(primitive1(A,Aprime,s,omega),-1)-upper.hp.prism(along)


def gamma4(A,B,C,Aprime,Bprime,Cprime,s,omega):
    source_value=full_source(A,B,C,Aprime,Bprime,Cprime,s,omega)
    raw=source_value+p.ds(phase4(A,B,C,Aprime,Bprime,Cprime,s,omega),s)
    def value(vertices):
        q=F(raw(vertices))
        if q.denominator!=1:
            raise ArithmeticError(f'auxiliary gamma4 is not integral: {q}')
        return q.numerator
    return p.Cochain(5,value)


def phi_without_A_source(A,B,C,s,omega):
    """Phi+V_P, for prisms with the A data pulled back constantly."""
    I=upper.hp.interval
    Ai,Bi,si,wi=I(A),I(B,True),I(s),I(omega)
    bi=p.binary(p.differential(Bi)+upper.primary(Ai,si,wi))
    fi=upper.fsharp(Ai,Bi,si,wi);f=upper.fsharp(A,B,s,omega)
    return (p.half(p.binary(p.E(C,omega)+upper.hE(p.binary(p.differential(C)),f)))
            -upper.hp.prism(p.theta(bi,fi,si,wi)))


def reduced_full_source(A,B,C,Aprime,Bprime,Cprime,s,omega):
    """Same prism as full_source when A,Aprime,s,omega are constant along I."""
    n=A.degree
    alpha=upper.hp.hD(p.binary(A),p.binary(Aprime),s)
    beta=upper.beta_sharp(A,B,Aprime,Bprime,s,omega)
    b,c=upper.curvature(A,B,C,s,omega)
    bp,cp=upper.curvature(Aprime,Bprime,Cprime,s,omega)
    theta=theta_phase(n).phase(*(as_cochain(x)for x in(b,c,bp,cp,s,omega)))
    return (phi_without_A_source(A,B,C,s,omega)+phi_without_A_source(Aprime,Bprime,Cprime,s,omega)
            -phi_without_A_source(A+Aprime,p.binary(B+Bprime+alpha),p.binary(C+Cprime+beta),s,omega)
            -p.Cochain(theta.degree,theta))
