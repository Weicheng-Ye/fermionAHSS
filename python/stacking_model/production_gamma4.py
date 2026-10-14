"""Full legal production stacking correction in physical degree k=4.

The rational source is resolved by the fixed shared-sign universal chain
contractor, not by solving a cochain equation on X. Coefficients are certified
on all three mixed small source cells. Symmetry is not yet asserted here.
"""
from fractions import Fraction as F
from functools import lru_cache
import upper_pair_source as auxiliary
import production_gamma4_comparison as comparison
import production_upper_binary_comparison as binary
import upper_phase_diagnostic
p,upper,rc=auxiliary.p,auxiliary.upper,auxiliary.rc1


def source(A,Aprime,s,omega):
    return (auxiliary.source(A,Aprime,s,omega)+comparison.A_phase(A,s,omega)
            +comparison.A_phase(Aprime,s,omega)-comparison.A_phase(A+Aprime,s,omega)
            +p.half(binary.affine_J0_source(A,Aprime,s,omega)))


@lru_cache(None)
def source_value(pair):
    if not any(pair[0].left) or not any(pair[0].right):
        return F(0)
    A,Aprime,s=auxiliary.from_model1(pair[0])
    return source(A,Aprime,s,auxiliary.v1.from_omega(pair[1]))(tuple(range(6)))


def evaluate(chain):
    return sum(c*source_value(x) for x,c in chain.items())


@lru_cache(None)
def certificate():
    lower,upper_basis=rc.total_basis(4),rc.total_basis(5)
    values=tuple(evaluate(rc.Gtot(b)) for b in upper_basis)
    rows=tuple(tuple(rc.total_boundary(b).get(a,0) for a in lower) for b in upper_basis)
    if rows!=((-2,0),(-2,0),(0,-2)):
        raise ArithmeticError('the fixed shared-sign source basis changed')
    coefficients=(F(1,8),F(-1,8))
    if any((sum(a*b for a,b in zip(row,coefficients))-v)%1 for row,v in zip(rows,values)):
        raise ArithmeticError('production source has an uncancelled mixed period')
    return dict(lower=lower,upper=upper_basis,rows=rows,values=values,coefficients=coefficients)


def primitive(A,Aprime,s,omega):
    if (A.degree,Aprime.degree,s.degree,omega.degree)!=(1,1,1,2):
        raise ValueError('production gamma4 expects integral degree-one A inputs')
    small=dict(zip(rc.total_basis(4),(F(1,8),F(-1,8))))
    def value(vertices):
        diag=auxiliary.v1.to_diag(A,Aprime,s,vertices)
        x=(rc.Borel(diag.signs,diag.left,diag.right),auxiliary.v1.to_omega(omega,vertices))
        return evaluate(rc.Htot(x))+sum(c*small.get(b,0) for b,c in rc.Ftot(x).items())
    return p.Cochain(4,value)


def lower_product(A,B,C,Aprime,Bprime,Cprime,s,omega):
    alpha=upper.hp.hD(p.binary(A),p.binary(Aprime),s)
    beta=upper.beta_sharp(A,B,Aprime,Bprime,s,omega)
    return A+Aprime,p.binary(B+Bprime+alpha),p.binary(C+Cprime+beta)


def phase4(A,B,C,Aprime,Bprime,Cprime,s,omega):
    """The legal production phase; both inputs satisfy the first three equations."""
    I=upper.hp.interval
    ramp=auxiliary.reduced_full_source(I(A),I(B,True),I(C,True),I(Aprime),I(Bprime,True),I(Cprime,True),I(s),I(omega))
    As,Bs,Cs=lower_product(A,B,C,Aprime,Bprime,Cprime,s,omega)
    Qdiff=(comparison.comparison_gauge(A,B,C,s,omega)+comparison.comparison_gauge(Aprime,Bprime,Cprime,s,omega)
           -comparison.comparison_gauge(As,Bs,Cs,s,omega))
    return (p.scale(primitive(A,Aprime,s,omega),-1)-upper.hp.prism(ramp)
            -p.half(binary.affine_J0_primitive(A,B,Aprime,Bprime,s,omega))-Qdiff)


def omega4(A,B,C,s,omega):
    return upper_phase_diagnostic.production_phase(1,A,B,C,s,omega)


def gamma4(A,B,C,Aprime,Bprime,Cprime,s,omega):
    As,Bs,Cs=lower_product(A,B,C,Aprime,Bprime,Cprime,s,omega)
    raw=(omega4(A,B,C,s,omega)+omega4(Aprime,Bprime,Cprime,s,omega)-omega4(As,Bs,Cs,s,omega)
         +p.ds(phase4(A,B,C,Aprime,Bprime,Cprime,s,omega),s))
    def value(vertices):
        q=F(raw(vertices))
        if q.denominator!=1:
            raise ArithmeticError(f'production gamma4 is not integral: {q}')
        return q.numerator
    return p.Cochain(5,value)


if __name__=='__main__':
    import time
    for b in rc.total_basis(5):
        t=time.time();v=evaluate(rc.Gtot(b))
        print(b,v,v%1,'seconds',time.time()-t,flush=True)
