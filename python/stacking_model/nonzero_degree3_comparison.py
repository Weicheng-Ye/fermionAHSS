"""Explicit k3 potential comparison on legal A,B and unrestricted C.

The nonzero rank is oriented by its sign. Odd rank contracts B together
with omega; even rank uses the fixed BC2 times K(F2,2) contraction. The
per-rank five-cell check protects the fixed R0/R1 suspension calibration.
"""
from functools import lru_cache
from fractions import Fraction as F
import off_shell_beta as off
from upper_phase_diagnostic import production_phase
p,hp,cm=off.p,off.hp,off.hp.cm


def raw_potential(A,B,C,s,omega):
    I=hp.interval
    ai,bi,ci,si,wi=I(A,True),I(B,True),I(C,True),I(s),I(omega)
    a,b=off.curvature(ai,bi,si,wi)
    c=p.binary(p.differential(ci)+off.raw_H(ai,bi,si,wi))
    return hp.prism(production_phase(1,a,b,c,si,wi))


def potential(A,B,C,s,omega):
    tau=hp.low.secondary(A,B,s,omega)
    t=p.binary(p.differential(C)+tau)
    return production_phase(0,A,B,C,s,omega)+p.half(p.cup(t,tau,2))


def difference(A,B,C,s,omega):
    """Closed Q/Z source Psi=Phi_J-Phi_hat-h(s t)."""
    tau=hp.low.secondary(A,B,s,omega)
    t=p.binary(p.differential(C)+tau)
    return potential(A,B,C,s,omega)-raw_potential(A,B,C,s,omega)-p.half(p.cup(s,t))


def omega_section(vertices,omega):
    n=len(vertices)-1
    matrix=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            value=omega((vertices[i],vertices[i+1],vertices[j+1]))
            if j>i+1:value+=omega((vertices[i],vertices[i+1],vertices[j]))
            matrix[i][j]=matrix[j][i]=value%2
    return cm.Diag('c2',tuple((0,tuple(row))for row in matrix))


EVEN_BASES=((0,(3,)),(0,(1,1)),(1,(2,)),(2,(1,)),(4,()))
EVEN_LOWER_BASES=((0,(2,)),(1,(1,)),(3,()))


@lru_cache(None)
def even_source(rank,pair):
    b,w=pair
    B=p.Cochain(1,lambda f:sum(b.edges[f[0]:f[1]])%2)
    omega=p.Cochain(2,lambda f:sum(w.rows[i][1][j]
                    for i in range(f[0],f[1])for j in range(f[1],f[2]))%2)
    A=p.Cochain(0,lambda f:rank)
    return difference(A,B,p.zero(2),p.zero(1),omega)(tuple(range(5)))


@lru_cache(None)
def even_inclusion(basis):
    degree,word=basis
    b=cm.Atom('c2',0,(1,)*degree)
    return cm.linear({(b,w):c for w,c in cm.G('c2',word).items()},cm.product_shuffle)


@lru_cache(None)
def even_certificate(rank):
    if rank<0 or rank%2:
        raise ValueError('the oriented even-rank certificate needs rank>=0 even')
    values=tuple(sum(c*even_source(rank,pair)for pair,c in even_inclusion(basis).items())
                 for basis in EVEN_BASES)
    # The first two cells are cycles. The next two have the same boundary
    # 2*(1,(1,)); the last has boundary 2*(3,()). Only these three periods
    # must vanish, not the value on every noncycle simplex.
    if any(F(v)%1 for v in (values[0],values[1],values[2]-values[3])):
        raise ArithmeticError(f'the R0/R1 rank comparison has a residual period: {values}')
    return values


@lru_cache(None)
def even_small_coefficients(rank):
    values=even_certificate(rank)
    return dict(zip(EVEN_LOWER_BASES,(F(0),(values[2]%1)/2,(values[4]%1)/2)))


@lru_cache(None)
def even_forward(pair):
    result={}
    for (b,w),c in cm.product_AW(pair).items():
        if cm._degen_indices(b):continue
        for word,d in cm.F(w).items():
            cm.add(result,{(cm._degree(b),word):c*d})
    return result


@lru_cache(None)
def even_homotopy(pair):
    result=dict(cm.product_homotopy(pair))
    for (b,w),c in cm.product_AW(pair).items():
        for v,d in cm.H(w).items():
            cm.add(result,cm.product_shuffle((b,v)),c*(-1)**cm._degree(b)*d)
    return result


def even_oriented_gauge(A,B,C,omega):
    """The rank is an arbitrary nonnegative even locally constant integer."""
    I=hp.interval
    cpart=hp.prism(difference(I(A),I(B),I(C,True),p.zero(1),I(omega)))
    def source_value(vertices):
        rank=A(vertices[:1])
        coefficients=even_small_coefficients(rank)
        pair=(cm.Atom('c2',0,tuple(B((a,b))for a,b in zip(vertices,vertices[1:]))),
              omega_section(vertices,omega))
        return (sum(c*even_source(rank,q)for q,c in even_homotopy(pair).items())
                +sum(c*coefficients.get(q,0)for q,c in even_forward(pair).items()))
    return cpart+p.Cochain(3,source_value)


def odd_oriented_gauge(A,B,C,omega):
    I=hp.interval
    bi=p.binary(I(B)+I(B,True))
    ci=p.binary(I(C)+I(C,True))
    wi=p.binary(p.differential(bi))
    return p.scale(hp.prism(difference(I(A),bi,ci,p.zero(1),wi)),-1)


@lru_cache(None)
def zero_rank_constructor():
    from a0_degree3 import DegreeThreeA0Stacking
    return DegreeThreeA0Stacking(hp.ROOT)


def gauge(A,B,C,s,omega):
    """q3 with Phi_J-Phi_hat=h(s t)+delta_s q3, modulo integers.

    Requires delta_s A=0 and delta B=QD(rho A) globally. C is arbitrary.
    Componentwise rank-zero, even, and odd cases are selected at vertices;
    legal A makes each case constant on each connected component.
    """
    if (A.degree,B.degree,C.degree)!=(0,1,2):
        raise ValueError('the k3 comparison expects cochain degrees0,1,2')
    I=hp.interval
    ell=p.Cochain(0,lambda f:int(A(f)<0))
    gi=I(ell,True)
    ai=p.Cochain(0,lambda f:(-1)**gi(f)*I(A)(f))
    si=p.binary(I(s)+p.differential(gi))
    orientation=hp.prism(difference(ai,I(B),I(C),si,I(omega)))
    absolute=p.Cochain(0,lambda f:abs(A(f)))
    even=even_oriented_gauge(absolute,B,C,omega)
    odd=odd_oriented_gauge(absolute,B,C,omega)
    zero_cache=[]
    def value(vertices):
        rank=A(vertices[:1])
        if not rank:
            if not zero_cache:
                from cochains import Cochain
                args=tuple(Cochain(c.degree,c)for c in (B,C,s,omega))
                zero_cache.append(zero_rank_constructor().closed_pair_gauge(*args))
            return zero_cache[0](vertices)
        oriented=odd if rank%2 else even
        return (-1)**int(rank<0)*oriented(vertices)-orientation(vertices)
    return p.Cochain(3,value)
