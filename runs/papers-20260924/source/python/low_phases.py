"""Danus phases in input degrees 0, 1 and 2, with fixed universal choices.

Only the explicitly displayed d_s D terms in the final assemblers are
omitted. Their omission leaves the integral differential unchanged.
All source contractions are universal; no cochain equation on the input
space is solved here. See low_calibration.py for the fixed finite selectors.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations

try:
    from . import phase_eval as p, chain_models as cm, source_primitive as sp
    from . import odd_comparison as oc, odd_primitive as op
except ImportError:
    import phase_eval as p
    import chain_models as cm
    import source_primitive as sp
    import odd_comparison as oc
    import odd_primitive as op


SOURCE_VALUES = tuple(map(Fraction, ('0', '3/4', '0', '3/4', '1/4')))
ODD_LAMBDAS = (0,)*8


def transported(left, right, s):
    """AW multiplication with transport of the sign-valued right factor."""
    degree = left.degree
    return p.Cochain(degree+right.degree,
        lambda z: left(z[:degree+1])*(-1)**s((z[0], z[degree]))*right(z[degree:]))


def pontryagin(omega):
    return p.cup(omega, omega, integral=True) + p.cup(
        omega, p.differential(omega), 1, integral=True)


def secondary(A, b, s, omega):
    """Literal current calibrated psi, including full integral carries."""
    n = A.degree
    a = p.binary(A)
    t = p.binary(p.divide(A-a, 2, 'A carry'))
    B = p.divide(p.differential(a), 2, 'binary Bockstein')
    e = p.binary(B)
    CB = p.binary(p.divide(B+e, 2, 'plus carry'))
    u,v,w = p.square(a,2),p.cup(omega,a),p.cup(s,e)
    FA = p.binary(p.zeta2(omega,a)+p.chi(a)
        +p.cup(u,v,n+1)+p.cup(u,w,n+1)+p.cup(v,w,n+1)
        +p.zeta1(s,e)+p.cup(p.cup(omega,s,1),e)
        +p.cup(s,u)+p.cup(p.cup(s,s),CB))
    F = p.binary(p.E(b,omega)+FA)
    g = p.binary(p.Q(t,2)+p.cup(omega,t)
        +p.cup(p.binary(p.differential(t)),p.cup(s,a),n)
        +p.zeta1(s,a)+p.cup(p.cup(omega,s,1),a)+p.cup(s,b))
    q = p.cup(omega,B,integral=True)+p.cup(B,B,n-1,integral=True)
    L = p.divide(q-p.ds(g,s),2,'secondary L')
    return p.binary(F+p.binary(L)+p.cup(p.cup(p.cup(s,s),s),a))


def boundary(b,s,omega):
    return p.theta(b,p.zero(b.degree+1),s,omega)


def R0(A,b,s,omega):
    if A.degree != 0 or b.degree != 1:
        raise ValueError('R0 expects degrees (0,1)')
    K = p.divide(A,2,'even input')
    t = p.binary(K)
    quarter = p.divide(K-t,2,'quarter input carry')
    Z = p.cup(t,omega)
    sb = p.cup(s,b)
    Nb = p.divide(Z+sb-p.binary(Z+sb),2,'binary sum carry')
    lam = p.binary(p.cup(quarter,omega,integral=True)+Nb)
    H = p.QD(b,s,omega)
    r = p.binary(p.divide(p.differential(omega),2,'background Bockstein'))
    z = p.cup(t,r)
    even = boundary(b,s,omega)-p.scale(p.cup(K,pontryagin(omega),integral=True),Fraction(1,8))
    even = even+p.half(p.binary(p.cup(H,z,2)+p.E(lam,omega)
            +p.cup(p.binary(H+z),p.binary(p.differential(lam)),2)))

    absolute = p.Cochain(0,lambda z: abs(A(z)))
    odd_t = p.binary(p.divide(absolute-p.binary(absolute),2,'odd input carry'))
    odd_lam = p.binary(p.cup(b,b)+p.cup(odd_t,p.Q(b,1)))
    dbcube = p.differential(p.cup(p.cup(b,b,integral=True),b,integral=True))
    oriented = p.scale(p.E(odd_lam,omega),Fraction(-1,2))
    oriented = oriented+p.Cochain(4,lambda z: Fraction(3+2*odd_t((z[0],)),8)*dbcube(z))
    ell = p.Cochain(0,lambda z:int(A(z)<0))
    gauge = p.Cochain(0,lambda z:z[0][-1]*ell(z))
    sI = p.binary(s+p.differential(gauge))
    AI = p.Cochain(0,lambda z:(-1)**gauge(z)*A(z))
    psiI = secondary(AI,b,sI,omega)
    residualI = p.E(psiI,omega)
    prism = p.prism_cochain(residualI)
    def evaluate(vertices):
        value = A((vertices[0],))
        if value % 2 == 0:
            return even(vertices)
        return (-1)**ell((vertices[0],))*oriented(vertices)+Fraction(prism(vertices)%2,2)
    return p.Cochain(4,evaluate)


# Homogeneous D_infinity chains.  Coordinate 2k-e identifies its Cayley tree
# with Z; actual reverse bar edges remain distinct chain generators.
def _coordinate(g): return 2*g[0]-g[1]
def _vertex(q): return ((q+(q%2))//2,q%2)
def _gdegen(g): return any(a==b for a,b in zip(g,g[1:]))
def _gbasis(g): return {} if _gdegen(g) else {g:1}


def _gboundary(g):
    out={}
    for i in range(len(g)) if len(g)>1 else ():
        cm.add(out,_gbasis(g[:i]+g[i+1:]),(-1)**i)
    return out


def _cone(root,chain):
    out={}
    for g,c in chain.items(): cm.add(out,_gbasis((root,)+g),c)
    return out


def _K(root,g):
    if len(g)==1:
        x,y=_coordinate(root),_coordinate(g[0]); step=1 if y>x else -1
        return {(_vertex(i),_vertex(i+step)):1 for i in range(x,y,step)}
    ends=set(g)
    if len(ends)!=2 or abs(_coordinate(g[0])-_coordinate(g[1]))!=1:
        raise ArithmeticError('Cayley contractor received a simplex outside the wedge')
    near=min(ends,key=lambda v:abs(_coordinate(v)-_coordinate(root)))
    return _gbasis((near,)+g)


@lru_cache(None)
def retract_group(g):
    if _gdegen(g): return {}
    if len(g)==1: return {g:1}
    return cm.linear(cm.linear(_gboundary(g),retract_group),lambda x:_K(g[0],x))


@lru_cache(None)
def homotopy_group(g):
    if len(g)==1 or _gdegen(g): return {}
    out=cm.add({g:1},retract_group(g),-1)
    cm.add(out,cm.linear(_gboundary(g),homotopy_group),-1)
    return _cone(g[0],out)


def _pair_degen(pair):
    g,w=pair
    return {i for i in range(len(g)-1) if g[i]==g[i+1]} & cm._degen_indices(w)


def _pair_basis(pair): return {} if _pair_degen(pair) else {pair:1}


def _pair_pull(pair,h,v):
    g,w=pair
    wr,factor=cm._pull(w,v)
    return (tuple(g[i] for i in h),wr),factor


def group_product_AW(pair):
    g,w=pair; n=len(g)-1; out={}
    for i in range(n+1):
        q,c=_pair_pull(pair,tuple(range(i+1)),tuple(range(i,n+1)))
        if not _gdegen(q[0]) and not cm._degen_indices(q[1]): cm.add(out,{q:c})
    return out


def group_product_shuffle(pair):
    g,w=pair; out={}
    if _gdegen(g) or cm._degen_indices(w): return out
    for h,v,c in cm._shuffles(len(g)-1,len(w.rows)):
        q,f=_pair_pull(pair,h,v); cm.add(out,_pair_basis(q),c*f)
    return out


def group_product_H(pair):
    out={}
    if _pair_degen(pair): return out
    for (h,v),c in cm._raw_h_terms(len(pair[0])-1).items():
        q,f=_pair_pull(pair,h,v); cm.add(out,_pair_basis(q),c*f)
    return out


def _omega_diagonal(vertices,omega):
    n=len(vertices)-1; M=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            M[i][j]=(omega((vertices[i],vertices[i+1],vertices[j+1]))
                     +omega((vertices[i],vertices[i+1],vertices[j])))%2
            M[j][i]=M[i][j]
    return cm.Diag('c2',tuple((0,tuple(row)) for row in M))


def to_pair1(vertices,A,s,omega):
    g=[(0,0)]
    for i in range(1,len(vertices)):
        g.append((A((vertices[0],vertices[i])),s((vertices[0],vertices[i]))))
    return tuple(g),_omega_diagonal(vertices,omega)


@lru_cache(None)
def from_pair1(pair):
    g,w=pair; n=len(g)-1
    def sv(ij):
        i,j=ij; return g[i][1]^g[j][1]
    def av(ij):
        i,j=ij; return (-1)**g[i][1]*(g[j][0]-g[i][0])
    def ov(ijk):
        i,j,k=ijk; return sum(w.rows[u][1][v] for u in range(i,j) for v in range(j,k))%2
    return p.make_source(1,n,sv,av,ov)


@lru_cache(None)
def odd_raw():
    def cached_fill(gamma):
        @lru_cache(None)
        def value(pair):
            return sum(gamma(oc.from_pair(term)) for term in oc.T(pair))%2
        return p.Cochain(5,lambda vertices:value(oc.to_pair(vertices)))
    return op.build_Uraw(op.build_odd(p.A1,p.s,p.omega),ODD_LAMBDAS,cached_fill)


def _background_source(w,n=1):
    return p.make_source(n,len(w.rows),lambda ij:0,lambda ij:0,
        lambda ijk:sum(w.rows[u][1][v] for u in range(ijk[0],ijk[1])
                     for v in range(ijk[1],ijk[2]))%2)


def _background_value(cochain,w): return cochain(_background_source(w))


@lru_cache(None)
def odd_base_normalization():
    raw=odd_raw()
    R=p.binary(p.divide(p.differential(p.omega),2,'omega Bockstein'))
    correction=p.half(p.cup(p.omega,R))
    cycle={(4,):-1,(1,2):1,(2,1):-1}
    chain=cm.linear(cycle,lambda w:cm.G('c2',w))
    period=sum(c*_background_value(raw,w) for w,c in chain.items())
    bit=2*(period%1)
    if bit.denominator!=1 or int(bit) not in (0,1):
        raise ArithmeticError('odd relative normalization is not binary')
    bit=int(bit)
    phi=raw-p.scale(correction,bit)
    vals=[]
    for word in ((1,2),(2,1)):
        vals.append(sum(c*_background_value(phi,w) for w,c in cm.G('c2',word).items())%1)
    B,C=vals
    ef={(3,):(C-B)/2,(1,1):-(B+C)/4}
    @lru_cache(None)
    def E0_diagonal(diag):
        return sum(c*_background_value(phi,w) for w,c in cm.H(diag).items())+sum(
            c*ef.get(w,0) for w,c in cm.F(diag).items())
    E0=p.Cochain(4,lambda vertices:E0_diagonal(_omega_diagonal(vertices,p.omega)))
    return raw-p.scale(correction,bit)-p.ds(E0,p.s),dict(eta=bit,period=str(period),B=str(B),C=str(C))


@lru_cache(None)
def source_phase1():
    data=p.source(p.A1,p.s,p.omega)
    return p.theta(data['q'],data['k'],p.s,p.omega)


def _evaluate_pair(cochain,pair):
    return (-1)**pair[0][0][1]*cochain(from_pair1(pair))


def _wedge_primitive(pair):
    g,w=pair
    # Constant group simplices pull back Ut|s=0=0, as do the r branch.
    if len(set(g))<2: return Fraction(0)
    edges=[((-1)**a[1]*(b[0]-a[0]),a[1]^b[1]) for a,b in zip(g,g[1:]) if a!=b]
    if all(a==0 for a,s in edges): return Fraction(0)
    if any(edge!=(1,1) for edge in edges):
        raise ArithmeticError('retracted simplex does not lie on either dihedral branch')
    return _evaluate_pair(odd_base_normalization()[0],pair)


@lru_cache(None)
def V1_pair(pair):
    """Memoize the actual A,s,omega input, independently of b,c/source IDs."""
    C=group_product_H(pair)
    projected={}
    for (g,w),c in group_product_AW(pair).items():
        for h,d in homotopy_group(g).items(): cm.add(C,group_product_shuffle((h,w)),c*d)
        for r,d in retract_group(g).items(): cm.add(projected,group_product_shuffle((r,w)),c*d)
    return sum(c*_evaluate_pair(source_phase1(),q) for q,c in C.items())+sum(
        c*_wedge_primitive(q) for q,c in projected.items())


def V1(vertices,A=None,s=None,omega=None):
    A,s,omega=A or p.A1,s or p.s,omega or p.omega
    return V1_pair(to_pair1(vertices,A,s,omega))


@lru_cache(None)
def V2_pair(pair):
    """The fixed source depends only on its two diagonal input simplices."""
    ef=sp.ef_values(SOURCE_VALUES)
    return p.evaluate_r(p.phase2()['phi'],sp.Htot(pair))+sum(
        coefficient*ef.get(key,0) for key,coefficient in sp.Ftot(pair).items())


def _registered(n,vertices,A,b,c,s,omega):
    return p.make_source(n,len(vertices)-1,
        lambda f:s(tuple(vertices[i] for i in f)),
        lambda f:A(tuple(vertices[i] for i in f)),
        lambda f:omega(tuple(vertices[i] for i in f)),
        lambda f:b(tuple(vertices[i] for i in f)),
        lambda f:c(tuple(vertices[i] for i in f)))


def R1_raw(A,b,s,omega):
    data=p.source(A,s,omega)
    q0,g,k=data['p'],data['g'],data['k']
    e0=p.cup(s,data['t']); de0=p.binary(p.differential(e0))
    bp=p.binary(b+e0); y=p.QD(bp,s,omega)
    sb=p.cup(s,b)
    N=p.binary(p.divide(g+sb-p.binary(g+sb),2,'splitting carry'))
    lam=p.binary(N+p.cup(b,e0,1)+p.cup(q0,e0,2)+p.cup(q0,de0,3)
        +p.cup(s,p.cup(b,e0,2)+p.cup(q0,e0,3)+p.cup(q0,de0,4)))
    dlam=p.binary(p.differential(lam)); psi=p.binary(y+k+dlam)
    pol=p.binary(p.cup(y,k,3)+p.cup(k,p.binary(p.differential(y)),4)+p.Q(k,1))
    bI=p.Cochain(2,lambda z:bp(z)*z[-1][-1])
    prism=p.prism_cochain(p.theta(p.binary(p.differential(bI)),p.QD(bI,s,omega),s,omega))
    core=p.half(p.binary(pol+p.E(lam,omega)+p.cup(psi,dlam,3)
        +p.E(p.Q(bp,1),omega)+p.cup(s,psi)))
    V=p.Cochain(5,lambda z:V1(z,A,s,omega))
    return core-prism-V


def build_phase(n,A,b,c,s,omega,*,epsilon=None,alpha=None):
    """Canonical low-degree integral d5 phase, before applying d_s.

    epsilon/alpha are normally loaded from the independently evaluated finite
    selectors. Explicit values are accepted only to reproduce those selectors.
    """
    if n not in (0,1,2): raise ValueError('low phase input must be 0, 1 or 2')
    if (A.degree,b.degree,c.degree,s.degree,omega.degree)!=(n,n+1,n+2,1,2):
        raise ValueError('invalid defining-system degrees')
    if n==0: return p.half(p.E(c,omega))+R0(A,b,s,omega)
    if epsilon is None or (n==2 and alpha is None):
        try: from .low_calibration import constants
        except ImportError: from low_calibration import constants
        bits=constants()
        epsilon=bits['epsilon'] if epsilon is None else epsilon
        alpha=bits['alpha'] if alpha is None else alpha
    if n==1:
        return p.half(p.E(c,omega))+R1_raw(A,b,s,omega)-p.scale(
            transported(pontryagin(omega),A,s),Fraction(epsilon,4))
    base=p.Cochain(6,lambda z:p.phase2()['Xi'](_registered(n,z,A,b,c,s,omega))
                    -V2_pair(p.to_diags(_registered(n,z,A,b,c,s,omega))))
    # The three sign-system factors multiply in order; the middle transport
    # occurs at the start of the second A factor and the last at the third.
    cube=transported(transported(A,A,s),A,s)
    return base+p.scale(transported(pontryagin(omega),A,s),Fraction(alpha,4))-p.scale(cube,Fraction(1,4))
