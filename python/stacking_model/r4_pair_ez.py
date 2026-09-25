"""Common-sign Borel plus equivariant fiber-EZ contraction.

Actual reduced K(Z,n) cocycle simplices are used. The two-fiber tensor
complex has no relative cells below degree 2n. Exact integer EZ operators
come from fermionAHSS chain_models.py under its MIT license/attribution.
"""
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
import r2_pair_chain as r2
cm,add,linear,sign=r2.cm,r2.add,r2.linear,r2.sign


@lru_cache(None)
def faces(q,n):
    return tuple(combinations(range(q+1),n+1))


@dataclass(frozen=True)
class KZ:
    n:int
    q:int
    values:tuple

    def __post_init__(self):
        if self.n<1 or self.q<0 or len(self.values)!=len(faces(self.q,self.n)):
            raise ValueError('invalid K(Z,n) simplex dimensions')
        if any(not isinstance(a,int) for a in self.values):
            raise ValueError('integral K(Z,n) simplex required')
        values=dict(zip(faces(self.q,self.n),self.values))
        for f in combinations(range(self.q+1),self.n+2):
            if sum(sign(i)*values[f[:i]+f[i+1:]] for i in range(len(f))):
                raise ValueError('K(Z,n) simplex must be a cocycle')


@dataclass(frozen=True)
class Borel:
    signs:tuple
    left:KZ
    right:KZ

    def __post_init__(self):
        if self.left.q!=self.right.q or self.left.n!=self.right.n or any(s not in (0,1) for s in self.signs):
            raise ValueError('invalid shared-sign Borel rectangle')


@dataclass(frozen=True)
class Small:
    signs:tuple
    left:KZ
    right:KZ


def negate(a):
    return KZ(a.n,a.q,tuple(-v for v in a.values))


@lru_cache(None)
def kz_pull(a,vertices):
    values=dict(zip(faces(a.q,a.n),a.values))
    return KZ(a.n,len(vertices)-1,tuple(0 if len(set(v))!=a.n+1 else values[v]
              for f in faces(len(vertices)-1,a.n)
              for v in (tuple(vertices[i] for i in f),)))


@lru_cache(None)
def kz_degens(a):
    out=set()
    for i in range(a.q):
        restriction=kz_pull(a,tuple(j for j in range(a.q+1) if j!=i))
        surjection=tuple(j if j<=i else j-1 for j in range(a.q+1))
        if kz_pull(restriction,surjection)==a:
            out.add(i)
    return out


def h_pull(x,vertices):
    invert=sum(x.signs[:vertices[0]])%2
    return Borel(tuple(sum(x.signs[a:b])%2 for a,b in zip(vertices,vertices[1:])),
                 negate(x.left) if invert else x.left,negate(x.right) if invert else x.right)


def v_pull(x,vertices):
    return Borel(x.signs,kz_pull(x.left,vertices),kz_pull(x.right,vertices))


def rect_pull(x,hv,vv):
    return v_pull(h_pull(x,hv),vv)


def degree(x):
    if isinstance(x,Borel):
        if len(x.signs)!=x.left.q:
            raise ValueError('diagonal expected')
        return x.left.q
    if isinstance(x,KZ):
        return x.q
    if isinstance(x,(cm.Atom,cm.Diag)):
        return cm._degree(x)
    return degree(x[0])


def degens(x):
    if isinstance(x,Borel):
        return {i for i,s in enumerate(x.signs) if not s}&kz_degens(x.left)&kz_degens(x.right)
    if isinstance(x,KZ):
        return kz_degens(x)
    if isinstance(x,(cm.Atom,cm.Diag)):
        return cm._degen_indices(x)
    return degens(x[0])&degens(x[1])


def basis(x):
    return {} if degens(x) else {x:1}


def pull(x,vertices):
    if isinstance(x,Borel):
        return rect_pull(x,vertices,vertices),1
    if isinstance(x,KZ):
        return kz_pull(x,vertices),1
    if isinstance(x,(cm.Atom,cm.Diag)):
        return cm._pull(x,vertices)
    a,c=pull(x[0],vertices);b,d=pull(x[1],vertices)
    return (a,b),c*d


def boundary(x):
    if degens(x):
        return {}
    n,out=degree(x),{}
    for i in range(n+1) if n else ():
        y,c=pull(x,tuple(j for j in range(n+1) if j!=i))
        add(out,basis(y),sign(i)*c)
    return out


def product_AW(pair):
    x,y=pair;n,out=degree(x),{}
    for p in range(n+1):
        a,c=pull(x,tuple(range(p+1)));b,d=pull(y,tuple(range(p,n+1)))
        if not degens(a) and not degens(b):
            add(out,{(a,b):c*d})
    return out


def product_shuffle(pair):
    x,y=pair;out={}
    if degens(x) or degens(y):
        return out
    for hv,vv,c in cm._shuffles(degree(x),degree(y)):
        a,ca=pull(x,hv);b,cb=pull(y,vv)
        add(out,basis((a,b)),c*ca*cb)
    return out


def product_homotopy(pair):
    if degens(pair):
        return {}
    out={}
    for (hv,vv),c in cm._h_terms(degree(pair)).items():
        a,ca=pull(pair[0],hv);b,cb=pull(pair[1],vv)
        add(out,basis((a,b)),c*ca*cb)
    return out


def rect_degenerate(x):
    return any(not s for s in x.signs) or bool(kz_degens(x.left)&kz_degens(x.right))


def rect_basis(x):
    return {} if rect_degenerate(x) else {x:1}


@lru_cache(None)
def FD(x):
    if degens(x):
        return {}
    n,out=degree(x),{}
    for p in range(n+1):
        add(out,rect_basis(rect_pull(x,tuple(range(p+1)),tuple(range(p,n+1)))))
    return out


@lru_cache(None)
def GD(x):
    if rect_degenerate(x):
        return {}
    out={}
    for h,v,c in cm._shuffles(len(x.signs),x.left.q):
        add(out,basis(rect_pull(x,h,v)),c)
    return out


@lru_cache(None)
def HD(x):
    if degens(x):
        return {}
    out={}
    for (h,v),c in cm._h_terms(degree(x)).items():
        add(out,basis(rect_pull(x,h,v)),c)
    return out


def F0(x):
    if rect_degenerate(x):
        return {}
    return {Small(x.signs,a,b):c for (a,b),c in product_AW((x.left,x.right)).items()}


def G0(small):
    return {Borel(small.signs,a,b):c for (a,b),c in product_shuffle((small.left,small.right)).items()}


def H0(x):
    if rect_degenerate(x):
        return {}
    return {Borel(x.signs,a,b):sign(len(x.signs))*c
            for (a,b),c in product_homotopy((x.left,x.right)).items()}


@lru_cache(None)
def F(x):
    return linear(FD(x),F0)


@lru_cache(None)
def G(small):
    return linear(G0(small),GD)


@lru_cache(None)
def H(x):
    return add(dict(HD(x)),linear(linear(FD(x),H0),GD))


def small_basis(x):
    return {} if any(not s for s in x.signs) or degens(x.left) or degens(x.right) else {x:1}


def small_boundary(x):
    p,out=len(x.signs),{}
    for i in range(p+1) if p else ():
        vertices=tuple(j for j in range(p+1) if j!=i)
        invert=sum(x.signs[:vertices[0]])%2
        y=Small(tuple(sum(x.signs[a:b])%2 for a,b in zip(vertices,vertices[1:])),
                negate(x.left) if invert else x.left,negate(x.right) if invert else x.right)
        add(out,small_basis(y),sign(i))
    for a,c in boundary(x.left).items():
        add(out,small_basis(Small(x.signs,a,x.right)),sign(p)*c)
    for b,c in boundary(x.right).items():
        add(out,small_basis(Small(x.signs,x.left,b)),sign(p+x.left.q)*c)
    return out


def Ftot(pair):
    out={}
    for (x,y),c in product_AW(pair).items():
        for a,d in F(x).items():
            for b,e in cm.F(y).items():
                add(out,{(a,b):c*d*e})
    return out


def Gtot(b):
    small,omega=b;out={}
    for x,c in G(small).items():
        for y,d in cm.G('c2',omega).items():
            add(out,product_shuffle((x,y)),c*d)
    return out


def Htot(pair):
    out=dict(product_homotopy(pair))
    for (x,y),c in product_AW(pair).items():
        for a,d in H(x).items():
            add(out,product_shuffle((a,y)),c*d)
        for a,d in linear(F(x),G).items():
            for b,e in cm.H(y).items():
                add(out,product_shuffle((a,b)),c*sign(degree(x))*d*e)
    return out
