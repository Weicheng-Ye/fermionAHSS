"""Exact integer contractions for iterated binary diagonal bar models.

Adapted from fermionAHSS/python/r3_chain.py and chain_models.py under their
MIT license and original attribution. Unlike the integral divided-power
specialization, the full nonzero differential of every inner binary small
complex is retained. No production file or calibration is modified.
"""
from dataclasses import dataclass
from functools import cached_property, lru_cache
from itertools import product

from h_tau_primitive import cm

add, linear, sign = cm.add, cm.linear, cm.sign


@dataclass(frozen=True)
class Simplex:
    level: int
    data: tuple

    @cached_property
    def _cached_hash(self):
        return hash((self.level,self.data))

    def __hash__(self):
        return self._cached_hash

    def __post_init__(self):
        if self.level < 1:
            raise ValueError('binary bar level must be positive')
        if self.level == 1:
            if any(x not in (0,1) for x in self.data):
                raise ValueError('binary edges required')
        elif any(not isinstance(x,Simplex) or x.level != self.level-1
                 or len(x.data) != len(self.data) for x in self.data):
            raise ValueError('each diagonal bar row must have the matching degree')


@dataclass(frozen=True)
class Rectangle:
    level: int
    q: int
    rows: tuple

    @cached_property
    def _cached_hash(self):
        return hash((self.level,self.q,self.rows))

    def __hash__(self):
        return self._cached_hash


@lru_cache(None)
def identity(level,q):
    return Simplex(level,(0,)*q if level==1 else (identity(level-1,q),)*q)


def multiply(x,y):
    if x.level!=y.level or len(x.data)!=len(y.data):
        raise ValueError('simplicial group multiplication needs equal dimensions')
    return Simplex(x.level,tuple(a^b for a,b in zip(x.data,y.data)) if x.level==1
                   else tuple(multiply(a,b) for a,b in zip(x.data,y.data)))


def is_identity(x):
    return not any(x.data) if x.level==1 else all(is_identity(row) for row in x.data)


def degree(x):
    if isinstance(x,Simplex):
        return len(x.data)
    if isinstance(x,(cm.Atom,cm.Diag)):
        return cm._degree(x)
    if not x:
        return 0
    ds={degree(y) for y in x}
    if len(ds)!=1:
        raise ValueError('simplicial product factors have unequal degrees')
    return ds.pop()


@lru_cache(None)
def degens(x):
    if isinstance(x,Simplex):
        if x.level==1:
            return frozenset(i for i,b in enumerate(x.data) if not b)
        return frozenset(i for i,row in enumerate(x.data)
                         if is_identity(row) and all(i in degens(y) for y in x.data))
    if isinstance(x,(cm.Atom,cm.Diag)):
        return frozenset(cm._degen_indices(x))
    result=set(range(degree(x)))
    for y in x:
        result.intersection_update(degens(y))
    return frozenset(result)


def basis(x):
    return {} if degens(x) else {x:1}


def h_pull(x,vertices):
    rows=[]
    for start,stop in zip(vertices,vertices[1:]):
        row=identity(x.level-1,x.q)
        for other in x.rows[start:stop]:
            row=multiply(row,other)
        rows.append(row)
    return Rectangle(x.level,x.q,tuple(rows))


def v_pull(x,vertices):
    return Rectangle(x.level,len(vertices)-1,tuple(pull(row,vertices)[0] for row in x.rows))


def rect_pull(x,hv,vv):
    return v_pull(h_pull(x,hv),vv)


@lru_cache(None)
def pull(x,vertices):
    if isinstance(x,Simplex):
        if x.level==1:
            return Simplex(1,tuple(sum(x.data[a:b])%2
                                  for a,b in zip(vertices,vertices[1:]))),1
        rect=rect_pull(Rectangle(x.level,degree(x),x.data),vertices,vertices)
        return Simplex(x.level,rect.rows),1
    if isinstance(x,(cm.Atom,cm.Diag)):
        return cm._pull(x,vertices)
    rows=[];coefficient=1
    for row in x:
        y,c=pull(row,vertices);rows.append(y);coefficient*=c
    return tuple(rows),coefficient


@lru_cache(None)
def boundary(x):
    if degens(x):
        return {}
    n,out=degree(x),{}
    for i in range(n+1) if n else ():
        y,c=pull(x,tuple(j for j in range(n+1) if j!=i))
        add(out,basis(y),sign(i)*c)
    return out


def AW(pair):
    x,y=pair;n,out=degree(x),{}
    if degree(y)!=n:
        raise ValueError('AW requires equal degrees')
    for p in range(n+1):
        a,c=pull(x,tuple(range(p+1)));b,d=pull(y,tuple(range(p,n+1)))
        if not degens(a) and not degens(b):
            add(out,{(a,b):c*d})
    return out


def shuffle(pair):
    x,y=pair;out={}
    if degens(x) or degens(y):
        return out
    for hv,vv,c in cm._shuffles(degree(x),degree(y)):
        a,ca=pull(x,hv);b,cb=pull(y,vv)
        add(out,basis((a,b)),c*ca*cb)
    return out


def product_H(pair):
    if degens(pair):
        return {}
    out={}
    for (hv,vv),c in cm._h_terms(degree(pair)).items():
        a,ca=pull(pair[0],hv);b,cb=pull(pair[1],vv)
        add(out,basis((a,b)),c*ca*cb)
    return out


def rect_degenerate(x):
    return any(is_identity(row) for row in x.rows) or any(
        all(i in degens(row) for row in x.rows) for i in range(x.q))


def rect_basis(x):
    return {} if rect_degenerate(x) else {x:1}


def h_boundary(x):
    p,out=len(x.rows),{}
    for i in range(p+1) if p else ():
        add(out,rect_basis(h_pull(x,tuple(j for j in range(p+1) if j!=i))),sign(i))
    return out


@lru_cache(None)
def FD(x):
    n,out=degree(x),{}
    if degens(x):
        return out
    rect=Rectangle(x.level,n,x.data)
    for p in range(n+1):
        add(out,rect_basis(rect_pull(rect,tuple(range(p+1)),tuple(range(p,n+1)))))
    return out


@lru_cache(None)
def GD(x):
    if rect_degenerate(x):
        return {}
    out={}
    for hv,vv,c in cm._shuffles(len(x.rows),x.q):
        add(out,basis(Simplex(x.level,rect_pull(x,hv,vv).rows)),c)
    return out


@lru_cache(None)
def HD(x):
    if degens(x):
        return {}
    n,out=degree(x),{};rect=Rectangle(x.level,n,x.data)
    for (hv,vv),c in cm._h_terms(n).items():
        add(out,basis(Simplex(x.level,rect_pull(rect,hv,vv).rows)),c)
    return out


def multi_AW(rows):
    if len(rows)<2:
        return {rows:1} if not rows or not degens(rows[0]) else {}
    out={}
    for (left,right),c in AW((rows[:-1],rows[-1])).items():
        for word,d in multi_AW(left).items():
            add(out,{word+(right,):c*d})
    return out


def multi_shuffle(word):
    if len(word)<2:
        return {word:1} if not word or not degens(word[0]) else {}
    out={}
    for left,c in multi_shuffle(word[:-1]).items():
        for (a,b),d in shuffle((left,word[-1])).items():
            add(out,{a+(b,):c*d})
    return out


def multi_H(rows):
    if len(rows)<2:
        return {}
    out={};pair=(rows[:-1],rows[-1])
    for (left,right),c in product_H(pair).items():
        add(out,{left+(right,):c})
    for (left,right),c in AW(pair).items():
        for a,d in multi_H(left).items():
            for (b,e),f in shuffle((a,right)).items():
                add(out,{b+(e,):c*d*f})
    return out


def rect_f0(x):
    return {word:c for word,c in multi_AW(x.rows).items()
            if not any(degree(a)==0 for a in word)}


def rect_g0(level,word):
    if not word:
        return {Rectangle(level,0,()):1}
    out={}
    for rows,c in multi_shuffle(word).items():
        add(out,rect_basis(Rectangle(level,degree(rows),rows)),c)
    return out


def rect_h0(x):
    out={}
    for rows,c in multi_H(x.rows).items():
        add(out,rect_basis(Rectangle(x.level,x.q+1,rows)),sign(len(x.rows))*c)
    return out


def bar_sign(word):
    p=len(word)
    return sign(p+sum((p-j-1)*degree(a) for j,a in enumerate(word)))


@lru_cache(None)
def FV(x):
    out,term={},{x:1}
    for _ in range(len(x.rows)+1):
        for y,c in term.items():
            for word,d in rect_f0(y).items():
                add(out,{word:c*d*bar_sign(word)})
        term={a:-c for a,c in linear(linear(term,rect_h0),h_boundary).items()}
        if not term:
            return out
    raise ArithmeticError('vertical projection perturbation did not terminate')


@lru_cache(None)
def HV(x):
    out,term={},{x:1}
    for _ in range(len(x.rows)+1):
        add(out,linear(term,rect_h0))
        term={a:-c for a,c in linear(linear(term,rect_h0),h_boundary).items()}
        if not term:
            return out
    raise ArithmeticError('vertical homotopy perturbation did not terminate')


@lru_cache(None)
def Fbar(x):
    return linear(FD(x),FV)


@lru_cache(None)
def Gbar(level,word):
    return linear({x:bar_sign(word)*c for x,c in rect_g0(level,word).items()},GD)


@lru_cache(None)
def Hbar(x):
    return add(dict(HD(x)),linear(linear(FD(x),HV),GD))


def unit_p(x):
    return F(x) if degree(x) else {}


def unit_Q(x):
    return add(basis(x),linear(unit_p(x),lambda b:G(x.level,b)),-1) if degree(x) else {}


def unit_u(x):
    return linear(linear(unit_Q(x),H),unit_Q)


@lru_cache(None)
def unit_h(x):
    return linear(linear(unit_u(x),boundary),unit_u)


def tensor_p(word):
    out={():1}
    for x in word:
        nxt={}
        for prefix,c in out.items():
            for letter,d in unit_p(x).items():
                add(nxt,{prefix+(letter,):c*d})
        out=nxt
    return out


def tensor_i(level,word):
    out={():1}
    for letter in word:
        nxt={}
        for prefix,c in out.items():
            for x,d in G(level-1,letter).items():
                add(nxt,{prefix+(x,):c*d})
        out=nxt
    return out


def tensor_h(word):
    out,prefix,prefix_degree={},{():1},0
    for j,x in enumerate(word):
        for left,c in prefix.items():
            for y,d in unit_h(x).items():
                add(out,{left+(y,)+word[j+1:]:sign(prefix_degree+1)*c*d})
        nxt={}
        for left,c in prefix.items():
            for y,d in linear(unit_p(x),lambda b:G(x.level,b)).items():
                add(nxt,{left+(y,):c*d})
        prefix=nxt;prefix_degree+=degree(x)+1
    return out


@lru_cache(None)
def atom_product(x,y):
    out={}
    for hv,vv,c in cm._shuffles(degree(x),degree(y)):
        z=multiply(pull(x,hv)[0],pull(y,vv)[0])
        if degree(z):
            add(out,basis(z),c)
    return out


def bar_delta(word):
    out={};internal_degree=0
    for j in range(len(word)-1):
        internal_degree+=degree(word[j])
        for x,c in atom_product(word[j],word[j+1]).items():
            add(out,{word[:j]+(x,)+word[j+2:]:sign(j+internal_degree)*c})
    return out


def bar_boundary(word):
    out=bar_delta(word);internal_degree=0
    for j,x in enumerate(word):
        for y,c in boundary(x).items():
            add(out,{word[:j]+(y,)+word[j+1:]:sign(j+1+internal_degree)*c})
        internal_degree+=degree(x)
    return out


@lru_cache(None)
def transfer_F(word):
    out,term={},{word:1}
    for _ in range(len(word)+1):
        add(out,linear(term,tensor_p))
        term={a:-c for a,c in linear(linear(term,tensor_h),bar_delta).items()}
        if not term:
            return out
    raise ArithmeticError('binary projection transfer did not terminate')


@lru_cache(None)
def transfer_G(level,word):
    out,term={},tensor_i(level,word)
    for _ in range(len(word)+1):
        add(out,term)
        term={a:-c for a,c in linear(linear(term,bar_delta),tensor_h).items()}
        if not term:
            return out
    raise ArithmeticError('binary inclusion transfer did not terminate')


@lru_cache(None)
def transfer_H(word):
    out,term={},{word:1}
    for _ in range(len(word)+1):
        add(out,linear(term,tensor_h))
        term={a:-c for a,c in linear(linear(term,tensor_h),bar_delta).items()}
        if not term:
            return out
    raise ArithmeticError('binary homotopy transfer did not terminate')


@lru_cache(None)
def small_boundary(level,word):
    if level==1:
        return {word-1:2} if word>0 and word%2==0 else {}
    # The inner binary differential is NOT zero: retain the full boundary.
    return linear(linear(transfer_G(level,word),bar_boundary),tensor_p)


@lru_cache(None)
def F(x):
    if x.level==1:
        return {} if degens(x) else {degree(x):1}
    return linear(Fbar(x),transfer_F)


@lru_cache(None)
def G(level,word):
    if level==1:
        return {Simplex(1,(1,)*word):1}
    return linear(transfer_G(level,word),lambda w:Gbar(level,w))


@lru_cache(None)
def H(x):
    if x.level==1:
        return {}
    return add(dict(Hbar(x)),linear(linear(Fbar(x),transfer_H),lambda w:Gbar(x.level,w)))


def small_degree(level,word):
    return word if level==1 else sum(small_degree(level-1,a)+1 for a in word)


@lru_cache(None)
def small_basis(level,n):
    if level==1:
        return (n,)
    if n==0:
        return ((),)
    return tuple((letter,)+tail for d in range(1,n)
                 for letter in small_basis(level-1,d)
                 for tail in small_basis(level,n-d-1))


def cochain_value(x,vertices):
    """Rectangular-sum cocycle of degree level on a diagonal bar simplex."""
    if x.level==1:
        return sum(x.data[vertices[0]:vertices[1]])%2
    return sum(cochain_value(row,vertices[1:])
               for row in x.data[vertices[0]:vertices[1]])%2


def coordinate(level,source,vertices):
    """Natural binary alternating section, recursively contracting each edge."""
    q=len(vertices)-1
    def normalized(indices):
        if len(set(indices))!=len(indices):
            return 0
        return source(tuple(vertices[i] for i in sorted(indices)))%2
    def build(m,value):
        if m==1:
            return Simplex(1,tuple(value((i,i+1)) for i in range(q)))
        rows=[]
        for i in range(q):
            def next_value(indices,i=i,value=value):
                def evaluate(v):
                    return 0 if len(set(v))!=len(v) else value(tuple(sorted(v)))
                return evaluate((i,)+indices)^evaluate((i+1,)+indices)
            rows.append(build(m-1,next_value))
        return Simplex(m,tuple(rows))
    return build(level,normalized)
