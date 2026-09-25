"""Legal A=0 gamma5,gamma6 from reduced binary cocycles and product EZ.

All chain operations are finite, and use the MIT-attributed interval-cut EZ
operators from fermionAHSS. The universal pair source has no nonintegral
small period; no source cochain is solved on X. Full off-shell compatibility
is a separate step from the legal-data identity implemented here.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations

from cochains import Cochain, binary_sum, cup, differential, signed_differential
from compatible_sector import E, Q, QD, alpha, integral, polarization
from lower_stacking import pullback_interval, prism
import h_tau_primitive as production
p, cm = production.p, production.cm


@lru_cache(None)
def faces(q, m):
    return tuple(combinations(range(q+1), m+1))


@dataclass(frozen=True)
class KB:
    m: int
    q: int
    values: tuple

    def __post_init__(self):
        if self.m < 1 or self.q < 0 or len(self.values) != len(faces(self.q,self.m)):
            raise ValueError('invalid reduced binary cocycle simplex')
        if any(value not in (0,1) for value in self.values):
            raise ValueError('binary cocycle values must be zero or one')
        table = dict(zip(faces(self.q,self.m),self.values))
        for face in combinations(range(self.q+1),self.m+2):
            if sum(table[face[:i]+face[i+1:]] for i in range(len(face)))%2:
                raise ValueError('the binary input simplex is not closed')


@lru_cache(None)
def kb_pull(x, vertices):
    table = dict(zip(faces(x.q,x.m),x.values))
    values = []
    for face in faces(len(vertices)-1,x.m):
        pulled = tuple(vertices[i] for i in face)
        values.append(0 if len(set(pulled)) != x.m+1 else table[pulled])
    return KB(x.m,len(vertices)-1,tuple(values))


@lru_cache(None)
def kb_degens(x):
    result = set()
    for i in range(x.q):
        face = kb_pull(x,tuple(j for j in range(x.q+1) if j != i))
        degeneracy = tuple(j if j <= i else j-1 for j in range(x.q+1))
        if kb_pull(face,degeneracy) == x:
            result.add(i)
    return result


def degree(x):
    return x.q if isinstance(x,KB) else cm._degree(x) if isinstance(x,cm.Diag) else degree(x[0])


def degens(x):
    if isinstance(x,KB):
        return kb_degens(x)
    if isinstance(x,cm.Diag):
        return cm._degen_indices(x)
    return degens(x[0]) & degens(x[1])


def pull(x, vertices):
    if isinstance(x,KB):
        return kb_pull(x,vertices),1
    if isinstance(x,cm.Diag):
        return cm._pull(x,vertices)
    a,c = pull(x[0],vertices)
    b,d = pull(x[1],vertices)
    return (a,b),c*d


def basis(x):
    return {} if degens(x) else {x:1}


def boundary(x):
    if degens(x):
        return {}
    out = {}
    n = degree(x)
    for i in range(n+1) if n else ():
        face,c = pull(x,tuple(j for j in range(n+1) if j != i))
        cm.add(out,basis(face),(-1)**i*c)
    return out


def AW(x):
    out = {}
    n = degree(x)
    for i in range(n+1):
        a,c = pull(x[0],tuple(range(i+1)))
        b,d = pull(x[1],tuple(range(i,n+1)))
        if not degens(a) and not degens(b):
            cm.add(out,{(a,b):c*d})
    return out


def shuffle(pair):
    a,b = pair
    if degens(a) or degens(b):
        return {}
    out = {}
    for left,right,c in cm._shuffles(degree(a),degree(b)):
        aa,ca = pull(a,left)
        bb,cb = pull(b,right)
        cm.add(out,basis((aa,bb)),c*ca*cb)
    return out


def product_H(pair):
    if degens(pair):
        return {}
    out = {}
    for (left,right),c in cm._h_terms(degree(pair)).items():
        a,ca = pull(pair[0],left)
        b,cb = pull(pair[1],right)
        cm.add(out,basis((a,b)),c*ca*cb)
    return out


@lru_cache(None)
def forward(x):
    if isinstance(x,(KB,cm.Diag)):
        return basis(x)
    out = {}
    for (a,b),c in AW(x).items():
        for aa,ca in forward(a).items():
            for bb,cb in forward(b).items():
                cm.add(out,{(aa,bb):c*ca*cb})
    return out


@lru_cache(None)
def inclusion(x):
    if isinstance(x,(KB,cm.Diag)):
        return basis(x)
    out = {}
    for a,ca in inclusion(x[0]).items():
        for b,cb in inclusion(x[1]).items():
            cm.add(out,shuffle((a,b)),ca*cb)
    return out


@lru_cache(None)
def homotopy(x):
    if isinstance(x,(KB,cm.Diag)):
        return {}
    out = dict(product_H(x))
    for (a,b),c in AW(x).items():
        for aa,ca in homotopy(a).items():
            cm.add(out,shuffle((aa,b)),c*ca)
        for aa,ca in cm.linear(forward(a),inclusion).items():
            for bb,cb in homotopy(b).items():
                cm.add(out,shuffle((aa,bb)),c*(-1)**degree(a)*ca*cb)
    return out


class HigherA0Stacking:
    """The full legal-data integral correction in physical degree k=5 or6."""

    def __init__(self,k):
        if k not in (5,6):
            raise ValueError('this reduced source constructor covers k=5,6')
        self.k,self.m = k,k-2
        if k == 5 and self.minimal_period()%1:
            raise ArithmeticError('the fixed minimal source period is nonintegral')

    def boundary_phase(self,b,s,omega):
        if (b.degree,s.degree,omega.degree) != (self.m,1,2):
            raise ValueError('inconsistent input degrees')
        bi = pullback_interval(b,True)
        si,wi = pullback_interval(s),pullback_interval(omega)
        raw = p.theta(differential(bi).mod2(),QD(bi,si,wi),si,wi)
        suspended = prism(Cochain(raw.degree,raw))
        if self.m == 3:
            return suspended-F(1,2)*E(cup(s,b),omega)
        return -suspended+F(1,2)*binary_sum(E(Q(b,1),omega),cup(s,QD(b,s,omega)))

    def omega(self,b,c,s,omega):
        return F(1,2)*E(c,omega)+self.boundary_phase(b,s,omega)

    def source(self,b,bp,s,omega):
        beta = alpha(b,bp,s)
        tau,taup = QD(b,s,omega),QD(bp,s,omega)
        phase = (self.boundary_phase(binary_sum(b,bp),s,omega)
                 -self.boundary_phase(b,s,omega)-self.boundary_phase(bp,s,omega))
        return phase+F(1,2)*binary_sum(polarization(tau,taup),E(beta,omega),
            polarization(binary_sum(tau,taup),differential(beta).mod2()))

    def inputs(self,simplex):
        (b,bp),(sb,w) = simplex
        def binary_cochain(x):
            table = dict(zip(faces(x.q,x.m),x.values))
            return Cochain(x.m,lambda face:table[face])
        s = Cochain(1,lambda face:sum(row[0] for row in sb.rows[face[0]:face[1]])%2)
        return binary_cochain(b),binary_cochain(bp),s,binary_cochain(w)

    @lru_cache(None)
    def evaluate(self,simplex):
        return self.source(*self.inputs(simplex))(tuple(range(degree(simplex)+1)))

    def chain_value(self,chain):
        return sum(c*self.evaluate(x) for x,c in chain.items())

    def section(self,b,bp,s,omega,vertices):
        q = len(vertices)-1
        def binary_model(c):
            return KB(c.degree,q,tuple(c(tuple(vertices[i] for i in f))%2
                                      for f in faces(q,c.degree)))
        sb = cm.Diag('zsign',tuple((s((vertices[i],vertices[i+1])),(0,)*q)
                                  for i in range(q)))
        return ((binary_model(b),binary_model(bp)),(sb,binary_model(omega)))

    @lru_cache(None)
    def minimal_period(self):
        m = self.m
        if m != 3:
            return F(0)
        generator = KB(m,m,(1,))
        point_s = cm.Diag('zsign',())
        point_w = KB(2,0,())
        return self.chain_value(inclusion(((generator,generator),(point_s,point_w))))

    def primitive(self,b,bp,s,omega):
        return Cochain(self.k,lambda vertices:self.chain_value(
            homotopy(self.section(b,bp,s,omega,vertices))))

    def phase(self,b,c,bp,cp,s,omega):
        beta = alpha(b,bp,s)
        return F(1,2)*binary_sum(polarization(c,cp),
            polarization(binary_sum(c,cp),beta))+self.primitive(b,bp,s,omega)

    def gamma(self,b,c,bp,cp,s,omega):
        total_b = binary_sum(b,bp)
        total_c = binary_sum(c,cp,alpha(b,bp,s))
        return integral(self.omega(b,c,s,omega)+self.omega(bp,cp,s,omega)
                        -self.omega(total_b,total_c,s,omega)
                        +signed_differential(self.phase(b,c,bp,cp,s,omega),s),
                        'legal A=0 gamma'+str(self.k))
