"""Higher one-fiber production/Theta source comparison via binary bars.

The source is evaluated with the actual production rational lift. The fixed
small-complex coefficients are independently certified on every source cell;
a failed certificate raises an error rather than assuming the source exact.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations

from cochains import Cochain,zero
from a0_high_gamma import HigherA0Stacking,p
import binary_iterated_bar as bar

TREE=('b',('s','omega'))


class HigherThetaComparison:
    def __init__(self,m):
        if m not in (3,4):
            raise ValueError('this comparison covers m=3,4')
        self.m=m
        self.rule=HigherA0Stacking(m+2)

    def source(self,b,s,omega):
        theta=p.theta(b,zero(self.m+1),s,omega)
        return self.rule.boundary_phase(b,s,omega)-Cochain(theta.degree,theta)

    def inputs(self,simplex):
        b,(s,w)=simplex
        return (Cochain(self.m,lambda f:bar.cochain_value(b,f)),
                Cochain(1,lambda f:sum(row[0] for row in s.rows[f[0]:f[1]])%2),
                Cochain(2,lambda f:bar.cochain_value(w,f)))

    @lru_cache(None)
    def source_value(self,simplex):
        if bar.degree(simplex)!=self.m+3:
            raise ValueError('wrong comparison source degree')
        vertices=tuple(range(self.m+4))
        b,s,w=self.inputs(simplex)
        bv=tuple(b(face) for face in combinations(vertices,self.m+1))
        # The source is literally zero on the zero-B axis, including all
        # contractible diagonal-bar simplices with zero rectangular cocycle.
        if not any(bv):
            return F(0)
        sv=tuple(s((i,i+1)) for i in range(self.m+3))
        wv=tuple(w(face) for face in combinations(vertices,3))
        return self.source_coordinates(bv,sv,wv)

    @lru_cache(None)
    def source_coordinates(self,bvalues,sedges,wvalues):
        """Cache the actual normalized cochains, independently of bar extras."""
        vertices=tuple(range(self.m+4))
        bt=dict(zip(combinations(vertices,self.m+1),bvalues))
        wt=dict(zip(combinations(vertices,3),wvalues))
        b=Cochain(self.m,lambda face:bt[face])
        s=Cochain(1,lambda face:sum(sedges[face[0]:face[1]])%2)
        w=Cochain(2,lambda face:wt[face])
        return self.source(b,s,w)(vertices)

    def chain_value(self,chain):
        return sum(c*self.source_value(x) for x,c in chain.items())

    def tensor(self,pair,left,right):
        out={}
        for a,c in left(pair[0]).items():
            for b,d in right(pair[1]).items():
                bar.add(out,{(a,b):c*d})
        return out

    @lru_cache(None)
    def forward(self,simplex,tree=TREE):
        if isinstance(tree,str):
            return ({} if bar.degens(simplex) else {bar.degree(simplex):1}) if tree=='s' else bar.F(simplex)
        return bar.linear(bar.AW(simplex),lambda pair:self.tensor(pair,
            lambda x:self.forward(x,tree[0]),lambda x:self.forward(x,tree[1])))

    @lru_cache(None)
    def inclusion(self,key,tree=TREE):
        if isinstance(tree,str):
            if tree=='s':
                return {bar.cm.Diag('zsign',tuple((1,(0,)*key) for _ in range(key))):1}
            return bar.G(self.m if tree=='b' else 2,key)
        return bar.linear(self.tensor(key,lambda b:self.inclusion(b,tree[0]),
                                      lambda b:self.inclusion(b,tree[1])),bar.shuffle)

    @lru_cache(None)
    def homotopy(self,simplex,tree=TREE):
        if isinstance(tree,str):
            return {} if tree=='s' else bar.H(simplex)
        out=dict(bar.product_H(simplex))
        for (a,b),c in bar.AW(simplex).items():
            tensor={(aa,b):d for aa,d in self.homotopy(a,tree[0]).items()}
            projected=bar.linear(self.forward(a,tree[0]),lambda q:self.inclusion(q,tree[0]))
            for aa,d in projected.items():
                for bb,e in self.homotopy(b,tree[1]).items():
                    bar.add(tensor,{(aa,bb):bar.sign(bar.degree(a))*d*e})
            bar.add(out,bar.linear(tensor,bar.shuffle),c)
        return out

    def small_basis(self,n):
        return tuple((b,(ds,w)) for db in range(self.m,n+1)
                     for b in bar.small_basis(self.m,db)
                     for ds in range(n-db+1)
                     for w in bar.small_basis(2,n-db-ds))

    @lru_cache(None)
    def small_boundary(self,key):
        b,(ds,w)=key
        out={(a,(ds,w)):c for a,c in bar.small_boundary(self.m,b).items()}
        db=bar.small_degree(self.m,b)
        if ds%2:
            bar.add(out,{(b,(ds-1,w)):-2*bar.sign(db)})
        for a,c in bar.small_boundary(2,w).items():
            bar.add(out,{(b,(ds,a)):bar.sign(db+ds)*c})
        return out

    @lru_cache(None)
    def certificate(self):
        lower,upper=self.small_basis(self.m+2),self.small_basis(self.m+3)
        rows=[];values=[]
        for key in upper:
            boundary=self.small_boundary(key)
            rows.append(tuple(boundary.get(b,0) for b in lower))
            values.append(self.chain_value(self.inclusion(key)))
        if len(lower)!=5 or len(upper)!=11:
            raise ArithmeticError('the prescribed binary comparison basis changed')
        # Complete universal certificates fix these coefficients once and
        # for all. No coefficients are fitted at runtime or on X.
        solution=(F(0),)*5 if self.m==3 else (F(1,4),F(0),F(0),F(1,4),F(0))
        if any((sum(a*b for a,b in zip(row,solution))-value)%1
               for row,value in zip(rows,values)):
            raise ArithmeticError('fixed comparison certificate failed its residual check')
        return dict(lower=lower,upper=upper,rows=tuple(rows),values=tuple(values),
                    coefficients=solution,denominator=4)

    def section(self,b,s,omega,vertices):
        n=len(vertices)-1
        sb=bar.cm.Diag('zsign',tuple((s((vertices[i],vertices[i+1])),(0,)*n)
                                    for i in range(n)))
        return bar.coordinate(self.m,b,vertices),(sb,bar.coordinate(2,omega,vertices))

    def primitive(self,b,s,omega):
        certificate=self.certificate()
        small=dict(zip(certificate['lower'],certificate['coefficients']))
        def value(vertices):
            simplex=self.section(b,s,omega,vertices)
            return self.chain_value(self.homotopy(simplex))+sum(
                c*small.get(key,F(0)) for key,c in self.forward(simplex).items())
        return Cochain(self.m+2,value)
