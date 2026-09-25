"""One degree-indexed API for the exact four-cochain differential/product.

The zero predicate tests an entire cochain on X.  Degree six additionally
uses the explicit finite section adapter; its only injected maps are
linear basis conversions.  Degree seven is an endpoint convention, and
its product is needed only on degree-six differential images.
"""
from dataclasses import dataclass

import off_shell_beta as off
import stacking_lower as lower
import all_cochain_upper as upper
from cochains import Cochain as LocalCochain
from coherent_low_commutative import (DegreeOneCommutativeStacking,
                                     DegreeTwoCommutativeStacking)
p=off.p


def cochain(c):
    return p.Cochain(c.degree,c)


def local(c):
    return LocalCochain(c.degree,c)


@dataclass(frozen=True)
class State:
    A: object
    B: object
    C: object
    D: object

    def __post_init__(self):
        k=self.k
        if (self.A.degree,self.B.degree,self.C.degree)!=(k-3,k-2,k-1):
            raise ValueError('expected degrees (k-3,k-2,k-1,k+1)')

    @property
    def k(self):return self.D.degree-1


class Stacking:
    def __init__(self,s,omega,is_zero,*,degree_six=None):
        if (s.degree,omega.degree)!=(1,2):
            raise ValueError('the twists have degrees one and two')
        self.s,self.omega=cochain(s),cochain(omega)
        self.is_zero,self.degree_six=is_zero,degree_six
        if not is_zero(p.binary(p.differential(self.s))):
            raise ValueError('s must be a binary cocycle')
        if not is_zero(p.binary(p.differential(self.omega))):
            raise ValueError('omega must be a binary cocycle')
        if degree_six is not None and any(not is_zero(p.binary(left-cochain(right)))
                for left,right in ((self.s,degree_six.s),(self.omega,degree_six.omega))):
            raise ValueError('the degree-six adapter must use these same twists')
        self.low2=DegreeTwoCommutativeStacking()
        self.low1=DegreeOneCommutativeStacking(successor=self.low2)

    def legal_pair(self,A,B):
        return all(self.is_zero(c)for c in off.curvature(A,B,self.s,self.omega))

    def triple(self,A,B,C):
        legal=self.legal_pair(A,B)
        pure=(self.is_zero(A)and self.is_zero(B)
              and self.is_zero(p.binary(p.differential(C))))
        full=legal and self.is_zero(p.binary(p.differential(C)+off.current_f(
                    off.LowerPair(A,B,True),self.s,self.omega)))
        return upper.Triple(A,B,C,legal,pure,full)

    def finite(self,x):
        if self.degree_six is None:
            raise ValueError('provide the fixed ProductionG6Section adapter for k6/k7')
        model=self.degree_six
        return tuple(model.vector(x.A)),model.bits(x.B),model.bits(x.C)

    @staticmethod
    def converted(x):
        return State(*(cochain(c)for c in(x.A,x.B,x.C,x.D)))

    def d(self,x):
        x=self.converted(x);k=x.k;s,w=self.s,self.omega
        if k not in range(-1,8):
            raise ValueError('d is provided in k=-1..6 and at the k7 endpoint')
        if k<=0:
            return State(p.zero(k-2),p.zero(k-1),p.zero(k),p.ds(x.D,s))
        if k==1:
            g=cochain(self.low1.g(*map(local,(x.C,s,w))))
            return State(p.zero(-1),p.zero(0),p.binary(p.differential(x.C)),p.ds(x.D,s)+g)
        if k==2:
            from a0_gamma import lower_d2
            b,c=lower_d2(*map(local,(x.B,x.C,s,w)))
            g=self.low2.g(*map(local,(x.B,x.C,s,w)),
                             b_closed=self.is_zero(p.binary(p.differential(x.B))))
            return State(p.zero(0),cochain(b),cochain(c),p.ds(x.D,s)+cochain(g))
        u=self.triple(x.A,x.B,x.C)
        a,b,c=upper.current_curvature(u,s,w)
        if k<=5:
            g=upper.g(u,s,w)
        elif k==6:
            value=self.finite(x)
            g=self.degree_six.integral_cochain(8,self.degree_six.correction(value))
        elif self.is_zero(x.A)and self.is_zero(x.B):
            from compatible_sector import pure_c_g
            g=pure_c_g(*map(local,(x.C,s,w)))
        else:
            g=p.zero(9)
        return State(a,b,c,p.ds(x.D,s)+cochain(g))

    def xtimes(self,x,y):
        x,y=self.converted(x),self.converted(y)
        if x.k!=y.k:
            raise ValueError('stacking inputs must have the same physical degree')
        k=x.k;s,w=self.s,self.omega
        if k not in range(-1,8):
            raise ValueError('stacking is provided through k6, with k7 on differential images')
        if k<=0:
            return State(x.A+y.A,p.binary(x.B+y.B),p.binary(x.C+y.C),x.D+y.D)
        if k==1:
            gamma=self.low1.gamma(*map(local,(x.C,y.C,s,w)))
            return State(p.zero(-2),p.zero(-1),p.binary(x.C+y.C),x.D+y.D+cochain(gamma))
        if k==2:
            from a0_gamma import beta2
            B=p.binary(x.B+y.B)
            beta=cochain(beta2(*map(local,(x.B,y.B,s))))
            gamma=self.low2.gamma(*map(local,(x.B,x.C,y.B,y.C,s,w)),
                b_closed=self.is_zero(p.binary(p.differential(x.B))),
                bp_closed=self.is_zero(p.binary(p.differential(y.B))),
                sum_closed=self.is_zero(p.binary(p.differential(B))))
            return State(p.zero(-1),B,p.binary(x.C+y.C+beta),x.D+y.D+cochain(gamma))
        A=x.A+y.A;B=p.binary(x.B+y.B+lower.alpha(x.A,y.A,s))
        if k==7:
            # This endpoint product is required only on legal images of d6.
            if not(self.legal_pair(x.A,x.B)and self.legal_pair(y.A,y.B)):
                raise ValueError('the k7 endpoint product is defined on d6 images')
            beta=lower.legal_beta(x.A,x.B,y.A,y.B,s,w)
            left,right=self.finite(x),self.finite(y)
            gamma=self.degree_six.integral_cochain(8,
                      self.degree_six.successor_gamma(left,right))
        else:
            u,v=self.triple(x.A,x.B,x.C),self.triple(y.A,y.B,y.C)
            total_legal=self.legal_pair(A,B)
            beta=lower.all_cochain_beta(off.LowerPair(x.A,x.B,u.lower_legal),
                        off.LowerPair(y.A,y.B,v.lower_legal),total_legal,s,w)
            C=p.binary(x.C+y.C+beta)
            total_pure=(self.is_zero(A)and self.is_zero(B)
                        and self.is_zero(p.binary(p.differential(C))))
            if k<=5:
                gamma=upper.gamma(u,v,total_legal,total_pure,s,w)
            else:
                left,right=self.finite(x),self.finite(y)
                gamma=self.degree_six.integral_cochain(7,self.degree_six.gamma(left,right))
        return State(A,B,p.binary(x.C+y.C+beta),x.D+y.D+cochain(gamma))
