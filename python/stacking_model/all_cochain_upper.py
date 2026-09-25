"""Branchwise transport to the original all-cochain g in k3, k4 and k5.

All successor corrections use the explicitly rephased legal product.
Input flags are global cochain conditions, as in stacking_lower.py.
"""
from dataclasses import dataclass
from itertools import combinations
import natural_upper as natural
import nonzero_offshell_comparison as comparison
import closed_ab_upper as closed
import stacking_lower as lower
import off_shell_beta as off
import pure_c_normalization as normalized
from fractions import Fraction as F
p,hp=off.p,off.hp


@dataclass(frozen=True)
class Triple:
    A: object
    B: object
    C: object
    lower_legal: bool
    pure_closed: bool = False
    full_legal: bool | None = None

    def __post_init__(self):
        if (self.B.degree,self.C.degree)!=(self.A.degree+1,self.A.degree+2):
            raise ValueError('inconsistent triple degrees')
        if self.pure_closed:
            object.__setattr__(self,'full_legal',True)
        if self.full_legal and not self.lower_legal:
            raise ValueError('a fully legal triple has legal lower pair')


def on_simplex(A,B,C,s,omega,vertices):
    def vanishes(c):return not any(c(f) for f in combinations(vertices,c.degree+1))
    legal=all(vanishes(c)for c in off.curvature(A,B,s,omega))
    pure=vanishes(A)and vanishes(B)and vanishes(p.binary(p.differential(C)))
    full=legal and vanishes(p.binary(p.differential(C)+off.current_f(
                off.LowerPair(A,B,True),s,omega)))
    return Triple(A,B,C,legal,pure,full)


def r(u,s):
    return (p.half(p.cup(s,u.C))if u.A.degree>=1 and u.pure_closed
            else p.zero(u.C.degree+1))


def K(u,s):
    return natural.integral(p.ds(r(u,s),s),'integral successor rephasing')


def current_curvature(u,s,omega):
    a,b=off.curvature(u.A,u.B,s,omega)
    pair=off.LowerPair(u.A,u.B,u.lower_legal)
    return a,b,p.binary(p.differential(u.C)+off.current_f(pair,s,omega))


def product(u,v,sum_legal,sum_pure_closed,s,omega):
    alpha=lower.alpha(u.A,v.A,s)
    beta=lower.all_cochain_beta(off.LowerPair(u.A,u.B,u.lower_legal),
             off.LowerPair(v.A,v.B,v.lower_legal),sum_legal,s,omega)
    return Triple(u.A+v.A,p.binary(u.B+v.B+alpha),p.binary(u.C+v.C+beta),sum_legal,sum_pure_closed,
                  True if u.full_legal and v.full_legal else None)


def coordinate(u,s,omega):
    e=off.coordinate_gauge(off.LowerPair(u.A,u.B,u.lower_legal),s,omega)
    return u.A,u.B,p.binary(u.C+e)


def integer_gauge(u,s,omega):
    if not u.lower_legal:
        if u.A.degree==0:
            return p.zero(4)
        return comparison.off_branch_integer_gauge(u.A,u.B,u.C,s,omega)
    if u.A.degree==0:
        import nonzero_degree3_comparison as degree3
        q=degree3.gauge(u.A,u.B,u.C,s,omega)
        raw=degree3.difference(u.A,u.B,u.C,s,omega)-p.ds(q,s)
        return natural.integral(raw,'degree-three closed-branch upper comparison')
    elif u.A.degree==1:
        q=comparison.legal_pair_gauge4(u.A,u.B,u.C,s,omega)
    elif u.A.degree==2:
        q=comparison.legal_pair_gauge5(u.A,u.B,u.C,s,omega)
    else:
        raise ValueError('the transported upper operation covers k3,k4,k5')
    t=current_curvature(u,s,omega)[2]
    raw=(comparison.legal_pair_potential(u.A,u.B,u.C,s,omega)
         -comparison.raw_potential(u.A,u.B,u.C,s,omega)-p.half(p.cup(s,t))-p.ds(q,s))
    if u.A.degree==2:
        cube=hp.low.transported(hp.low.transported(u.A,u.A,s),u.A,s)
        raw=raw+p.scale(cube,F(8,3))
    return natural.integral(raw,'closed-branch upper comparison')


def legal_gamma(u,v,total,s,omega):
    """Rephased prescribed legal product, also used as the successor."""
    if u.A.degree==0:
        import closed_ab_degree3 as degree3
        answer=degree3.gamma(u.A,u.B,u.C,v.A,v.B,v.C,s,omega)
        import degree3_legal_recalibration as recalibration
        correction=recalibration.gamma_change(u.A,u.B,u.C,v.A,v.B,v.C,s,omega)
        return p.Cochain(answer.degree,answer)+p.Cochain(correction.degree,correction)
    return normalized.gamma(u.A,u.B,u.C,v.A,v.B,v.C,s,omega)-K(u,s)-K(v,s)+K(total,s)


def g(u,s,omega):
    """The prescribed J/G upper differential, without a change of D."""
    if u.A.degree not in (0,1,2):
        raise ValueError('this direct differential covers k3,k4,k5')
    if u.lower_legal:
        if u.A.degree==0:
            import nonzero_degree3_comparison as degree3
            t=current_curvature(u,s,omega)[2]
            return natural.integral(p.ds(degree3.potential(u.A,u.B,u.C,s,omega),s)
                                    -p.half(p.E(t,omega)),'degree-three J')
        return closed.J(u.A,u.B,u.C,s,omega)
    from upper_phase_diagnostic import production_phase
    I=hp.interval
    A,B,C,ss,ww=I(u.A,True),I(u.B,True),I(u.C,True),I(s),I(omega)
    a,b=off.curvature(A,B,ss,ww)
    c=p.binary(p.differential(C)+off.raw_H(A,B,ss,ww))
    T=p.ds(production_phase(u.A.degree+1,a,b,c,ss,ww),ss)
    return natural.integral(p.scale(hp.prism(T),-1),'current upper differential')


def gamma(u,v,sum_legal,sum_pure_closed,s,omega):
    if u.A.degree not in(0,1,2):raise ValueError('all-cochain upper transport covers k3,k4,k5')
    total=product(u,v,sum_legal,sum_pure_closed,s,omega)
    change=p.zero(u.A.degree+4)-K(u,s)-K(v,s)+K(total,s)
    if u.lower_legal and v.lower_legal:
        t,tp,ts=[current_curvature(w,s,omega)[2]for w in(u,v,total)]
        # All three images are closed pure-C triples. The entire half-lift
        # difference is an integer carry, before taking any coboundary.
        carry=natural.integral(p.half(p.cup(s,t))+p.half(p.cup(s,tp))-p.half(p.cup(s,ts)),
                               'pure-C successor carry')
        if u.A.degree==0:
            import closed_ab_degree3 as degree3
            answer=degree3.gamma(u.A,u.B,u.C,v.A,v.B,v.C,s,omega)
            base=p.Cochain(answer.degree,answer)
            if u.full_legal is None or v.full_legal is None:
                raise ValueError('k3 needs global full_legal flags; use on_simplex or Stacking.triple')
            if u.full_legal and v.full_legal:
                import degree3_legal_recalibration as recalibration
                correction=recalibration.gamma_change(u.A,u.B,u.C,v.A,v.B,v.C,s,omega)
                base=base+p.Cochain(correction.degree,correction)
        else:
            base=normalized.closed_gamma(u.A,u.B,u.C,v.A,v.B,v.C,s,omega)
        return base-carry+change
    au,bu,cu=coordinate(u,s,omega);av,bv,cv=coordinate(v,s,omega)
    # On this branch, the transformed current product equals the raw
    # natural product exactly; no legal-beta comparison is needed.
    return (natural.gamma(au,bu,cu,av,bv,cv,s,omega)+integer_gauge(u,s,omega)
            +integer_gauge(v,s,omega)-integer_gauge(total,s,omega)+change)
