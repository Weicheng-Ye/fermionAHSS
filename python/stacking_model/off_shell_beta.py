"""Prism transport of calibrated legal beta to arbitrary lower cochains.

The fixed comparison h with delta h=Tau+H supplies a local extension q.
The global input flags implement exactly the existing piecewise differential.
This module never changes a production operation or solves a cochain on X.
The successor legal beta is an explicit supplied universal formula.
"""
from dataclasses import dataclass
from itertools import combinations

import h_tau_primitive as hp
from h_tau_primitive import p, low


@dataclass(frozen=True)
class LowerPair:
    A: object
    B: object
    legal: bool

    def __post_init__(self):
        if self.A.degree < 0 or self.B.degree != self.A.degree+1:
            raise ValueError('expected A of degree n>=0 and B of degree n+1')


def primary(A,s,omega):
    return p.QD(p.binary(A),s,omega)


def curvature(A,B,s,omega):
    return (p.ds(A,s),p.binary(p.differential(B)+primary(A,s,omega)))


def alpha(A,Aprime,s):
    return hp.hD(p.binary(A),p.binary(Aprime),s)


def raw_H(A,B,s,omega):
    ai,bi=hp.interval(A,True),hp.interval(B,True)
    si,wi=hp.interval(s),hp.interval(omega)
    outa,outb=curvature(ai,bi,si,wi)
    return p.binary(hp.prism(low.secondary(outa,outb,si,wi)))


def local_comparison_extension(A,B,s,omega):
    """q equals h on locally legal simplex data, and zero otherwise."""
    n=A.degree
    if n==0:
        return p.zero(2)
    a,b=curvature(A,B,s,omega)
    h=hp.comparison_gauge(A,B,s,omega)
    def value(vertices):
        if any(a(face) for face in combinations(vertices,n+2)):
            return 0
        if any(b(face) for face in combinations(vertices,n+3)):
            return 0
        return h(vertices)
    return p.Cochain(n+2,value)


def natural_f(A,B,s,omega):
    """A natural auxiliary f, equal to the calibrated Tau on legal inputs."""
    return p.binary(raw_H(A,B,s,omega)
                    +p.differential(local_comparison_extension(A,B,s,omega)))


def coordinate_gauge(pair,s,omega):
    return (p.zero(pair.A.degree+2) if pair.legal else
            local_comparison_extension(pair.A,pair.B,s,omega))


def current_f(pair,s,omega):
    if pair.legal:
        return low.secondary(pair.A,pair.B,s,omega)
    return raw_H(pair.A,pair.B,s,omega)


def _closed_pair_residual(A,B,Aprime,Bprime,s,omega,next_legal_beta):
    correction=alpha(A,Aprime,s)
    Asum,Bsum=A+Aprime,p.binary(B+Bprime+correction)
    a,b=curvature(A,B,s,omega)
    ap,bp=curvature(Aprime,Bprime,s,omega)
    return p.binary(natural_f(Asum,Bsum,s,omega)+natural_f(A,B,s,omega)
                    +natural_f(Aprime,Bprime,s,omega)
                    +next_legal_beta(a,b,ap,bp,s,omega))


def raw_beta(A,B,Aprime,Bprime,s,omega,next_legal_beta):
    si,wi=hp.interval(s),hp.interval(omega)
    residual=_closed_pair_residual(
        hp.interval(A,True),hp.interval(B,True),
        hp.interval(Aprime,True),hp.interval(Bprime,True),
        si,wi,next_legal_beta)
    return p.binary(hp.prism(residual))


def beta(pair,pairprime,sum_legal,s,omega,legal_beta,next_legal_beta):
    """Exact beta for the existing piecewise f and explicit legal formulas.

    All flags mean global legality. legal_beta and next_legal_beta are
    the calibrated universal formulas at n and n+1, respectively. They
    must be normalized and obey delta beta=Delta Tau on their domain.
    For the off-shell branch this formula evaluates next_legal_beta only
    on curvature pairs, which are always legal.
    """
    A,B,Ap,Bp=pair.A,pair.B,pairprime.A,pairprime.B
    if A.degree!=Ap.degree:
        raise ValueError('both inputs must have the same degree')
    if pair.legal and pairprime.legal:
        if not sum_legal:
            raise ValueError('two legal inputs necessarily have a legal lower sum')
        return legal_beta(A,B,Ap,Bp,s,omega)
    total=LowerPair(A+Ap,p.binary(B+Bp+alpha(A,Ap,s)),sum_legal)
    return p.binary(raw_beta(A,B,Ap,Bp,s,omega,next_legal_beta)
                    +coordinate_gauge(pair,s,omega)
                    +coordinate_gauge(pairprime,s,omega)
                    +coordinate_gauge(total,s,omega))
