"""Production legal gamma5 with the matched mixed beta correction.

The old beta5 is changed by rho(A) cup rho(Aprime). The residual universal
source has zero periods only after the two explicit B transgressions.
"""
from functools import lru_cache
from fractions import Fraction as F
import upper_pair_source as aux
import production_gamma5_comparison as comparison
from universal_pair6 import SourcePrimitive6
p = aux.p


def mixed(A, Ap):
    return p.binary(p.cup(p.binary(A), p.binary(Ap)))


def transgression_source(A, Ap, s, omega):
    return p.half(p.binary(p.cup(p.binary(A), aux.upper.primary(Ap,s,omega))
                          +p.cup(p.binary(Ap), aux.upper.primary(A,s,omega))))


def transgression(A, B, Ap, Bp):
    return p.half(p.binary(p.cup(p.binary(A), Bp)+p.cup(p.binary(Ap), B)))


def source(A, Ap, s, omega):
    # aux.source uses the current, matched beta5 supplied by beta_sharp.
    return (aux.source(A, Ap, s, omega)+comparison.A_phase(A,s,omega)
            +comparison.A_phase(Ap,s,omega)-comparison.A_phase(A+Ap,s,omega)
            +comparison.affine_B_source(A,Ap,s,omega)
            -transgression_source(A,Ap,s,omega))


FIXED_COEFFICIENTS = (F('7/8'),F('3/8'),F('1/4'),F('1/4'),F('1/2'),F('1/2'),F('1/2'),F('3/4'),F('1/2'),F('3/4'),F('3/4'),F('3/4'),F('1/2'),F('1/2'),F('1/4'),F('0'),F('1/4'),F('0'),F('0'),F('3/4'),F('3/4'),F('3/8'),F('0'),F('1/4'),F('1/4'),F('3/8'),F('0'),F('0'),F('0'),F('0'),F('1/4'),F('7/8'),F('0'),F('1/4'),F('3/4'),F('3/4'),F('3/4'),F('0'),F('1/8'),F('5/8'),F('0'),F('0'),F('0'),F('5/8'),F('0'),F('0'),F('0'),F('0'),)


class ProductionSourcePrimitive(SourcePrimitive6):
    def coefficients(self):
        # Matched source periods and pivot values are recorded in the certificate.
        return FIXED_COEFFICIENTS


@lru_cache(None)
def constructor():
    return ProductionSourcePrimitive(source)


def phase(A,B,C,Ap,Bp,Cp,s,omega):
    if (A.degree,B.degree,C.degree,Ap.degree,Bp.degree,Cp.degree)!=(2,3,4,2,3,4):
        raise ValueError('production gamma5 expects two degree-(2,3,4) triples')
    upper=aux.upper;I=upper.hp.interval
    alpha=upper.hp.hD(p.binary(A),p.binary(Ap),s)
    Bsum=p.binary(B+Bp+alpha)
    along=aux.reduced_full_source(I(A),I(B,True),I(C,True),I(Ap),I(Bp,True),I(Cp,True),I(s),I(omega))
    gauge=(comparison.comparison_gauge(A,B,s,omega)+comparison.comparison_gauge(Ap,Bp,s,omega)
           -comparison.comparison_gauge(A+Ap,Bsum,s,omega))
    return (p.scale(constructor().primitive(A,Ap,s,omega),-1)-upper.hp.prism(along)
            -comparison.affine_B_primitive(A,B,Ap,Bp,s,omega)-gauge
            -transgression(A,B,Ap,Bp))


def gamma(A,B,C,Ap,Bp,Cp,s,omega):
    from upper_phase_diagnostic import production_phase
    upper=aux.upper
    alpha=upper.hp.hD(p.binary(A),p.binary(Ap),s)
    beta=upper.beta_sharp(A,B,Ap,Bp,s,omega)
    Bs,Cs=p.binary(B+Bp+alpha),p.binary(C+Cp+beta)
    raw=(production_phase(2,A,B,C,s,omega)+production_phase(2,Ap,Bp,Cp,s,omega)
         -production_phase(2,A+Ap,Bs,Cs,s,omega)+p.ds(phase(A,B,C,Ap,Bp,Cp,s,omega),s))
    def value(vertices):
        result=F(raw(vertices))
        if result.denominator!=1:
            raise ArithmeticError(f'production legal gamma5 is not integral: {result}')
        return result.numerator
    return p.Cochain(6,value)
