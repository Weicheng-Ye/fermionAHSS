"""Explicit production gamma6 on full legal defining systems.

This constructs the complete two- and three-primary correction. Its
universal source is reduced by the integral H7=0 certificate. It does not
claim the nonzero-A all-cochain extension or its exchange homotopy.
"""
from fractions import Fraction as F
from functools import lru_cache
import upper_pair_source as aux
import production_gamma6_comparison as comparison
import production_upper_binary_comparison as binary
from universal_pair7 import SourcePrimitive7
p = aux.p


def source(A, Ap, s, omega):
    return (aux.source(A, Ap, s, omega)+comparison.A_phase(A, s, omega)
            +comparison.A_phase(Ap, s, omega)-comparison.A_phase(A+Ap, s, omega)
            +p.half(binary.affine_J0_source(A, Ap, s, omega)))


FIXED_COEFFICIENTS = (F(19,24),F(7,24),F(11,12),F(5,12),F(0),F(11,12),
                      F(11,12),F(3,8),F(3,8),F(0),F(0),F(0))


class ProductionSourcePrimitive(SourcePrimitive7):
    def coefficients(self):
        # Ten actual source evaluations and the integral H7 certificate
        # determine these fixed coefficients; see the JSON certificate.
        return FIXED_COEFFICIENTS


@lru_cache(None)
def constructor():
    return ProductionSourcePrimitive(source)


def phase(A, B, C, Ap, Bp, Cp, s, omega):
    if (A.degree, B.degree, C.degree, Ap.degree, Bp.degree, Cp.degree) != (3,4,5,3,4,5):
        raise ValueError('production gamma6 expects two degree-(3,4,5) triples')
    upper = aux.upper
    interval = upper.hp.interval
    alpha = upper.hp.hD(p.binary(A), p.binary(Ap), s)
    Bsum = p.binary(B+Bp+alpha)
    along = aux.reduced_full_source(interval(A), interval(B, True), interval(C, True),
                            interval(Ap), interval(Bp, True), interval(Cp, True),
                            interval(s), interval(omega))
    gauge_difference = (comparison.comparison_gauge(A, B, s, omega)
                         +comparison.comparison_gauge(Ap, Bp, s, omega)
                         -comparison.comparison_gauge(A+Ap, Bsum, s, omega))
    return (p.scale(constructor().primitive(A, Ap, s, omega), -1)
            -upper.hp.prism(along)
            -p.half(binary.affine_J0_primitive(A, B, Ap, Bp, s, omega))
            -gauge_difference)


def gamma(A, B, C, Ap, Bp, Cp, s, omega):
    from upper_phase_diagnostic import production_phase
    upper = aux.upper
    alpha = upper.hp.hD(p.binary(A), p.binary(Ap), s)
    beta = upper.beta_sharp(A, B, Ap, Bp, s, omega)
    total_b, total_c = p.binary(B+Bp+alpha), p.binary(C+Cp+beta)
    raw = (production_phase(3, A, B, C, s, omega)
           +production_phase(3, Ap, Bp, Cp, s, omega)
           -production_phase(3, A+Ap, total_b, total_c, s, omega)
           +p.ds(phase(A, B, C, Ap, Bp, Cp, s, omega), s))
    def value(vertices):
        result = F(raw(vertices))
        if result.denominator != 1:
            raise ArithmeticError(f'production legal gamma6 is not integral: {result}')
        return result.numerator
    return p.Cochain(7, value)
