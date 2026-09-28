"""D-layer terms of the native transferred stacking model in degrees 3 to 6.

The D layer of a state in package degree k is corrected by g_k in the
differential and by gamma_k in the product (doc/all_cochain_differential.md,
Section 3). Degrees 3 to 5 use the prescribed all-cochain transport of
all_cochain_upper on every branch; degree 6, the cutoff, uses J_6 and the
legal gamma_6 on the legal lower locus. Two exact simplifications apply to
the legal branch of every degree from 4 to 6:

- The production correction is closed_ab_upper.gamma plus the pure-C
  normalization delta_s(phase_correction) of pure_c_normalization, whose
  explicit integral primitive integer_phase_correction shows it to be an
  integral coboundary. A D-gauge with that primitive carries one product to
  the other, so the stacking group on gauge classes is the same, while the
  second evaluation of the phase at A=B=0 is saved. The native product uses
  closed_ab_upper.gamma directly, with the half-lift carry of the three C
  residuals and the K change of the pure-C rephasing, which are not
  coboundaries in general and are kept. Native D representatives therefore
  differ from the finite-section reference model's by D-coboundaries.
- Outside the legal lower locus in degree 6 the note uses the section
  retraction R_6 of a finite free cochain model. The native model never
  evaluates that branch: native curvature stops at the first nonzero layer,
  so g_6 is requested only on lower-legal triples, and products and gauge
  actions receive flat states and differential images, whose lower pairs
  are legal. Such a request is reported as unresolved.

The stacking correction of two degree-six states whose A layers are not
both zero evaluates the universal pair source of production_gamma6. Each of
its universal terms contains three V3 source contractions on eight-vertex
universal simplices, and the pair contraction itself has thousands of terms
per output simplex, so this correction costs orders of magnitude more than
any degree-five term. It is therefore refused, leaving the degree
unresolved, unless FERMIONAHSS_DEGREE_SIX_A_STACKING=1 is set; the
differential J_6 with nonzero A and every term with A=0 are evaluated.

Copyright (c) 2026 koAHSS contributors; MIT license.
"""
import os

import all_cochain_upper as upper
import closed_ab_upper as closed
import natural_upper as natural

p = upper.p


A_STACKING = os.environ.get('FERMIONAHSS_DEGREE_SIX_A_STACKING', '0').strip() == '1'


class NativeDegreeSixLimit(RuntimeError):
    """A degree-six request the native model leaves unresolved."""


class SectionBranchRequired(NativeDegreeSixLimit):
    """A degree-six D-layer value outside the legal lower locus was requested."""


class PairSourceLimit(NativeDegreeSixLimit):
    """The stacking correction of degree-six states with a nonzero A layer."""


def _check(k, triple, name):
    if not isinstance(triple, upper.Triple) or triple.A.degree != k - 3:
        raise ValueError(f'the native degree-{k} {name} expects triples of degrees '
                         f'({k - 3},{k - 2},{k - 1})')


def g(k, triple, s, omega):
    """g_k(A,B,C) of (A13) in degree k: J or the prism G below the cutoff, J_6 at it."""
    if k in (3, 4, 5):
        return upper.g(triple, s, omega)
    if k != 6:
        raise ValueError('the native D-layer differential covers degrees 3 to 6')
    _check(k, triple, 'differential')
    if not triple.lower_legal:
        raise SectionBranchRequired(
            'the degree-six D-layer differential outside the legal lower locus '
            'needs the section retraction of the finite cochain model, which the '
            'native resolution model does not evaluate')
    return closed.J(triple.A, triple.B, triple.C, s, omega)


def legal_gamma(u, v, sum_legal, sum_pure_closed, s, omega):
    """gamma_k of two lower-legal triples in degrees 4 to 6, without the pure-C normalization.

    This is the legal branch of all_cochain_upper.gamma with
    closed_ab_upper.gamma in place of the normalized correction; the C
    residuals of the two inputs may be nonzero.
    """
    total = upper.product(u, v, sum_legal, sum_pure_closed, s, omega)
    change = p.zero(u.A.degree + 4) - upper.K(u, s) - upper.K(v, s) + upper.K(total, s)
    t, tp, ts = [upper.current_curvature(w, s, omega)[2] for w in (u, v, total)]
    # All three images are closed pure-C triples; the entire half-lift
    # difference is an integer carry before taking any coboundary.
    carry = natural.integral(p.half(p.cup(s, t)) + p.half(p.cup(s, tp))
                             - p.half(p.cup(s, ts)), 'pure-C successor carry')
    base = closed.gamma(u.A, u.B, u.C, v.A, v.B, v.C, s, omega)
    return base - carry + change


def gamma(k, u, v, sum_legal, sum_pure_closed, a_zero, s, omega):
    """gamma_k of two triples in degree k; ``a_zero`` says whether both A layers vanish on R."""
    if k == 3:
        return upper.gamma(u, v, sum_legal, sum_pure_closed, s, omega)
    if k not in (4, 5, 6):
        raise ValueError('the native stacking correction covers degrees 3 to 6')
    _check(k, u, 'stacking correction')
    _check(k, v, 'stacking correction')
    if not (u.lower_legal and v.lower_legal):
        if k == 6:
            raise SectionBranchRequired(
                'the degree-six stacking correction outside the legal lower locus '
                'needs the section retraction of the finite cochain model, which the '
                'native resolution model does not evaluate')
        return upper.gamma(u, v, sum_legal, sum_pure_closed, s, omega)
    if k == 6 and not a_zero and not A_STACKING:
        raise PairSourceLimit(
            'the degree-six stacking correction of states with a nonzero A layer '
            'evaluates the universal pair source of production_gamma6, whose terms '
            'nest the V3 contraction and which is not practical to evaluate; set '
            'FERMIONAHSS_DEGREE_SIX_A_STACKING=1 to evaluate it regardless of its '
            'running time')
    return legal_gamma(u, v, sum_legal, sum_pure_closed, s, omega)
