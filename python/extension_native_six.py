"""Degree-six D-layer terms of the native transferred stacking model.

Package degree six is the cutoff of the all-cochain differential
(doc/all_cochain_differential.md, Section 3). Its lower three components
are the same formulas as in degrees three to five: delta_s A, delta B+P_6(A)
and delta C+f_6(A,B), with f_6 the calibrated Tau on the legal lower locus
L_6 and the prism H_6 outside it. The D-layer terms differ from degrees
three to five only outside L_6:

- g_6 is J_6 on L_6, the integral cochain of (A5) evaluated by the fixed
  degree-three phase Omega_6 with the residual correction, for arbitrary C;
- gamma_6 on two lower-legal triples is the legal production stacking
  correction (closed_ab_upper.gamma with production_gamma6.phase), with the
  common pure-C normalization and the K rephasing of the successor, exactly
  as the reference finite-section model evaluates it on legal data;
- outside L_6 the note uses the section retraction R_6 of a finite free
  cochain model. The native model never evaluates that branch: native
  curvature stops at the first nonzero layer, so g_6 is requested only on
  lower-legal triples, and products and gauge actions receive flat states
  and differential images, whose lower pairs are legal. A request outside
  L_6 is reported as unresolved instead of being answered by a retraction
  that the note does not define on the resolution.

The stacking correction of two states whose A layers are not both zero
evaluates the universal pair source of production_gamma6. Each of its
universal terms contains three V3 source contractions on eight-vertex
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
import pure_c_normalization as normalized

p = upper.p


A_STACKING = os.environ.get('FERMIONAHSS_DEGREE_SIX_A_STACKING', '0').strip() == '1'


class NativeDegreeSixLimit(RuntimeError):
    """A degree-six request the native model leaves unresolved."""


class SectionBranchRequired(NativeDegreeSixLimit):
    """A degree-six D-layer value outside the legal lower locus was requested."""


class PairSourceLimit(NativeDegreeSixLimit):
    """The stacking correction of degree-six states with a nonzero A layer."""


def _check(triple, name):
    if not isinstance(triple, upper.Triple) or triple.A.degree != 3:
        raise ValueError(f'the native degree-six {name} expects triples of degrees (3,4,5)')


def g(triple, s, omega):
    """g_6(A,B,C) of (A13) on the legal lower locus: the J_6 branch.

    C is arbitrary; on a full defining system this is T_3 = delta_s Omega_6.
    """
    _check(triple, 'differential')
    if not triple.lower_legal:
        raise SectionBranchRequired(
            'the degree-six D-layer differential outside the legal lower locus '
            'needs the section retraction of the finite cochain model, which the '
            'native resolution model does not evaluate')
    return closed.J(triple.A, triple.B, triple.C, s, omega)


def gamma(u, v, sum_legal, sum_pure_closed, a_zero, s, omega):
    """gamma_6 of two lower-legal triples, with the pure-C successor rephasing.

    This is the legal branch of all_cochain_upper.gamma in package degree
    six: the rephased production correction closed_gamma, minus the
    half-lift carry of the three C residuals, plus the K change of the
    pure-C rephasing. The C residuals may be nonzero. ``a_zero`` records
    whether both A layers vanish on the resolution.
    """
    _check(u, 'stacking correction')
    _check(v, 'stacking correction')
    if not (u.lower_legal and v.lower_legal):
        raise SectionBranchRequired(
            'the degree-six stacking correction outside the legal lower locus '
            'needs the section retraction of the finite cochain model, which the '
            'native resolution model does not evaluate')
    if not a_zero and not A_STACKING:
        raise PairSourceLimit(
            'the degree-six stacking correction of states with a nonzero A layer '
            'evaluates the universal pair source of production_gamma6, whose terms '
            'nest the V3 contraction and which is not practical to evaluate; set '
            'FERMIONAHSS_DEGREE_SIX_A_STACKING=1 to evaluate it regardless of its '
            'running time')
    total = upper.product(u, v, sum_legal, sum_pure_closed, s, omega)
    change = p.zero(7) - upper.K(u, s) - upper.K(v, s) + upper.K(total, s)
    t, tp, ts = [upper.current_curvature(w, s, omega)[2] for w in (u, v, total)]
    # All three images are closed pure-C triples; the entire half-lift
    # difference is an integer carry before taking any coboundary.
    carry = natural.integral(p.half(p.cup(s, t)) + p.half(p.cup(s, tp))
                             - p.half(p.cup(s, ts)), 'pure-C successor carry')
    base = normalized.closed_gamma(u.A, u.B, u.C, v.A, v.B, v.C, s, omega)
    return base - carry + change
