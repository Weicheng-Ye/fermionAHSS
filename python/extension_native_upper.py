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

Two exact identities dispatch the A=0 sector of degrees 5 and 6 to the
direct formulas of a0_high_gamma and compatible_sector (doc/extensions.md):

- For a lower-legal triple with A=0 on R, the D-layer curvature
  g_k(0,B,C) = closed_ab_upper.J(0,B,C) equals, as a cochain,
  delta_s[Omega_a0(B,C) + h(t cup_{k-1} Q_D(B))] - h(E(t)) with
  Omega_a0 = HigherA0Stacking(k).omega and t the C residual, and equals
  compatible_sector.pure_c_g(C) when B=0 as well: production_phase(n,0,B,C)
  is Omega_a0(B,C) for closed B because the source splitting of A=0 vanishes,
  f^sharp(0,B) = Q_D(B), and every comparison gauge and universal primitive
  vanishes at A=0. No universal simplex is registered on this branch.
- For two fully legal triples with A=B=0 on R, the product correction
  closed_ab_upper.gamma equals compatible_sector.pure_c_gamma(C,C') plus
  delta_s of the explicit integral cochain
  -I(pure_c_gamma(Cl,C'l)) + I(lift(pol(delta(Cl),delta(C'l)))) - [s cup C][s cup C']
  (I the right prism, l the interval coordinate), an integral coboundary
  absorbed by a D-gauge; the half-lift carry and the K change are kept.
  With a nonzero B the difference of the two corrections is delta_s of a
  rational cochain whose integrality is only tested, so that sector keeps
  closed_ab_upper.gamma.

The stacking correction of two degree-six states with two nonzero A layers
evaluates the universal pair source of production_gamma6. Each of its
universal terms contains three V3 source contractions on eight-vertex
universal simplices, and the pair contraction itself has thousands of terms
per output simplex, so this correction costs orders of magnitude more than
any degree-five term. It is therefore refused, leaving the degree
unresolved, unless FERMIONAHSS_DEGREE_SIX_A_STACKING=1 is set; the
differential J_6 with nonzero A, every term with A=0 and the correction of
two states of which at most one has a nonzero A layer are evaluated: with
one fiber zero the pair primitive vanishes identically (every contraction
term has a zero fiber and the relative small basis has no such word), which
universal_pair7 short-circuits.

Copyright (c) 2026 koAHSS contributors; MIT license.
"""
import os

import all_cochain_upper as upper
import closed_ab_upper as closed
import compatible_sector as sector
import natural_upper as natural
from a0_high_gamma import HigherA0Stacking
from cochains import Cochain as LocalCochain

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


_A0 = {}


def _a0(k):
    """The one A=0 stacking instance per degree; it memoizes its universal values."""
    if k not in _A0:
        _A0[k] = HigherA0Stacking(k)
    return _A0[k]


def _local(cochain):
    """A phase_eval cochain as the cochain type of the stacking formulas."""
    return LocalCochain(cochain.degree, cochain)


def a0_curvature(k, triple, s, omega, b_zero):
    """g_k(0,B,C) of a lower-legal triple with A=0 on R, k=5,6, as an exact cochain.

    Equal to closed_ab_upper.J(0,B,C): pure_c_g(C) when B=0 as well, otherwise
    delta_s[Omega_a0(B,C) + h(t cup_{k-1} Q_D(B))] - h(E(t)).
    """
    if b_zero:
        return p.Cochain(k + 2, sector.pure_c_g(_local(triple.C), _local(s), _local(omega)))
    tau = p.QD(triple.B, s, omega)
    t = p.binary(p.differential(triple.C) + tau)
    omega_a0 = _a0(k).omega(_local(triple.B), _local(triple.C), _local(s), _local(omega))
    potential = p.Cochain(k + 1, omega_a0) + p.half(p.binary(p.cup(t, tau, k - 1)))
    return natural.integral(p.ds(potential, s) - p.half(p.E(t, omega)), 'A=0 D-layer curvature')


def g(k, triple, s, omega, a_zero=False, b_zero=False):
    """g_k(A,B,C) of (A13) in degree k: J or the prism G below the cutoff, J_6 at it.

    ``a_zero`` and ``b_zero`` say whether the A and B layers vanish on R; a
    lower-legal triple with A=0 in degree 5 or 6 uses the direct A=0 formulas.
    """
    if k in (5, 6) and a_zero and triple.lower_legal:
        _check(k, triple, 'differential')
        return a0_curvature(k, triple, s, omega, b_zero)
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


def legal_gamma(u, v, sum_legal, sum_pure_closed, s, omega, pure_sector=False):
    """gamma_k of two lower-legal triples in degrees 4 to 6, without the pure-C normalization.

    This is the legal branch of all_cochain_upper.gamma with
    closed_ab_upper.gamma in place of the normalized correction; the C
    residuals of the two inputs may be nonzero. With ``pure_sector`` (both
    triples fully legal with A=B=0 on R, degrees 5 and 6) the correction is
    compatible_sector.pure_c_gamma, which differs from closed_ab_upper.gamma
    by the integral coboundary of the module docstring.
    """
    total = upper.product(u, v, sum_legal, sum_pure_closed, s, omega)
    change = p.zero(u.A.degree + 4) - upper.K(u, s) - upper.K(v, s) + upper.K(total, s)
    t, tp, ts = [upper.current_curvature(w, s, omega)[2] for w in (u, v, total)]
    # All three images are closed pure-C triples; the entire half-lift
    # difference is an integer carry before taking any coboundary.
    carry = natural.integral(p.half(p.cup(s, t)) + p.half(p.cup(s, tp))
                             - p.half(p.cup(s, ts)), 'pure-C successor carry')
    if pure_sector:
        base = p.Cochain(u.A.degree + 4, sector.pure_c_gamma(
            _local(u.C), _local(v.C), _local(s), _local(omega)))
    else:
        base = closed.gamma(u.A, u.B, u.C, v.A, v.B, v.C, s, omega)
    return base - carry + change


def gamma(k, u, v, sum_legal, sum_pure_closed, a_zero, s, omega, b_zero=False):
    """gamma_k of two triples in degree k.

    ``a_zero`` is the pair of flags saying whether the A layer of ``u``,
    respectively of ``v``, vanishes on R (a single flag stands for both), and
    ``b_zero`` says whether both B layers vanish on R.
    """
    if isinstance(a_zero, bool):
        a_zero = (a_zero, a_zero)
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
    if k == 6 and not any(a_zero) and not A_STACKING:
        raise PairSourceLimit(
            'the degree-six stacking correction of two states with nonzero A layers '
            'evaluates the universal pair source of production_gamma6, whose terms '
            'nest the V3 contraction and which is not practical to evaluate; set '
            'FERMIONAHSS_DEGREE_SIX_A_STACKING=1 to evaluate it regardless of its '
            'running time')
    pure_sector = (k in (5, 6) and all(a_zero) and b_zero
                   and bool(u.full_legal) and bool(v.full_legal))
    return legal_gamma(u, v, sum_legal, sum_pure_closed, s, omega, pure_sector)
