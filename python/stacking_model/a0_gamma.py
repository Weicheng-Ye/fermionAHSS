"""Integral A=0 stacking at k=2 for arbitrary B,C,D.

The successor gamma at k=3 is supplied on legal lower data, which includes
every differential image of a k=2 input. This does not implement the full
all-cochain k=3 product. Fixed normalized prisms and exact rational lifts
are used; no cochain equation on X is solved.

The degree-one Theta specialization follows the MIT-licensed production
formulas in fermionAHSS/python/phase_eval.py; its calibrated chi_1 is zero.
"""
from fractions import Fraction as F

from cochains import (Cochain, binary_sum, cup, differential,
                      signed_differential, word_op, zero)
from compatible_sector import E, QD, integral
from lower_stacking import pullback_interval, prism
from paper_stacking import repaired_paper_corrections


def theta_degree_one(b, s, omega):
    """The production Theta(b,0), as its exact real lift, for closed b in C1."""
    if (b.degree, s.degree, omega.degree) != (1, 1, 2):
        raise ValueError("expected degrees (1,1,2)")
    b2 = cup(b, b)
    wb, sb2 = cup(omega, b), cup(s, b2)
    h = binary_sum(word_op("123134", omega, omega, b, b),
                   cup(wb, sb2, 2),
                   word_op("123243", s, s, b2, b2),
                   cup(cup(omega, s, 1), b2),
                   cup(cup(s, s), b2))
    return (F(1, 2)*h + F(3, 4)*cup(omega, b2)
            + F(1, 4)*cup(b2, b2))


def omega3(b, c, s, omega):
    """Exact production Omega3 at A=0, with closed degree-one B."""
    return F(1, 2)*E(c, omega) + theta_degree_one(b, s, omega)


def _paper_c_gauge(c, s):
    return F(1, 2)*binary_sum(cup(s, c),
                              cup(c, differential(c).mod2(), 2))


def _eighth_cube(b):
    return F(1, 8)*cup(cup(b, b, integral=True), b, integral=True)


def source_polarization_phase(b, bp, s):
    """Primitive of the polarization of h(s B^3+s^2 B^2)."""
    lam = cup(b, bp, 1)
    return (F(1, 4)*cup(s, cup(b, bp, integral=True), integral=True)
            + F(1, 2)*binary_sum(
                cup(s, binary_sum(cup(lam, b), cup(bp, lam))),
                cup(cup(s, s), lam)))


def phase3_production(b, c, bp, cp, s, omega):
    """Legal k3 phase in the production Omega gauge.

    Requires delta B=delta B'=0, delta C=QD(B), delta C'=QD(B').
    The formula can be evaluated diagnostically elsewhere, but its phase
    identity is asserted only on this domain.
    """
    correction = repaired_paper_corrections(3, b, c, bp, cp, s, omega)
    total_b = binary_sum(b, bp)
    total_c = binary_sum(c, cp, correction.beta)
    return (correction.phase - _paper_c_gauge(total_c, s)
            + _paper_c_gauge(c, s) + _paper_c_gauge(cp, s)
            + _eighth_cube(total_b) - _eighth_cube(b) - _eighth_cube(bp)
            + source_polarization_phase(b, bp, s))


def lower_d2(b, c, s, omega):
    """The first two nonzero output coordinates at k2, on arbitrary inputs."""
    if (b.degree, c.degree) != (0, 1):
        raise ValueError("expected k2 degrees B=0,C=1")
    return (differential(b).mod2(),
            binary_sum(differential(c).mod2(), QD(b, s, omega)))


def beta2(b, bp, s):
    return binary_sum(cup(differential(b).mod2(), bp), cup(s, cup(b, bp)))


def beta3_closed(b, bp, s):
    return binary_sum(cup(b, bp), cup(s, cup(b, bp, 1)))


def _interval_image2(b, c, s, omega):
    si, wi = pullback_interval(s), pullback_interval(omega)
    bi, ci = lower_d2(pullback_interval(b, True),
                      pullback_interval(c, True), si, wi)
    return bi, ci, si, wi


def raw_phi2(b, c, s, omega):
    """One-prism rational potential, valid for every k2 lower input."""
    bi, ci, si, wi = _interval_image2(b, c, s, omega)
    return prism(omega3(bi, ci, si, wi))


def phi2(b, c, s, omega, *, b_closed):
    """The current piecewise differential's potential.

    b_closed records global closedness on X, not merely on an evaluated
    simplex. Its two branches agree modulo integral cochains.
    """
    if not b_closed:
        return raw_phi2(b, c, s, omega)
    tau = cup(omega, b)
    defect = binary_sum(differential(c).mod2(), tau)
    return F(1, 2)*E(c, omega) + F(1, 2)*cup(defect, tau, 1)


def g2(b, c, s, omega, *, b_closed):
    """Exact current piecewise g2, evaluated with its rational potential."""
    db, dc = lower_d2(b, c, s, omega)
    return integral(signed_differential(phi2(b, c, s, omega,
                                             b_closed=b_closed), s)
                    - omega3(db, dc, s, omega), "g2")


def phase2_all_cochains(b, c, bp, cp, s, omega):
    """Explicit phase coordinate used to make gamma2 integral."""
    bi, ci, si, wi = _interval_image2(b, c, s, omega)
    bpi, cpi, _, _ = _interval_image2(bp, cp, s, omega)
    # Stacking the two ramps and ramping their stack differ in C by eta.
    # The nonclosed-B term delta B cup B' in beta2 is responsible for it.
    eta = binary_sum(beta2(pullback_interval(b, True),
                            pullback_interval(bp, True), si),
                     pullback_interval(beta2(b, bp, s), True))
    total_ci = binary_sum(ci, cpi, beta3_closed(bi, bpi, si))
    return -prism(phase3_production(bi, ci, bpi, cpi, si, wi)
                  + phase3_C_transport(total_ci, eta, wi))


def phase3_C_transport(c, lam, omega):
    """h(E lambda+H(C,delta lambda)), transporting C to C+delta lambda."""
    dl = differential(lam).mod2()
    return F(1, 2)*binary_sum(E(lam, omega), cup(c, dl, 1),
                              cup(differential(c).mod2(), dl, 2))


def gamma3_legal(b, c, bp, cp, s, omega):
    """Integral successor correction on full legal k3 lower data."""
    total_b = binary_sum(b, bp)
    total_c = binary_sum(c, cp, beta3_closed(b, bp, s))
    return integral(omega3(b, c, s, omega) + omega3(bp, cp, s, omega)
                    - omega3(total_b, total_c, s, omega)
                    + signed_differential(
                        phase3_production(b, c, bp, cp, s, omega), s),
                    "gamma3 on legal lower data")


def gamma2(b, c, bp, cp, s, omega, *, b_closed, bp_closed, sum_closed):
    """Integral k2 correction, with strict compatibility on all cochains.

    The successor product is gamma3_legal. The three flags are global
    closedness of B, B', and B+B', exactly as for the original piecewise d.
    D,D' enter only additively and can be arbitrary integral cochains.
    """
    total_b = binary_sum(b, bp)
    total_c = binary_sum(c, cp, beta2(b, bp, s))
    db, dc = lower_d2(b, c, s, omega)
    dbp, dcp = lower_d2(bp, cp, s, omega)
    change = (phi2(b, c, s, omega, b_closed=b_closed)
              + phi2(bp, cp, s, omega, b_closed=bp_closed)
              - phi2(total_b, total_c, s, omega, b_closed=sum_closed))
    return integral(change + signed_differential(
        phase2_all_cochains(b, c, bp, cp, s, omega), s)
        + phase3_production(db, dc, dbp, dcp, s, omega), "gamma2")
