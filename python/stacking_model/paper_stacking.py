"""Phase-coordinate stacking of arXiv:2310.19058v2, for A=0 and k=1,2,3.

The correction ``phase`` has degree k and is valued in Q/Z; it is NOT the
degree-(k+1) integral correction gamma for a Bockstein coordinate D.
Formulas may be evaluated off shell for diagnostics, but the paper proves
compatibility only on its defining data. No all-cochain chain-map claim is
made here. A generic-s k=3 transcription check currently fails; the regression
test preserves the witness. That part is experimental, not verified stacking.
Source: https://arxiv.org/html/2310.19058v2
"""
from dataclasses import dataclass
from fractions import Fraction as F

from cochains import Cochain, binary_sum, cup, differential, signed_differential, zero


@dataclass(frozen=True)
class PaperCorrections:
    alpha: Cochain
    beta: Cochain
    phase: Cochain


def _degrees(k, b, c, bp, cp, s, omega):
    if k not in (1, 2, 3):
        raise NotImplementedError("the paper supplies A=0 formulas only for k=1,2,3")
    actual = (b.degree, c.degree, bp.degree, cp.degree, s.degree, omega.degree)
    if actual != (k - 2, k - 1, k - 2, k - 1, 1, 2):
        raise ValueError("expected deg(B)=k-2, deg(C)=k-1, deg(s)=1, deg(omega)=2")


def _triple(a, b, c):
    return cup(cup(a, b), c)


def paper_corrections(k, b, c, bp, cp, s, omega):
    """Return alpha=0, beta, and additive rational phase e (modulo integers).

    B,C,B',C',s,omega are reduced to binary cochains. A is absent in this API.
    The returned rational lift is deliberately NOT reduced modulo one: its
    integral coboundary carries are needed when converting to D coordinates.
    For k=3 and nonzero s this is an experimental transcription with a known
    consistency failure; see test_generic_s_transcription_mismatch_is_exposed.
    """
    _degrees(k, b, c, bp, cp, s, omega)
    b, c, bp, cp, s, omega = (a.mod2() for a in (b, c, bp, cp, s, omega))
    if k == 1:
        return PaperCorrections(zero(-1), zero(0), zero(1))
    if k == 2:
        m = _triple(s, b, bp)
        epsilon = binary_sum(cup(m, binary_sum(c, cp)), cup(c, cp),
                             cup(c, differential(cp).mod2(), 1))
        phase = (F(1, 2) * epsilon + F(3, 4) * _triple(omega, b, bp)
                 + F(1, 2) * cup(omega, m, 1))
        return PaperCorrections(zero(0), m, phase)

    h = cup(b, bp, 1)
    total_b = binary_sum(b, bp)
    m = binary_sum(cup(b, bp), cup(s, h))
    total_c = binary_sum(c, cp, m)
    dc, dcp = differential(c).mod2(), differential(cp).mod2()
    epsilon = binary_sum(
        cup(c, cp, 1), cup(dc, c, 2), cup(cp, dc, 2), cup(dcp, cp, 2),
        cup(binary_sum(dc, dcp), total_c, 2), cup(m, total_c, 1),
        cup(m, differential(total_c).mod2(), 2))
    theta0 = F(1, 2) * _triple(b, h, bp) + F(1, 4) * _triple(b, bp, bp)
    theta_omega = (F(1, 2) * cup(cup(omega, total_b), cup(b, bp), 2)
                   + F(3, 4) * cup(omega, h))
    theta_s_half = binary_sum(
        _triple(s, s, h), _triple(s, cup(cup(s, b, 1), bp, 1), h),
        _triple(b, cup(s, bp, 1), h), _triple(cup(s, b, 1), total_b, h),
        _triple(cup(s, b, 1), h, b), _triple(s, h, bp))
    # Add the two integer-valued binary lifts before multiplying by 1/4;
    # total_b and every individual cup expression are reduced modulo two.
    theta_s = (F(1, 2) * theta_s_half + F(1, 4) * _triple(s, bp, b)
               + F(1, 4) * _triple(s, h, total_b) + F(5, 8) * _triple(s, b, bp))
    theta_mixed = F(1, 2) * binary_sum(
        cup(cup(omega, s, 1), h), cup(cup(omega, cup(s, h), 2), total_b),
        cup(cup(omega, cup(s, bp), 2), b),
        cup(cup(omega, cup(s, b), 2), binary_sum(bp, h)))
    phase = F(1, 2) * epsilon + theta0 + theta_omega + theta_s + theta_mixed
    return PaperCorrections(zero(1), m, phase)


def antiunitary_phase_repair(b, bp, s):
    """Derived local degree-three correction chi, for closed B,B',s.

    On (0123), write x=(B01,B12,B23), y=(B'01,B'12,B'23),
    and t=(s01,s12). All sums below are integer sums of binary lifts.
    This representative was solved exactly over Z/4 from the universal
    F2^3 cochain consistency equations, with chi=0 when s, B, or B' is zero.
    It is an independent repair in these cup conventions, not an equation
    attributed to the paper. No off-shell or higher-degree assertion follows.
    """
    if (b.degree, bp.degree, s.degree) != (1, 1, 1):
        raise ValueError("B, B', and s must have degree one")
    b, bp, s = b.mod2(), bp.mod2(), s.mod2()

    def evaluate(t):
        edges = tuple((t[j], t[j + 1]) for j in range(3))
        x1, x2, x3 = (b(e) for e in edges)
        y1, y2, y3 = (bp(e) for e in edges)
        t1, t2 = s(edges[0]), s(edges[1])
        quarter = t2 * x1 * y3 + t1 * (x2 + x3) * y2 * y3
        half = (t2 * y3 * (y1 * (x1 + x2 + x3) + x1 * (x2 + x3 + y2 + t1))
                + t1 * x3 * y2 * (y1 + x1))
        return F(1, 4) * quarter + F(1, 2) * half

    return Cochain(3, evaluate)


def repaired_paper_corrections(k, b, c, bp, cp, s, omega):
    """Paper lower layers and phase, with the derived k=3 antiunitary repair.

    The correction has been checked on legal data only. This function is not
    an all-cochain differential homomorphism and does not cover nonzero A.
    """
    original = paper_corrections(k, b, c, bp, cp, s, omega)
    if k != 3:
        return original
    return PaperCorrections(original.alpha, original.beta,
                            original.phase + antiunitary_phase_repair(b, bp, s))


def lower_obstruction(k, b, s, omega):
    """Paper's equation dC=f(B) on closed B, at A=0."""
    if k == 1:
        return zero(1)
    if k == 2:
        return cup(omega, b)
    if k == 3:
        return binary_sum(cup(omega, b), _triple(s, b, b))
    raise NotImplementedError("k must be 1, 2, or 3")


def obstruction_phase(k, c, s, omega):
    """The paper's exact obstruction gauge: Eqs.12,18,75 as additive Q/Z."""
    c, s, omega = (a.mod2() for a in (c, s, omega))
    if c.degree != k - 1 or s.degree != 1 or omega.degree != 2:
        raise ValueError("incorrect cochain degrees")
    if k == 1:
        return F(1, 2) * cup(omega, c)
    dc = differential(c).mod2()
    if k == 2:
        return F(1, 2) * binary_sum(cup(omega, c), cup(dc, dc, 1),
                                  cup(dc, c), cup(s, dc))
    if k != 3:
        raise NotImplementedError("k must be 1, 2, or 3")
    sign = binary_sum(cup(omega, c), cup(c, c), cup(c, dc, 1),
                      differential(binary_sum(cup(s, c), cup(c, dc, 2))).mod2())

    def exceptional_sign(t):
        a, b, d, e, f = t
        return (omega((a, b, e)) * dc((b, d, e, f))
                + dc((a, b, d, f)) * dc((a, d, e, f))) % 2

    def quarter(t):
        a, b, d, e, f = t
        return dc((a, b, d, e)) * (1 - dc((a, b, d, f)))

    return (F(1, 2) * binary_sum(sign, Cochain(4, exceptional_sign))
            - F(1, 4) * Cochain(4, quarter))


def phase_consistency_residual(k, b, c, bp, cp, s, omega):
    """ds(e)-Omega(total)+Omega(left)+Omega(right), valued in Q/Z.

    The paper claims vanishing on legal data. This implementation verifies it
    for k=1,2 and k=3 with s=0; a generic-s k=3 witness currently fails.
    """
    correction = paper_corrections(k, b, c, bp, cp, s, omega)
    total_c = binary_sum(c, cp, correction.beta)
    return (signed_differential(correction.phase, s)
            - obstruction_phase(k, total_c, s, omega)
            + obstruction_phase(k, c, s, omega)
            + obstruction_phase(k, cp, s, omega)).mod1()


def beta_swap_gauge(k, b, bp):
    """For closed B,B': beta+beta_swapped=d(gauge); k=2 is symmetric.

    At k=3 the gauge changes C. The upper phase swap therefore also requires
    the corresponding C-gauge transport, not just subtracting a coboundary.
    """
    if k in (1, 2):
        return zero(k - 2)
    if k == 3:
        return cup(b.mod2(), bp.mod2(), 1)
    raise NotImplementedError("k must be 1, 2, or 3")


def phase_swap_gauge_k2(b, c, bp, cp):
    """sigma with e2-e2_swapped = ds(sigma) mod 1, for legal k=2 data.

    B,B' are closed degree-zero cochains and dC=omega*B, dC'=omega*B'.
    beta is already symmetric here, so no additional C-gauge transport occurs.
    """
    if (b.degree, c.degree, bp.degree, cp.degree) != (0, 1, 0, 1):
        raise ValueError("expected k=2 data")
    return F(1, 2) * binary_sum(cup(c.mod2(), cp.mod2(), 1), cup(c.mod2(), bp.mod2()))


def phase_gauge_transport_k3(c, lam, s, omega):
    """K(C,lambda) with ds K=Omega(C+d lambda)-Omega(C), modulo one.

    This changes the degree-two C representative. The expression uses the
    obstruction gauge in Eq.75 and is meaningful without choosing a basis.
    """
    if (c.degree, lam.degree, s.degree, omega.degree) != (2, 1, 1, 2):
        raise ValueError("expected degrees (2,1,1,2)")
    c, lam, s, omega = (a.mod2() for a in (c, lam, s, omega))
    dc, dl = differential(c).mod2(), differential(lam).mod2()
    return F(1, 2) * binary_sum(cup(lam, dl), cup(omega, lam), cup(c, dl, 1),
                                cup(dc, dl, 2), cup(s, dl), cup(dl, dc, 2))


def phase_swap_gauge_k3(b, c, bp, cp, s):
    """Explicit sigma2 implementing prime exchange of the repaired k=3 law.

    Set lambda=B cup1 B', N=C+C'+beta. On legal data, modulo integers,
      e(B',C';B,C)-e(B,C;B',C')
        = K(N,lambda) + ds sigma2.
    The C output changes by d lambda at the same time. The formula is natural
    in local simplex values and contains no global cohomology-basis choice.
    """
    if (b.degree, c.degree, bp.degree, cp.degree, s.degree) != (1, 2, 1, 2, 1):
        raise ValueError("expected k=3 data")
    b, c, bp, cp, s = (a.mod2() for a in (b, c, bp, cp, s))

    def evaluate(t):
        x1, x2 = b(t[:2]), b(t[1:])
        y1, y2 = bp(t[:2]), bp(t[1:])
        st = s(t[:2])
        cv, cpv = c(t), cp(t)
        complex_part = F(1, 2) * (cv * cpv + (cv + cpv) * (x1 * y2 + x2 * y1))
        majorana_part = (
            F(1, 8) * y2 * (x1 + st * x2)
            + F(1, 4) * (3 * st * x2 * y1 + x1 * x2 * (y1 + 3 * y2) + st * x1 * y2)
            + F(1, 2) * (st * x1 * x2 * (y1 + y2) + st * x1 * y1 * y2
                         + (1 - st) * x2 * y1 * y2))
        return complex_part + majorana_part

    return Cochain(2, evaluate)
