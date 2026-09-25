"""Explicit A=0 lower stacking for the existing piecewise f.

Finite prisms give the off-shell coordinate correction. A caller supplies
the global closedness of B on its complex; a single evaluated simplex is
not enough to select that branch. No cochain equation is solved on X.
"""
from dataclasses import dataclass

from cochains import Cochain, binary_sum, differential
from compatible_sector import QD, alpha


def pullback_interval(c, multiply=False):
    """Normalized pullback to X times I, optionally times its last coordinate."""
    def value(vertices):
        base = tuple(v[0] for v in vertices)
        if any(a == b for a, b in zip(base, base[1:])):
            return 0
        result = c(base)
        return result * vertices[-1][1] if multiply else result
    return Cochain(c.degree, value)


def prism(c):
    """Signed right prism; callers reduce binary outputs explicitly."""
    def value(vertices):
        result = 0
        for j in range(len(vertices)):
            simplex = (tuple((v, 0) for v in vertices[:j+1])
                       + tuple((v, 1) for v in vertices[j:]))
            result += (-1)**j * c(simplex)
        return result
    return Cochain(c.degree-1, value)


def prism_gauge(B, s, omega):
    """K=I QD(B ell); delta K=H(0,B)+QD B by prism Stokes."""
    return prism(QD(pullback_interval(B, True), pullback_interval(s),
                    pullback_interval(omega))).mod2()


def raw_H_zero_A(B, s, omega):
    """The A=0 restriction of the original raw H; k2 has its stated convention."""
    if B.degree == 0:
        return QD(B, s, omega)
    return prism(QD(differential(pullback_interval(B, True)).mod2(),
                    pullback_interval(s), pullback_interval(omega))).mod2()


@dataclass(frozen=True)
class GlobalBinaryCochain:
    cochain: Cochain
    closed: bool


def off_shell_coordinate_gauge(B, s, omega):
    if B.closed or B.cochain.degree == 0:
        from cochains import zero
        return zero(B.cochain.degree+1)
    return prism_gauge(B.cochain, s, omega)


def f_zero_A(B, s, omega):
    return QD(B.cochain, s, omega) if B.closed else raw_H_zero_A(B.cochain, s, omega)


def beta_zero_A(B, Bprime, sum_closed, s, omega):
    """The exact third-component correction for arbitrary B,B' at A=A'=0.

    sum_closed must record global closedness of B+B'. Together with the two
    input flags this implements exactly the old piecewise f's branch rule.
    """
    total = GlobalBinaryCochain(binary_sum(B.cochain, Bprime.cochain), sum_closed)
    return binary_sum(alpha(B.cochain, Bprime.cochain, s),
                      off_shell_coordinate_gauge(B, s, omega),
                      off_shell_coordinate_gauge(Bprime, s, omega),
                      off_shell_coordinate_gauge(total, s, omega))
