"""Strict all-cochain stacking in the sector A=B=0, and the primary alpha.

This is not an implementation of the full four-layer stacking problem.
See ALL_COCHAIN_OBSTRUCTION.md for an obstruction for the existing cutoff d.
"""
from dataclasses import dataclass
from fractions import Fraction

from cochains import (Cochain, binary_sum, cup, differential,
                      signed_differential)


def Q(x, j):
    """Cochain Steenrod square including its coboundary term."""
    r = x.degree
    return binary_sum(cup(x, x, r-j),
                      cup(x, differential(x).mod2(), r-j+1))


def E(x, omega):
    return binary_sum(Q(x, 2), cup(omega, x))


def QD(x, s, omega):
    return binary_sum(E(x, omega), cup(s, Q(x, 1)))


def polarization(x, y, j=2):
    """h^j_r(x,y), with delta h + h(delta x,delta y) = Delta Q^j."""
    if x.degree != y.degree:
        raise ValueError("polarization inputs must have the same degree")
    r = x.degree
    return binary_sum(cup(x, y, r-j+1),
                      cup(differential(x).mod2(), y, r-j+2))


def alpha(A, Aprime, s):
    """The primary correction only; no nonzero-A beta or gamma is supplied."""
    a, aprime = A.mod2(), Aprime.mod2()
    return binary_sum(polarization(a, aprime),
                      cup(s, polarization(a, aprime, 1)))


def integral(cochain, label):
    """Require exact integrality on evaluation, never truncate a fraction."""
    def evaluate(simplex):
        value = Fraction(cochain(simplex))
        if value.denominator != 1:
            raise ArithmeticError(f"{label} is not integral: {value}")
        return value.numerator
    return Cochain(cochain.degree, evaluate, label)


def pure_c_g(C, s, omega):
    """The A=B=0 restriction of the existing all-cochain differential."""
    dc = differential(C).mod2()
    return integral(Fraction(1, 2) * (
        signed_differential(E(C, omega), s) - E(dc, omega)), "g")


def pure_c_gamma(C, Cprime, s, omega):
    """Integral gamma satisfying strict compatibility on arbitrary cochains.

    The successor polarization is essential. Cochain products inside E and
    polarization are binary; the five summands here are added over Q.
    """
    dc, dcp = differential(C).mod2(), differential(Cprime).mod2()
    total = binary_sum(C, Cprime)
    return integral(Fraction(1, 2) * (
        E(C, omega) + E(Cprime, omega) - E(total, omega)
        + signed_differential(polarization(C, Cprime), s)
        + polarization(dc, dcp)), "gamma")


def closed_c_swap_primitive(C, Cprime, s):
    """Integral K with gamma(C,C')-gamma(C',C)=delta_s K for CLOSED C,C'.

    This primitive is defined only on closed inputs. The caller must check
    closedness on its complex; integrality is additionally checked lazily.
    """
    pointwise = cup(C, Cprime, C.degree)
    return integral(Fraction(1, 2) * (
        polarization(C, Cprime) - polarization(Cprime, C)
        - signed_differential(pointwise, s)), "closed-input swap primitive")


def swap_homotopy(C, Cprime, s):
    """Integral K on arbitrary inputs with gamma-gamma^op=ds K+K(delta,delta)."""
    dc, dcp = differential(C).mod2(), differential(Cprime).mod2()
    return integral(Fraction(1, 2) * (
        polarization(C, Cprime) - polarization(Cprime, C)
        - signed_differential(cup(C, Cprime, C.degree), s)
        + cup(dc, dcp, dc.degree)), "swap homotopy")


@dataclass(frozen=True)
class PureCState:
    """The tuple (0,0,C,D) in degree k=C.degree+1."""
    C: Cochain
    D: Cochain

    def __post_init__(self):
        if self.C.degree < 0 or self.D.degree != self.C.degree + 2:
            raise ValueError("require C of degree k-1>=0 and D of degree k+1")


def d(state, s, omega):
    return PureCState(differential(state.C).mod2(),
                      signed_differential(state.D, s)
                      + pure_c_g(state.C, s, omega))


def stack(left, right, s, omega):
    return PureCState(binary_sum(left.C, right.C),
                      left.D + right.D
                      + pure_c_gamma(left.C, right.C, s, omega))
