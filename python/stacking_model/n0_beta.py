"""Calibrated legal-pair beta at input degree n=0 (physical degree k=3).

This supplies the formerly missing universal A-pair source primitive.
It does not supply beta on arbitrary nonclosed A,B or an integral gamma.
The inputs obey delta_s A=delta_s A'=0 and delta B=QD(rho A), likewise
for B'. All divisions are exact and evaluated lazily.
"""
from dataclasses import dataclass
from fractions import Fraction

from cochains import Cochain, binary_sum, cup, differential, signed_differential
from compatible_sector import QD, alpha as h_D, integral


@dataclass(frozen=True)
class Source0:
    a: Cochain
    t: Cochain
    g: Cochain
    P: Cochain
    k0: Cochain


def source0(A, s, omega):
    """The degree-zero specialization of the fixed secondary source split.

    For a signed integral 0-cocycle, g(012)=t(0)*omega(012)+a*s(01)*s(12).
    Thus the source depends only on A modulo four, including negative A.
    """
    if (A.degree, s.degree, omega.degree) != (0, 1, 2):
        raise ValueError("source0 expects degrees (0,1,2)")
    a = A.mod2()
    t = integral(Fraction(1, 2) * (A-a), "A carry").mod2()
    g = binary_sum(cup(t, omega), cup(cup(a, s), s))
    P = cup(omega, a)
    numerator = -signed_differential(g, s) + cup(s, P, integral=True)
    k0 = binary_sum(integral(Fraction(1, 2) * numerator, "source L0"),
                    cup(cup(cup(s, s), s), a))
    return Source0(a, t, g, P, k0)


def source_primitive(A, Aprime, s, omega):
    """Explicit normalized universal V with delta V=S(A,A').

    On (012), a,t are the two bits of A(0) modulo four, b,u those of A'(0),
    x=s(01), y=s(12), and w=omega(012). The six-monomial polynomial is
    a*u*x*y + w*(t*u + a*b*(x+y) + (t*b+a*u)*x*y), modulo two.
    """
    if (A.degree, Aprime.degree, s.degree, omega.degree) != (0, 0, 1, 2):
        raise ValueError("source_primitive expects degrees (0,0,1,2)")

    def evaluate(face):
        value, valueprime = Fraction(A(face[:1])), Fraction(Aprime(face[:1]))
        if value.denominator != 1 or valueprime.denominator != 1:
            raise ArithmeticError("A and A' must be integral")
        left, right = value.numerator % 4, valueprime.numerator % 4
        a, t, b, u = left % 2, left // 2, right % 2, right // 2
        x, y, w = s(face[:2]) % 2, s(face[1:]) % 2, omega(face) % 2
        return (a*u*x*y + w*(t*u + a*b*(x+y) + (t*b+a*u)*x*y)) % 2

    return Cochain(2, evaluate, "V0")


def splitting_carry(source, B, s):
    """Canonical carry lambda=([g]+[sB]-[g+sB])/2, reduced modulo two."""
    summand = source.g + cup(s, B)
    return integral(Fraction(1, 2) * (summand-summand.mod2()), "lambda").mod2()


def source_polarization(A, Aprime, s, omega):
    """The exact A-pair source S whose primitive is source_primitive."""
    left, right, total = (source0(x, s, omega) for x in (A, Aprime, A+Aprime))
    correction = h_D(A, Aprime, s)
    return binary_sum(total.k0, left.k0, right.k0, QD(correction, s, omega),
                      h_D(left.P, right.P, s),
                      h_D(binary_sum(left.P, right.P),
                          differential(correction).mod2(), s))


def beta_legal(A, B, Aprime, Bprime, s, omega):
    """Full beta on calibrated legal n=0 lower pairs, including nonzero A.

    The domain is delta_s A=0 and delta B=QD(rho A), for each input.
    Closedness is a global condition that the caller checks on its complex.
    The returned degree-two binary cochain obeys delta beta=Delta tau.
    """
    if (A.degree, B.degree, Aprime.degree, Bprime.degree,
            s.degree, omega.degree) != (0, 1, 0, 1, 1, 2):
        raise ValueError("beta_legal expects degrees (0,1,0,1,1,2)")
    correction = h_D(A, Aprime, s)
    Bsum = binary_sum(B, Bprime, correction)
    left, right, total = (source0(x, s, omega) for x in (A, Aprime, A+Aprime))
    return binary_sum(h_D(B, Bprime, s),
                      h_D(binary_sum(B, Bprime), correction, s),
                      splitting_carry(left, B, s),
                      splitting_carry(right, Bprime, s),
                      splitting_carry(total, Bsum, s),
                      source_primitive(A, Aprime, s, omega))
