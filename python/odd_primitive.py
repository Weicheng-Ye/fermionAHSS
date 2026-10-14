"""Current odd-branch source primitive, with its specified gauge prism.

This implements current r1.md Section 5, not the older source polynomial
in the imported Appendix B. The finite detector values and Fill operator
are supplied by the verification driver; no value is inferred from KO.
"""
from fractions import Fraction

try:
    from . import phase_eval as p
except ImportError:
    import phase_eval as p


def _integral(c, name):
    def evaluate(vertices):
        value = c(vertices)
        if getattr(value, "denominator", 1) != 1:
            raise ArithmeticError(f"{name} is not integral: {value}")
        return int(value)
    return p.Cochain(c.degree, evaluate)


def build_odd(A, s, omega):
    """Return Phi, Vg, kappa, z and Bs for the literal odd branch A=~s.

    Input vertices are (source_id, index, *existing_prism_levels).
    The new gauge level is appended, and all basic input cochains ignore
    it. Every actual prism simplex begins at its bottom first vertex,
    hence no unprinted transport multiplier is introduced in its sum.
    """
    if (A.degree, s.degree, omega.degree) != (1, 1, 2):
        raise ValueError("odd branch expects degrees (1,1,2)")

    def checked_A(vertices):
        av, sv = A(vertices), s(vertices)
        if av != sv or sv not in (0, 1):
            raise ValueError("the odd branch requires the literal input A=~s")
        return av

    A_odd = p.Cochain(1, checked_A)
    base_source = p.source1(A_odd, s, omega)
    Phi = p.theta(base_source["q"], base_source["k"], s, omega)

    tau = p.Cochain(0, lambda vertices: vertices[0][-1])
    sI = p.binary(s + p.differential(tau))
    AI = sI
    gauge_source = p.source1(AI, sI, omega)
    PhiI = p.theta(gauge_source["q"], gauge_source["k"], sI, omega)
    II = _integral(p.ds(PhiI, sI), "odd gauge source differential")
    KPhi = p.prism_cochain(PhiI)
    KI = _integral(p.prism_cochain(II), "odd gauge integral prism")
    Vg = p.scale(KPhi, Fraction(-1, 2))
    kappa = p.binary(KI)

    R = p.binary(p.divide(p.differential(omega), 2,
                          "odd background Bockstein"))
    z = (p.cup(p.cup(s, s), R), p.cup(omega, R), p.square(R, 2))
    Bs = tuple(p.binary(p.divide(p.ds(zi, s), 2,
                                "odd twisted binary Bockstein"))
               for zi in z)
    return dict(A=A_odd, s=s, omega=omega, source=base_source,
                Phi=Phi, Vg=Vg, kappa=kappa, z=z, Bs=Bs,
                R=R, sI=sI, AI=AI, gauge_source=gauge_source,
                PhiI=PhiI, II=II, KI=KI)


def _coefficients(lambdas):
    values = tuple(lambdas)
    if len(values) != 8 or any(value not in (0, 1) for value in values):
        raise ValueError("eight explicitly evaluated binary detector values required")
    l1, l2, l3, l4, l5, l6, l7, l8 = values
    if l8 or (l5 + l6 + l7) % 2 or (l1 + l2 + l4) % 2:
        raise ArithmeticError("current odd-source detector relations failed")
    if l1 or l4:
        raise ArithmeticError("current proved odd-source obstruction bits are nonzero")
    return (l3, (l7 + l5) % 2, l5)


def build_gamma(data, lambdas):
    """The current zero-obstruction residual after its three Bocksteins."""
    alpha = _coefficients(lambdas)
    result = data["kappa"]
    for bit, Bs in zip(alpha, data["Bs"]):
        if bit:
            result = result + Bs
    return p.binary(result)


def build_Uraw(data, lambdas, fill):
    """Assemble Uraw from actual detector bits and the printed finite Fill.

    fill receives the degree-six binary cochain gamma and must return
    the degree-five binary cochain prescribed by Appendix B.6.
    Current mu1=mu2=0 is checked against the supplied detector values.
    """
    alpha = _coefficients(lambdas)
    gamma = build_gamma(data, lambdas)
    filled = fill(gamma)
    if filled.degree != 5:
        raise ValueError("odd Fill must return a degree-five cochain")
    Uraw = data["Vg"] + p.half(filled)
    for bit, z in zip(alpha, data["z"]):
        if bit:
            Uraw = Uraw + p.scale(z, Fraction(1, 4))
    return Uraw
