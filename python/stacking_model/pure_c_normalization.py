"""One exact pure-C representative for the legal upper products k4--k6.

The correction is an explicit integral coboundary.  It normalizes the
phase before taking a prism; it never tests a global pure-C condition.
"""
import closed_ab_upper as closed
from natural_upper import integral
p,upper=closed.p,closed.upper


def phase_correction(A,B,C,Ap,Bp,Cp,s,omega):
    n=A.degree
    pure=closed.phase(p.zero(n),p.zero(n+1),C,
                      p.zero(n),p.zero(n+1),Cp,s,omega)
    return p.half(upper.hE(C,Cp))-pure


def phase(A,B,C,Ap,Bp,Cp,s,omega):
    return (closed.phase(A,B,C,Ap,Bp,Cp,s,omega)
            +phase_correction(A,B,C,Ap,Bp,Cp,s,omega))


def correction(A,B,C,Ap,Bp,Cp,s,omega):
    return integral(p.ds(phase_correction(A,B,C,Ap,Bp,Cp,s,omega),s),
                    'pure-C normalization')


def integer_phase_correction(A,B,C,Ap,Bp,Cp,s,omega):
    """An explicit integral primitive of ``correction``.

    chi is closed modulo one on unrestricted C,Cprime.  The normalized
    prism therefore gives chi-delta_s I chi integral, with boundary
    delta_s chi.  This proves that the normalization preserves the
    integral equivalence class, not only the obstruction equation.
    """
    I=upper.hp.interval
    chi=phase_correction(A,B,C,Ap,Bp,Cp,s,omega)
    along=phase_correction(I(A),I(B),I(C,True),I(Ap),I(Bp),I(Cp,True),I(s),I(omega))
    return integral(chi-p.ds(upper.hp.prism(along),s),'normalization primitive')


def gamma(A,B,C,Ap,Bp,Cp,s,omega):
    if A.degree==1:
        from production_gamma4 import gamma4 as formula
    elif A.degree==2:
        from production_gamma5 import gamma as formula
    elif A.degree==3:
        from production_gamma6 import gamma as formula
    else:
        raise ValueError('production normalization covers k4,k5,k6')
    return (formula(A,B,C,Ap,Bp,Cp,s,omega)
            +correction(A,B,C,Ap,Bp,Cp,s,omega))


def closed_gamma(A,B,C,Ap,Bp,Cp,s,omega):
    return (closed.gamma(A,B,C,Ap,Bp,Cp,s,omega)
            +correction(A,B,C,Ap,Bp,Cp,s,omega))
