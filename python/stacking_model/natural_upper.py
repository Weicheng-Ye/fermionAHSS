"""Natural nested-prism upper stacking, with prescribed legal successor.

This is an auxiliary all-cochain construction in physical degrees3,4,5.
The separate transport module restores the calibrated piecewise g.
"""
from fractions import Fraction as F
import off_shell_beta as off
import stacking_lower as lower
from upper_phase_diagnostic import production_phase
p,hp=off.p,off.hp


def curvature(A,B,C,s,omega):
    a,b=off.curvature(A,B,s,omega)
    return a,b,p.binary(p.differential(C)+off.natural_f(A,B,s,omega))


def product(A,B,C,Ap,Bp,Cp,s,omega):
    alpha=off.alpha(A,Ap,s)
    beta=off.raw_beta(A,B,Ap,Bp,s,omega,lower.legal_beta)
    return A+Ap,p.binary(B+Bp+alpha),p.binary(C+Cp+beta)


def successor_gamma(A,B,C,Ap,Bp,Cp,s,omega):
    if A.degree not in (1,2,3):
        raise ValueError('this auxiliary construction has successors k4,k5,k6')
    from pure_c_normalization import gamma
    return gamma(A,B,C,Ap,Bp,Cp,s,omega)


def integral(cochain,label):
    def value(vertices):
        answer=F(cochain(vertices))
        if answer.denominator!=1:raise ArithmeticError(f'{label} is not integral: {answer}')
        return answer.numerator
    return p.Cochain(cochain.degree,value)


def ghat(A,B,C,s,omega):
    if A.degree not in (0,1,2):raise ValueError('the auxiliary upper operation covers k3,k4,k5')
    I=hp.interval
    a,b,c=curvature(I(A,True),I(B,True),I(C,True),I(s),I(omega))
    T=p.ds(production_phase(A.degree+1,a,b,c,I(s),I(omega)),I(s))
    return integral(p.scale(hp.prism(T),-1),'natural upper differential')


def residual(A,B,C,Ap,Bp,Cp,s,omega):
    As,Bs,Cs=product(A,B,C,Ap,Bp,Cp,s,omega)
    u=curvature(A,B,C,s,omega);v=curvature(Ap,Bp,Cp,s,omega)
    return (ghat(A,B,C,s,omega)+ghat(Ap,Bp,Cp,s,omega)-ghat(As,Bs,Cs,s,omega)
            +successor_gamma(*u,*v,s,omega))


def gamma(A,B,C,Ap,Bp,Cp,s,omega):
    I=hp.interval
    source=residual(I(A,True),I(B,True),I(C,True),I(Ap,True),I(Bp,True),I(Cp,True),I(s),I(omega))
    return integral(hp.prism(source),'natural upper stacking')
