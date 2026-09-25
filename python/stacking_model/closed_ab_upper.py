"""Strict upper stacking on legal (A,B), with completely arbitrary C,D.

This uses the actual production J branch in physical degrees4,5,6.
The all-cochain A,B extension is a separate construction.
"""
from fractions import Fraction as F
import closed_a_upper as upper
from upper_phase_diagnostic import production_phase
p=upper.p


def phase(A,B,C,Ap,Bp,Cp,s,omega):
    if A.degree==1:
        from production_gamma4 import phase4 as formula
    elif A.degree==2:
        from production_gamma5 import phase as formula
    elif A.degree==3:
        from production_gamma6 import phase as formula
    else:
        raise ValueError('this closed-AB construction covers k4,k5,k6')
    return formula(A,B,C,Ap,Bp,Cp,s,omega)


def curvature(A,B,C,s,omega):
    return p.binary(p.differential(C)+upper.fsharp(A,B,s,omega))


def potential(A,B,C,s,omega):
    n=A.degree;tau=upper.fsharp(A,B,s,omega)
    t=p.binary(p.differential(C)+tau)
    return production_phase(n,A,B,C,s,omega)+p.half(p.binary(p.cup(t,tau,n+2)))


def J(A,B,C,s,omega):
    t=curvature(A,B,C,s,omega)
    raw=p.ds(potential(A,B,C,s,omega),s)-p.half(p.E(t,omega))
    def value(vertices):
        result=F(raw(vertices))
        if result.denominator!=1:raise ArithmeticError(f'J is not integral: {result}')
        return result.numerator
    return p.Cochain(A.degree+5,value)


def gamma(A,B,C,Ap,Bp,Cp,s,omega):
    alpha=upper.hp.hD(p.binary(A),p.binary(Ap),s)
    beta=upper.beta_sharp(A,B,Ap,Bp,s,omega)
    As,Bs,Cs=A+Ap,p.binary(B+Bp+alpha),p.binary(C+Cp+beta)
    t,tp=curvature(A,B,C,s,omega),curvature(Ap,Bp,Cp,s,omega)
    # This plus sign fixes the exact chosen pure-C successor, even though
    # its opposite sign is equivalent modulo integers.
    successor_phase=p.half(upper.hE(t,tp))
    raw=(potential(A,B,C,s,omega)+potential(Ap,Bp,Cp,s,omega)
         -potential(As,Bs,Cs,s,omega)+p.ds(phase(A,B,C,Ap,Bp,Cp,s,omega),s)
         +successor_phase)
    def value(vertices):
        result=F(raw(vertices))
        if result.denominator!=1:raise ArithmeticError(f'closed-AB gamma is not integral: {result}')
        return result.numerator
    return p.Cochain(A.degree+4,value)
