"""Closed legal pair-phase changes selected by the degree-three product.

The zero-rank selector is constant on every connected component of legal
signed rank cochains.  This correction is used only if both entire input
triples are legal; it is deliberately not an arbitrary-C natural formula.
"""
from cochains import Cochain, signed_differential
import gamma3_odd_symmetry as odd_exchange
from compatible_sector import integral
import closed_ab_degree3 as previous
import paper_commutative_production as commutative


def phase_change(A,B,C,Ap,Bp,Cp,s,omega):
    old=previous.phase(A,B,C,Ap,Bp,Cp,s,omega)
    new=commutative.phase(B,C,Bp,Cp,s,omega)
    odd_change=odd_exchange.phase_correction(A,B,Ap,Bp)
    return Cochain(3,lambda face:(new(face)-old(face))
                   if A(face[:1])==0 and Ap(face[:1])==0 else odd_change(face))


def gamma_change(A,B,C,Ap,Bp,Cp,s,omega):
    return integral(signed_differential(phase_change(A,B,C,Ap,Bp,Cp,s,omega),s),
                    'closed legal degree-three recalibration')


def exchange(A,B,C,Ap,Bp,Cp,s,omega):
    """Selected full-legal lambda,K,sigma,L,M on arbitrary components."""
    import gamma3_even_symmetry as even_exchange
    args=tuple(Cochain(c.degree,c)for c in(A,B,C,Ap,Bp,Cp,s,omega))
    A,B,C,Ap,Bp,Cp,s,omega=args
    data=(commutative.exchange(B,C,Bp,Cp,s,omega),
          even_exchange.exchange(*args),odd_exchange.exchange(*args))
    def branch(face):
        a,ap=A(face[:1]),Ap(face[:1])
        return 0 if a==ap==0 else 2 if a%2 or ap%2 else 1
    return tuple(Cochain(degree,lambda face,j=j:data[branch(face)][j](face))
                 for j,degree in enumerate((1,3,2,4,3)))
