"""Exact binary algebra in the n=3 production/auxiliary phase comparison.

This records a binary identity only. A production top-phase comparison is
not yet completed by this identity.
"""
import closed_a_upper as upper
p=upper.p


def comparison_terms(B,k,lam,s,omega):
    m=B.degree+1
    P=p.binary(p.differential(B));y=p.QD(B,s,omega)
    t=p.binary(p.differential(lam));r=p.binary(p.Q(P,1)+p.cup(s,P))
    tau=p.binary(y+k+t)
    V0=p.binary(p.Q(B,1)+p.cup(s,B));Z0=p.binary(V0+P)
    L=p.binary(lam+Z0)
    pol=p.binary(p.cup(y,k,m)+p.cup(k,p.binary(p.differential(y)),m+1)+p.Q(k,1))
    production=p.binary(pol+p.E(lam,omega)+p.cup(p.binary(y+k),t,m)
                        +p.E(p.Q(B,1),omega)+p.cup(s,tau))
    theta_pair_C=p.binary(upper.hE(k,y)+upper.hE(p.binary(k+y),r))
    transport=p.binary(p.E(L,omega)+upper.hE(p.binary(k+y+r),p.binary(t+r)))
    residual=p.binary(production+p.Q(tau,1)+theta_pair_C+transport)
    J=p.binary(p.Q(y,1)+p.cup(s,y)+p.E(p.cup(s,B),omega))
    Aonly=p.binary(p.cup(s,k)+p.Q(r,1)+p.binary(p.differential(k))
                   +p.E(P,omega)+upper.hE(p.Q(P,1),p.cup(s,P)))
    primitive=p.binary(p.Q(lam,1)+p.cup(s,lam)+upper.hE(lam,Z0)
                       +p.cup(t,r,m+1)+p.cup(p.binary(y+k),t,m+1)
                       +upper.hE(V0,P)+upper.hE(p.Q(B,1),p.cup(s,B)))
    return residual,J,Aonly,primitive


def J0(B,s,omega):
    y=p.QD(B,s,omega)
    return p.binary(p.Q(y,1)+p.cup(s,y)+p.E(p.cup(s,B),omega))


def j0_homotopy(x,y,s,omega):
    """All-cochain additivity homotopy: delta j+j(delta x,delta y)=Delta J0."""
    dx,dy=p.binary(p.differential(x)),p.binary(p.differential(y))
    h=upper.hp.hD(x,y,s);dh=p.binary(p.differential(h))
    t=upper.hp.hD(dx,dy,s)
    Lx,Ly=p.QD(x,s,omega),p.QD(y,s,omega)
    z,v=p.binary(Lx+Ly),p.binary(dh+t)
    K=lambda a,b:p.cup(a,b,a.degree)
    composition=p.binary(K(Lx,Ly)+p.Q(h,1)+K(z,v)+K(dh,t))
    return p.binary(composition+p.cup(s,h)+upper.hE(p.cup(s,x),p.cup(s,y)))


def affine_J0_source(A,Aprime,s,omega):
    """The A-only residue of J0 at B+B'+alpha, on legal B,B'."""
    alpha=upper.hp.hD(p.binary(A),p.binary(Aprime),s)
    P,Pp=upper.primary(A,s,omega),upper.primary(Aprime,s,omega)
    return p.binary(J0(alpha,s,omega)+j0_homotopy(P,Pp,s,omega)
                    +j0_homotopy(p.binary(P+Pp),p.binary(p.differential(alpha)),s,omega))


def affine_J0_primitive(A,B,Aprime,Bprime,s,omega):
    alpha=upper.hp.hD(p.binary(A),p.binary(Aprime),s)
    return p.binary(j0_homotopy(B,Bprime,s,omega)
                    +j0_homotopy(p.binary(B+Bprime),alpha,s,omega))
