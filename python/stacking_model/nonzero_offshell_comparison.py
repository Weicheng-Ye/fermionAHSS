"""Exact nonzero-A comparison gauges in physical degrees four and five.

The lower bridge allows signed-closed A and unrestricted B. The upper
potential comparison is calibrated on delta B=P(A), with arbitrary C.
Degree three is not covered by these fixed source contractors.
"""
from functools import lru_cache
from fractions import Fraction as F
from itertools import combinations
import off_shell_beta as off
import closed_a_upper as sharp
hp,p=off.hp,off.p


def r(cochain,s):
    return p.binary(p.Q(cochain,1)+p.cup(s,cochain))


def lower_gauge(A,B,s,omega):
    """L=h+q_loc+P cup_top b, for signed-closed A in degrees 1..3."""
    if A.degree not in (1,2,3):
        raise ValueError("the fixed A-only comparison covers degrees 1..3")
    b=p.binary(p.differential(B)+off.primary(A,s,omega))
    P=off.primary(A,s,omega)
    h=hp.comparison_gauge(A,B,s,omega)
    return p.binary(h+off.local_comparison_extension(A,B,s,omega)
                    +p.cup(P,b,P.degree))


def lower_comparison_residual(A,B,s,omega):
    b=p.binary(p.differential(B)+off.primary(A,s,omega))
    return p.binary(off.natural_f(A,B,s,omega)+sharp.fsharp(A,B,s,omega)
                    +r(b,s)+p.differential(lower_gauge(A,B,s,omega)))


def natural_curvature(A,B,C,s,omega):
    a,b=off.curvature(A,B,s,omega)
    c=p.binary(p.differential(C)+off.natural_f(A,B,s,omega))
    return a,b,c


def raw_potential(A,B,C,s,omega):
    """Phi_hat=I Omega_next(Fhat((A,B,C) ell)); n=1,2 here."""
    if A.degree not in (1,2):
        raise ValueError('the production successor potential is used in n=1,2')
    from upper_phase_diagnostic import production_phase
    I=hp.interval
    a,b,c=natural_curvature(I(A,True),I(B,True),I(C,True),I(s),I(omega))
    return hp.prism(production_phase(A.degree+1,a,b,c,I(s),I(omega)))


@lru_cache(None)
def _binary_comparison(m):
    from higher_theta_comparison import HigherThetaComparison
    return HigherThetaComparison(m)


@lru_cache(None)
def _theta_pair(m):
    from theta_pair_phase import ThetaPairPhase
    return ThetaPairPhase(m)


def _adapt(c):
    from cochains import Cochain
    return Cochain(c.degree,c)


def output_phase_gauge(A,B,C,s,omega):
    """K with Omega_next(Fhat)=−Theta(sharp curvature)+delta_s K, mod1.

    A must be signed closed. B,C are arbitrary. Both existing universal
    helpers are evaluated only on their actual closed/legal domains.
    """
    if A.degree not in (1,2):
        raise ValueError('this production comparison covers n=1,2')
    from production_gamma6_comparison import C_transport
    b,c=sharp.curvature(A,B,C,s,omega)
    rb=r(b,s)
    L=lower_gauge(A,B,s,omega)
    qb=_binary_comparison(b.degree).primitive(*map(_adapt,(b,s,omega)))
    pair=_theta_pair(b.degree).phase(*map(_adapt,(b,c,b,c,s,omega)))
    return (p.Cochain(qb.degree,qb)+p.half(sharp.hE(c,rb))
            +C_transport(p.binary(c+rb),L,omega)-p.Cochain(pair.degree,pair))


def closed_comparison(A,B,C,s,omega):
    """Closed Q/Z phase on signed-closed A and completely unrestricted B,C."""
    return (raw_potential(A,B,C,s,omega)-sharp.phi(A,B,C,s,omega)
            -output_phase_gauge(A,B,C,s,omega))


def A_source(A,s,omega):
    """Explicit closed A-only source before the degree-specific calibration."""
    return closed_comparison(A,p.zero(A.degree+1),p.zero(A.degree+2),s,omega)


def BC_comparison_gauge(A,B,C,s,omega):
    """delta q=closed_comparison(A,B,C)−S_A, modulo integers."""
    I=hp.interval
    return hp.prism(closed_comparison(I(A),I(B,True),I(C,True),I(s),I(omega)))


def degree_five_raw_source(A,s,omega):
    """The original eight-cell k5 source, before the matched transgressions."""
    if A.degree!=2:
        raise ValueError('degree-five calibration expects integral A2')
    from production_gamma5_comparison import A_phase
    k0=hp.uniform_source(A,s,omega)['k0']
    return A_phase(A,s,omega)-A_source(A,s,omega)-p.half(p.cup(s,k0))


def signed_cube(A,s):
    return hp.low.transported(hp.low.transported(A,A,s),A,s)


def degree_five_final_source(A,s,omega):
    """Corrected exact source, including the two B transgressions and cube."""
    P=off.primary(A,s,omega)
    return (degree_five_raw_source(A,s,omega)-p.half(p.cup(p.binary(A),P))
            -p.half(p.cup(p.cup(s,s),P))+p.scale(signed_cube(A,s),F(8,3)))


@lru_cache(None)
def primary_suspension_comparison():
    """Fixed V with delta V=h[Q1 QD(b)+E(s b)] for closed binary b3."""
    from fractions import Fraction as F
    from cochains import Cochain
    from higher_theta_comparison import HigherThetaComparison

    class PrimarySuspensionComparison(HigherThetaComparison):
        def source(self,b,s,omega):
            op=p.half(p.binary(p.Q(p.QD(b,s,omega),1)+p.E(p.cup(s,b),omega)))
            return Cochain(op.degree,op)

        @lru_cache(None)
        def certificate(self):
            lower,upper=self.small_basis(5),self.small_basis(6)
            rows=tuple(tuple(self.small_boundary(b).get(a,0)for a in lower)for b in upper)
            values=tuple(self.chain_value(self.inclusion(b))for b in upper)
            coefficients=(F(1,4),F(0),F(0),F(1,4),F(0))
            if len(lower)!=5 or len(upper)!=11:
                raise ArithmeticError('the primary suspension source basis changed')
            if any((sum(a*b for a,b in zip(row,coefficients))-v)%1
                   for row,v in zip(rows,values)):
                raise ArithmeticError('the primary suspension certificate failed')
            return dict(lower=lower,upper=upper,rows=rows,values=values,
                        coefficients=coefficients,denominator=4)

    return PrimarySuspensionComparison(3)


def degree_four_final_source(A,s,omega):
    """The remaining k4 A-only source after the h(s t) rephasing."""
    if A.degree!=1:
        raise ValueError('degree-four calibration expects integral A1')
    from production_gamma4_comparison import A_phase
    P=off.primary(A,s,omega)
    V=primary_suspension_comparison().primitive(*map(_adapt,(P,s,omega)))
    k0=hp.uniform_source(A,s,omega)['k0']
    return (A_phase(A,s,omega)-A_source(A,s,omega)+p.Cochain(V.degree,V)
            -p.half(p.cup(s,k0)))


K4_SMALL_COEFFICIENTS={(1,(2,)):F(0),(2,(1,)):F(1,4),(4,()):F(1,4)}


@lru_cache(None)
def _k4_source_value(pair):
    A,s,omega=hp.from_pair(pair,1)
    if all(not A((i,i+1))for i in range(5)):
        return 0
    # The group contractor uses a homogeneous global frame; production
    # cochains return values in the first vertex's local frame.
    return (-1)**pair[0][0][1]*degree_four_final_source(A,s,omega)(tuple(range(6)))


def _k4_small_value(group,word):
    if hp.low._gdegen(group) or len(group)<2:
        return 0
    edges=tuple(((-1)**a[1]*(b[0]-a[0]),a[1]^b[1])
                for a,b in zip(group,group[1:]))
    if all(a==0 for a,s in edges):
        return 0
    if any(edge!=(1,1)for edge in edges):
        raise ArithmeticError('the retracted source is outside the dihedral wedge')
    return (-1)**group[0][1]*K4_SMALL_COEFFICIENTS.get((len(group)-1,word),0)


@lru_cache(None)
def _k4_primitive_value(pair):
    large=sum(c*_k4_source_value(x)for x,c in hp.group_homotopy(pair).items())
    small=0
    for (g,w),coefficient in hp.low.group_product_AW(pair).items():
        for group,c in hp.low.retract_group(g).items():
            for word,d in hp.cm.F(w).items():
                small+=coefficient*c*d*_k4_small_value(group,word)
    return large+small


def degree_four_calibration_primitive(A,s,omega):
    """q_A=S H+e F for the fully evaluated five-cell k4 A-only source."""
    if A.degree!=1:
        raise ValueError('degree-four calibration expects integral A1')
    return p.Cochain(4,lambda vertices:
                     _k4_primitive_value(hp.low.to_pair1(vertices,A,s,omega)))


def primary_suspension_gauge(B,s,omega):
    """q_U with U(B)=V(delta B)+delta_s q_U, for unrestricted B2."""
    I=hp.interval
    Bi,si,wi=I(B,True),I(s),I(omega)
    U=p.half(p.binary(p.Q(p.QD(Bi,si,wi),1)+p.E(p.cup(si,Bi),wi)))
    V=primary_suspension_comparison().primitive(*map(_adapt,
                        (p.binary(p.differential(Bi)),si,wi)))
    return hp.prism(U-p.Cochain(V.degree,V))


def legal_pair_potential(A,B,C,s,omega):
    from upper_phase_diagnostic import production_phase
    tau=hp.low.secondary(A,B,s,omega)
    t=p.binary(p.differential(C)+tau)
    return production_phase(A.degree,A,B,C,s,omega)+p.half(p.cup(t,tau,A.degree+2))


def legal_pair_gauge4(A,B,C,s,omega):
    """q4 with Phi_J−Phi_hat=h(s t)+delta_s q4 modulo one.

    Domain: A1 signed closed and delta B2=P(A); C3 is arbitrary. This
    compares the actual prescribed J4 with the natural suspended model.
    """
    if A.degree!=1:
        raise ValueError('this calibration gauge is for physical degree four')
    from production_gamma4_comparison import comparison_gauge
    I=hp.interval
    def difference(a,b,c,ss,ww):
        tau=hp.low.secondary(a,b,ss,ww)
        t=p.binary(p.differential(c)+tau)
        return legal_pair_potential(a,b,c,ss,ww)-raw_potential(a,b,c,ss,ww)-p.half(p.cup(ss,t))
    cpart=hp.prism(difference(I(A),I(B),I(C,True),I(s),I(omega)))
    lam=hp.splitting_carry(A,B,s,omega)
    return (cpart+degree_four_calibration_primitive(A,s,omega)
            +primary_suspension_gauge(B,s,omega)-p.half(p.cup(s,lam))
            +comparison_gauge(A,B,p.zero(3),s,omega)
            -BC_comparison_gauge(A,B,p.zero(3),s,omega))


K5_SMALL_COEFFICIENTS=dict(zip(hp.basis2(5),(F(0),F(3,8),F(1,8),F(0))))


@lru_cache(None)
def _k5_source_value(pair):
    A,s,omega=hp.from_pair(pair,2)
    if all(not A(face)for face in combinations(range(7),3)):
        return F(0)
    return degree_five_final_source(A,s,omega)(tuple(range(7)))


@lru_cache(None)
def _k5_primitive_value(pair):
    return (sum(c*_k5_source_value(x)for x,c in hp.sp.Htot(pair).items())
            +sum(c*K5_SMALL_COEFFICIENTS.get(key,F(0))for key,c in hp.sp.Ftot(pair).items()))


def degree_five_calibration_primitive(A,s,omega):
    """q_A=S H+e F with the complete corrected eight-cell k5 certificate."""
    if A.degree!=2:
        raise ValueError('degree-five calibration expects integral A2')
    return p.Cochain(5,lambda vertices:_k5_primitive_value(hp.to_pair(A,s,omega,vertices)))


def legal_pair_gauge5(A,B,C,s,omega):
    """q5 with Phi_J−Phi_hat=h(s t)−8A^3/3+delta_s q5, modulo one.

    A2 is signed closed, delta B3=P(A), and C4 is arbitrary. The signed
    rational cubic term is exactly closed and has zero effect on g.
    """
    if A.degree!=2:
        raise ValueError('this calibration gauge is for physical degree five')
    from production_gamma5_comparison import comparison_gauge
    I=hp.interval
    def difference(a,b,c,ss,ww):
        tau=hp.low.secondary(a,b,ss,ww)
        t=p.binary(p.differential(c)+tau)
        return legal_pair_potential(a,b,c,ss,ww)-raw_potential(a,b,c,ss,ww)-p.half(p.cup(ss,t))
    cpart=hp.prism(difference(I(A),I(B),I(C,True),I(s),I(omega)))
    lam=hp.splitting_carry(A,B,s,omega)
    return (cpart+degree_five_calibration_primitive(A,s,omega)
            +p.half(p.cup(p.binary(A),B))+p.half(p.cup(p.cup(s,s),B))
            -p.half(p.cup(s,lam))+comparison_gauge(A,B,s,omega)
            -BC_comparison_gauge(A,B,p.zero(4),s,omega))


def integral_C_transport(c,L,s,omega):
    """Exact integer carry for shifting an arbitrary C cochain by delta L.

    Production Omega's other terms cancel. No universal A phase is
    evaluated in this operation.
    """
    from production_gamma6_comparison import C_transport
    shifted=p.binary(c+p.differential(L))
    raw=p.half(p.E(shifted,omega))-p.half(p.E(c,omega))-p.ds(C_transport(c,L,omega),s)
    def value(vertices):
        answer=F(raw(vertices))
        if answer.denominator!=1:
            raise ArithmeticError('the C-transport carry is not integral')
        return answer.numerator
    return p.Cochain(c.degree+2,value)


def off_branch_integer_gauge(A,B,C,s,omega):
    """Z_off with G_current=ghat(A,B,C+q_loc)+delta_s Z_off."""
    I=hp.interval
    ai,bi,ci,si,wi=I(A,True),I(B,True),I(C,True),I(s),I(omega)
    q=off.local_comparison_extension(A,B,s,omega)
    M=p.binary(I(q,True)+off.local_comparison_extension(ai,bi,si,wi))
    c_raw=p.binary(p.differential(ci)+off.raw_H(ai,bi,si,wi))
    return p.scale(hp.prism(integral_C_transport(c_raw,M,si,wi)),-1)


def natural_input_C_shift(A,B,C,kappa,s,omega):
    """Integer T with ghat(A,B,C+delta kappa)−ghat(A,B,C)=delta_s T."""
    I=hp.interval
    si,wi=I(s),I(omega)
    _,_,c=natural_curvature(I(A,True),I(B,True),I(C,True),si,wi)
    dk=p.binary(p.differential(kappa))
    along=integral_C_transport(c,I(dk,True),si,wi)
    return (hp.prism(along)+p.half(p.E(dk,omega))
            -p.ds(p.half(p.E(kappa,omega)),s))
