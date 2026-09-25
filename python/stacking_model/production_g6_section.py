"""Concrete production operations for the finite general g6 section engine.

Only linear finite-basis conversions are injected. Every nonlinear lower
operation, phase, integral carry, and legal-cylinder prism is selected here
from the fixed formulas, including the common pure-C normalization and K
rephasing. The g6 correction itself retains the original J on legal AB.
"""
from functools import lru_cache

import g6_general_section as engine
import lower_raw_comparison as comparison
import pure_c_normalization as normalized
import closed_ab_upper as closed
from g6_repair import integral_vector,f2_apply,int_apply
from upper_phase_diagnostic import production_phase
p,off,lower=comparison.p,comparison.off,comparison.lower


class ProductionG6Section:
    """An explicit finite-model g6 and current-product Gamma6/Gamma7.

    integral_cochain(q,values), binary_cochain(q,bits), bits(cochain), and
    vector(cochain) are inverse evaluations in a fixed finite cochain basis.
    They must also implement its normalized simplicial pullbacks. These
    linear coordinate maps are the only user-supplied callbacks.
    """
    def __init__(self,coordinates,dim_b,dim_c,delta_b,delta_c,s,omega,
                 integral_cochain,binary_cochain,bits,vector):
        self.s,self.omega=s,omega
        self.integral_cochain,self.binary_cochain=integral_cochain,binary_cochain
        self.bits,self.vector=bits,vector
        self.coordinates=coordinates
        self.delta_b,self.delta_c=tuple(delta_b),tuple(delta_c)
        self.zero=((0,)*coordinates.dimension,0,0)
        self.zero_gamma=integral_vector(vector(p.zero(7)),'zero Gamma')
        self.zero_top=integral_vector(vector(p.zero(8)),'zero T')
        self.model=engine.FiniteCurrentSection(coordinates,dim_b,dim_c,delta_b,delta_c,
                         self.primary,self.current_f,self.coordinate_gauge)
        self.natural=engine.NaturalSectionStacking(self.model,self.alpha,self.natural_beta,
                         self.T,self.J,self.legal_gamma,self.closed_gamma,
                         self.pure_successor,self.transports)
        self.current=engine.CurrentSectionStacking(self.natural,self.current_beta,
                         self.local_pair_gauge,self.C_transport)

    @lru_cache(None)
    def cochains(self,value):
        A,B,C=value
        return self.integral_cochain(3,A),self.binary_cochain(4,B),self.binary_cochain(5,C)

    @lru_cache(None)
    def primary(self,A):
        return self.bits(off.primary(self.integral_cochain(3,A),self.s,self.omega))

    def on_locus(self,A,B):
        return (not any(int_apply(self.coordinates.delta,A))
                and not(f2_apply(self.delta_b,B)^self.primary(A)))

    def lower_pair(self,A,B):
        return off.LowerPair(self.integral_cochain(3,A),self.binary_cochain(4,B),self.on_locus(A,B))

    @lru_cache(None)
    def current_f(self,A,B):
        return self.bits(off.current_f(self.lower_pair(A,B),self.s,self.omega))

    @lru_cache(None)
    def coordinate_gauge(self,A,B):
        return self.bits(off.coordinate_gauge(self.lower_pair(A,B),self.s,self.omega))

    @lru_cache(None)
    def alpha(self,A,Ap):
        if not any(A)or not any(Ap):return 0
        return self.bits(off.alpha(self.integral_cochain(3,A),self.integral_cochain(3,Ap),self.s))

    @lru_cache(None)
    def natural_beta(self,A,B,Ap,Bp):
        if (not any(A)and not B)or(not any(Ap)and not Bp):return 0
        return self.bits(comparison.natural_beta(self.integral_cochain(3,A),self.binary_cochain(4,B),
                         self.integral_cochain(3,Ap),self.binary_cochain(4,Bp),self.s,self.omega))

    @lru_cache(None)
    def current_beta(self,A,B,Ap,Bp):
        if (not any(A)and not B)or(not any(Ap)and not Bp):return 0
        As=tuple(a+b for a,b in zip(A,Ap));Bs=B^Bp^self.alpha(A,Ap)
        return self.bits(lower.all_cochain_beta(self.lower_pair(A,B),self.lower_pair(Ap,Bp),
                         self.on_locus(As,Bs),self.s,self.omega))

    @lru_cache(None)
    def local_pair_gauge(self,A,B,Ap,Bp):
        if (not any(A)and not B)or(not any(Ap)and not Bp):return 0
        return self.bits(comparison.local_extension(self.integral_cochain(3,A),self.binary_cochain(4,B),
                         self.integral_cochain(3,Ap),self.binary_cochain(4,Bp),self.s,self.omega))

    def calibrated_total(self,u,v):
        A,B,C=u;Ap,Bp,Cp=v
        return (tuple(a+b for a,b in zip(A,Ap)),B^Bp^self.alpha(A,Ap),
                C^Cp^self.current_beta(A,B,Ap,Bp))

    @lru_cache(None)
    def K6(self,u):
        A,B,C=u
        if any(A)or B or f2_apply(self.delta_c,C):return self.zero_gamma
        c=self.binary_cochain(5,C)
        return integral_vector(self.vector(p.ds(p.half(p.cup(self.s,c)),self.s)),'K6')

    @lru_cache(None)
    def T(self,u):
        if u==self.zero:return self.zero_top
        if self.model.lower(u)!=self.model.zero_image:
            raise ValueError('production T requires a full legal defining system')
        return integral_vector(self.vector(p.ds(production_phase(3,*self.cochains(u),self.s,self.omega),self.s)),'T6')

    @lru_cache(None)
    def J(self,u):
        if not self.on_locus(u[0],u[1]):
            raise ValueError('production J requires legal AB')
        return integral_vector(self.vector(closed.J(*self.cochains(u),self.s,self.omega)),'J6')

    @lru_cache(None)
    def legal_gamma(self,u,v):
        if u==self.zero or v==self.zero:return self.zero_gamma
        total=self.calibrated_total(u,v)
        base=integral_vector(self.vector(normalized.gamma(*self.cochains(u),*self.cochains(v),self.s,self.omega)),
                             'normalized legal Gamma6')
        return engine._add(base,engine._neg(self.K6(u)),engine._neg(self.K6(v)),self.K6(total))

    @lru_cache(None)
    def closed_gamma(self,u,v):
        if u==self.zero or v==self.zero:return self.zero_gamma
        total=self.calibrated_total(u,v)
        t,tp,ts=(self.binary_cochain(6,self.model.lower(x)[2])for x in(u,v,total))
        carry=p.half(p.cup(self.s,t))+p.half(p.cup(self.s,tp))-p.half(p.cup(self.s,ts))
        base=normalized.closed_gamma(*self.cochains(u),*self.cochains(v),self.s,self.omega)-carry
        return engine._add(integral_vector(self.vector(base),'rephased closed-AB Gamma6'),
                           engine._neg(self.K6(u)),engine._neg(self.K6(v)),self.K6(total))

    @lru_cache(None)
    def pure_successor(self,c,cp):
        from compatible_sector import pure_c_gamma
        from cochains import Cochain
        adapt=lambda x:Cochain(x.degree,x)
        C,Cp=self.binary_cochain(6,c),self.binary_cochain(6,cp)
        base=pure_c_gamma(*map(adapt,(C,Cp,self.s,self.omega)))
        K=lambda z:p.ds(p.half(p.cup(self.s,z)),self.s)
        value=p.Cochain(8,base)-K(C)-K(Cp)+K(p.binary(C+Cp))
        return integral_vector(self.vector(value),'rephased pure-C Gamma7')

    def transports(self,*values):
        result=engine.legal_cylinder_transports(*(self.cochains(v)for v in values),self.s,self.omega)
        return tuple(integral_vector(self.vector(c),'legal-cylinder T prism')for c in result)

    def C_transport(self,value,lam):
        correction=engine.C_shift_transport(self.cochains(value),self.binary_cochain(4,lam),self.s,self.omega)
        return integral_vector(self.vector(correction),'current product C transport')

    def correction(self,u):return self.current.correction(u)
    def gamma(self,u,v):return self.current.gamma(u,v)
    def successor_gamma(self,w,wp):return self.current.successor_gamma(w,wp)
    def stack(self,u,v):return self.current.stack(u,v)
