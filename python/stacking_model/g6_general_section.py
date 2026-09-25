"""Finite section and two-cylinder g6 construction without associativity.

Finite cochains use (integral A tuple, binary B bits, binary C bits).
Injected maps are evaluations of the explicit cochain formulas already
implemented in this scratch package. The section solves only lower
defining equations. Its exhaustive search is finite, but exponential.
See G6_GENERAL_SECTION.md for the contracts and exact transport signs.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import product

from g6_repair import (f2_apply, f2_solve, f2_kernel, xor_columns,
                       integral_vector, int_apply)


def _add(*vectors):
    if not vectors:
        return ()
    if any(len(v)!=len(vectors[0])for v in vectors):
        raise ValueError('incompatible integral cochain dimensions')
    return tuple(sum(entries)for entries in zip(*vectors))


def _neg(vector):
    return tuple(-x for x in vector)


def _injective_solve(matrix,target,ncols):
    """Exact solution in an injective integer column complement."""
    if len(matrix)!=len(target) or any(len(row)!=ncols for row in matrix):
        raise ValueError('incompatible integral equation dimensions')
    rows=[list(map(Fraction,row))+[Fraction(value)]for row,value in zip(matrix,target)]
    for col in range(ncols):
        pivot=next((i for i in range(col,len(rows))if rows[i][col]),None)
        if pivot is None:
            raise ValueError('the declared integral complement is not injective')
        rows[col],rows[pivot]=rows[pivot],rows[col]
        leading=rows[col][col]
        rows[col]=[x/leading for x in rows[col]]
        for i in range(len(rows)):
            if i!=col and rows[i][col]:
                scale=rows[i][col]
                rows[i]=[a-scale*b for a,b in zip(rows[i],rows[col])]
    if any(not any(row[:ncols])and row[-1]for row in rows):
        return None
    solution=tuple(rows[i][-1]for i in range(ncols))
    if any(x.denominator!=1 for x in solution):
        return None
    return tuple(x.numerator for x in solution)


class FiniteCurrentSection:
    """Explicit section of F_current, conjugated to natural C coordinates.

    coordinates is IntegralCycleCoordinates for delta_s on integral A3.
    primary(A) and current_f(A,B) are literal P and globally branched f6.
    gauge(A,B) is the current-to-natural C shift, zero on the legal AB locus.
    All three callbacks return binary coordinate vectors.
    """
    def __init__(self,coordinates,dim_b,dim_c,delta_b,delta_c,
                 primary,current_f,gauge):
        self.coordinates=coordinates
        self.dim_b,self.dim_c=dim_b,dim_c
        self.delta_b,self.delta_c=tuple(delta_b),tuple(delta_c)
        self.primary,self.current_f,self.gauge=primary,current_f,gauge
        if len(self.delta_b)!=dim_c:
            raise ValueError('delta B must target the C-coordinate space')
        for rows,dimension in ((self.delta_b,dim_b),(self.delta_c,dim_c)):
            if any(not isinstance(row,int)or row<0 or row>>dimension for row in rows):
                raise ValueError('binary differential has incompatible dimensions')
        if any(f2_apply(self.delta_c,f2_apply(self.delta_b,1<<j))for j in range(dim_b)):
            raise ValueError('binary differentials do not square to zero')
        self.closed_b=f2_kernel(self.delta_b,dim_b)
        self.zero=((0,)*coordinates.dimension,0,0)
        self.zero_image=((0,)*len(coordinates.delta),0,0)
        if self.lower(self.zero)!=self.zero_image or gauge(self.zero[0],0):
            raise ValueError('the prescribed lower operations must be pointed')

    @staticmethod
    def _bits(value,dimension,name):
        if not isinstance(value,int)or value<0 or value>>dimension:
            raise ValueError(f'{name} is outside its binary coordinate space')
        return value

    def validate(self,value):
        A,B,C=value
        A=integral_vector(A,'A')
        if len(A)!=self.coordinates.dimension:
            raise ValueError('A has incompatible integral coordinates')
        return A,self._bits(B,self.dim_b,'B'),self._bits(C,self.dim_c,'C')

    def lower(self,value):
        A,B,C=self.validate(value)
        P=self._bits(self.primary(A),self.dim_c,'P')
        f=self._bits(self.current_f(A,B),len(self.delta_c),'f')
        return (int_apply(self.coordinates.delta,A),f2_apply(self.delta_b,B)^P,
                f2_apply(self.delta_c,C)^f)

    def conjugate(self,value):
        """The same involution converts current to natural and back."""
        A,B,C=self.validate(value)
        e=self._bits(self.gauge(A,B),self.dim_c,'C-coordinate gauge')
        return A,B,C^e

    def natural_lower(self,value):
        return self.lower(self.conjugate(value))

    @lru_cache(None)
    def current_section(self,image):
        a,b,c=image
        a=integral_vector(a,'output A')
        b=self._bits(b,self.dim_c,'output B')
        c=self._bits(c,len(self.delta_c),'output C')
        coordinates=self.coordinates
        rank=coordinates.cycle_rank
        complement=tuple(tuple(sum(row[j]*coordinates.basis[j][col]
                              for j in range(coordinates.dimension))
                              for col in range(rank,coordinates.dimension))
                         for row in coordinates.delta)
        fixed=_injective_solve(complement,a,coordinates.dimension-rank)
        if fixed is None:
            raise ValueError('the integral output is not a coboundary')
        # Tau and raw H depend on A only modulo four. The global legality
        # predicate is unchanged by adding four integral cycles. No such
        # periodicity is assumed for the comparison cochain gauge.
        for residues in product(range(4),repeat=rank):
            A=int_apply(coordinates.basis,residues+fixed)
            rhs=b^self._bits(self.primary(A),self.dim_c,'P')
            B0=f2_solve(self.delta_b,rhs,self.dim_b)
            if B0 is None:
                continue
            for bits in range(1<<len(self.closed_b)):
                B=B0^xor_columns(self.closed_b,bits)
                f=self._bits(self.current_f(A,B),len(self.delta_c),'f')
                C=f2_solve(self.delta_c,c^f,self.dim_c)
                if C is not None:
                    answer=A,B,C
                    if self.lower(answer)!=(a,b,c):
                        raise ArithmeticError('the section failed its defining equations')
                    return answer
        raise ValueError('the triple is not in the image of the lower differential')

    def section(self,image):
        return self.conjugate(self.current_section(image))


def natural_product(left,right,s,omega):
    """The explicit natural lower product, including the local lambda3."""
    import lower_raw_comparison as comparison
    p=comparison.p
    A,B,C=left;Ap,Bp,Cp=right
    alpha=comparison.off.alpha(A,Ap,s)
    beta=comparison.natural_beta(A,B,Ap,Bp,s,omega)
    return A+Ap,p.binary(B+Bp+alpha),p.binary(C+Cp+beta)


def natural_left_division(left,total,s,omega):
    """Unique triangular solution right of natural_product(left,right)=total."""
    import lower_raw_comparison as comparison
    p=comparison.p
    A,B,C=left;At,Bt,Ct=total
    Ap=At-A
    Bp=p.binary(Bt+B+comparison.off.alpha(A,Ap,s))
    Cp=p.binary(Ct+C+comparison.natural_beta(A,B,Ap,Bp,s,omega))
    return Ap,Bp,Cp


def legal_cylinder_transports(reference,referencep,total_reference,z,zp,kappa,r,s,omega):
    """Two explicit integral prisms; no associator or unknown primitive."""
    import lower_raw_comparison as comparison
    from upper_phase_diagnostic import production_phase
    p,hp=comparison.p,comparison.hp
    I=hp.interval;si,wi=I(s),I(omega)
    const=lambda triple:tuple(I(c)for c in triple)
    ramp=lambda triple:tuple(I(c,True)for c in triple)
    mul=lambda u,v:natural_product(u,v,si,wi)
    div=lambda u,v:natural_left_division(u,v,si,wi)
    at,bt=ramp(reference),ramp(referencep)
    path1=div(mul(at,bt),mul(mul(at,const(z)),mul(bt,const(zp))))
    ct=ramp(total_reference)
    path2=div(ct,mul(mul(ct,const(kappa)),const(r)))
    def transport(path):
        T=p.ds(production_phase(3,*path,si,wi),si)
        def integer(vertices):
            value=Fraction(T(vertices))
            if value.denominator!=1:
                raise ArithmeticError('the legal cylinder T is not integral')
            return value.numerator
        return hp.prism(p.Cochain(8,integer))
    return transport(path1),transport(path2)


def C_shift_transport(value,lam,s,omega):
    """Integral L with g(C+delta lam)-g(C)=delta_s L on legal AB.

    C may be arbitrary for J; on a full defining system this is the same
    transport for T. Off the legal AB locus, apply this to the legal
    retracted triple instead. The lower section is independent of C shifts
    by a coboundary, so its retract undergoes exactly the same C shift.
    """
    import closed_a_upper as upper
    p=upper.p
    A,B,C=value;shift=p.binary(p.differential(lam))
    gauge=p.half(p.binary(p.E(lam,omega)+upper.hE(C,shift)))
    # The production R3 is independent of C; its exact lift cancels.
    raw=p.half(p.E(p.binary(C+shift),omega)-p.E(C,omega))-p.ds(gauge,s)
    def integer(vertices):
        result=Fraction(raw(vertices))
        if result.denominator!=1:
            raise ArithmeticError('the C-coboundary transport is not integral')
        return result.numerator
    return p.Cochain(7,integer)


class NaturalSectionStacking:
    """Finite-coordinate section identity for the explicit natural lower loop.

    alpha,beta are evaluations of off.alpha and comparison.natural_beta.
    T and J, legal_gamma, closed_gamma, pure_successor are the proved
    production formulas. transports evaluates legal_cylinder_transports
    in the finite basis. These are fixed operations, not primitives sought
    by this class. It does not assume or test associativity.
    """
    def __init__(self,section,alpha,beta,T,J,legal_gamma,closed_gamma,
                 pure_successor,transports):
        self.model=section;self.alpha,self.beta=alpha,beta
        self.T,self.J=T,J
        self.legal_gamma,self.closed_gamma=legal_gamma,closed_gamma
        self.pure_successor,self.transports=pure_successor,transports
        self.zero=section.zero;self.zero_image=section.zero_image
        self.zero_gamma=integral_vector(legal_gamma(self.zero,self.zero),'legal gamma')
        self.zero_top=integral_vector(T(self.zero),'T')
        if any(self.zero_gamma)or any(self.zero_top):
            raise ValueError('upper operations must be pointed')

    def lower(self,value):
        return self.model.natural_lower(value)

    def stack(self,left,right):
        A,B,C=self.model.validate(left);Ap,Bp,Cp=self.model.validate(right)
        if not any(A)and not B:
            return Ap,Bp,C^Cp
        if not any(Ap)and not Bp:
            return A,B,C^Cp
        return self.model.validate((_add(A,Ap),B^Bp^self.alpha(A,Ap),
                                    C^Cp^self.beta(A,B,Ap,Bp)))

    def divide(self,left,total):
        A,B,C=self.model.validate(left);At,Bt,Ct=self.model.validate(total)
        if A==At and B==Bt:
            return self.zero[0],0,C^Ct
        Ap=_add(At,_neg(A));Bp=Bt^B^self.alpha(A,Ap)
        answer=self.model.validate((Ap,Bp,Ct^C^self.beta(A,B,Ap,Bp)))
        if self.stack(left,answer)!=total:
            raise ArithmeticError('triangular left division failed')
        return answer

    def retract(self,value):
        z=self.divide(self.model.section(self.lower(value)),value)
        if self.lower(z)!=self.zero_image:
            raise ArithmeticError('the lower homomorphism failed in section division')
        return z

    @staticmethod
    def pure_image(image):
        return not any(image[0])and not image[1]

    def j(self,image):
        return (integral_vector(self.J(self.model.section(image)),'J')
                if self.pure_image(image)else self.zero_top)

    def q(self,value):
        image=self.lower(value)
        return (integral_vector(self.closed_gamma(self.model.section(image),self.retract(value)),
                                'closed-AB section gamma')
                if self.pure_image(image)else self.zero_gamma)

    def correction(self,value):
        if self.pure_image(self.lower(value)):
            return integral_vector(self.J(value),'J')
        return integral_vector(self.T(self.retract(value)),'T')

    def pair_data(self,left,right):
        w,wp=self.lower(left),self.lower(right)
        a,b=self.model.section(w),self.model.section(wp)
        z,zp=self.retract(left),self.retract(right)
        q=self.stack(a,b);total=self.stack(left,right);wt=self.lower(q)
        if self.lower(total)!=wt:
            raise ArithmeticError('the explicit lower product is not compatible with F')
        reference=self.model.section(wt)
        kappa=self.divide(reference,q);r=self.divide(q,total)
        if self.lower(kappa)!=self.zero_image or self.lower(r)!=self.zero_image:
            raise ArithmeticError('the divided comparison triples are not legal')
        return w,wp,wt,a,b,reference,z,zp,kappa,r,total

    def gamma(self,left,right):
        if self.pure_image(self.lower(left))and self.pure_image(self.lower(right)):
            return integral_vector(self.closed_gamma(left,right),'closed-AB gamma')
        _,_,_,a,b,reference,z,zp,kappa,r,total=self.pair_data(left,right)
        if z==self.zero and zp==self.zero:
            # Both inputs are section representatives: r=0, the first
            # cylinder is zero and the second is the constant kappa.
            # The two legal Gamma terms are normalized unit values.
            if r!=self.zero:
                raise ArithmeticError('reference pair division failed')
            return integral_vector(_add(_neg(self.q(left)),_neg(self.q(right)),self.q(total)),
                                   'reference-pair section gamma6')
        first,second=self.transports(a,b,reference,z,zp,kappa,r)
        return integral_vector(_add(self.legal_gamma(z,zp),self.legal_gamma(kappa,r),
                                    _neg(first),_neg(second),_neg(self.q(left)),
                                    _neg(self.q(right)),self.q(total)),'section gamma6')

    def successor_gamma(self,w,wp):
        if self.pure_image(w)and self.pure_image(wp):
            return integral_vector(self.pure_successor(w[2],wp[2]),'pure-C gamma7')
        q=self.stack(self.model.section(w),self.model.section(wp))
        wt=self.lower(q);kappa=self.retract(q)
        return integral_vector(_add(self.T(kappa),self.j(wt),_neg(self.j(w)),_neg(self.j(wp))),
                               'image gamma7')


class CurrentSectionStacking:
    """Transport the natural section formula to the existing lower product.

    current_beta is the existing globally branched beta6. local_pair_gauge
    evaluates comparison.local_extension in degree-four bits. transport
    evaluates C_shift_transport on the specified legal-AB triple. Thus all
    callback equations are supplied above or in the existing modules.
    """
    def __init__(self,natural,current_beta,local_pair_gauge,transport):
        self.natural=natural;self.model=natural.model
        self.current_beta=current_beta
        self.local_pair_gauge=local_pair_gauge;self.transport=transport

    def stack(self,left,right):
        A,B,C=self.model.validate(left);Ap,Bp,Cp=self.model.validate(right)
        return self.model.validate((_add(A,Ap),B^Bp^self.natural.alpha(A,Ap),
                                    C^Cp^self.current_beta(A,B,Ap,Bp)))

    def correction(self,value):
        return self.natural.correction(self.model.conjugate(value))

    def gamma(self,left,right):
        total=self.stack(left,right)
        lnat,rnat=self.model.conjugate(left),self.model.conjugate(right)
        n=self.natural
        if n.pure_image(self.model.lower(left))and n.pure_image(self.model.lower(right)):
            lam=0
        else:
            lam=self.local_pair_gauge(left[0],left[1],right[0],right[1])
        lam=self.model._bits(lam,self.model.dim_b,'local pair gauge')
        shifted=total[0],total[1],total[2]^f2_apply(self.model.delta_b,lam)
        if self.model.conjugate(n.stack(lnat,rnat))!=shifted:
            raise ArithmeticError('the prescribed current/natural product comparison failed')
        base=(total if n.pure_image(self.model.lower(total))
              else n.retract(self.model.conjugate(total)))
        correction=self.transport(base,lam) if lam else n.zero_gamma
        return integral_vector(_add(n.gamma(lnat,rnat),correction),'current section gamma6')

    def successor_gamma(self,w,wp):
        return self.natural.successor_gamma(w,wp)
