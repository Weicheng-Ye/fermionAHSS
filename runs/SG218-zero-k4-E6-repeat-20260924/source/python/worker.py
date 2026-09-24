#!/usr/bin/env python3
"""Exact batch protocol used by GAP's natural d5 adapter (stdin/stdout)."""
from fractions import Fraction
from itertools import combinations
import json
import os
import sys

if sys.version_info < (3, 10):
    raise SystemExit('natural d5 requires Python 3.10 or newer')

try:
    from .runtime_cache import (install_cochain_cache_limit, cochain_cache_limit,
                                install_chain_cache_limit, chain_cache_limit)
except ImportError:
    from runtime_cache import (install_cochain_cache_limit, cochain_cache_limit,
                               install_chain_cache_limit, chain_cache_limit)

# Set the evaluation policy before phase_eval constructs its input cochains.
# A zero capacity disables memoization; all arithmetic remains exact.
COCHAIN_CACHE_ENTRIES = install_cochain_cache_limit(
    int(os.environ.get('KOAHSS_COCHAIN_CACHE_ENTRIES', '256')))
CHAIN_CACHE_ENTRIES = install_chain_cache_limit(
    int(os.environ.get('KOAHSS_CHAIN_CACHE_ENTRIES', '256')))

try:
    from . import phase_eval as pe
except ImportError:
    import phase_eval as pe


def current_psi(A, b, s, omega):
    """Independent direct calibrated lower formula, including degree zero."""
    n=A.degree
    a=pe.binary(A)
    t=pe.binary(pe.divide(A-a,2,'input carry'))
    B=pe.divide(pe.differential(a),2,'binary input Bockstein')
    e=pe.binary(B)
    carry=pe.binary(pe.divide(B+e,2,'plus carry'))
    u,v,w=pe.square(a,2),pe.cup(omega,a),pe.cup(s,e)
    z2=pe.zeta2(omega,a) if n else pe.zero(3)
    z1=pe.zeta1(s,a) if n else pe.zero(2)
    F=pe.binary(pe.E(b,omega)+z2+pe.chi(a)
                +pe.cup(u,v,n+1)+pe.cup(u,w,n+1)+pe.cup(v,w,n+1)
                +pe.zeta1(s,e)+pe.cup(pe.cup(omega,s,1),e)
                +pe.cup(s,u)+pe.cup(pe.cup(s,s),carry))
    G=pe.binary(pe.Q(t,2)+pe.cup(omega,t)
                +pe.cup(pe.binary(pe.differential(t)),pe.cup(s,a),n)
                +z1+pe.cup(pe.cup(omega,s,1),a)+pe.cup(s,b))
    q=pe.cup(omega,B,integral=True)+pe.cup(B,B,n-1,integral=True)
    L=pe.divide(q-pe.ds(G,s),2,'secondary L')
    psi=pe.binary(F+pe.binary(L)+pe.cup(pe.cup(pe.cup(s,s),s),a))
    return pe.binary(u+v+w),psi


def import_sample(n, sample):
    if not isinstance(sample,dict) or set(sample)!={'s','A','omega','b','c'}:
        raise ValueError('each d5 sample must contain exactly s,A,omega,b,c')
    dim=n+4
    tables={}
    for name,degree in [('s',1),('A',n),('omega',2),('b',n+1),('c',n+2)]:
        faces=tuple(combinations(range(dim+1),degree+1))
        values=sample[name]
        if not isinstance(values,list) or len(values)!=len(faces) or any(type(v) is not int for v in values):
            raise ValueError('invalid integral face data for '+name)
        if name!='A' and any(v not in (0,1) for v in values):
            raise ValueError('nonbinary '+name+' cochain')
        tables[name]=dict(zip(faces,values))
    return pe.make_source(n,dim,*(tables[name].__getitem__ for name in ['s','A','omega','b','c']))


def check_on_faces(cochain, vertices, modulus, name):
    for indices in combinations(range(len(vertices)),cochain.degree+1):
        value=cochain(tuple(vertices[i] for i in indices))
        if (value % modulus if modulus else value):
            raise ArithmeticError('natural d5 '+name+' failed on face '+str(indices))


def evaluate(request):
    if not isinstance(request,dict) or request.get('schema')!=1:
        raise ValueError('unsupported natural d5 protocol')
    n=request.get('degree')
    if type(n) is not int or n not in range(4):
        raise ValueError('natural d5 supports only degrees 0,1,2,3')
    samples=request.get('samples')
    if not isinstance(samples,list):
        raise ValueError('natural d5 samples must be a list')
    A=pe.input_cochain('A',n)
    b=pe.input_cochain('b',n+1)
    c=pe.input_cochain('c',n+2)
    s,omega=pe.s,pe.omega
    primary,psi=current_psi(A,b,s,omega)
    identities=[(pe.differential(s),2,'ds=0'),
                (pe.differential(omega),2,'domega=0'),
                (pe.ds(A,s),0,'d_s A=0'),
                (pe.differential(b)-primary,2,'db=Da'),
                (pe.differential(c)-psi,2,'dc=psi')]
    phase=None
    output=[]
    for sample in samples:
        vertices=import_sample(n,sample)
        for cochain,modulus,name in identities:
            check_on_faces(cochain,vertices,modulus,name)
        if not any(sample['A']) and not any(sample['b']) and not any(sample['c']):
            value=Fraction(0) # Every displayed formula term vanishes literally.
        else:
            if phase is None:
                if n<3:
                    try:
                        from .low_phases import build_phase
                    except ImportError:
                        from low_phases import build_phase
                    phase=build_phase(n,A,b,c,s,omega)
                else:
                    try:
                        from .high_phase import build_phase
                    except ImportError:
                        from high_phase import build_phase
                    phase=build_phase(A,b,c,s,omega)
            value=Fraction(phase(vertices))
        output.append([value.numerator,value.denominator])
    return dict(status='computed',schema=1,degree=n,phases=output,
                audit=dict(reference='Danus chi7_tail epsilon100 eta101 mu0',
                           defining_systems_checked=len(samples),
                           cochain_cache_entries=cochain_cache_limit(),
                           chain_cache_entries=chain_cache_limit(),
                           arithmetic='exact integer and rational',
                           local_R_solve=False,old_correction_applied=False))


if __name__=='__main__':
    try:
        result=evaluate(json.load(sys.stdin))
    except Exception as error:
        print(json.dumps(dict(status='error',type=type(error).__name__,reason=str(error))))
        sys.exit(1)
    print(json.dumps(result,separators=(',',':')))
