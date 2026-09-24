"""Current degree-three Danus d5 phase, with its fixed suspension selectors.

Only the final d_s D3 term is omitted: its further differential is zero.
The degree-three source primitive is evaluated by fixed universal chain
operators, never by solving a residual equation on the supplied space.
"""
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from itertools import combinations
import json

try:
    from . import phase_eval as p, low_phases as low, r3_source as rs, r3_chain as rc
    from .mod3_power import tertiary_three_primary_phase
except ImportError:
    import phase_eval as p, low_phases as low, r3_source as rs, r3_chain as rc
    from mod3_power import tertiary_three_primary_phase

SOURCE_DATA = Path(__file__).with_name('r3_source.json')


@lru_cache(None)
def source_coefficients():
    if not SOURCE_DATA.exists():
        raise RuntimeError('fixed universal R3 source table is missing')
    return rs.ef_from_data(json.loads(SOURCE_DATA.read_text()))


def corrections(A,s,omega):
    a=p.binary(A)
    ss=p.cup(s,s)
    return dict(L=p.scale(low.transported(low.pontryagin(omega),A,s),Fraction(1,4)),
                N=p.half(p.cup(ss,p.square(a,2))),
                O=p.half(p.cup(omega,p.square(a,2))),
                M=p.half(p.cup(p.cup(ss,omega),a)))


def build_phase(A,b,c,s,omega,*,coefficients=None,ef=None):
    if (A.degree,b.degree,c.degree,s.degree,omega.degree)!=(3,4,5,1,2):
        raise ValueError('degree-three d5 expects degrees (3,4,5,1,2)')
    if coefficients is None:
        try: from .high_calibration import constants
        except ImportError: from high_calibration import constants
        coefficients=constants()
    if ef is None: ef=source_coefficients()
    data=p.source(A,s,omega)
    g,k=data['g'],data['k']
    y=p.QD(b,s,omega)
    summand=g+p.cup(s,b)
    lam=p.binary(p.divide(summand-p.binary(summand),2,'degree-three splitting carry'))
    dlam=p.binary(p.differential(lam))
    psi=p.binary(y+k+dlam)
    core=p.half(p.binary(p.cup(y,k,5)+p.cup(k,p.binary(p.differential(y)),6)
        +p.Q(k,1)+p.E(lam,omega)+p.cup(p.binary(y+k),dlam,5)))
    bI=p.Cochain(4,lambda z:b(z)*z[-1][-1])
    prism=p.prism_cochain(p.theta(p.binary(p.differential(bI)),p.QD(bI,s,omega),s,omega))
    @lru_cache(None)
    def source_value(pair):
        return (rs.evaluate_r(rs.phi(),rc.Htot(pair))
                +sum(weight*ef.get(basis,0) for basis,weight in rc.Ftot(pair).items()))
    def source_on_vertices(z):
        # V3 is literally relative to A=0, including arbitrary backgrounds.
        if not any(A(face) for face in combinations(z,4)):return Fraction(0)
        registered=low._registered(3,z,A,b,c,s,omega)
        return source_value(rs.to_diags(registered))
    V=p.Cochain(7,source_on_vertices)
    phase=p.half(p.E(c,omega))+core-prism-V+p.half(p.binary(
        p.E(p.Q(b,1),omega)+p.cup(s,psi)))
    terms=corrections(A,s,omega)
    for name,coefficient in [('L',coefficients['c4']),
            ('N',coefficients['cN']+coefficients['epsilon_c']),
            ('O',coefficients['cO']+coefficients['epsilon_c']),
            ('M',coefficients['cM'])]:
        phase=phase+p.scale(terms[name],coefficient)
    return phase+tertiary_three_primary_phase(A,s)
