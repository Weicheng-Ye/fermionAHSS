"""Finite R2-to-R3 suspension selectors from r3_suspension.md, Section 3.

The historical V2fin occurs here, before the current lower cubic correction.
This command evaluates the actual comparison maps and records their periods.
"""
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
import hashlib
import json

try:
    from . import phase_eval as p, low_phases as low, low_calibration as lc
    from . import source_primitive as sp, r3_source as rs, r3_chain as rc, chain_models as cm
    from .high_phase import source_coefficients, SOURCE_DATA
except ImportError:
    import phase_eval as p, low_phases as low, low_calibration as lc
    import source_primitive as sp, r3_source as rs, r3_chain as rc, chain_models as cm
    from high_phase import source_coefficients, SOURCE_DATA

DATA=Path(__file__).with_name('high_calibration.json')


def provenance():
    folder=Path(__file__).resolve().parent
    paths=[folder/name for name in ('high_calibration.py','high_phase.py',
           'low_phases.py','low_calibration.json','source_primitive.py')]
    paths.append(SOURCE_DATA)
    return {f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in paths}


def integer(value,label):
    value=Fraction(value)
    if value.denominator!=1: raise ArithmeticError(f'{label} is not integral: {value}')
    return int(value)


def suspended_v3(vertices,ef):
    result=Fraction(0)
    for z,sign in p.prism(vertices):
        def av(f):
            i,j,k,l=f
            return p.A2((z[i],z[j],z[k]))*(z[l][-1]-z[k][-1])
        source=p.make_source(3,len(z)-1,
            lambda f:p.s(tuple(z[i] for i in f)),av,
            lambda f:p.omega(tuple(z[i] for i in f)))
        result+=sign*rs.V3(source,ef) # kappa7=+I
    return result


CYCLES={'F':{((0,3),()):1},'Y':{((0,1),(1,1)):1},
        'W':{((2,2),()):1},'V':{((0,2),(1,)):1},
        'ZT':{((1,1),(2,)):1,((2,1),(1,)):1}}
PREPARED=Path(__file__).with_name('high_calibration_periods')


def preparation_hashes():
    folder=Path(__file__).resolve().parent
    names=('high_calibration.py','r3_source.py','r3_chain.py','phase_eval.py',
           'chain_models.py','exterior_bar.py','cochain_tools.py',
           'source_primitive.py','low_phases.py','low_calibration.json')
    out={name:hashlib.sha256((folder/name).read_bytes()).hexdigest() for name in names}
    path=p.PACKAGE/'data/chi-calibrated-degree7-anf.g'
    out[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    return out


def _save_checkpoint(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix('.tmp')
    temp.write_text(json.dumps(data,indent=2)+'\n')
    temp.replace(path)


def _clear_period_sources():
    # These are the two registration caches used in this isolated driver.
    # Their values contain source ids, so clear both before recycling ids.
    p.from_diags.cache_clear()
    rs.clear_source_memos()


def normalized_lower(cycle):
    out={}
    for b,c in cycle.items():
        for pair,d in sp.Gtot(b).items():
            # Cochain classification factors through j2 r2 strictly.
            projected=p.to_diags(p.from_diags(pair))
            cm.add(out,cm._basis(projected),c*d)
    return out


def normalized_suspension(lower):
    out={}
    for pair,coefficient in lower.items():
        vertices=p.from_diags(pair)
        for z,sign in p.prism(vertices):
            def av(f):
                i,j,k,l=f
                return p.A2((z[i],z[j],z[k]))*(z[l][-1]-z[k][-1])
            source=p.make_source(3,len(z)-1,
                lambda f:p.s(tuple(z[i] for i in f)),av,
                lambda f:p.omega(tuple(z[i] for i in f)))
            cm.add(out,rc.basis(rs.to_diags(source)),coefficient*sign)
    return out


def prepare_cycle(name,progress=None):
    def report(event,**fields):
        if progress:progress(event,cycle=name,**fields)
    if name not in CYCLES:raise ValueError('unknown suspension detector')
    path=PREPARED/(name+'.json')
    hashes=preparation_hashes()
    previous=json.loads(path.read_text()) if path.exists() else None
    if previous is not None:
        if previous.get('source_hashes')!=hashes or previous.get('cycle')!=name:
            raise RuntimeError('stale suspension checkpoint: '+str(path))
        if previous.get('status')=='prepared':
            report('cycle_reused')
            return previous
    lower=normalized_lower(CYCLES[name])
    report('normalized_lower',terms=len(lower))
    ef2=sp.ef_values(low.SOURCE_VALUES)
    alpha,beta=lc.constants()['alpha'],lc.constants()['beta']
    A,s,omega=p.A2,p.s,p.omega
    L=p.scale(low.transported(low.pontryagin(omega),A,s),Fraction(1,4))
    M=p.half(p.cup(p.cup(p.cup(s,s),omega),p.binary(A)))
    cube=low.transported(low.transported(A,A,s),A,s)
    old=Fraction(0); orientation=Fraction(0)
    for index,(pair,coefficient) in enumerate(lower.items()):
        z=p.from_diags(pair)
        old+=coefficient*(sp.V2(z,ef2)-alpha*L(z)-beta*M(z))
        if name=='F':orientation+=coefficient*cube(z)
        if name=='Y':orientation+=coefficient*L(z)
        if (index+1)%20==0:report('lower_progress',done=index+1,total=len(lower))
    suspended=normalized_suspension(lower)
    report('normalized_suspension',terms=len(suspended),V2fin=str(old))
    H={};F={}
    for index,(pair,coefficient) in enumerate(suspended.items()):
        cm.add(H,rc.Htot(pair),coefficient)
        cm.add(F,rc.Ftot(pair),coefficient)
        if (index+1)%10==0:report('contractor_progress',done=index+1,
                                  total=len(suspended),H_terms=len(H),F_terms=len(F))
    report('contractor_aggregated',H_terms=len(H),F_terms=len(F))
    # Ordering is fixed by the printed chain algorithms. Save its digest
    # before accepting a resumed partial phase sum.
    items=tuple(H.items())
    digest=hashlib.sha256(repr(items).encode()).hexdigest()
    data=dict(schema=1,status='partial',cycle=name,source_hashes=hashes,
        lower_terms=len(lower),suspended_terms=len(suspended),
        H_terms=len(items),H_digest=digest,H_processed=0,H_value='0',
        F_values=[[rs._key(b),c] for b,c in F.items()],
        V2fin=str(old),orientation=str(orientation))
    if previous is not None:
        if any(previous.get(k)!=data[k] for k in
               ('H_terms','H_digest','F_values','V2fin','orientation')):
            raise RuntimeError('resumed suspension chain differs from checkpoint')
        data['H_processed']=previous['H_processed']
        data['H_value']=previous['H_value']
    _save_checkpoint(path,data)
    _clear_period_sources()
    value=Fraction(data['H_value']);start=data['H_processed']
    source=rs.phi()
    for index in range(start,len(items)):
        pair,coefficient=items[index]
        value+=coefficient*source(rs.from_diags(pair))
        if (index+1)%500==0 or index+1==len(items):
            data.update(H_processed=index+1,H_value=str(value))
            _save_checkpoint(path,data)
            report('source_H_progress',done=index+1,total=len(items),value=str(value))
            _clear_period_sources()
    data.update(status='prepared',H_processed=len(items),H_value=str(value))
    _save_checkpoint(path,data)
    report('cycle_prepared',H_value=str(value),V2fin=str(old),orientation=str(orientation))
    _clear_period_sources()
    return data


def evaluate(progress=None):
    prepared={name:prepare_cycle(name,progress) for name in CYCLES}
    ef=source_coefficients()
    periods={};rows=[]
    for name,data in prepared.items():
        small=sum(c*ef.get(rs._basis_key(b),0) for b,c in data['F_values'])
        upper=Fraction(data['H_value'])+small
        value=upper-Fraction(data['V2fin'])
        periods[name]=value%1
        rows.append(dict(cycle=name,suspended_V3_H=data['H_value'],
            suspended_V3_F=str(small),suspended_V3=str(upper),
            V2fin=data['V2fin'],difference=str(value)))
        if progress:progress('cycle_period',**rows[-1],phase=str(value%1))
    aF=integer(Fraction(prepared['F']['orientation']),'Euler cube orientation')
    jL=integer(4*(Fraction(prepared['Y']['orientation'])%1),'Pontryagin orientation')%4
    if aF not in (-1,1) or jL not in (1,3):
        raise ArithmeticError(f'comparison orientations are not units: {aF},{jL}')
    xi=(aF*periods['F'])%1
    if xi not in (Fraction(1,4),Fraction(3,4)):
        raise ArithmeticError(f'quarter-cubic suspension audit failed: xi={xi}')
    return dict(status='evaluated',schema=1,reference='chi7_tail/epsilon100/eta101/R2sharp',
        source_hashes=provenance(),aF=aF,jL=jL,xi=str(xi),
        c4=jL*integer(4*periods['Y'],'c4')%4,
        cN=integer(2*periods['W'],'cN')%2,
        cO=integer(2*periods['V'],'cO')%2,
        cM=integer(2*periods['ZT'],'cM')%2,
        epsilon_c=integer(2*((xi-Fraction(1,4))%1),'epsilon_c'),
        periods={name:str(value) for name,value in periods.items()},rows=rows)


@lru_cache(None)
def constants():
    if not DATA.exists():raise RuntimeError('fixed R3 suspension calibration is missing')
    data=json.loads(DATA.read_text())
    if (data.get('status')!='evaluated' or data.get('schema')!=1 or
        data.get('reference')!='chi7_tail/epsilon100/eta101/R2sharp' or
        data.get('source_hashes')!=provenance()):
        raise RuntimeError('R3 suspension calibration does not match current sources')
    if data.get('c4') not in range(4) or any(data.get(k) not in (0,1)
        for k in ('cN','cO','cM','epsilon_c')):
        raise RuntimeError('invalid R3 suspension coefficients')
    return data


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=DATA)
    parser.add_argument('--prepare-only',action='store_true')
    parser.add_argument('--cycle',choices=tuple(CYCLES))
    args=parser.parse_args()
    report=lambda event,**kw:print(json.dumps(dict(event=event,**kw)),flush=True)
    if args.prepare_only:
        for name in ([args.cycle] if args.cycle else CYCLES):prepare_cycle(name,report)
        raise SystemExit(0)
    if args.cycle:parser.error('--cycle requires --prepare-only')
    result=evaluate(report)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True)
