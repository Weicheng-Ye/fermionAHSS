"""Fixed finite epsilon and alpha selectors for the current Danus reference."""
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
import json
try:
    from . import low_phases as l, phase_eval as p, chain_models as cm
    from . import source_primitive as sp
except ImportError:
    import low_phases as l, phase_eval as p, chain_models as cm, source_primitive as sp

DATA=Path(__file__).with_name('low_calibration.json')


def _integer(value,label):
    value=Fraction(value)
    if value.denominator!=1: raise ArithmeticError(f'{label} is not integral: {value}')
    return int(value)


def _background_cycle(): return cm.G('c2',(1,1))


def cycle(group_edge):
    out={}
    for w,c in _background_cycle().items():
        cm.add(out,l.group_product_shuffle((((0,0),group_edge),w)),c)
    return out


def _suspended_v2(vertices):
    total=Fraction(0); ef=sp.ef_values(l.SOURCE_VALUES)
    for z,sign in p.prism(vertices):
        def av(f):
            i,j,k=f
            return p.A1((z[i],z[j]))*(z[k][-1]-z[j][-1])
        source=p.make_source(2,len(z)-1,
            lambda f:p.s(tuple(z[i] for i in f)),av,
            lambda f:p.omega(tuple(z[i] for i in f)))
        total-=sign*sp.V2(source,ef)
    return total


def evaluate(progress=None):
    def report(event,**data):
        if progress:progress(event,**data)
    P=l.pontryagin(p.omega)
    background=sum(c*l._background_value(P,w) for w,c in _background_cycle().items())
    j=_integer(background,'Pontryagin cycle')%4
    if j not in (1,3):raise ArithmeticError('Pontryagin cycle is not a generator')
    rank=Fraction(0); rank_rows=[]
    for i,(pair,c) in enumerate(cycle((2,0)).items()):
        z=l.from_pair1(pair)
        # On the binary oriented K circle, all half-valued Reven terms
        # vanish and Rraw differs from -V by an exact phase Sq1(k)/2.
        v=l.V1(z)
        K=p.divide(p.A1,2,'rank even input')
        even=Fraction(-l.transported(P,K,p.s)(z),8)
        value=-v-even
        rank+=c*value
        rank_rows.append(dict(index=i,coefficient=c,V=str(v),Reven=str(even),difference=str(value)))
        report('rank_term',**rank_rows[-1])
    rank*=j
    epsilon=_integer(2*(rank%1),'epsilon')
    if epsilon not in (0,1):raise ArithmeticError('rank discrepancy is not order two')
    report('epsilon',value=epsilon,period=str(rank),j=j)
    yperiod=Fraction(0); yrows=[]
    for i,(pair,c) in enumerate(cycle((1,1)).items()):
        z=l.from_pair1(pair)
        vs=_suspended_v2(z)
        v1=l.V1(z)
        l1=Fraction(l.transported(P,p.A1,p.s)(z),4)
        # e0=s*t is zero on the literal odd branch, hence J_A=0.
        val=vs-v1-epsilon*l1
        yperiod+=c*val
        yrows.append(dict(index=i,coefficient=c,suspended_V2=str(vs),V1=str(v1),L1=str(l1),difference=str(val)))
        report('alpha_term',**yrows[-1])
    alpha=j*_integer(4*(yperiod%1),'alpha quarter')%4
    report('alpha',value=alpha,period=str(yperiod))
    return dict(status='evaluated',reference='chi7_tail/epsilon100/eta101',
        epsilon=epsilon,alpha=alpha,beta=0,j=j,
        rank_period=str(rank),alpha_period=str(yperiod),
        odd_base_normalization=l.odd_base_normalization()[1],
        rank_rows=rank_rows,alpha_rows=yrows)


@lru_cache(None)
def constants():
    if not DATA.exists():
        raise RuntimeError('Danus normalization data missing; run low_calibration.py first')
    data=json.loads(DATA.read_text())
    if data.get('status')!='evaluated' or data.get('reference')!='chi7_tail/epsilon100/eta101':
        raise RuntimeError('Danus normalization data have the wrong convention')
    if data.get('epsilon') not in (0,1) or data.get('alpha') not in (0,1,2,3) or data.get('beta')!=0:
        raise RuntimeError('invalid Danus normalization coefficients')
    return data


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=DATA)
    args=parser.parse_args()
    result=evaluate(lambda event,**kw:print(json.dumps(dict(event=event,**kw)),flush=True))
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if not k.endswith('_rows')}),flush=True)
