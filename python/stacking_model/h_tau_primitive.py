"""Fixed-source comparison of raw H and calibrated Tau.

Research implementation: calibration failures are reported, never replaced by
a local solve. Source formulas are the production uniform splitting in all
needed degrees. Existing universal contractors are used without modification.
"""
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import os
import sys

ROOT = Path(os.environ.get('FERMIONAHSS_ROOT', str(Path(__file__).resolve().parents[2])))
sys.path.insert(0, str(ROOT / 'python'))
import phase_eval as p
import low_phases as low
import chain_models as cm
import source_primitive as sp
import r3_chain as rc
import r3_source as r3s
import exterior_bar as eb
from g6_repair import f2_solve


def hD(x, y, s):
    dx = p.binary(p.differential(x))
    return p.binary(p.cup(x,y,x.degree-1)+p.cup(dx,y,x.degree)
                    +p.cup(s,p.binary(p.cup(x,y,x.degree)+p.cup(dx,y,x.degree+1))))


def uniform_source(A, s, omega):
    """The n-independent k0,g specialization of phase_eval.source."""
    n = A.degree
    a = p.binary(A)
    t = p.binary(p.divide(A-a,2,'A carry'))
    B = p.divide(p.differential(a),2,'binary Bockstein')
    e = p.binary(B)
    CB = p.binary(p.divide(B+e,2,'plus carry'))
    u,v,w = p.square(a,2),p.cup(omega,a),p.cup(s,e)
    primary = p.binary(u+v+w)
    FA = p.binary(p.zeta2(omega,a)+p.chi(a)
                  +p.cup(u,v,n+1)+p.cup(u,w,n+1)+p.cup(v,w,n+1)
                  +p.zeta1(s,e)+p.cup(p.cup(omega,s,1),e)
                  +p.cup(s,u)+p.cup(p.cup(s,s),CB))
    g = p.binary(p.Q(t,2)+p.cup(omega,t)
                 +p.cup(p.binary(p.differential(t)),p.cup(s,a),n)
                 +p.zeta1(s,a)+p.cup(p.cup(omega,s,1),a))
    q = p.cup(omega,B,integral=True)+p.cup(B,B,n-1,integral=True)
    LA = p.divide(q-p.ds(g,s)+p.cup(s,primary,integral=True),2,'source LA')
    k0 = p.binary(FA+p.binary(LA)+p.cup(p.cup(p.cup(s,s),s),a))
    return dict(k0=k0,g=g,p=primary)


def interval(c, scaled=False):
    def value(vertices):
        base = tuple(v[0] for v in vertices)
        if any(a == b for a,b in zip(base,base[1:])):
            return 0
        result = c(base)
        return result*vertices[-1][1] if scaled else result
    return p.Cochain(c.degree,value)


def prism(c):
    return p.Cochain(c.degree-1,lambda vertices: sum(
        (-1)**j*c(tuple((v,0) for v in vertices[:j+1])+tuple((v,1) for v in vertices[j:]))
        for j in range(len(vertices))))


def closed_source(A,s,omega):
    """Z_n with Tau+H=Z_n+delta L_n on legal lower pairs."""
    source = uniform_source(A,s,omega)
    ai,si,wi = interval(A,True),interval(s),interval(omega)
    AI = p.ds(ai,si)
    W = p.QD(p.binary(ai),si,wi)
    suspended = uniform_source(AI,si,wi)['k0']
    P = source['p']
    return p.binary(source['k0']+prism(p.binary(suspended+p.QD(W,si,wi)))+hD(P,P,s))


def splitting_carry(A,B,s,omega):
    total = uniform_source(A,s,omega)['g']+p.cup(s,B)
    return p.binary(p.divide(total-p.binary(total),2,'splitting carry'))


def known_gauge(A,B,s,omega):
    ai,bi,si,wi = interval(A,True),interval(B,True),interval(s),interval(omega)
    AI,V = p.ds(ai,si),p.binary(p.differential(bi))
    W = p.QD(p.binary(ai),si,wi)
    BI = p.binary(V+W)
    return p.binary(splitting_carry(A,B,s,omega)+prism(p.QD(bi,si,wi))
                    +prism(hD(V,W,si))+prism(splitting_carry(AI,BI,si,wi)))


def from_pair(pair,n):
    signed,omega = pair
    if n==2:
        def sv(f):
            return sum(row[0] for row in signed.rows[f[0]:f[1]])%2
        def av(f):
            i,j,k=f
            return sum((-1)**sv((i,r))*signed.rows[r][1][t]
                       for r in range(i,j) for t in range(j,k))
    elif n==3:
        def sv(f):
            return sum(row.sigma for row in signed.rows[f[0]:f[1]])%2
        def av(f):
            i,j,k,l=f
            return sum((-1)**sv((i,r))*signed.rows[r].matrix[t][v]
                       for r in range(i,j) for t in range(j,k) for v in range(k,l))
    else:
        g=signed
        def sv(f):return g[f[0]][1]^g[f[1]][1]
        def av(f):return (-1)**g[f[0]][1]*(g[f[1]][0]-g[f[0]][0])
    def ov(f):
        i,j,k=f
        return sum(omega.rows[r][1][t] for r in range(i,j) for t in range(j,k))%2
    return p.Cochain(n,av),p.Cochain(1,sv),p.Cochain(2,ov)


@lru_cache(None)
def source_value(pair,n):
    return closed_source(*from_pair(pair,n))(tuple(range(n+4)))


def evaluate(chain,n):
    return sum(c*source_value(pair,n) for pair,c in chain.items())%2


def to_pair(A,s,omega,vertices):
    n=A.degree
    registered=p.make_source(n,len(vertices)-1,
        lambda f:s(tuple(vertices[i] for i in f)),
        lambda f:A(tuple(vertices[i] for i in f)),
        lambda f:omega(tuple(vertices[i] for i in f)))
    if n==1:
        return low.to_pair1(registered,p.A1,p.s,p.omega)
    if n==2:
        return p.to_diags(registered)
    if n==3:
        return r3s.to_diags(registered)
    raise ValueError('source contraction exists only in degrees one through three')


def group_homotopy(pair):
    """Dihedral wedge retraction followed by the omega contraction."""
    result=dict(low.group_product_H(pair))
    for (g,w),coefficient in low.group_product_AW(pair).items():
        for h,c in low.homotopy_group(g).items():
            cm.add(result,low.group_product_shuffle((h,w)),coefficient*c)
        for r,c in low.retract_group(g).items():
            for h,d in cm.H(w).items():
                cm.add(result,low.group_product_shuffle((r,h)),coefficient*(-1)**(len(g)-1)*c*d)
    return result


@lru_cache(None)
def wedge_certificate():
    values=[]
    bases=[]
    for degree in range(1,5):
        for word in rc.omega_words(4-degree):
            g=tuple((0,0) if i%2==0 else (1,1) for i in range(degree+1))
            chain={}
            for w,c in cm.G('c2',word).items():
                cm.add(chain,low.group_product_shuffle((g,w)),c)
            values.append(evaluate(chain,1))
            bases.append((degree,word))
    return dict(n=1,basis_upper=bases,values=values,coefficients=[] if not any(values) else None)


def basis2(degree):
    return tuple(((p,k),w) for k in range(1,degree//2+1)
                 for p in range(degree-2*k+1) for w in rc.omega_words(degree-p-2*k))


def boundary2(basis):
    (p,k),w=basis
    result={(a,w):c for a,c in eb.small_boundary(p,k).items()}
    for v,c in rc.omega_boundary(w).items():
        cm.add(result,{((p,k),v):(-1)**p*c})
    return result


@lru_cache(None)
def certificate(n):
    if n==1:
        return wedge_certificate()
    if n not in (2,3):
        raise ValueError('small source matrix covers n=1,2,3')
    basisfun,boundary,G=(basis2,boundary2,sp.Gtot) if n==2 else (rc.total_basis,rc.total_boundary,rc.Gtot)
    lower,upper=basisfun(n+2),basisfun(n+3)
    indices={b:i for i,b in enumerate(lower)}
    rows=[]
    values=[]
    for index,b in enumerate(upper):
        rows.append(sum((c%2)<<indices[a] for a,c in boundary(b).items() if a in indices))
        chain=G(b)
        value=evaluate(chain,n)
        values.append(value)
    rhs=sum(value<<i for i,value in enumerate(values))
    solution=f2_solve(rows,rhs,len(lower))
    return dict(n=n,basis_lower=lower,basis_upper=upper,rows=rows,values=values,
                coefficients=None if solution is None else [(solution>>i)&1 for i in range(len(lower))])


def source_primitive(A,s,omega):
    """v_n=Z_n H; its entire fixed small restriction vanishes."""
    n=A.degree
    if (n,s.degree,omega.degree) not in ((1,1,2),(2,1,2),(3,1,2)):
        raise ValueError('source primitive expects degrees (n,1,2), n=1,2,3')
    cert=certificate(n)
    if any(cert['values']):
        raise ArithmeticError('the fixed comparison source has a nonzero small value')
    H=group_homotopy if n==1 else sp.Htot if n==2 else rc.Htot
    return p.Cochain(n+2,lambda vertices:evaluate(H(to_pair(A,s,omega,vertices)),n))


def comparison_gauge(A,B,s,omega):
    """h with delta h=Tau_n(A;B)+H_(n+3)(A,B), on legal pairs only."""
    if A.degree==0:
        return p.zero(2)
    if B.degree!=A.degree+1:
        raise ValueError('B must have degree deg(A)+1')
    return p.binary(known_gauge(A,B,s,omega)+source_primitive(A,s,omega))


if __name__=='__main__':
    import argparse,json
    parser=argparse.ArgumentParser()
    parser.add_argument('--degree',type=int,choices=(1,2,3),required=True)
    parser.add_argument('--report',type=Path)
    args=parser.parse_args()
    report=json.dumps(certificate(args.degree),indent=2)+'\n'
    print(report,end='',flush=True)
    if args.report:
        args.report.write_text(report)
