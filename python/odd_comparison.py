"""The exact Appendix B.6 filling map and current suspension beta selector."""
from fractions import Fraction
from functools import lru_cache

try:
    from . import chain_models as cm, phase_eval as pe, source_primitive as sp, odd_primitive as op
except ImportError:
    import chain_models as cm
    import phase_eval as pe
    import source_primitive as sp
    import odd_primitive as op


def _mod2(chain):
    return {a: 1 for a, c in chain.items() if c % 2}


@lru_cache(None)
def from_pair(pair):
    group, background = pair
    n = len(group.edges)
    def sv(indices):
        r, t = indices
        return sum(group.edges[r:t]) % 2
    def ov(indices):
        r, l, t = indices
        return sum(background.rows[i][1][j]
                   for i in range(r,l) for j in range(l,t)) % 2
    return pe.make_source(1, n, sv, sv, ov)


def to_pair(vertices):
    n = len(vertices)-1
    group = cm.Atom('c2', 0, tuple(pe.s((vertices[i], vertices[i+1]))
                                     for i in range(n)))
    matrix = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            matrix[i][j] = (pe.omega((vertices[i], vertices[i+1], vertices[j+1]))
                             +pe.omega((vertices[i], vertices[i+1], vertices[j]))) % 2
            matrix[j][i] = matrix[i][j]
    return group, cm.Diag('c2', tuple((0,tuple(row)) for row in matrix))


def F(pair):
    out = {}
    for (a,b), c in cm.product_AW(pair).items():
        for w,d in cm.F(b).items():
            cm.add(out, {(len(a.edges),w):c*d})
    return out


@lru_cache(None)
def G(basis):
    p,w = basis
    x = cm.Atom('c2',0,(1,)*p)
    out = {}
    for y,c in cm.G('c2',w).items():
        cm.add(out,cm.product_shuffle((x,y)),c)
    return out


def H(pair):
    out = dict(cm.product_homotopy(pair))
    for (a,b),c in cm.product_AW(pair).items():
        for bh,d in cm.H(b).items():
            cm.add(out,cm.product_shuffle((a,bh)),(-1)**len(a.edges)*c*d)
    return out


HSM = {(3,):(1,2), (5,):(1,4), (1,3):(1,1,2), (3,1):(2,1,1)}


def T(pair):
    out = H(pair)
    for (p,w),c in F(pair).items():
        if w in HSM:
            cm.add(out,G((p,HSM[w])),c)
    return _mod2(out)


def Fill(gamma):
    return pe.Cochain(5, lambda z: sum(gamma(from_pair(pair))
                                      for pair in T(to_pair(z))) % 2)


CYCLES = ({(6,()):1}, {(4,(1,)):1}, {(3,(2,)):1}, {(2,(1,1)):1},
          {(1,(4,)):1}, {(1,(1,2)):1,(1,(2,1)):1},
          {(0,(2,2)):1}, {(0,(1,1,1)):1})


def odd_values(progress=None):
    data = op.build_odd(pe.A1,pe.s,pe.omega)
    values = []
    for i,cycle in enumerate(CYCLES):
        chain = _mod2(cm.linear(cycle,G))
        value = sum(data['kappa'](from_pair(pair)) for pair in chain) % 2
        values.append(value)
        if progress:
            progress('odd_source_value', index=i+1,terms=len(chain),value=value)
    # This checks the source-obstruction vanishing rather than assuming it.
    op.build_gamma(data,values)
    return values


def tor_chain():
    return cm.add(dict(G((2,(2,)))),G((3,(1,))),-1)


def odd_period(values, progress=None):
    data = op.build_odd(pe.A1,pe.s,pe.omega)
    raw = op.build_Uraw(data,values,Fill)
    total = Fraction(0)
    rows = []
    for i,(pair,c) in enumerate(tor_chain().items()):
        value = raw(from_pair(pair))
        total += c*value
        rows.append(dict(index=i,coefficient=c,value=str(value)))
        if progress:
            progress('odd_tor_term',**rows[-1])
    return total, rows


def suspended_source_chain():
    out = {}
    for pair,c in tor_chain().items():
        vertices = from_pair(pair)
        for prism,sign in pe.prism(vertices):
            n = len(prism)-1
            def av(indices):
                i,j,k = indices
                return pe.A1((prism[i],prism[j]))*(prism[k][-1]-prism[j][-1])
            source = pe.make_source(2,n,
                lambda indices:pe.s(tuple(prism[i] for i in indices)),av,
                lambda indices:pe.omega(tuple(prism[i] for i in indices)))
            # kappa_6 is minus the raw right prism.
            cm.add(out,{source:-c*sign})
    return out


def suspended_v2(values, progress=None, part=None):
    ef = sp.ef_values(values)
    chain = suspended_source_chain()
    total, rows = Fraction(0), []
    for i,(vertices,c) in enumerate(chain.items()):
        if part is not None and part != i:
            continue
        value = sp.V2(vertices,ef)
        total += c*value
        rows.append(dict(index=i,coefficient=c,value=str(value)))
        if progress:
            progress('suspended_v2_term',total_terms=len(chain),**rows[-1])
    return total, rows
