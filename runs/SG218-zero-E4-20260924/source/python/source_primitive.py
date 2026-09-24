"""Finite V2 contractor prescribed in current r2.md, with no local solves."""
from fractions import Fraction
from functools import lru_cache
try:
    from . import chain_models as cm, phase_eval as pe
except ImportError:
    import chain_models as cm
    import phase_eval as pe


def _tensor_maps(pair, left, right):
    x, y = pair
    out = {}
    for a, c in left(x).items():
        for b, d in right(y).items():
            cm.add(out, {(a, b): c*d})
    return out


def Ftot(pair):
    return cm.linear(cm.product_AW(pair), lambda xy: _tensor_maps(xy, cm.F, cm.F))


@lru_cache(None)
def Gtot(basis):
    u, w = basis
    return cm.linear(_tensor_maps((u, w), lambda a: cm.G('zsign', a),
                                   lambda b: cm.G('c2', b)), cm.product_shuffle)


def Htot(pair):
    out = dict(cm.product_homotopy(pair))
    for (x, y), coefficient in cm.product_AW(pair).items():
        tens = {(h, y): c for h, c in cm.H(x).items()}
        gf = cm.linear(cm.F(x), lambda a: cm.G('zsign', a))
        for a, c in gf.items():
            for b, d in cm.H(y).items():
                cm.add(tens, {(a, b): (-1)**len(x.rows)*c*d})
        cm.add(out, cm.linear(tens, cm.product_shuffle), coefficient)
    return out


SOURCE_BASES = (((2,1),(2,)), ((0,1),(1,2)), ((0,1),(2,1)),
                ((3,2),()), ((1,2),(1,)))


def source_values(progress=None):
    values = []
    phi = pe.phase2()['phi']
    for basis in SOURCE_BASES:
        chain = Gtot(basis)
        if progress:
            progress('source_chain', basis=basis, terms=len(chain))
        value = Fraction(0)
        for index, (pair, coefficient) in enumerate(chain.items()):
            value += coefficient*phi(pe.from_diags(pair))
            if progress and (index+1) % 100 == 0:
                progress('source_progress', basis=basis, evaluated=index+1,
                         total=len(chain))
        values.append(value % 1)
        if progress:
            progress('source_value', basis=basis, raw=str(value),
                     phase=str(value % 1))
    return values


def ef_values(values):
    C, G, H, A, B = map(Fraction, values)
    return {((2,1),(1,)): -C/2,
            ((0,1),(3,)): (H-G)/2,
            ((0,1),(1,1)): -(G+H)/4,
            ((2,2),()): -A/2,
            ((0,2),(1,)): -B/2}


def V2(vertices, ef):
    pair = pe.to_diags(vertices)
    return (pe.evaluate_r(pe.phase2()['phi'], Htot(pair))
            +sum(coefficient*ef.get(key, 0) for key, coefficient in Ftot(pair).items()))
