"""Explicit shared-sign n=4 source-pair primitive with arbitrary omega.

The relative fiber tensor complex starts in degree eight, so its degree-seven
source projection is zero. Only common-sign Borel and equivariant fiber EZ
operators occur, with no source solve or calibration.
"""
from functools import lru_cache
from pathlib import Path
import sys

sys.path.insert(0, str(Path(str(Path(__file__).resolve().parents[1]))))
import phase_eval as p
import r4_pair_ez as rc
from h_tau_primitive import uniform_source
from v3_pair_untwisted import from_omega, to_omega


def from_diag(diag):
    def sv(f):
        return sum(diag.signs[f[0]:f[1]]) % 2
    def component(second):
        simplex=diag.right if second else diag.left
        values=dict(zip(rc.faces(simplex.q,4),simplex.values))
        return p.Cochain(4,lambda f: (-1)**sum(diag.signs[:f[0]])*values[f])
    return component(False),component(True),p.Cochain(1,sv)


def to_diag(A,Aprime,s,vertices):
    q=len(vertices)-1
    signs=tuple(s((vertices[i],vertices[i+1])) for i in range(q))
    def fiber(source):
        return rc.KZ(4,q,tuple((-1)**sum(signs[:f[0]])*source(tuple(vertices[i] for i in f))
                               for f in rc.faces(q,4)))
    return rc.Borel(signs,fiber(A),fiber(Aprime))


def hD(x, y, s):
    dx = p.binary(p.differential(x))
    return p.binary(p.cup(x, y, x.degree - 1) + p.cup(dx, y, x.degree)
                    + p.cup(s, p.binary(p.cup(x, y, x.degree)
                                         + p.cup(dx, y, x.degree + 1))))


def source_pair(A, Aprime, s, omega):
    a, ap = p.binary(A), p.binary(Aprime)
    alpha = hD(a, ap, s)
    primary, pp = p.QD(a, s, omega), p.QD(ap, s, omega)
    return p.binary(uniform_source(A + Aprime, s, omega)['k0']
                    + uniform_source(A, s, omega)['k0'] + uniform_source(Aprime, s, omega)['k0']
                    + p.QD(alpha, s, omega) + hD(primary, pp, s)
                    + hD(p.binary(primary + pp), p.binary(p.differential(alpha)), s))


@lru_cache(None)
def source_value(pair):
    if rc.degree(pair) != 7:
        raise ValueError('source pair has degree seven')
    if not any(pair[0].left.values) or not any(pair[0].right.values):
        return 0
    A, Aprime, s = from_diag(pair[0])
    return source_pair(A, Aprime, s, from_omega(pair[1]))(tuple(range(8)))


def evaluate_source(chain):
    return sum(source_value(pair) for pair, coefficient in chain.items() if coefficient % 2) % 2


def primitive(A, Aprime, s, omega):
    """V4=S H; the relative tensor target has no cells below degree eight."""
    if (A.degree, Aprime.degree, s.degree, omega.degree) != (4, 4, 1, 2):
        raise ValueError('V4 expects degrees (4,4,1,2)')
    def value(vertices):
        pair=(to_diag(A,Aprime,s,vertices),to_omega(omega,vertices))
        return evaluate_source(rc.Htot(pair))
    return p.Cochain(6,value)


def legal_beta(A, B, Aprime, Bprime, s, omega):
    """Degree-seven beta on legal first-two-layer data; no C arguments."""
    alpha = hD(p.binary(A), p.binary(Aprime), s)
    Asum = A + Aprime
    Bsum = p.binary(B + Bprime + alpha)
    def splitting_carry(a, b):
        g = uniform_source(a, s, omega)['g']
        summand = g + p.cup(s, b)
        return p.binary(p.divide(summand - p.binary(summand), 2,
                                 'stacking secondary splitting carry'))
    base = p.binary(hD(B, Bprime, s) + hD(p.binary(B + Bprime), alpha, s)
                    + splitting_carry(A, B) + splitting_carry(Aprime, Bprime)
                    + splitting_carry(Asum, Bsum))
    return p.binary(base + primitive(A, Aprime, s, omega))


