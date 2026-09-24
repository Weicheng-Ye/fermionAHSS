"""Current chi7_tail phases on finite universal input simplices.

Exact rational values are canonical real lifts, not division in R/Z.
This module is verification scaffolding, not a production operation.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
import ast
import json
from pathlib import Path
import re

PACKAGE = Path(__file__).resolve().parents[1]
ROOT = PACKAGE.parent
try:
    from .cochain_tools import (Cochain, cup, differential, q2, square, zero,
                               zeta1, zeta2, word_op, cut_terms)
except ImportError:
    from cochain_tools import (Cochain, cup, differential, q2, square, zero,
                               zeta1, zeta2, word_op, cut_terms)


def scale(c, q):
    return Cochain(c.degree, lambda z: q * c(z))


def half(c):
    return scale(c, Fraction(1, 2))


def binary(c):
    return c.mod2()


def ds(c, s):
    return differential(c) - scale(cup(s, c, integral=True), 2)


def divide(c, denominator, name='even quotient'):
    def value(z):
        numerator = c(z)
        if numerator % denominator:
            raise ArithmeticError(f'{name}: {numerator} is not divisible by {denominator}')
        return numerator // denominator
    return Cochain(c.degree, value)


def Q(c, degree):
    return binary(cup(c, c, c.degree-degree)
                  + cup(c, binary(differential(c)), c.degree-degree+1))


def E(c, omega):
    return binary(Q(c, 2) + cup(omega, c))


def QD(c, s, omega):
    return binary(E(c, omega) + cup(s, Q(c, 1)))


@lru_cache(None)
def chi_anf(n):
    if n not in range(8):
        raise ValueError('the calibrated chi family ends at degree seven')
    # Consume the very same certified ANF that GAP evaluates. No note files,
    # external checkout, compiler or expanded word lists are runtime inputs.
    text = (PACKAGE/'data/chi-calibrated-degree7-anf.g').read_text()
    pattern = (r'rec\(degree := '+str(n)+r', bitsPerFace := (\d+),\s*'
               r'wordCount := \d+, cutCount := \d+,\s*'
               r'faces := (\[.*?\]),\s*monomials := (\[.*?\])\s*\)')
    match = re.search(pattern,text,re.S)
    if match is None:
        raise ValueError('missing or malformed packaged chi ANF degree '+str(n))
    width = int(match[1])
    faces = tuple(tuple(i-1 for i in face) for face in ast.literal_eval(match[2]))
    if faces != tuple(combinations(range(n+4),n+1)):
        raise ValueError('packaged chi face ordering changed')
    masks = []
    for code in ast.literal_eval(match[3]):
        mask = 0
        while code:
            face = code & ((1 << width)-1)
            if not 1 <= face <= len(faces):
                raise ValueError('invalid packaged chi face code')
            mask |= 1 << (face-1)
            code >>= width
        masks.append(mask)
    return faces,tuple(masks)


def chi(c):
    faces, masks = chi_anf(c.degree)
    @lru_cache(None)
    def evaluate_active(active):
        # Different universal simplices often restrict the input to the
        # same binary face pattern. The ANF depends only on that pattern.
        return sum((active & mask) == mask for mask in masks) % 2
    def value(z):
        active = sum(c(tuple(z[j] for j in face)) << i for i, face in enumerate(faces))
        return evaluate_active(active)
    return Cochain(c.degree+3, value)


def source(A, s, omega):
    """The A-only source splitting, n=1,2,3, with epsilon100."""
    n = A.degree
    if n not in (1, 2, 3):
        raise ValueError('source splitting is only implemented in degrees1,2,3')
    a = binary(A)
    t = binary(divide(A-a, 2, 'A carry'))
    B = divide(differential(a), 2, 'binary Bockstein')
    e = binary(B)
    CB = binary(divide(B+e, 2, 'plus carry'))
    u, v, w = square(a, 2), cup(omega, a), cup(s, e)
    p = binary(u+v+w)
    FA = binary(zeta2(omega, a)+chi(a)
                +cup(u, v, n+1)+cup(u, w, n+1)+cup(v, w, n+1)
                +zeta1(s, e)+cup(cup(omega, s, 1), e)
                +cup(s, u)+cup(cup(s, s), CB))
    g = binary(Q(t, 2)+cup(omega, t)
               +cup(binary(differential(t)), cup(s, a), n)
               +zeta1(s, a)+cup(cup(omega, s, 1), a))
    qint = cup(omega, B, integral=True)+cup(B, B, n-1, integral=True)
    LA = divide(qint-ds(g, s)+cup(s, p, integral=True), 2, 'source LA')
    k0 = binary(FA+binary(LA)+cup(cup(cup(s, s), s), a))
    result = dict(a=a, t=t, B=B, e=e, CB=CB, p=p, q0=p,
                  g=g, k0=k0, k=k0, q=p, FA=FA, LA=LA)
    if n == 1:
        e0 = cup(s, t)
        de0 = binary(differential(e0))
        result['q'] = cup(binary(omega+cup(s, s)), a)
        result['k'] = binary(k0+QD(e0, s, omega)+cup(de0, p, 2)
                             +cup(s, cup(de0, p, 3)))
    return result


source1 = source


def theta(q, v, s, omega):
    """Exact prescribed boundary phase with current eta101, real lift."""
    m = q.degree
    B = divide(differential(q), 2, 'Theta binary Bockstein')
    e = binary(B)
    CB = binary(divide(B+e, 2, 'Theta plus carry'))
    u, wp, se = square(q, 2), cup(omega, q), cup(s, e)
    H = binary(zeta2(omega, q)+chi(q)
               +cup(u, wp, m+1)+cup(u, se, m+1)+cup(wp, se, m+1)
               +zeta1(s, e)+cup(cup(omega, s, 1), e)
               +cup(s, u)+cup(cup(s, s), CB))
    qint = cup(omega, B, integral=True)+cup(B, B, m-1, integral=True)
    Z = binary(cup(s, square(q, 2))+cup(omega, square(q, 1)))
    return half(binary(E(v, omega)+H))+scale(qint, Fraction(1, 4))+half(Z)


# Vertices (source_id, original_index, *prism_levels). All input tables use
# the first two coordinates; prism cochains may inspect the last coordinate.
SOURCES = []
SOURCE_KEYS = {}


def make_source(n, dim, s_value, A_value, omega_value, b_value=None, c_value=None):
    fields = {}
    for name, degree, fun in [('s', 1, s_value), ('A', n, A_value),
                               ('omega', 2, omega_value),
                               ('b', n+1, b_value), ('c', n+2, c_value)]:
        fields[name] = {indices: fun(indices) if fun else 0
                        for indices in combinations(range(dim+1), degree+1)}
    key = (n, dim, tuple(tuple(table.items()) for table in fields.values()))
    if key in SOURCE_KEYS:
        return tuple((SOURCE_KEYS[key], j) for j in range(dim+1))
    source_id = len(SOURCES)
    SOURCES.append((n, dim, fields))
    SOURCE_KEYS[key] = source_id
    return tuple((source_id, j) for j in range(dim+1))


def input_cochain(name, degree):
    def value(vertices):
        source_id = vertices[0][0]
        if any(v[0] != source_id for v in vertices):
            raise ValueError('mixed source simplices')
        indices = tuple(v[1] for v in vertices)
        if len(set(indices)) != len(indices):
            return 0
        return SOURCES[source_id][2][name][indices]
    return Cochain(degree, value)


s = input_cochain('s', 1)
omega = input_cochain('omega', 2)
A1 = input_cochain('A', 1)
A2 = input_cochain('A', 2)
A3 = input_cochain('A', 3)
b2 = input_cochain('b', 3)
c2 = input_cochain('c', 4)


def prism(vertices):
    """Raw right prism; its boundary is (top-bottom) - prism boundary."""
    n = len(vertices)-1
    for j in range(n+1):
        yield ((tuple(v+(0,) for v in vertices[:j+1])
                +tuple(v+(1,) for v in vertices[j:])), (-1)**j)


def prism_cochain(c):
    return Cochain(c.degree-1,
                   lambda z: sum(sign*c(v) for v, sign in prism(z)))


@lru_cache(None)
def phase2():
    data = source(A2, s, omega)
    p, k, g = data['p'], data['k'], data['g']
    y = QD(b2, s, omega)
    sb = cup(s, b2)
    lam = binary(divide(g+sb-binary(g+sb), 2, 'splitting carry'))
    dlam = binary(differential(lam))
    psi = binary(y+k+dlam)
    pol = binary(cup(y, k, 4)+cup(k, binary(differential(y)), 5)+Q(k, 1))
    correction = binary(E(lam, omega)+cup(binary(y+k), dlam, 4))
    UA = binary(zeta1(s, p)+cup(s, k)+cup(cup(omega, s, 1), p))
    bI = Cochain(3, lambda z: b2(z)*z[-1][-1])
    free = prism_cochain(theta(binary(differential(bI)), QD(bI, s, omega), s, omega))
    Xi = (half(E(c2, omega))+half(pol)+free+half(correction)
          -half(E(sb, omega))+half(UA))
    return dict(**data, psi=psi, Xi=Xi, phi=theta(p, k, s, omega),
                lambda_=lam, free=free, pol=pol, correction=correction, UA=UA)


@lru_cache(None)
def from_diags(pair):
    """Rectangular-sum r on signed and mod-two diagonal bar simplices."""
    signed, binary_diag = pair
    n = len(signed.rows)
    if len(binary_diag.rows) != n:
        raise ValueError('diagonal product dimensions differ')

    def sv(indices):
        i, j = indices
        return sum(row[0] for row in signed.rows[i:j]) % 2

    def av(indices):
        r, l, t = indices
        return sum((-1)**sv((r, i))*signed.rows[i][1][j]
                   for i in range(r, l) for j in range(l, t))

    def ov(indices):
        r, l, t = indices
        return sum(binary_diag.rows[i][1][j]
                   for i in range(r, l) for j in range(l, t)) % 2
    return make_source(2, n, sv, av, ov)


def to_diags(vertices):
    """Signed skew section j (not a replacement for its nontrivial r)."""
    try:
        from .chain_models import Diag
    except ImportError:
        from chain_models import Diag
    n = len(vertices)-1
    matrix, omega_matrix = [[0]*n for _ in range(n)], [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            matrix[i][j] = (A2((vertices[i], vertices[i+1], vertices[j+1]))
                             - A2((vertices[i], vertices[i+1], vertices[j])))
            matrix[j][i] = -(-1)**s((vertices[i], vertices[j]))*matrix[i][j]
            omega_matrix[i][j] = (omega((vertices[i], vertices[i+1], vertices[j+1]))
                                   + omega((vertices[i], vertices[i+1], vertices[j]))) % 2
            omega_matrix[j][i] = omega_matrix[i][j]
    return (Diag('zsign', tuple((s((vertices[i], vertices[i+1])), tuple(matrix[i]))
                               for i in range(n))),
            Diag('c2', tuple((0, tuple(row)) for row in omega_matrix)))


def evaluate_r(c, chain):
    return sum(coefficient*c(from_diags(pair)) for pair, coefficient in chain.items())
