"""Fixed pair contraction for the n=3 binary source, with s=0 and arbitrary omega.

This uses the existing universal r3_chain contractions on both K(Z,3)
factors. It does not choose a primitive on an input space.
"""
from functools import lru_cache
from itertools import combinations
import sys
from pathlib import Path

REPO = Path(str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(REPO / 'python'))
import phase_eval as p
import r3_chain as rc


def from_diag(diag):
    def value(indices):
        i, j, k, l = indices
        return sum(diag.rows[r].matrix[t][v]
                   for r in range(i, j) for t in range(j, k) for v in range(k, l))
    return p.Cochain(3, value)


def to_diag(A, vertices):
    n = len(vertices) - 1
    def abar(indices):
        if len(set(indices)) != 4:
            return 0
        ordered = tuple(sorted(indices))
        sign = (-1) ** sum(indices[i] > indices[j]
                          for i in range(4) for j in range(i + 1, 4))
        return sign * A(tuple(vertices[i] for i in ordered))
    rows = []
    for r in range(n):
        def z(i, j, k):
            return abar((r, i, j, k)) - abar((r + 1, i, j, k))
        matrix = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                matrix[i][j] = z(i, i + 1, j + 1) - z(i, i + 1, j)
                matrix[j][i] = -matrix[i][j]
        rows.append(rc.U2(0, tuple(map(tuple, matrix))))
    return rc.Diag3(tuple(rows))


def h(x, y):
    return p.binary(p.cup(x, y, x.degree - 1)
                    + p.cup(p.binary(p.differential(x)), y, x.degree))


def source_pair(A, Aprime, omega=None):
    s, w = p.zero(1), omega if omega is not None else p.zero(2)
    a, ap = p.binary(A), p.binary(Aprime)
    alpha = h(a, ap)
    primary, pp = p.QD(a, s, w), p.QD(ap, s, w)
    return p.binary(p.source(A + Aprime, s, w)['k0']
                    + p.source(A, s, w)['k0'] + p.source(Aprime, s, w)['k0']
                    + p.QD(alpha, s, w) + h(primary, pp)
                    + h(p.binary(primary + pp), p.binary(p.differential(alpha))))


@lru_cache(None)
def source_value(pair):
    degree = rc.degree(pair)
    if degree != 6:
        raise ValueError('source pair has degree six')
    return source_pair(from_diag(pair[0]), from_diag(pair[1]))(tuple(range(7)))


def pair_projection(pair):
    out = {}
    for (x, y), coefficient in rc.product_AW(pair).items():
        for left, a in rc.F(x).items():
            for right, b in rc.F(y).items():
                rc.add(out, {(left, right): coefficient * a * b})
    return out


@lru_cache(None)
def pair_inclusion(words):
    out = {}
    for x, a in rc.G(words[0]).items():
        for y, b in rc.G(words[1]).items():
            rc.add(out, rc.product_shuffle((x, y)), a * b)
    return out


def pair_homotopy(pair):
    out = dict(rc.product_homotopy(pair))
    for (x, y), coefficient in rc.product_AW(pair).items():
        for a, weight in rc.H(x).items():
            rc.add(out, rc.product_shuffle((a, y)), coefficient * weight)
        for a, weight in rc.linear(rc.F(x), rc.G).items():
            for b, other in rc.H(y).items():
                rc.add(out, rc.product_shuffle((a, b)),
                       coefficient * (-1) ** rc.degree(x) * weight * other)
    return out


def evaluate_source(chain):
    return sum(coefficient * source_value(pair) for pair, coefficient in chain.items()) % 2


def primitive(A, Aprime):
    return p.Cochain(5, lambda vertices: evaluate_source(
        pair_homotopy((to_diag(A, vertices), to_diag(Aprime, vertices)))))



def from_omega(diag):
    def value(indices):
        i, j, k = indices
        return sum(diag.rows[r][1][t]
                   for r in range(i, j) for t in range(j, k)) % 2
    return p.Cochain(2, value)


def to_omega(omega, vertices):
    n = len(vertices) - 1
    matrix = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            second = 0 if j == i + 1 else omega((vertices[i], vertices[i + 1], vertices[j]))
            matrix[i][j] = (omega((vertices[i], vertices[i + 1], vertices[j + 1]))
                            + second) % 2
            matrix[j][i] = matrix[i][j]
    return rc.cm.Diag('c2', tuple((0, tuple(row)) for row in matrix))


@lru_cache(None)
def source_background_value(triple):
    pair, omega = triple
    if rc.degree(triple) != 6:
        raise ValueError('source pair has degree six')
    return source_pair(from_diag(pair[0]), from_diag(pair[1]),
                       from_omega(omega))(tuple(range(7)))


def background_homotopy(triple):
    out = dict(rc.product_homotopy(triple))
    for (pair, omega), coefficient in rc.product_AW(triple).items():
        for a, weight in pair_homotopy(pair).items():
            rc.add(out, rc.product_shuffle((a, omega)), coefficient * weight)
        for a, weight in rc.linear(pair_projection(pair), pair_inclusion).items():
            for b, other in rc.cm.H(omega).items():
                rc.add(out, rc.product_shuffle((a, b)),
                       coefficient * (-1) ** rc.degree(pair) * weight * other)
    return out


def primitive_with_omega(A, Aprime, omega):
    def value(vertices):
        triple = ((to_diag(A, vertices), to_diag(Aprime, vertices)),
                  to_omega(omega, vertices))
        return sum(coefficient * source_background_value(simplex)
                   for simplex, coefficient in background_homotopy(triple).items()) % 2
    return p.Cochain(5, value)


if __name__ == '__main__':
    import json
    word = ((0, 1),)
    chain = pair_inclusion((word, word))
    print(json.dumps({'fundamental_pair_chain_terms': len(chain),
                      'fundamental_pair_period': evaluate_source(chain)}), flush=True)
