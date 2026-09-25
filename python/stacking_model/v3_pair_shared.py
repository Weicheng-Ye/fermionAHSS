"""Explicit shared-sign n=3 source-pair primitive with arbitrary omega.

The only solve is a fixed universal 12-by-2 binary matrix, independent of
X and its input cochains. All universal chains use r3_pair_chain's shared
sign component, so the two local systems are never separated by EZ.
"""
from functools import lru_cache
from pathlib import Path
import sys

sys.path.insert(0, str(Path(str(Path(__file__).resolve().parents[1]))))
import phase_eval as p
import r3_pair_chain as rc
from g6_repair import f2_solve
from v3_pair_untwisted import from_omega, to_omega


def from_diag(diag):
    def sign_value(indices):
        i, j = indices
        return sum(a.sigma for a in diag.rows[i:j]) % 2
    def component(second):
        def value(indices):
            i, j, k, l = indices
            return sum((-1) ** sign_value((i, r)) * matrix[t][v]
                       for r in range(i, j)
                       for matrix in ((diag.rows[r].matrix_prime if second
                                       else diag.rows[r].matrix),)
                       for t in range(j, k) for v in range(k, l))
        return p.Cochain(3, value)
    return component(False), component(True), p.Cochain(1, sign_value)


def to_diag(A, Aprime, s, vertices):
    n = len(vertices) - 1
    def sv(i, j):
        return 0 if i == j else s((vertices[i], vertices[j]))
    def matrix_for(source, r):
        def abar(indices):
            if len(set(indices)) != 4:
                return 0
            ordered = tuple(sorted(indices))
            sign = (-1) ** sum(indices[i] > indices[j]
                              for i in range(4) for j in range(i + 1, 4))
            return sign * (-1) ** sv(0, ordered[0]) * source(
                tuple(vertices[i] for i in ordered))
        def z(i, j, k):
            return (-1) ** sv(0, r) * (abar((r, i, j, k)) - abar((r + 1, i, j, k)))
        matrix = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                matrix[i][j] = z(i, i + 1, j + 1) - z(i, i + 1, j)
                matrix[j][i] = -matrix[i][j]
        return tuple(map(tuple, matrix))
    return rc.Diag3(tuple(rc.U2(sv(r, r + 1), matrix_for(A, r), matrix_for(Aprime, r))
                          for r in range(n)))


def hD(x, y, s):
    dx = p.binary(p.differential(x))
    return p.binary(p.cup(x, y, x.degree - 1) + p.cup(dx, y, x.degree)
                    + p.cup(s, p.binary(p.cup(x, y, x.degree)
                                         + p.cup(dx, y, x.degree + 1))))


def source_pair(A, Aprime, s, omega):
    a, ap = p.binary(A), p.binary(Aprime)
    alpha = hD(a, ap, s)
    primary, pp = p.QD(a, s, omega), p.QD(ap, s, omega)
    return p.binary(p.source(A + Aprime, s, omega)['k0']
                    + p.source(A, s, omega)['k0'] + p.source(Aprime, s, omega)['k0']
                    + p.QD(alpha, s, omega) + hD(primary, pp, s)
                    + hD(p.binary(primary + pp), p.binary(p.differential(alpha)), s))


@lru_cache(None)
def source_value(pair):
    if rc.degree(pair) != 6:
        raise ValueError('source pair has degree six')
    if (not any(any(any(row) for row in a.matrix) for a in pair[0].rows)
            or not any(any(any(row) for row in a.matrix_prime) for a in pair[0].rows)):
        return 0
    A, Aprime, s = from_diag(pair[0])
    return source_pair(A, Aprime, s, from_omega(pair[1]))(tuple(range(7)))


def evaluate_source(chain):
    return sum(coefficient * source_value(pair) for pair, coefficient in chain.items()) % 2


@lru_cache(None)
def fixed_certificate():
    basis5, basis6 = rc.total_basis(5), rc.total_basis(6)
    indices = {basis: i for i, basis in enumerate(basis5)}
    rows, values = [], []
    for basis in basis6:
        boundary = rc.total_boundary(basis)
        rows.append(sum((coefficient % 2) << indices[face]
                        for face, coefficient in boundary.items() if face in indices))
        values.append(evaluate_source(rc.Gtot(basis)))
    rhs = sum(value << i for i, value in enumerate(values))
    solution = f2_solve(tuple(rows), rhs, len(basis5))
    if solution is None:
        raise ArithmeticError('fixed shared-source pair primitive does not exist')
    return dict(basis5=basis5, basis6=basis6, rows=tuple(rows), values=tuple(values),
                coefficients=tuple((solution >> i) & 1 for i in range(len(basis5))))


@lru_cache(None)
def primitive_value(pair):
    """Reuse the fixed universal value across independently built cochains."""
    certificate = fixed_certificate()
    coefficients = dict(zip(certificate['basis5'], certificate['coefficients']))
    if (not any(any(any(row) for row in a.matrix) for a in pair[0].rows)
            or not any(any(any(row) for row in a.matrix_prime) for a in pair[0].rows)):
        return 0
    return (evaluate_source(rc.Htot(pair))
            +sum(coefficient*coefficients.get(basis, 0)
                 for basis, coefficient in rc.Ftot(pair).items())) % 2


def primitive(A, Aprime, s, omega):
    if (A.degree, Aprime.degree, s.degree, omega.degree) != (3, 3, 1, 2):
        raise ValueError('V3 expects degrees (3,3,1,2)')
    def value(vertices):
        pair = (to_diag(A, Aprime, s, vertices), to_omega(omega, vertices))
        return primitive_value(pair)
    return p.Cochain(5, value)



def legal_beta(A, B, Aprime, Bprime, s, omega):
    """Degree-six beta on legal first-two-layer data; no C arguments."""
    alpha = hD(p.binary(A), p.binary(Aprime), s)
    Asum = A + Aprime
    Bsum = p.binary(B + Bprime + alpha)
    def splitting_carry(a, b):
        g = p.source(a, s, omega)['g']
        summand = g + p.cup(s, b)
        return p.binary(p.divide(summand - p.binary(summand), 2,
                                 'stacking secondary splitting carry'))
    base = p.binary(hD(B, Bprime, s) + hD(p.binary(B + Bprime), alpha, s)
                    + splitting_carry(A, B) + splitting_carry(Aprime, Bprime)
                    + splitting_carry(Asum, Bsum))
    return p.binary(base + primitive(A, Aprime, s, omega))


if __name__ == '__main__':
    import json
    print(json.dumps(fixed_certificate(), indent=2), flush=True)
