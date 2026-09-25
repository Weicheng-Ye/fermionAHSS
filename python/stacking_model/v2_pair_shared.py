"""Shared-sign degree-two source-pair contraction and fixed calibration."""
from functools import lru_cache
from pathlib import Path
import sys

sys.path.insert(0, str(Path(str(Path(__file__).resolve().parents[1]))))
import phase_eval as p
import r2_pair_chain as rc
from g6_repair import f2_solve
from v3_pair_shared import source_pair, hD
from v3_pair_untwisted import from_omega, to_omega


def from_diag(diag):
    def sign_value(indices):
        i, j = indices
        return sum(a.sigma for a in diag.rows[i:j]) % 2
    def component(second):
        def value(indices):
            i, j, k = indices
            return sum((-1) ** sign_value((i, r)) * edges[t]
                       for r in range(i, j)
                       for edges in ((diag.rows[r].edges_prime if second
                                      else diag.rows[r].edges),)
                       for t in range(j, k))
        return p.Cochain(2, value)
    return component(False), component(True), p.Cochain(1, sign_value)


def to_diag(A, Aprime, s, vertices):
    n = len(vertices) - 1
    def matrix_for(source):
        matrix = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                second = 0 if j == i + 1 else source((vertices[i], vertices[i + 1], vertices[j]))
                matrix[i][j] = source((vertices[i], vertices[i + 1], vertices[j + 1])) - second
                matrix[j][i] = -(-1) ** s((vertices[i], vertices[j])) * matrix[i][j]
        return matrix
    matrix, matrix_prime = matrix_for(A), matrix_for(Aprime)
    return rc.Diag2(tuple(rc.U1(s((vertices[i], vertices[i + 1])), tuple(matrix[i]),
                               tuple(matrix_prime[i])) for i in range(n)))


@lru_cache(None)
def source_value(pair):
    if rc.degree(pair) != 5:
        raise ValueError('source pair has degree five')
    A, Aprime, s = from_diag(pair[0])
    return source_pair(A, Aprime, s, from_omega(pair[1]))(tuple(range(6)))


def evaluate_source(chain):
    return sum(coefficient * source_value(pair) for pair, coefficient in chain.items()) % 2


@lru_cache(None)
def fixed_certificate():
    basis4, basis5 = rc.total_basis(4), rc.total_basis(5)
    indices = {basis: i for i, basis in enumerate(basis4)}
    rows, values = [], []
    for basis in basis5:
        boundary = rc.total_boundary(basis)
        rows.append(sum((coefficient % 2) << indices[face]
                        for face, coefficient in boundary.items() if face in indices))
        values.append(evaluate_source(rc.Gtot(basis)))
    rhs = sum(value << i for i, value in enumerate(values))
    solution = f2_solve(tuple(rows), rhs, len(basis4))
    return dict(basis4=basis4, basis5=basis5, rows=tuple(rows), values=tuple(values),
                coefficients=None if solution is None else
                tuple((solution >> i) & 1 for i in range(len(basis4))))


def primitive(A, Aprime, s, omega):
    if (A.degree, Aprime.degree, s.degree, omega.degree) != (2, 2, 1, 2):
        raise ValueError('V2 expects degrees (2,2,1,2)')
    certificate = fixed_certificate()
    if certificate['coefficients'] is None:
        raise ArithmeticError('fixed shared-source V2 has a nonzero small obstruction')
    coefficients = dict(zip(certificate['basis4'], certificate['coefficients']))
    def value(vertices):
        pair = (to_diag(A, Aprime, s, vertices), to_omega(omega, vertices))
        return (evaluate_source(rc.Htot(pair))
                + sum(coefficient * coefficients.get(basis, 0)
                      for basis, coefficient in rc.Ftot(pair).items())) % 2
    return p.Cochain(4, value)


def legal_beta(A, B, Aprime, Bprime, s, omega):
    """Degree-five beta in the matched production-upper calibration.

    The closed term rho(A) cup rho(Aprime) fixes the upper source periods.
    It changes no secondary boundary identity and vanishes on either A=0
    axis. It must also be used by the preceding off-shell successor prism.
    """
    alpha = hD(p.binary(A), p.binary(Aprime), s)
    Asum = A + Aprime
    Bsum = p.binary(B + Bprime + alpha)
    def splitting_carry(a, b):
        summand = p.source(a, s, omega)['g'] + p.cup(s, b)
        return p.binary(p.divide(summand - p.binary(summand), 2,
                                 'stacking secondary splitting carry'))
    base = p.binary(hD(B, Bprime, s) + hD(p.binary(B + Bprime), alpha, s)
                    + splitting_carry(A, B) + splitting_carry(Aprime, Bprime)
                    + splitting_carry(Asum, Bsum))
    return p.binary(base + primitive(A, Aprime, s, omega)
                    +p.cup(p.binary(A), p.binary(Aprime)))


if __name__ == '__main__':
    import json
    certificate = fixed_certificate()
    print(json.dumps(certificate, indent=2), flush=True)
