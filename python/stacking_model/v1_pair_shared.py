"""Shared-sign degree-one source-pair contraction and fixed calibration."""
from functools import lru_cache
from pathlib import Path
import sys

sys.path.insert(0, str(Path(str(Path(__file__).resolve().parents[1]))))
import phase_eval as p
import r1_pair_chain as rc
from g6_repair import f2_solve
from v3_pair_shared import source_pair, hD
from v3_pair_untwisted import from_omega, to_omega


def from_diag(diag):
    def sign_value(indices):
        i, j = indices
        return sum(diag.signs[i:j]) % 2
    def component(second):
        edges = diag.right if second else diag.left
        return p.Cochain(1, lambda f: (-1) ** sum(diag.signs[:f[0]])
                         * sum(edges[f[0]:f[1]]))
    return component(False), component(True), p.Cochain(1, sign_value)


def to_diag(A, Aprime, s, vertices):
    n = len(vertices) - 1
    signs = tuple(s((vertices[i], vertices[i + 1])) for i in range(n))
    def edges(source):
        return tuple((-1) ** sum(signs[:i]) * source((vertices[i], vertices[i + 1]))
                     for i in range(n))
    return rc.Borel(signs, edges(A), edges(Aprime))


@lru_cache(None)
def source_value(pair):
    if rc.degree(pair) != 4:
        raise ValueError('source pair has degree four')
    A, Aprime, s = from_diag(pair[0])
    return source_pair(A, Aprime, s, from_omega(pair[1]))(tuple(range(5)))


def evaluate_source(chain):
    return sum(coefficient * source_value(pair) for pair, coefficient in chain.items()) % 2


@lru_cache(None)
def fixed_certificate():
    basis3, basis4 = rc.total_basis(3), rc.total_basis(4)
    indices = {basis: i for i, basis in enumerate(basis3)}
    rows, values = [], []
    for basis in basis4:
        boundary = rc.total_boundary(basis)
        rows.append(sum((coefficient % 2) << indices[face]
                        for face, coefficient in boundary.items() if face in indices))
        values.append(evaluate_source(rc.Gtot(basis)))
    rhs = sum(value << i for i, value in enumerate(values))
    solution = f2_solve(tuple(rows), rhs, len(basis3))
    return dict(basis3=basis3, basis4=basis4, rows=tuple(rows), values=tuple(values),
                coefficients=None if solution is None else
                tuple((solution >> i) & 1 for i in range(len(basis3))))


def primitive(A, Aprime, s, omega):
    if (A.degree, Aprime.degree, s.degree, omega.degree) != (1, 1, 1, 2):
        raise ValueError('V1 expects degrees (1,1,1,2)')
    certificate = fixed_certificate()
    if certificate['coefficients'] is None:
        raise ArithmeticError('fixed shared-source V1 has a nonzero small obstruction')
    coefficients = dict(zip(certificate['basis3'], certificate['coefficients']))
    def value(vertices):
        pair = (to_diag(A, Aprime, s, vertices), to_omega(omega, vertices))
        return (evaluate_source(rc.Htot(pair))
                + sum(coefficient * coefficients.get(basis, 0)
                      for basis, coefficient in rc.Ftot(pair).items())) % 2
    return p.Cochain(3, value)


def legal_beta(A, B, Aprime, Bprime, s, omega):
    """Degree-four beta on legal first-two-layer data, independent of C."""
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
    return p.binary(base + primitive(A, Aprime, s, omega))


if __name__ == '__main__':
    import json
    certificate = fixed_certificate()
    print(json.dumps(certificate, indent=2), flush=True)
