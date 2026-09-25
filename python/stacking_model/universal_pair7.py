"""Explicit primitive for a normalized closed signed pair source of degree7.

The universal input fibers are K(Z,3), with one common sign and omega.
Its relative degree-seven homology is zero integrally. Only ten fixed
source values enter the canonical Smith extension; no equation is solved
on the user's cochain complex.
"""
from fractions import Fraction as F
from functools import lru_cache
import v3_pair_shared as pair
from h_tau_primitive import r3s


def unimodular_inverse(matrix):
    n = len(matrix)
    a = [[F(x) for x in row]+[F(i == j) for j in range(n)]
         for i, row in enumerate(matrix)]
    for j in range(n):
        i = next(i for i in range(j, n) if a[i][j])
        a[j], a[i] = a[i], a[j]
        scale = a[j][j]
        a[j] = [x/scale for x in a[j]]
        for i in range(n):
            if i != j and a[i][j]:
                scale = a[i][j]
                a[i] = [x-scale*y for x, y in zip(a[i], a[j])]
    result = [row[n:] for row in a]
    if any(x.denominator != 1 for row in result for x in row):
        raise ArithmeticError('universal Smith transformation is not unimodular')
    return tuple(tuple(int(x) for x in row) for row in result)


@lru_cache(None)
def fixed_complex():
    lo, mid, hi = [pair.rc.total_basis(q) for q in (6, 7, 8)]
    matrix = tuple(tuple(pair.rc.total_boundary(b).get(a, 0) for a in lo) for b in mid)
    diagonal, u, v = r3s.smith([list(row) for row in matrix])
    rank = sum(bool(diagonal[i][i]) for i in range(min(len(mid), len(lo))))
    inverse_u = unimodular_inverse(u)
    boundary8 = tuple(tuple(pair.rc.total_boundary(b).get(a, 0) for b in hi) for a in mid)
    transformed = [[sum(inverse_u[j][i]*boundary8[j][k] for j in range(len(mid)))
                    for k in range(len(hi))] for i in range(len(mid))]
    if any(x for row in transformed[:rank] for x in row):
        raise ArithmeticError('the universal boundary does not land in the kernel')
    kernel_diagonal, _, _ = r3s.smith(transformed[rank:])
    if any(abs(kernel_diagonal[i][i]) != 1 for i in range(len(mid)-rank)):
        raise ArithmeticError('the relative universal degree-seven homology is nonzero')
    used = tuple(sorted({j for row in u[:rank] for j, x in enumerate(row) if x}))
    return dict(lower=lo, upper=mid, matrix=matrix, diagonal=diagonal,
                u=u, v=v, rank=rank, used=used,
                dimensions=(len(lo), len(mid), len(hi)),
                kernel_diagonal=tuple(kernel_diagonal[i][i] for i in range(len(mid)-rank)))


class SourcePrimitive7:
    """Source must be normalized, degree7, and signed-closed modulo integers."""
    def __init__(self, source):
        self.source = source

    @lru_cache(None)
    def source_value(self, simplex):
        A, Ap, s = pair.from_diag(simplex[0])
        if (not any(any(any(row) for row in a.matrix) for a in simplex[0].rows)
                or not any(any(any(row) for row in a.matrix_prime) for a in simplex[0].rows)):
            return F(0)
        value = self.source(A, Ap, s, pair.from_omega(simplex[1]))
        if value.degree != 7:
            raise ValueError('the universal source must have degree seven')
        return F(value(tuple(range(8))))

    def evaluate(self, chain):
        return sum(c*self.source_value(simplex) for simplex, c in chain.items())

    @lru_cache(None)
    def coefficients(self):
        data = fixed_complex()
        values = {j: self.evaluate(pair.rc.Gtot(data['upper'][j])) for j in data['used']}
        coordinates = [F(0)]*len(data['lower'])
        for i in range(data['rank']):
            rhs = sum(c*values[j] for j, c in enumerate(data['u'][i]) if c)
            coordinates[i] = (rhs % 1)/data['diagonal'][i][i]
        return tuple(sum(c*x for c, x in zip(row, coordinates)) % 1 for row in data['v'])

    def primitive(self, A, Ap, s, omega):
        if (A.degree, Ap.degree, s.degree, omega.degree) != (3, 3, 1, 2):
            raise ValueError('expected degrees (3,3,1,2)')
        small = dict(zip(fixed_complex()['lower'], self.coefficients()))
        def value(vertices):
            simplex = (pair.to_diag(A, Ap, s, vertices), pair.to_omega(omega, vertices))
            return self.evaluate(pair.rc.Htot(simplex))+sum(
                c*small.get(b, F(0)) for b, c in pair.rc.Ftot(simplex).items())
        return pair.p.Cochain(6, value)
