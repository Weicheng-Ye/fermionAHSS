"""Fixed signed degree-six pair-source primitive after checking all periods."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import json
import v2_pair_shared as pair
from h_tau_primitive import r3s
from universal_pair7 import unimodular_inverse


@lru_cache(None)
def fixed_complex():
    lo, mid, hi = [pair.rc.total_basis(q) for q in (5, 6, 7)]
    matrix = tuple(tuple(pair.rc.total_boundary(b).get(a, 0) for a in lo) for b in mid)
    diagonal, u, v = r3s.smith([list(row) for row in matrix])
    rank = sum(bool(diagonal[i][i]) for i in range(min(len(mid), len(lo))))
    inverse_u = unimodular_inverse(u)
    boundary = tuple(tuple(pair.rc.total_boundary(b).get(a, 0) for b in hi) for a in mid)
    transformed = [[sum(inverse_u[j][i]*boundary[j][k] for j in range(len(mid)))
                    for k in range(len(hi))] for i in range(len(mid))]
    if any(x for row in transformed[:rank] for x in row):
        raise ArithmeticError('the universal boundary does not land in the kernel')
    d, u2, _ = r3s.smith(transformed[rank:])
    r2 = sum(bool(d[i][i]) for i in range(min(len(d), len(d[0]))))
    u2_inverse = unimodular_inverse(u2)
    indices = [i for i in range(len(mid)-rank) if i >= r2 or abs(d[i][i]) != 1]
    cycles = tuple(tuple(sum(u[rank+j][i]*u2_inverse[j][h]
                        for j in range(len(mid)-rank)) for i in range(len(mid)))
                   for h in indices)
    orders = tuple(d[i][i] if i < r2 else 0 for i in indices)
    if orders != (2, 2, 0, 0):
        raise ArithmeticError('unexpected universal degree-six homology')
    used = tuple(sorted({j for row in u[:rank] for j, x in enumerate(row) if x}))
    return dict(lower=lo, upper=mid, matrix=matrix, diagonal=diagonal,
                u=u, v=v, rank=rank, used=used, cycles=cycles, orders=orders,
                dimensions=(len(lo), len(mid), len(hi)))


class SourcePrimitive6:
    """The source is closed modulo one; all four homology periods must vanish."""
    def __init__(self, source):
        self.source = source

    @lru_cache(None)
    def source_value(self, simplex):
        A, Ap, s = pair.from_diag(simplex[0])
        vertices = tuple(range(7))
        if (not any(A(f) for f in combinations(vertices, 3))
                or not any(Ap(f) for f in combinations(vertices, 3))):
            return F(0)
        value = self.source(A, Ap, s, pair.from_omega(simplex[1]))
        if value.degree != 6:
            raise ValueError('the source must have degree six')
        return F(value(vertices))

    def evaluate(self, chain):
        return sum(c*self.source_value(simplex) for simplex, c in chain.items())

    @lru_cache(None)
    def small_value(self, index):
        return self.evaluate(pair.rc.Gtot(fixed_complex()['upper'][index]))

    @lru_cache(None)
    def coefficients(self):
        data = fixed_complex()
        periods = tuple(sum(c*self.small_value(j) for j, c in enumerate(row) if c)
                        for row in data['cycles'])
        if any(x % 1 for x in periods):
            raise ArithmeticError(f'the universal source has nonzero periods: {periods}')
        coordinates = [F(0)]*len(data['lower'])
        for i in range(data['rank']):
            rhs = sum(c*self.small_value(j) for j, c in enumerate(data['u'][i]) if c)
            coordinates[i] = (rhs % 1)/data['diagonal'][i][i]
        return tuple(sum(c*x for c, x in zip(row, coordinates)) % 1 for row in data['v'])

    def primitive(self, A, Ap, s, omega):
        if (A.degree, Ap.degree, s.degree, omega.degree) != (2, 2, 1, 2):
            raise ValueError('expected degrees (2,2,1,2)')
        small = dict(zip(fixed_complex()['lower'], self.coefficients()))
        def value(vertices):
            simplex = (pair.to_diag(A, Ap, s, vertices), pair.to_omega(omega, vertices))
            return self.evaluate(pair.rc.Htot(simplex))+sum(
                c*small.get(b, F(0)) for b, c in pair.rc.Ftot(simplex).items())
        return pair.p.Cochain(5, value)
