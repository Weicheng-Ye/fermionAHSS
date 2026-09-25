"""Nonlocal finite-model repair of the degree-six upper-cutoff extension.

The injected P, tau, T and J are the fixed production operations, not
unknown callbacks to be solved. All projections and defining cochain choices
are computed here. This module does not implement a full stacking product.

Binary cochains are nonnegative bit-vectors; binary matrices are tuples of
bit-vector rows. Integral matrices are tuples of tuple rows, acting on
column vectors. No third-party dependency is used.
"""
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache


def parity(x):
    return x.bit_count() & 1


def f2_apply(rows, vector):
    return sum(parity(row & vector) << i for i, row in enumerate(rows))


def f2_rows(columns, dimension):
    return tuple(sum(((column >> row) & 1) << j
                     for j, column in enumerate(columns))
                 for row in range(dimension))


def f2_rref(rows, ncols):
    data = list(rows)
    pivot_columns = []
    position = 0
    for col in range(ncols):
        found = next((i for i in range(position, len(data))
                      if (data[i] >> col) & 1), None)
        if found is None:
            continue
        data[position], data[found] = data[found], data[position]
        for i in range(len(data)):
            if i != position and ((data[i] >> col) & 1):
                data[i] ^= data[position]
        pivot_columns.append(col)
        position += 1
    return tuple(data), tuple(pivot_columns)


def f2_solve(rows, target, ncols):
    if any(row < 0 or row >> ncols for row in rows):
        raise ValueError('binary matrix row is outside its declared dimension')
    if target < 0 or target >> len(rows):
        raise ValueError('binary target is outside its declared dimension')
    augmented = tuple(row | (((target >> i) & 1) << ncols)
                      for i, row in enumerate(rows))
    reduced, pivots = f2_rref(augmented, ncols)
    mask = (1 << ncols) - 1
    if any((row & mask) == 0 and ((row >> ncols) & 1)
           for row in reduced):
        return None
    # Every free variable is zero; RREF supplies the pivot values.
    return sum(((reduced[i] >> ncols) & 1) << col
               for i, col in enumerate(pivots))


def f2_kernel(rows, ncols):
    reduced, pivots = f2_rref(rows, ncols)
    free = [j for j in range(ncols) if j not in pivots]
    result = []
    for col in free:
        vector = 1 << col
        for i, pivot in enumerate(pivots):
            vector |= ((reduced[i] >> col) & 1) << pivot
        result.append(vector)
    return tuple(result)


def xor_columns(columns, coefficients):
    value = 0
    for i, column in enumerate(columns):
        if (coefficients >> i) & 1:
            value ^= column
    return value


class BinarySubspace:
    """A subspace with a fully computed complement and exact coordinates."""
    def __init__(self, dimension, generators):
        self.dimension = dimension
        selected = []
        for vector in generators:
            if vector < 0 or vector >> dimension:
                raise ValueError('binary generator is outside its dimension')
            if f2_solve(f2_rows(selected, dimension), vector,
                        len(selected)) is None:
                selected.append(vector)
        self.basis = tuple(selected)
        for j in range(dimension):
            vector = 1 << j
            if f2_solve(f2_rows(selected, dimension), vector,
                        len(selected)) is None:
                selected.append(vector)
        self.full_basis = tuple(selected)
        rows = f2_rows(selected, dimension)
        # Inverse columns are the coordinates of the standard vectors.
        inverse_columns = tuple(f2_solve(rows, 1 << j, dimension)
                                for j in range(dimension))
        if any(value is None for value in inverse_columns):
            raise AssertionError('basis completion failed')
        self.inverse_rows = f2_rows(inverse_columns, dimension)

    def coordinates(self, vector):
        if vector < 0 or vector >> self.dimension:
            raise ValueError('binary vector is outside its dimension')
        return f2_apply(self.inverse_rows, vector)

    def project(self, vector):
        return xor_columns(self.basis, self.coordinates(vector))

    def quotient(self, vector):
        return self.coordinates(vector) >> len(self.basis)


def integral_vector(values, name):
    result = []
    for value in values:
        if isinstance(value, int):
            result.append(value)
        elif isinstance(value, Fraction) and value.denominator == 1:
            result.append(value.numerator)
        else:
            raise ArithmeticError(f'{name} returned a nonintegral cochain')
    return tuple(result)


def int_apply(matrix, vector):
    return tuple(sum(a * b for a, b in zip(row, vector))
                 for row in matrix)


def int_product(left, right):
    if not left:
        return ()
    ncols = len(right[0]) if right else 0
    return tuple(tuple(sum(left[i][j] * right[j][k]
                           for j in range(len(right)))
                       for k in range(ncols))
                 for i in range(len(left)))


def rational_rank(matrix, ncols):
    data = [[Fraction(value) for value in row] for row in matrix]
    position = 0
    for col in range(ncols):
        found = next((i for i in range(position, len(data)) if data[i][col]), None)
        if found is None:
            continue
        data[position], data[found] = data[found], data[position]
        scale = data[position][col]
        data[position] = [value / scale for value in data[position]]
        for i in range(len(data)):
            if i != position and data[i][col]:
                coefficient = data[i][col]
                data[i] = [a - coefficient * b
                           for a, b in zip(data[i], data[position])]
        position += 1
    return position


@dataclass(frozen=True)
class IntegralCycleCoordinates:
    """Certified integral cycle projection from a unimodular decomposition.

    basis has cycle columns first, followed by complement columns. Its
    integer inverse is checked in both orders. Injectivity of delta on the
    complement proves that the first block is the entire integral kernel.
    """
    delta: tuple
    basis: tuple
    inverse: tuple
    cycle_rank: int

    def __post_init__(self):
        n = len(self.basis)
        if not 0 <= self.cycle_rank <= n:
            raise ValueError('invalid cycle rank')
        matrices = (self.basis, self.inverse)
        if any(len(matrix) != n or any(len(row) != n for row in matrix)
               for matrix in matrices):
            raise ValueError('basis and inverse must be square matrices')
        if any(len(row) != n for row in self.delta):
            raise ValueError('integral delta has incompatible dimensions')
        if any(not isinstance(x, int)
               for matrix in (self.delta,) + matrices
               for row in matrix for x in row):
            raise ValueError('cycle coordinates must be integral')
        identity = tuple(tuple(int(i == j) for j in range(n)) for i in range(n))
        if (int_product(self.basis, self.inverse) != identity
                or int_product(self.inverse, self.basis) != identity):
            raise ValueError('basis inverse is not an integer inverse')
        transformed = int_product(self.delta, self.basis)
        if any(value for row in transformed for value in row[:self.cycle_rank]):
            raise ValueError('declared cycle columns are not cycles')
        complement = tuple(row[self.cycle_rank:] for row in transformed)
        if rational_rank(complement, n - self.cycle_rank) != n - self.cycle_rank:
            raise ValueError('declared cycle block does not span the kernel')

    @property
    def dimension(self):
        return len(self.basis)

    def coordinates(self, vector):
        if len(vector) != self.dimension or any(not isinstance(x, int) for x in vector):
            raise ValueError('A has incompatible integral coordinates')
        return int_apply(self.inverse, vector)[:self.cycle_rank]

    def from_cycles(self, coefficients):
        if len(coefficients) != self.cycle_rank:
            raise ValueError('wrong number of integral cycle coordinates')
        return int_apply(self.basis, tuple(coefficients)
                         + (0,) * (self.dimension - self.cycle_rank))

    def project(self, vector):
        return self.from_cycles(self.coordinates(vector))



class IntegerParityKernel:
    """Integral lattice ker(M:Z^n->F2^r), with a fixed coset section.

    A free-column basis vector is e_f + sum M_rref[p,f] e_p; each
    pivot contributes 2e_p. Subtracting the parity remainder at pivots
    gives an exact lattice decomposition without a Smith algorithm.
    """
    def __init__(self, rows, dimension):
        if any(row < 0 or row >> dimension for row in rows):
            raise ValueError('parity matrix has incompatible dimensions')
        self.dimension = dimension
        reduced, pivots = f2_rref(rows, dimension)
        self.rows = reduced[:len(pivots)]
        self.pivots = pivots
        self.free = tuple(j for j in range(dimension) if j not in pivots)
        columns = []
        for col in self.free:
            vector = [int(i == col) for i in range(dimension)]
            for row, pivot in zip(self.rows, self.pivots):
                vector[pivot] = (row >> col) & 1
            columns.append(tuple(vector))
        for pivot in self.pivots:
            columns.append(tuple(2 * int(i == pivot) for i in range(dimension)))
        self.columns = tuple(columns)

    def apply(self, coordinates):
        if len(coordinates) != self.dimension:
            raise ValueError('parity lattice coordinate dimension mismatch')
        return tuple(sum(column[i] * value
                         for column, value in zip(self.columns, coordinates))
                     for i in range(self.dimension))

    def split(self, vector):
        """Return (kernel coordinates, remainder); input=K(coords)+remainder."""
        if len(vector) != self.dimension or any(not isinstance(x, int) for x in vector):
            raise ValueError('parity lattice input must be an integral vector')
        residue = sum((value % 2) << i for i, value in enumerate(vector))
        remainder = [0] * self.dimension
        for row, pivot in zip(self.rows, self.pivots):
            remainder[pivot] = parity(row & residue)
        kernel_vector = tuple(value - carry for value, carry in zip(vector, remainder))
        free_values = tuple(kernel_vector[i] for i in self.free)
        pivot_values = []
        for row, pivot in zip(self.rows, self.pivots):
            numerator = kernel_vector[pivot] - sum(((row >> col) & 1) * value
                                                   for col, value in zip(self.free, free_values))
            if numerator % 2:
                raise AssertionError('parity lattice remainder was not removed')
            pivot_values.append(numerator // 2)
        coordinates = free_values + tuple(pivot_values)
        if tuple(a + b for a, b in zip(self.apply(coordinates), remainder)) != tuple(vector):
            raise AssertionError('parity lattice decomposition failed')
        return coordinates, tuple(remainder)

    def project(self, vector):
        coordinates, _ = self.split(vector)
        return self.apply(coordinates)


class G6Repair:
    """Explicit finite-model retraction and repaired upper correction.

    delta_b has target dimension dim_c. delta_c has source dimension dim_c.
    qd_closed is QD on closed degree-four binary cochains. tau is only
    evaluated on valid first-two-layer data; tertiary is only evaluated on
    full defining systems. j_legal is the existing J6 on the larger L6
    domain, where C is arbitrary. delta_top is delta_s:C^8->C^9.

    The operation contracts are the production identities: QD is additive
    on cohomology; tau changes by QD under closed B shifts; induced Tau is
    additive on primary-admissible A modulo the complete B indeterminacy;
    and tau depends on a closed integral A only modulo four. The returned
    defining equations are checked directly before use; those checks do not
    independently establish all universal operation identities.
    """
    def __init__(self, a_coordinates, dim_b, dim_c, delta_b, delta_c,
                 delta_top, primary, tau, qd_closed, tertiary, j_legal):
        self.a_coordinates = a_coordinates
        self.dim_b, self.dim_c = dim_b, dim_c
        self.delta_b, self.delta_c = tuple(delta_b), tuple(delta_c)
        self.delta_top = tuple(tuple(row) for row in delta_top)
        self.primary, self.tau = primary, tau
        self.qd_closed = qd_closed
        self.tertiary, self.j_legal = tertiary, j_legal
        if len(self.delta_b) != dim_c:
            raise ValueError('delta B must target the degree-five C space')
        if any(row < 0 or row >> dim_b for row in self.delta_b):
            raise ValueError('delta B row has incompatible dimensions')
        if any(row < 0 or row >> dim_c for row in self.delta_c):
            raise ValueError('delta C row has incompatible dimensions')
        if any(f2_apply(self.delta_c, f2_apply(self.delta_b, 1 << i))
               for i in range(dim_b)):
            raise ValueError('binary coboundaries do not square to zero')
        self.closed_b = f2_kernel(self.delta_b, dim_b)
        self.closed_c = BinarySubspace(dim_c, f2_kernel(self.delta_c, dim_c))
        boundary_columns = tuple(f2_apply(self.delta_c, 1 << i)
                                 for i in range(dim_c))
        self.boundaries = BinarySubspace(len(self.delta_c), boundary_columns)
        qd_columns = tuple(self.boundaries.quotient(qd_closed(b))
                           for b in self.closed_b)
        quotient_dim = len(self.delta_c) - len(self.boundaries.basis)
        self.qd_rows = f2_rows(qd_columns, quotient_dim)
        coefficient_kernel = f2_kernel(self.qd_rows, len(self.closed_b))
        self.allowed_b = BinarySubspace(
            dim_b, (xor_columns(self.closed_b, z) for z in coefficient_kernel))
        zero_a = (0,) * a_coordinates.dimension
        if primary(zero_a) or tau(zero_a, 0):
            raise ValueError('production lower operations must be pointed')
        zero_j = integral_vector(j_legal(zero_a, 0, 0), 'J')
        if any(zero_j):
            raise ValueError('the existing J correction must be pointed')
        self.top_dimension = len(zero_j)
        if any(len(row) != self.top_dimension for row in self.delta_top):
            raise ValueError('top differential has incompatible dimensions')
        # First admissibility lattice: primary obstruction modulo im delta B.
        primary_boundaries = BinarySubspace(
            dim_c, (f2_apply(self.delta_b, 1 << i) for i in range(dim_b)))
        rank_a = a_coordinates.cycle_rank
        primary_columns = []
        for i in range(rank_a):
            coefficients = tuple(int(i == j) for j in range(rank_a))
            cycle = a_coordinates.from_cycles(coefficients)
            primary_columns.append(primary_boundaries.quotient(primary(cycle)))
        primary_rows = f2_rows(primary_columns, dim_c - len(primary_boundaries.basis))
        self.primary_lattice = IntegerParityKernel(primary_rows, rank_a)
        # Second admissibility lattice: Tau modulo its complete B indeterminacy.
        indeterminacy = BinarySubspace(
            len(self.delta_c), boundary_columns + tuple(qd_closed(b) for b in self.closed_b))
        secondary_columns = []
        for coefficients in self.primary_lattice.columns:
            cycle = a_coordinates.from_cycles(coefficients)
            b0 = f2_solve(self.delta_b, primary(cycle), dim_b)
            if b0 is None:
                raise ArithmeticError('primary operation is not additive on cohomology')
            secondary_columns.append(indeterminacy.quotient(tau(cycle, b0)))
        secondary_rows = f2_rows(secondary_columns,
                                len(self.delta_c) - len(indeterminacy.basis))
        self.secondary_lattice = IntegerParityKernel(secondary_rows, rank_a)

    def in_locus(self, a, b):
        self.a_coordinates.coordinates(a)
        if b < 0 or b >> self.dim_b:
            raise ValueError('B is outside its binary dimension')
        return (not any(int_apply(self.a_coordinates.delta, a))
                and f2_apply(self.delta_b, b) == self.primary(a))

    @lru_cache(None)
    def reference_b(self, residue):
        """Return B0 for this residue, or None if no B,C exist."""
        a = self.a_coordinates.from_cycles(residue)
        b = f2_solve(self.delta_b, self.primary(a), self.dim_b)
        if b is None:
            return None
        t = self.tau(a, b)
        coefficients = f2_solve(self.qd_rows, self.boundaries.quotient(t),
                                len(self.closed_b))
        if coefficients is None:
            return None
        b ^= xor_columns(self.closed_b, coefficients)
        if f2_solve(self.delta_c, self.tau(a, b), self.dim_c) is None:
            raise ArithmeticError('Tau defining-cochain variation contract failed')
        return b

    def retract(self, a, b, c):
        a = tuple(a)
        if b < 0 or b >> self.dim_b or c < 0 or c >> self.dim_c:
            raise ValueError('binary input is outside its dimension')
        coordinates = self.a_coordinates.coordinates(a)
        primary_coordinates, _ = self.primary_lattice.split(coordinates)
        secondary_coordinates, _ = self.secondary_lattice.split(primary_coordinates)
        coefficients = self.primary_lattice.apply(
            self.secondary_lattice.apply(secondary_coordinates))
        residue = tuple(value % 4 for value in coefficients)
        b0 = self.reference_b(residue)
        if b0 is None:
            raise ArithmeticError('secondary operation is not additive modulo indeterminacy')
        a_star = self.a_coordinates.from_cycles(coefficients)
        b_star = b0 ^ self.allowed_b.project(b ^ b0)
        if not self.in_locus(a_star, b_star):
            raise ArithmeticError('retraction failed the primary equation')
        t = self.tau(a_star, b_star)
        c0 = f2_solve(self.delta_c, t, self.dim_c)
        if c0 is None:
            raise ArithmeticError('retraction failed Tau modulo-four/variation contract')
        c_star = c0 ^ self.closed_c.project(c ^ c0)
        if f2_apply(self.delta_c, c_star) != t:
            raise ArithmeticError('retraction failed the secondary equation')
        return a_star, b_star, c_star

    def correction(self, a, b, c):
        a = tuple(a)
        if c < 0 or c >> self.dim_c:
            raise ValueError('C is outside its binary dimension')
        if self.in_locus(a, b):
            result = integral_vector(self.j_legal(a, b, c), 'J')
            if len(result) != self.top_dimension:
                raise ValueError('J returned incompatible top coordinates')
            return result
        legal = self.retract(a, b, c)
        result = integral_vector(self.tertiary(*legal), 'T')
        if len(result) != self.top_dimension:
            raise ValueError('T returned incompatible top coordinates')
        if any(int_apply(self.delta_top, result)):
            raise ArithmeticError('T on a full defining system is not closed')
        return result
