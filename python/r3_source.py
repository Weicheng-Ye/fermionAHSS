"""The fixed universal degree-three source primitive from r3_source.md.

Only this finite universal 46-by-94 matrix is reduced.  No matrix or
primitive is selected from the cochains of a supplied user space.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
import hashlib
import json
import gc
from pathlib import Path

try:
    from . import r3_chain as rc, chain_models as cm, phase_eval as pe
except ImportError:
    import r3_chain as rc
    import chain_models as cm
    import phase_eval as pe


def _u2_value(a, indices):
    i, j, k = indices
    return sum(a.matrix[r][t] for r in range(i, j) for t in range(j, k))


@lru_cache(None)
def from_diags(pair):
    signed, omega = pair
    n = len(signed.rows)
    if len(omega.rows) != n:
        raise ValueError("D3 source product needs equal dimensions")

    def sv(indices):
        i, j = indices
        return sum(a.sigma for a in signed.rows[i:j]) % 2

    def av(indices):
        i, j, k, l = indices
        return sum((-1)**sv((i, r)) * _u2_value(signed.rows[r], (j, k, l))
                   for r in range(i, j))

    def ov(indices):
        i, j, k = indices
        return sum(omega.rows[r][1][t] for r in range(i, j)
                   for t in range(j, k)) % 2

    return pe.make_source(3, n, sv, av, ov)


def to_diags(vertices):
    n = len(vertices) - 1
    A = pe.A3
    s = pe.s

    def abar(indices):
        if len(set(indices)) != len(indices):
            return 0
        ordered = tuple(sorted(indices))
        permutation_sign = (-1)**sum(indices[i] > indices[j]
                                      for i in range(4) for j in range(i + 1, 4))
        return permutation_sign * (-1)**s((vertices[0], vertices[ordered[0]])) * A(
            tuple(vertices[i] for i in ordered))

    rows = []
    for r in range(n):
        def z(i, j, k):
            return (-1)**s((vertices[0], vertices[r])) * (
                abar((r, i, j, k)) - abar((r + 1, i, j, k)))
        matrix = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                matrix[i][j] = z(i, i + 1, j + 1) - z(i, i + 1, j)
                matrix[j][i] = -matrix[i][j]
        rows.append(rc.U2(s((vertices[r], vertices[r + 1])),
                          tuple(map(tuple, matrix))))
    omega_matrix = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            omega_matrix[i][j] = (pe.omega((vertices[i], vertices[i + 1], vertices[j + 1]))
                                  + pe.omega((vertices[i], vertices[i + 1], vertices[j]))) % 2
            omega_matrix[j][i] = omega_matrix[i][j]
    return (rc.Diag3(tuple(rows)),
            cm.Diag('c2', tuple((0, tuple(row)) for row in omega_matrix)))


@lru_cache(None)
def phi():
    data = pe.source(pe.A3, pe.s, pe.omega)
    return pe.theta(data['p'], data['k'], pe.s, pe.omega)


def evaluate_r(cochain, chain):
    return sum(c * cochain(from_diags(pair)) for pair, c in chain.items())


def smith(matrix):
    """Deterministic Euclidean Smith recipe in r3_source.md §8.

    Returns D,U,V with U * matrix * V = D. Remainders are nonnegative;
    a nonzero remainder is immediately promoted to the pivot position.
    """
    a = [list(row) for row in matrix]
    m, n = len(a), len(a[0]) if a else 0
    u = [[int(i == j) for j in range(m)] for i in range(m)]
    v = [[int(i == j) for j in range(n)] for i in range(n)]

    def swap_rows(i, j):
        a[i], a[j] = a[j], a[i]
        u[i], u[j] = u[j], u[i]

    def swap_cols(i, j):
        for row in a:
            row[i], row[j] = row[j], row[i]
        for row in v:
            row[i], row[j] = row[j], row[i]

    def row_add(i, j, scale):
        a[i] = [x + scale * y for x, y in zip(a[i], a[j])]
        u[i] = [x + scale * y for x, y in zip(u[i], u[j])]

    def col_add(i, j, scale):
        for row in a:
            row[i] += scale * row[j]
        for row in v:
            row[i] += scale * row[j]

    for k in range(min(m, n)):
        entries = [(abs(a[i][j]), i, j) for i in range(k, m)
                   for j in range(k, n) if a[i][j]]
        if not entries:
            break
        _, i, j = min(entries)
        swap_rows(i, k)
        swap_cols(j, k)
        while True:
            if a[k][k] < 0:
                a[k] = [-x for x in a[k]]
                u[k] = [-x for x in u[k]]
            restart = False
            for i in range(k + 1, m):
                quotient, remainder = divmod(a[i][k], a[k][k])
                if quotient:
                    row_add(i, k, -quotient)
                if remainder:
                    swap_rows(i, k)
                    restart = True
                    break
            if restart:
                continue
            for j in range(k + 1, n):
                quotient, remainder = divmod(a[k][j], a[k][k])
                if quotient:
                    col_add(j, k, -quotient)
                if remainder:
                    swap_cols(j, k)
                    restart = True
                    break
            if restart:
                continue
            offender = next(((i, j) for i in range(k + 1, m)
                             for j in range(k + 1, n) if a[i][j] % a[k][k]), None)
            if offender is None:
                break
            row_add(k, offender[0], 1)
    return a, u, v


def _matmul(a, b):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def _rowmul(row, matrix):
    return [sum(x * y for x, y in zip(row, col)) for col in zip(*matrix)]


def _key(basis):
    word, omega = basis
    return [list(map(list, word)), list(omega)]


def _basis_key(value):
    word, omega = value
    return tuple(map(tuple, word)), tuple(omega)


def provenance():
    directory = Path(__file__).resolve().parent
    paths = [directory / name for name in
             ('r3_chain.py', 'r3_source.py', 'chain_models.py', 'exterior_bar.py',
              'phase_eval.py', 'cochain_tools.py')]
    paths.append(pe.PACKAGE / 'data/chi-calibrated-degree7-anf.g')
    return {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}


def clear_source_memos():
    """Release one universal period's finite simplex memo tables.

    SOURCE ids can be recycled only after EVERY cochain evaluator cache and
    this module's source-registration cache have been emptied. This helper
    is called solely by the isolated universal precompute command.
    """
    for obj in gc.get_objects():
        if isinstance(obj, pe.Cochain):
            obj.evaluate.cache_clear()
    from_diags.cache_clear()
    pe.SOURCES.clear()
    pe.SOURCE_KEYS.clear()
    gc.collect()


def precompute(progress=None):
    b7, b8 = rc.total_basis(7), rc.total_basis(8)
    positions = {b: i for i, b in enumerate(b7)}
    matrix = [[0] * len(b8) for _ in b7]
    for j, b in enumerate(b8):
        for t, coefficient in rc.total_boundary(b).items():
            if any(k for _, k in t[0]):
                matrix[positions[t]][j] += coefficient
    diagonal, u, v = smith(matrix)
    if _matmul(_matmul(u, matrix), v) != diagonal:
        raise ArithmeticError('universal Smith transformation failed')
    if any(diagonal[i][j] for i in range(len(b7)) for j in range(len(b8)) if i != j):
        raise ArithmeticError('Smith output is not diagonal')
    rank = sum(diagonal[i][i] != 0 for i in range(len(b7)))
    f = []
    source_phi = phi()
    for i, b in enumerate(b8):
        chain = rc.Gtot(b)
        if progress:
            progress('source_chain', index=i, total=len(b8), terms=len(chain))
        value = Fraction(evaluate_r(source_phi, chain)) % 1
        if (4 * value).denominator != 1:
            raise ArithmeticError('Theta5 has ceased to be quarter-valued')
        f.append(value)
        if progress:
            progress('source_value', index=i, value=str(value))
        clear_source_memos()
    transformed = [value % 1 for value in _rowmul(f, v)]
    if any(transformed[rank:]):
        raise ArithmeticError('nonzero universal R3 source cycle obstruction')
    ep = [Fraction(0)] * len(b7)
    for i in range(rank):
        d = diagonal[i][i]
        power = 0
        while d % 2 == 0:
            d //= 2
            power += 1
        ni = int(4 * transformed[i])
        mi = (pow(d, -1, 4) * ni) % 4
        ep[i] = Fraction(mi, 2**(power + 2))
    ef = _rowmul(ep, u)
    if any((a - b) % 1 for a, b in zip(_rowmul(ef, matrix), f)):
        raise ArithmeticError('fixed dyadic universal phase division failed')
    return dict(schema=1, source_hashes=provenance(),
                degree7_basis=list(map(_key, b7)), degree8_basis=list(map(_key, b8)),
                boundary=matrix, smith_diagonal=[diagonal[i][i] for i in range(rank)],
                smith_u=u, smith_v=v, source_values=list(map(str, f)),
                source_cycle_values=list(map(str, transformed[rank:])),
                ef_values=list(map(str, ef)))


def ef_from_data(data, check_hashes=True):
    if data.get('schema') != 1:
        raise ValueError('unsupported universal R3 source data')
    if check_hashes and data['source_hashes'] != provenance():
        raise ValueError('universal R3 source data hashes do not match current code')
    basis7 = tuple(map(_basis_key, data['degree7_basis']))
    if basis7 != rc.total_basis(7):
        raise ValueError('universal R3 source basis convention changed')
    values = tuple(map(Fraction, data['ef_values']))
    if len(values) != len(basis7):
        raise ValueError('universal R3 source data must supply every degree-seven coefficient')
    if any(value.denominator & (value.denominator - 1) for value in values):
        raise ValueError('universal R3 source data must be dyadic')
    return dict(zip(basis7, values))


def V3(vertices, ef):
    pair = to_diags(vertices)
    return (evaluate_r(phi(), rc.Htot(pair))
            + sum(c * ef.get(b, 0) for b, c in rc.Ftot(pair).items()))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    def progress(event, **fields):
        print(json.dumps(dict(event=event, **fields)), flush=True)
    data = precompute(progress)
    args.output.write_text(json.dumps(data, indent=2) + '\n')
    print(json.dumps(dict(event='completed', smith=data['smith_diagonal'],
                          ef=data['ef_values'])), flush=True)
