"""Universal source values shared by the page and extension workers.

Three universal functions of the fixed formulas are evaluated on registered
universal simplices and kept in the store of universal_values.py:

- the theta values of the degree-three source, one per pair of
  r3_source.to_diags (an r3_chain.Diag3 with skew U2 rows and a binary
  omega diagonal); every V3 sum reads them through r3_source.evaluate_r;
- the theta-pair sources of a0_high_gamma.HigherA0Stacking and its
  ThetaPairPhase subclass, one value per class, degree and universal tensor
  simplex, shared by every instance;
- the degree-two source low_phases.V2_pair.

The first two are keyed by an exact, invertible string encoding of their
universal simplex (a few hundred bytes each) instead of the nested key
objects, which the tables would otherwise retain in memory. The theta key
keeps the integer labels of the pair modulo four only: the source splitting
reads its integer cochain A through a = A mod 2 and the carry
t = ((A - a)/2) mod 2, every later cochain of the source and of Theta is a
function of those binary cochains, and the reconstruction from_diags reads
only the strict upper triangles, as signed rectangular sums, so the value is
a function of (sigma, upper labels mod 4, W) (doc/universal_value_growth.md,
Theorem 4). The theta table is the only one kept out of the bundled values:
one dense degree-three pair alone has tens of thousands of terms.

Copyright (c) 2026 koAHSS contributors; MIT license.
"""
from array import array
import base64
import functools

from universal_values import CompactTable, UniversalTable


def _pack(values):
    values = list(values)
    for code in ('b', 'i'):
        try:
            return code + base64.b64encode(array(code, values).tobytes()).decode('ascii')
        except OverflowError:
            continue
    raise OverflowError('universal simplex entries exceed 32 bits')


def _unpack(text):
    data = array(text[0])
    data.frombytes(base64.b64decode(text[1:]))
    return list(data)


def _bits(values):
    return format(sum(int(v) << i for i, v in enumerate(values)), 'x')


def _unbits(text, count):
    number = int(text, 16)
    return tuple((number >> i) & 1 for i in range(count))


# Names of the tables kept in the store only, never in data/universal-values.json.
UNBUNDLED = frozenset(['r3_source.phi_value'])


def compact_diag3_pair(pair):
    """Reduced string key of an (r3_chain.Diag3, chain_models.Diag) pair.

    The key holds each row's sigma, the strict upper triangle of its skew
    matrix modulo four and the strict upper triangle of the binary omega
    diagonal: the theta value of the pair is a function of these (module
    docstring). Two pairs with the same key have the same value; the pair
    that expand_diag3_pair rebuilds from a key is one of them.
    """
    diag, omega = pair
    n = len(diag.rows)
    if omega.mode != 'c2' or len(omega.rows) != n:
        raise ValueError('a degree-three source pair needs a binary omega diagonal of the same size')
    flat = []
    for row in diag.rows:
        flat.append(row.sigma % 2)
        flat.extend(row.matrix[i][j] % 4 for i in range(n) for j in range(i + 1, n))
    for i in range(n):
        line = omega.rows[i][1]
        flat.extend(line[j] % 2 for j in range(i + 1, n))
    return f'R3:{n}:{_pack(flat)}'


def expand_diag3_pair(text):
    """The representative pair of a reduced key: labels in 0..3, symmetric omega."""
    import chain_models
    import r3_chain
    tag, n, data = text.split(':', 2)
    if tag != 'R3':
        raise ValueError('not a degree-three source pair key')
    n = int(n)
    values = iter(_unpack(data))
    rows = []
    for _ in range(n):
        sigma = next(values)
        matrix = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                matrix[i][j] = next(values)
                matrix[j][i] = -matrix[i][j]
        rows.append(r3_chain.U2(sigma, tuple(map(tuple, matrix))))
    omega = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            omega[i][j] = omega[j][i] = next(values)
    return (r3_chain.Diag3(tuple(rows)),
            chain_models.Diag('c2', tuple((0, tuple(line)) for line in omega)))


def compact_theta_pair(key):
    """Exact string key of (class name, m, ((b, bp), (sign diagonal, omega)))."""
    name, m, ((b, bp), (sign, w)) = key
    if b.q != bp.q or b.q != w.q or len(sign.rows) != b.q or b.m != m or bp.m != m:
        raise ValueError('inconsistent theta-pair simplex')
    return ':'.join(('TP', name, str(m), str(b.q), _bits(b.values), _bits(bp.values),
                     _bits(row[0] for row in sign.rows), _bits(w.values)))


def expand_theta_pair(text):
    import a0_high_gamma
    import chain_models
    tag, name, m, q, b, bp, signs, w = text.split(':')
    if tag != 'TP':
        raise ValueError('not a theta-pair key')
    m, q = int(m), int(q)
    count = len(a0_high_gamma.faces(q, m))
    simplex = ((a0_high_gamma.KB(m, q, _unbits(b, count)), a0_high_gamma.KB(m, q, _unbits(bp, count))),
               (chain_models.Diag('zsign', tuple((s, (0,) * q) for s in _unbits(signs, q))),
                a0_high_gamma.KB(2, q, _unbits(w, len(a0_high_gamma.faces(q, 2))))))
    return (name, m, simplex)


def install_theta_values(tables, wrap=lambda function: function):
    """Route every V3 sum through a table of theta values per universal pair."""
    import r3_source
    table = CompactTable('r3_source.phi_value',
                         wrap(lambda pair: r3_source.phi()(r3_source.from_diags(pair))),
                         compact_diag3_pair, expand_diag3_pair)
    original = r3_source.evaluate_r

    @functools.wraps(original)
    def evaluate_r(cochain, chain):
        if cochain is r3_source.phi():
            return sum(c * table(pair) for pair, c in chain.items())
        return original(cochain, chain)
    r3_source.evaluate_r = evaluate_r
    tables[table.name] = table
    return table


def _theta_instance(name, m):
    if name == 'ThetaPairPhase':
        from theta_pair_phase import ThetaPairPhase
        return ThetaPairPhase(m)
    if name == 'HigherA0Stacking':
        import a0_high_gamma
        return a0_high_gamma.HigherA0Stacking(m + 2)
    raise KeyError(f'unknown theta-pair source class {name}')


def install_theta_pair_values(tables, wrap=lambda function: function):
    """One table of theta-pair source values for all HigherA0Stacking instances."""
    import a0_high_gamma
    cls = a0_high_gamma.HigherA0Stacking
    original = cls.evaluate
    raw = getattr(original, '__wrapped__', original)
    instances = {}

    def compute(key):
        name, m, simplex = key
        instance = instances.get((name, m))
        if instance is None:
            instance = instances.setdefault((name, m), _theta_instance(name, m))
        return raw(instance, simplex)
    table = CompactTable('a0_high_gamma.HigherA0Stacking.evaluate', wrap(compute),
                         compact_theta_pair, expand_theta_pair)

    def evaluate(self, simplex):
        key = (type(self).__qualname__, self.m, simplex)
        instances.setdefault(key[:2], self)
        return table(key)
    cls.evaluate = evaluate
    tables[table.name] = table
    return table


def install_v2_values(tables, wrap=lambda function: function):
    """The degree-two source, keyed like V1 by its pair of diagonal simplices."""
    import low_phases
    table = UniversalTable('low_phases.V2_pair', wrap(low_phases.V2_pair))
    low_phases.V2_pair = table
    tables[table.name] = table
    return table
