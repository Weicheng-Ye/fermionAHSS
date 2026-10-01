"""Compile the degree-six unary coefficient table data/unary-gamma6-coefficients.json.

The 22 constants q_b of unary_gamma6, one per relative unary small cell b
of r3_chain.total_basis(6), are

    q_b = S(H_P Delta_* G b) + c(F_P Delta_* G b)  modulo Z_(2),

with G the unary contraction map of r3_chain, Delta_* the diagonal
U2(sigma,M) -> U2(sigma,M,M), (F_P, H_P) the pair contraction of
r3_pair_chain, c = production_gamma6.FIXED_COEFFICIENTS on the pair small
cells and S the pair source production_gamma6.source, zero on a simplex
with a zero A fiber. Each cell is evaluated exactly. The source is split
literally as S = -V(A) - V(A') + V(A+A') + L(A,A'), with V the
A-only source primitive of closed_a_upper and L the rest of the source;
each V value is taken on the cocycle section of the reconstructed unary
simplex (j composed with r), never on the raw universal simplex.

Every computed cell is reported as one JSON line (also appended to the
--record file); --collect reads such lines back. The table is written only
when all 22 cells are known, every cell was computed with the current
formula sources, equal cells agree, every entry has a denominator dividing
48 and, with --reference, every entry equals the reference table. Run it
after any change to the formula sources, which invalidates the table.

usage: python3 python/generate_unary_coefficients.py [--jobs N]
       python3 python/generate_unary_coefficients.py --cells 3 11 --record FILE --no-write
       python3 python/generate_unary_coefficients.py --collect FILE ... [--reference TABLE]

The universal-value cache is not used unless --cache names a store
directory; the bundled universal values are used unless
FERMIONAHSS_BUNDLED_VALUES=0.

Copyright (c) 2026 koAHSS contributors; MIT license.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import gc
import json
import multiprocessing
import os
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import extension_transfer  # noqa: E402,F401  (installs the worker's exact policy)
import extension_acceleration as acceleration  # noqa: E402
import universal_values as universal  # noqa: E402
import unary_gamma6  # noqa: E402
import r3_chain as unary  # noqa: E402
import r3_pair_chain as pair  # noqa: E402
import production_gamma6 as gamma  # noqa: E402
import v3_pair_shared  # noqa: E402
from universal_pair7 import one_fiber_zero  # noqa: E402

COMPILER = 'python/generate_unary_coefficients.py'
DEFINITION = ('q_b = (Delta^* P_S)(G b) mod Z_(2); light-transport-proof (20), '
              'light-completion-tertiary (4.3)-(4.5)')
p = gamma.p
aux = gamma.aux
upper = aux.upper


def open_store(cache=None):
    """Attach the universal-value store: the bundled values and, if given, a cache."""
    previous = os.environ.get('FERMIONAHSS_CACHE_DIR')
    os.environ['FERMIONAHSS_CACHE_DIR'] = '' if cache is None else str(cache)
    try:
        acceleration.persist()
    finally:
        if previous is None:
            del os.environ['FERMIONAHSS_CACHE_DIR']
        else:
            os.environ['FERMIONAHSS_CACHE_DIR'] = previous


def diagonal(x):
    """Delta_* of a unary cell: the same matrix in both fibers."""
    return pair.Diag3(tuple(pair.U2(a.sigma, a.matrix, a.matrix) for a in x.rows))


def local_source(A, Ap, s, omega):
    """The pair source without its A-only primitives: S + V(A) + V(A') - V(A+A')."""
    n = A.degree
    B = p.zero(n+1)
    alpha = upper.hp.hD(p.binary(A), p.binary(Ap), s)
    beta = upper.beta_sharp(A, B, Ap, B, s, omega)
    P, Pp = upper.primary(A, s, omega), upper.primary(Ap, s, omega)
    k, kp = upper.fsharp(A, B, s, omega), upper.fsharp(Ap, B, s, omega)
    theta = aux.theta_phase(n).phase(*(aux.as_cochain(x) for x in (P, k, Pp, kp, s, omega)))
    return (p.scale(aux.phi_without_A_source(A+Ap, alpha, beta, s, omega), -1)
            - p.Cochain(theta.degree, theta)
            + gamma.comparison.A_phase(A, s, omega) + gamma.comparison.A_phase(Ap, s, omega)
            - gamma.comparison.A_phase(A+Ap, s, omega)
            + p.half(gamma.binary.affine_J0_source(A, Ap, s, omega)))


def local_value(simplex, memo):
    """L on one universal pair simplex, with the universal evaluation policy."""
    try:
        return memo[simplex]
    except KeyError:
        pass
    with acceleration.universal():
        A, Ap, s = v3_pair_shared.from_diag(simplex[0])
        value = local_source(A, Ap, s, v3_pair_shared.from_omega(simplex[1]))(tuple(range(8)))
    memo[simplex] = value
    return value


def source_on_chain(chain):
    """S on a pair chain, split into the A-only primitives and the rest."""
    filtered = {x: c for x, c in chain.items() if not one_fiber_zero(x)}
    projected = {}
    for (x, w), c in filtered.items():
        for choice, factor in (('left', -1), ('right', -1), ('sum', 1)):
            rows = []
            for a in x.rows:
                if choice == 'left':
                    matrix = a.matrix
                elif choice == 'right':
                    matrix = a.matrix_prime
                else:
                    matrix = tuple(tuple(l+r for l, r in zip(lrow, rrow))
                                   for lrow, rrow in zip(a.matrix, a.matrix_prime))
                rows.append(unary.U2(a.sigma, matrix))
            # The cocycle section of the reconstructed unary simplex (j o r).
            A0, s0, w0 = upper.hp.from_pair((unary.Diag3(tuple(rows)), w), 3)
            canonical = upper.hp.to_pair(A0, s0, w0, tuple(range(8)))
            unary.add(projected, {canonical: c*factor})
    value = sum(c*upper._high_source_value(x) for x, c in projected.items())
    memo = {}
    for x, c in filtered.items():
        value += c*local_value(x, memo)
    return value


def cell_value(index):
    """The exact rational (Delta^* P_S)(G b) of the cell b of the given index."""
    cell = unary_gamma6.unary_cells()[index]
    hchain, fchain = {}, {}
    for (x, w), c in unary.Gtot(cell).items():
        d = (diagonal(x), w)
        pair.add(hchain, pair.Htot(d), c)
        pair.add(fchain, pair.Ftot(d), c)
    small = dict(zip(pair.total_basis(6), gamma.FIXED_COEFFICIENTS))
    return source_on_chain(hchain) + sum(c*small.get(b, Fraction(0)) for b, c in fchain.items())


def compile_cell(index, cache=None):
    """One cell as a record line: exact value, value modulo one, provenance."""
    open_store(cache)
    gc.disable()
    started = time.monotonic()
    try:
        value = cell_value(index)
    finally:
        acceleration.flush()
        gc.enable()
        gc.collect()
    classes = universal._key_classes()
    return {'index': index, 'cell': universal._encode(unary_gamma6.unary_cells()[index], classes),
            'value': str(value), 'mod1': str(value % 1),
            'seconds': round(time.monotonic()-started, 1),
            'sourceHash': unary_gamma6.source_hash(), 'compiler': COMPILER}


def read_records(paths):
    records = []
    for path in paths:
        for line in Path(path).read_text().splitlines():
            try:
                record = json.loads(line)
            except ValueError:
                continue
            if isinstance(record, dict) and 'index' in record and 'mod1' in record:
                records.append(record)
    return records


def reference_values(path):
    """q_mod1 and, where recorded, the exact values of a reference table."""
    data = json.loads(Path(path).read_text())
    q = [Fraction(v) for v in data['q_mod1']]
    exact = {}
    for key, entry in data.get('raw', {}).items():
        exact[int(key)] = Fraction(entry['value'] if isinstance(entry, dict) else entry)
    return q, exact


def assemble(records, reference=None):
    """The table record, or SystemExit naming every reason it cannot be written."""
    problems = []
    provenance = unary_gamma6.source_hash()
    cells = unary_gamma6.encoded_cells()
    values = {}
    for record in records:
        index = record['index']
        if not isinstance(index, int) or not 0 <= index < len(cells):
            problems.append(f'record with an invalid index {index!r}')
            continue
        if record.get('sourceHash') != provenance:
            problems.append(f'cell {index} was computed from other formula sources')
            continue
        if record.get('cell') != cells[index]:
            problems.append(f'cell {index} was computed on another cell')
            continue
        value = Fraction(record['value'])
        if Fraction(record['mod1']) != value % 1:
            problems.append(f'cell {index} records an inconsistent value modulo one')
            continue
        if index in values and values[index] != value:
            problems.append(f'cell {index} has conflicting values {values[index]} and {value}')
            continue
        values[index] = value
    missing = [i for i in range(len(cells)) if i not in values]
    if missing:
        problems.append('uncomputed cells: ' + ','.join(map(str, missing)))
    for index, value in sorted(values.items()):
        if unary_gamma6.DENOMINATOR % (value % 1).denominator:
            problems.append(f'cell {index} has a denominator not dividing '
                            f'{unary_gamma6.DENOMINATOR}: {value % 1}')
    cross = None
    if reference is not None and not missing:
        q, exact = reference_values(reference)
        differ = [i for i in range(len(cells)) if q[i] != values[i] % 1]
        if differ:
            problems.append('cells differing from the reference modulo one: '
                            + ', '.join(f'{i} ({values[i] % 1} != {q[i]})' for i in differ))
        raw_differ = [i for i in sorted(exact) if exact[i] != values[i]]
        if raw_differ:
            print('note: exact values differing from the reference by integers: '
                  + ', '.join(f'{i} ({values[i]} != {exact[i]})' for i in raw_differ),
                  file=sys.stderr)
        try:
            shown = str(Path(reference).resolve().relative_to(ROOT))
        except ValueError:
            shown = str(reference)
        cross = {'reference': shown, 'agree': len(cells)-len(differ),
                 'agreeExact': len(exact)-len(raw_differ)}
    if problems:
        raise SystemExit('the unary coefficient table is not written:\n  ' + '\n  '.join(problems))
    table = {'schema': unary_gamma6.SCHEMA, 'definition': DEFINITION, 'basis': unary_gamma6.BASIS,
             'cells': cells, 'q_mod1': [str(values[i] % 1) for i in range(len(cells))],
             'raw': {str(i): str(values[i]) for i in range(len(cells))},
             'pairCoefficients': [str(c) for c in gamma.FIXED_COEFFICIENTS],
             'sourceHash': provenance, 'compiler': COMPILER}
    if cross is not None:
        table['crossCheck'] = cross
    unary_gamma6.check_table(table)
    return table


def write_table(table, path):
    """One key per line; the cells one per line."""
    lines = []
    for key, value in table.items():
        if key == 'cells':
            text = '[\n' + ',\n'.join(json.dumps(c, separators=(',', ':')) for c in value) + '\n]'
        else:
            text = json.dumps(value)
        lines.append(json.dumps(key) + ': ' + text)
    Path(path).write_text('{\n' + ',\n'.join(lines) + '\n}\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--cells', nargs='*', type=int, metavar='INDEX',
                        help='cells to compute (default: all 22, none with --collect)')
    parser.add_argument('--jobs', type=int, default=1, help='parallel processes, one cell each')
    parser.add_argument('--record', metavar='FILE', help='append every computed cell to FILE')
    parser.add_argument('--collect', nargs='*', default=[], metavar='FILE',
                        help='read computed cells from record files')
    parser.add_argument('--reference', metavar='TABLE',
                        help='a table whose q_mod1 every entry must equal')
    parser.add_argument('--cache', metavar='DIR', help='universal-value store directory')
    parser.add_argument('--output', default=str(unary_gamma6.TABLE), metavar='FILE',
                        help='the table file to write')
    parser.add_argument('--no-write', action='store_true', help='compute and check only')
    args = parser.parse_args()
    indices = args.cells
    if indices is None:
        indices = [] if args.collect else list(range(unary_gamma6.CELLS))
    if any(not 0 <= i < unary_gamma6.CELLS for i in indices):
        parser.error(f'cell indices lie in 0..{unary_gamma6.CELLS-1}')
    records = read_records(args.collect)
    started = time.time()

    def report(record):
        line = json.dumps(record, separators=(',', ':'))
        print(line, flush=True)
        if args.record:
            with open(args.record, 'a') as out:
                out.write(line + '\n')
        records.append(record)

    if args.jobs > 1 and len(indices) > 1:
        context = multiprocessing.get_context('spawn')
        with ProcessPoolExecutor(max_workers=args.jobs, mp_context=context,
                                 max_tasks_per_child=1) as pool:
            for record in pool.map(compile_cell, indices, [args.cache]*len(indices)):
                report(record)
    else:
        for index in indices:
            report(compile_cell(index, args.cache))
    print(f'{len(indices)} cells computed in {time.time()-started:.0f} s', file=sys.stderr)
    if args.no_write:
        return
    table = assemble(records, args.reference)
    write_table(table, args.output)
    print(f'wrote {args.output}', file=sys.stderr)


if __name__ == '__main__':
    main()
