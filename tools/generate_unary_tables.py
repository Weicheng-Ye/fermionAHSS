#!/usr/bin/env python3
"""Compile the one-fiber coefficient tables data/unary-degree{5,6}-coefficients.json.

The constants of degree k (docs/power_relations.md, (P4)) are

    c_b = (Delta^* P_k)(G b)  modulo Z,     b in the small cells of degree k,

with G the one-fiber contraction map (source_primitive.Gtot for k=5, r3_chain.Gtot
for k=6), Delta the diagonal of the universal models and P_k the two-fiber
primitive with the fixed coefficients of production_gamma{k}. Every cell is
evaluated exactly with the two-fiber contraction; the result is written only
when every cell is known and has a denominator dividing 48.

usage: python3 tools/generate_unary_tables.py --degree 5 [--cells 0 1 ...]
"""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'python'))
sys.path.insert(0, str(ROOT / 'python' / 'stacking_model'))
import unary_tables  # noqa: E402

DENOMINATOR = 48


def cell_value(k, index):
    """The exact (Delta^* P_k)(G b) of the cell b of the given index."""
    cell = unary_tables.unary_cells(k)[index]
    if k == 5:
        import source_primitive as unary
        import r2_pair_chain as rc
        import production_gamma5 as gamma
        from universal_pair6 import one_fiber_zero
        prim = gamma.constructor()
        small = dict(zip(rc.total_basis(5), prim.coefficients()))

        def diagonal(x):
            return rc.Diag2(tuple(rc.U1(row[0], tuple(row[1]), tuple(row[1])) for row in x.rows))
        Htot, Ftot = rc.Htot, rc.Ftot
    elif k == 6:
        import r3_chain as unary
        import r3_pair_chain as rc
        import production_gamma6 as gamma
        from universal_pair7 import one_fiber_zero
        prim = gamma.constructor()
        small = dict(zip(rc.total_basis(6), gamma.FIXED_COEFFICIENTS))

        def diagonal(x):
            return rc.Diag3(tuple(rc.U2(a.sigma, a.matrix, a.matrix) for a in x.rows))
        Htot, Ftot = rc.Htot, rc.Ftot
    else:
        raise ValueError('unary tables exist for degrees five and six')
    total = Fraction(0)
    for (x, w), c in unary.Gtot(cell).items():
        d = (diagonal(x), w)
        if one_fiber_zero(d):
            continue
        total += c * (sum(c2 * prim.source_value(y) for y, c2 in Htot(d).items())
                      + sum(c2 * small.get(b, Fraction(0)) for b, c2 in Ftot(d).items()))
    return total


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--degree', type=int, choices=(5, 6), required=True)
    parser.add_argument('--cells', type=int, nargs='*')
    parser.add_argument('--output')
    args = parser.parse_args()
    k = args.degree
    cells = unary_tables.unary_cells(k)
    indices = args.cells if args.cells else list(range(len(cells)))
    values = {}
    for index in indices:
        started = time.monotonic()
        values[index] = cell_value(k, index)
        print(json.dumps({'degree': k, 'index': index, 'cell': repr(cells[index]), 'value': str(values[index]),
                          'mod1': str(values[index] % 1), 'seconds': round(time.monotonic() - started, 1)}), flush=True)
    if set(values) != set(range(len(cells))):
        print('partial run: the table is written only when every cell is known', file=sys.stderr)
        return
    for index, value in values.items():
        if DENOMINATOR % (value % 1).denominator:
            raise SystemExit(f'cell {index} has a denominator not dividing {DENOMINATOR}: {value}')
    table = {'schema': unary_tables.SCHEMA,
             'definition': f'c_b = (Delta^* P_{k})(G b) mod Z; docs/power_relations.md (P4)',
             'basis': unary_tables.TABLES[k]['basis'],
             'cells': unary_tables.encoded_cells(k),
             'q_mod1': [str(values[i] % 1) for i in range(len(cells))],
             'raw': {str(i): str(values[i]) for i in range(len(cells))},
             'pairCoefficients': [str(c) for c in unary_tables.pair_coefficients(k)],
             'sourceHash': unary_tables.source_hash(),
             'compiler': 'tools/generate_unary_tables.py'}
    output = Path(args.output) if args.output else ROOT / 'data' / unary_tables.TABLES[k]['file']
    output.write_text(json.dumps(table, indent=1) + '\n')
    print('wrote', output)


if __name__ == '__main__':
    main()
