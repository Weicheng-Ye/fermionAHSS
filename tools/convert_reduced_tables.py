#!/usr/bin/env python3
"""Build the shipped tables of the reduced diagonal source from the shipped pair sources.

For every diagonal pair key Delta(tau) of production_gamma{k}.source_value in
data/universal-values/ the reduced source of unary_tables is, exactly,

    reduced(tau) = pair(Delta tau) + 2 V_n(key tau),

with V_n the one-fiber source primitive (low_phases.V2_pair for k = 5,
closed_a_upper._high_source_value for k = 6) at the re-sectioned simplex
(docs/power_relations.md, Section 12). The values whose leaf is shipped are
converted here without any evaluation; with --compute-missing the missing
leaves are evaluated (one V_n leaf each) and appended to the leaf table as
well. The manifest's table sizes are updated.

usage: python3 tools/convert_reduced_tables.py [--compute-missing] [--degree 5|6]
"""
import argparse
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault('FERMIONAHSS_CACHE_DIR', '')
sys.path.insert(0, str(ROOT / 'python'))
sys.path.insert(0, str(ROOT / 'python' / 'stacking_model'))
import universal_values as universal  # noqa: E402
import unary_tables  # noqa: E402
import chain_models as cm  # noqa: E402
import r2_pair_chain  # noqa: E402,F401
import r3_pair_chain  # noqa: E402,F401
import r3_chain  # noqa: E402

BUNDLE = universal.BUNDLED
LEAF = {5: ('low_phases.V2_pair', 'low_phases', 'V2_pair'),
        6: ('closed_a_upper._high_source_value', 'closed_a_upper', '_high_source_value')}


def rows(name, classes):
    path = BUNDLE / (name + '.jsonl')
    if not path.exists():
        return []
    out = []
    for line in path.read_text().splitlines():
        try:
            key, value = json.loads(line)
        except ValueError:
            continue
        out.append((universal._decode(key, classes), universal._decode(value, classes)))
    return out


def one_fiber(k, pair_key):
    """The raw one-fiber simplex tau of a diagonal pair key Delta(tau), or None."""
    x, w = pair_key
    if k == 5:
        if any(tuple(u.edges) != tuple(u.edges_prime) for u in x.rows):
            return None
        return (cm.Diag('zsign', tuple((u.sigma, tuple(u.edges)) for u in x.rows)), w)
    if any(u.matrix != u.matrix_prime for u in x.rows):
        return None
    return (r3_chain.Diag3(tuple(r3_chain.U2(u.sigma, u.matrix) for u in x.rows)), w)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--compute-missing', action='store_true')
    parser.add_argument('--degree', type=int, choices=(5, 6), nargs='*', default=[5, 6])
    args = parser.parse_args()
    classes = universal._key_classes()
    manifest_path = BUNDLE / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    if manifest.get('provenance') != universal.source_provenance():
        raise SystemExit('the bundle was made for other formula sources')
    for k in args.degree:
        leaf_name, module_name, function_name = LEAF[k]
        pair = rows(f'production_gamma{k}.source_value', classes)
        leaves = dict(rows(leaf_name, classes))
        reduced_name = f'unary_tables.reduced_source_value_{k}'
        known = {}
        for key, value in rows(reduced_name, classes):
            known[json.dumps(universal._encode(key, classes), sort_keys=True)] = value
        converted = missing = computed = 0
        started = time.time()
        reduced_lines, leaf_lines = [], []
        for pair_key, pair_value in pair:
            tau = one_fiber(k, pair_key)
            if tau is None:
                continue
            text = json.dumps(universal._encode(tau, classes), sort_keys=True)
            if text in known:
                continue
            normal = unary_tables._normal_form(k, tau)
            leaf = leaves.get(normal)
            if leaf is None:
                if not args.compute_missing:
                    missing += 1
                    continue
                module = __import__(module_name)
                leaf = getattr(module, function_name)(normal)
                leaves[normal] = leaf
                leaf_lines.append(json.dumps([universal._encode(normal, classes), universal._encode(leaf, classes)],
                                             separators=(',', ':')))
                computed += 1
                if computed % 10 == 0:
                    print(f'  degree {k}: {computed} leaves computed, {time.time() - started:.0f} s', flush=True)
            value = pair_value + 2 * leaf
            known[text] = value
            reduced_lines.append(json.dumps([universal._encode(tau, classes), universal._encode(value, classes)],
                                            separators=(',', ':')))
            converted += 1
        for name, lines in ((reduced_name, reduced_lines), (leaf_name, leaf_lines)):
            if lines:
                with open(BUNDLE / (name + '.jsonl'), 'a') as out:
                    out.write('\n'.join(lines) + '\n')
            if name == reduced_name:
                manifest['tables'][name] = len(known)
            elif lines:
                manifest['tables'][name] = len(leaves)
        print(f'degree {k}: {len(pair)} pair values, {converted} converted ({computed} leaves computed), '
              f'{missing} left without a shipped leaf, {time.time() - started:.0f} s', flush=True)
    manifest_path.write_text(json.dumps(manifest, indent=1) + '\n')


if __name__ == '__main__':
    main()
