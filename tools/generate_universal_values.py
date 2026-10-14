#!/usr/bin/env python3
"""Recompute the universal values shipped in data/universal-values/.

The workers load this bundle (one JSON-lines file per table, read when the
table is first used) before their cache when its recorded provenance equals
the hash of the current formula sources. Every key of the shipped bundle,
and every key of the cache directories given with --add, is evaluated again
with the current formulas; the bundle is rewritten with the current
provenance. Run it after any change to the formula sources, whose old
values the workers then ignore. Keys whose new value differs from a value
found in the inputs are reported.

usage: python3 tools/generate_universal_values.py [--add CACHE ...]

A CACHE is a store directory of the cache directory (one JSON-lines file per
table), a bundle directory, or a legacy single-file bundle.
"""
import argparse
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
# Evaluate every value afresh: neither a cache nor the shipped bundle is consulted.
os.environ['FERMIONAHSS_CACHE_DIR'] = ''
sys.path.insert(0, str(ROOT / 'python'))
sys.path.insert(0, str(ROOT / 'python' / 'stacking_model'))
import universal_values as universal  # noqa: E402
BUNDLE = universal.BUNDLED
universal.BUNDLED = ROOT / 'data' / 'no-bundle-while-generating'
import stacking_worker  # noqa: E402,F401  (installs the worker's exact policy)
import acceleration  # noqa: E402
import universal_sources  # noqa: E402


def read(path, classes, keys, previous):
    """Read a bundle directory, a store directory, or a legacy single-file bundle."""
    path = Path(path)
    if path.is_dir():
        tables = {}
        for file in sorted(path.glob('*.jsonl')):
            entries = tables[file.stem] = []
            for line in file.read_text().splitlines():
                try:
                    entries.append(tuple(json.loads(line)))
                except ValueError:
                    continue
    else:
        tables = json.loads(path.read_text()).get('tables', {})
    for name, entries in tables.items():
        for key, value in entries:
            text = json.dumps(key, sort_keys=True)
            keys.setdefault(name, {})[text] = universal._decode(key, classes)
            previous.setdefault((name, text), set()).add(json.dumps(value, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--add', nargs='*', default=[], metavar='CACHE',
                        help='store or bundle directories, or legacy bundle files, whose keys are added')
    args = parser.parse_args()
    acceleration.persist()
    classes = universal._key_classes()
    keys, previous = {}, {}
    for path in ([BUNDLE] if BUNDLE.exists() else []) + args.add:
        read(path, classes, keys, previous)
    started, changed, sizes = time.time(), 0, {}
    BUNDLE.mkdir(parents=True, exist_ok=True)
    # Every registered table except the per-term theta table is written,
    # empty when no input lists its keys. Source values first: they evaluate
    # most of the V1 values they need.
    for name in sorted(acceleration._tables, reverse=True):
        if name in universal_sources.UNBUNDLED:
            continue
        table = acceleration._tables[name]
        rows = []
        for text in sorted(keys.get(name, {})):
            value = json.dumps(universal._encode(table.recompute(keys[name][text]), classes),
                               sort_keys=True)
            changed += any(old != value for old in previous.get((name, text), ()))
            rows.append('[' + text + ',' + value + ']')
        (BUNDLE / (name + '.jsonl')).write_text('\n'.join(rows) + ('\n' if rows else ''))
        sizes[name] = len(rows)
    (BUNDLE / 'manifest.json').write_text(json.dumps(
        {'schema': 2, 'provenance': universal.source_provenance(), 'tables': sizes}, indent=1) + '\n')
    print(f'{sum(sizes.values())} values in {time.time() - started:.0f} s; '
          f'{changed} differ from the inputs; wrote {BUNDLE}')


if __name__ == '__main__':
    main()
