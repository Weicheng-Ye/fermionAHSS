"""Recompute the universal values shipped in data/universal-values.json.

The koFull and page workers load this file before their cache when the
recorded provenance equals the hash of the current formula sources, so a
first run starts with these values. Every key listed in the file, and every key of the
cache files given with --add, is evaluated again with the current formulas;
the file is rewritten with the current provenance and one sorted entry per
line. Run it after any change to the formula sources, whose old values the
worker then ignores. Keys whose new value differs from a value found in the
inputs are reported.

usage: python3 python/generate_universal_values.py [--add CACHE.json ...]
"""
import argparse
import json
import os
from pathlib import Path
import sys
import time

# Evaluate every value afresh: neither a cache nor the old file is consulted.
os.environ['FERMIONAHSS_CACHE_DIR'] = ''
os.environ['FERMIONAHSS_BUNDLED_VALUES'] = '0'
sys.path.insert(0, str(Path(__file__).resolve().parent))
import extension_transfer  # noqa: E402,F401  (installs the worker's exact policy)
import extension_acceleration as acceleration  # noqa: E402
import universal_values as universal  # noqa: E402


def read(path, classes, keys, previous):
    data = json.loads(Path(path).read_text())
    for name, entries in data.get('tables', {}).items():
        for key, value in entries:
            text = json.dumps(key, sort_keys=True)
            keys.setdefault(name, {})[text] = universal._decode(key, classes)
            previous.setdefault((name, text), set()).add(json.dumps(value, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--add', nargs='*', default=[], metavar='CACHE.json',
                        help='cache files whose keys are added')
    args = parser.parse_args()
    acceleration.persist()
    classes = universal._key_classes()
    keys, previous = {}, {}
    for path in ([universal.BUNDLED] if universal.BUNDLED.exists() else []) + args.add:
        read(path, classes, keys, previous)
    started, changed, tables = time.time(), 0, {}
    # Source values first: they evaluate most of the V1 values they need.
    for name in sorted(keys, reverse=True):
        table = acceleration._tables[name]
        rows = tables[name] = []
        for text in sorted(keys[name]):
            value = json.dumps(universal._encode(table(keys[name][text]), classes),
                               sort_keys=True)
            changed += any(old != value for old in previous[(name, text)])
            rows.append('[' + text + ',' + value + ']')
    body = ',\n'.join(json.dumps(name) + ':[\n' + ',\n'.join(rows) + '\n]'
                      for name, rows in sorted(tables.items()))
    universal.BUNDLED.write_text(
        '{"schema":1,"provenance":' + json.dumps(universal.source_provenance())
        + ',"tables":{\n' + body + '\n}}\n')
    print(f'{sum(len(r) for r in tables.values())} values in {time.time() - started:.0f} s; '
          f'{changed} differ from the inputs; wrote {universal.BUNDLED}')


if __name__ == '__main__':
    main()
