#!/usr/bin/env python3
"""Convert a single-file universal-value bundle into the per-table directory layout.

usage: python3 tools/split_universal_bundle.py data/universal-values.json data/universal-values
The directory holds manifest.json (schema 2, provenance, table sizes) and one
<table>.jsonl file per table; python/universal_values.py reads a table's file
only when the table is first used.
"""
import json, sys
from pathlib import Path
source, target = Path(sys.argv[1]), Path(sys.argv[2])
data = json.loads(source.read_text())
if data.get('schema') != 1:
    raise SystemExit('expected a schema-1 bundle')
target.mkdir(parents=True, exist_ok=True)
sizes = {}
for name, rows in data['tables'].items():
    with open(target / (name + '.jsonl'), 'w') as out:
        for key, value in rows:
            out.write(json.dumps([key, value], separators=(',', ':')) + '\n')
    sizes[name] = len(rows)
(target / 'manifest.json').write_text(json.dumps({'schema': 2, 'provenance': data['provenance'], 'tables': sizes}, indent=1) + '\n')
print('wrote', target, sizes)
