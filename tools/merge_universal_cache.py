#!/usr/bin/env python3
"""Merge the universal values of a cache store into the shipped bundle.

usage: python3 tools/merge_universal_cache.py <store directory>...
A store directory is cache/universal-values-<hash>/ (one JSON-lines file per
table) written by the workers; only stores whose hash is the current
provenance are merged, and a key already in the bundle keeps its value. The
per-term theta table (universal_sources.UNBUNDLED) is never bundled.
"""
import json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'python')); sys.path.insert(0, str(ROOT / 'python' / 'stacking_model'))
import universal_values as universal  # noqa: E402
import universal_sources  # noqa: E402
provenance = universal.source_provenance()
bundle = universal.BUNDLED
manifest_path = bundle / 'manifest.json'
manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {'schema': 2, 'provenance': provenance, 'tables': {}}
if manifest.get('provenance') != provenance:
    raise SystemExit('the bundle was made for other formula sources; regenerate it with tools/generate_universal_values.py')
added = {}
for store in sys.argv[1:]:
    store = Path(store)
    if not store.name.endswith(provenance[:16]):
        print('skipping', store, '(other provenance)'); continue
    for file in sorted(store.glob('*.jsonl')):
        name = file.stem
        if name in universal_sources.UNBUNDLED:
            continue
        target = bundle / file.name
        known = set()
        if target.exists():
            for line in target.read_text().splitlines():
                try: known.add(json.dumps(json.loads(line)[0], sort_keys=True))
                except (ValueError, IndexError): continue
        new_lines = []
        for line in file.read_text().splitlines():
            try: key, value = json.loads(line)
            except (ValueError, TypeError): continue
            text = json.dumps(key, sort_keys=True)
            if text in known: continue
            known.add(text); new_lines.append(json.dumps([key, value], separators=(',', ':')))
        if new_lines:
            with open(target, 'a') as out:
                out.write('\n'.join(new_lines) + '\n')
        added[name] = added.get(name, 0) + len(new_lines)
        manifest['tables'][name] = len(known)
manifest_path.write_text(json.dumps(manifest, indent=1) + '\n')
print('added', added, 'bundle sizes', manifest['tables'])
