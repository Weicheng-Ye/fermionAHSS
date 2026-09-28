"""Universal values kept across worker processes.

A universal function depends only on its key, never on the group, the
resolution or the process. Its values are shipped in
data/universal-values.json and kept in a cache file of the cache directory,
both keyed by the hash of every formula source; a file with another hash is
ignored. The koFull worker keeps the tables listed in
extension_acceleration.py (the degree-four and degree-six pair sources, V1,
V3 and the n=3 legal-beta source), the page worker low_phases.V1_pair.

Copyright (c) 2026 koAHSS contributors; MIT license.
"""
from fractions import Fraction
import dataclasses
import functools
import hashlib
import importlib
import json
import os
from pathlib import Path
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
# Frozen dataclasses that occur in universal keys; nothing else is decoded.
KEY_CLASSES = (('chain_models', 'Diag'), ('r1_pair_chain', 'Borel'),
               ('r1_pair_chain_signed', 'Borel'), ('r3_chain', 'Diag3'),
               ('r3_chain', 'U2'), ('r3_pair_chain', 'Diag3'), ('r3_pair_chain', 'U2'),
               ('a0_high_gamma', 'KB'))
# Runtime and test modules do not determine universal values.
_NOT_SOURCES = {'extension_transfer.py', 'extension_worker.py',
                'extension_degree_six.py', 'extension_native_six.py',
                'extension_acceleration.py',
                'worker.py', 'runtime_cache.py', 'generate_universal_values.py',
                'universal_values.py'}
# Values computed with the current sources and shipped with the package; see
# generate_universal_values.py. A file with other provenance is ignored.
BUNDLED = ROOT / 'data' / 'universal-values.json'


def bundled_values_enabled():
    return os.environ.get('FERMIONAHSS_BUNDLED_VALUES', '1').strip() != '0'


def cache_directory():
    configured = os.environ.get('FERMIONAHSS_CACHE_DIR')
    if configured is not None:
        return Path(configured) if configured.strip() else None
    base = os.environ.get('XDG_CACHE_HOME') or str(Path.home() / '.cache')
    return Path(base) / 'fermionAHSS'


def source_provenance():
    """Hash of every formula source, calibration table and packaged ANF."""
    digest = hashlib.sha256()
    folder = ROOT / 'python'
    paths = [path for pattern in ('*.py', '*.json')
             for path in folder.glob(pattern)
             if path.name not in _NOT_SOURCES and not path.name.startswith('test_')]
    paths += list((folder / 'stacking_model').glob('*.py'))
    paths += [folder / 'stacking_model' / 'provenance.json',
              ROOT / 'data' / 'chi-calibrated-degree7-anf.g']
    for path in sorted(paths):
        digest.update(str(path.relative_to(ROOT)).encode() + b'\0')
        digest.update(hashlib.sha256(path.read_bytes()).digest())
    return digest.hexdigest()


def _key_classes():
    classes = {}
    for module_name, name in KEY_CLASSES:
        module = sys.modules.get(module_name)
        if module is None:
            try:
                module = importlib.import_module(module_name)
            except ImportError:
                continue
        classes[module_name + '.' + name] = getattr(module, name)
    return classes


class _Unsupported(Exception):
    pass


def _encode(value, classes):
    if value is None or isinstance(value, (bool, str)):
        return value
    if isinstance(value, int):
        return value
    if isinstance(value, Fraction):
        return {'f': f'{value.numerator}/{value.denominator}'}
    if isinstance(value, tuple):
        return {'t': [_encode(v, classes) for v in value]}
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        name = type(value).__module__ + '.' + type(value).__qualname__
        if classes.get(name) is not type(value):
            raise _Unsupported(name)
        return {'d': name, 'v': [_encode(getattr(value, field.name), classes)
                                 for field in dataclasses.fields(value)]}
    raise _Unsupported(type(value).__name__)


def _decode(data, classes):
    if not isinstance(data, dict):
        return data
    if 'f' in data:
        return Fraction(data['f'])
    if 't' in data:
        return tuple(_decode(v, classes) for v in data['t'])
    if 'd' in data:
        return classes[data['d']](*(_decode(v, classes) for v in data['v']))
    raise ValueError('malformed universal cache entry')


class UniversalTable:
    """Persistent values of one pure universal function of a hashable key."""
    def __init__(self, name, function):
        self.name = name
        self.function = function
        self.values = {}
        self.pending = {}
        functools.update_wrapper(self, function)

    def __call__(self, key):
        try:
            return self.values[key]
        except KeyError:
            pass
        answer = self.function(key)
        self.values[key] = answer
        self.pending[key] = answer
        return answer


class UniversalStore:
    def __init__(self, directory, provenance, tables):
        self.path = (None if directory is None else
                     directory / f'universal-values-{provenance[:16]}.json')
        self.provenance = provenance
        self.tables = tables
        self.classes = _key_classes()

    def _read(self):
        try:
            data = json.loads(self.path.read_text())
        except (OSError, ValueError):
            return {}
        if data.get('schema') != 1 or data.get('provenance') != self.provenance:
            return {}
        return data.get('tables', {})

    def _insert(self, tables):
        for name, entries in tables.items():
            table = self.tables.get(name)
            if table is None:
                continue
            for key, value in entries:
                try:
                    table.values.setdefault(_decode(key, self.classes),
                                            _decode(value, self.classes))
                except (KeyError, TypeError, ValueError):
                    continue

    def load(self):
        if self.path is not None:
            self._insert(self._read())

    def load_bundled(self, path):
        """Add shipped values; they are never written to the cache file."""
        try:
            data = json.loads(path.read_text())
        except (OSError, ValueError):
            return
        if data.get('schema') == 1 and data.get('provenance') == self.provenance:
            self._insert(data.get('tables', {}))

    def flush(self):
        """Merge new values into the cache file; failures only skip persistence.

        Tables of other workers found in the file are written back unchanged.
        """
        if self.path is None or not any(t.pending for t in self.tables.values()):
            return
        try:
            stored = self._read()
            for name, table in self.tables.items():
                entries = stored.setdefault(name, [])
                known = {json.dumps(key, sort_keys=True) for key, _ in entries}
                for key, value in table.pending.items():
                    try:
                        encoded = _encode(key, self.classes)
                        answer = _encode(value, self.classes)
                    except _Unsupported:
                        continue
                    text = json.dumps(encoded, sort_keys=True)
                    if text not in known:
                        known.add(text)
                        entries.append([encoded, answer])
            self.path.parent.mkdir(parents=True, exist_ok=True)
            handle, temporary = tempfile.mkstemp(dir=self.path.parent,
                                                 prefix='.universal-', suffix='.tmp')
            with os.fdopen(handle, 'w') as out:
                json.dump(dict(schema=1, provenance=self.provenance, tables=stored), out,
                          separators=(',', ':'))
            os.replace(temporary, self.path)
        except OSError:
            self.path = None
            return
        for table in self.tables.values():
            table.pending.clear()


def open_store(tables):
    """The store of these tables: bundled values first, then the cache file."""
    store = UniversalStore(cache_directory(), source_provenance(), tables)
    if bundled_values_enabled():
        store.load_bundled(BUNDLED)
    store.load()
    return store
