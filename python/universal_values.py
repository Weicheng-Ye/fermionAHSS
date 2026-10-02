"""Universal values kept across worker processes.

A universal function depends only on its key, never on the group, the
resolution or the process. Its values are shipped in
data/universal-values.json and kept in the cache directory, in one
append-only file of JSON lines per table under a directory named by the hash
of every formula source; files with another hash are ignored. The koFull
worker keeps the tables listed in extension_acceleration.py and
universal_sources.py (the pair sources of degrees four to six, V1, V2, V3,
the theta values of the degree-three source and the theta-pair sources), the
page worker V1, V2 and the theta values.

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

ROOT = Path(__file__).resolve().parents[1]
# Frozen dataclasses that occur in universal keys; nothing else is decoded.
KEY_CLASSES = (('chain_models', 'Diag'), ('r1_pair_chain', 'Borel'),
               ('r1_pair_chain_signed', 'Borel'), ('r3_chain', 'Diag3'),
               ('r3_chain', 'U2'), ('r3_pair_chain', 'Diag3'), ('r3_pair_chain', 'U2'),
               ('r2_pair_chain', 'Diag2'), ('r2_pair_chain', 'U1'), ('a0_high_gamma', 'KB'))
# Runtime and test modules do not determine universal values.
_NOT_SOURCES = {'extension_transfer.py', 'extension_worker.py',
                'extension_degree_six.py', 'extension_native_upper.py',
                'extension_three_local.py',
                'extension_acceleration.py', 'universal_sources.py',
                'worker.py', 'runtime_cache.py', 'generate_universal_values.py',
                'universal_values.py', 'unary_gamma6.py', 'generate_unary_coefficients.py',
                'extension_light.py'}
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
        self._values = {}
        self._stored_entries = []
        self.pending = {}
        functools.update_wrapper(self, function)

    @property
    def values(self):
        # Decode a table only when it is used. Low-degree calculations do
        # not need the large dataclass keys of the higher universal phases.
        # Preserve insertion order: existing values, bundled values, then
        # cache entries, with the first value of a key winning as before.
        if not self._stored_entries:
            return self._values
        entries, self._stored_entries = self._stored_entries, []
        for rows, classes in entries:
            for key, value in rows:
                try:
                    self._values.setdefault(_decode(key, classes),
                                            _decode(value, classes))
                except (KeyError, TypeError, ValueError):
                    continue
        return self._values

    def __call__(self, key):
        try:
            return self.values[key]
        except KeyError:
            pass
        answer = self.function(key)
        self.values[key] = answer
        self.pending[key] = answer
        return answer

    def recompute(self, key):
        """Evaluate the function afresh on a stored key, bypassing the table."""
        return self.function(key)


class CompactTable(UniversalTable):
    """A table whose stored keys are exact, invertible string encodings.

    ``compact`` maps a key object to its string and ``expand`` inverts it;
    the key objects themselves are not retained. A string key is accepted
    directly, so persisted keys need no decoding.
    """
    def __init__(self, name, function, compact, expand):
        super().__init__(name, function)
        self.compact = compact
        self.expand = expand

    def __call__(self, key):
        short = key if isinstance(key, str) else self.compact(key)
        try:
            return self.values[short]
        except KeyError:
            pass
        answer = self.function(self.expand(short) if isinstance(key, str) else key)
        self.values[short] = answer
        self.pending[short] = answer
        return answer

    def recompute(self, key):
        return self.function(self.expand(key) if isinstance(key, str) else key)


class UniversalStore:
    """One append-only file of JSON lines per table, under the provenance directory."""
    def __init__(self, directory, provenance, tables):
        self.directory = (None if directory is None else
                          directory / f'universal-values-{provenance[:16]}')
        self.provenance = provenance
        self.tables = tables
        self.classes = _key_classes()

    def _table_path(self, name):
        return self.directory / (name + '.jsonl')

    def _read_lines(self, path):
        try:
            text = path.read_text()
        except OSError:
            return []
        entries = []
        for line in text.splitlines():
            try:
                key, value = json.loads(line)
            except ValueError:
                continue
            entries.append((key, value))
        return entries

    def _insert(self, tables):
        for name, entries in tables.items():
            table = self.tables.get(name)
            if table is None:
                continue
            table._stored_entries.append((entries, self.classes))

    def load(self):
        if self.directory is None:
            return
        for name in self.tables:
            path = self._table_path(name)
            if path.exists():
                self._insert({name: self._read_lines(path)})

    def load_bundled(self, path):
        """Add shipped values; they are never written to the cache file."""
        try:
            data = json.loads(path.read_text())
        except (OSError, ValueError):
            return
        if data.get('schema') == 1 and data.get('provenance') == self.provenance:
            self._insert(data.get('tables', {}))

    def flush(self):
        """Append new values to the tables' files; failures only skip persistence.

        Appending whole lines keeps the cost proportional to the new values;
        a key appended twice by two processes is read once, the first value
        winning, and a torn line is skipped when read.
        """
        if self.directory is None or not any(t.pending for t in self.tables.values()):
            return
        try:
            self.directory.mkdir(parents=True, exist_ok=True)
            for name, table in self.tables.items():
                if not table.pending:
                    continue
                lines = []
                for key, value in table.pending.items():
                    try:
                        encoded = _encode(key, self.classes)
                        answer = _encode(value, self.classes)
                    except _Unsupported:
                        continue
                    lines.append(json.dumps([encoded, answer], separators=(',', ':')))
                if lines:
                    with open(self._table_path(name), 'a') as out:
                        out.write('\n'.join(lines) + '\n')
        except OSError:
            self.directory = None
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
