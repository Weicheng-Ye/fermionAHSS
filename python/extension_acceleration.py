"""Exact execution policy for the extension transfer worker.

The interval-cut engine itself lives in cochain_tools.py and phase_eval.py,
shared with the page worker. This module, installed by extension_transfer.py
alone, adds policies that change how the koFull worker evaluates the
unchanged formulas; every cochain value is exactly that of the original
evaluator.

1. Structural zeros: zero() returns one canonical zero per degree.  A word
   with a zero factor, a sum or difference with zero, a scalar multiple, mod
   two reduction, coboundary, quotient, interval lift, prism or transport of
   zero vanishes by its definition and returns that canonical zero.
2. Identity memoization: pure cochain builders return the same cochain for
   the same input objects, so a term shared by several products (for
   example the potential of one operand) is evaluated once.
3. Persistent universal values: values of universal source functions are
   kept across processes, keyed by the hashes of every formula source.

Copyright (c) 2026 koAHSS contributors; MIT license.
"""
from collections import OrderedDict
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

import cochain_tools
import phase_eval

ROOT = Path(__file__).resolve().parents[1]
_installed = False
_Cochain = cochain_tools.Cochain


def is_zero(cochain):
    """True only for a cochain known to vanish identically by construction."""
    return getattr(cochain, 'structural_zero', False)


# 1. Structural zeros --------------------------------------------------------

_ZEROS = {}
_original_zero = cochain_tools.zero


def zero(degree):
    answer = _ZEROS.get(degree)
    if answer is None:
        answer = _original_zero(degree)
        answer.structural_zero = True
        _ZEROS[degree] = answer
    return answer


def _word_rule(function):
    @functools.wraps(function)
    def word_op(word, *cochains, integral_index=None):
        # Every interval-cut term multiplies one value of each factor.
        if any(is_zero(c) for c in cochains):
            return zero(cochain_tools.word_degree(word, cochains))
        return function(word, *cochains, integral_index=integral_index)
    return word_op


def _check_degrees(a, b):
    if a.degree != b.degree:
        raise AssertionError('cochain degrees differ')


_original_add = _Cochain.__add__
_original_sub = _Cochain.__sub__
_original_mod2 = _Cochain.mod2


def _add(self, other):
    _check_degrees(self, other)
    if is_zero(other):
        return self
    if is_zero(self):
        return other
    return _original_add(self, other)


def _sub(self, other):
    _check_degrees(self, other)
    if is_zero(other):
        return self
    return _original_sub(self, other)


def _mod2(self):
    if is_zero(self):
        return self
    return _original_mod2(self)


_original_differential = cochain_tools.differential


def _differential(a):
    if is_zero(a):
        return zero(a.degree + 1)
    return _original_differential(a)


def _zero_preserving(function, degree_shift=0):
    """Return zero when the (linear) cochain argument vanishes."""
    @functools.wraps(function)
    def wrapper(c, *args, **kwargs):
        if is_zero(c):
            return zero(c.degree + degree_shift)
        return function(c, *args, **kwargs)
    return wrapper


def _scale_rule(function):
    @functools.wraps(function)
    def scale(c, q):
        if is_zero(c) or q == 0:
            return zero(c.degree)
        return function(c, q)
    return scale


def _transported_rule(function):
    @functools.wraps(function)
    def transported(left, right, s):
        if is_zero(left) or is_zero(right):
            return zero(left.degree + right.degree)
        return function(left, right, s)
    return transported


# 2. Identity memoization -------------------------------------------------------

class _IdentityMemo:
    """Bounded memo of a pure builder, keyed by its argument objects.

    Cochains compare by identity, so a hit requires the very same input
    cochains; keys hold strong references and cannot be recycled.
    """
    def __init__(self, function, limit):
        self.function = function
        self.limit = limit
        self.table = OrderedDict()
        functools.update_wrapper(self, function)

    def __call__(self, *args, **kwargs):
        key = (args, tuple(sorted(kwargs.items()))) if kwargs else args
        try:
            answer = self.table[key]
        except KeyError:
            answer = self.function(*args, **kwargs)
            self.table[key] = answer
            if len(self.table) > self.limit:
                self.table.popitem(last=False)
        except TypeError:
            return self.function(*args, **kwargs)
        else:
            self.table.move_to_end(key)
        if isinstance(answer, dict):
            return dict(answer)
        return answer


# Pure builders of the k=3..5 upper formulas and their lower-layer inputs.
MEMOIZED = {
    'closed_ab_upper': ('phase', 'curvature', 'potential', 'J'),
    'closed_a_upper': ('fsharp', 'beta_sharp', 'hE', 'primary', 'curvature'),
    'production_gamma4': ('source', 'primitive', 'lower_product', 'phase4'),
    'production_gamma4_comparison': (
        'source', 'shifted_lambda', 'fshift', 'curvature_gauge', 'source_prism',
        'reduced_phi_shift', 'closed_difference', 'shifted_A_phase',
        'reduced_A_phase', 'A_phase', 'shifted_comparison_gauge', 'comparison_gauge'),
    'upper_phase_diagnostic': ('production_phase',),
    'h_tau_primitive': ('uniform_source',),
    'off_shell_beta': ('primary', 'curvature', 'alpha', 'raw_H', 'natural_f',
                       'coordinate_gauge', 'current_f'),
    'stacking_lower': ('alpha', 'legal_beta', 'all_cochain_beta'),
    'all_cochain_upper': ('r', 'K', 'current_curvature', 'coordinate', 'integer_gauge'),
    'natural_upper': ('curvature', 'ghat'),
    'pure_c_normalization': ('phase_correction', 'correction'),
    'production_upper_binary_comparison': ('J0', 'affine_J0_source', 'affine_J0_primitive'),
    'upper_pair_source': ('source', 'full_source', 'reduced_full_source'),
    'phase_eval': ('source', 'source_splitting', 'theta', 'E', 'QD', 'ds', 'hD', 'integral'),
    'cochain_tools': ('Q', 'interval_pullback', 'right_prism'),
}
# Enough for the terms shared by consecutive products; larger tables only
# retain memory (measured on the degree-four and degree-six workers).
MEMO_LIMIT = 32


# 3. Persistent universal values -------------------------------------------------

UNIVERSAL = (('production_gamma4', 'source_value'), ('low_phases', 'V1_pair'))
# Frozen dataclasses that occur in universal keys; nothing else is decoded.
KEY_CLASSES = (('chain_models', 'Diag'), ('r1_pair_chain', 'Borel'),
               ('r1_pair_chain_signed', 'Borel'))
# Runtime and test modules do not determine universal values.
_NOT_SOURCES = {'extension_transfer.py', 'extension_worker.py',
                'extension_degree_six.py', 'extension_acceleration.py',
                'worker.py', 'runtime_cache.py'}


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

    def load(self):
        if self.path is None:
            return
        for name, entries in self._read().items():
            table = self.tables.get(name)
            if table is None:
                continue
            for key, value in entries:
                try:
                    table.values.setdefault(_decode(key, self.classes),
                                            _decode(value, self.classes))
                except (KeyError, TypeError, ValueError):
                    continue

    def flush(self):
        """Merge new values into the cache file; failures only skip persistence."""
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


_store = None


def flush():
    if _store is not None:
        _store.flush()


# Installation --------------------------------------------------------------------

def _rebind(original, replacement):
    """Replace every module-level reference, including imported names.

    This module keeps its own references to the originals it wraps.
    """
    this = sys.modules[__name__]
    for module in list(sys.modules.values()):
        namespace = getattr(module, '__dict__', None)
        if not namespace or module is this:
            continue
        for name, value in list(namespace.items()):
            if value is original:
                namespace[name] = replacement


def persist():
    """Keep universal values in the cache directory across worker processes."""
    global _store
    if not _installed or _store is not None:
        return
    directory = cache_directory()
    _store = UniversalStore(directory, source_provenance() if directory else '',
                            _tables)
    _store.load()


_tables = {}


def install():
    """Install the exact evaluation policy once, before any request."""
    global _installed
    if _installed:
        return
    modules = {name: importlib.import_module(name) for name in
               set(MEMOIZED) | {name for name, _ in UNIVERSAL}}
    # 1. Structural zeros.
    _Cochain.__add__ = _add
    _Cochain.__sub__ = _sub
    _Cochain.mod2 = _mod2
    _rebind(_original_zero, zero)
    _rebind(_original_differential, _differential)
    _rebind(cochain_tools.word_op, _word_rule(cochain_tools.word_op))
    _rebind(phase_eval.scale, _scale_rule(phase_eval.scale))
    _rebind(phase_eval.divide, _zero_preserving(phase_eval.divide))
    _rebind(cochain_tools.interval_pullback, _zero_preserving(cochain_tools.interval_pullback))
    _rebind(cochain_tools.right_prism,
            _zero_preserving(cochain_tools.right_prism, degree_shift=-1))
    _rebind(phase_eval.integral, _zero_preserving(phase_eval.integral))
    low = modules['low_phases']
    _rebind(low.transported, _transported_rule(low.transported))
    # 2. Identity memoization, wrapping the zero-aware versions above.
    for module_name, names in MEMOIZED.items():
        module = modules[module_name]
        for name in names:
            _rebind(getattr(module, name), _IdentityMemo(getattr(module, name), MEMO_LIMIT))
    # 3. Universal values; persist() attaches the cross-process store.
    for module_name, name in UNIVERSAL:
        original = getattr(modules[module_name], name)
        table = _tables[module_name + '.' + name] = UniversalTable(
            module_name + '.' + name, original)
        _rebind(original, table)
    _installed = True
