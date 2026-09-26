"""Exact execution policy for the extension transfer worker.

The calibrated formula sources, and the modules whose hashes certify the
bundled calibration tables, stay byte-identical.  Like runtime_cache.py,
this module only changes how the worker evaluates those unchanged formulas;
every cochain value is exactly the value of the original evaluator.  It is
installed by extension_transfer.py alone, never by the page worker.

1. Evaluation: interval-cut words read their faces through precomputed
   C-level getters, calling a cochain enters its memo table without an
   intermediate Python frame, and memo tables are created without
   functools.update_wrapper.
2. Structural zeros: zero() returns one canonical zero per degree.  A word
   with a zero factor, a sum or difference with zero, a scalar multiple, mod
   two reduction, coboundary, quotient, interval lift, prism or transport of
   zero vanishes by its definition and returns that canonical zero.
3. Identity memoization: pure cochain builders return the same cochain for
   the same input objects, so a term shared by several products (for
   example the potential of one operand) is evaluated once.
4. Persistent universal values: values of universal source functions are
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
import operator
import os
from pathlib import Path
import sys
import tempfile

import cochain_tools
import phase_eval

ROOT = Path(__file__).resolve().parents[1]
_installed = False


def is_zero(cochain):
    """True only for a cochain known to vanish identically by construction."""
    return getattr(cochain, 'structural_zero', False)


# 1. Evaluation -------------------------------------------------------------

_Cochain = cochain_tools.Cochain
_lru_cache_wrapper = getattr(functools, '_lru_cache_wrapper', None)
_CacheInfo = getattr(functools, '_CacheInfo', None)


def _memo(evaluate):
    # The same unbounded C memo table that lru_cache(None) builds, without
    # copying metadata; a replaced cochain_tools.lru_cache policy still wins.
    if (cochain_tools.lru_cache is functools.lru_cache
            and _lru_cache_wrapper is not None and _CacheInfo is not None):
        wrapper = _lru_cache_wrapper(evaluate, None, False, _CacheInfo)
        wrapper.__wrapped__ = evaluate
        return wrapper
    return cochain_tools.lru_cache(None)(evaluate)


def _cochain_init(self, degree, evaluate, name=""):
    self.degree = degree
    self.evaluate = _memo(evaluate)
    self.name = name


def _getter(face):
    """Return a C-level function taking a simplex to its face tuple."""
    if len(face) == 1:
        return operator.itemgetter(slice(face[0], face[0] + 1))
    return operator.itemgetter(*face)


def _word_op(word, *cochains, integral_index=None):
    degree = sum(c.degree for c in cochains) - (len(word) - len(cochains))
    # Every interval-cut term multiplies one value of each factor.
    if any(is_zero(c) for c in cochains):
        return zero(degree)
    terms = cochain_tools.cut_terms(word, tuple(c.degree for c in cochains),
                                    integral_index)
    prepared = tuple((weight, tuple((c, _getter(face))
                                    for c, face in zip(cochains, faces)))
                     for faces, weight in terms)
    integral = integral_index is not None

    def evaluate(simplex):
        if type(simplex) is not tuple:
            simplex = tuple(simplex)
        total = 0
        for weight, factors in prepared:
            value = weight
            for c, get in factors:
                value *= c(get(simplex))
                if not value:
                    break
            total += value
        return total if integral else total % 2

    return _Cochain(degree, evaluate)


def _chi(c):
    faces, masks = phase_eval.chi_anf(c.degree)
    getters = tuple(_getter(face) for face in faces)

    @functools.lru_cache(None)
    def evaluate_active(active):
        return sum((active & mask) == mask for mask in masks) % 2

    def value(z):
        if type(z) is not tuple:
            z = tuple(z)
        return evaluate_active(sum(c(get(z)) << i for i, get in enumerate(getters)))
    return _Cochain(c.degree + 3, value)


# 2. Structural zeros --------------------------------------------------------

_ZEROS = {}
_original_zero = cochain_tools.zero


def zero(degree):
    answer = _ZEROS.get(degree)
    if answer is None:
        answer = _original_zero(degree)
        answer.structural_zero = True
        _ZEROS[degree] = answer
    return answer


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


# 3. Identity memoization -------------------------------------------------------

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
    'h_tau_primitive': ('hD', 'uniform_source', 'interval', 'prism'),
    'off_shell_beta': ('primary', 'curvature', 'alpha', 'raw_H', 'natural_f',
                       'coordinate_gauge', 'current_f'),
    'stacking_lower': ('alpha', 'legal_beta', 'all_cochain_beta'),
    'all_cochain_upper': ('r', 'K', 'current_curvature', 'coordinate', 'integer_gauge'),
    'natural_upper': ('curvature', 'integral', 'ghat'),
    'pure_c_normalization': ('phase_correction', 'correction'),
    'production_upper_binary_comparison': ('J0', 'affine_J0_source', 'affine_J0_primitive'),
    'upper_pair_source': ('source', 'full_source', 'reduced_full_source'),
    'phase_eval': ('source', 'theta', 'Q', 'E', 'QD', 'ds'),
}
# Enough for the terms shared by consecutive products; larger tables only
# retain memory (measured on the degree-four and degree-six workers).
MEMO_LIMIT = 32


# 4. Persistent universal values -------------------------------------------------

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
    # 1. Evaluation.
    _Cochain.__init__ = _cochain_init
    _Cochain.__call__ = property(operator.attrgetter('evaluate'))
    _rebind(cochain_tools.word_op, _word_op)
    _rebind(phase_eval.chi, _chi)
    # 2. Structural zeros.
    _Cochain.__add__ = _add
    _Cochain.__sub__ = _sub
    _Cochain.mod2 = _mod2
    _rebind(_original_zero, zero)
    _rebind(_original_differential, _differential)
    _rebind(phase_eval.scale, _scale_rule(phase_eval.scale))
    _rebind(phase_eval.divide, _zero_preserving(phase_eval.divide))
    hp = modules['h_tau_primitive']
    _rebind(hp.interval, _zero_preserving(hp.interval))
    _rebind(hp.prism, _zero_preserving(hp.prism, degree_shift=-1))
    natural = modules['natural_upper']
    _rebind(natural.integral, _zero_preserving(natural.integral))
    low = modules['low_phases']
    _rebind(low.transported, _transported_rule(low.transported))
    # 3. Identity memoization, wrapping the zero-aware versions above.
    for module_name, names in MEMOIZED.items():
        module = modules[module_name]
        for name in names:
            _rebind(getattr(module, name), _IdentityMemo(getattr(module, name), MEMO_LIMIT))
    # 4. Universal values; persist() attaches the cross-process store.
    for module_name, name in UNIVERSAL:
        original = getattr(modules[module_name], name)
        table = _tables[module_name + '.' + name] = UniversalTable(
            module_name + '.' + name, original)
        _rebind(original, table)
    _installed = True
