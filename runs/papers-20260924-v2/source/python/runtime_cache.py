"""Bound disposable cochain and universal-chain memo tables in the GAP worker.

This execution policy is separate from the mathematical source modules whose
hashes certify the fixed universal tables.  It changes only functools cache
capacity: evicted values are recomputed by the same exact evaluator.  It does
not discard registered source simplices, recycle source IDs, or alter any
chain operator, lift, coefficient, or normalization.
"""
from functools import lru_cache
from importlib import import_module
import gc

try:
    from . import cochain_tools
except ImportError:
    import cochain_tools


_configured_entries = None


def cochain_cache_limit():
    """Return the active policy, including diagnostic-process reconfiguration."""
    return _configured_entries


def install_cochain_cache_limit(entries=256):
    """Install the single-threaded worker policy before evaluating a request.

    Cochain constructors resolve ``lru_cache`` in cochain_tools at call time.
    Replace that module-local decorator, leaving functools itself unchanged.
    Existing Cochains (when the worker is imported into a diagnostic process)
    must also be bounded.  No cached result is part of the mathematical state.
    """
    global _configured_entries
    if type(entries) is not int or entries < 0:
        raise ValueError('cochain cache entries must be a nonnegative integer')
    if _configured_entries == entries:
        return entries

    def bounded_cache(maxsize=128, typed=False):
        return lru_cache(entries if maxsize is None else maxsize, typed=typed)

    replacements = [(obj, lru_cache(entries)(obj.evaluate.__wrapped__))
                    for obj in gc.get_objects()
                    if isinstance(obj, cochain_tools.Cochain)]
    # All fallible decorator calls finish before the active policy changes.
    cochain_tools.lru_cache = bounded_cache
    for obj, evaluate in replacements:
        obj.evaluate = evaluate
    _configured_entries = entries
    return entries


# These three modules contain pure integral chain operators.  Registration
# factories, source coefficient loaders, and the source-ID registry are
# deliberately outside this list.
_CHAIN_CACHE_MODULES = ('chain_models', 'r3_chain', 'exterior_bar')
_configured_chain_entries = None


def _chain_modules():
    return [import_module('.' + name, __package__) if __package__ else
            import_module(name) for name in _CHAIN_CACHE_MODULES]


def _chain_functions(module):
    for name, value in tuple(vars(module).items()):
        if (getattr(value, '__module__', None) == module.__name__
                and hasattr(value, 'cache_info')
                and hasattr(value, '__wrapped__')):
            yield name, value


def chain_cache_limit():
    """Return the active universal-chain memo capacity, or None if unconfigured."""
    return _configured_chain_entries


def install_chain_cache_limit(entries=256):
    """Bound pure chain-operator caches without altering mathematical sources.

    Call before a worker request.  Each module-owned functools cache is
    replaced by a bounded wrapper around precisely its original evaluator.
    Recursive calls resolve the newly bounded module globals.  Eviction only
    loses recomputable chain dictionaries: registered sources, source IDs,
    factory caches, and calibrated coefficient data are never changed.

    No old wrappers are retained here, since that would retain their former
    unbounded caches.  Reconfiguration is intended for a single-threaded
    worker or diagnostic process, not concurrent evaluators.
    """
    global _configured_chain_entries
    if type(entries) is not int or entries < 0:
        raise ValueError('chain cache entries must be a nonnegative integer')
    if _configured_chain_entries == entries:
        return entries
    replacements = []
    for module in _chain_modules():
        for name, value in _chain_functions(module):
            typed = value.cache_parameters()['typed']
            wrapped = lru_cache(entries, typed=typed)(value.__wrapped__)
            replacements.append((module, name, wrapped))
    # A failed import/decorator must leave every binding and policy unchanged.
    for module, name, wrapped in replacements:
        setattr(module, name, wrapped)
    _configured_chain_entries = entries
    return entries


def chain_cache_statistics():
    """Report the actual bounded wrappers for diagnostics and regression checks."""
    return {module.__name__.rsplit('.', 1)[-1] + '.' + name:
            value.cache_info()._asdict()
            for module in _chain_modules()
            for name, value in _chain_functions(module)}
