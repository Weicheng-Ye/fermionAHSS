"""The degree-six unary coefficient table of the light A-over-D residue.

The unary replacement of the degree-six diagonal phase needs 22 fixed
constants q_b, one per relative unary small cell b of
r3_chain.total_basis(6):

    q_b = (Delta^* P_S)(G b)  modulo Z_(2),

with S the degree-seven pair source of production_gamma6, P_S = S H_P + c F_P
its fixed primitive on the pair contraction (F_P, G_P, H_P) of
r3_pair_chain with the coefficients c of production_gamma6.FIXED_COEFFICIENTS,
Delta the diagonal of the universal cocycle spaces and G the unary
contraction map of r3_chain. They are compiled once by
generate_unary_coefficients.py and shipped in
data/unary-gamma6-coefficients.json. The table is valid only for the formula
sources it was compiled from: coefficients() refuses it, raising
UnaryTableUnavailable, unless the file exists, has schema 1, 22 non-null
entries on the cells of r3_chain.total_basis(6), the pair coefficients of
production_gamma6.FIXED_COEFFICIENTS, denominators dividing 48, and the
source hash of the current formula sources. FERMIONAHSS_UNARY_TABLE names
another table file.

The runtime evaluates the unary primitive

    P^(A) = S_Delta(H_U j(A)) + q(F j(A))

(light-transport-proof (21)-(22)): S_Delta is the diagonal pull-back of the
pair source, expanded as -2 V3(A) + V3(2A) - Phi_0(2A,alpha,beta)
- Theta_pair + 2 R_A(A) - R_A(2A) + h(affine J0(A,A)), with each V3 value
taken on the cocycle section of its unary simplex; (F, G, H_U) is the unary
contraction of r3_chain and j its closed-cocycle section. diagonal_phase is
production_gamma6.phase on two equal triples with its initial summand
-P_S(A,A) replaced by -P^(A); it never evaluates the upper two-A pair
contractor.

Copyright (c) 2026 koAHSS contributors; MIT license.
"""
from fractions import Fraction
from functools import lru_cache
import json
import os
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
ROOT = HERE.parent
# The table validator is also a standalone API. Install the vendored formula
# import path even when no transfer worker has been imported by the caller.
from extension_worker import api as _stacking_api  # noqa: E402,F401

TABLE = ROOT / 'data' / 'unary-gamma6-coefficients.json'
SCHEMA = 1
CELLS = 22
# All 22 constants lie in (1/48)Z; an entry with another denominator is refused.
DENOMINATOR = 48
BASIS = 'r3_chain.total_basis(6)'


class UnaryTableUnavailable(Exception):
    """The coefficient table is missing or does not match the current sources."""
    def __init__(self, reason):
        super().__init__(reason)
        self.reason = reason


def table_path(path=None):
    """The table file: path, else FERMIONAHSS_UNARY_TABLE, else the bundled file."""
    if path is not None:
        return Path(path)
    configured = os.environ.get('FERMIONAHSS_UNARY_TABLE')
    if configured is not None and configured.strip():
        return Path(configured)
    return TABLE


@lru_cache(None)
def source_hash():
    """The hash of the formula sources, computed once per process."""
    import universal_values
    return universal_values.source_provenance()


@lru_cache(None)
def unary_cells():
    """The 22 relative unary small cells of degree six, in table order."""
    import r3_chain
    return tuple(r3_chain.total_basis(6))


def encoded_cells():
    import universal_values
    classes = universal_values._key_classes()
    return [universal_values._encode(cell, classes) for cell in unary_cells()]


def pair_coefficients():
    import production_gamma6
    return tuple(production_gamma6.FIXED_COEFFICIENTS)


def _fraction(value, what):
    if not isinstance(value, str):
        raise UnaryTableUnavailable(f'{what} is not an exact rational string')
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError):
        raise UnaryTableUnavailable(f'{what} is not an exact rational: {value!r}') from None


def check_table(data, provenance=None):
    """The constants {cell: q_b} of a decoded table, or UnaryTableUnavailable."""
    if not isinstance(data, dict) or data.get('schema') != SCHEMA:
        raise UnaryTableUnavailable('the unary coefficient table has an unknown schema')
    if data.get('basis') != BASIS:
        raise UnaryTableUnavailable('the unary coefficient table names another basis')
    values = data.get('q_mod1')
    if not isinstance(values, list) or len(values) != CELLS or any(v is None for v in values):
        raise UnaryTableUnavailable(f'the unary coefficient table needs {CELLS} non-null entries')
    q = [_fraction(v, f'entry {i}') for i, v in enumerate(values)]
    for i, value in enumerate(q):
        if not 0 <= value < 1:
            raise UnaryTableUnavailable(f'entry {i} is not reduced modulo one: {value}')
        if DENOMINATOR % value.denominator:
            raise UnaryTableUnavailable(
                f'entry {i} has a denominator not dividing {DENOMINATOR}: {value}')
    raw = data.get('raw')
    if raw is not None:
        if not isinstance(raw, dict) or set(raw) != {str(i) for i in range(CELLS)}:
            raise UnaryTableUnavailable('the exact values of the table do not cover its cells')
        for i, value in enumerate(q):
            if _fraction(raw[str(i)], f'exact value {i}') % 1 != value:
                raise UnaryTableUnavailable(f'entry {i} is not its exact value modulo one')
    if data.get('cells') != encoded_cells():
        raise UnaryTableUnavailable(
            'the unary coefficient table was compiled on other cells than ' + BASIS)
    pair = data.get('pairCoefficients')
    if (not isinstance(pair, list)
            or [_fraction(v, 'pair coefficient') for v in pair] != list(pair_coefficients())):
        raise UnaryTableUnavailable(
            'the unary coefficient table was compiled with other pair coefficients')
    expected = source_hash() if provenance is None else provenance
    if data.get('sourceHash') != expected:
        raise UnaryTableUnavailable(
            'the unary coefficient table was compiled from other formula sources; '
            'regenerate it with python3 python/generate_unary_coefficients.py')
    return dict(zip(unary_cells(), q))


@lru_cache(None)
def _load(path, mtime, size):
    try:
        data = json.loads(Path(path).read_text())
    except OSError:
        raise UnaryTableUnavailable(f'the unary coefficient table {path} is missing') from None
    except ValueError:
        raise UnaryTableUnavailable(f'the unary coefficient table {path} is not JSON') from None
    return check_table(data)


def coefficients(path=None):
    """The checked constants q_b, keyed by the cells of r3_chain.total_basis(6)."""
    path = table_path(path)
    try:
        status = path.stat()
    except OSError:
        raise UnaryTableUnavailable(f'the unary coefficient table {path} is missing') from None
    return dict(_load(str(path.resolve()), status.st_mtime_ns, status.st_size))


def table_status(path=None):
    """Whether the table is valid, with the refusal reason when it is not."""
    try:
        coefficients(path)
    except UnaryTableUnavailable as refusal:
        return {'status': 'invalid', 'reason': refusal.reason, 'path': str(table_path(path))}
    return {'status': 'valid', 'path': str(table_path(path))}


def diagonal_source(A, s, omega):
    """S_Delta(A) on one universal simplex: the pair source of production_gamma6 at (A, A)."""
    import production_gamma6 as formula
    aux = formula.aux
    upper = aux.upper
    p = formula.p
    B = p.zero(4)
    alpha = upper.hp.hD(p.binary(A), p.binary(A), s)
    beta = upper.beta_sharp(A, B, A, B, s, omega)
    P = upper.primary(A, s, omega)
    kappa = upper.fsharp(A, B, s, omega)
    theta = aux.theta_phase(3).phase(*(aux.as_cochain(x) for x in (P, kappa, P, kappa, s, omega)))
    local = (p.scale(aux.phi_without_A_source(A + A, alpha, beta, s, omega), -1)
             - p.Cochain(theta.degree, theta)
             + p.scale(formula.comparison.A_phase(A, s, omega), 2)
             - formula.comparison.A_phase(A + A, s, omega)
             + p.half(formula.binary.affine_J0_source(A, A, s, omega)))
    return (p.scale(upper.source_primitive(A, s, omega), -2)
            + upper.source_primitive(A + A, s, omega) + local)


@lru_cache(4096)
def source_value(simplex):
    """S_Delta on a raw unary universal simplex; zero on the zero A fiber."""
    import h_tau_primitive as hp
    if not any(x for a in simplex[0].rows for row in a.matrix for x in row):
        return Fraction(0)
    A, s, omega = hp.from_pair(simplex, 3)
    return Fraction(diagonal_source(A, s, omega)(tuple(range(8))))


def primitive(A, s, omega):
    """The unary degree-six primitive P^(A) of a closed signed degree-three A."""
    import h_tau_primitive as hp
    import r3_chain as unary
    import production_gamma6 as formula
    if (A.degree, s.degree, omega.degree) != (3, 1, 2):
        raise ValueError('the unary primitive expects degrees (3,1,2)')
    small = coefficients()

    def value(vertices):
        simplex = hp.to_pair(A, s, omega, vertices)
        return (sum(c * source_value(x) for x, c in unary.Htot(simplex).items())
                + sum(c * small.get(b, Fraction(0)) for b, c in unary.Ftot(simplex).items()))
    return formula.p.Cochain(6, value)


def diagonal_phase(A, B, C, s, omega):
    """production_gamma6.phase(x, x) with -P_S(A,A) replaced by -P^(A)."""
    import production_gamma6 as formula
    aux = formula.aux
    upper = aux.upper
    p = formula.p
    interval = upper.hp.interval
    alpha = upper.hp.hD(p.binary(A), p.binary(A), s)
    Bsum = p.binary(B + B + alpha)
    along = aux.reduced_full_source(interval(A), interval(B, True), interval(C, True),
                                    interval(A), interval(B, True), interval(C, True),
                                    interval(s), interval(omega))
    gauge = (formula.comparison.comparison_gauge(A, B, s, omega)
             + formula.comparison.comparison_gauge(A, B, s, omega)
             - formula.comparison.comparison_gauge(A + A, Bsum, s, omega))
    return (p.scale(primitive(A, s, omega), -1) - upper.hp.prism(along)
            - p.half(formula.binary.affine_J0_primitive(A, B, A, B, s, omega)) - gauge)
