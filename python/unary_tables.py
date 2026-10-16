"""One-fiber (unary) primitives of the diagonal pair phases in degrees five and six.

The relation of an A generator multiplies its flat lift by itself: the only
two-fiber universal primitive it meets is the pair primitive P_k(A, A) on the
diagonal (docs/power_relations.md). On the one-fiber universal model
(M_n; F, G, H) of V_n, n = k-3, with its closed-cocycle section j,

    P^(A) = S_Delta(H j(A)) + sum_b c_b <F j(A), b>,     c_b = (Delta^* P_k)(G b),

where S_Delta(A) = S_k(A, A) is the diagonal pair source and b runs over the
small cells of degree k. The contraction identity x = dHx + Hdx + GFx gives
P_k(A, A) = P^(A) + delta_s W + N with W rational and N integral, so the
relation changes only by a pure-D gauge. The constants c_b are compiled once
by tools/generate_unary_tables.py from the two-fiber primitive and shipped
in data/; a table is accepted only for the formula sources it was compiled
from. diagonal_phase(k, A, B, C) is production_gamma{k}.phase(x, x) with its
initial summand -P_k(A, A) replaced by -P^(A).

Degree four keeps the two-fiber primitive: its pair model (the shared-sign
degree-one pair) is small and never dominated a measurement.

The undoubled term. The diagonal source S_k(A_tau, A_tau) of every chain
term contains -V_P(A_tau) - V_P(A_tau) = -2 V_n(tau), the one-fiber source
primitive at the chain term itself. Summed over the chain H j(A) with the
chain's coefficients, these terms give exactly zero: H H = 0 and F H = 0
hold for the contraction, and the source phi_n vanishes on every simplex of
H applied to a re-sectioned chain term (verified term by term on random
closed cocycles at the production dimension, degrees five and six;
docs/power_relations.md). primitive() therefore evaluates the reduced
source S_k(A,A) + 2 V_n(A), which costs one V_n leaf per chain term instead
of two, through the persistent tables reduced_source_value_{5,6}; a stored
pair-source value converts exactly when the V_n value of the term is stored
as well. Only the H part of (P4) changes; the constants c_b and every result
are unchanged.

Copyright (c) 2026 koAHSS contributors; MIT license.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
ROOT = HERE.parent
from worker_api import api as _stacking_api  # noqa: E402,F401  (installs the formula import path)

SCHEMA = 1
DENOMINATOR = 48
TABLES = {
    5: dict(file='unary-degree5-coefficients.json', basis='h_tau_primitive.basis2(5)', cells=4),
    6: dict(file='unary-gamma6-coefficients.json', basis='r3_chain.total_basis(6)', cells=22),
}


class UnaryTableUnavailable(Exception):
    """The coefficient table is missing or does not match the current sources."""
    def __init__(self, reason):
        super().__init__(reason)
        self.reason = reason


@lru_cache(None)
def source_hash():
    """The hash of the formula sources, computed once per process."""
    import universal_values
    return universal_values.source_provenance()


@lru_cache(None)
def unary_cells(k):
    """The relative one-fiber small cells of degree k, in table order."""
    if k == 5:
        import h_tau_primitive as hp
        return tuple(hp.basis2(5))
    if k == 6:
        import r3_chain
        return tuple(r3_chain.total_basis(6))
    raise ValueError('unary tables exist for degrees five and six')


def encoded_cells(k):
    import universal_values
    classes = universal_values._key_classes()
    return [universal_values._encode(cell, classes) for cell in unary_cells(k)]


def pair_coefficients(k):
    if k == 5:
        import production_gamma5
        return tuple(production_gamma5.FIXED_COEFFICIENTS)
    import production_gamma6
    return tuple(production_gamma6.FIXED_COEFFICIENTS)


def _fraction(value, what):
    if not isinstance(value, str):
        raise UnaryTableUnavailable(f'{what} is not an exact rational string')
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError):
        raise UnaryTableUnavailable(f'{what} is not an exact rational: {value!r}') from None


def check_table(k, data, provenance=None):
    """The constants {cell: c_b} of a decoded table of degree k, or UnaryTableUnavailable."""
    spec = TABLES[k]
    if not isinstance(data, dict) or data.get('schema') != SCHEMA:
        raise UnaryTableUnavailable('the unary coefficient table has an unknown schema')
    if data.get('basis') != spec['basis']:
        raise UnaryTableUnavailable('the unary coefficient table names another basis')
    values = data.get('q_mod1')
    if not isinstance(values, list) or len(values) != spec['cells'] or any(v is None for v in values):
        raise UnaryTableUnavailable(f"the unary coefficient table needs {spec['cells']} non-null entries")
    q = [_fraction(v, f'entry {i}') for i, v in enumerate(values)]
    for i, value in enumerate(q):
        if not 0 <= value < 1:
            raise UnaryTableUnavailable(f'entry {i} is not reduced modulo one: {value}')
        if DENOMINATOR % value.denominator:
            raise UnaryTableUnavailable(
                f'entry {i} has a denominator not dividing {DENOMINATOR}: {value}')
    raw = data.get('raw')
    if raw is not None:
        if not isinstance(raw, dict) or set(raw) != {str(i) for i in range(spec['cells'])}:
            raise UnaryTableUnavailable('the exact values of the table do not cover its cells')
        for i, value in enumerate(q):
            if _fraction(raw[str(i)], f'exact value {i}') % 1 != value:
                raise UnaryTableUnavailable(f'entry {i} is not its exact value modulo one')
    if data.get('cells') != encoded_cells(k):
        raise UnaryTableUnavailable(
            'the unary coefficient table was compiled on other cells than ' + spec['basis'])
    pair = data.get('pairCoefficients')
    if (not isinstance(pair, list)
            or [_fraction(v, 'pair coefficient') for v in pair] != list(pair_coefficients(k))):
        raise UnaryTableUnavailable(
            'the unary coefficient table was compiled with other pair coefficients')
    expected = source_hash() if provenance is None else provenance
    if data.get('sourceHash') != expected:
        raise UnaryTableUnavailable(
            'the unary coefficient table was compiled from other formula sources; '
            'regenerate it with python3 tools/generate_unary_tables.py')
    return dict(zip(unary_cells(k), q))


@lru_cache(None)
def coefficients(k):
    """The checked constants c_b of degree k, keyed by the small cells."""
    path = ROOT / 'data' / TABLES[k]['file']
    try:
        data = json.loads(path.read_text())
    except OSError:
        raise UnaryTableUnavailable(f'the unary coefficient table {path} is missing') from None
    except ValueError:
        raise UnaryTableUnavailable(f'the unary coefficient table {path} is not JSON') from None
    return check_table(k, data)


def table_status(k):
    """Whether the table of degree k is valid, with the refusal reason when it is not."""
    try:
        coefficients(k)
    except UnaryTableUnavailable as refusal:
        return {'status': 'invalid', 'reason': refusal.reason, 'degree': k}
    return {'status': 'valid', 'degree': k}


# ---------------------------------------------------------------------------
# The diagonal source and the one-fiber primitive.

def diagonal_source(k, A, s, omega):
    """S_Delta(A) = S_k(A, A) on one universal simplex."""
    if k == 5:
        import production_gamma5
        return production_gamma5.source(A, A, s, omega)
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


def reduced_diagonal_source(k, A, s, omega):
    """S_k(A, A) + 2 V_P(A): the diagonal pair source without its undoubled V_n term."""
    if k == 5:
        import production_gamma5 as formula
        aux = formula.aux
        upper = aux.upper
        comparison = formula.comparison
        p = formula.p
        B = p.zero(3)
        alpha = upper.hp.hD(p.binary(A), p.binary(A), s)
        beta = upper.beta_sharp(A, B, A, B, s, omega)
        P = upper.primary(A, s, omega)
        kappa = upper.fsharp(A, B, s, omega)
        theta = aux.theta_phase(2).phase(*(aux.as_cochain(x) for x in (P, kappa, P, kappa, s, omega)))
        # upper.phi(2A, alpha, beta) keeps the doubled term -V_P(2A) and its prism.
        return (p.scale(upper.phi(A + A, alpha, beta, s, omega), -1) - p.Cochain(theta.degree, theta)
                + p.scale(comparison.A_phase(A, s, omega), 2) - comparison.A_phase(A + A, s, omega)
                + comparison.affine_B_source(A, A, s, omega) - formula.transgression_source(A, A, s, omega))
    if k != 6:
        raise ValueError('the reduced diagonal source covers degrees five and six')
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
    return upper.source_primitive(A + A, s, omega) + local


def _reduced_source_value(k, simplex):
    """The reduced diagonal source on a raw one-fiber universal simplex of degree k+1; zero on a zero fiber."""
    from fractions import Fraction as F
    d = diagonal_simplex(k, simplex)
    if k == 5:
        import v2_pair_shared as pair
        if not any(any(u.edges) for u in d[0].rows):
            return F(0)
    else:
        import v3_pair_shared as pair
        if not any(x for u in d[0].rows for row in u.matrix for x in row):
            return F(0)
    A, _, s = pair.from_diag(d[0])
    value = reduced_diagonal_source(k, A, s, pair.from_omega(d[1]))
    return F(value(tuple(range(k + 2))))


def reduced_source_value_5(simplex):
    return _reduced_source_value(5, simplex)


def reduced_source_value_6(simplex):
    return _reduced_source_value(6, simplex)


def _normal_form(k, simplex):
    """The key under which the one-fiber source primitive V_n stores its value at this simplex."""
    if k == 5:
        import phase_eval as pe
        return pe.to_diags(pe.from_diags(simplex))
    import r3_source as rs
    return rs.to_diags(rs.from_diags(simplex))


def _stored(table, key):
    """The stored value of a persistent table at key, or None (never evaluates)."""
    values = getattr(table, 'values', None)
    if isinstance(values, dict):
        return values.get(key)
    return None


def diagonal_simplex(k, simplex):
    """The two-fiber universal simplex Delta(tau) of a one-fiber simplex tau: both fibers carry its labels."""
    x, w = simplex
    if k == 5:
        import r2_pair_chain as rc
        return (rc.Diag2(tuple(rc.U1(row[0], tuple(row[1]), tuple(row[1])) for row in x.rows)), w)
    import r3_pair_chain as rc
    return (rc.Diag3(tuple(rc.U2(a.sigma, a.matrix, a.matrix) for a in x.rows)), w)


def _source_value(k, simplex):
    """S_Delta on a raw one-fiber universal simplex of degree k+1.

    It is the two-fiber source S_k(A, A) on the diagonal simplex Delta(tau),
    read through the persistent table of production_gamma{k}.source_value
    (zero on a zero fiber): the values shipped and cached for the pair
    primitive serve the one-fiber primitive, and conversely.
    """
    if k == 5:
        import production_gamma5 as formula
    else:
        import production_gamma6 as formula
    return formula.constructor().source_value(diagonal_simplex(k, simplex))


def source_value_5(simplex):
    return _source_value(5, simplex)


def source_value_6(simplex):
    return _source_value(6, simplex)


def _chains(k):
    if k == 5:
        import source_primitive as unary
    else:
        import r3_chain as unary
    return unary.Htot, unary.Ftot


def primitive(k, A, s, omega):
    """The one-fiber primitive P^(A) of a closed signed degree-(k-3) A, a cochain of degree k."""
    import h_tau_primitive as hp
    p = _stacking_api.p
    n = k - 3
    if (A.degree, s.degree, omega.degree) != (n, 1, 2):
        raise ValueError(f'the unary primitive of degree {k} expects degrees ({n},1,2)')
    small = coefficients(k)
    Htot, Ftot = _chains(k)
    reduced = {5: reduced_source_value_5, 6: reduced_source_value_6}[k]
    pair_table = {5: source_value_5, 6: source_value_6}[k]
    if k == 5:
        import production_gamma5 as formula
        import low_phases
        leaf_table = low_phases.V2_pair
    else:
        import production_gamma6 as formula
        import closed_a_upper
        leaf_table = closed_a_upper._high_source_value
    pair_store = formula.constructor().source_value

    def value_of(x):
        # The reduced value: stored (the shipped table covers the shipped
        # pair values), converted exactly from a stored pair value and the
        # V_n value of the term (one leaf at the undoubled label when it is
        # not stored, cheaper than the doubled leaf of the reduced source),
        # or evaluated.
        found = _stored(reduced, x)
        if found is not None:
            return found
        pair_value = _stored(pair_store, diagonal_simplex(k, x))
        if pair_value is not None:
            converted = pair_value + 2 * leaf_table(_normal_form(k, x))
            put = getattr(reduced, 'put', None)
            if put is not None:
                put(x, converted)
            return converted
        return reduced(x)

    def value(vertices):
        if not any(A(tuple(vertices[i] for i in f)) for f in combinations(range(len(vertices)), n + 1)):
            return Fraction(0)
        simplex = hp.to_pair(A, s, omega, vertices)
        return (sum(c * value_of(x) for x, c in Htot(simplex).items())
                + sum(c * small.get(b, Fraction(0)) for b, c in Ftot(simplex).items()))
    return p.Cochain(k, value)


def diagonal_phase(k, A, B, C, s, omega):
    """production_gamma{k}.phase(x, x) with -P_k(A,A) replaced by -P^(A)."""
    if k == 5:
        import production_gamma5 as formula
        aux = formula.aux
        upper = aux.upper
        comparison = formula.comparison
        p = formula.p
        I = upper.hp.interval
        alpha = upper.hp.hD(p.binary(A), p.binary(A), s)
        Bsum = p.binary(B + B + alpha)
        along = aux.reduced_full_source(I(A), I(B, True), I(C, True), I(A), I(B, True), I(C, True), I(s), I(omega))
        g1 = comparison.comparison_gauge(A, B, s, omega)
        gauge = g1 + g1 - comparison.comparison_gauge(A + A, Bsum, s, omega)
        return (p.scale(primitive(5, A, s, omega), -1) - upper.hp.prism(along)
                - comparison.affine_B_primitive(A, B, A, B, s, omega) - gauge
                - formula.transgression(A, B, A, B))
    if k != 6:
        raise ValueError('the unary diagonal phase covers degrees five and six')
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
    return (p.scale(primitive(6, A, s, omega), -1) - upper.hp.prism(along)
            - p.half(formula.binary.affine_J0_primitive(A, B, A, B, s, omega)) - gauge)
