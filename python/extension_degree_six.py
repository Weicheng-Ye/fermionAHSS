"""Bounded finite-section adapter for the selected degree-six formulas.

Only exact linear coordinates are supplied by GAP. The nonlinear operations
come from the bundled production model, whose files the worker verifies before
calling this module. Copyright (c) 2026 koAHSS contributors; MIT license.
"""
from functools import lru_cache


class DegreeSixResourceLimit(RuntimeError):
    """The finite section would exceed its declared search bound."""


def configure_degree_six(bar, payload):
    """Attach the production adapter, preserving the bar basis and twists.

    Matrices act on columns here. ``basis`` has the integral kernel columns
    first; ``deltaB`` and ``deltaC`` are lists of packed binary matrix rows.
    ``maxSectionCandidates`` bounds the worst-case number of (A,B) candidates
    for an uncached section request, rather than precomputing all states.
    """
    import four_cochain_stacking as api
    from g6_repair import IntegralCycleCoordinates
    from production_g6_section import ProductionG6Section

    required = {"delta", "basis", "inverse", "cycleRank", "dimB", "dimC",
                "deltaB", "deltaC"}
    if not isinstance(payload, dict) or not required <= payload.keys():
        raise ValueError("incomplete degree-six linear coordinates")
    if payload.keys() - required - {"maxSectionCandidates"}:
        raise ValueError("unknown degree-six coordinate option")
    dimensions = {q: bar.dimension(q) for q in range(3, 7)}
    for key, expected in (("dimB", dimensions[4]), ("dimC", dimensions[5])):
        if type(payload[key]) is not int or payload[key] != expected:
            raise ValueError("degree-six " + key + " does not match the bar basis")
    rank = payload["cycleRank"]
    if type(rank) is not int or not 0 <= rank <= dimensions[3]:
        raise ValueError("invalid degree-six integral cycle rank")
    limit = payload.get("maxSectionCandidates", 4096)
    if type(limit) is not int or limit < 1:
        raise ValueError("maxSectionCandidates must be a positive integer")

    def integer_matrix(key, rows, columns):
        value = payload[key]
        if (not isinstance(value, (list, tuple)) or len(value) != rows
                or any(not isinstance(row, (list, tuple)) or len(row) != columns
                       or any(type(x) is not int for x in row) for row in value)):
            raise ValueError("invalid degree-six " + key + " matrix")
        return tuple(tuple(row) for row in value)

    delta = integer_matrix("delta", dimensions[4], dimensions[3])
    basis = integer_matrix("basis", dimensions[3], dimensions[3])
    inverse = integer_matrix("inverse", dimensions[3], dimensions[3])

    def binary_matrix(key, rows, columns):
        value = payload[key]
        if (not isinstance(value, (list, tuple)) or len(value) != rows
                or any(type(row) is not int or row < 0 or row >> columns
                       for row in value)):
            raise ValueError("invalid packed degree-six " + key + " matrix")
        return tuple(value)

    delta_b = binary_matrix("deltaB", dimensions[5], dimensions[4])
    delta_c = binary_matrix("deltaC", dimensions[6], dimensions[5])

    # Check that the supplied linear data use this exact normalized basis.
    # In particular, a correct abstract kernel in another resolution is unsafe.
    for q, supplied, signed in ((3, delta, True), (4, delta_b, False),
                                (5, delta_c, False)):
        for j in range(dimensions[q]):
            unit = [int(i == j) for i in range(dimensions[q])]
            cochain = bar.cochain(q, unit)
            differential = api.p.ds(cochain, bar.s) if signed else api.p.differential(cochain)
            actual = bar.vector(differential)
            expected = ([row[j] for row in supplied] if signed else
                        [(row >> j) & 1 for row in supplied])
            if not signed:
                actual = [x % 2 for x in actual]
            if actual != expected:
                raise ValueError("degree-six differential differs from the bar basis")

    coordinates = IntegralCycleCoordinates(delta, basis, inverse, rank)

    def binary_cochain(q, value):
        dimension = bar.dimension(q)
        if type(value) is not int or value < 0 or value >> dimension:
            raise ValueError("invalid degree-six binary cochain coordinates")
        return bar.cochain(q, [(value >> j) & 1 for j in range(dimension)])

    def bits(cochain):
        return sum((value % 2) << j for j, value in enumerate(bar.vector(cochain)))

    section = ProductionG6Section(
        coordinates, dimensions[4], dimensions[5], delta_b, delta_c,
        bar.s, bar.w, bar.cochain, binary_cochain, bits,
        lambda cochain: tuple(bar.vector(cochain)))
    # Replace the model's unbounded memoization, rather than retaining a
    # second unbounded copy beneath this worker-owned cache.
    original_section = section.model.current_section.__wrapped__
    search_exponent = 2 * rank + len(section.model.closed_b)

    @lru_cache(maxsize=256)
    def bounded_section(image):
        # The original lexicographic search chooses exactly this pointed value.
        # Its defining equations were checked by the model constructor.
        if image == section.model.zero_image:
            return section.model.zero
        if search_exponent >= limit.bit_length():
            raise DegreeSixResourceLimit(
                "degree-six section requires a search bound of 2^"
                + str(search_exponent) + "; maxSectionCandidates=" + str(limit))
        return original_section(section.model, image)

    section.model.current_section = bounded_section
    section.max_section_candidates = limit
    section.section_search_exponent = search_exponent
    # Mutate the worker only after every linear identity has been checked.
    bar.rule = api.Stacking(bar.s, bar.w, bar.is_zero, degree_six=section)
    bar.degree_six = section
    return section
