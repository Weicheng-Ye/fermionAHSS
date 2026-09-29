"""Exact worker regression on the C2 comparison and its protocol contracts.

The C2 native and normalized-bar coordinates coincide; the separate GAP
transport tests exercise sparse nonidentity comparisons. These tests do not
assert completeness of the transferred gauge relation.
Copyright (c) 2026 koAHSS contributors; MIT license.
"""
import contextlib
import io
import json
import os
import unittest
from unittest.mock import patch

from extension_worker import FiniteBar
from extension_degree_six import configure_degree_six
import extension_transfer
from extension_transfer import TransferredModel, TransferResourceLimit, degrees, integer
import extension_native_upper as native_module


def c2_models(sign=0, omega=0, top=7):
    """The C2 models with resolution degrees 0..top; degree six needs top=8."""
    bar = FiniteBar([[0, 1], [1, 0]], [sign], [omega])
    setup = dict(schema=1, multiplication=bar.mul, ranks=[1] * (top + 1),
        ordinary=[[[0 if n % 2 == 0 else 2]] for n in range(top)],
        signed=[[[(-2 if n % 2 == 0 else 0) if sign else
                  (0 if n % 2 == 0 else 2)]] for n in range(top)],
        g=[[[[1, 1, list(bar.simplices(n)[0])]]] for n in range(top + 1)],
        s=[sign], omega=[omega])
    calls = []
    def transport(kind, degree, vertices):
        calls.append((kind, degree, vertices))
        return dict(status="computed", terms=[[0, 1, 1]] if kind == "f" else [])
    # Construction must validate twists natively and never construct Stacking.
    with patch("extension_worker.api.Stacking", side_effect=AssertionError("dense constructor")):
        native = TransferredModel(setup, transport)
    return native, bar, calls


def state(a, b, c, d):
    return dict(A=a, B=b, C=c, D=d)


def same_class(bar, k, actual, expected):
    """Equal A, B and C, and D layers differing by a signed coboundary on the bar.

    The native product omits the pure-C normalization of the production
    correction, an integral coboundary, so its D representatives differ from
    the complete-bar model's by D-gauges.
    """
    from extension_worker import api
    if any(actual[f] != expected[f] for f in "ABC"):
        return False
    differences = [a - e for a, e in zip(actual["D"], expected["D"])]
    if not any(differences):
        return True
    generators = []
    for j in range(bar.dimension(k)):
        unit = [int(i == j) for i in range(bar.dimension(k))]
        generators.append(bar.vector(api.p.ds(bar.cochain(k, unit), bar.s)))
    if len(differences) != 1:
        raise NotImplementedError("the coboundary test is written for one-cell bars")
    modulus = generators[0][0] if generators else 0
    return modulus != 0 and differences[0] % modulus == 0


def random_cochains(seed, dimension):
    """Random binary/integral cochains on the standard simplex of a dimension."""
    import random
    from itertools import combinations
    from extension_worker import api
    p = api.p
    rng = random.Random(seed)

    def faces(degree):
        return list(combinations(range(dimension + 1), degree + 1))

    def binary(degree):
        table = {f: rng.randrange(2) for f in faces(degree)}
        return p.Cochain(degree, lambda f: table[tuple(f)])

    def closed_binary(degree):
        return p.binary(p.differential(binary(degree - 1)))

    def integral(degree, low=-2, high=3):
        table = {f: rng.randrange(low, high) for f in faces(degree)}
        return p.Cochain(degree, lambda f: table[tuple(f)])

    def twists(sign, omega):
        return (closed_binary(1) if sign else p.zero(1), closed_binary(2) if omega else p.zero(2))
    return dict(binary=binary, closed_binary=closed_binary, integral=integral, twists=twists, p=p)


def lazy_c2_models(sign=0, omega=0, top=7):
    """The C2 models of c2_models, with the chains g requested from the transport."""
    bar = FiniteBar([[0, 1], [1, 0]], [sign], [omega])
    eager, _, _ = c2_models(sign, omega, top)
    setup = dict(schema=1, multiplication=bar.mul, ranks=[1] * (top + 1),
        ordinary=eager.matrices[False], signed=eager.matrices[True],
        gMode="lazy", s=[sign], omega=[omega])
    calls = []
    def transport(kind, degree, vertices):
        calls.append((kind, degree, vertices))
        if kind == "g":
            return dict(status="computed", terms=[[1, 1, list(bar.simplices(degree)[0])]])
        return dict(status="computed", terms=[[0, 1, 1]] if kind == "f" else [])
    return TransferredModel(setup, transport), eager, calls


class ExtensionTransferTests(unittest.TestCase):
    def test_twists_do_not_request_bar_values_during_setup(self):
        _, _, calls = c2_models(1, 1)
        self.assertEqual(calls, [])

    def test_lazy_chains_are_requested_once_on_first_use(self):
        native, eager, calls = lazy_c2_models(1, 1)
        self.assertEqual(calls, [])
        # The zero test stops at the first nonzero pairing; a projection
        # requests every chain of its degree once and keeps it.
        self.assertFalse(native.is_zero(native.lift(2, [1], False), True))
        self.assertEqual([c for c in calls if c[0] == "g"], [("g", 2, 0)])
        self.assertEqual(native.project(native.lift(3, [3], True), True), [3])
        self.assertEqual(native.project(native.lift(3, [5], True), True), [5])
        self.assertEqual([c for c in calls if c[0] == "g"], [("g", 2, 0), ("g", 3, 0)])
        # The same values as with the chains supplied in full.
        for k, data in ((3, state([0], [1], [0], [3])), (4, state([1], [0], [1], [2]))):
            self.assertEqual(native.kappa(k, data), eager.kappa(k, data))
        with self.assertRaises(ValueError):
            native.chain(2, 1)

    def test_projection_is_linear_and_uses_checked_identities(self):
        native, _, calls = lazy_c2_models(1, 1)
        p = native.lift(3, [2], True).__class__
        def direct(cochain, weight):
            return [sum(t[weight] * cochain(tuple(t[2])) for t in native.chain(cochain.degree, j))
                    for j in range(native.ranks[cochain.degree])]
        a, b = native.lift(3, [2], True), native.lift(3, [5], True)
        other = p(3, lambda t: 7)
        total = a + b - other
        self.assertEqual(native.project(total, True), [integer(v) for v in direct(total, 1)])
        self.assertEqual(native.project(total, True), [2 + 5 - 7])
        # Each lift pairs through the identity f g = 1, checked once per degree.
        self.assertTrue(native._identities[("lift", 3, 1, True)])
        # A pairing is kept per cochain: the same sum is not paired twice.
        count = len(calls)
        self.assertEqual(native.project(total, True), [0])
        self.assertEqual(len(calls), count)
        binary = (native.lift(2, [1], False) + native.lift(2, [1], True)).mod2()
        self.assertEqual(native.project(binary, False), [integer(v) % 2 for v in direct(binary, 0)])
        image = native.homotopy(native.lift(4, [3], True), True)
        self.assertEqual(native.project(image, True), [integer(v) for v in direct(image, 1)])

    def test_projection_evaluates_when_an_identity_fails(self):
        native, _, _ = lazy_c2_models()
        # A comparison with f g = 2 on degree three: the lift identity fails,
        # and the pairing is evaluated on the chains instead.
        transport = native.transport
        native.transport = lambda kind, degree, vertices: (
            dict(status="computed", terms=[[0, 2, 2]]) if kind == "f"
            else transport(kind, degree, vertices))
        native._transport.clear()
        self.assertEqual(native.project(native.lift(3, [3], True), True), [6])
        self.assertFalse(native._identities[("lift", 3, 1, True)])

    def test_lazy_chain_refusal_is_a_resource_limit(self):
        native, _, _ = lazy_c2_models()
        native.transport = lambda kind, degree, vertices: dict(
            status="unresolved", code="support-budget",
            reason="sparse comparison exceeds the per-degree support budget")
        with self.assertRaises(TransferResourceLimit):
            native.project(native.lift(2, [1], True), True)

    def test_curvature_prefix_stops_after_first_obstruction(self):
        native, bar, _ = c2_models(1, 0)
        data = state([1], [1], [1], [-2])
        complete = bar.export(bar.rule.d(bar.state(3, data)))
        self.assertEqual(complete, state([-2], [0], [1], [-3]))
        self.assertEqual(native.kappa(3, data), state([-2], [0], [0], [0]))
        data = state([0], [1], [0], [0])
        complete = bar.export(bar.rule.d(bar.state(3, data)))
        self.assertNotEqual(complete["C"], [0])
        self.assertNotEqual(complete["D"], [0])
        self.assertEqual(native.kappa(3, data), dict(complete, D=[0]))

    def test_full_flat_curvature_matches_complete_bar(self):
        native, bar, _ = c2_models()
        data = state([1], [0], [1], [7])
        self.assertEqual(native.kappa(3, data),
                         bar.export(bar.rule.d(bar.state(3, data))))

    def test_product_and_division_keep_integral_d_carries(self):
        native, bar, _ = c2_models()
        left = state([0], [1], [0], [7])
        right = state([-1], [1], [0], [-4])
        actual, gauge = native.product(3, left, right)
        expected = bar.export(bar.rule.xtimes(bar.state(3, left), bar.state(3, right)))
        self.assertEqual(actual, expected)
        self.assertEqual(native.values(gauge), state([], [0], [0], [0]))
        self.assertEqual(native.divide_left(3, left, actual), right)
        self.assertEqual(native.product(3, dict(left, D=[18]), dict(right, D=[-9]))[0],
                         dict(actual, D=[actual["D"][0] + 6]))

    def test_action_uses_the_full_bar_boundary_of_a_nonflat_gauge(self):
        native, bar, _ = c2_models(1, 1)
        gauge = state([], [1], [0], [3])
        canonical = state([0], [0], [1], [0])
        boundary = bar.rule.d(bar.state(2, gauge))
        self.assertNotEqual(bar.export(boundary)["C"], [0])
        expected = bar.export(bar.rule.xtimes(boundary, bar.state(3, canonical)))
        self.assertEqual(native.act(3, gauge, canonical)[0], expected)

    def test_degree_four_action_on_a_nonclosed_signed_gauge_a(self):
        native, bar, _ = c2_models(1, 0)
        gauge = state([1], [1], [1], [-2])
        canonical = native.zero(4)
        actual, _ = native.act(4, gauge, canonical)
        expected = bar.export(bar.rule.xtimes(bar.rule.d(bar.state(3, gauge)),
                                             bar.state(4, canonical)))
        self.assertTrue(same_class(bar, 4, actual, expected), (actual, expected))
        self.assertEqual({f: actual[f] for f in "ABC"}, {f: v for f, v in state([-2], [0], [1], [-3]).items() if f != "D"})

    def test_degree_five_pure_c_square_matches_complete_bar(self):
        native, bar, _ = c2_models()
        data = state([0], [0], [1], [0])
        self.assertEqual(native.kappa(5, data), native.zero(6))
        actual, _ = native.product(5, data, data)
        expected = bar.export(bar.rule.xtimes(bar.state(5, data), bar.state(5, data)))
        self.assertTrue(same_class(bar, 5, actual, expected), (actual, expected))
        self.assertEqual({f: actual[f] for f in "ABC"}, {f: [0] for f in "ABC"})

    def test_pure_c_products_use_the_direct_formula_up_to_a_coboundary(self):
        # A=B=0 products of fully legal degree-five states are corrected by
        # pure_c_gamma; the reference model's product differs by a D-coboundary.
        import compatible_sector
        for sign, omega in ((1, 0), (0, 1), (1, 1)):
            native, bar, _ = c2_models(sign, omega)
            states = [state([0], [0], [c], [d]) for c, d in ((1, 0), (1, 2), (1, -1))]
            states = [x for x in states if native.kappa(5, x) == native.zero(6)]
            if not states:
                continue
            pairs = [(left, right) for left in states for right in states]
            expected = [bar.export(bar.rule.xtimes(bar.state(5, left), bar.state(5, right)))
                        for left, right in pairs]
            with patch.object(native_module.closed, "gamma",
                              side_effect=AssertionError("production gamma on the pure sector")):
                for (left, right), reference in zip(pairs, expected):
                    actual, _ = native.product(5, left, right)
                    self.assertTrue(same_class(bar, 5, actual, reference), (sign, omega, actual, reference))
                    self.assertEqual(native.divide_left(5, left, actual), right)

    def test_a0_curvature_equals_the_production_j_exactly(self):
        # Theorem: g_k(0,B,C) = closed_ab_upper.J(0,B,C) as cochains for closed B.
        from fractions import Fraction
        import all_cochain_upper as upper
        for k in (5, 6):
            for sign, omega in ((0, 0), (1, 0), (0, 1), (1, 1)):
                for b_zero in (True, False):
                    tools = random_cochains(10 * k + 2 * sign + omega + (1 if b_zero else 0), k + 2)
                    p = tools["p"]
                    s, w = tools["twists"](sign, omega)
                    B = p.zero(k - 2) if b_zero else tools["closed_binary"](k - 2)
                    C = tools["binary"](k - 1)
                    triple = upper.Triple(p.zero(k - 3), B, C, True)
                    simplex = tuple(range(k + 3))
                    direct = native_module.a0_curvature(k, triple, s, w, b_zero)
                    production = native_module.closed.J(triple.A, B, C, s, w)
                    self.assertEqual(Fraction(direct(simplex)), Fraction(production(simplex)),
                                     (k, sign, omega, b_zero))

    def test_a0_states_never_evaluate_the_production_phase(self):
        native, _, _ = c2_models(1, 0, top=8)
        with patch.object(native_module.closed, "J", side_effect=AssertionError("production J at A=0")):
            self.assertEqual(native.kappa(6, state([0], [1], [0], [0])), native.zero(7))
            self.assertEqual(sorted(native.kappa(5, state([0], [0], [1], [0]))), ["A", "B", "C", "D"])

    def test_closed_prism_equals_the_prism_of_theta(self):
        # I Theta(delta(b l), Q_D(b l)) for closed b, in closed form.
        from fractions import Fraction
        import a0_high_gamma
        from cochains import Cochain as LocalCochain
        from compatible_sector import QD
        from lower_stacking import prism, pullback_interval
        for m in (3, 4):
            for sign, omega in ((0, 0), (1, 0), (0, 1), (1, 1)):
                tools = random_cochains(100 + 10 * m + 2 * sign + omega, m + 3)
                p = tools["p"]
                s, w = tools["twists"](sign, omega)
                b = tools["closed_binary"](m)
                simplex = tuple(range(m + 4))
                bi = pullback_interval(LocalCochain(m, b), True)
                si, wi = pullback_interval(LocalCochain(1, s)), pullback_interval(LocalCochain(2, w))
                raw = p.theta(a0_high_gamma.differential(bi).mod2(), QD(bi, si, wi), si, wi)
                expected = prism(LocalCochain(raw.degree, raw))(simplex)
                actual = a0_high_gamma.closed_prism(LocalCochain(m, b), LocalCochain(1, s),
                                                    LocalCochain(2, w))(simplex)
                self.assertEqual(Fraction(actual), Fraction(expected), (m, sign, omega))

    def test_pair_source_prisms_of_pulled_back_data_vanish(self):
        # phi(A,0,0) = -V_P(A): the prism of Theta on pulled-back data is zero.
        from fractions import Fraction
        import closed_a_upper as upper
        for sign, omega in ((0, 0), (1, 0), (0, 1), (1, 1)):
            tools = random_cochains(200 + 2 * sign + omega, 6)
            p = tools["p"]
            s, w = tools["twists"](sign, omega)
            A = p.ds(tools["integral"](1), s)
            simplex = tuple(range(7))
            phi = upper.phi(A, p.zero(3), p.zero(4), s, w)
            self.assertEqual(Fraction(phi(simplex)), -Fraction(upper.source_primitive(A, s, w)(simplex)),
                             (sign, omega))

    def test_pair_primitive_vanishes_with_one_zero_fiber(self):
        from fractions import Fraction
        import production_gamma5
        import production_gamma6
        for module, degree, dimension in ((production_gamma5, 2, 6), (production_gamma6, 3, 7)):
            for sign, omega in ((0, 0), (1, 1)):
                tools = random_cochains(300 + degree + sign, dimension)
                p = tools["p"]
                s, w = tools["twists"](sign, omega)
                A = p.ds(tools["integral"](degree - 1), s)
                primitive = module.constructor().primitive(A, p.zero(degree), s, w)
                self.assertEqual(Fraction(primitive(tuple(range(dimension + 1)))), 0)
                primitive = module.constructor().primitive(p.zero(degree), A, s, w)
                self.assertEqual(Fraction(primitive(tuple(range(dimension + 1)))), 0)
        # The degree-six correction of one nonzero A layer is admitted by the guard.
        import all_cochain_upper as upper
        p = upper.p
        legal = upper.Triple(p.zero(3), p.zero(4), p.zero(5), True)
        with patch.object(native_module, "A_STACKING", False):
            correction = native_module.gamma(6, legal, legal, True, False, (False, True), p.zero(1), p.zero(2))
        self.assertEqual(correction.degree, 7)

    def test_debug_phi_values_use_requested_vertices_and_are_immutable_in_cache(self):
        native, _, _ = c2_models()
        data = state([2], [0], [0], [9])
        samples = {f: [[i % 2 for i in range(n + 1)]]
                   for f, n in zip("ABCD", degrees(3))}
        request = dict(operation="phi_values", degree=3, state=data, simplices=samples)
        result = native.calculate(request)
        self.assertEqual(result["state"], data)
        result["state"]["A"][0] = 99
        self.assertEqual(native.calculate(request)["state"], data)

    def test_zero_predicate_pairs_with_resolution_chains_only(self):
        native, _, _ = c2_models()
        with patch.object(native, "simplices", side_effect=AssertionError("bar enumeration")):
            self.assertTrue(native.is_zero(native.lift(2, [0], True), False))
            self.assertFalse(native.is_zero(native.lift(2, [3], True), False))
            self.assertFalse(native.is_zero(native.lift(2, [1], False), True))
            self.assertTrue(native.is_zero(native.lift(2, [2], True), True))

    def test_bar_value_export_keeps_its_enumeration_budget(self):
        native, _, _ = c2_models()
        # A three-element multiplication table gives 2^n normalized simplices.
        native.mul = [[(a + b) % 3 for b in range(3)] for a in range(3)]
        native.size = 3
        native.max_debug_simplices = 1
        with self.assertRaises(TransferResourceLimit):
            native.simplices(2)

    def test_label_vertices_are_normalized_by_gap(self):
        bar = FiniteBar([[0, 1], [1, 0]], [0], [0])
        calls = []
        def transport(kind, degree, vertices):
            calls.append(tuple(vertices))
            return dict(status="computed", terms=[[0, 1, 1]] if kind == "f" else [])
        setup = dict(schema=1, vertexMode="labels", ranks=[1] * 8,
            ordinary=[[[0 if n % 2 == 0 else 2]] for n in range(7)],
            signed=[[[0 if n % 2 == 0 else 2]] for n in range(7)],
            g=[[[[1, 1, list(bar.simplices(n)[0])]]] for n in range(8)], s=[0], omega=[0])
        native = TransferredModel(setup, transport)
        self.assertIsNone(native.normalize((5, 5, 2)))
        self.assertEqual(native.normalize((7, 3)), (7, 3))
        self.assertEqual(native.lift(1, [1], True)((7, 3)), 1)
        self.assertEqual(calls, [(7, 3)])
        with self.assertRaises(ValueError):
            native.simplices(1)

    def test_degree_two_products_and_actions_match_complete_bar(self):
        for sign, omega in ((0, 0), (1, 0), (0, 1), (1, 1)):
            native, bar, _ = c2_models(sign, omega)
            for left, right in ((state([], [1], [0], [3]), state([], [1], [1], [-2])),
                                (state([], [0], [1], [0]), state([], [1], [0], [5]))):
                with self.subTest(sign=sign, omega=omega, left=left, right=right):
                    expected = bar.export(bar.rule.xtimes(bar.state(2, left), bar.state(2, right)))
                    self.assertEqual(native.product(2, left, right)[0], expected)
                    # Native curvature is exact through its first nonzero
                    # layer and zero-filled afterwards.
                    complete = bar.export(bar.rule.d(bar.state(2, left)))
                    first = next((i for i, f in enumerate("ABCD") if any(complete[f])), 3)
                    expected = {f: complete[f] if i <= first else [0] * len(complete[f])
                                for i, f in enumerate("ABCD")}
                    self.assertEqual(native.kappa(2, left), expected)
            gauge = state([], [], [1], [4])
            canonical = state([], [1], [0], [0])
            boundary = bar.rule.d(bar.state(1, gauge))
            expected = bar.export(bar.rule.xtimes(boundary, bar.state(2, canonical)))
            self.assertEqual(native.act(2, gauge, canonical)[0], expected)

    def test_invalid_state_coordinates_are_rejected(self):
        native, _, _ = c2_models()
        for bad in ([], [0, 0], [1.5], [True]):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                native.kappa(3, state([0], [0], [0], bad))
        with self.assertRaises(ValueError):
            native.kappa(3, state([0], [2], [0], [0]))

    def test_full_product_rejects_nonflat_inputs_but_partial_layer_is_available(self):
        native, _, _ = c2_models(1, 0)
        data = state([1], [0], [0], [0])
        request = dict(operation="xtimes", degree=3, state=data, other=native.zero(3))
        with self.assertRaises(ValueError):
            native.calculate(request)
        self.assertEqual(native.calculate(dict(request, upto=0))["state"]["A"], [1])

    def test_nonlinear_product_cache_preserves_distinct_d_carries(self):
        native, bar, _ = c2_models()
        left = state([0], [1], [0], [7])
        right = state([-1], [1], [0], [-4])
        request = dict(operation="xtimes", degree=3, state=left, other=right)
        first = native.calculate(request)["state"]
        count = len(native._cores)
        shifted = dict(request, state=dict(left, D=[18]), other=dict(right, D=[-9]))
        with patch.object(native, "product", side_effect=AssertionError("nonlinear cache miss")):
            actual = native.calculate(shifted)["state"]
        self.assertEqual(len(native._cores), count)
        self.assertEqual(actual, dict(first, D=[first["D"][0] + 6]))
        self.assertEqual(actual, bar.export(bar.rule.xtimes(
            bar.state(3, shifted["state"]), bar.state(3, shifted["other"]))))
        first["D"][0] = 999
        self.assertEqual(native.calculate(shifted)["state"], actual)
        self.assertIs(native.phi(3, left)[2], native.phi(3, shifted["state"])[2])

    def test_signed_curvature_cache_and_action_cache_preserve_d_carries(self):
        native, _, _ = c2_models(1, 1)
        first = native.kappa(3, state([0], [0], [0], [9]))
        count = len(native._kappas)
        self.assertEqual(first["D"], [-18])
        self.assertEqual(native.kappa(3, state([0], [0], [0], [-3]))["D"], [6])
        self.assertEqual(len(native._kappas), count)
        native, bar, _ = c2_models()
        gauge = state([], [1], [0], [3])
        canonical = state([0], [0], [1], [7])
        request = dict(operation="act", degree=3, gauge=gauge, state=canonical)
        first = native.calculate(request)["state"]
        shifted = dict(request, gauge=dict(gauge, D=[-2]), state=dict(canonical, D=[11]))
        with patch.object(native, "act", side_effect=AssertionError("nonlinear cache miss")):
            actual = native.calculate(shifted)["state"]
        self.assertEqual(actual["D"], [first["D"][0] - 6])
        expected = bar.export(bar.rule.xtimes(bar.rule.d(bar.state(2, shifted["gauge"])),
                                             bar.state(3, shifted["state"])))
        self.assertEqual(actual, expected)

    def test_degree_six_curvature_stops_before_the_d_layer_off_the_legal_locus(self):
        # A nonclosed A (sign 0: delta_s A = 2A) reports its first obstruction
        # only; the degree-six D-layer formula is never requested there.
        native, _, _ = c2_models(0, 1, top=8)
        with patch.object(native_module, "g", side_effect=AssertionError("D-layer formula off the legal locus")):
            self.assertEqual(native.kappa(6, state([1], [1], [1], [5])),
                             state([2], [0], [0], [0]))
        request = dict(operation="xtimes", degree=6, state=state([1], [0], [0], [0]),
                       other=native.zero(6))
        with self.assertRaises(ValueError):
            native.calculate(request)
        self.assertEqual(native.calculate(dict(request, upto=0))["state"]["A"], [1])
        with self.assertRaises(ValueError):
            native.calculate(dict(operation="d", degree=7, state=native.zero(6)))

    def test_degree_six_section_branch_is_reported_as_unresolved(self):
        import all_cochain_upper as upper
        p = upper.p
        triple = upper.Triple(p.zero(3), p.zero(4), p.zero(5), False)
        s, omega = p.zero(1), p.zero(2)
        with self.assertRaises(native_module.SectionBranchRequired):
            native_module.g(6, triple, s, omega)
        with self.assertRaises(native_module.SectionBranchRequired):
            native_module.gamma(6, triple, triple, False, False, True, s, omega)
        with self.assertRaises(ValueError):
            native_module.g(6, upper.Triple(p.zero(2), p.zero(3), p.zero(4), True), s, omega)
        # A nonzero A layer in a degree-six stacking correction is refused
        # before any formula is evaluated, unless the override is set.
        legal = upper.Triple(p.zero(3), p.zero(4), p.zero(5), True)
        with patch.object(native_module, "A_STACKING", False), \
             self.assertRaises(native_module.PairSourceLimit):
            native_module.gamma(6, legal, legal, True, False, False, s, omega)
        self.assertTrue(issubclass(native_module.PairSourceLimit, native_module.NativeDegreeSixLimit))
        # The worker loop reports the branch as unresolved, not as an error.
        native, _, _ = c2_models(1, 0, top=8)
        lines = [json.dumps(dict(operation="setup")),
                 json.dumps(dict(operation="d", degree=6, state=native.zero(6)))]
        out = io.StringIO()
        with patch.object(extension_transfer, "TransferredModel", return_value=native), \
             patch.object(native, "calculate", side_effect=native_module.SectionBranchRequired("section")), \
             patch.object(extension_transfer.acceleration, "persist"), \
             patch.object(extension_transfer.acceleration, "flush"), \
             patch.object(extension_transfer.sys, "stdin", io.StringIO("\n".join(lines) + "\n")), \
             contextlib.redirect_stdout(out):
            extension_transfer.serve()
        answers = [json.loads(line) for line in out.getvalue().splitlines() if line.strip()]
        self.assertEqual(answers[0]["status"], "computed")
        self.assertEqual(answers[1]["status"], "unresolved")
        self.assertEqual(answers[1]["exception"], "SectionBranchRequired")

    @unittest.skipUnless(os.environ.get("FERMIONAHSS_SLOW_TRANSFER_TESTS") == "1",
                         "degree-six formulas on the complete C2 bar; opt in explicitly")
    def test_degree_six_operations_match_the_complete_bar_section_model(self):
        from test_extension_degree_six import c2_data
        # The omega twist alone has no nonzero flat state among the candidates
        # (its C=1 state has a nonzero D curvature and its A layer is not
        # closed), so the sign twist is compared with and without omega. The
        # curvature of nonzero-A states (J_6 with nonzero A, minutes each) is
        # compared for the sign twist alone; with omega the A=0 candidates
        # exercise the product and the action.
        for sign, omega in ((1, 0), (1, 1)):
            native, bar, _ = c2_models(sign, omega, top=8)
            configure_degree_six(bar, c2_data(signed=bool(sign)))
            candidates = ((0, 0, 1), (1, 0, 0), (1, 1, 1), (2, 0, 1)) if omega == 0 else ((0, 0, 1), (0, 1, 0), (0, 1, 1))
            flat = [data for a, b, c in candidates
                    for data in (state([a], [b], [c], [0]),)
                    if not any(any(v) for v in native.kappa(6, data).values())]
            self.assertTrue(flat)
            for data in flat:
                complete = bar.export(bar.rule.d(bar.state(6, data)))
                self.assertEqual(complete, native.zero(7))
            # Two nonzero A layers in a degree-six product are refused by the
            # guard on the nested pair source; with one nonzero A layer the
            # pair primitive vanishes and the product is compared with the
            # reference model. Products and actions of zero-A states are
            # compared on the flat states with a zero A layer (none besides
            # zero for the omega twist without sign).
            a_zero = [data for data in flat if data["A"] == [0]]
            nonzero = [data for data in flat if data["A"] != [0]]
            for data in nonzero:
                with self.assertRaises(native_module.PairSourceLimit):
                    native.product(6, data, data)
            for data in nonzero[:1]:
                left = state([0], [0], [0], [3])
                actual, _ = native.product(6, left, data)
                expected = bar.export(bar.rule.xtimes(bar.state(6, left), bar.state(6, data)))
                self.assertTrue(same_class(bar, 6, actual, expected), (actual, expected))
            if not a_zero:
                continue
            left = dict(a_zero[0], D=[3])
            right = dict(a_zero[-1], D=[-1])
            actual, _ = native.product(6, left, right)
            expected = bar.export(bar.rule.xtimes(bar.state(6, left), bar.state(6, right)))
            self.assertTrue(same_class(bar, 6, actual, expected), (actual, expected))
            self.assertEqual(native.divide_left(6, left, actual), right)
            gauge = state([0], [1], [0], [2])
            actual, _ = native.act(6, gauge, a_zero[0])
            expected = bar.export(bar.rule.xtimes(bar.rule.d(bar.state(5, gauge)),
                                                 bar.state(6, a_zero[0])))
            self.assertTrue(same_class(bar, 6, actual, expected), (actual, expected))

    def test_three_primary_phase_in_input_degree_two(self):
        # (2/3) lift(P^1_s rho_3 A) with P^1 the cube modulo three in degree two;
        # its failure of additivity is integral, and it enters the degree-two
        # production phase.
        from fractions import Fraction
        import mod3_power
        import low_phases
        tools = random_cochains(500, 6)
        p = tools["p"]
        simplex = tuple(range(7))
        for sign in (0, 1):
            s, w = tools["twists"](sign, 0)
            A = p.ds(tools["integral"](1), s)
            Ap = p.ds(tools["integral"](1), s)
            phase = mod3_power.tertiary_three_primary_phase(A, s)
            self.assertEqual(phase.degree, 6)
            power = mod3_power.reduced_power_1(A, s)
            self.assertEqual(Fraction(phase(simplex)), Fraction(2 * power(simplex), 3))
            self.assertIn(power(simplex), (0, 1, 2))
            # The phase is defined only modulo integers up to the universal pair
            # primitive: its failure of additivity is a third-integral cochain.
            total = mod3_power.tertiary_three_primary_phase(A + Ap, s)
            other = mod3_power.tertiary_three_primary_phase(Ap, s)
            self.assertIn(Fraction(phase(simplex) + other(simplex) - total(simplex)).denominator, (1, 3))
        self.assertEqual(mod3_power.tertiary_three_primary_phase(p.zero(1), p.zero(1)).degree, 5)
        # A cocycle with a nonzero cube modulo three: the constant 1 on 2-faces.
        A = p.Cochain(2, lambda face: 1)
        s = p.zero(1)
        self.assertEqual(Fraction(mod3_power.tertiary_three_primary_phase(A, s)(simplex)), Fraction(2, 3))
        # The degree-two production phase contains the term.
        b, c, w = p.zero(3), p.zero(4), p.zero(2)
        full = low_phases.build_phase(2, A, b, c, s, w)
        with patch.object(mod3_power, "tertiary_three_primary_phase", lambda A, s: p.zero(A.degree + 4)):
            without = low_phases.build_phase(2, A, b, c, s, w)
        self.assertEqual(Fraction(full(simplex)) - Fraction(without(simplex)), Fraction(2, 3))

    def test_layer_limited_operations_agree_with_their_full_layers(self):
        # A layer-limited action, division and differential equal the full
        # ones through the requested layer and are zero above it.
        native, _, _ = c2_models(1, 1)
        canonical = state([1], [0], [1], [-3])
        gauge = state([0], [1], [0], [2])
        full, _ = native.act(4, gauge, canonical)
        for upto in (0, 1, 2):
            partial, _ = native.act(4, gauge, canonical, upto)
            for layer, f in enumerate("ABCD"):
                self.assertEqual(partial[f], full[f] if layer <= upto else [0] * len(full[f]), (upto, f))
        left = state([0], [1], [0], [7])
        total = state([1], [1], [1], [-4])
        quotient = native.divide_left(4, left, total)
        for upto in (1, 2):
            partial = native.divide_left(4, left, total, upto=upto)
            for layer, f in enumerate("ABCD"):
                self.assertEqual(partial[f], quotient[f] if layer <= upto else [0] * len(quotient[f]), (upto, f))
        curvature = native.kappa(4, state([1], [1], [1], [5]))
        limited = native.kappa(4, state([1], [1], [1], [5]), 1)
        self.assertEqual(limited["A"], curvature["A"])
        self.assertEqual(limited["B"], curvature["B"])
        # Requests carry the layer index; the answers above it are zero.
        request = dict(operation="act", degree=4, state=canonical, gauge=gauge, upto=1)
        answer = native.calculate(request)["state"]
        self.assertEqual({f: answer[f] for f in "AB"}, {f: full[f] for f in "AB"})
        self.assertEqual(answer["D"], [0] * len(full["D"]))
        answer = native.calculate(dict(operation="d", degree=4, state=state([1], [1], [1], [5]), upto=1))["state"]
        self.assertEqual(answer["C"], [0] * len(answer["C"]))
        answer = native.calculate(dict(operation="divide_left", degree=4, state=left, other=total, upto=2))["state"]
        self.assertEqual({f: answer[f] for f in "ABC"}, {f: quotient[f] for f in "ABC"})
        with self.assertRaises(ValueError):
            native.calculate(dict(operation="act", degree=4, state=canonical, gauge=gauge, upto=4))

    def test_cached_left_division_restores_total_minus_left_d(self):
        native, _, _ = c2_models()
        left = state([0], [1], [0], [7])
        right = state([-1], [1], [0], [-4])
        total = native.calculate(dict(operation="xtimes", degree=3,
                                      state=left, other=right))["state"]
        request = dict(operation="divide_left", degree=3, state=left, other=total)
        self.assertEqual(native.calculate(request)["state"], right)
        shifted = dict(request, state=dict(left, D=[18]),
                       other=dict(total, D=[total["D"][0] + 6]))
        with patch.object(native, "divide_left", side_effect=AssertionError("nonlinear cache miss")):
            actual = native.calculate(shifted)["state"]
        self.assertEqual(actual, dict(right, D=[-9]))
        self.assertEqual(native.calculate(dict(operation="xtimes", degree=3,
            state=shifted["state"], other=actual))["state"], shifted["other"])

    def test_three_local_operations_on_the_c2_setup(self):
        # The two-layer model of the prime three (degree five, B=C=0): curvature,
        # product, division and gauge action, and their worker requests.
        native, bar, _ = c2_models()
        local = native.three_local
        A = state([1], [0], [0], [0])
        self.assertEqual(local.kappa(5, A), native.zero(6))
        product = local.product(5, A, A)
        # gamma(A,A) = -2 (A A A + A A A) = -4 A^3 on the one-cell resolution.
        self.assertEqual(product, state([2], [0], [0], [-4]))
        self.assertEqual(local.divide_left(5, A, product), A)
        # A gauge (u,w) of degree four acts through its boundary (delta u, delta w):
        # delta u = 2u in degree one and delta w = 2w in degree five of the
        # untwisted C2 resolution, and the D layer receives gamma(delta u, A).
        gauge = state([1], [0], [0], [2])
        self.assertEqual(local.act(5, gauge, A), state([3], [0], [0], [4 - 12]))
        self.assertEqual(native.calculate(dict(operation="xtimes", degree=5, state=A,
            other=A, prime=3))["state"], product)
        self.assertEqual(native.calculate(dict(operation="divide_left", degree=5, state=A,
            other=product, prime=3))["state"], A)
        self.assertEqual(native.calculate(dict(operation="act", degree=5, state=A,
            gauge=gauge, prime=3))["state"], state([3], [0], [0], [-8]))
        self.assertEqual(native.calculate(dict(operation="d", degree=5, state=product,
            prime=3))["state"], native.zero(6))
        self.assertEqual(native.calculate(dict(operation="d", degree=4, state=gauge,
            prime=3))["state"], state([2], [0], [0], [4]))
        # The complete model is untouched by the localized requests.
        self.assertEqual(native.kappa(5, A)["A"], [0])
        with self.assertRaises(ValueError):
            local.product(5, state([1], [1], [0], [0]), A)
        with self.assertRaises(ValueError):
            native.calculate(dict(operation="d", degree=5, state=A, prime=5))
        with self.assertRaises(ValueError):
            native.calculate(dict(operation="xtimes", degree=4, state=state([1], [0], [0], [0]),
                other=state([1], [0], [0], [0]), prime=3))

    def test_three_local_correction_is_the_polarization_up_to_a_cup_one_coboundary(self):
        # gamma(A,A') = -2 (AAA' + AA'A') differs from the polarization of
        # Omega'(A) = (2/3) AAA by (2/3) delta Xi with
        # Xi = 2 A(A u1 A') + (A u1 A')A + 2 (A u1 A')A' + A'(A u1 A'), the
        # Hirsch identity A'A - AA' = delta(A u1 A') applied word by word.
        from fractions import Fraction
        from low_phases import transported
        tools = random_cochains(11, 6)
        p = tools["p"]
        s = p.zero(1)
        simplex = tuple(range(7))
        for trial in range(3):
            A = p.ds(tools["integral"](1), s)
            Ap = p.ds(tools["integral"](1), s)
            word = lambda x, y, z: transported(transported(x, y, s), z, s)
            polarization = p.scale(word(A, A, A) + word(Ap, Ap, Ap) - word(A + Ap, A + Ap, A + Ap),
                                   Fraction(2, 3))
            gamma = p.scale(word(A, A, Ap) + word(A, Ap, Ap), -2)
            self.assertEqual(Fraction(gamma(simplex)).denominator, 1)
            c1 = p.cup(A, Ap, 1, integral=True)
            face = tuple(range(5))
            self.assertEqual((transported(Ap, A, s) - transported(A, Ap, s))(face), p.ds(c1, s)(face))
            xi = (p.scale(transported(A, c1, s), 2) + transported(c1, A, s)
                  + p.scale(transported(c1, Ap, s), 2) + transported(Ap, c1, s))
            self.assertEqual(Fraction(gamma(simplex)) - Fraction(polarization(simplex)),
                             Fraction(2, 3) * Fraction(p.ds(xi, s)(simplex)))

    def test_degree_six_three_local_primitives(self):
        # The natural primitives of the degree-three reduced power modulo three,
        # for the untwisted and a twisted sign system: delta phi is the cross
        # effect of P^1 on cocycles, delta chi is P^1 of a coboundary, and the
        # brackets of the degree-six stacking correction and gauge boundary are
        # divisible by three.
        from fractions import Fraction
        from itertools import combinations
        import mod3_power as m
        tools = random_cochains(2026, 8)
        for sign in (0, 1):
            s, _ = tools["twists"](sign, 0)
            for trial in range(2):
                a = m.signed_coboundary(tools["integral"](2), s)
                ap = m.signed_coboundary(tools["integral"](2), s)
                u = tools["integral"](2)
                phi = m.cross_effect_primitive(a, ap, s)
                dphi = m.signed_coboundary(phi, s)
                chi = m.coboundary_primitive(u, s)
                dchi = m.signed_coboundary(chi, s)
                du = m.signed_coboundary(u, s)
                power_du = m.reduced_power_1(du, s)
                powers = [m.reduced_power_1(x, s) for x in (a, ap, a + ap)]
                for z in combinations(range(9), 8):
                    cross = (powers[2](z) - powers[0](z) - powers[1](z)) % 3
                    self.assertEqual(dphi(z) % 3, cross)
                    self.assertEqual(Fraction(2 * (powers[0](z) + powers[1](z) - powers[2](z) + dphi(z)), 3).denominator, 1)
                    self.assertEqual(dchi(z) % 3, power_du(z))
                    self.assertEqual(Fraction(2 * (power_du(z) - dchi(z)), 3).denominator, 1)
        self.assertEqual(len(m.cyclic_diagonal(3, 6)), 448)
        self.assertEqual(len(m.diagonal_terms(3, 6, (3, 3, 3), rotation=(0, 2, 1))), 48)

    def test_degree_six_three_local_operations_on_the_twisted_c2_setup(self):
        # On the sign-twisted C2 resolution A=[1] is a degree-three cocycle;
        # products, division, gauge boundaries and actions of the degree-six
        # two-layer model are flat and consistent, and the requests reach it.
        native, bar, _ = c2_models(sign=1, top=8)
        local = native.three_local
        flat = lambda k, data: not any(any(v) for v in local.kappa(k, data).values())
        A = state([1], [0], [0], [3])
        self.assertTrue(flat(6, A))
        product = local.product(6, A, A)
        self.assertEqual(product["A"], [2])
        self.assertTrue(flat(6, product))
        self.assertEqual(local.divide_left(6, A, product), A)
        gauge = state([1], [0], [0], [1])
        boundary = local.boundary(6, gauge)
        self.assertEqual(boundary["A"], native.coboundary(2, [1], True))
        self.assertTrue(flat(6, boundary))
        acted = local.act(6, gauge, A)
        self.assertTrue(flat(6, acted))
        self.assertEqual(acted, local.product(6, boundary, A))
        self.assertEqual(native.calculate(dict(operation="d", degree=5, state=gauge, prime=3,
            role="gauge"))["state"], boundary)
        self.assertEqual(native.calculate(dict(operation="xtimes", degree=6, state=A, other=A,
            prime=3))["state"], product)
        self.assertEqual(native.calculate(dict(operation="act", degree=6, state=A, gauge=gauge,
            prime=3))["state"], acted)
        self.assertEqual(native.calculate(dict(operation="divide_left", degree=6, state=A,
            other=product, prime=3))["state"], A)
        self.assertEqual(native.calculate(dict(operation="d", degree=6, state=product,
            prime=3))["state"], native.zero(7))


if __name__ == "__main__":
    unittest.main()
