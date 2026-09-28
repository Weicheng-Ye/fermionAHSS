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
from extension_transfer import TransferredModel, TransferResourceLimit, degrees
import extension_native_six as six


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


class ExtensionTransferTests(unittest.TestCase):
    def test_twists_do_not_request_bar_values_during_setup(self):
        _, _, calls = c2_models(1, 1)
        self.assertEqual(calls, [])

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
        self.assertEqual(actual, expected)
        self.assertEqual(actual, state([-2], [0], [1], [-3]))

    def test_degree_five_pure_c_square_matches_complete_bar(self):
        native, bar, _ = c2_models()
        data = state([0], [0], [1], [0])
        self.assertEqual(native.kappa(5, data), native.zero(6))
        actual, _ = native.product(5, data, data)
        expected = bar.export(bar.rule.xtimes(bar.state(5, data), bar.state(5, data)))
        self.assertEqual(actual, expected)
        self.assertEqual(actual, native.zero(5))

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
        with patch.object(six, "g", side_effect=AssertionError("D-layer formula off the legal locus")):
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
        with self.assertRaises(six.SectionBranchRequired):
            six.g(triple, s, omega)
        with self.assertRaises(six.SectionBranchRequired):
            six.gamma(triple, triple, False, False, True, s, omega)
        with self.assertRaises(ValueError):
            six.g(upper.Triple(p.zero(2), p.zero(3), p.zero(4), True), s, omega)
        # A nonzero A layer in a degree-six stacking correction is refused
        # before any formula is evaluated, unless the override is set.
        legal = upper.Triple(p.zero(3), p.zero(4), p.zero(5), True)
        with patch.object(six, "A_STACKING", False), self.assertRaises(six.PairSourceLimit):
            six.gamma(legal, legal, True, False, False, s, omega)
        self.assertTrue(issubclass(six.PairSourceLimit, six.NativeDegreeSixLimit))
        # The worker loop reports the branch as unresolved, not as an error.
        native, _, _ = c2_models(1, 0, top=8)
        lines = [json.dumps(dict(operation="setup")),
                 json.dumps(dict(operation="d", degree=6, state=native.zero(6)))]
        out = io.StringIO()
        with patch.object(extension_transfer, "TransferredModel", return_value=native), \
             patch.object(native, "calculate", side_effect=six.SectionBranchRequired("section")), \
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
        for sign, omega in ((1, 0), (0, 1)):
            native, bar, _ = c2_models(sign, omega, top=8)
            configure_degree_six(bar, c2_data(signed=bool(sign)))
            flat = [data for a, b, c in ((0, 0, 1), (1, 0, 0), (1, 1, 1), (2, 0, 1))
                    for data in (state([a], [b], [c], [0]),)
                    if not any(any(v) for v in native.kappa(6, data).values())]
            self.assertTrue(flat)
            for data in flat:
                complete = bar.export(bar.rule.d(bar.state(6, data)))
                self.assertEqual(complete, native.zero(7))
            left = dict(flat[0], D=[3])
            right = dict(flat[-1], D=[-1])
            actual, _ = native.product(6, left, right)
            expected = bar.export(bar.rule.xtimes(bar.state(6, left), bar.state(6, right)))
            self.assertEqual(actual, expected)
            self.assertEqual(native.divide_left(6, left, actual), right)
            gauge = state([1], [1], [0], [2])
            actual, _ = native.act(6, gauge, flat[0])
            expected = bar.export(bar.rule.xtimes(bar.rule.d(bar.state(5, gauge)),
                                                 bar.state(6, flat[0])))
            self.assertEqual(actual, expected)

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


if __name__ == "__main__":
    unittest.main()
