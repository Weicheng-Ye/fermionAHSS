"""Exact worker regression on the C2 comparison and its protocol contracts.

The C2 native and normalized-bar coordinates coincide; the separate GAP
transport tests exercise sparse nonidentity comparisons. These tests do not
assert completeness of the transferred gauge relation.
Copyright (c) 2026 koAHSS contributors; MIT license.
"""
import unittest
from unittest.mock import patch

from extension_worker import FiniteBar
from extension_transfer import TransferredModel, TransferResourceLimit, degrees


def c2_models(sign=0, omega=0):
    bar = FiniteBar([[0, 1], [1, 0]], [sign], [omega])
    setup = dict(schema=1, multiplication=bar.mul, ranks=[1] * 8,
        ordinary=[[[0 if n % 2 == 0 else 2]] for n in range(7)],
        signed=[[[(-2 if n % 2 == 0 else 0) if sign else
                  (0 if n % 2 == 0 else 2)]] for n in range(7)],
        g=[[[[1, 1, list(bar.simplices(n)[0])]]] for n in range(8)],
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

    def test_global_zero_predicate_rejects_a_budget_before_sampling(self):
        native, _, _ = c2_models()
        # A three-element multiplication table gives 2^n normalized simplices.
        native.mul = [[(a + b) % 3 for b in range(3)] for a in range(3)]
        native.size = 3
        native.max_flags = 1
        with self.assertRaises(TransferResourceLimit):
            native.simplices(2)

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
