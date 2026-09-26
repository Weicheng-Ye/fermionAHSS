"""Bounded prerequisites for transferring the bundled stacking operations.

These check literal complete-bar identities, not naturality, completeness of
transferred gauges, or an implementation of the proposed transferred model.
The larger off-shell fourth-layer checks are opt-in because their universal
prisms are expensive: FERMIONAHSS_SLOW_TRANSFER_TESTS=1 python3 -m unittest
discover -s python -p test_extension_transfer_contracts.py -v

Copyright (c) 2026 koAHSS contributors; MIT license.
"""
import os
import unittest

from extension_worker import FiniteBar


def c2(sign, omega):
    return FiniteBar([[0, 1], [1, 0]], [sign], [omega])


def state(a, b, c, d):
    return dict(A=a, B=b, C=c, D=d)


class TransferFormulaContractTests(unittest.TestCase):
    def assert_square_zero(self, bar, degree, data, expected_boundary):
        boundary = bar.rule.d(bar.state(degree, data))
        self.assertEqual(bar.export(boundary), expected_boundary)
        second = bar.export(bar.rule.d(boundary))
        self.assertEqual(second, {
            name: [0] * bar.dimension(n)
            for name, n in zip("ABCD", (degree - 1, degree, degree + 1,
                                        degree + 3))
        })

    def assert_multiplicative(self, bar, degree, left, right):
        x, y = bar.state(degree, left), bar.state(degree, right)
        actual = bar.export(bar.rule.d(bar.rule.xtimes(x, y)))
        expected = bar.export(bar.rule.xtimes(bar.rule.d(x), bar.rule.d(y)))
        self.assertEqual(actual, expected)
        return actual

    def test_degree_two_square_zero_retains_nonzero_signed_carries(self):
        self.assert_square_zero(c2(0, 0), 2,
            state([], [0], [0], [-2]),
            state([0], [0], [0], [-4]))
        self.assert_square_zero(c2(1, 1), 2,
            state([], [1], [0], [0]),
            state([0], [0], [1], [0]))

    def test_degree_three_square_zero_on_nonclosed_integral_a(self):
        self.assert_square_zero(c2(1, 0), 3,
            state([1], [1], [1], [-2]),
            state([-2], [0], [1], [-3]))

    def test_degree_four_square_zero_on_legal_pair_with_nonzero_c_residual(self):
        self.assert_square_zero(c2(1, 0), 4,
            state([0], [1], [1], [-2]),
            state([0], [0], [1], [0]))

    def test_degree_two_multiplicativity_with_three_twists(self):
        cases = (
            (0, 0, state([], [0], [1], [2]), state([], [1], [0], [-1]), 2),
            (1, 0, state([], [1], [0], [1]), state([], [1], [0], [1]), 0),
            (1, 1, state([], [1], [0], [1]), state([], [1], [0], [2]), 0),
        )
        for sign, omega, left, right, expected_d in cases:
            with self.subTest(sign=sign, omega=omega):
                actual = self.assert_multiplicative(c2(sign, omega), 2, left, right)
                self.assertEqual(actual["D"], [expected_d])

    def test_degree_three_multiplicativity_with_mixed_integral_a(self):
        self.assert_multiplicative(c2(0, 0), 3,
            state([0], [1], [0], [1]),
            state([-1], [1], [0], [-1]))

    def test_zero_is_a_literal_two_sided_unit_in_degree_three(self):
        bar = c2(0, 0)
        data = state([-1], [1], [0], [-3])
        x = bar.state(3, data)
        zero = bar.state(3, state([0], [0], [0], [0]))
        self.assertEqual(bar.export(bar.rule.d(zero)),
                         state([0], [0], [0], [0]))
        self.assertEqual(bar.export(bar.rule.xtimes(zero, x)), data)
        self.assertEqual(bar.export(bar.rule.xtimes(x, zero)), data)

    @unittest.skipUnless(os.environ.get("FERMIONAHSS_SLOW_TRANSFER_TESTS") == "1",
                         "off-shell degree-four universal prism; opt in explicitly")
    def test_degree_four_square_zero_on_nonclosed_integral_a(self):
        self.assert_square_zero(c2(0, 0), 4,
            state([-1], [0], [0], [0]),
            state([-2], [0], [0], [-2]))


if __name__ == "__main__":
    unittest.main()
