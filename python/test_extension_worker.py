"""The nonlinear cache must preserve exact integral D carries."""
import unittest
from extension_worker import FiniteBar


class ExtensionWorkerTests(unittest.TestCase):
    def test_invalid_d_coordinates_are_rejected_before_caching(self):
        bar = FiniteBar([[0,1],[1,0]], [0], [0])
        zero = dict(A=[], B=[], C=[], D=[0])
        for bad in ([], [0,0], [0.5], [0.0], [True]):
            for side in ("state", "other"):
                request = dict(operation="xtimes", degree=0, state=zero, other=zero)
                request[side] = dict(zero,D=bad)
                with self.assertRaises(ValueError):
                    bar.calculate(request)

    def test_stacking_shares_nonlinear_part_without_dropping_carries(self):
        bar = FiniteBar([[0,1],[1,0]], [0], [0])
        x = dict(A=[1], B=[0], C=[1], D=[7])
        y = dict(A=[0], B=[1], C=[0], D=[-4])
        actual = bar.calculate(dict(operation="xtimes", degree=3, state=x, other=y))
        expected = bar.export(bar.rule.xtimes(bar.state(3,x), bar.state(3,y)))
        self.assertEqual(actual, expected)
        before = len(bar._answers)
        shifted = bar.calculate(dict(operation="xtimes", degree=3,
                                     state=dict(x,D=[18]), other=dict(y,D=[-9])))
        self.assertEqual(shifted["D"], [actual["D"][0]+6])
        self.assertEqual(len(bar._answers), before)
        actual["A"][0] = 999
        self.assertEqual(bar.calculate(dict(operation="xtimes", degree=3,
                                          state=x, other=y))["A"], [1])

    def test_signed_differential_carry_is_exact(self):
        bar = FiniteBar([[0,1],[1,0]], [1], [1])
        state = dict(A=[0], B=[0], C=[0], D=[9])
        actual = bar.calculate(dict(operation="d", degree=3, state=state))
        expected = bar.export(bar.rule.d(bar.state(3,state)))
        self.assertEqual(actual, expected)
        self.assertEqual(actual["D"], [-18])
        shifted = bar.calculate(dict(operation="d", degree=3, state=dict(state,D=[-3])))
        self.assertEqual(shifted["D"], [6])
        self.assertEqual(len(bar._answers), 1)


if __name__ == "__main__":
    unittest.main()
