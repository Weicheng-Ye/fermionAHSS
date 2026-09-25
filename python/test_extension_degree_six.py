"""Exact finite-bar checks of the bundled degree-six section adapter."""
import unittest

from extension_worker import FiniteBar
from extension_degree_six import configure_degree_six, DegreeSixResourceLimit


def c2_data(signed=False, limit=4096):
    return dict(delta=[[0 if signed else 2]], basis=[[1]], inverse=[[1]],
                cycleRank=int(signed), dimB=1, dimC=1, deltaB=[0], deltaC=[0],
                maxSectionCandidates=limit)


class DegreeSixTests(unittest.TestCase):
    def test_nonzero_off_shell_A_and_successor(self):
        bar = FiniteBar([[0, 1], [1, 0]], [0], [1])
        section = configure_degree_six(bar, c2_data())

        def state(a, b, c, d):
            return bar.state(6, {name: [value] for name, value in zip("ABCD", (a, b, c, d))})

        x, y = state(1, 0, 0, 0), state(1, 1, 0, 0)
        dx, dy = bar.rule.d(x), bar.rule.d(y)
        total = bar.rule.xtimes(x, y)
        self.assertEqual(bar.export(total), dict(A=[2], B=[0], C=[1], D=[0]))
        left = bar.export(bar.rule.d(total))
        right = bar.export(bar.rule.xtimes(dx, dy))
        self.assertEqual(left, right)
        self.assertEqual(left, dict(A=[4], B=[0], C=[0], D=[1]))
        self.assertEqual(bar.export(bar.rule.d(dx)), dict(A=[0], B=[0], C=[0], D=[0]))
        self.assertEqual(section.section_search_exponent, 1)

    def test_zero_dimensional_positive_cochains(self):
        bar = FiniteBar([[0]], [], [])
        data = dict(delta=[], basis=[], inverse=[], cycleRank=0, dimB=0, dimC=0,
                    deltaB=[], deltaC=[])
        section = configure_degree_six(bar, data)
        empty = dict(A=[], B=[], C=[], D=[])
        x = bar.state(6, empty)
        self.assertEqual(bar.export(bar.rule.d(x)), empty)
        self.assertEqual(bar.export(bar.rule.xtimes(x, x)), empty)
        self.assertEqual(section.model.current_section(section.model.zero_image), ((), 0, 0))

    def test_wrong_basis_differential_is_rejected(self):
        bar = FiniteBar([[0, 1], [1, 0]], [0], [0])
        with self.assertRaisesRegex(ValueError, "differs from the bar basis"):
            configure_degree_six(bar, c2_data(signed=True))

    def test_section_limit_preserves_pointed_zero(self):
        bar = FiniteBar([[0, 1], [1, 0]], [1], [1])
        section = configure_degree_six(bar, c2_data(signed=True, limit=1))
        self.assertEqual(section.model.current_section(section.model.zero_image), ((0,), 0, 0))
        with self.assertRaises(DegreeSixResourceLimit):
            section.model.current_section(((0,), 1, 0))

    def test_malformed_linear_coordinates_are_rejected(self):
        bar = FiniteBar([[0, 1], [1, 0]], [0], [0])
        for change, message in ((dict(extra=1), "unknown"),
                                (dict(deltaB=[2]), "packed"),
                                (dict(inverse=[[2]]), "integer inverse"),
                                (dict(maxSectionCandidates=True), "positive integer")):
            data = c2_data()
            data.update(change)
            with self.subTest(change=change), self.assertRaisesRegex(ValueError, message):
                configure_degree_six(bar, data)


if __name__ == "__main__":
    unittest.main()
