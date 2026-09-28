"""The extension worker's evaluation policy returns the original values.

Bounded universal memos, the in-process theta tables and the recycled source
registry change where values are kept, never which values are computed.
Copyright (c) 2026 koAHSS contributors; MIT license.
"""
import os
import random
import unittest

os.environ.setdefault("FERMIONAHSS_CACHE_DIR", "")
import extension_transfer  # noqa: E402,F401  (installs the policy)
import extension_acceleration as acceleration  # noqa: E402
import phase_eval  # noqa: E402
import r3_chain  # noqa: E402
import r3_source  # noqa: E402
import chain_models  # noqa: E402


def random_pair(seed, n=8):
    """A universal simplex of the degree-eight source phase: nine vertices."""
    rng = random.Random(seed)
    rows = []
    for _ in range(n):
        matrix = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                matrix[i][j] = rng.randrange(-2, 3)
                matrix[j][i] = -matrix[i][j]
        rows.append(r3_chain.U2(rng.randrange(2), tuple(map(tuple, matrix))))
    omega = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            omega[i][j] = omega[j][i] = rng.randrange(2)
    return (r3_chain.Diag3(tuple(rows)),
            chain_models.Diag('c2', tuple((0, tuple(row)) for row in omega)))


class AccelerationPolicyTests(unittest.TestCase):
    def test_theta_values_are_tabled_per_pair_and_survive_recycling(self):
        table = acceleration._local_tables['r3_source.phi_value']
        original_evaluate_r = r3_source.evaluate_r.__wrapped__
        pairs = [random_pair(seed) for seed in (1, 2, 3)]
        chain = {pairs[0]: 2, pairs[1]: -1, pairs[2]: 1}
        expected = original_evaluate_r(r3_source.phi(), chain)
        self.assertEqual(r3_source.evaluate_r(r3_source.phi(), chain), expected)
        self.assertTrue(all(table.compact(pair) in table.values for pair in pairs))
        self.assertEqual(len({table.compact(pair) for pair in pairs}), 3)
        registered = len(phase_eval.SOURCES)
        self.assertGreater(registered, 0)
        # Recycling empties the registry and every memo, and the same values
        # come back from fresh registrations.
        saved = acceleration.SOURCE_RECYCLE_THRESHOLD
        try:
            acceleration.SOURCE_RECYCLE_THRESHOLD = 0
            self.assertTrue(acceleration.recycle_sources())
        finally:
            acceleration.SOURCE_RECYCLE_THRESHOLD = saved
        self.assertEqual(len(phase_eval.SOURCES), 0)
        self.assertEqual(r3_source.evaluate_r(r3_source.phi(), chain), expected)
        self.assertEqual(original_evaluate_r(r3_source.phi(), chain), expected)
        self.assertEqual(table(pairs[0]) * 2 + table(pairs[1]) * -1 + table(pairs[2]), expected)
        self.assertFalse(acceleration.recycle_sources())

    def test_universal_evaluations_construct_bounded_memos(self):
        with acceleration.universal():
            bounded = phase_eval.Cochain(1, lambda t: 1)
        plain = phase_eval.Cochain(1, lambda t: 1)
        self.assertEqual(bounded.evaluate.cache_parameters()['maxsize'],
                         acceleration.UNIVERSAL_MEMO_LIMIT)
        # Ordinary cochains keep the unbounded memo table of cochain_tools,
        # the C wrapper without cache_parameters.
        self.assertFalse(hasattr(plain.evaluate, 'cache_parameters'))
        self.assertEqual(acceleration._universal_depth, 0)

    def test_theta_pair_source_values_are_shared_between_instances(self):
        from theta_pair_phase import ThetaPairPhase
        import a0_high_gamma
        table = acceleration._local_tables['a0_high_gamma.HigherA0Stacking.evaluate']
        first, second = ThetaPairPhase(4), ThetaPairPhase(4)
        # The source has degree seven: a degree-seven tensor simplex of two
        # binary 4-cocycles, a sign diagonal and a binary 2-cocycle.
        q, rng = 7, random.Random(11)

        def cocycle(m):
            lower = {face: rng.randrange(2) for face in a0_high_gamma.faces(q, m - 1)}
            return a0_high_gamma.KB(m, q, tuple(
                sum(lower[face[:i] + face[i + 1:]] for i in range(m + 1)) % 2
                for face in a0_high_gamma.faces(q, m)))
        b, zero = cocycle(4), a0_high_gamma.KB(4, q, (0,) * len(a0_high_gamma.faces(q, 4)))
        s = chain_models.Diag('zsign', tuple((rng.randrange(2), (0,) * q) for _ in range(q)))
        w = cocycle(2)
        simplex = ((b, zero), (s, w))
        before = len(table.values)
        value = first.evaluate(simplex)
        self.assertEqual(second.evaluate(simplex), value)
        self.assertEqual(len(table.values), before + 1)
        self.assertEqual(value, acceleration._local_tables[
            'a0_high_gamma.HigherA0Stacking.evaluate'].function(('ThetaPairPhase', 4, simplex)))


if __name__ == "__main__":
    unittest.main()
