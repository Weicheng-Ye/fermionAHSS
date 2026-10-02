"""The universal values shipped in data/universal-values.json are current."""
import json
from fractions import Fraction
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import extension_transfer  # noqa: E402,F401  (installs the worker's exact policy)
import extension_acceleration as acceleration  # noqa: E402
import universal_values as universal  # noqa: E402
import universal_sources  # noqa: E402


class BundledUniversalValues(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads(universal.BUNDLED.read_text())

    def test_provenance_matches_the_formula_sources(self):
        # The worker ignores a file computed with other sources.
        self.assertEqual(self.data['schema'], 1)
        self.assertEqual(self.data['provenance'], universal.source_provenance(),
                         'regenerate with python3 python/generate_universal_values.py')

    def test_entries_are_sorted_and_unique(self):
        self.assertEqual(sorted(self.data['tables']),
                         sorted(set(acceleration._tables) - universal_sources.UNBUNDLED))
        for entries in self.data['tables'].values():
            keys = [json.dumps(key, sort_keys=True) for key, _ in entries]
            self.assertEqual(keys, sorted(set(keys)))

    def test_sampled_values_equal_a_fresh_evaluation(self):
        classes = universal._key_classes()
        for name, entries in sorted(self.data['tables'].items()):
            table = acceleration._tables[name]
            values = [(key, universal._decode(value, classes)) for key, value in entries]
            nonzero = [(key, value) for key, value in values if value != 0]
            for key, value in nonzero[::max(1, len(nonzero) // 3)][:3]:
                self.assertEqual(table.recompute(universal._decode(key, classes)), value)


class DeferredUniversalValues(unittest.TestCase):
    def test_unused_tables_are_not_decoded_and_first_value_wins(self):
        first = universal.UniversalTable('first', lambda key: Fraction(key, 4))
        unused = universal.UniversalTable('unused', lambda key: 1 / 0)
        store = universal.UniversalStore(None, 'test', {'first': first, 'unused': unused})
        first.values[2] = Fraction(1, 2)
        store._insert({'first': [(1, {'f': '1/4'}), (2, {'f': '9/4'})],
                       'unused': [(1, {'d': 'unknown.Class', 'v': []})]})
        store._insert({'first': [(1, {'f': '7/4'}), (3, {'f': '3/4'})]})
        with patch.object(universal, '_decode', wraps=universal._decode) as decode:
            self.assertEqual(first(1), Fraction(1, 4))
            self.assertEqual(first.values, {1: Fraction(1, 4), 2: Fraction(1, 2),
                                            3: Fraction(3, 4)})
            count = decode.call_count
            self.assertEqual(first(3), Fraction(3, 4))
            self.assertEqual(decode.call_count, count)
        self.assertEqual(first.pending, {})
        self.assertTrue(unused._stored_entries)
        self.assertEqual(first(4), 1)
        self.assertEqual(first.pending, {4: 1})


if __name__ == '__main__':
    unittest.main()
