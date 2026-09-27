"""The universal values shipped in data/universal-values.json are current."""
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import extension_transfer  # noqa: E402,F401  (installs the worker's exact policy)
import extension_acceleration as acceleration  # noqa: E402


class BundledUniversalValues(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads(acceleration.BUNDLED.read_text())

    def test_provenance_matches_the_formula_sources(self):
        # The worker ignores a file computed with other sources.
        self.assertEqual(self.data['schema'], 1)
        self.assertEqual(self.data['provenance'], acceleration.source_provenance(),
                         'regenerate with python3 python/generate_universal_values.py')

    def test_entries_are_sorted_and_unique(self):
        self.assertEqual(sorted(self.data['tables']), sorted(acceleration._tables))
        for entries in self.data['tables'].values():
            keys = [json.dumps(key, sort_keys=True) for key, _ in entries]
            self.assertEqual(keys, sorted(set(keys)))

    def test_sampled_values_equal_a_fresh_evaluation(self):
        classes = acceleration._key_classes()
        for name, entries in sorted(self.data['tables'].items()):
            function = acceleration._tables[name].function
            values = [(key, acceleration._decode(value, classes)) for key, value in entries]
            nonzero = [(key, value) for key, value in values if value != 0]
            for key, value in nonzero[::len(nonzero) // 3][:3]:
                self.assertEqual(function(acceleration._decode(key, classes)), value)


if __name__ == '__main__':
    unittest.main()
