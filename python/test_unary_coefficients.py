"""Validation of the compiled degree-six unary constants.

Copyright (c) 2026 koAHSS contributors; MIT license.
"""
import copy
from fractions import Fraction
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import unary_gamma6 as unary


class UnaryCoefficientTests(unittest.TestCase):
    def setUp(self):
        self.table = json.loads(unary.TABLE.read_text())

    def test_bundled_table_matches_formula_sources_and_exact_values(self):
        values = unary.check_table(self.table)
        self.assertEqual(len(values), 22)
        for i, cell in enumerate(unary.unary_cells()):
            self.assertEqual(values[cell], Fraction(self.table['raw'][str(i)]) % 1)

    def test_validator_imports_without_a_worker(self):
        subprocess.run([sys.executable, '-c',
            'import unary_gamma6; assert unary_gamma6.table_status()["status"] == "valid"'],
            cwd=Path(__file__).parent, check=True, capture_output=True, text=True)

    def test_invalid_tables_are_refused(self):
        mutations = [lambda d: d.update(schema=2),
                     lambda d: d.update(sourceHash='stale'),
                     lambda d: d.update(basis='another basis'),
                     lambda d: d['q_mod1'].pop(),
                     lambda d: d['q_mod1'].__setitem__(0, None),
                     lambda d: d['q_mod1'].__setitem__(0, '1/5'),
                     lambda d: d['q_mod1'].__setitem__(0, '1'),
                     lambda d: d['raw'].__setitem__('0', '1/2'),
                     lambda d: d['cells'].reverse(),
                     lambda d: d['pairCoefficients'].__setitem__(0, '0')]
        for mutate in mutations:
            data = copy.deepcopy(self.table)
            mutate(data)
            with self.assertRaises(unary.UnaryTableUnavailable):
                unary.check_table(data)

    def test_missing_and_malformed_files_are_refused(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'table.json'
            self.assertEqual(unary.table_status(path)['status'], 'invalid')
            path.write_text('{')
            self.assertEqual(unary.table_status(path)['status'], 'invalid')
            path.write_text(json.dumps(self.table))
            self.assertEqual(unary.table_status(path)['status'], 'valid')

    @unittest.skipUnless(os.environ.get('FERMIONAHSS_SLOW_LIGHT_TESTS') == '1',
                         'opt-in reconstruction from the pair contractor')
    def test_recompile_two_cells_from_the_pair_contractor(self):
        from generate_unary_coefficients import compile_cell
        for i in (14, 18):
            record = compile_cell(i, cache=os.environ.get('FERMIONAHSS_CACHE_DIR') or None)
            self.assertEqual(Fraction(record['value']), Fraction(self.table['raw'][str(i)]))


if __name__ == '__main__':
    unittest.main()
