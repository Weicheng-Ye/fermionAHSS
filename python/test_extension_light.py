"""Exact transport and failure contracts of the light residue evaluator.

Copyright (c) 2026 koAHSS contributors; MIT license.
"""
from fractions import Fraction
import itertools
import json
import os
from pathlib import Path
import unittest
from unittest.mock import patch

from extension_light import LightEvaluator, p
from extension_transfer import TransferredModel, TransferResourceLimit
from test_extension_transfer import c2_models, lazy_c2_models


class LightTransportTests(unittest.TestCase):
    def test_primitive_and_chain_map_on_nonidentity_c4_comparison(self):
        fixture = json.loads((Path(__file__).parent / 'test_fixtures' / 'light_c4.json').read_text())
        records = {(r['kind'], r['degree'], tuple(r['vertices']) if isinstance(r['vertices'], list)
                    else r['vertices']): r['answer'] for r in fixture['transport']}
        calls = []
        def transport(kind, degree, vertices):
            calls.append((kind, degree, vertices))
            return records[kind, degree, vertices]
        model = TransferredModel(fixture['setup'], transport)
        light = LightEvaluator(model)
        def value(vertices):
            normalized = model.normalize(vertices)
            return 0 if normalized is None else sum((i + 2) * v for i, v in enumerate(normalized))
        for n in (1, 2):
            r = p.Cochain(n, value)
            source = p.ds(r, model.s)
            vector = model.project(r, True)
            self.assertEqual(model.coboundary(n, vector, True), light.integer(source))
            primitive = light.prim(source, vector, True)
            trusted = light._prim_normalized(source, vector, True)
            for tail in itertools.product(range(4), repeat=n + 1):
                simplex = (0,) + tail
                if model.normalize(simplex) is not None:
                    self.assertEqual(p.ds(primitive, model.s)(simplex), source(simplex))
                    self.assertEqual(trusted(simplex[:-1]), primitive(simplex[:-1]))
                    self.assertEqual(p.ds(trusted, model.s)(simplex), source(simplex))
            half = p.scale(r, Fraction(1, 2))
            self.assertEqual(light.rational_ds(half), light.rational(p.ds(half, model.s)))
            self.assertEqual(light._rational_ds_normalized(half), light.rational_ds(half))
        self.assertTrue(any(k == 'h' and records[k, n, v]['terms'] for k, n, v in calls))

    def test_eager_and_lazy_rational_pairing(self):
        for sign in (0, 1):
            for make in (c2_models, lazy_c2_models):
                model, *_ = make(sign, 0, top=8)
                light = LightEvaluator(model)
                c = p.scale(model.lift(3, [5], True), Fraction(1, 3))
                self.assertEqual(light.rational(c), [Fraction(5, 3)])
                self.assertEqual(light.rational_ds(c),
                                 light.rational(p.ds(c, model.s)))
                answer = model.calculate(dict(operation='light', degree=6,
                    task='atom_curvature', data=dict(b=[0], c=[0])))
                self.assertEqual(answer['result'], dict(J=[0]))

    def test_primitive_equation_with_twisted_coefficients(self):
        for sign, n in ((0, 1), (1, 2)):
            model, bar, _ = c2_models(sign)
            light = LightEvaluator(model)
            r = model.lift(n, [3], True)
            source = p.ds(r, model.s)
            primitive = light.prim(source, [3], True)
            for simplex in bar.simplices(n + 1):
                self.assertEqual(p.ds(primitive, model.s)(simplex), source(simplex))
            binary = light.prim(p.binary(source), [1])
            self.assertEqual(light.bin(p.differential(binary)), light.bin(source))

    def test_chain_map_rejects_unnormalized_cochain(self):
        model, *_ = c2_models()
        light = LightEvaluator(model)
        # The face (0,0) of the normalized simplex (0,1,0) is omitted
        # by the normalized boundary. The constant cochain does not vanish.
        with self.assertRaisesRegex(ArithmeticError, 'not normalized'):
            light.rational_ds(p.Cochain(1, lambda f: Fraction(1, 2)))

    def test_primitive_rejects_unnormalized_source(self):
        model, *_ = c2_models()
        primitive = LightEvaluator(model).prim(p.Cochain(2, lambda f: 1), [0])
        with self.assertRaisesRegex(ArithmeticError, 'not normalized'):
            primitive((0, 1))

    def test_production_tasks_do_not_repeat_normalization_checks(self):
        model, *_ = c2_models(top=8)
        with patch.dict(os.environ, FERMIONAHSS_LIGHT_NORMALIZATION_CHECKS='0'):
            light = LightEvaluator(model)
        model.light = light
        with patch.object(light, 'check_normalized', side_effect=AssertionError('repeated check')), \
                patch.object(light, 'normalized_degeneracies', side_effect=AssertionError('repeated check')):
            # The production request constructs a primitive and pairs an
            # integral coboundary of a rational potential in degrees 4--6.
            for k in (4, 5, 6):
                answer = model.calculate(dict(operation='light', degree=k,
                    task='atom_curvature', data=dict(b=[0], c=[0])))
                self.assertEqual(answer['result'], dict(J=[0]))
            primitive = light._prim_normalized(model.lift(2, [2], True), [0], True)
            self.assertEqual(primitive((0, 1)), 0)

    def test_normalized_rational_coboundary_needs_no_successor_chains(self):
        for sign, expected in ((0, Fraction(10, 3)), (1, Fraction(0))):
            model, _, calls = lazy_c2_models(sign, top=8)
            with patch.dict(os.environ, FERMIONAHSS_LIGHT_NORMALIZATION_CHECKS='0'):
                light = LightEvaluator(model)
            c = p.scale(model.lift(3, [5], True), Fraction(1, 3))
            self.assertEqual(light._rational_ds_normalized(c), [expected])
            self.assertEqual([call for call in calls if call[0] == 'g'], [('g', 3, 0)])

    def test_audit_mode_checks_trusted_paths_and_generic_helpers_stay_checked(self):
        for audit in ('0', '1'):
            model, *_ = c2_models()
            with patch.dict(os.environ, FERMIONAHSS_LIGHT_NORMALIZATION_CHECKS=audit):
                light = LightEvaluator(model)
            potential = p.Cochain(1, lambda f: Fraction(1, 2))
            source = p.Cochain(2, lambda f: 1)
            with self.assertRaisesRegex(ArithmeticError, 'not normalized'):
                light.rational_ds(potential)
            with self.assertRaisesRegex(ArithmeticError, 'not normalized'):
                light.prim(source, [0])((0, 1))
            if audit == '1':
                with self.assertRaisesRegex(ArithmeticError, 'not normalized'):
                    light._rational_ds_normalized(potential)
                with self.assertRaisesRegex(ArithmeticError, 'not normalized'):
                    light._prim_normalized(source, [0])((0, 1))

    def test_zero_primitive_needs_no_comparison_and_preserves_formula_zero(self):
        model, *_ = c2_models()
        light = LightEvaluator(model)
        with patch.object(model, 'terms', side_effect=AssertionError('unneeded comparison')):
            for signed in (False, True):
                primitive = light.prim(p.zero(3), [0], signed)
                self.assertEqual(primitive((0, 1, 0)), 0)
                # Higher formulas must see the same zero as native H(0),
                # so they do not build costly universal pair contractors.
                self.assertTrue(getattr(primitive, 'structural_zero', False))
                self.assertEqual(light.integer(p.ds(primitive, model.s)), [0])

    def test_primary_sources_vanish_on_degeneracies(self):
        for sign in (0, 1):
            for omega in (0, 1):
                model, *_ = c2_models(sign, omega)
                light = LightEvaluator(model)
                for n in (1, 2, 3):
                    source = p.QD(light.lift(n, [1]), model.s, model.omega)
                    light.normalized_degeneracies(source, tuple(i % 2 for i in range(n + 2)))

    def test_exact_and_localized_residues(self):
        self.assertEqual(LightEvaluator.exact([Fraction(-4), 3], 'test'), [-4, 3])
        self.assertEqual(LightEvaluator.modular([Fraction(2, 3)], 8, 2, 'test'), [6])
        with self.assertRaises(ArithmeticError):
            LightEvaluator.exact([Fraction(1, 2)], 'test')
        with self.assertRaises(ArithmeticError):
            LightEvaluator.modular([Fraction(1, 2)], 8, 2, 'test')

    def test_worker_failures_remain_in_band_and_stream_recovers(self):
        model, *_ = c2_models()
        light = LightEvaluator(model)
        request = dict(task='qd', degree=3, data=dict(b=[0]))
        for error, field in ((TransferResourceLimit('bounded'), 'refused'),
                             (MemoryError('bounded'), 'refused'),
                             (ArithmeticError('bad identity'), 'error')):
            with patch.dict(light.TASKS, qd=lambda *args, error=error: (_ for _ in ()).throw(error)):
                answer = light.calculate(request)
                self.assertIn(field, answer)
                if field == 'refused':
                    self.assertEqual(answer[field]['code'], 'resource-limit')
            self.assertEqual(light.calculate(request)['result'], dict(QD=[0]))


if __name__ == '__main__':
    unittest.main()
