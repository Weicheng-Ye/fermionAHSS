"""Checks of the identities in doc/light-transport-reduction.md.

These check the formulas and their runtime reductions against the original
marked cochains. Finite checks supplement the proofs in the note.
Copyright (c) 2026 koAHSS contributors; MIT license.
"""
from fractions import Fraction
from itertools import combinations, product
import json
from math import prod
from pathlib import Path
import random
import unittest
from unittest.mock import patch

from extension_light import LightEvaluator, p
from extension_transfer import TransferredModel
from test_extension_transfer import lazy_c2_models


def c4_comparison(reorient=False):
    fixture = json.loads((Path(__file__).parent / 'test_fixtures' /
                          'light_c4.json').read_text())
    records = {(r['kind'], r['degree'], tuple(r['vertices'])
                if isinstance(r['vertices'], list) else r['vertices']): r['answer']
               for r in fixture['transport']}
    if reorient:
        # Change each odd-degree native basis by -1. This is the same
        # integral retraction, with nonzero carries for canonical 0/1 input
        # vectors; f and g change by inverse basis maps and H is unchanged.
        for name in ('ordinary', 'signed'):
            fixture['setup'][name] = [[[-x for x in row] for row in matrix]
                                       for matrix in fixture['setup'][name]]
        for (kind, degree, _), answer in records.items():
            if degree % 2 and kind in ('f', 'g'):
                for term in answer['terms']:
                    if kind == 'f':
                        term[1], term[2] = -term[1], -term[2]
                    else:
                        term[0], term[1] = -term[0], -term[1]
    calls = []

    def transport(kind, degree, vertices):
        calls.append((kind, degree, vertices))
        return records[kind, degree, vertices]

    return TransferredModel(fixture['setup'], transport), calls


def homotopy(model, cochain):
    return p.Cochain(cochain.degree - 1, lambda face: sum(
        t[0] * cochain(tuple(t[2])) for t in model.terms('h', face)) % 2)


def simplex_fields(seed, vertices=7):
    rng = random.Random(seed)

    def arbitrary(degree):
        values = {f: rng.randrange(2)
                  for f in combinations(range(vertices), degree + 1)}
        return p.Cochain(degree, lambda f: values.get(tuple(f), 0))

    def closed(degree):
        return p.binary(p.differential(arbitrary(degree - 1)))

    return arbitrary, closed


def extension_exponent(orders, relation):
    """Exponent of (H + Z a)/(2a-t), by its finite normal-form group law."""
    exponent = 1
    for a in range(2):
        for h in product(*(range(n) for n in orders)):
            b, value = 0, (0,) * len(orders)
            for order in range(1, 2 * prod(orders) + 1):
                carry, b = divmod(b + a, 2)
                value = tuple((x + y + carry * t) % n
                              for x, y, t, n in zip(value, h, relation, orders))
                if b == 0 and not any(value):
                    exponent = max(exponent, order)
                    break
            else:
                raise AssertionError('normal-form group element has no order')
    return exponent


class LightTransportReductionTests(unittest.TestCase):
    def test_runtime_b_carry_preserves_the_original_primitive(self):
        for reorient in (False, True):
            model, calls = c4_comparison(reorient)
            light = LightEvaluator(model)
            for n in (1, 2):
                b = light.lift(n, [1])
                C = light.lift(n + 1, [1])
                P = p.binary(p.differential(light.lift(n, [1])))
                old = light._prim_normalized(p.binary(p.hD(b, b, model.s) + C + P), [1])
                new = light.b_primitive(n + 2, [1], P, [1], C)
                for face in product(range(4), repeat=n + 1):
                    self.assertEqual(new(face), old(face), (n, face, reorient))

    def test_b_gauge_projection_uses_the_lower_degree_carry(self):
        for k in (2, 3, 4):
            model, calls = c4_comparison(reorient=True)
            light = LightEvaluator(model)
            data = dict(b=[1], Cref=[[1]], gauge={}, want=['xC', 'target'])
            new = light.task_gauge(k, data)
            self.assertFalse(any(kind == 'g' and n >= k - 1 for kind, n, _ in calls))
            light._reduce_transport = False
            self.assertEqual(new, light.task_gauge(k, data))

    def test_strict_pairings_do_not_request_comparison_chains(self):
        model, calls = c4_comparison()
        model.strict_identities = True
        light = LightEvaluator(model)
        with patch.object(model, 'chain', side_effect=AssertionError('unneeded g')):
            self.assertEqual(light.bin(light.lift(2, [1])), [1])
            self.assertEqual(light.integer(light.lift(2, [3], True)), [3])
            h = light._homotopy(p.Cochain(3, lambda f: 1))
            self.assertEqual(light.bin(h), [0])
            self.assertEqual(light.bin(light.lift(2, [1]) + h), [1])

    def test_primary_d_gauge_keeps_the_exact_original_homotopy(self):
        model, _ = c4_comparison(reorient=True)
        light = LightEvaluator(model)
        b, y = light.lift(2, [1]), light.lift(1, [1])
        P = p.QD(y, model.s, model.omega)
        C = light.lift(3, [1])
        old = light._prim_normalized(p.binary(p.hD(b, b, model.s) + C + P), [1])
        # A deliberately different valid comparison K=HP+delta L+Lambda z
        # exercises both Pi K and delta H K, rather than a zero comparison.
        L = p.Cochain(1, lambda f: 0 if model.normalize(f) is None else
                      int(model.normalize(f)[1] in (1, 2)))
        K = p.binary(model.homotopy(P, False) + p.differential(L) + light.lift(2, [1]))
        model.primary_comparison = True
        with patch.object(light, 'primary_comparison', return_value=(K, None)):
            new = light.b_primitive(4, [1], P, [1], C, primary_gauge=[1])
            for face in product(range(4), repeat=3):
                self.assertEqual(new(face), old(face), face)

    def test_runtime_pure_c_curvature_and_marked_square_class(self):
        for k, sign, omega in product(range(1, 7), range(2), range(2)):
            model, _, calls = lazy_c2_models(sign, omega, top=8)
            light = LightEvaluator(model)
            new = light.task_c_mark(k, {'c': [1]})
            self.assertTrue(all(n <= k + 1 for kind, n, _ in calls if kind == 'g'))
            light._reduce_transport = False
            old = light.task_c_mark(k, {'c': [1]})
            self.assertEqual(new['J'], old['J'], (k, sign, omega))
            difference = new['gamma'][0] - old['gamma'][0]
            incoming = model.matrices[True][k][0][0]
            self.assertEqual(difference % incoming if incoming else difference, 0,
                             (k, sign, omega))

    def test_bockstein_carry_and_old_primitive_on_nonidentity_comparison(self):
        model, calls = c4_comparison(reorient=True)
        light = LightEvaluator(model)
        nonzero_carries = 0
        for n in (1, 2):
            lifted = model.lift(n, [1], True)
            b = p.binary(lifted)
            carry = p.divide(lifted - b, 2, 'Bockstein comparison carry')
            binary_carry = p.binary(carry)
            q_r = [x // 2 for x in model.coboundary(n, [1], True)]
            lifted_q = model.lift(n + 1, q_r, True)
            bockstein = p.divide(p.ds(b, model.s), 2, 'binary Bockstein')
            x = p.hD(b, b, model.s)
            reference = p.binary(lifted_q)
            source = p.binary(x + reference)

            # With this reference r_new=0, so the new primitive is the carry.
            start = len(calls)
            for tail in product(range(4), repeat=n):
                face = (0,) + tail
                nonzero_carries += int(carry(face) != 0)
            self.assertFalse(any(kind == 'h' for kind, *_ in calls[start:]))

            old_vector = light.bin(binary_carry)
            old = light._prim_normalized(source, old_vector)
            correction = p.binary(p.differential(homotopy(model, binary_carry)))
            for tail in product(range(4), repeat=n + 1):
                face = (0,) + tail
                self.assertEqual(bockstein(face),
                                 (lifted_q - p.ds(carry, model.s))(face))
                self.assertEqual(x(face), bockstein(face) % 2)
                self.assertEqual(p.differential(binary_carry)(face) % 2,
                                 source(face))
                self.assertEqual(old(face[:-1]),
                                 (binary_carry(face[:-1]) + correction(face[:-1])) % 2)
        self.assertGreater(nonzero_carries, 0)

    def test_integral_relation_primitive_is_a_lift(self):
        model, calls = c4_comparison()
        light = LightEvaluator(model)
        for n in (0, 2):
            u = [3]
            a_r = [x // 2 for x in model.coboundary(n, u, True)]
            a = model.lift(n + 1, a_r, True)
            new = model.lift(n, u, True)
            old = light._prim_normalized(p.scale(a, 2), u, True)
            start = len(calls)
            for tail in product(range(4), repeat=n):
                new((0,) + tail)
            self.assertFalse(any(kind == 'h' for kind, *_ in calls[start:]))
            for tail in product(range(4), repeat=n + 1):
                face = (0,) + tail
                self.assertEqual(p.ds(new, model.s)(face), 2 * a(face))
                self.assertEqual(new(face[:-1]), old(face[:-1]))

    def test_primary_comparison_with_nontrivial_coherent_cup_changes(self):
        # Work on a simplex with Lambda=id but genuinely different higher
        # cups. Arbitrary bilinear homotopies define the second system by
        # equation (8); this exercises all four terms of K_D, including the
        # comparison of the nested sign times Sq1 operation.
        for seed in range(8):
            arbitrary, closed = simplex_fields(5300 + seed)
            r, s, w = arbitrary(0), closed(1), closed(2)

            def K(i, a, b):
                if i < 0:
                    return p.zero(a.degree + b.degree - i - 1)
                return p.cup(r, p.cup(a, b, i + 1))

            def native_cup(i, a, b):
                if i < 0:
                    return p.zero(a.degree + b.degree - i)
                return p.binary(p.cup(a, b, i) + p.differential(K(i, a, b))
                                + K(i, p.differential(a), b)
                                + K(i, a, p.differential(b))
                                + K(i - 1, a, b) + K(i - 1, b, a))

            for n in (1, 2, 3):
                b = closed(n)
                sq1 = native_cup(n - 1, b, b)
                native_source = p.binary(native_cup(n - 2, b, b)
                                         + native_cup(0, w, b)
                                         + native_cup(0, s, sq1))
                comparison = p.binary(K(n - 2, b, b) + K(0, w, b)
                                      + K(0, s, sq1) + p.cup(s, K(n - 1, b, b)))
                residual = p.binary(p.differential(comparison)
                                    + p.QD(b, s, w) + native_source)
                for face in combinations(range(7), n + 3):
                    self.assertEqual(residual(face), 0, (seed, n, face))

    def test_degree_four_canonical_a_and_leading_square(self):
        for seed in range(12):
            _, closed = simplex_fields(8700 + seed)
            s, w = closed(1), closed(2)
            unit = p.Cochain(0, lambda f: 1)
            u = p.scale(unit, -1)
            for face in combinations(range(7), 2):
                self.assertEqual(p.ds(u, s)(face), 2 * s(face))
            for face in combinations(range(7), 3):
                self.assertEqual(p.ds(s, s)(face), 0)
                self.assertEqual(p.hD(s, s, s)(face), 0)
                self.assertEqual(p.QD(unit, s, w)(face), w(face))

    def test_c_coboundary_change_has_an_integral_d_carry(self):
        for seed, n in product(range(8), (2, 3)):
            arbitrary, closed = simplex_fields(6400 + seed, vertices=8)
            s, w = closed(1), closed(2)
            c, eta = arbitrary(n), arbitrary(n - 1)
            delta_eta = p.binary(p.differential(eta))
            changed = p.binary(c + delta_eta)
            v = p.binary(p.E(eta, w) + p.polarization(c, delta_eta, 2))
            potential_difference = p.half(p.E(changed, w) - p.E(c, w))
            carry = p.integral(potential_difference - p.half(p.ds(v, s)),
                               'C coboundary comparison carry')
            for face in combinations(range(8), n + 3):
                self.assertEqual(carry(face).denominator, 1)
            residual = p.ds(carry, s) - p.ds(potential_difference, s)
            for face in combinations(range(8), n + 4):
                self.assertEqual(residual(face), 0, (seed, n, face))

    def test_pure_c_marked_class_without_successor_comparison_chains(self):
        for k, sign, omega in product((3, 4), range(2), range(2)):
            with self.subTest(k=k, sign=sign, omega=omega):
                model, _, calls = lazy_c2_models(sign, omega, top=8)
                light = LightEvaluator(model)
                c = light.lift(k - 1, [1])
                e = light.integer(p.E(c, model.omega))
                delta_e = model.coboundary(k + 1, e, True)
                self.assertTrue(all(x % 2 == 0 for x in delta_e))
                curvature = [x // 2 for x in delta_e]
                y = [(x // 2) % 2 for x in
                     model.coboundary(k - 1, [1], False)]
                delta_y = model.coboundary(k, y, True)
                self.assertTrue(all(x % 2 == 0 for x in delta_y))
                gamma = [x + y // 2 for x, y in zip(e, delta_y)]
                self.assertTrue(all(n <= k + 1 for kind, n, _ in calls
                                    if kind == 'g'))

                old = light.task_c_mark(k, {'c': [1]})
                self.assertEqual(curvature, old['J'])
                difference = old['gamma'][0] - gamma[0]
                incoming = model.matrices[True][k][0][0]
                if incoming:
                    self.assertEqual(difference % incoming, 0)
                else:
                    self.assertEqual(difference, 0)
                # Equal curvature gives the same native D completion;
                # a coboundary difference then gives the same marked class.

    def test_absorption_is_a_filtered_shear_of_a_nonsplit_lower_group(self):
        # H=Z/4{c}+Z/2{b}, D=<2c>. The extension below C is nonsplit.
        elements = list(product(range(4), range(2)))

        def add(x, y):
            return ((x[0] + y[0]) % 4, (x[1] + y[1]) % 2)

        for lc, lb in ((1, 0), (0, 1), (1, 1)):
            character = lambda h: (lc * h[0] + lb * h[1]) % 2
            for d in ((0, 0), (2, 0)):
                def shear(h):
                    return add(h, d if character(h) else (0, 0))
                for h in elements:
                    self.assertEqual(shear(shear(h)), h)
                    self.assertEqual((shear(h)[0] % 2, shear(h)[1]),
                                     (h[0] % 2, h[1]))
                    for j in elements:
                        self.assertEqual(shear(add(h, j)), add(shear(h), shear(j)))
                    if character(h):
                        self.assertEqual(shear(add(h, d)), h)

    def test_absorption_hypotheses_cannot_be_omitted(self):
        self.assertEqual(extension_exponent((2,), (0,)), 2)
        self.assertEqual(extension_exponent((2,), (1,)), 4)
        self.assertEqual(extension_exponent((2, 4), (1, 0)), 4)
        self.assertEqual(extension_exponent((2, 4), (1, 1)), 8)

    def test_division_retains_the_residue_only_at_the_required_precision(self):
        for q in (2, 4, 8, 16):
            for m in (2, 4, 8):
                for z in range(-35, 36):
                    numerator = q * z
                    self.assertEqual((numerator % (m * q)) // q, z % m)
        self.assertNotEqual((4 % 4) // 4, 1 % 2)

    def test_pairing_a_coboundary_with_a_target_character(self):
        model, _ = c4_comparison()
        light = LightEvaluator(model)
        for n in (1, 2):
            def value(face):
                normalized = model.normalize(face)
                if normalized is None:
                    return Fraction(0)
                return Fraction(sum((i + 2) * x for i, x in enumerate(normalized)), 4)
            phase = p.Cochain(n, value)
            for ell in (-3, 0, 2):
                projected = light.rational(p.ds(phase, model.s))[0] * ell
                boundary = model.matrices[True][n][0][0] * ell
                paired = sum(boundary * t[1] * phase(tuple(t[2]))
                             for t in model.chain(n, 0))
                self.assertEqual(projected, paired)


if __name__ == '__main__':
    unittest.main()
