"""Exact normalization certificates for the light-task formula constructors.

The ANF and reduced-power checks are symbolic. The section checks exercise
normalized contractors, complementing the proof in doc/light-normalization.md.
Copyright (c) 2026 koAHSS contributors; MIT license.
"""
import unittest

from extension_worker import api  # Installs the stacking import path.
import phase_eval as p
import mod3_power
import chain_models as cm
import low_phases as low
import source_primitive as sp
import r3_chain as r3
import h_tau_primitive as hp
import v1_pair_shared as pair1
import v2_pair_shared as pair2
import v3_pair_shared as pair3
import r1_pair_chain_signed as pair1_signed
import a0_high_gamma as a0



def chi_certificate():
    rows = []
    for n in range(8):
        faces, masks = p.chi_anf(n)
        residuals = []
        for d in range(n + 3):
            variables, substitution = {}, []
            for face in faces:
                image = tuple(i if i <= d else i - 1 for i in face)
                if len(set(image)) != len(image):
                    substitution.append(None)
                else:
                    if image not in variables:
                        variables[image] = len(variables)
                    substitution.append(1 << variables[image])
            result = set()
            for monomial in masks:
                image = 0
                while monomial:
                    bit = monomial & -monomial
                    monomial -= bit
                    factor = substitution[bit.bit_length() - 1]
                    if factor is None:
                        break
                    image |= factor
                else:
                    if image in result:
                        result.remove(image)
                    else:
                        result.add(image)
            residuals.append(len(result))
        assert not any(residuals), (n, residuals)
        rows.append(dict(degree=n, monomials=len(masks),
                         degeneracy_residual_monomials=residuals))
    return rows


def p1_certificate():
    rows = []
    for n in (2, 3):
        terms = mod3_power.reduced_power_terms(n)
        counts = [sum(not any(j in face and j + 1 in face for face in factors)
                      for factors, _ in terms) for j in range(n + 4)]
        assert not any(counts), (n, counts)
        rows.append(dict(input_degree=n, terms=len(terms),
                         terms_without_degenerate_factor=counts))
    return rows


def arbitrary(n, salt):
    def value(face):
        if any(a == b for a, b in zip(face, face[1:])):
            return 0
        return sum((i + salt) * (v + 1) + salt * i * i
                   for i, v in enumerate(face)) % 11 - 5
    return p.Cochain(n, value)


def background(salt):
    eps = lambda v: (v * v + v // 2 + salt * v) % 2
    s = p.Cochain(1, lambda f: (eps(f[0]) + eps(f[1])) % 2)
    w = p.binary(p.differential(arbitrary(1, salt + 10)))
    return s, w


def section_checks():
    counts = dict(unary_sections=0, binary_pair_sections=0,
                  rational_degree_one_pair_sections=0, theta_sections=0)
    for salt in (1, 2):
        s, w = background(salt)
        for n, pair in ((1, pair1), (2, pair2), (3, pair3)):
            A, Ap = arbitrary(n, salt + 20), arbitrary(n, salt + 30)
            for q in (n + 2, n + 3, n + 4):
                base = tuple(range(q))
                for d in range(q):
                    vertices = base[:d] + (base[d],) + base[d:]
                    unary = hp.to_pair(A, s, w, vertices)
                    if n == 1:
                        assert d in low._pair_degen(unary)
                        assert not low.group_product_AW(unary)
                        assert not low.group_product_H(unary)
                        assert not hp.group_homotopy(unary)
                    elif n == 2:
                        assert d in cm._degen_indices(unary)
                        assert not sp.Ftot(unary) and not sp.Htot(unary)
                    else:
                        assert d in r3.degens(unary)
                        assert not r3.Ftot(unary) and not r3.Htot(unary)
                    counts['unary_sections'] += 1
                    x = (pair.to_diag(A, Ap, s, vertices),
                         pair.to_omega(w, vertices))
                    assert d in pair.rc.degens(x)
                    assert not pair.rc.Ftot(x) and not pair.rc.Htot(x)
                    counts['binary_pair_sections'] += 1
                    if n == 1:
                        diag = x[0]
                        signed = (pair1_signed.Borel(diag.signs, diag.left, diag.right), x[1])
                        assert d in pair1_signed.degens(signed)
                        assert not pair1_signed.Ftot(signed) and not pair1_signed.Htot(signed)
                        counts['rational_degree_one_pair_sections'] += 1
        for m in (3, 4, 5):
            # Avoid constructor source-period calculations: only inspect the
            # section and its normalized chain operators.
            rule = object.__new__(a0.HigherA0Stacking)
            rule.m, rule.k = m, m + 2
            b = p.binary(p.differential(arbitrary(m - 1, salt + 40)))
            bp = p.binary(p.differential(arbitrary(m - 1, salt + 50)))
            q, base = m + 2, tuple(range(m + 2))
            for d in range(q):
                vertices = base[:d] + (base[d],) + base[d:]
                x = rule.section(b, bp, s, w, vertices)
                assert d in a0.degens(x)
                assert not a0.forward(x) and not a0.homotopy(x)
                counts['theta_sections'] += 1
    return counts


class LightNormalizationTests(unittest.TestCase):
    def test_chi_degeneracies_are_zero_polynomials(self):
        self.assertEqual(len(chi_certificate()), 8)

    def test_reduced_power_terms_have_degenerate_factors(self):
        self.assertEqual(len(p1_certificate()), 2)

    def test_universal_sections_preserve_degeneracies(self):
        self.assertEqual(section_checks(), dict(unary_sections=90,
            binary_pair_sections=90, rational_degree_one_pair_sections=24,
            theta_sections=36))


if __name__ == '__main__':
    unittest.main()
