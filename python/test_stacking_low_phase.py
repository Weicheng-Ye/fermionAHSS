"""Exhaustive local certificates for the bounded production stacking table.

Associativity and exchange are checked modulo an explicit universal binary
edge coboundary. Taking the signed differential then gives an integral
boundary witness for the corresponding D-coordinate defect.
"""
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import unittest


PAYLOAD = json.loads((Path(__file__).resolve().parents[1]
    / "data/stacking-low-phase.json").read_text())
PHASES = [Fraction(*pair) for pair in PAYLOAD["values"]]


def phase(b, c, bp, cp, s):
    index = 0
    for bit in (b, bp, *c, *cp, *s):
        index = 2 * index + bit
    return PHASES[index]


def composed(b, c, bp, cp, s):
    return ((b + bp) % 2,
            tuple((x + y + z*b*bp) % 2 for x, y, z in zip(c, cp, s)))


def binary_primitive(values, rank):
    """Solve delta lambda=values on B(F2^rank), returning all edge values."""
    pivots = {}
    for (g, h), value in values.items():
        if (2 * value).denominator != 1:
            raise AssertionError("the phase defect is not a binary half lift")
        mask = (1 << g) ^ (1 << h) ^ (1 << (g ^ h))
        rhs = int(2 * value) % 2
        while mask:
            pivot = mask.bit_length() - 1
            if pivot in pivots:
                mask ^= pivots[pivot][0]
                rhs ^= pivots[pivot][1]
            else:
                pivots[pivot] = (mask, rhs)
                break
        if not mask and rhs:
            raise AssertionError("phase defect has no universal binary primitive")
    solution = [0] * (1 << rank)
    for pivot in sorted(pivots):
        mask, rhs = pivots[pivot]
        solution[pivot] = (rhs + sum(solution[j]
            for j in range(pivot) if mask & (1 << j))) % 2
    for (g, h), value in values.items():
        assert (solution[g] + solution[h] + solution[g ^ h]) % 2 == int(2 * value) % 2
    return solution


def bits(value, rank):
    return tuple((value >> j) & 1 for j in range(rank))


class LowStackingPhaseTests(unittest.TestCase):
    def test_payload_integrity_and_normalization(self):
        encoded = json.dumps(PAYLOAD["values"], separators=(",", ":")).encode()
        self.assertEqual(hashlib.sha256(encoded).hexdigest(),
            "4ae828fdfa4c5ebcc93b3919c4f7ee9615d356f432dda908256130f1b456b11b")
        self.assertEqual(PAYLOAD["values_sha256"], hashlib.sha256(encoded).hexdigest())
        self.assertEqual(len(PHASES), 256)
        self.assertEqual(len(PAYLOAD["source_sha256"]), 9)
        for b, x, y, s, t in itertools.product(range(2), repeat=5):
            self.assertEqual(phase(b, (x,y), 0, (0,0), (s,t)), 0)
            self.assertEqual(phase(0, (0,0), b, (x,y), (s,t)), 0)

    def test_all_signed_integral_carries(self):
        # Every tetrahedron's three one-cocycles are determined by adjacent edges.
        for b, bp in itertools.product(range(2), repeat=2):
            for edge_bits in itertools.product(range(2), repeat=9):
                c, cp, s = edge_bits[:3], edge_bits[3:6], edge_bits[6:]
                def faces(v):
                    return ((v[1],v[2]), ((v[0]+v[1])%2,v[2]),
                            (v[0],(v[1]+v[2])%2), (v[0],v[1]))
                gamma = sum(((-1)**s[0] if j == 0 else (-1)**j)
                    * phase(b, faces(c)[j], bp, faces(cp)[j], faces(s)[j])
                    for j in range(4))
                self.assertEqual(gamma.denominator, 1)

    def test_associativity_has_universal_boundary_witness(self):
        # All 8 constant-B triples and 256 pairs of universal edges.
        for b, bp, bpp in itertools.product(range(2), repeat=3):
            values = {}
            for g, h in itertools.product(range(16), repeat=2):
                x, y = bits(g, 4), bits(h, 4)
                c, cp, cpp, s = ((x[j], y[j]) for j in range(4))
                b01, c01 = composed(b,c,bp,cp,s)
                b12, c12 = composed(bp,cp,bpp,cpp,s)
                values[g,h] = (phase(b,c,bp,cp,s)+phase(b01,c01,bpp,cpp,s)
                    -phase(bp,cp,bpp,cpp,s)-phase(b,c,b12,c12,s))
            self.assertEqual(binary_primitive(values, 4)[0], 0)

    def test_exchange_has_universal_boundary_witness(self):
        for b, bp in itertools.product(range(2), repeat=2):
            values = {}
            for g, h in itertools.product(range(8), repeat=2):
                x, y = bits(g, 3), bits(h, 3)
                c, cp, s = ((x[j], y[j]) for j in range(3))
                values[g,h] = phase(b,c,bp,cp,s)-phase(bp,cp,b,c,s)
            self.assertEqual(binary_primitive(values, 3)[0], 0)


if __name__ == "__main__":
    unittest.main()
