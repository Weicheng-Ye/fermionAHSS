"""Exact cochains on ordered simplices and normalized interval-cut cups.

The interval-cut engine is adapted from fermionAHSS/python/cochain_tools.py.
Copyright (c) 2026 koAHSS contributors; distributed under the MIT license.
Values are integers or fractions. Binary reduction and rational lifts are
explicit; ``cup`` returns a binary value unless ``integral=True`` is supplied.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache


@lru_cache(None)
def cut_terms(word, degrees, integral_index=None):
    labels = tuple(int(c) - 1 for c in word)
    degree = sum(degrees) - (len(labels) - len(degrees))
    if degree < 0 or any(d < 0 for d in degrees):
        return ()
    last = {c: max(j for j, z in enumerate(labels) if z == c) for c in set(labels)}
    result = Counter()

    def visit(j, previous, faces, cuts):
        c = labels[j]
        endpoints = range(previous, degree + 1) if j < len(labels) - 1 else (degree,)
        for endpoint in endpoints:
            interval = frozenset(range(previous, endpoint + 1))
            if faces[c] & interval:
                continue
            union = faces[c] | interval
            if len(union) > degrees[c] + 1 or (j == last[c] and len(union) != degrees[c] + 1):
                continue
            changed = list(faces)
            changed[c] = union
            extended = cuts + (endpoint,)
            if j + 1 != len(labels):
                visit(j + 1, endpoint, changed, extended)
                continue
            key = tuple(tuple(sorted(face)) for face in changed)
            sign = 1
            if integral_index is not None:
                p, q = degrees
                inner = [j != last[c] for j, c in enumerate(labels)]
                lengths = [extended[j + 1] - extended[j] + inner[j] for j in range(len(labels))]
                epsilon = sum(extended[j + 1] for j in range(len(labels)) if inner[j])
                epsilon += sum(lengths[j] * lengths[k] for j in range(len(labels))
                               for k in range(j + 1, len(labels)) if labels[j] > labels[k])
                epsilon += integral_index * (p + q) + integral_index * (integral_index - 1) // 2
                sign = (-1) ** epsilon
            result[key] += sign

    visit(0, 0, [frozenset() for _ in degrees], (0,))
    if integral_index is None:
        return tuple((faces, value % 2) for faces, value in result.items() if value % 2)
    return tuple((faces, value) for faces, value in result.items() if value)


class Cochain:
    """A lazy, exact cochain whose values use the first-vertex trivialization."""

    def __init__(self, degree, evaluate, name=""):
        self.degree = degree
        self.evaluate = lru_cache(None)(evaluate)
        self.name = name

    def __call__(self, simplex):
        simplex = tuple(simplex)
        if self.degree < 0:
            return 0
        if len(simplex) != self.degree + 1:
            raise ValueError("simplex size does not match cochain degree")
        return self.evaluate(simplex)

    def __add__(self, other):
        if self.degree != other.degree:
            raise ValueError("cochain degrees differ")
        return Cochain(self.degree, lambda t: self(t) + other(t))

    def __sub__(self, other):
        return self + (-other)

    def __neg__(self):
        return Cochain(self.degree, lambda t: -self(t))

    def __mul__(self, scalar):
        if not isinstance(scalar, (int, Fraction)):
            return NotImplemented
        return Cochain(self.degree, lambda t: scalar * self(t))

    __rmul__ = __mul__

    def mod2(self):
        return Cochain(self.degree, lambda t: self(t) % 2)

    def mod1(self):
        return Cochain(self.degree, lambda t: self(t) % 1)


def zero(degree):
    return Cochain(degree, lambda t: 0)


def binary_sum(*cochains):
    if not cochains:
        raise ValueError("binary_sum requires at least one cochain")
    result = cochains[0]
    for cochain in cochains[1:]:
        result = result + cochain
    return result.mod2()


def word_op(word, *cochains, integral_index=None):
    degree = sum(c.degree for c in cochains) - (len(word) - len(cochains))
    terms = cut_terms(word, tuple(c.degree for c in cochains), integral_index)

    def evaluate(simplex):
        total = 0
        for faces, weight in terms:
            value = weight
            for c, face in zip(cochains, faces):
                value *= c(tuple(simplex[j] for j in face))
                if not value:
                    break
            total += value
        return total if integral_index is not None else total % 2

    return Cochain(degree, evaluate)


def cup(a, b, index=0, integral=False):
    if index < 0:
        return zero(a.degree + b.degree - index)
    word = "".join(str(1 + j % 2) for j in range(index + 2))
    return word_op(word, a, b, integral_index=index if integral else None)


def differential(a):
    if a.degree < 0:
        return zero(a.degree + 1)
    return Cochain(a.degree + 1, lambda t: sum((-1) ** j * a(t[:j] + t[j + 1:])
                                              for j in range(len(t))))


def signed_differential(a, s):
    """Coboundary with transport (-1)^s along the first edge; requires ds=0."""
    if s.degree != 1:
        raise ValueError("s must have degree one")
    if a.degree < 0:
        return zero(a.degree + 1)
    return differential(a) - 2 * cup(s.mod2(), a, integral=True)
