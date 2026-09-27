"""Exact cochains on ordered simplices and normalized interval-cut cups.

The interval-cut engine is the package's python/cochain_tools.py; this module
wraps it in a validating Cochain class used by the degree-three formulas.
Copyright (c) 2026 koAHSS contributors; distributed under the MIT license.
Values are integers or fractions. Binary reduction and rational lifts are
explicit; ``cup`` returns a binary value unless ``integral=True`` is supplied.
"""
from fractions import Fraction
from functools import lru_cache

from cochain_tools import (coboundary_evaluator, cup_word, cut_terms,
                           word_degree, word_evaluator)


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
    return Cochain(word_degree(word, cochains),
                   word_evaluator(word, cochains, integral_index))


def cup(a, b, index=0, integral=False):
    if index < 0:
        return zero(a.degree + b.degree - index)
    return word_op(cup_word(index), a, b, integral_index=index if integral else None)


def differential(a):
    if a.degree < 0:
        return zero(a.degree + 1)
    return Cochain(a.degree + 1, coboundary_evaluator(a))


def signed_differential(a, s):
    """Coboundary with transport (-1)^s along the first edge; requires ds=0."""
    if s.degree != 1:
        raise ValueError("s must have degree one")
    if a.degree < 0:
        return zero(a.degree + 1)
    return differential(a) - 2 * cup(s.mod2(), a, integral=True)
