"""Exact interval-cut cochain operations and the lexicographic certificate.

All arithmetic is integral or Boolean; no external packages are required.
This module is the single implementation of the interval-cut engine: the
vendored stacking_model/cochains.py wraps these functions in its validating
Cochain class, and the page and extension workers both evaluate through it.
"""
from collections import Counter
import functools
from functools import lru_cache
from itertools import combinations
import json
from operator import attrgetter, itemgetter
import time


@lru_cache(None)
def cut_terms(word, degrees, integral_index=None):
    labels = tuple(int(c) - 1 for c in word)
    degree = sum(degrees) - (len(labels) - len(degrees))
    # Cochains of negative degree vanish, so their words have no terms.
    if degree < 0 or any(d < 0 for d in degrees):
        return ()
    last = {c: max(j for j, z in enumerate(labels) if z == c) for c in set(labels)}
    result = Counter()

    def visit(j, previous, faces, cuts):
        c = labels[j]
        for endpoint in (range(previous, degree + 1) if j < len(labels) - 1 else (degree,)):
            interval = frozenset(range(previous, endpoint + 1))
            # Repeated within-label vertices cannot satisfy the dimension sum.
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


_lru_cache_wrapper = getattr(functools, '_lru_cache_wrapper', None)
_CacheInfo = getattr(functools, '_CacheInfo', None)


def _memoize(evaluate):
    # The same unbounded C memo table that lru_cache(None) builds, without
    # copying metadata. A replaced module-level lru_cache (runtime_cache.py)
    # is honoured.
    if lru_cache is functools.lru_cache and _lru_cache_wrapper is not None:
        wrapper = _lru_cache_wrapper(evaluate, None, False, _CacheInfo)
        wrapper.__wrapped__ = evaluate
        return wrapper
    return lru_cache(None)(evaluate)


class Cochain:
    def __init__(self, degree, evaluate, name=""):
        self.degree = degree
        self.evaluate = _memoize(evaluate)
        self.name = name

    # Calling a cochain enters its memo table without a Python frame.
    __call__ = property(attrgetter('evaluate'))

    def __add__(self, other):
        assert self.degree == other.degree
        return Cochain(self.degree, lambda t: self(t) + other(t))

    def __sub__(self, other):
        assert self.degree == other.degree
        return Cochain(self.degree, lambda t: self(t) - other(t))

    def mod2(self):
        return Cochain(self.degree, lambda t: self(t) % 2)


def zero(degree):
    return Cochain(degree, lambda t: 0)


def face_getter(face):
    """A C-level map from a simplex to its face tuple."""
    if len(face) == 1:
        return itemgetter(slice(face[0], face[0] + 1))
    return itemgetter(*face)


def word_degree(word, cochains):
    return sum(c.degree for c in cochains) - (len(word) - len(cochains))


def word_evaluator(word, cochains, integral_index=None):
    """Evaluation function of an interval-cut word on the given cochains."""
    terms = cut_terms(word, tuple(c.degree for c in cochains), integral_index)
    prepared = tuple((weight, tuple((c, face_getter(face))
                                    for c, face in zip(cochains, faces)))
                     for faces, weight in terms)
    integral = integral_index is not None

    def evaluate(simplex):
        if type(simplex) is not tuple:
            simplex = tuple(simplex)
        total = 0
        for weight, factors in prepared:
            value = weight
            for c, get in factors:
                value *= c(get(simplex))
                if not value:
                    break
            total += value
        return total if integral else total % 2

    return evaluate


def word_op(word, *cochains, integral_index=None):
    return Cochain(word_degree(word, cochains),
                   word_evaluator(word, cochains, integral_index))


def cup_word(index):
    """The interval-cut word of the cup-i product."""
    return "".join(str(1 + j % 2) for j in range(index + 2))


def cup(a, b, index=0, integral=False):
    if index < 0:
        return zero(a.degree + b.degree - index)
    return word_op(cup_word(index), a, b, integral_index=index if integral else None)


def coboundary_evaluator(a):
    """Evaluation function of the ordinary coboundary of a."""
    return lambda t: sum((-1) ** j * a(t[:j] + t[j + 1:]) for j in range(len(t)))


def differential(a):
    return Cochain(a.degree + 1, coboundary_evaluator(a))


def interval_pullback(c, scaled=False):
    """Normalized pullback to X times the interval, optionally times its last coordinate."""
    def value(vertices):
        base = tuple(v[0] for v in vertices)
        if any(a == b for a, b in zip(base, base[1:])):
            return 0
        result = c(base)
        return result * vertices[-1][1] if scaled else result
    return Cochain(c.degree, value)


def right_prism(c):
    """Signed right prism over the vertex pairs (v, 0) and (v, 1)."""
    return Cochain(c.degree - 1, lambda vertices: sum(
        (-1) ** j * c(tuple((v, 0) for v in vertices[:j + 1]) + tuple((v, 1) for v in vertices[j:]))
        for j in range(len(vertices))))


def Q(c, degree):
    """Cochain Steenrod square Q^degree, including its coboundary term."""
    return (cup(c, c, c.degree - degree)
            + cup(c, differential(c).mod2(), c.degree - degree + 1)).mod2()


def square(a, k):
    return cup(a, a, a.degree - k)


def q2(a):
    return Q(a, 2)


def zeta1(x, y):
    word = "1232" + "".join(str(4 - j % 2) for j in range(y.degree))
    return word_op(word, x, x, y, y)


def zeta2(x, y):
    word = "1231" + "".join(str(3 + j % 2) for j in range(y.degree + 1))
    return word_op(word, x, x, y, y)


def shuffled_cycle(base_degree, interval_degree=8):
    """EZ shuffle of [g|...|g] and the oriented interval simplex; mod-two signs."""
    for base_steps in combinations(range(base_degree + interval_degree), base_degree):
        selected = set(base_steps)
        g = level = 0
        vertices = [(g, level)]
        for step in range(base_degree + interval_degree):
            if step in selected:
                g ^= 1
            else:
                level += 1
            vertices.append((g, level))
        yield tuple(vertices)


class Polynomial(frozenset):
    """Boolean ANF over F2: each integer monomial is a variable-support mask."""
    def __xor__(self, other):
        return Polynomial(self.symmetric_difference(other))

    def __and__(self, other):
        out = set()
        for a in self:
            for b in other:
                value = a | b
                if value in out:
                    out.remove(value)
                else:
                    out.add(value)
        return Polynomial(out)


def primary_lex_certificate(n=9):
    """Prove first support of the CLOSED cochain Sq2 Sq1 is suspended x.

    Only two 1-bits occur in the proposed minimum. We certify that Q is
    identically zero below its first 1-bit, evaluate the two intervening
    codes, then evaluate the minimum. This proves the global lex minimum.
    """
    r = n + 3
    coordinates = list(combinations(range(1, r + 1), n))
    index = {f: j for j, f in enumerate(coordinates)}
    first_face = (1,) + tuple(range(5, r + 1))
    first_index = index[first_face]
    exponent = len(coordinates) - first_index - 1
    minimum = 2 ** exponent + 2

    def polynomial_Q(fixed_prefix_zero):
        @lru_cache(None)
        def input_value(face):
            if face[0] == 0:
                j = index[face[1:]]
                return Polynomial() if j < fixed_prefix_zero else Polynomial([1 << j])
            result = Polynomial()
            for j in range(len(face)):
                result ^= input_value((0,) + face[:j] + face[j + 1:])
            return result

        def polynomial_cup(p, q, i, left, right):
            terms = cut_terms("".join(str(1 + j % 2) for j in range(i + 2)), (p, q))

            @lru_cache(None)
            def evaluate(face):
                result = Polynomial()
                for (u, v), weight in terms:
                    result ^= left(tuple(face[j] for j in u)) & right(tuple(face[j] for j in v))
                return result
            return evaluate

        e = polynomial_cup(n, n, n - 1, input_value, input_value)
        Q = polynomial_cup(n + 1, n + 1, n - 1, e, e)
        return Q(tuple(range(r + 1)))

    remainder = polynomial_Q(first_index + 1)
    assert not remainder, "Q has support earlier than the proposed leading bit"

    def Q_on_code(code):
        @lru_cache(None)
        def input_value(face):
            if face[0] == 0:
                j = index[face[1:]]
                return (code >> (len(coordinates) - j - 1)) & 1
            return sum(input_value((0,) + face[:j] + face[j + 1:])
                       for j in range(len(face))) % 2
        a = Cochain(n, input_value)
        return square(square(a, 1), 2)(tuple(range(r + 1)))

    evaluations = [Q_on_code(2 ** exponent + j) for j in range(3)]
    assert evaluations == [0, 0, 1], evaluations
    return {"n": n, "independent_coordinates": len(coordinates),
            "zero_prefix_length": first_index + 1,
            "restricted_ANF_monomials": len(remainder),
            "evaluations_at_leading_bit_plus_0_1_2": evaluations,
            "first_nonzero_simplex_code": minimum,
            "conclusion": "lexicographically least Adem primitive is zero on this simplex"}


