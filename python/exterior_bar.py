"""Exact crossed exterior bar comparison from Appendix E.5--E.6.

This is a bounded mathematical verification prototype, not production GAP
code.  The coefficient algebra has basis ``1,s,e,es``, with s*e=-e*s,
s*s=1 and e*e=0.  Reduced bar letters are ``t,e,es``, where t=s-1;
in particular t*t=-2*t.  A word is a tuple of these letters.

Public small-complex API (all values are sparse integer dictionaries):

* differential(word): sign-coefficient bar boundary, keyed by words;
* small_boundary(p, k): boundary in the small complex, keyed by (p,k);
* f_small(word), g_small(p,k), h_small(word): comparison and homotopy;
* word_degree(word), words_through(max_degree): degree and finite inventory.

For independent verification, free_q_boundary, free_p_boundary, h_q, h_p,
f_full, g_full and k_full act on one free basis element.  Q basis elements
are (coefficient, word); P basis elements are (coefficient,p,k).  The
homotopy k_full is graded E-linear: K(a*q)=(-1)^|a| a*K(q).

After tensoring with the right sign module, coefficients 1,s,e,es act as
1,-1,0,0, respectively.  No floating-point arithmetic or external package
is used.  The private cached dictionaries are never exposed or mutated.
"""

from functools import lru_cache
from typing import Dict, Hashable, Iterator, Mapping, Tuple, TypeVar


Word = Tuple[str, ...]
QBasis = Tuple[str, Word]
PBasis = Tuple[str, int, int]
SmallBasis = Tuple[int, int]
Key = TypeVar("Key", bound=Hashable)
COEFFICIENTS = ("1", "s", "e", "es")
LETTERS = ("t", "e", "es")
_COEFFICIENT_DEGREE = {"1": 0, "s": 0, "e": 1, "es": 1}
_LETTER_DEGREE = {"t": 0, "e": 1, "es": 1}
_COEFFICIENT_BITS = {"1": (0, 0), "s": (0, 1), "e": (1, 0), "es": (1, 1)}
_BITS_COEFFICIENT = {bits: name for name, bits in _COEFFICIENT_BITS.items()}
_SIGN_ACTION = {"1": 1, "s": -1, "e": 0, "es": 0}
_LETTER_EXPANSION = {"t": {"s": 1, "1": -1}, "e": {"e": 1}, "es": {"es": 1}}


def _add(target: Dict[Key, int], key: Key, value: int) -> None:
    """Accumulate signed coefficients without retaining zero terms."""
    updated = target.get(key, 0) + value
    if updated:
        target[key] = updated
    else:
        target.pop(key, None)


def _extend(target: Dict[Key, int], source: Mapping[Key, int], scale: int = 1) -> None:
    for key, value in source.items():
        _add(target, key, scale * value)


def _validate_coefficient(coefficient: str) -> None:
    if coefficient not in _COEFFICIENT_DEGREE:
        raise ValueError(f"Unknown exterior coefficient: {coefficient!r}")


def _validate_word(word: Word) -> None:
    if not isinstance(word, tuple) or any(letter not in LETTERS for letter in word):
        raise ValueError("A bar word must be a tuple with letters 't', 'e', or 'es'.")


def _validate_indices(p: int, k: int) -> None:
    if not isinstance(p, int) or not isinstance(k, int) or p < 0 or k < 0:
        raise ValueError("The small-complex indices p and k must be nonnegative integers.")


def coefficient_product(left: str, right: str) -> Dict[str, int]:
    """Multiply two coefficient basis elements, including the reflection sign."""
    _validate_coefficient(left)
    _validate_coefficient(right)
    left_e, left_s = _COEFFICIENT_BITS[left]
    right_e, right_s = _COEFFICIENT_BITS[right]
    if left_e and right_e:
        return {}
    name = _BITS_COEFFICIENT[(left_e + right_e, left_s ^ right_s)]
    return {name: (-1) ** (left_s * right_e)}


def letter_product(left: str, right: str) -> Dict[str, int]:
    """Multiply augmentation-ideal letters and re-expand in t,e,es."""
    _validate_word((left, right))
    expanded: Dict[str, int] = {}
    for a, coefficient_a in _LETTER_EXPANSION[left].items():
        for b, coefficient_b in _LETTER_EXPANSION[right].items():
            _extend(expanded, coefficient_product(a, b), coefficient_a * coefficient_b)
    # The product lies in the augmentation ideal: its 1 coefficient is
    # the negative of its s coefficient.  Converting back must retain t.
    if expanded.get("1", 0) + expanded.get("s", 0):
        raise ArithmeticError("An augmentation-ideal product acquired nonzero augmentation.")
    return {letter: expanded[name] for letter, name in (("t", "s"), ("e", "e"), ("es", "es"))
            if expanded.get(name, 0)}


def word_degree(word: Word) -> int:
    """Return the total suspended degree of a reduced bar word."""
    _validate_word(word)
    return sum(1 + _LETTER_DEGREE[letter] for letter in word)


def words_through(max_degree: int) -> Iterator[Word]:
    """Enumerate every word through a finite degree, in degree/lexical order."""
    _validate_indices(max_degree, 0)
    by_degree = [[()]] + [[] for _ in range(max_degree)]
    for degree in range(1, max_degree + 1):
        for letter in LETTERS:
            remainder = degree - 1 - _LETTER_DEGREE[letter]
            if remainder >= 0:
                by_degree[degree].extend((letter,) + tail for tail in by_degree[remainder])
        by_degree[degree].sort()
    for words in by_degree:
        yield from words


def _left_q(coefficient: str, chain: Mapping[QBasis, int], scale: int = 1) -> Dict[QBasis, int]:
    result: Dict[QBasis, int] = {}
    for (right, word), value in chain.items():
        for merged, sign in coefficient_product(coefficient, right).items():
            _add(result, (merged, word), scale * value * sign)
    return result


def _left_p(coefficient: str, chain: Mapping[PBasis, int], scale: int = 1) -> Dict[PBasis, int]:
    result: Dict[PBasis, int] = {}
    for (right, p, k), value in chain.items():
        for merged, sign in coefficient_product(coefficient, right).items():
            _add(result, (merged, p, k), scale * value * sign)
    return result


@lru_cache(maxsize=None)
def _q_boundary(word: Word) -> Dict[QBasis, int]:
    result: Dict[QBasis, int] = {}
    if not word:
        return result
    for coefficient, value in _LETTER_EXPANSION[word[0]].items():
        _add(result, (coefficient, word[1:]), -value)
    prefix_degree = 0
    for index in range(len(word) - 1):
        prefix_degree += _LETTER_DEGREE[word[index]]
        sign = (-1) ** (index + prefix_degree)
        for letter, value in letter_product(word[index], word[index + 1]).items():
            merged = word[:index] + (letter,) + word[index + 2:]
            _add(result, ("1", merged), sign * value)
    return result


def free_q_boundary(basis: QBasis) -> Dict[QBasis, int]:
    """Free-bar differential, with d(a*q)=(-1)^|a| a*dq."""
    coefficient, word = basis
    _validate_coefficient(coefficient)
    _validate_word(word)
    return _left_q(coefficient, _q_boundary(word), (-1) ** _COEFFICIENT_DEGREE[coefficient])


def free_p_boundary(basis: PBasis) -> Dict[PBasis, int]:
    """Appendix E.5 free-resolution differential."""
    coefficient, p, k = basis
    _validate_coefficient(coefficient)
    _validate_indices(p, k)
    result: Dict[PBasis, int] = {}
    if k:
        _add(result, ("e", p, k - 1), (-1) ** p)
    if p:
        _add(result, ("s", p - 1, k), 1)
        _add(result, ("1", p - 1, k), (-1) ** (p + k))
    return _left_p(coefficient, result, (-1) ** _COEFFICIENT_DEGREE[coefficient])


def h_p(basis: PBasis) -> Dict[PBasis, int]:
    """The specified underlying integral augmentation contraction of P."""
    coefficient, p, k = basis
    _validate_coefficient(coefficient)
    _validate_indices(p, k)
    if coefficient == "s" and k == 0:
        return {("1", p + 1, 0): 1}
    if coefficient == "e":
        return {("1", p, k + 1): (-1) ** p}
    if coefficient == "es":
        return {("s", p, k + 1): (-1) ** (p + 1)}
    return {}


def h_q(basis: QBasis) -> Dict[QBasis, int]:
    """The specified integral augmentation contraction, -[a-epsilon(a)|w]."""
    coefficient, word = basis
    _validate_coefficient(coefficient)
    _validate_word(word)
    if coefficient == "1":
        return {}
    letter = "t" if coefficient == "s" else coefficient
    return {("1", (letter,) + word): -1}


@lru_cache(maxsize=None)
def _f_basis(word: Word) -> Dict[PBasis, int]:
    if not word:
        return {("1", 0, 0): 1}
    result: Dict[PBasis, int] = {}
    for boundary_basis, boundary_value in _q_boundary(word).items():
        for image_basis, image_value in f_full(boundary_basis).items():
            _extend(result, h_p(image_basis), boundary_value * image_value)
    return result


def f_full(basis: QBasis) -> Dict[PBasis, int]:
    """The recursively specified degree-zero E-linear map Q to P."""
    coefficient, word = basis
    _validate_coefficient(coefficient)
    _validate_word(word)
    return _left_p(coefficient, _f_basis(word))


@lru_cache(maxsize=None)
def _g_basis(p: int, k: int) -> Dict[QBasis, int]:
    if p == k == 0:
        return {("1", ()): 1}
    result: Dict[QBasis, int] = {}
    for boundary_basis, boundary_value in free_p_boundary(("1", p, k)).items():
        for image_basis, image_value in g_full(boundary_basis).items():
            _extend(result, h_q(image_basis), boundary_value * image_value)
    return result


def g_full(basis: PBasis) -> Dict[QBasis, int]:
    """The recursively specified degree-zero E-linear map P to Q."""
    coefficient, p, k = basis
    _validate_coefficient(coefficient)
    _validate_indices(p, k)
    return _left_q(coefficient, _g_basis(p, k))


@lru_cache(maxsize=None)
def _k_basis(word: Word) -> Dict[QBasis, int]:
    residual: Dict[QBasis, int] = {("1", word): 1}
    for basis, value in _f_basis(word).items():
        _extend(residual, g_full(basis), -value)
    for basis, value in _q_boundary(word).items():
        _extend(residual, k_full(basis), -value)
    result: Dict[QBasis, int] = {}
    for basis, value in residual.items():
        _extend(result, h_q(basis), value)
    return result


def k_full(basis: QBasis) -> Dict[QBasis, int]:
    """The degree-one E-linear comparison homotopy, including its Koszul sign."""
    coefficient, word = basis
    _validate_coefficient(coefficient)
    _validate_word(word)
    return _left_q(coefficient, _k_basis(word), (-1) ** _COEFFICIENT_DEGREE[coefficient])


def _project_q(chain: Mapping[QBasis, int]) -> Dict[Word, int]:
    result: Dict[Word, int] = {}
    for (coefficient, word), value in chain.items():
        _add(result, word, _SIGN_ACTION[coefficient] * value)
    return result


def _project_p(chain: Mapping[PBasis, int]) -> Dict[SmallBasis, int]:
    result: Dict[SmallBasis, int] = {}
    for (coefficient, p, k), value in chain.items():
        _add(result, (p, k), _SIGN_ACTION[coefficient] * value)
    return result


def differential(word: Word) -> Dict[Word, int]:
    """Boundary of one word after tensoring with the right sign module."""
    _validate_word(word)
    return _project_q(_q_boundary(word))


def small_boundary(p: int, k: int) -> Dict[SmallBasis, int]:
    """Small sign boundary: (-1+(-1)^(p+k))*u_(p-1,k)."""
    _validate_indices(p, k)
    if p and (p + k) % 2:
        return {(p - 1, k): -2}
    return {}


def f_small(word: Word) -> Dict[SmallBasis, int]:
    """Bar-to-small comparison after applying the sign character."""
    _validate_word(word)
    return _project_p(_f_basis(word))


def g_small(p: int, k: int) -> Dict[Word, int]:
    """Small-to-bar comparison after applying the sign character."""
    _validate_indices(p, k)
    return _project_q(_g_basis(p, k))


def h_small(word: Word) -> Dict[Word, int]:
    """Comparison homotopy with d*h+h*d=1-g*f on the sign bar complex."""
    _validate_word(word)
    return _project_q(_k_basis(word))
