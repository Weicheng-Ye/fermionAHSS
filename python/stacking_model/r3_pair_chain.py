"""Shared-sign two-fiber U2 bar contraction.

Adapted from fermionAHSS/python/r3_chain.py, under the repository MIT
license and original attribution. Existing calibration files are unchanged.

Finite crossed U2 bar contraction for the current Danus degree-three source.

This is r3_source.md §§6--8, including both perturbation series.  A small
letter is (component, divided_power_index); (1,0) denotes sigma-1.  No
multiplicativity or equivariance of the U2 contraction is assumed.
"""
from dataclasses import dataclass
from functools import lru_cache
from itertools import product

try:
    from . import chain_models as cm
except ImportError:
    import chain_models as cm

add, linear, sign = cm.add, cm.linear, cm.sign


@dataclass(frozen=True)
class U2:
    sigma: int
    matrix: tuple
    matrix_prime: tuple

    def __post_init__(self):
        n = len(self.matrix)
        if self.sigma not in (0, 1) or len(self.matrix_prime) != n:
            raise ValueError("shared U2 needs equal matrix dimensions and binary sigma")
        for matrix in (self.matrix, self.matrix_prime):
            if any(len(row) != n for row in matrix):
                raise ValueError("shared U2 matrices must be square")
            if any(matrix[i][i] or any(matrix[j][i] != -matrix[i][j]
                                      for j in range(n)) for i in range(n)):
                raise ValueError("shared U2 matrices must be integrally skew")


@dataclass(frozen=True)
class Diag3:
    rows: tuple

    def __post_init__(self):
        if any(len(a.matrix) != len(self.rows) for a in self.rows):
            raise ValueError("degree-n D3 needs n rows of degree n")


@dataclass(frozen=True)
class Rect3:
    q: int
    rows: tuple


def zero_u2(n, sigma=0):
    return U2(sigma, ((0,) * n,) * n, ((0,) * n,) * n)


def zero_diag(n):
    return Diag3((zero_u2(n),) * n)


def _inner(a, second=False):
    matrix = a.matrix_prime if second else a.matrix
    return cm.Diag('zsign', tuple((0, row) for row in matrix))


def _outer_pair(x, y, sigma=0):
    if len(x.rows) != len(y.rows):
        raise ValueError("unsigned U2 pair dimensions differ")
    n = len(x.rows)
    def matrix(diag):
        if any(s for s, _ in diag.rows):
            raise ArithmeticError("unsigned U2 factor acquired a sign")
        return tuple(tuple(diag.rows[i][1][j] if i < j else
                           -diag.rows[j][1][i] if i > j else 0
                           for j in range(n)) for i in range(n))
    return U2(sigma, matrix(x), matrix(y))


def _u2_pull(a, vertices):
    x, c = cm._pull(_inner(a), vertices)
    y, d = cm._pull(_inner(a, True), vertices)
    if c != 1 or d != 1:
        raise ArithmeticError("unsigned U2 pull acquired transport")
    return _outer_pair(x, y, a.sigma)


def _u2_degens(a):
    return cm._degen_indices(_inner(a)) & cm._degen_indices(_inner(a, True))


def _u2_identity(a):
    return not a.sigma and not any(any(row) for row in a.matrix + a.matrix_prime)


def _u2_mul(a, b):
    if len(a.matrix) != len(b.matrix):
        raise ValueError("simplicial multiplication needs equal dimensions")
    def combine(left, right):
        return tuple(tuple(x + sign(a.sigma) * y for x, y in zip(r, t))
                     for r, t in zip(left, right))
    return U2(a.sigma ^ b.sigma, combine(a.matrix, b.matrix),
              combine(a.matrix_prime, b.matrix_prime))


def _h_pull(x, vertices):
    coefficient = sign(sum(a.sigma for a in x.rows[:vertices[0]]))
    rows = []
    for start, stop in zip(vertices, vertices[1:]):
        a = zero_u2(x.q)
        for b in x.rows[start:stop]:
            a = _u2_mul(a, b)
        rows.append(a)
    return Rect3(x.q, tuple(rows)), coefficient


def _v_pull(x, vertices):
    return Rect3(len(vertices) - 1, tuple(_u2_pull(a, vertices) for a in x.rows))


def _rect_pull(x, hvertices, vvertices):
    y, coefficient = _h_pull(x, hvertices)
    return _v_pull(y, vvertices), coefficient


def _rect(x):
    return Rect3(len(x.rows), x.rows)


def _as_diag(x):
    if x.q != len(x.rows):
        raise ValueError("rectangle is not diagonal")
    return Diag3(x.rows)


def _diag_degens(x):
    return {i for i, a in enumerate(x.rows)
            if _u2_identity(a) and all(i in _u2_degens(b) for b in x.rows)}


def _rect_degenerate(x):
    if any(_u2_identity(a) for a in x.rows):
        return True
    return any(all(i in _u2_degens(a) for a in x.rows) for i in range(x.q))


def degree(x):
    if isinstance(x, Diag3):
        return len(x.rows)
    if isinstance(x, U2):
        return len(x.matrix)
    if isinstance(x, (cm.Atom, cm.Diag)):
        return cm._degree(x)
    if not x:
        return 0
    degrees = {degree(a) for a in x}
    if len(degrees) != 1:
        raise ValueError("product factors have unequal dimensions")
    return degrees.pop()


def degens(x):
    if isinstance(x, Diag3):
        return _diag_degens(x)
    if isinstance(x, U2):
        return _u2_degens(x)
    if isinstance(x, (cm.Atom, cm.Diag)):
        return cm._degen_indices(x)
    answer = set(range(degree(x)))
    for a in x:
        answer &= degens(a)
    return answer


def basis(x):
    return {} if degens(x) else {x: 1}


def pull(x, vertices):
    if isinstance(x, Diag3):
        y, coefficient = _rect_pull(_rect(x), vertices, vertices)
        return _as_diag(y), coefficient
    if isinstance(x, U2):
        return _u2_pull(x, vertices), 1
    if isinstance(x, (cm.Atom, cm.Diag)):
        return cm._pull(x, vertices)
    out, coefficient = [], 1
    for a in x:
        y, c = pull(a, vertices)
        out.append(y)
        coefficient *= c
    return tuple(out), coefficient


def boundary(x):
    n, out = degree(x), {}
    if degens(x):
        return out
    for i in range(n + 1) if n else ():
        y, c = pull(x, tuple(j for j in range(n + 1) if j != i))
        add(out, basis(y), sign(i) * c)
    return out


def product_AW(pair):
    x, y = pair
    n, out = degree(x), {}
    if degree(y) != n:
        raise ValueError("AW needs equal dimensions")
    for p in range(n + 1):
        a, c = pull(x, tuple(range(p + 1)))
        b, d = pull(y, tuple(range(p, n + 1)))
        if not degens(a) and not degens(b):
            add(out, {(a, b): c * d})
    return out


def product_shuffle(pair):
    x, y = pair
    out = {}
    if degens(x) or degens(y):
        return out
    for h, v, c in cm._shuffles(degree(x), degree(y)):
        a, ca = pull(x, h)
        b, cb = pull(y, v)
        add(out, basis((a, b)), c * ca * cb)
    return out


def product_homotopy(pair):
    x, y = pair
    n, out = degree(x), {}
    if degree(y) != n:
        raise ValueError("product homotopy needs equal dimensions")
    if degens(pair):
        return out
    for (h, v), c in cm._h_terms(n).items():
        a, ca = pull(x, h)
        b, cb = pull(y, v)
        add(out, basis((a, b)), c * ca * cb)
    return out


def _h_boundary(x):
    p, out = len(x.rows), {}
    for i in range(p + 1) if p else ():
        y, c = _h_pull(x, tuple(j for j in range(p + 1) if j != i))
        if not _rect_degenerate(y):
            add(out, {y: sign(i) * c})
    return out


def _v_boundary(x):
    out = {}
    for i in range(x.q + 1) if x.q else ():
        y = _v_pull(x, tuple(j for j in range(x.q + 1) if j != i))
        if not _rect_degenerate(y):
            add(out, {y: sign(i)})
    return out


def _rect_boundary(x):
    return add(_h_boundary(x), _v_boundary(x), sign(len(x.rows)))


@lru_cache(None)
def _F_D(x):
    if _diag_degens(x):
        return {}
    n, out = len(x.rows), {}
    for p in range(n + 1):
        y, c = _rect_pull(_rect(x), tuple(range(p + 1)), tuple(range(p, n + 1)))
        if not _rect_degenerate(y):
            add(out, {y: c})
    return out


@lru_cache(None)
def _G_D(x):
    if _rect_degenerate(x):
        return {}
    out = {}
    for h, v, c in cm._shuffles(len(x.rows), x.q):
        y, d = _rect_pull(x, h, v)
        y = _as_diag(y)
        if not _diag_degens(y):
            add(out, {y: c * d})
    return out


@lru_cache(None)
def _H_D(x):
    if _diag_degens(x):
        return {}
    out = {}
    for (h, v), c in cm._h_terms(len(x.rows)).items():
        y, d = _rect_pull(_rect(x), h, v)
        y = _as_diag(y)
        if not _diag_degens(y):
            add(out, {y: c * d})
    return out


def _multi_AW(factors):
    if len(factors) < 2:
        return {factors: 1} if not factors or not degens(factors[0]) else {}
    out = {}
    for (left, right), c in product_AW((factors[:-1], factors[-1])).items():
        for word, d in _multi_AW(left).items():
            add(out, {word + (right,): c * d})
    return out


def _multi_shuffle(word):
    if len(word) < 2:
        return {word: 1} if not word or not degens(word[0]) else {}
    out = {}
    for left, c in _multi_shuffle(word[:-1]).items():
        for (a, b), d in product_shuffle((left, word[-1])).items():
            add(out, {a + (b,): c * d})
    return out


def _multi_homotopy(factors):
    if len(factors) < 2:
        return {}
    out = {}
    pair = (factors[:-1], factors[-1])
    for (left, right), c in product_homotopy(pair).items():
        add(out, {left + (right,): c})
    for (left, right), c in product_AW(pair).items():
        for a, d in _multi_homotopy(left).items():
            for (b, e), f in product_shuffle((a, right)).items():
                add(out, {b + (e,): c * d * f})
    return out


def _rect_f0(x):
    return {word: c for word, c in _multi_AW(x.rows).items()
            if not any(not a.matrix and not a.sigma for a in word)}


def _rect_g0(word):
    if not word:
        return {Rect3(0, ()): 1}
    out = {}
    for factors, c in _multi_shuffle(word).items():
        y = Rect3(degree(factors), factors)
        if not _rect_degenerate(y):
            add(out, {y: c})
    return out


def _rect_h0(x):
    out = {}
    for factors, c in _multi_homotopy(x.rows).items():
        y = Rect3(x.q + 1, factors)
        if not _rect_degenerate(y):
            add(out, {y: sign(len(x.rows)) * c})
    return out


def _bar_sign(word):
    p = len(word)
    return sign(p + sum((p - j - 1) * len(a.matrix) for j, a in enumerate(word)))


@lru_cache(None)
def _F_V(x):
    out, term = {}, {x: 1}
    for _ in range(len(x.rows) + 1):
        for y, c in term.items():
            for word, d in _rect_f0(y).items():
                add(out, {word: c * d * _bar_sign(word)})
        term = {a: -c for a, c in linear(linear(term, _rect_h0), _h_boundary).items()}
        if not term:
            return out
    raise AssertionError("vertical perturbation failed to terminate")


@lru_cache(None)
def _H_V(x):
    out, term = {}, {x: 1}
    for _ in range(len(x.rows) + 1):
        add(out, linear(term, _rect_h0))
        term = {a: -c for a, c in linear(linear(term, _rect_h0), _h_boundary).items()}
        if not term:
            return out
    raise AssertionError("vertical homotopy perturbation failed to terminate")


def _G_V(word):
    return {x: _bar_sign(word) * c for x, c in _rect_g0(word).items()}


@lru_cache(None)
def _F_bar(x):
    return linear(_F_D(x), _F_V)


@lru_cache(None)
def _G_bar(word):
    return linear(_G_V(word), _G_D)


@lru_cache(None)
def _H_bar(x):
    return add(dict(_H_D(x)), linear(linear(_F_D(x), _H_V), _G_D))


@lru_cache(None)
def unit_p(a):
    if not a.matrix:
        return {(1, 0, 0): 1} if a.sigma else {}
    out = {}
    for (x, y), coefficient in product_AW((_inner(a), _inner(a, True))).items():
        for (p, k), left in cm.F(x).items():
            for (q, l), right in cm.F(y).items():
                if p or q or not k + l:
                    raise ArithmeticError("unsigned U2 pair projection left divided powers")
                add(out, {(a.sigma, k, l): coefficient * left * right})
    return out


@lru_cache(None)
def unit_i(letter):
    sigma, k, l = letter
    if k + l == 0:
        if sigma != 1:
            raise ValueError("only sigma-1 is an augmentation letter")
        return {zero_u2(0, 1): 1}
    out = {}
    for x, left in cm.G('zsign', (0, k)).items():
        for y, right in cm.G('zsign', (0, l)).items():
            for (a, b), coefficient in product_shuffle((x, y)).items():
                add(out, basis(_outer_pair(a, b, sigma)), coefficient * left * right)
    return out


def _unit_Q(a):
    return add(basis(a), linear(unit_p(a), unit_i), -1) if a.matrix else {}


@lru_cache(None)
def _unit_raw(a):
    if not a.matrix:
        return {}
    pair = (_inner(a), _inner(a, True))
    chain = dict(product_homotopy(pair))
    for (x, y), coefficient in product_AW(pair).items():
        for left, weight in cm.H(x).items():
            add(chain, product_shuffle((left, y)), coefficient * weight)
        for left, weight in cm.linear(cm.F(x), lambda w: cm.G('zsign', w)).items():
            for right, other in cm.H(y).items():
                add(chain, product_shuffle((left, right)),
                    coefficient * sign(degree(x)) * weight * other)
    out = {}
    for (x, y), coefficient in chain.items():
        add(out, basis(_outer_pair(x, y, a.sigma)), coefficient)
    return out


def _unit_u(a):
    return linear(linear(_unit_Q(a), _unit_raw), _unit_Q)


@lru_cache(None)
def unit_h(a):
    return linear(linear(_unit_u(a), boundary), _unit_u)


def _tensor_p(word):
    out = {(): 1}
    for a in word:
        nxt = {}
        for prefix, c in out.items():
            for letter, d in unit_p(a).items():
                add(nxt, {prefix + (letter,): c * d})
        out = nxt
    return out


def _tensor_i(word):
    out = {(): 1}
    for letter in word:
        nxt = {}
        for prefix, c in out.items():
            for a, d in unit_i(letter).items():
                add(nxt, {prefix + (a,): c * d})
        out = nxt
    return out


def _tensor_h(word):
    out, prefix, prefix_degree = {}, {(): 1}, 0
    for j, a in enumerate(word):
        for left, c in prefix.items():
            for b, d in unit_h(a).items():
                add(out, {left + (b,) + word[j + 1:]: sign(prefix_degree + 1) * c * d})
        nxt = {}
        for left, c in prefix.items():
            for b, d in linear(unit_p(a), unit_i).items():
                add(nxt, {left + (b,): c * d})
        prefix = nxt
        prefix_degree += len(a.matrix) + 1
    return out


def _raw_expansion(a):
    out = {a: 1}
    if not a.matrix and a.sigma:
        out[zero_u2(0)] = -1
    return out


@lru_cache(None)
def _atom_product(a, b):
    out = {}
    for left, c in _raw_expansion(a).items():
        for right, d in _raw_expansion(b).items():
            for h, v, e in cm._shuffles(len(left.matrix), len(right.matrix)):
                z = _u2_mul(_u2_pull(left, h), _u2_pull(right, v))
                if (z.matrix or z.sigma) and not degens(z):
                    add(out, {z: c * d * e})
    return out


def _bar_delta(word):
    out = {}
    if word and not word[0].matrix and word[0].sigma:
        add(out, {word[1:]: 2})
    internal_degree = 0
    for j in range(len(word) - 1):
        internal_degree += len(word[j].matrix)
        for a, c in _atom_product(word[j], word[j + 1]).items():
            add(out, {word[:j] + (a,) + word[j + 2:]: sign(j + internal_degree) * c})
    return out


def bar_boundary(word):
    out = _bar_delta(word)
    internal_degree = 0
    for j, a in enumerate(word):
        for b, c in boundary(a).items():
            add(out, {word[:j] + (b,) + word[j + 1:]: sign(j + 1 + internal_degree) * c})
        internal_degree += len(a.matrix)
    return out


@lru_cache(None)
def transfer_F(word):
    out, term = {}, {word: 1}
    for _ in range(len(word) + 1):
        add(out, linear(term, _tensor_p))
        term = {a: -c for a, c in linear(linear(term, _tensor_h), _bar_delta).items()}
        if not term:
            return out
    raise AssertionError("U2 projection perturbation failed to terminate")


@lru_cache(None)
def transfer_G(word):
    out, term = {}, _tensor_i(word)
    for _ in range(len(word) + 1):
        add(out, term)
        term = {a: -c for a, c in linear(linear(term, _bar_delta), _tensor_h).items()}
        if not term:
            return out
    raise AssertionError("U2 inclusion perturbation failed to terminate")


@lru_cache(None)
def transfer_H(word):
    out, term = {}, {word: 1}
    for _ in range(len(word) + 1):
        add(out, linear(term, _tensor_h))
        term = {a: -c for a, c in linear(linear(term, _tensor_h), _bar_delta).items()}
        if not term:
            return out
    raise AssertionError("U2 homotopy perturbation failed to terminate")


@lru_cache(None)
def small_boundary(word):
    return linear(linear(transfer_G(word), _bar_delta), _tensor_p)


@lru_cache(None)
def F(x):
    return linear(_F_bar(x), transfer_F)


@lru_cache(None)
def G(word):
    return linear(transfer_G(word), _G_bar)


@lru_cache(None)
def H(x):
    return add(dict(_H_bar(x)), linear(linear(_F_bar(x), transfer_H), _G_bar))


def word_degree(word):
    return sum(2 * (k + l) + 1 for sigma, k, l in word)


@lru_cache(None)
def words(n):
    if n == 0:
        return ((),)
    out = []
    for total in range((n - 1) // 2 + 1):
        letters = (((1, 0, 0),) if total == 0 else
                   tuple((sigma, k, total-k) for sigma in (0, 1)
                         for k in range(total + 1)))
        for letter in letters:
            for tail in words(n - 2 * total - 1):
                out.append((letter,) + tail)
    return tuple(sorted(out, key=lambda w: (len(w), w)))


@lru_cache(None)
def omega_words(n):
    if n == 0:
        return ((),)
    return tuple(sorted(((k,) + tail for k in range(1, n)
                         for tail in omega_words(n - k - 1)), key=lambda w: (len(w), w)))


def total_basis(n, relative=True):
    out = []
    for p in range(n + 1):
        for word in words(p):
            if relative and (not any(k for _, k, _ in word) or not any(l for _, _, l in word)):
                continue
            for omega in omega_words(n - p):
                out.append((word, omega))
    return tuple(sorted(out, key=lambda b: (word_degree(b[0]), len(b[0]), b[0],
                                           sum(k + 1 for k in b[1]), len(b[1]), b[1])))


def omega_boundary(word):
    atoms = tuple(cm.Atom('c2', 0, (1,) * n) for n in word)
    return {tuple(len(a.edges) for a in w): c for w, c in cm._bar_boundary(atoms).items()}


def total_boundary(b):
    word, omega = b
    out = {(w, omega): c for w, c in small_boundary(word).items()}
    for w, c in omega_boundary(omega).items():
        add(out, {(word, w): sign(word_degree(word)) * c})
    return out


def Ftot(pair):
    out = {}
    for (x, y), c in product_AW(pair).items():
        for a, d in F(x).items():
            for b, e in cm.F(y).items():
                add(out, {(a, b): c * d * e})
    return out


@lru_cache(None)
def Gtot(b):
    word, omega = b
    out = {}
    for x, c in G(word).items():
        for y, d in cm.G('c2', omega).items():
            add(out, product_shuffle((x, y)), c * d)
    return out


def Htot(pair):
    out = dict(product_homotopy(pair))
    for (x, y), c in product_AW(pair).items():
        for a, d in H(x).items():
            add(out, product_shuffle((a, y)), c * d)
        for a, d in linear(F(x), G).items():
            for b, e in cm.H(y).items():
                add(out, product_shuffle((a, b)), c * sign(degree(x)) * d * e)
    return out
