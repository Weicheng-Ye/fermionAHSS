"""Exact sparse chain operators printed in the current Danus source notes.

No local cochain solver occurs here. All chains are normalized, integral,
and represented by dictionaries with zero coefficients removed.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations, product


def add(out, source, scale=1):
    for key, coefficient in source.items():
        value = out.get(key, 0) + scale * coefficient
        if value:
            out[key] = value
        else:
            out.pop(key, None)
    return out


def linear(chain, operation):
    answer = {}
    for key, coefficient in chain.items():
        add(answer, operation(key), coefficient)
    return answer


def sign(n):
    return -1 if n % 2 else 1


@dataclass(frozen=True)
class Diag:
    mode: str
    rows: tuple

    def __post_init__(self):
        if self.mode not in ("zsign", "c2"):
            raise ValueError("mode must be zsign or c2")
        n = len(self.rows)
        if any(s not in (0, 1) or len(e) != n for s, e in self.rows):
            raise ValueError("a diagonal simplex must have n rows of n edges")
        if self.mode == "c2" and any(s or any(e not in (0, 1) for e in v)
                                     for s, v in self.rows):
            raise ValueError("c2 rows have sigma zero and binary edges")


@dataclass(frozen=True)
class Rect:
    mode: str
    q: int
    rows: tuple


@dataclass(frozen=True)
class Atom:
    mode: str
    sigma: int
    edges: tuple


def zero_diag(mode, n):
    return Diag(mode, tuple((0, (0,) * n) for _ in range(n)))


def _sum_edges(mode, entries):
    value = sum(entries)
    return value % 2 if mode == "c2" else value


def _row_mul(mode, left, right):
    s, a = left
    t, b = right
    return (s ^ t, tuple(_sum_edges(mode, (x, sign(s) * y))
                         for x, y in zip(a, b)))


def _h_pull(x, vertices):
    """Horizontal pullback, together with its initial-frame transport."""
    rows = x.rows
    factor = sign(sum(s for s, _ in rows[:vertices[0]])) if x.mode == "zsign" else 1
    out = []
    for start, stop in zip(vertices, vertices[1:]):
        row = (0, (0,) * x.q)
        for entry in rows[start:stop]:
            row = _row_mul(x.mode, row, entry)
        out.append(row)
    return Rect(x.mode, x.q, tuple(out)), factor


def _v_pull(x, vertices):
    rows = tuple((s, tuple(_sum_edges(x.mode, edges[a:b])
                           for a, b in zip(vertices, vertices[1:])))
                 for s, edges in x.rows)
    return Rect(x.mode, len(vertices) - 1, rows)


def _rect_pull(x, hvertices, vvertices):
    y, factor = _h_pull(x, hvertices)
    return _v_pull(y, vvertices), factor


def _diag_rect(x):
    return Rect(x.mode, len(x.rows), x.rows)


def _as_diag(x):
    if len(x.rows) != x.q:
        raise ValueError("rectangle is not diagonal")
    return Diag(x.mode, x.rows)


def diag_face(x, i):
    n = len(x.rows)
    if not 0 <= i <= n or n == 0:
        raise ValueError("invalid face")
    vertices = tuple(j for j in range(n + 1) if j != i)
    return _as_diag(_rect_pull(_diag_rect(x), vertices, vertices)[0])


def diag_degeneracy(x, i):
    n = len(x.rows)
    if not 0 <= i <= n:
        raise ValueError("invalid degeneracy")
    vertices = tuple(range(i + 1)) + tuple(range(i, n + 1))
    return _as_diag(_rect_pull(_diag_rect(x), vertices, vertices)[0])


def _diag_degeneracies(x):
    n = len(x.rows)
    return {i for i, (s, edges) in enumerate(x.rows)
            if not s and not any(edges)
            and all(not row[1][i] for row in x.rows)}


def diag_boundary(x):
    if _diag_degeneracies(x):
        return {}
    n = len(x.rows)
    answer = {}
    for i in range(n + 1) if n else ():
        vertices = tuple(j for j in range(n + 1) if j != i)
        y, factor = _rect_pull(_diag_rect(x), vertices, vertices)
        y = _as_diag(y)
        if not _diag_degeneracies(y):
            add(answer, {y: sign(i) * factor})
    return answer


def _rect_degenerate(x):
    if any(not s and not any(e) for s, e in x.rows):
        return True
    return any(all(not e[i] for _, e in x.rows) for i in range(x.q))


def _h_boundary(x):
    p = len(x.rows)
    out = {}
    for i in range(p + 1) if p else ():
        v = tuple(j for j in range(p + 1) if i != j)
        y, factor = _h_pull(x, v)
        if not _rect_degenerate(y):
            add(out, {y: sign(i) * factor})
    return out


def _v_boundary(x):
    out = {}
    for i in range(x.q + 1) if x.q else ():
        v = tuple(j for j in range(x.q + 1) if i != j)
        y = _v_pull(x, v)
        if not _rect_degenerate(y):
            add(out, {y: sign(i)})
    return out


def _rect_boundary(x):
    return add(_h_boundary(x), _v_boundary(x), sign(len(x.rows)))


# Universal product homotopy: keys are pairs of nondecreasing vertex maps.
def _u_degenerate(pair):
    a, b = pair
    return any(a[i] == a[i + 1] and b[i] == b[i + 1]
               for i in range(len(a) - 1))


def _u_chain(pair):
    return {} if _u_degenerate(pair) else {pair: 1}


@lru_cache(None)
def _shuffles(p, q):
    result = []
    for positions in combinations(range(p + q), p):
        hs = set(positions)
        a, b = [0], [0]
        inversion = 0
        for j in range(p + q):
            if j in hs:
                inversion += b[-1]
                a.append(a[-1] + 1)
                b.append(b[-1])
            else:
                a.append(a[-1])
                b.append(b[-1] + 1)
        result.append((tuple(a), tuple(b), sign(inversion)))
    return tuple(result)


def _u_boundary(pair):
    n = len(pair[0]) - 1
    out = {}
    for i in range(n + 1) if n else ():
        y = tuple(v[:i] + v[i + 1:] for v in pair)
        add(out, _u_chain(y), sign(i))
    return out


def _u_gf(pair):
    a, b = pair
    n = len(a) - 1
    out = {}
    for p in range(n + 1):
        # A factor of a tensor may be degenerate independently.
        left, right = a[:p + 1], b[p:]
        if any(x == y for x, y in zip(left, left[1:])) or any(
                x == y for x, y in zip(right, right[1:])):
            continue
        for h, v, coefficient in _shuffles(p, n - p):
            y = (tuple(left[j] for j in h), tuple(right[j] for j in v))
            add(out, _u_chain(y), coefficient)
    return out


def _u_Q(pair):
    return add(_u_chain(pair), _u_gf(pair), -1)


def _u_apply_terms(pair, terms):
    out = {}
    for (h, v), coefficient in terms.items():
        y = (tuple(pair[0][j] for j in h),
             tuple(pair[1][j] for j in v))
        add(out, _u_chain(y), coefficient)
    return out


@lru_cache(None)
def _raw_h_terms(n):
    if n == 0:
        return {}
    identity = tuple(range(n + 1))
    cycle = _u_Q((identity, identity))
    for i in range(n + 1):
        face = tuple(j for j in identity if j != i)
        add(cycle, _u_apply_terms((face, face), _raw_h_terms(n - 1)), -sign(i))
    out = {}
    for (a, b), coefficient in cycle.items():
        add(out, _u_chain(((0,) + a, (0,) + b)), coefficient)
    return out


def _u_raw(pair):
    return _u_apply_terms(pair, _raw_h_terms(len(pair[0]) - 1))


def _u_u(pair):
    return linear(linear(_u_Q(pair), _u_raw), _u_Q)


@lru_cache(None)
def _h_terms(n):
    identity = tuple(range(n + 1))
    return linear(linear(_u_u((identity, identity)), _u_boundary), _u_u)


def _degree(x):
    if isinstance(x, Diag):
        return len(x.rows)
    if isinstance(x, Atom):
        return len(x.edges)
    if not x:
        return 0
    degrees = {_degree(y) for y in x}
    if len(degrees) != 1:
        raise ValueError("product factors have different degrees")
    return degrees.pop()


def _degen_indices(x):
    if isinstance(x, Diag):
        return _diag_degeneracies(x)
    if isinstance(x, Atom):
        return {i for i, v in enumerate(x.edges) if v == 0}
    n = _degree(x)
    answer = set(range(n))
    for factor in x:
        answer &= _degen_indices(factor)
    return answer


def _basis(x):
    return {} if _degen_indices(x) else {x: 1}


def _pull(x, vertices):
    if isinstance(x, Diag):
        y, factor = _rect_pull(_diag_rect(x), vertices, vertices)
        return _as_diag(y), factor
    if isinstance(x, Atom):
        return Atom(x.mode, x.sigma, tuple(
            _sum_edges(x.mode, x.edges[a:b])
            for a, b in zip(vertices, vertices[1:]))), 1
    out, coefficient = [], 1
    for factor in x:
        y, c = _pull(factor, vertices)
        out.append(y)
        coefficient *= c
    return tuple(out), coefficient


def _boundary(x):
    n = _degree(x)
    out = {}
    if _degen_indices(x):
        return out
    for i in range(n + 1) if n else ():
        y, coefficient = _pull(x, tuple(j for j in range(n + 1) if i != j))
        add(out, _basis(y), sign(i) * coefficient)
    return out


def product_AW(pair):
    x, y = pair
    n = _degree(x)
    if _degree(y) != n:
        raise ValueError("AW expects equal dimensions")
    out = {}
    for p in range(n + 1):
        a, ca = _pull(x, tuple(range(p + 1)))
        b, cb = _pull(y, tuple(range(p, n + 1)))
        if not _degen_indices(a) and not _degen_indices(b):
            add(out, {(a, b): ca * cb})
    return out


def product_shuffle(pair):
    x, y = pair
    if _degen_indices(x) or _degen_indices(y):
        return {}
    out = {}
    for h, v, coefficient in _shuffles(_degree(x), _degree(y)):
        a, ca = _pull(x, h)
        b, cb = _pull(y, v)
        add(out, _basis((a, b)), coefficient * ca * cb)
    return out


def product_homotopy(pair):
    x, y = pair
    n = _degree(x)
    if _degree(y) != n:
        raise ValueError("product homotopy expects equal dimensions")
    if _degen_indices(pair):
        return {}
    out = {}
    for (h, v), coefficient in _h_terms(n).items():
        a, ca = _pull(x, h)
        b, cb = _pull(y, v)
        add(out, _basis((a, b)), coefficient * ca * cb)
    return out


@lru_cache(None)
def _F_D(x):
    if _diag_degeneracies(x):
        return {}
    n = len(x.rows)
    out = {}
    for p in range(n + 1):
        y, factor = _rect_pull(_diag_rect(x), tuple(range(p + 1)),
                               tuple(range(p, n + 1)))
        if not _rect_degenerate(y):
            add(out, {y: factor})
    return out


@lru_cache(None)
def _G_D(x):
    if _rect_degenerate(x):
        return {}
    out = {}
    for h, v, coefficient in _shuffles(len(x.rows), x.q):
        y, factor = _rect_pull(x, h, v)
        y = _as_diag(y)
        if not _diag_degeneracies(y):
            add(out, {y: coefficient * factor})
    return out


@lru_cache(None)
def _H_D(x):
    if _diag_degeneracies(x):
        return {}
    out = {}
    for (h, v), coefficient in _h_terms(len(x.rows)).items():
        y, factor = _rect_pull(_diag_rect(x), h, v)
        y = _as_diag(y)
        if not _diag_degeneracies(y):
            add(out, {y: coefficient * factor})
    return out


def _multi_AW(factors):
    if not factors:
        return {(): 1}
    if len(factors) == 1:
        return {} if _degen_indices(factors[0]) else {factors: 1}
    out = {}
    for (left, right), coefficient in product_AW((factors[:-1], factors[-1])).items():
        for word, c in _multi_AW(left).items():
            add(out, {word + (right,): coefficient * c})
    return out


def _multi_shuffle(word):
    if not word:
        return {(): 1}
    if len(word) == 1:
        return {} if _degen_indices(word[0]) else {word: 1}
    out = {}
    for left, coefficient in _multi_shuffle(word[:-1]).items():
        for (a, b), c in product_shuffle((left, word[-1])).items():
            add(out, {a + (b,): coefficient * c})
    return out


def _multi_homotopy(factors):
    if len(factors) < 2:
        return {}
    out = {}
    pair = (factors[:-1], factors[-1])
    for (left, right), coefficient in product_homotopy(pair).items():
        add(out, {left + (right,): coefficient})
    for (left, right), coefficient in product_AW(pair).items():
        for a, c in _multi_homotopy(left).items():
            for (b, d), e in product_shuffle((a, right)).items():
                add(out, {b + (d,): coefficient * c * e})
    return out


def _rect_f0(x):
    factors = tuple(Atom(x.mode, s, e) for s, e in x.rows)
    out = {}
    for word, coefficient in _multi_AW(factors).items():
        if not any(not a.edges and not a.sigma for a in word):
            add(out, {word: coefficient})
    return out


def _rect_g0(word):
    if not word:
        raise ValueError("empty bar requires mode; handled by caller")
    out = {}
    for factors, coefficient in _multi_shuffle(word).items():
        y = Rect(word[0].mode, _degree(factors),
                 tuple((a.sigma, a.edges) for a in factors))
        if not _rect_degenerate(y):
            add(out, {y: coefficient})
    return out


def _rect_h0(x):
    factors = tuple(Atom(x.mode, s, e) for s, e in x.rows)
    out = {}
    for factors, coefficient in _multi_homotopy(factors).items():
        y = Rect(x.mode, x.q + 1, tuple((a.sigma, a.edges) for a in factors))
        if not _rect_degenerate(y):
            add(out, {y: sign(len(x.rows)) * coefficient})
    return out


def _bar_sign(word):
    p = len(word)
    return sign(p + sum((p - j - 1) * len(a.edges) for j, a in enumerate(word)))


@lru_cache(None)
def _F_V(x):
    answer = {}
    term = {x: 1}
    # Each application decreases horizontal length.
    for _ in range(len(x.rows) + 1):
        for y, coefficient in term.items():
            for word, c in _rect_f0(y).items():
                add(answer, {word: coefficient * c * _bar_sign(word)})
        term = linear(linear(term, _rect_h0), _h_boundary)
        term = {a: -c for a, c in term.items()}
        if not term:
            return answer
    raise AssertionError("vertical perturbation failed to terminate")


@lru_cache(None)
def _H_V(x):
    answer = {}
    term = {x: 1}
    for _ in range(len(x.rows) + 1):
        add(answer, linear(term, _rect_h0))
        term = linear(linear(term, _rect_h0), _h_boundary)
        term = {a: -c for a, c in term.items()}
        if not term:
            return answer
    raise AssertionError("vertical perturbation failed to terminate")


def _G_V(mode, word):
    if not word:
        return {Rect(mode, 0, ()): 1}
    return {x: _bar_sign(word) * c for x, c in _rect_g0(word).items()}


@lru_cache(None)
def _F_bar(x):
    return linear(_F_D(x), _F_V)


@lru_cache(None)
def _G_bar(mode, word):
    return linear(_G_V(mode, word), _G_D)


@lru_cache(None)
def _H_bar(x):
    return add(dict(_H_D(x)), linear(linear(_F_D(x), _H_V), _G_D))


def _atom_boundary(a):
    return _boundary(a)


def _atom_p(a):
    if not a.edges:
        return {"t": 1} if a.sigma else {}
    if len(a.edges) == 1 and a.edges[0]:
        return {"es" if a.sigma else "e": a.edges[0]}
    return {}


def _atom_i(letter):
    return Atom("zsign", int(letter in ("t", "es")),
                () if letter == "t" else (1,))


def _unit_Q(a):
    if not a.edges:
        return {}
    return add(_basis(a), {_atom_i(k): v for k, v in _atom_p(a).items()}, -1)


def _unit_raw(a):
    n = len(a.edges)
    if not n or _degen_indices(a):
        return {}
    vertices = [0]
    for edge in a.edges:
        vertices.append(vertices[-1] + edge)
    start, stop = vertices[-2:]
    coefficient = sign(n) * (1 if stop >= start else -1)
    out = {}
    for j in range(min(start, stop), max(start, stop)):
        new = vertices[:-1] + [j, j + 1]
        b = Atom(a.mode, a.sigma, tuple(y - x for x, y in zip(new, new[1:])))
        add(out, _basis(b), coefficient)
    return out


def _unit_u(a):
    return linear(linear(_unit_Q(a), _unit_raw), _unit_Q)


@lru_cache(None)
def _unit_h(a):
    return linear(linear(_unit_u(a), _atom_boundary), _unit_u)


def _tensor_p(word):
    out = {(): 1}
    for a in word:
        next_out = {}
        for prefix, c in out.items():
            for letter, d in _atom_p(a).items():
                add(next_out, {prefix + (letter,): c * d})
        out = next_out
    return out


def _tensor_i(word):
    return {tuple(_atom_i(letter) for letter in word): 1}


def _tensor_h(word):
    out = {}
    prefix = {(): 1}
    prefix_degree = 0
    for j, a in enumerate(word):
        for left, c in prefix.items():
            for b, d in _unit_h(a).items():
                add(out, {left + (b,) + word[j + 1:]:
                          sign(prefix_degree + 1) * c * d})
        next_prefix = {}
        for left, c in prefix.items():
            for letter, d in _atom_p(a).items():
                add(next_prefix, {left + (_atom_i(letter),): c * d})
        prefix = next_prefix
        prefix_degree += len(a.edges) + 1
    return out


def _raw_atom_expansion(a):
    out = {a: 1}
    if not a.edges and a.sigma:
        out[Atom(a.mode, 0, ())] = -1
    return out


@lru_cache(None)
def _atom_product(a, b):
    out = {}
    for left, c in _raw_atom_expansion(a).items():
        for right, d in _raw_atom_expansion(b).items():
            for h, v, coefficient in _shuffles(len(left.edges), len(right.edges)):
                x, _ = _pull(left, h)
                y, _ = _pull(right, v)
                sigma, edges = _row_mul(a.mode, (x.sigma, x.edges),
                                       (y.sigma, y.edges))
                z = Atom(a.mode, sigma, edges)
                if (z.edges or z.sigma) and not _degen_indices(z):
                    add(out, {z: c * d * coefficient})
    return out


def _bar_delta(word):
    out = {}
    if word and not word[0].edges and word[0].sigma:
        # Minus the signed action of sigma-1.
        add(out, {word[1:]: 2})
    degree = 0
    for j in range(len(word) - 1):
        degree += len(word[j].edges)
        for a, c in _atom_product(word[j], word[j + 1]).items():
            add(out, {word[:j] + (a,) + word[j + 2:]: sign(j + degree) * c})
    return out


def _bar_boundary(word):
    out = _bar_delta(word)
    degree = 0
    for j, a in enumerate(word):
        for b, c in _atom_boundary(a).items():
            add(out, {word[:j] + (b,) + word[j + 1:]:
                      sign(j + 1 + degree) * c})
        degree += len(a.edges)
    return out


@lru_cache(None)
def _transfer_H(word):
    out, term = {}, {word: 1}
    for _ in range(len(word) + 1):
        add(out, linear(term, _tensor_h))
        term = linear(linear(term, _tensor_h), _bar_delta)
        term = {a: -c for a, c in term.items()}
        if not term:
            return out
    raise AssertionError("unit-edge perturbation failed to terminate")


@lru_cache(None)
def _transfer_I(word):
    out, term = {}, _tensor_i(word)
    for _ in range(len(word) + 1):
        add(out, term)
        term = linear(linear(term, _bar_delta), _tensor_h)
        term = {a: -c for a, c in term.items()}
        if not term:
            return out
    raise AssertionError("inclusion perturbation failed to terminate")


def _exterior():
    try:
        from . import exterior_bar
    except ImportError:
        import exterior_bar
    return exterior_bar


@lru_cache(None)
def F(x):
    raw = _F_bar(x)
    if x.mode == "c2":
        return {tuple(len(a.edges) for a in word): c for word, c in raw.items()}
    return linear(linear(raw, _tensor_p), _exterior().f_small)


@lru_cache(None)
def G(mode, smallbasis):
    if mode == "c2":
        word = tuple(Atom("c2", 0, (1,) * n) for n in smallbasis)
        return _G_bar(mode, word)
    words = _exterior().g_small(*smallbasis)
    return linear(linear(words, _transfer_I), lambda word: _G_bar(mode, word))


@lru_cache(None)
def H(x):
    out = dict(_H_bar(x))
    if x.mode == "c2":
        return out
    raw = _F_bar(x)
    add(out, linear(linear(raw, _transfer_H),
                    lambda word: _G_bar(x.mode, word)))
    extra = linear(linear(raw, _tensor_p), _exterior().h_small)
    add(out, linear(linear(extra, _transfer_I),
                    lambda word: _G_bar(x.mode, word)))
    return out
