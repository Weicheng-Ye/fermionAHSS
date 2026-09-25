"""Shared-sign Borel contraction for two signed integral degree-one cocycles.

Uses the exact integer circle/EZ homotopies of fermionAHSS (MIT attribution)
and a finite horizontal perturbation. Signed coefficient chains are used for rational sources; the two fibers
have one common inversion action.
"""
from dataclasses import dataclass
from functools import lru_cache
import r2_pair_chain as r2
cm, add, linear, sign = r2.cm, r2.add, r2.linear, r2.sign


@dataclass(frozen=True)
class Borel:
    signs: tuple
    left: tuple
    right: tuple

    def __post_init__(self):
        if len(self.left) != len(self.right) or any(s not in (0, 1) for s in self.signs):
            raise ValueError('invalid Borel rectangle')


def h_pull(x, vertices):
    frame = sign(sum(x.signs[:vertices[0]]))
    return Borel(tuple(sum(x.signs[a:b]) % 2 for a, b in zip(vertices, vertices[1:])),
                 tuple(frame * a for a in x.left), tuple(frame * a for a in x.right))


def v_pull(x, vertices):
    return Borel(x.signs, tuple(sum(x.left[a:b]) for a, b in zip(vertices, vertices[1:])),
                 tuple(sum(x.right[a:b]) for a, b in zip(vertices, vertices[1:])))


def rect_pull(x, hv, vv):
    return v_pull(h_pull(x, hv), vv)


def degree(x):
    if isinstance(x, Borel):
        if len(x.signs) != len(x.left):
            raise ValueError('diagonal expected')
        return len(x.signs)
    if isinstance(x, (cm.Atom, cm.Diag)):
        return cm._degree(x)
    return degree(x[0])


def degens(x):
    if isinstance(x, Borel):
        return {i for i, (s, a, b) in enumerate(zip(x.signs, x.left, x.right))
                if not s and not a and not b}
    if isinstance(x, (cm.Atom, cm.Diag)):
        return cm._degen_indices(x)
    return degens(x[0]) & degens(x[1])


def basis(x):
    return {} if degens(x) else {x: 1}


def pull(x, vertices):
    if isinstance(x, Borel):
        return rect_pull(x, vertices, vertices), sign(sum(x.signs[:vertices[0]]))
    if isinstance(x, (cm.Atom, cm.Diag)):
        return cm._pull(x, vertices)
    a, c = pull(x[0], vertices)
    b, d = pull(x[1], vertices)
    return (a, b), c * d


def boundary(x):
    if degens(x):
        return {}
    n, out = degree(x), {}
    for i in range(n + 1) if n else ():
        y, c = pull(x, tuple(j for j in range(n + 1) if j != i))
        add(out, basis(y), sign(i) * c)
    return out


def product_AW(pair):
    x, y = pair
    n, out = degree(x), {}
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
    for hv, vv, c in cm._shuffles(degree(x), degree(y)):
        a, ca = pull(x, hv)
        b, cb = pull(y, vv)
        add(out, basis((a, b)), c * ca * cb)
    return out


def product_homotopy(pair):
    if degens(pair):
        return {}
    out = {}
    for (hv, vv), c in cm._h_terms(degree(pair)).items():
        a, ca = pull(pair[0], hv)
        b, cb = pull(pair[1], vv)
        add(out, basis((a, b)), c * ca * cb)
    return out


def rect_degenerate(x):
    return any(not s for s in x.signs) or any(not a and not b for a, b in zip(x.left, x.right))


def rect_basis(x):
    return {} if rect_degenerate(x) else {x: 1}


def h_boundary(x):
    p, out = len(x.signs), {}
    for i in range(p + 1) if p else ():
        vertices = tuple(j for j in range(p + 1) if j != i)
        y = h_pull(x, vertices)
        add(out, rect_basis(y), sign(i + sum(x.signs[:vertices[0]])))
    return out


def v_boundary(x):
    q, out = len(x.left), {}
    for i in range(q + 1) if q else ():
        y = v_pull(x, tuple(j for j in range(q + 1) if j != i))
        add(out, rect_basis(y), sign(i))
    return out


def rect_boundary(x):
    return add(h_boundary(x), v_boundary(x), sign(len(x.signs)))


@lru_cache(None)
def FD(x):
    if degens(x):
        return {}
    n, out = degree(x), {}
    for p in range(n + 1):
        add(out, rect_basis(rect_pull(x, tuple(range(p + 1)), tuple(range(p, n + 1)))))
    return out


@lru_cache(None)
def GD(x):
    if rect_degenerate(x):
        return {}
    out = {}
    for h, v, c in cm._shuffles(len(x.signs), len(x.left)):
        add(out, basis(rect_pull(x, h, v)), c * sign(sum(x.signs[:h[0]])))
    return out


@lru_cache(None)
def HD(x):
    if degens(x):
        return {}
    out = {}
    for (h, v), c in cm._h_terms(degree(x)).items():
        add(out, basis(rect_pull(x, h, v)), c * sign(sum(x.signs[:h[0]])))
    return out


def f0(x):
    if rect_degenerate(x):
        return {}
    if not x.left:
        return {(len(x.signs), 0, 0): 1}
    return {(len(x.signs), k, l): c for (_, k, l), c in
            r2.unit_p(r2.U1(0, x.left, x.right)).items()}


def g0(small):
    p, k, l = small
    if k + l == 0:
        return {Borel((1,) * p, (), ()): 1}
    return {Borel((1,) * p, a.edges, a.edges_prime): c
            for a, c in r2.unit_i((0, k, l)).items()}


def h0(x):
    return {Borel(x.signs, a.edges, a.edges_prime): sign(len(x.signs)) * c
            for a, c in r2.unit_h(r2.U1(0, x.left, x.right)).items()}


@lru_cache(None)
def FV(x):
    out, term = {}, {x: 1}
    for _ in range(len(x.signs) + 1):
        add(out, linear(term, f0))
        term = {a: -c for a, c in linear(linear(term, h0), h_boundary).items()}
        if not term:
            return out
    raise ArithmeticError('Borel projection series did not terminate')


@lru_cache(None)
def GV(small):
    out, term = {}, g0(small)
    for _ in range(small[0] + 1):
        add(out, term)
        term = {a: -c for a, c in linear(linear(term, h_boundary), h0).items()}
        if not term:
            return out
    raise ArithmeticError('Borel inclusion series did not terminate')


@lru_cache(None)
def HV(x):
    out, term = {}, {x: 1}
    for _ in range(len(x.signs) + 1):
        add(out, linear(term, h0))
        term = {a: -c for a, c in linear(linear(term, h0), h_boundary).items()}
        if not term:
            return out
    raise ArithmeticError('Borel homotopy series did not terminate')


@lru_cache(None)
def F(x):
    return linear(FD(x), FV)


@lru_cache(None)
def G(small):
    return linear(GV(small), GD)


@lru_cache(None)
def H(x):
    return add(dict(HD(x)), linear(linear(FD(x), HV), GD))


@lru_cache(None)
def small_boundary(small):
    return linear(linear(GV(small), rect_boundary), FV)


def word_degree(small):
    return sum(small)


omega_words, omega_boundary = r2.omega_words, r2.omega_boundary


def total_basis(n, relative=True):
    return tuple(((p, k, l), w) for p in range(n + 1) for k in (0, 1) for l in (0, 1)
                 if not relative or k == l == 1
                 for w in omega_words(n - p - k - l) if n >= p + k + l)


def total_boundary(b):
    small, omega = b
    out = {(a, omega): c for a, c in small_boundary(small).items()}
    for w, c in omega_boundary(omega).items():
        add(out, {(small, w): sign(word_degree(small)) * c})
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
    small, omega = b
    out = {}
    for x, c in G(small).items():
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
