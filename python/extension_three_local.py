"""Three-local stacking of the rows q=0 and q=-4 on the transferred resolution.

Away from the prime two the binary rows of the window vanish, and the
three-local ko spectrum in the window is the two-stage tower of the rows
q=0 (layer A) and q=-4 (layer D) with k-invariant 2 beta_3 P^1 rho_3. In
package degree five the input degree of A is two, where P^1 rho_3 A is the
reduction of the transported cube A cup A cup A: the rational potential
(2/3) lift(P^1_s rho_3 A) differs from Omega'(A) = (2/3) A cup A cup A by an
integral cochain, which is an isomorphism of stacking models. In the model
of Omega' the curvature of a cocycle A vanishes, the stacking correction is
the integral polynomial

    gamma(A, A') = -2 (A cup A cup A' + A cup A' cup A'),

whose difference from the polarization of Omega' is (2/3) delta_s Xi with
the cup-one expression Xi = 2 A(A cup_1 A') + (A cup_1 A')A + 2 (A cup_1 A')A'
+ A'(A cup_1 A'), and the boundary state of a gauge (u, w) is
(delta_s u, delta_s w). States are the four-layer records of the transferred
model with B = C = 0; every operation is evaluated on the resolution through
the comparison lift and projection of that model, so associativity and
commutativity hold up to signed coboundaries of the D layer exactly as in
the complete model.

The relation of a generator of order 3^a therefore reads, with 3^a A =
delta_s u and the local flat lift (A, 0),

    x_D = 2 * 3^(a-1) A cup A cup A + (2/3) delta_s [sum_j Xi(jA, A) - u cup delta_s u cup delta_s u]

modulo 3^a H_D; the last term is a Bockstein class and is computed, not
dropped. Only degree five is supported: in input degrees zero and one P^1
vanishes (the tower splits), and in input degree three the primitives of
the nineteen-term reduced power are not implemented.
Copyright (c) 2026 koAHSS contributors; MIT license.
"""
from extension_worker import api, BoundedCache
from low_phases import transported

p = api.p
FIELDS = "ABCD"
PRIME = 3


def degrees(k):
    return k - 3, k - 2, k - 1, k + 1


class ThreeLocalModel:
    """Two-layer stacking operations for three-primary classes of degree five."""

    def __init__(self, transferred):
        self.model = transferred
        self._gammas = BoundedCache(512)

    def supports(self, k):
        return k == 5

    def check_state(self, k, data):
        state = self.model.check_state(k, data)
        if any(state["B"]) or any(state["C"]):
            raise ValueError("three-local states have no binary layers")
        return state

    def zero(self, k):
        return self.model.zero(k)

    def _word(self, x, y, z):
        s = self.model.s
        return transported(transported(x, y, s), z, s)

    def gamma(self, k, a, b):
        """The D layer of the product of the flat states (a, 0) and (b, 0)."""
        if not self.supports(k):
            raise ValueError("the three-local stacking correction is implemented in degree five only")
        a, b = tuple(a), tuple(b)
        if not any(a) or not any(b):
            return [0] * self.model.dimension(k + 1)
        key = k, a, b
        if key not in self._gammas:
            n = degrees(k)[0]
            A = self.model.lift(n, a, True)
            Ap = self.model.lift(n, b, True)
            correction = p.scale(self._word(A, A, Ap) + self._word(A, Ap, Ap), -2)
            self._gammas[key] = tuple(self.model.project(correction, True))
        return list(self._gammas[key])

    def kappa(self, k, data, upto=3):
        """The curvature (delta_s A, 0, 0, delta_s D), stopping at the first nonzero layer."""
        state = self.check_state(k, data)
        answer = self.zero(k + 1)
        if not self.supports(k):
            raise ValueError("the three-local model covers degree five only")
        for layer in range(upto + 1):
            f = FIELDS[layer]
            if layer in (0, 3):
                answer[f] = self.model.coboundary(degrees(k)[layer], state[f], True)
            if any(answer[f]):
                break
        return answer

    def require_flat(self, k, data, upto=3):
        if upto >= 0 and any(any(v) for v in self.kappa(k, data, upto).values()):
            raise ValueError("three-local operation requires flat inputs through its requested layers")

    def product(self, k, left, right, upto=3):
        left, right = self.check_state(k, left), self.check_state(k, right)
        answer = self.zero(k)
        answer["A"] = [x + y for x, y in zip(left["A"], right["A"])]
        if upto == 3:
            answer["D"] = [x + y + g for x, y, g in
                           zip(left["D"], right["D"], self.gamma(k, left["A"], right["A"]))]
        return answer

    def act(self, k, gauge, canonical, upto=3):
        gauge = self.check_state(k - 1, gauge)
        canonical = self.check_state(k, canonical)
        boundary = self.zero(k)
        boundary["A"] = self.model.coboundary(degrees(k - 1)[0], gauge["A"], True)
        boundary["D"] = self.model.coboundary(degrees(k - 1)[3], gauge["D"], True)
        return self.product(k, boundary, canonical, upto)

    def divide_left(self, k, left, total, upto=3):
        left, total = self.check_state(k, left), self.check_state(k, total)
        right = self.zero(k)
        right["A"] = [t - x for x, t in zip(left["A"], total["A"])]
        if upto == 3:
            right["D"] = [t - x - g for x, t, g in
                          zip(left["D"], total["D"], self.gamma(k, left["A"], right["A"]))]
        return right

    def calculate(self, request):
        operation, k = request["operation"], request["degree"]
        if request.get("prime") != PRIME:
            raise ValueError("the transferred model localizes at the prime three only")
        if operation != "d" and not self.supports(k):
            raise ValueError("three-local stacking operations cover degree five only")
        upto = request.get("upto", 3)
        if type(upto) is not int or not 0 <= upto <= 3:
            raise ValueError("upto must be a layer index from zero through three")
        required = 3 if upto == 3 else upto - 1
        if operation == "d":
            if k not in (4, 5):
                raise ValueError("three-local curvature covers degree-five states and their gauges")
            if k == 4:
                # A gauge of a degree-five state: P^1 vanishes in input degree
                # one, so its boundary state is (delta_s u, 0, 0, delta_s w).
                state = self.check_state(k, request["state"])
                answer = self.zero(k + 1)
                for layer in (0, 3):
                    if layer <= upto:
                        answer[FIELDS[layer]] = self.model.coboundary(
                            degrees(k)[layer], state[FIELDS[layer]], True)
                return {"state": answer}
            return {"state": self.kappa(k, request["state"], upto)}
        if operation == "xtimes":
            self.require_flat(k, request["state"], required)
            self.require_flat(k, request["other"], required)
            return {"state": self.product(k, request["state"], request["other"], upto)}
        if operation == "act":
            self.require_flat(k, request["state"], required)
            return {"state": self.act(k, request["gauge"], request["state"], upto)}
        if operation == "divide_left":
            self.require_flat(k, request["state"], required)
            self.require_flat(k, request["other"], required)
            result = self.divide_left(k, request["state"], request["other"], upto)
            self.require_flat(k, result, required)
            product = self.product(k, request["state"], result, upto)
            if any(tuple(product[f]) != tuple(request["other"][f]) for f in FIELDS[:upto + 1]):
                raise ArithmeticError("three-local left division failed its exact product equality")
            return {"state": result}
        raise ValueError("unknown three-local extension transfer operation")
