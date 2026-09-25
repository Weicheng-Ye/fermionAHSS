"""Exact normalized finite-bar transport for the bundled four-cochain model.

Line-delimited JSON protocol. A process retains universal formula caches across
all lift, stacking and boundary calculations for one fixed group and twists.
Copyright (c) 2026 koAHSS contributors; MIT license.
"""
from functools import lru_cache
from collections import OrderedDict
import hashlib
import itertools
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "python" / "stacking_model"
# Reject an accidentally modified formula copy; calibration is not selected
# from a sibling checkout or from ambient Python paths.
manifest = json.loads((SOURCE / "provenance.json").read_text())
for name, metadata in manifest["files"].items():
    if hashlib.sha256((SOURCE / name).read_bytes()).hexdigest() != metadata["bundled_sha256"]:
        raise RuntimeError("stacking source checksum mismatch: " + name)
os.environ["FERMIONAHSS_ROOT"] = str(ROOT)
sys.path.insert(0, str(SOURCE))
import four_cochain_stacking as api


class BoundedCache(OrderedDict):
    def __init__(self, limit):
        super().__init__()
        self.limit = limit

    def __setitem__(self, key, value):
        super().__setitem__(key, value)
        self.move_to_end(key)
        while len(self) > self.limit:
            self.popitem(last=False)


class FiniteBar:
    def __init__(self, multiplication, sign, omega):
        self.mul = multiplication
        self.size = len(multiplication)
        if not self.size or any(len(row) != self.size for row in multiplication):
            raise ValueError("invalid multiplication table")
        self.inverse = [next(h for h in range(self.size) if multiplication[g][h] == 0)
                        for g in range(self.size)]
        self._cochains = BoundedCache(256)
        self._states = BoundedCache(128)
        self._answers = BoundedCache(256)
        self.s = self.cochain(1, sign)
        self.w = self.cochain(2, omega)
        self.rule = api.Stacking(self.s, self.w, self.is_zero)

    def dimension(self, n):
        return 0 if n < 0 else (self.size - 1) ** n

    @lru_cache(None)
    def simplices(self, n):
        if n < 0:
            return ()
        answer = []
        for increments in itertools.product(range(1, self.size), repeat=n):
            vertices = [0]
            for g in increments:
                vertices.append(self.mul[vertices[-1]][g])
            answer.append(tuple(vertices))
        return tuple(answer)

    def cochain(self, n, vector):
        if len(vector) != self.dimension(n) or any(type(x) is not int for x in vector):
            raise ValueError("invalid integral bar cochain")
        key = n, tuple(vector)
        if key in self._cochains:
            return self._cochains[key]
        vector = tuple(vector)
        def value(vertices):
            if n < 0:
                return 0
            index = 0
            for a, b in zip(vertices, vertices[1:]):
                g = self.mul[self.inverse[a]][b]
                if g == 0:
                    return 0
                index = index * (self.size - 1) + g - 1
            return vector[index]
        result = api.p.Cochain(n, value)
        self._cochains[key] = result
        return result

    def vector(self, cochain):
        answer = [cochain(sigma) for sigma in self.simplices(cochain.degree)]
        if any(x.denominator != 1 if hasattr(x, "denominator") else type(x) is not int
               for x in answer):
            raise ValueError("nonintegral output from the exact stacking model")
        return [int(x) for x in answer]

    def is_zero(self, cochain):
        return all(cochain(sigma) == 0 for sigma in self.simplices(cochain.degree))

    def state(self, k, data):
        cochains = tuple(self.cochain(n, data[f])
            for f, n in zip("ABCD", (k-3, k-2, k-1, k+1)))
        if any(x not in (0,1) for f in "BC" for x in data[f]):
            raise ValueError("B and C cochains must be binary")
        key = k, tuple(tuple(data[f]) for f in "ABCD")
        if key not in self._states:
            self._states[key] = api.State(*cochains)
        return self._states[key]

    def export(self, state):
        return {f: self.vector(getattr(state, f)) for f in "ABCD"}

    def calculate(self, request):
        k = request["degree"]
        if not -1 <= k <= 7:
            raise ValueError("this stacking bridge covers degrees -1 through 6 and the degree-7 endpoint")
        if k>=6 and not hasattr(self, "degree_six"):
            raise ValueError("degree-six finite section has not been configured")
        self.cochain(k+1, request["state"]["D"])
        if request["operation"] == "xtimes":
            self.cochain(k+1, request["other"]["D"])
        # D occurs only as delta_s D in d and D+Dprime in xtimes.
        # Share the expensive nonlinear A/B/C evaluation across all D carries.
        core = dict(request)
        core["state"] = dict(request["state"], D=[0] * self.dimension(k+1))
        if request["operation"] == "xtimes":
            core["other"] = dict(request["other"], D=[0] * self.dimension(k+1))
        key = json.dumps(core, sort_keys=True, separators=(",", ":"))
        if key not in self._answers:
            x = self.state(k, core["state"])
            if request["operation"] == "d":
                result = self.rule.d(x)
            elif request["operation"] == "xtimes":
                result = self.rule.xtimes(x, self.state(k, core["other"]))
            else:
                raise ValueError("unknown extension operation")
            self._answers[key] = self.export(result)
        answer = {f: list(values) for f, values in self._answers[key].items()}
        if request["operation"] == "d":
            carry = self.vector(api.p.ds(self.cochain(k+1, request["state"]["D"]), self.s))
        else:
            carry = [a+b for a,b in zip(request["state"]["D"], request["other"]["D"])]
        answer["D"] = [a+b for a,b in zip(answer["D"], carry)]
        return answer


def serve():
    model = None
    for line in sys.stdin:
        try:
            request = json.loads(line)
            if request["operation"] == "setup":
                model = FiniteBar(request["multiplication"], request["s"], request["omega"])
                if "degreeSix" in request:
                    from extension_degree_six import configure_degree_six
                    configure_degree_six(model, request["degreeSix"])
                answer = {"status": "computed"}
            elif model is None:
                raise ValueError("finite-bar model is not initialized")
            else:
                answer = {"status": "computed", "state": model.calculate(request)}
        except Exception as exc:
            status = "unresolved" if type(exc).__name__ == "DegreeSixResourceLimit" else "error"
            answer = {"status": status, "exception": type(exc).__name__, "reason": str(exc)}
        print(json.dumps(answer, separators=(",", ":")), flush=True)


if __name__ == "__main__":
    serve()
