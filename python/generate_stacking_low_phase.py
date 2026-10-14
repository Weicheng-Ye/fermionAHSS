"""Regenerate the finite, exact closed-input degree-two stacking phase.

The production runtime uses the checked-in payload, not this source workspace.
Pass the sibling fermionAHSS_stacking directory explicitly. The table contains
local cochain values, never expected classification groups. All 256 binary
triangle configurations are evaluated in the selected xtimes normalization.

Copyright (c) 2026 koAHSS contributors. Distributed under the MIT license.
"""
import argparse
import hashlib
import importlib
import itertools
import json
from pathlib import Path
import sys


SOURCES = (
    "cochains.py", "compatible_sector.py", "lower_stacking.py",
    "paper_stacking.py", "a0_gamma.py", "a0_degree2.py",
    "coherent_low_commutative.py", "paper_commutative_production.py",
    "four_cochain_stacking.py",
)


def generate(source):
    source = Path(source).resolve()
    sys.path.insert(0, str(source))
    cochains = importlib.import_module("cochains")
    selected = importlib.import_module("coherent_low_commutative")
    rule = selected.DegreeTwoCommutativeStacking()
    values = []

    def closed_one(x, y):
        potentials = (0, x, (x + y) % 2)
        return cochains.Cochain(1, lambda face:
            (potentials[face[0]] + potentials[face[1]]) % 2)

    for b, bp, x, y, xp, yp, s, t in itertools.product(range(2), repeat=8):
        phase = rule.phase(
            cochains.Cochain(0, lambda face, value=b: value), closed_one(x, y),
            cochains.Cochain(0, lambda face, value=bp: value), closed_one(xp, yp),
            closed_one(s, t), cochains.zero(2))((0, 1, 2))
        values.append([phase.numerator, phase.denominator])
    canonical = json.dumps(values, separators=(",", ":")).encode()
    return {
        "schema": 1,
        "calibration": "coherent-low-commutative/closed-B-C/omega-zero/v1",
        "domain": "B,Bprime are closed degree-zero; C,Cprime,s are binary one-cocycles; omega=0",
        "bit_order": ["B", "Bprime", "C01", "C12", "Cprime01", "Cprime12", "s01", "s12"],
        "values": values,
        "values_sha256": hashlib.sha256(canonical).hexdigest(),
        "source_sha256": {name: hashlib.sha256((source / name).read_bytes()).hexdigest()
                          for name in SOURCES},
        "copyright": "Copyright (c) 2026 koAHSS contributors; MIT license",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = json.dumps(generate(args.source), indent=2) + "\n"
    if args.output:
        args.output.write_text(payload)
    else:
        print(payload, end="")
