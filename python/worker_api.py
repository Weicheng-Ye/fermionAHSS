"""Exact normalized finite-bar transport for the bundled four-cochain model.

Line-delimited JSON protocol. A process retains universal formula caches across
all lift, stacking and boundary calculations for one fixed group and twists.
Copyright (c) 2026 koAHSS contributors; MIT license.
"""
from functools import lru_cache
from collections import OrderedDict
import gc
import itertools
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "python" / "stacking_model"
# The formulas are always loaded from this package, never from a sibling
# checkout or ambient Python paths. Their bundled hashes are not verified.
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
