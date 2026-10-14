# Status of this folder

`fermionAHSS_new` is the rewrite of `../fermionAHSS` (the twisted connective
real K-theory AHSS through E6 and the stacking extension problem). It is
complete and verified as of 2026-10-14, version 1.0.0. The user-facing
description is `README.md`; the mathematics is `docs/`; where each formula
is evaluated, the shortcuts and the caches are `docs/implementation.md`.
This file is the working status for whoever (person or agent) touches the
folder next.

## What was done

- Page differentials are components of the generalized stacking
  differential, evaluated by one page-mode Python worker
  (`gap/page_operations.gi`, `python/light_rows.py`): no separate secondary
  or tertiary operation files.
- The A-over-D relations of degrees five and six use the one-fiber
  primitives of `docs/power_relations.md` (`python/unary_tables.py`, tables
  in `data/`).
- Results of `koFull` are cached per degree, universal values and
  scaffolding per process and on disk; `FERMIONAHSS_CACHE_DIR` is the only
  environment variable (location of the cache; empty disables it).
- The Python formula modules (`python/*.py` except the worker modules listed
  in `_NOT_SOURCES` of `python/universal_values.py`, and
  `python/stacking_model/`) are byte-identical to the previous package.
  Their provenance hash (`8547700966fc5a5f…`) is recorded in
  `data/universal-values/manifest.json` and in the unary tables; any edit
  to a formula source invalidates them, and they must then be regenerated
  with `tools/generate_universal_values.py` and
  `tools/generate_unary_tables.py`. Worker caches can be folded into the
  shipped bundle with `tools/merge_universal_cache.py`.
- The GAP result cache is keyed by a hash of `gap/`, `python/`,
  `python/stacking_model/` and the formula data, so stored results are
  ignored automatically after any edit.

## Verification

`bench/` is the comparison harness against the previous package and is
not part of the package proper: it is kept in this working copy only
(ignored by git), nothing in `gap/`, `python/` or `tst/` uses it, and its
Slurm scripts hardcode this account's paths and the LPTMC cluster
(partition `normalx`). `bench/README.md` lists the drivers;
`bench/write_report.py` regenerates `bench/REPORT.md` from result
directories (`bench/results/` is ignored by git).

`bench/REPORT.md` is the record. Both packages ran alternately on the same
node in fresh processes with cold caches, and every result (E6 tables,
invariants, relation matrices, filtrations, relation rows) is identical:

| comparison | cases |
| --- | --- |
| finite groups, degrees 3 to 6, cold / warm / warm-compute passes | 33 |
| 3d space groups at degree 4, orientation twist and pin twist | 38 |
| page differentials and degree results of 13 finite groups, every twist of a basis, degrees 1 to 5 | 345 |
| page differentials of the same groups through degree 4 against the stage-one snapshot of the port (`bench/results/stage1-snapshot`) | 276 |

Cold wall time: 141 s against 596 s on the finite-group cases, 3398 s
against 4810 s on the space groups. The only slower case is space group
230 with the pin twist (93 s against 79 s): the new model's B-layer source
needs lifts on more faces of the comparison chains than the old secondary
operation did. Trivial cases pay about 0.3 s for starting the page worker.
`tst/run.g` (smoke tests) passes.

## Things learned the hard way

- Four exact shortcuts in the transport were needed to reach parity on
  space groups (encoded transport answers cached by anchored labels, the
  structural zero of an even lift, sign characters of affine matrices from
  generator characters, a contraction memo of 16384 entries); see
  `docs/implementation.md`. When a space-group case is slow, count the
  worker's transport callbacks per kind before touching any formula: the
  cost has always been glue (contractions of the HAP resolution, encoding),
  never the formulas.
- Never edit the package while a comparison campaign runs: the result
  cache key changes and the processes read the files as they start.
- The previous package is `../fermionAHSS` at commit cf1e2f3 and must stay
  unmodified; it is the reference for every comparison.

## Open items

- At degree six the stacking correction of two states that both have a
  nonzero A layer is refused (`A_STACKING` in `python/native_upper.py`: it
  would evaluate the universal pair source of `production_gamma6`, whose
  terms nest the $V_3$ contraction), so a relation that needs it is
  reported unresolved; the previous package had the same limit. Squares of
  A generators do not need it (one-fiber primitives).
- Space group 230 with the pin twist is 17% slower than before; closing
  that gap means a different choice of defining cochains in the model, not
  more caching.
