# fermionAHSS

An exact GAP package for the five rows `q = -4,...,0` of the twisted
connective real K-theory Atiyah–Hirzebruch spectral sequence through E6 and
for the stacking extension problem that assembles E6 into the classification
group of fermionic symmetry-protected topological phases, for finite groups
and for 3d space groups.

This is the rewrite of the previous `fermionAHSS` package. It computes the
same groups and presentations with one model: every page differential and
every relation of the extension problem is a component of the generalized
differential or of the product of the four-cochain stacking model
(`docs/differential.md`, `docs/product_*.md`). There are no separate
secondary or tertiary cohomology operations, no configuration switches, and
every expensive intermediate object is cached on disk.

## Requirements

GAP 4.12 or newer with the packages polycyclic, HAP and json; Python 3.10 or
newer on `PATH` as `python3` (standard library only). Space groups need the
local SpaceGroupCohomology package (`batch/koahss_groups.g`).

## Use

```gap
Read("/path/to/fermionAHSS_new/load.g");
full := koFull(CyclicGroup(4), [1], 0, 6);;   # group, s, omega, k
full.invariants;                              # [ 4 ]
koAHSSDisplay(full);                          # the E6 line and the group
```

- `koAHSS_batch(G or R, s, omega, k[, n][, rec(details:=true)])`: the pages
  of every degree through `k` (`n` pages from E2; one E6 table without `n`).
- `koAHSS(...)`: the line `p+q = k-3` of those pages.
- `koFull(G or R, s, omega, k)` or `koFull(detailedE6Result)`: the group of
  degree `k` (`invariants`, `degreeResult` with the presentation).
- `koFull_batch(...)`: every degree from -1 through `k`.

`s` and `omega` are binary cocycle vectors in the resolution's basis (scalar
`0` for the zero class). `k` is the package degree `p+q+3` in `[-1..6]`; the
physical spatial dimension is `k-1`. A supplied HAP resolution (including one
of an infinite group) is used as is and must have length at least `k+3`.
`batch/koahss_groups.g` builds a space or point group with its resolution,
and `batch/koahss_spacegroup_twists.g` its orientation and pin twists.

## What is computed, and how

1. E2 is the twisted cohomology of the resolution. The primary arrows
   (`Dbar`, `D`, `Dtilde`) are cup-i words on the resolution itself
   (`gap/hap_space.gi`, `gap/binary_tensors.gi`).
2. The higher arrows `Tau`, `Psi`, `T` are the C and D components of the
   generalized differential on the defining systems `(A,B,0,0)`, `(0,b,c,0)`
   and `(A,B,C,0)`: the defining cochains are solved on the resolution, the
   page worker (`python/stacking_worker.py`, `python/light_rows.py`) forms
   the bar cochains `Lambda r + H z` with the group-bar comparison, evaluates
   the components and pairs them back along the comparison chains
   (`gap/page_operations.gi`). The page engine (`gap/pages.gi`) skips every
   arrow with a trivial source or target; `Psi` uses the closed form of the
   component at `A=0`, and `T_0` the native Pontryagin-square Bockstein.
3. The extension problem (`gap/extension.gi`, `gap/relations.gi`,
   `gap/light_rows.gi`) measures the relation of each E6 generator from the
   same worker: target-layer shortcuts, primary-operation rows, light rows
   (residues of transported defining data) and the three-local rows, with
   the squares of A generators evaluated through the one-fiber primitives of
   `docs/power_relations.md` in degrees five and six (`python/unary_tables.py`).

The formulas are those of `docs/`; the Python formula modules in
`python/` and `python/stacking_model/` are shared with the previous package
byte for byte, so the shipped universal-value bundle and coefficient tables
(`data/`) remain valid (their provenance hash is checked on load).

## Caches

Everything persistent lives in `cache/` (or in `FERMIONAHSS_CACHE_DIR` when
set; an empty value disables persistence):

- universal values (`universal-values-<hash>/`), appended by every worker,
  on top of the shipped bundle `data/universal-values/` (one file per table,
  read when the table is first used);
- the scaffolding of the universal models (`scaffolding-<hash>/`), the
  certificates and Smith reductions built once per formula version;
- degree results (`result-<key>.json`): a computed degree is returned from
  the cache when the package sources, the resolution, the twists and the
  degree are the same (`gap/result_cache.gi`).

Delete the directory to recompute everything.

## Layout

| Directory | Contents |
| --- | --- |
| `gap/` | the GAP package: `api.gi` (entry points), `pages.gi` (page engine), `cochain_space.gi` and `hap_space.gi` (cohomology and primary operations on a resolution), `page_operations.gi` (Tau, Psi, T from the differential), `bar_transport.gi` and `cell_transport.gi` (comparison with the bar), `worker_model.gi` (the Python worker), `primary_transport.gi`, `linear_algebra.gi`, `flat_lifts.gi`, `gauge_compare.gi`, `relations.gi`, `light_rows.gi`, `extension.gi`, `result_cache.gi`, `display.gi` |
| `python/` | the exact cochain engine and the stacking formulas (`stacking_model/`), the worker (`stacking_worker.py`, `light_rows.py`, `acceleration.py`), the universal-value store and tables |
| `data/` | the chi ANF table, the universal-value bundle, the one-fiber coefficient tables, the paper fixtures |
| `docs/` | the mathematics: every formula the package evaluates |
| `tools/` | generators of the shipped tables |
| `batch/` | group and twist constructors for space and point groups |
| `tst/`, `examples/` | tests and examples |

## Tests

```sh
gap -q --quitonbreak -T tst/run.g
```
