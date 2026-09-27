# Shared cochain kernel: verification on 2026-09-27

Base revision: `f9eb8b6`, with the changes below uncommitted in the working
tree. The earlier [acceleration record](extension-acceleration-20260926.md)
is retained as historical evidence.

## Runtime changes

- One Python interval-cut engine. [cochain_tools.py](../../python/cochain_tools.py)
  holds the cut terms, the word evaluator with C-level face getters,
  memoized cochain calls without an intermediate Python frame, the
  coboundary, the interval pullback and right prism, and `Q`. The vendored
  [stacking_model/cochains.py](../../python/stacking_model/cochains.py)
  keeps its validating `Cochain` class but evaluates words and coboundaries
  through these functions; its copy of the engine was removed. The copy's
  negative-degree guard now applies to the shared cut terms: a word with a
  factor of negative degree has no terms.
- One set of shared helpers in [phase_eval.py](../../python/phase_eval.py):
  `polarization`, `hD`, `integral`, `chi` with precomputed face getters, and
  `source_splitting`, the degree-independent A-only splitting that
  `phase_eval.source` and `h_tau_primitive.uniform_source` computed
  separately. `compatible_sector` (Q, E, QD, polarization, alpha, integral),
  `lower_stacking` (interval pullback, prism), `h_tau_primitive` (hD,
  uniform source, interval, prism), `v3_pair_shared.hD`,
  `natural_upper.integral` and the lift in `upper_phase_diagnostic` call
  them.
- [extension_acceleration.py](../../python/extension_acceleration.py) no
  longer replaces the evaluator. It keeps structural zeros, identity
  memoization and the universal-value store, applied to the shared
  functions.
- Two source-hash checks were removed: the `stacking_model/provenance.json`
  manifest check in `extension_worker.py`, and the check of the recorded
  source hashes in `r3_source.ef_from_data`. The recorded hashes remain in
  those files as a record of the precompute. `high_calibration.json` is
  still checked against its sources, none of which changed.
- GAP reads each worker reply up to its newline
  (`KOAHSS_ExtensionReadWorkerLine`). `ReadLine` on a pty stream returns
  the characters that have arrived, and the worker's `print` writes a reply
  and its newline separately. On a loaded machine the newline could arrive
  late, and the next request then read it as an empty reply, failing with
  `syntax error at line 1 near:`. Both pty workers
  (`extension_transfer.gi`, `extension_bar.gi`) were affected on
  `f9eb8b6`.

Related functions that are not identical stay separate:
`phase_eval.source` adds the degree-one epsilon100 correction to the shared
splitting; `theta` and the pair phase in `theta_pair_phase.py` are distinct
operations; the stacking `differential` returns zero below degree zero; and
the GAP page engine is kept, since its Theta lift can differ from the Python
lift by an integral cochain. No calibration table or numerical payload
changed.

## Executed checks

GAP checks used fresh processes with `--quitonbreak` and an empty
`FERMIONAHSS_CACHE_DIR` unless stated. Most checks on the separate merged
worktree ran while other jobs loaded the machine (load average 20-60 on 32
threads), so only the side-by-side timings compare speed. "Old" is
`f9eb8b6`; the reference results are the recorded ones named below.

| Check | Observed result |
| --- | --- |
| `python3 -m unittest discover -s python -p 'test_*.py'`, old and new trees | Both: 33 tests OK, one skipped |
| Direct value comparison, old and new code: 25 nontrivial universal `source_value` keys, 25 `V1_pair` keys, and `Q`, `E`, `QD`, polarization, alpha, interval pullback with prism, signed differential and `hD` of the stacking modules on sample simplices | Identical |
| Worker replay of 54 degree-four and 7 degree-three requests recorded from the `29c045b` worker in `koFull(CyclicGroup(2),[1],[1],4)` | All answers identical to the recorded ones |
| `koFull(CyclicGroup(2),[1],[1],4)`, old and new trees side by side on an idle machine | `[[],[2],[],[2],[2],[16]]`; complete result summaries (relation matrices, Smith bases, lower coordinates, filtration stages, layer cochains, flat lifts, E6 table) identical to the `29c045b` baseline; empty cache 337.2 s old, 321.2 s new (peak RSS 2.44 and 2.35 GB); warm cache 36.0 s old, 33.8 s new |
| `LoadPackage("fermionAHSS"); Assert(0,TestPackage("fermionAHSS"));` | Passed, 206 s; the test files of `tst/smoke.tst` also passed one by one on the merged worktree |
| `examples/c2.g`, `examples/twisted_c2.g` | Output identical to the old tree |
| `examples/resolution_extensions.g` | Passed: C2 `[[0],[],[2,2],[2,2],[0,8]]`, C4 `[[],[2],[2],[2,4],[2]]`, C8 `[[0],[],[2,8],[2,2],[0,2,16]]`. On the loaded machine the old tree failed with `syntax error at line 1 near:`; run through a relay that forwards whole lines, it gave the same output as the new tree |
| C2 degree three through a relay that splits every worker line into three pieces 50 ms apart | Old tree: `syntax error at line 1 near:`; new tree: relations computed, invariants `[0,8]` |
| `examples/extension_papers.g` | 20 matches, 0 unresolved, 0 mismatches in 4 distinct calculations; output identical to the run of the earlier record |
| Eight space-group `koAHSS` cases at k=4 through E6 (SG008, SG031, SG207, SG208, SG214 with `w2+w1^2`; SG083, SG184, SG209 untwisted), run with `batch/koahss_worker.g` | Every event identical to `runs/papers-20260924-v3` apart from `runtime_ms`: page tables, T evaluations, secondary formula checks |
| SG008 `koAHSS`, old and new trees side by side | Identical events; 335.5 s wall and 328 s CPU old, 277.3 s wall and 271 s CPU new. GAP page times were unchanged; the saving is in the Python T-phase worker, which `runtime_ms` does not include |

No calibration table, numerical payload or recorded hash was changed.
`certified_ko` remains false and no new completeness or naturality claim is
made.
