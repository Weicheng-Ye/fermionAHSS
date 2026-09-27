# Bundled universal values and memo speed-ups: verification on 2026-09-27

Base revision: `34dcfcc`, with the changes below uncommitted in the working
tree. The earlier [worker speed-up record](worker-speedups-20260927.md) is
retained as historical evidence.

## Runtime changes

- [data/universal-values.json](../../data/universal-values.json) ships the
  1,129 pair-source and 818 `V1_pair` values that the C2 and Z4 examples
  need, with the hash of the formula sources that produced them. The
  `koFull` worker loads them before its cache when the hash matches the
  current sources and never writes them to the cache file;
  `FERMIONAHSS_BUNDLED_VALUES=0` ignores them.
  [generate_universal_values.py](../../python/generate_universal_values.py)
  recomputes every listed value with the current formulas (and the keys of
  any cache file given with `--add`), and `test_universal_values.py` checks
  the hash and a sample of values.
- The page T worker ends quietly when GAP closes its pty, instead of printing
  an `OSError` traceback into every `koAHSS` log. In serve mode it closes the
  descriptors inherited from GAP: a profiler's gzip pipe held open by the
  worker made `UnprofileLineByLine` wait forever.
- `KOAHSS_BoundedMemo` replaces GAP dictionaries in the natural-cochain and
  chi memo tables. It keeps the same 256 entries and forgets the oldest one
  first, but keeps its keys sorted: GAP's `RemoveDictionary` scanned the
  whole dictionary on each eviction, comparing lists of group elements. The
  per-simplex degeneracy and vertex checks use plain loops.
- `cochain_tools.LinearCochain`: a sum, difference or binary reduction of a
  sum extends one term list instead of nesting one memoized closure per
  operation. Terms are evaluated in the same order with the same exact
  arithmetic.

No formula, calibration table or cache limit changed.

## Measurement behind the GAP change

A line-by-line profile (`ProfileLineByLine`, wall time) of the SG008
`koAHSS` case attributed about 27% of GAP time to `dict.gi` line 113, the
linear `PositionProperty` of `RemoveDictionary`, reached from the memo
evictions; `ForAny`/`ForAll` closures took about 10% and group-element
matrix products about 7%. These are shares of a profiled run, not times.

## Executed checks

GAP checks used fresh processes with `--quitonbreak` and an empty
`FERMIONAHSS_CACHE_DIR`. "Old" is a `git archive` copy of `34dcfcc`;
side-by-side timings ran both trees at the same time on an otherwise idle
AMD EPYC 8124P. "Cold" runs set `FERMIONAHSS_BUNDLED_VALUES=0`.

| Check | Observed result |
| --- | --- |
| Flattened sums against the old kernel: 400 random expressions of sums, differences and binary reductions, 10,000 values | Identical values and value types |
| `KOAHSS_BoundedMemo` against the previous dictionary with its FIFO ring: 18,000 random lookups with limits 1, 7 and 256 | Identical results; 5.5 µs instead of 8.7 µs per lookup with evictions (keys of three 4×4 rational matrices) |
| `generate_universal_values.py` on the keys of the C2, Z4 and example caches | 1,947 values in 231 s, each equal to the value computed by the old kernel |
| Replay of the 54 recorded degree-four worker requests, cold, collector off | All answers identical; 161.7 s old, 160.0 s new; peak RSS 2.35 GB old, 1.84 GB new |
| One SG008 page T batch | Identical reply; 40.1 s old, 35.5 s new |
| Start the page worker, then `UnprofileLineByLine()` and quit | New tree exits after 10 s; the old tree waits until killed |
| `koFull(CyclicGroup(2),[1],[1],4)`, side by side | Complete result summaries identical to the `29c045b` baseline in all five runs; empty cache 164.1 s old, 156.2 s new cold, 24.0 s new with the bundled values; warm 24.8 s old, 23.6 s new; peak RSS 2.35, 1.84 and 0.42 GB, warm 0.53 and 0.42 GB |
| `koFull(Group((1,2,3,4)),[1],0,4)`, side by side | `[[],[2],[2],[2,4],[2],[4]]`; complete summaries identical in all three runs; 131.2 s old, 126.1 s new cold, 77.7 s with the bundled values; peak RSS 2.59, 2.04 and 1.49 GB |
| `examples/extension_papers.g`, empty cache, side by side | 20 matches, 0 unresolved, 0 mismatches for both trees; output identical to the earlier runs; 276.2 s old, 63.8 s new |
| SG008 with `w2+w1^2`, `koAHSS` k=4 through E6, side by side | Events identical between the trees and to `runs/papers-20260924-v3` apart from `runtime_ms`; wall 109.2 s old, 91.0 s new; GAP page time 32.9 s old, 20.8 s new; the worker traceback appears only in the old log |
| SG083 with `w2+w1^2`, same | Identical events, also to the campaign; wall 954.3 s old, 574.0 s new; GAP page time 681.0 s old, 329.2 s new (E5 209.5 and 145.8 s, E6 452.0 and 169.9 s) |
| Eight space-group `koAHSS` cases (SG008, SG031, SG207, SG208, SG214 with `w2+w1^2`; SG083, SG184, SG209 untwisted) | Every event identical to the campaign apart from `runtime_ms`; no tracebacks in the logs; no page worker left running |
| `LoadPackage("fermionAHSS"); Assert(0,TestPackage("fermionAHSS"));` | Passed |
| `examples/c2.g`, `examples/twisted_c2.g`, `examples/resolution_extensions.g` | Output identical to the pre-change runs |
| `python3 -m unittest discover -s python -p 'test_*.py'` | 36 tests passed, including the 3 new ones, one skipped as before; 292 s |

`certified_ko` remains false and no new completeness or naturality claim is
made.
