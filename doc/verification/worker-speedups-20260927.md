# Worker and word-evaluator speed-ups: verification on 2026-09-27

Base revision: `6e69ab5`, with the changes below uncommitted in the working
tree. The earlier [shared-kernel record](shared-kernel-20260927.md) is
retained as historical evidence.

## Runtime changes

- The Python workers (`worker.py`, `extension_transfer.py`,
  `extension_worker.py`) call `gc.disable()` in their `__main__`. Their heaps
  are long-lived memo tables with almost no reference cycles, so cyclic
  collection only rescanned them. Importing the modules leaves collection on.
- GAP starts the page T worker once per session (`worker.py --serve`) and
  sends every batch as one JSON line with an increasing id, which the reply
  echoes; a reply to an earlier, interrupted request is skipped. Previously
  each batch started a new process, which repeated the one-time set-up behind
  V1 (the wedge normalization and the odd-primitive values) and discarded its
  universal values.
- `KOAHSS_NaturalWordEvaluator(word,degrees,inputs[,cupIndex])` checks an
  interval-cut word and collects its cuts once, and returns the function
  evaluated on each simplex. `ctx.cup`, `ctx.zeta` and the primary scalar
  expression build their evaluators when the cochain is built;
  `koAHSSNaturalWordValue` and `koAHSSNaturalZetaValue` are wrappers with
  unchanged results. Previously every simplex evaluation checked the word
  again and rebuilt a string key for the pattern cache.
- The newline-safe reply reader is now `KOAHSS_ReadWorkerLine` in
  `natural_tertiary.gi`, formerly `KOAHSS_ExtensionReadWorkerLine` in
  `extension_bar.gi`, since the page worker uses it too.

No formula, calibration table, cache limit or numerical payload changed.

## Executed checks

GAP checks used fresh processes with `--quitonbreak` and an empty
`FERMIONAHSS_CACHE_DIR` unless stated. "Old" is a `git archive` copy of
`6e69ab5`; side-by-side timings ran both trees at the same time on an
otherwise idle AMD EPYC 8124P.

| Check | Observed result |
| --- | --- |
| Old and new word evaluators on random input: 4,000 binary words, 2,500 integral cup-i words (each also in binary mode), 800 zeta words; 462 of the shapes uncached, 1,028 zero shortcuts | The old function, the new wrapper and the new evaluator agree in every case |
| Page worker line protocol: ids, a malformed and an unsupported request, single-shot mode | Ids echoed, errors returned as replies, process keeps serving; single-shot mode unchanged |
| `koFull(CyclicGroup(2),[1],[1],4)`, old and new side by side | `[[],[2],[],[2],[2],[16]]`; complete result summaries identical to the `29c045b` baseline for both trees, cold and warm; empty cache 307.1 s old, 171.3 s new; warm 33.8 s old, 25.2 s new; peak RSS 2.35 GB for both |
| SG008 with `w2+w1^2`, `koAHSS` k=4 through E6 with `batch/koahss_worker.g`, side by side | Every event identical between the trees and to `runs/papers-20260924-v3` apart from `runtime_ms`; wall 159.0 s old, 109.0 s new; GAP page time 40.3 s old, 33.0 s new |
| SG083 with `w2+w1^2`, same | Identical events, also to the campaign; wall 1089.3 s old, 958.5 s new; GAP page time 802.1 s old, 683.1 s new |
| `LoadPackage("fermionAHSS"); Assert(0,TestPackage("fermionAHSS"));` | Passed |
| `examples/c2.g`, `examples/twisted_c2.g`, `examples/resolution_extensions.g` | Output identical to the pre-change runs |
| `examples/extension_papers.g` | 20 matches, 0 unresolved, 0 mismatches; output identical to the earlier run |
| Eight space-group `koAHSS` cases (SG008, SG031, SG207, SG208, SG214 with `w2+w1^2`; SG083, SG184, SG209 untwisted) | Every event identical to `runs/papers-20260924-v3` apart from `runtime_ms` |
| `python3 -m unittest discover -s python -p 'test_extension*.py'` | 29 tests passed, one skipped as before |
| Page workers after their GAP sessions ended | None left running |

`certified_ko` remains false and no new completeness or naturality claim is
made.
