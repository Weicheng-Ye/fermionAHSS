# Extension acceleration: verification on 2026-09-26

Base revision: `29c045b`, with the changes below uncommitted in the working
tree. The earlier [native-default record](native-extension-default-20260926.md)
is retained as historical evidence.

## Runtime changes

- The finite quotient audit was removed: `KOAHSS_ExtensionFiniteAudit`,
  `gap/extension_audit.gi` and `tst/extension_audit.tst` no longer exist,
  and there is no optional fallback. Commutativity and associativity of
  stacking on gauge classes are assumed; completed native degree results
  record `abelianQuotientAssumed=true`. Its non-audit assertions (C2 degree
  three is `[0,8]`, with computed canonical comparisons) moved to
  `tst/extensions.tst`.
- The worker's left division stores its final product as the core of the
  verifying product instead of evaluating the same product twice.
- `reflect` evaluates its D terms on the cached embedding of the native
  A, B, C coordinates, so their formulas are shared with other uses of the
  same state.
- `python/extension_acceleration.py` installs an exact evaluation policy in
  the extension worker only: C-level interval-cut face getters and cochain
  calls, canonical structural zeros, identity memoization of pure builders
  (32 recent entries each), and persistent universal source values
  (`production_gamma4.source_value`, `low_phases.V1_pair`) in
  `FERMIONAHSS_CACHE_DIR` (default `$XDG_CACHE_HOME/fermionAHSS` or
  `~/.cache/fermionAHSS`), keyed by the hashes of all formula sources.

No formula source, calibration table or source hash changed. The 58
`stacking_model` hashes, `high_calibration.constants()` and
`r3_source.ef_from_data` all validate against the unchanged files.

## Executed checks

All GAP checks used fresh processes with `--quitonbreak`. Timings are wall
times on the same machine (AMD EPYC 8124P); the old tree is a `git archive`
copy of `29c045b`.

| Check | Observed result |
| --- | --- |
| `koFull(CyclicGroup(2),[1],[1],4)`, old tree | `[[],[2],[],[2],[2],[16]]`; 1881.2 s |
| Same call, new tree, empty cache directory | Identical result; 324.6 s; peak RSS 2.44 GB |
| Same call, new tree, cache from the previous run | Identical result; 34.7 s; peak RSS 0.56 GB |
| Result comparison, degrees -1 through 4 | Identical invariants, relation matrices, Smith bases, relation vectors and lower coordinates, filtration stages, layer cochains, flat lifts and E6 table (cold and warm runs) |
| Worker replay of the 54 recorded degree-four requests before the audit | All 54 answers identical to the old worker; 327.4 s with an empty cache, 32.9 s warm, against about 709 s |
| `python3 -m unittest discover -s python -p 'test_extension*.py'` | 29 tests, one slow-contract skip as before; 269.9 s (old tree: 361.6 s) |
| `LoadPackage("fermionAHSS"); Assert(0,TestPackage("fermionAHSS"));` | Passed, exit 0 |
| Test files of `tst/smoke.tst`, timed one by one | All passed; 253 s in total (old tree with the audit file: 283 s) |
| `examples/c2.g`, `examples/twisted_c2.g` | Passed; unchanged page output |
| `examples/resolution_extensions.g` | Passed: C2, signed C4 and C8 through package degree 3 |
| `examples/extension_papers.g` | 20 matches, 0 unresolved, 0 mismatches in 4 distinct calculations; 335.9 s from an empty cache |

The cache directory was a fresh scratch directory for every check above;
the default user cache was not written. `certified_ko` remains false and no
new completeness or naturality claim is made.
