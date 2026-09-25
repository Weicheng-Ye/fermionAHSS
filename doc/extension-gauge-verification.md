# Staged gauge reduction: verification on 2026-09-26

The extension reducer now searches preceding-degree gauges in the order
`(0,0,0,d)`, `(0,0,c,d)`, `(0,b,c,d)`, `(a,b,c,d)`. Earlier components are
literally zero at each stage. It solves the remaining defining equations
together and accepts a result only after verifying

\[
\mathrm d g\mathbin{\times}h=x
\]

in all four cochain components. Here `x` is the original stacked state and
`h` is a product of the original stored lower-layer lifts. This retains
both the common marked basis and the lower carries introduced while
removing earlier components.

The E6 projection of a residual in D proposes its marked coordinates.
That projection alone is not an extension witness: the complete equation
above must still be solved. When projected coordinates are unavailable,
the reducer tries at most 32 finite marked lower normal forms. Its bounded
searches can still return unresolved; they do not prove nonequivalence.
See [the extension contract](extensions.md#complete-flat-representatives-and-boundary-comparisons).

## Previously failing C2h relation

The original point-group run left the degree-four C relation unresolved
for C2h with `s=w1`, `omega=w2+w1^2`. Its native HAP twist coordinates were
`s=[1,0]`, `omega=[1,1,0]`, with seed 1 and resolution length 7. The
abstract group is C2 × C2. This package degree corresponds to spatial
dimension 3.

The regression reuses that run's fixed full C and D lifts, reconstructs
the same complete normalized bar model, verifies both lifts are flat, and
recomputes `C xtimes C`. The ordinary integer solve using the D lift and
ordinary signed coboundaries fails, reproducing the reported limitation.
The new reducer then finds

\[
C\mathbin{\times}C=\mathrm d(0,0,c,d)\mathbin{\times}D.
\]

It checks both lower coordinates 0 and 1. Stage D fails for both; stage
CD exhausts all 16 leading C choices for coordinate 0, then proves
coordinate 1. The winning gauge has four nonzero C entries and 54 nonzero
D entries. The exact comparison, flatness checks, rejected attempts,
original twist certificate and source hashes are stored in
[extension-gauge-results-20260926.json](../data/extension-gauge-results-20260926.json).
Its `runtimeMs` is GAP CPU time for reduction, not wall time or Python
worker CPU time.

A separate fresh E6 calculation, using one shared resolution and backend
for page computation and projection, sends the stacked D cochain to
native coordinates `[23,0,0,13,0,-6]`, E2 cohomology coordinates `[1,1,0]`,
and E6 D coordinates `[1]`. This independently checks the coordinate
proposal used by the production fallback.

This verifies the saved C relation. An optional fresh complete C2h
calculation was stopped after approximately 15 minutes, while the Python
stacking worker was still computing. It produced no full classification;
that check is recorded as incomplete, not passed or mathematically
unresolved. Both test processes were terminated. The focused relation and
E6 projection checks above completed successfully.

The portable replay is:

```sh
gap -q --quitonbreak examples/c2h_gauge.g
```

This manually run example recomputes the stacking product and proves the
saved proposed coordinate 1. It does not use the historical run directory
or insert an expected group into the solver.
The example passed a GAP parse-only check; the equivalent saved-lift
computation described above was executed separately.

## Bounded regression tests

Fresh GAP checks passed for:

- Each possible winning stage: D, CD, BCD and ABCD.
- Preferred cohomology directions followed by an unpreferred ordinary
  kernel direction needed for a successful comparison.
- An earlier fixed C contribution with a nontrivial D carry, ensuring the
  final comparison uses the original target and the same full C lift.
- A wrong proposed D coordinate, finite enumeration limits, and a free
  lower generator without projected coordinates remaining unresolved.
- Preserved diagnostics for unresolved relation queries.
- The full package suite, including the existing C2 full-state extension
  and display/API checks.
- `examples/c2.g` and `examples/twisted_c2.g`, each in a fresh GAP process.

The files exercising the new behavior are
[extension_equivalence.tst](../tst/extension_equivalence.tst) and
[extension_gauge_reduction.tst](../tst/extension_gauge_reduction.tst).

## Existing batch runs

Existing point-group runs retain their frozen source snapshots. Editing
the working package does not change their queued jobs, prior results or
`summary.csv`. A fresh run directory is required to use this fix; see
[the point-group runner](../batch/POINT_GROUPS.md). No original run artifact
was overwritten by these checks.

All results remain within the bounded five-row stacking model, with
`certified_ko=false`.
