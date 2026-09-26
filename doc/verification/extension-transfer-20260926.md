# Resolution-transfer prerequisites: verification on 2026-09-26

Base revision: `760baf1`. The working tree was clean before this revision.
This record covers the Smith cache, comparison preflight and model
contract tests described in [the plan review](../transfer.md). It does
not record a completed transferred extension calculation.

## Executed checks

All GAP tests below used fresh processes with `--quitonbreak`.

| Check | Observed result |
| --- | --- |
| `Test("tst/integer_equations.tst")` | Passed: 325 nonempty affine cases against the uncached Smith solver, inconsistent right-hand sides, empty dimensions, mutation isolation, cache reuse and eviction |
| `Test("tst/extension_transfer.tst")` | Passed: C2/C4 group-ring retractions through cochain degree 7, trivial/signed cochain comparisons, failed support/term budgets, short resolution, missing transport, and the rank-one contractible-pair counterexample |
| Extension lift, equivalence and gauge-reduction test files | Passed in one fresh integration process |
| `LoadPackage("fermionAHSS"); Assert(0,TestPackage("fermionAHSS"));` | Passed, process exit 0; elapsed wall time 57.036 s |
| `gap -q --quitonbreak examples/c2.g` | Passed, exit 0; E2–E6 tables printed |
| `gap -q --quitonbreak examples/twisted_c2.g` | Passed, exit 0; E2–E6 tables printed |
| Default Python transfer-contract tests | Six passed, one explicitly skipped; 6.422 s |
| Python transfer-contract tests with slow option | All seven passed; 51.813 s |
| `git diff --check` | Passed |
| New documentation links and bundled stacking source SHA-256 checks | Passed; all source checksums match `python/stacking_model/provenance.json` |

The Python commands were:

```sh
python3 -m unittest discover -s python -p test_extension_transfer_contracts.py -v
FERMIONAHSS_SLOW_TRANSFER_TESTS=1 python3 -m unittest discover -s python -p test_extension_transfer_contracts.py -v
```

These tests evaluate complete C2 bar cochains. They cover selected
square-zero cases in physical degrees 2–4, including nonclosed integral
A, exact signed integral carries and a legal pair with nonzero C
residual; multiplicativity cases in physical degrees 2 and 3; and a
literal two-sided zero unit in degree 3. The opt-in test evaluates the
expensive degree-four off-shell prism. These are finite cases, not a proof
of any universal identity.

An additional temporary probe completed 18 C2 square-zero cases in
physical degrees 2–4, for `(s,omega)=(0,0),(1,0),(1,1)`. A longer twisted,
nonflat degree-three multiplicativity probe was interrupted after several
minutes. It is recorded as incomplete, not passed or failed; the later
multiplicativity cases in that probe were not run. The permanent tests
above are the portable regression evidence.

The complete package test log ended with:

```text
#I  No errors detected while testing package fermionahss version 0.1.0
#I  using the test file `/Users/victor/.gap/pkg/fermionAHSS/tst/smoke.tst'
FERMIONAHSS_PACKAGE_REVIEW_EXIT=0; elapsed_seconds=57.036
```

## Integer-solve measurement

Construct a length-seven HAP resolution of C4, with `s=[1]`, `omega=0`,
then the existing complete-bar extension model in package degree 4.
Its signed `model.matrix(5,true)` has 243 rows and 729 columns. Using
the first matrix row as b, solve `x*M=b`, then `x*M=2*b` in the same
process. Both exact equations and equality of the returned homogeneous
bases were checked.

| Measurement | GAP CPU milliseconds |
| --- | ---: |
| Build the bar matrix | 2375 |
| First solve, including Smith preparation | 3558 |
| Second solve with cached preparation and a different RHS | 4 |

Matrix rank: 182. Homogeneous rank: 61. The first solve is not a timing
of the previous checkout; it measures preparation plus solving in the
new implementation. These numbers do not include nonlinear products,
Python worker time, or a full `koFull` computation.

## Additional sparse comparison checks

The preflight also passed on the resolutions below. Times are GAP CPU
milliseconds; k3, k4 and k5 reused one transport in that order. These
are g-support/retraction checks only, with no nonlinear stacking, face/H
closure or P1/P2 evaluation.

| HAP group constructor | k3 | k4 | k5 | Distinct degree-six g simplices |
| --- | ---: | ---: | ---: | ---: |
| `AbelianGroup(IsPermGroup,[2,2])` | 3 | 7 | 16 | 64 |
| `DihedralGroup(8)` | 31 | 144 | 767 | 680 |
| `DihedralGroup(IsPermGroup,8)` | 35 | 204 | 1246 | 957 |
| `SymmetricGroup(4)` | 184 | 2170 | Uncomputed | 5432 |

The S4 process was stopped after six seconds of wall time during its k5
audit. Its k3/k4 checks completed; no k5 result or refusal is inferred.

## Mathematical review and remaining acceptance work

Independent review checked the conditional flatness/reflection arguments,
the normalized-homotopy construction, the actual HAP contractible-pair
counterexample and the algebraic gauge-completeness countermodel. The
latter refutes sufficiency of the plan's restricted hypotheses; it does
not refute the calibrated formulas or prove transfer impossible.

Code review corrected two preflight-test issues before the passing runs:
unprefixed comments in GAP test output, and the need for one spare HAP
contraction degree. A separate reviewer checked the Smith cache's exact
solution convention, mutation isolation and empty-module behavior.

Still required before production transfer:

- A gauge-completeness/coherence proof for the actual formulas, or an
  independently justified filtered isomorphism argument.
- A comparison applicable to rejected resolutions, or an explicit
  narrower domain; raw H side conditions must not be inferred from
  successful retraction checks.
- Complete support closure and correct P1/P2 branch tests, including
  nonvacuous degree-five cases.
- The lazy nonlinear worker, separate gauge action, consumer hooks,
  explicit-resolution entry point, and per-degree model selection.
- Whole-bar verification of transferred product/gauge certificates and
  the proposed C2, Z4T and C2h classification comparisons.

No Z4 whole-bar nonlinear transfer comparison, degree-five specialized
test campaign, new T1–T3 direct-operation campaign, or full Z4T extension
timing was run. Production formulas, calibration JSON, source hashes,
historical verification records and batch snapshots were not changed.
The newly added record does not supersede the older extension evidence.
