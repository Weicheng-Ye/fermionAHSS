# Native-resolution extensions: verification on 2026-09-26

Base revision: `23370ce`. This implementation follows the initial
[mathematical review](../transfer.md) and its separate
[prerequisite record](extension-transfer-20260926.md). The runtime API
and its completion criteria are described in
[resolution extensions](../resolution-extensions.md).

## Executed checks

All GAP checks used fresh processes with `--quitonbreak`.

| Check | Observed result |
| --- | --- |
| `tst/extension_transfer.tst` | Passed: exact retraction, resource refusals, normalized SDR identities, native worker integration, full-bar C3/C4 reflection, and actual explicit-resolution fallback |
| `tst/extension_transfer_certify.tst` | Passed: independent reference relations and finite audit; missing export and reference resource refusal stay unresolved; an export changing the marked leading representative raises an error |
| `tst/extension_transfer_hooks.tst` and `tst/api.tst` | Passed: supplied-resolution identity, model options, action/division hooks, native/refusal fallback, and preservation of shared reference-model caching |
| `LoadPackage("fermionAHSS"); Assert(0,TestPackage("fermionAHSS"));` | Passed, exit 0; 68.067 s, including the new transfer tests |
| `examples/c2.g` and `examples/twisted_c2.g` | Passed: E2–E6 tables |
| `examples/resolution_extensions.g` | Passed: native/reference equality through package degree 3 for C2 and signed C4, without fallback |
| C2 paper fixture calculations with transfer selected | All four calculations completed; 20 labeled matches, zero mismatches and zero unresolved results; runner exit caveat below |
| Python `test_extension_transfer*.py` | 21 tests: 20 passed, one explicit slow-contract skip; 34.379 s |
| Bundled formula source checks | All 58 SHA-256 hashes match `python/stacking_model/provenance.json` |
| Diff and documentation checks | `git diff --check` and local links in changed documentation passed |

The Python command was:

```sh
python3 -m unittest discover -s python -p 'test_extension_transfer*.py' -v
```

The skipped contract test is the previously executed expensive
degree-four off-shell prism case; its passing run belongs to the
prerequisite record. It was not rerun as part of this command. The 14 new
worker regressions include signed degree-four action on a nonclosed
gauge, a degree-five pure-C square compared with the full bar, exact
D-affine cache carries for products/actions/division, native twist
validation, curvature prefix semantics, and explicit global-flag limits.

## Nonidentity comparison and full-result checks

The normalized homotopy is checked over the integral group ring on all
40 normalized C4 simplices in degrees 0–3. The test verifies
`boundary*h + h*boundary = 1-g*f`, `f*h=0`, `h*h=0`, and `h*g=0`;
some tested homotopy chains are nonzero. Signed cochain behavior is
tested separately. These bounded tests supplement the algebraic
normalization argument; they do not prove all-degree identities.

The C3 native state `(A,B,C,D)=(0,0,1,0)` has native square `(0,0,0,1)`.
The test exports the full corrected embedding, checks its complete-bar
flatness, and verifies literally
`Phi(x) xtimes Phi(x) = d(reflectionGauge) xtimes Phi(nativeProduct)`.
It also checks exact native division. An additional temporary probe
verified the corresponding full-bar action equation for a nonflat
degree-two gauge with `B=C=1`. These C3 reflection gauges happened to be
zero.

A further permanent test on untwisted C4 has native
`x=(A,B,C,D)=(0,1,0,0)` and native square `(0,0,0,6)`. Its reflection
gauge has `C=[0,1,1]` and nonzero integral D coordinates. The test verifies
complete-bar flatness of both exported native states and the full bar
reflection equation above. This exercises nonzero binary and integral
homotopy corrections in an actual nonlinear product. The fresh targeted
test passed after this case was added; the full package run above
preceded this final test-only addition.

The portable full-result example observed:

| Resolution and twists | Invariants in package degrees -1,0,1,2,3 |
| --- | --- |
| Standard C2, `s=omega=0` | `[[0],[],[2,2],[2,2],[0,8]]` |
| Standard C4, `s=[1], omega=0` | `[[],[2],[2],[2,4],[2]]` |

Both select the native model at degree 3. C2 passes the literal bar-basis
isomorphism check. C4 is a proper retract, and its result additionally
passes the independent complete-bar certificate on the same marked
generators. The original supplied resolution object is retained.

The contractible-summand C2 resolution from the initial review is also
tested through `koFull(R,0,0,3,rec(extensionModel:="transfer"))`. Its
strict retraction fails, the reference path is selected with the exact
failure reason, and the result agrees with the standard C2 example.
Thus a valid but ineligible resolution is not mistaken for an eligible
native model.

## C2 fixture results

The [saved structured results](../../data/extension-paper-results-20260926-transfer.json)
retain all 20 comparisons, selected full lifts, relation matrices, native
power witnesses and finite multiplication tables. They include `[0,8]`
in untwisted package degree 3 and `[16]` in doubly twisted package
degree 4. Expected groups were compared only after each calculation
completed.

The fixture runner was strengthened while this long process was
executing. The running function emitted its complete results using the
earlier report format, then GAP encountered a trailing read error from
the changed file position (`Variable: 'h' must have a value`). Therefore
this process is **not recorded as a clean script exit**. The completed
JSON report was independently checked: all 20 actual invariant lists
equal their fixtures, and the nontrivial degree-three and degree-four
power witnesses all have `certificateLevel="transfer-R"` and verified
equalities. Their finite audits also passed. Separate fresh calculations
explicitly asserted native selection without fallback for the other
seven distinct degree-three/four cases. The final runner parses in a
fresh process; it now additionally asserts no unresolved results and no
higher-degree fallback. That revised runner was not rerun end to end.

## Review and limits

Independent reviews covered the normalized chain homotopy, lazy formula
orchestration, signed coefficient frames, D-affine caches, native gauge
action, supplied-resolution API, reference fallback, and completion
certificate. Review corrected the need to verify marked leading
representatives in the bar certificate, retained the reference worker
across degrees, preserved legacy zero/free-layer behavior on reference
setup refusal, and strengthened the transfer fixture runner against
unresolved or fallback-only successes.

No general gauge-completeness theorem is claimed. Proper retracts need
independent complete-bar certification; failing native searches retry
the reference engine. Global predicates are exhaustive exact tests in
bounded lower degrees, rather than the plan's proposed specialized
flag optimization. Infinite groups, multiple degree-zero generators,
and failed strict retractions are outside the native path.

No complete signed-C4 degree-four timing, full degree-five classification
campaign, or new Tau/Psi/T direct-operation campaign was run. The
calibrated formulas, numerical payloads and their source hashes are
unchanged. This record does not establish the plan's projected runtime
improvements or an all-degree naturality theorem. All results retain
`certified_ko=false`.
