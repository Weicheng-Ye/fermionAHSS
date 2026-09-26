# Native extension default: verification on 2026-09-26

Base revision: `c2961a2`; the working tree was clean at the start.
This revision implements the requested assumption of gauge completeness.
It supersedes the runtime certification policy in the earlier
[resolution-extension record](resolution-extensions-20260926.md), which
is retained as historical evidence.

## Runtime change

`koFull(groupOrResolution,s,omega,k)` and `koFull(detailedE6)` now use
native resolution coordinates in degrees 3–5. The options argument,
complete-bar completion certificate, certification module, exceptional
C2 certification bypass, and automatic reference fallback were removed.
Completed native results record `gaugeCompletenessAssumed=true`; this
is an assumption rather than a newly proved theorem.

Exact native flatness, gauge-action equations, relation measurements,
finite quotient audits, and sparse comparison identities remain checked.
The formulas and their calibration were not modified. Sparse bar
transport and exact lower-degree branch predicates remain part of the
native formulas; neither constructs the complete-bar extension model.
Explicit complete-bar comparisons remain available only as developer
reference helpers and regression tests, outside the `koFull` call path.

Degree 6 has no native implementation and returns unresolved without
constructing a complete-bar model. Native setup and search failures also
remain unresolved without a bar retry. The supported resolution and
resource requirements are listed in
[resolution extensions](../resolution-extensions.md).

## Executed checks

All GAP checks below used fresh processes with `--quitonbreak`.

| Check | Observed result |
| --- | --- |
| `tst/extension_transfer_hooks.tst` | Passed: native C8 completion with the complete-bar factory replaced by an error; exact supplied-resolution identity; no certification hook; worker closure; options rejected; setup failures and degree 6 do not fall back |
| `tst/api.tst` and `tst/extensions.tst` | Passed |
| `LoadPackage("fermionAHSS"); Assert(0,TestPackage("fermionAHSS"));` | Passed, exit 0; 68.835 seconds, including all native-default regressions |
| `examples/resolution_extensions.g` | Passed: C2, signed C4 and C8 through package degree 3 using the plain API |
| `examples/c2.g` and `examples/twisted_c2.g` | Passed: unchanged E2–E6 page examples |
| C7 through package degree 4 | Passed: plain `koFull(CyclicGroup(7),0,0,4)` returns computed; degree 3 is `[0,7]`, degree 4 is `[]`, both using the native model |
| Fixture example parsing | Both `extension_papers.g` and its historical alias parse in a fresh loaded GAP process |
| Formula integrity | All 58 bundled SHA-256 hashes match `python/stacking_model/provenance.json` |
| Independent code review | No actionable runtime issues found; worker cleanup, removed API forms and absence of complete-bar calls checked |

The supplied-resolution example produced:

| Resolution and twists | Invariants in package degrees -1,0,1,2,3 |
| --- | --- |
| C2, `s=omega=0` | `[[0],[],[2,2],[2,2],[0,8]]` |
| C4, `s=[1], omega=0` | `[[],[2],[2],[2,4],[2]]` |
| C8, `s=omega=0` | `[[0],[],[2,8],[2,2],[0,2,16]]` |

The C8 regression calls the real native solver after replacing
`KOAHSS_ExtensionBarModel` with a function that raises an error. The full
calculation and native finite quotient audit succeed. This directly
checks that the former 8192-coordinate complete-bar gate is not consulted:
the corresponding degree-five bar module would have `7^5=16807`
coordinates, while the supplied cyclic resolution has rank one.

The prior 20-fixture calculation is historical evidence; that expensive
calculation was not rerun in this revision. The native nonlinear worker
is unchanged. No new universal completeness or all-degree naturality
claim is made, and `certified_ko` remains false.
