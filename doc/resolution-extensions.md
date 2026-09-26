# Extensions on a supplied resolution

`koAHSS` and `koFull` now accept an integral HAP resolution directly.
The resolution object and the coordinates of its twists are retained:

```gap
R := ResolutionFiniteGroup(CyclicGroup(4),6);;
full := koFull(R,[1],0,3,rec(extensionModel:="transfer"));;
full.invariants;
full.degreeResults[5].modelSelection;  # package degree 3
```

The same option works with the finite-group constructor and a retained
detailed E6 calculation:

```gap
full := koFull(CyclicGroup(2),0,0,3,rec(extensionModel:="transfer"));;
ahss := koAHSS(R,[1],0,3,rec(details:=true));;
full := koFull(ahss,rec(extensionModel:="transfer"));;
```

`extensionModel` accepts `"bar"` and `"transfer"`. The default remains
`"bar"`, preserving existing calls. The supplied resolution must have
the boundary, group action, integral contracting homotopy and length
required by the page calculation; through E6 the length requirement is
`max(3,k+3)`. No resolution is reconstructed when R is supplied.

This release moves the extension states, defining-cochain solves, gauge
searches and integer matrices to R in package degrees 3–5 when its fixed
comparison passes the retraction check. Nonlinear formulas still use
bar simplices, evaluated lazily through the comparison. There is no
complete-bar matrix in the native search. Degrees -1 through 2 retain
their existing adapters, and degree 6 retains the complete-bar section.

## Completion and fallback

The distinction between a homotopy equivalence and a strict retraction,
and between flatness transfer and gauge completeness, matters here. The
[initial mathematical review](transfer.md) explains both issues. This
implementation does not assume the missing general gauge-completeness
theorem.

A completed native presentation is accepted in either of these cases:

1. The comparison is checked to be the literal C2 normalized-bar basis
   isomorphism in every required degree. The native state coordinates are
   then the reference coordinates themselves.
2. The presentation passes an independent complete-bar certification.
   Every stored native lift is exported using the full corrected
   embedding, including its D correction. Its earlier layers must be
   zero and its leading layer must equal the transported marked E6
   representative. The reference engine checks flatness, recomputes every
   torsion power and ordered lower product, solves the literal complete
   gauge comparison, and performs its existing finite quotient audit.

The second route uses the same acceptance criteria as the original
complete-bar engine on the native-selected generators. It is not a claim
that every bar gauge comes from a native gauge. It retains genuine bar
witnesses even when the native gauge search is incomplete.

If native setup, search, or completion certification remains unresolved,
`koFull` retries that degree with the reference bar engine. Its result
retains the unsuccessful native attempt. If both paths reach their
limits, the degree remains unresolved. Failed exact identities remain
errors; they are not treated as resource limits or split extensions.

| Result field | Meaning |
| --- | --- |
| `full.extensionModel` | Requested model |
| `degree.modelSelection` | Requested and selected model, fallback flag, and reason when applicable |
| `degree.transferAttempt` | Retained unresolved native attempt when the reference path was selected |
| `degree.barCertification` | Independent complete-bar lifts, relation comparisons and finite audit, when required |
| `degree.certificateLevel` | `"transfer-R"` for the checked C2 isomorphism, or `"complete-bar-certified-transfer"` after reference certification |

Native `canonicalComparison` records verify `target = act(gauge,canonical)`
on R. They contain `certificateLevel="transfer-R"` and have no literal
bar `boundary` field. The independent bar certificate retains the usual
`target = d(gauge) xtimes canonical` witnesses. Code consuming both kinds
of record must check for `boundary` before reading it. `searchComplete`
in a native search refers to the native affine family examined, not all
bar gauges.

## Exact transport and worker

The sparse comparison first checks \(fg=1\) over the integral group
ring on every basis generator through cochain degree k+2. This is
stronger than checking only trivial and sign characters. The current
native path also requires one degree-zero generator and a finite group.
An arbitrary supplied HAP resolution can therefore enter the API even
when it is ineligible for native transfer: its reference computation is
used with an explicit fallback reason. Infinite and multiple-orbit
resolutions are not silently treated as strict retractions.

Given the checked retraction, the model constructs the normalized
homotopy on chains:

\[
q=1-gf,\qquad u=qhq,\qquad h'=u\partial u.
\]

This supplies the annihilation and square-zero side conditions, rather
than assuming they hold for the original HAP comparison. The proof and
cochain signs are in [the review](transfer.md). All integral arithmetic
is exact; B and C are reduced modulo two before subsequent lifts.

The embedding of a native state is triangular:

\[
\Phi(w)_A=\Lambda w_A,\qquad
\Phi(w)_\ell=\Lambda w_\ell-HN_\ell(\Phi(w)_{<\ell}).
\]

Here \(\Lambda=f^*\), \(H=(h')^*\), and \(N_\ell\) is the fixed
nonlinear differential term in layer \(\ell\). Products reflect the
actual bar product of the embedded states back to R. Gauge actions
reflect `d(Phi(e)) xtimes Phi(canonical)`; the curvature of a native
gauge is never substituted for that full bar boundary.

The Python worker imports the existing checksum-verified formulas. Its
GAP callbacks request sparse f or normalized-homotopy chains on demand.
The top-degree formulas are projected only along g-support during native
operations. Whole-bar materialization occurs when requesting the separate
reference certificate.

To avoid ambiguous branch choices, the first implementation evaluates
global zero predicates exactly on all normalized simplices in the
required lower degree. It does not infer them from a sample or from
g-support. This conservative choice includes P1/P2 automatically but is
more expensive than the plan's proposed specialized flag tables. A test
exceeding its budget returns unresolved before sampling.

Native curvature returns a correctly shaped four-layer tuple whose
components are exact through the first nonzero obstruction. Later
components are zero-filled and are not obstruction values. The defining
cochain solver advances only after preceding components vanish. Full
products require flat inputs; triangular division evaluates only the
layers needed at each step and checks its final equation and flatness.

## Limits and performance

The preflight and sparse homotopy have explicit term budgets. The worker
limits each complete lower-degree flag test to 8192 simplices and bounds
its caches. These are implementation bounds, not mathematical claims of
nonexistence. The initial preflight counts g-support; it does not claim
to have built the entire face/homotopy closure in advance.

Native solving reduces matrix sizes substantially: for a cyclic C4
resolution the cochain rank is one in each degree. The calibrated
per-simplex formulas remain costly, however, and the conservative bar
certificate still incurs complete-bar work and its existing limits.
This release therefore does not deliver a universally bar-free solver
or the plan's projected end-to-end speedup for larger groups.

The portable [resolution example](../examples/resolution_extensions.g)
compares complete native/reference results for C2 and signed C4 through
package degree 3. The [paper comparison wrapper](../examples/extension_papers_transfer.g)
runs the existing C2 fixtures with transfer selected. Executed outcomes
and uncompleted checks are recorded separately in
[the implementation verification record](verification/resolution-extensions-20260926.md).
All results retain `certified_ko=false` and the existing five-row scope.
