# Extensions on a supplied resolution

`koFull` solves higher extensions in the same integral HAP resolution used
for the AHSS. It accepts a supplied resolution directly and retains the
resolution object and the coordinates of its twists:

```gap
R := ResolutionFiniteGroup(CyclicGroup(4),6);;
full := koFull(R,[1],0,3);;
full.invariants;
full.degreeResults[5].modelId;  # package degree 3
```

The finite-group constructor and a retained detailed E6 calculation use
the same native-resolution implementation:

```gap
full := koFull(CyclicGroup(2),0,0,3);;
ahss := koAHSS(R,[1],0,3,rec(details:=true));;
full := koFull(ahss);;
```

There is no `extensionModel` option. The supplied resolution must have
the boundary, group action, integral contracting homotopy and length
required by the page calculation; through E6 the length requirement is
`max(3,k+3)`. No resolution is reconstructed when R is supplied.

Extension states, defining-cochain solves, gauge searches and integer
matrices use R in package degrees 3–5 when its fixed comparison passes
the retraction check. Nonlinear formulas still use bar simplices,
evaluated lazily through sparse comparison chains. No complete-bar model
is constructed for extension certification or fallback. Degrees -1 through
2 retain their existing adapters. Degree-six extensions are explicitly
unresolved because the native degree-six formula is not implemented.

## Completion and the gauge-completeness assumption

The native calculation assumes gauge completeness: two flat bar states
represented by native states are bar-gauge equivalent exactly when their
native states are equivalent under the transferred gauge action. This is
an assumption of the computation, not a theorem established by the runtime
checks. The [initial mathematical review](transfer.md) records the gap in
the original justification. Full results and completed native degree results record
`gaugeCompletenessAssumed=true`.

The engine retains the following exact checks on R:

- Flatness of every chosen full lift, with its defining equations.
- Each measured stacking relation against the ordered lower product,
  with a native gauge satisfying `target = act(gauge,canonical)`.

Stacking is also assumed to be commutative and associative on gauge classes,
so the relations determine the group and no finite multiplication table is
audited. Completed native degree results record `abelianQuotientAssumed=true`.
Free quotient coordinates split in the intended abelian abutment category.

If setup or search reaches its bounds, the degree stays unresolved.
There is no retry with a complete-bar engine. Failed exact identities remain
errors; they are not treated as split extensions.

| Result field | Meaning |
| --- | --- |
| `full.gaugeCompletenessAssumed` | `true`; identifies the assumption used for native completion |
| `degree.modelId` | Identifier of the native model when setup supplied one |
| `degree.gaugeCompletenessAssumed` | `true` for a completed native higher calculation |
| `degree.certificateLevel` | `"transfer-R"` for a completed native calculation |
| `degree.abelianQuotientAssumed` | `true`; the group law on gauge classes is assumed commutative and associative |

Native `canonicalComparison` records verify `target = act(gauge,canonical)`
on R. They contain `certificateLevel="transfer-R"` and have no literal
bar `boundary` field. `searchComplete` refers to the native affine family
examined. There are no public model-selection, failed-transfer or
independent bar-certificate records.

## Exact transport and worker

The sparse comparison first checks \(fg=1\) over the integral group
ring on every basis generator through cochain degree k+2. This is
stronger than checking only trivial and sign characters. The current
native path also requires one degree-zero generator and a finite group.
A supplied HAP resolution may therefore be accepted by the AHSS API while
its higher extension calculation remains unresolved. Infinite groups,
multiple degree-zero generators and comparisons failing the strict
retraction identity are outside the current native extension domain.

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
gauge is never substituted for that full bar boundary. These evaluations
use sparse requested simplices, not a complete-bar coordinate array.

The Python worker imports the existing checksum-verified formulas. Its
GAP callbacks request sparse f or normalized-homotopy chains on demand.
The top-degree formulas are projected only along g-support during native
operations.

To avoid ambiguous branch choices, the implementation evaluates global
zero predicates exactly on all normalized simplices in the required lower
degree. It does not infer them from a sample or from g-support. This
conservative choice includes P1/P2 automatically but is more expensive
than the plan's proposed specialized flag tables. These lower-degree
predicate checks are distinct from constructing a complete bar model with
its cochain arrays and matrices. A predicate exceeding its budget returns
unresolved before sampling.

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

Native solving reduces matrix sizes substantially: cyclic resolutions
can have rank one in each degree even when the complete bar is too large.
The calibrated per-simplex formulas and lower-degree predicate checks
remain costly. Removing complete-bar certification removes that size gate;
it does not remove the native engine's separate search and transport bounds.

The portable [resolution example](../examples/resolution_extensions.g)
exercises C2, signed C4 and C8 extensions on supplied resolutions. Current
checks, including the removal of the finite audit and the exact worker
evaluation policy, are recorded in
[the acceleration verification record](verification/extension-acceleration-20260926.md).
The [native-default record](verification/native-extension-default-20260926.md)
is retained as historical evidence. Historical transfer and reference
comparisons remain in
[the initial implementation verification record](verification/resolution-extensions-20260926.md).
All results retain `certified_ko=false` and the existing five-row scope.
