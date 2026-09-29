# Extensions on a supplied resolution

`koFull` and `koFull_batch` solve higher extensions in the same integral
HAP resolution used for the AHSS. They accept a supplied resolution
directly and retain the resolution object and the coordinates of its
twists:

```gap
R := ResolutionFiniteGroup(CyclicGroup(4),6);;
full := koFull(R,[1],0,3);;
full.invariants;
full.degreeResult.modelId;          # package degree 3
batch := koFull_batch(R,[1],0,3);;
batch.degreeResults[5].modelId;     # package degree 3
```

The finite-group constructor and a retained detailed E6 calculation use
the same native-resolution implementation:

```gap
full := koFull(CyclicGroup(2),0,0,3);;
ahss := koAHSS_batch(R,[1],0,3,rec(details:=true));;
full := koFull(ahss);;
```

There is no model-selection option. The supplied resolution must have
the boundary, group action, integral contracting homotopy and length
required by the page calculation; through E6 the length requirement is
`max(3,k+3)`. No resolution is reconstructed when R is supplied.

Extension states, defining-cochain solves, gauge searches and integer
matrices use R in package degrees 1–6, for finite and infinite groups and
any number of degree-zero generators. Nonlinear formulas are evaluated on
simplices of a comparison complex, lazily through sparse comparison
chains: the normalized group bar when it retracts onto R, and otherwise
the cell complex described below. No complete-bar model
is constructed for extension certification or fallback. Degrees -1 and 0
have only D and need no model. In degree six, the cutoff of the
all-cochain differential, the D-layer terms \(J_6\) and the legal
\(\gamma_6\) are evaluated on the legal lower locus only, which is the
only locus the native operations reach; a request outside it leaves the
degree unresolved (see [extensions.md](extensions.md)).

## Completion and the gauge-completeness assumption

The native calculation assumes gauge completeness: two flat bar states
represented by native states are bar-gauge equivalent exactly when their
native states are equivalent under the transferred gauge action. This is
an assumption of the computation, not a theorem established by the runtime
checks; the [transfer note](transfer.md) gives a counterexample to its
derivation from those checks alone. `koFull_batch` results and completed
native degree results record `gaugeCompletenessAssumed=true`.

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
| `batch.gaugeCompletenessAssumed` | `true`; identifies the assumption used for native completion |
| `degree.modelId` | Identifier of the native model when setup supplied one |
| `degree.gaugeCompletenessAssumed` | `true` for a completed native higher calculation |
| `degree.certificateLevel` | `"transfer-R"` for a completed native calculation |
| `degree.abelianQuotientAssumed` | `true`; the group law on gauge classes is assumed commutative and associative |

Native `canonicalComparison` records verify `target = act(gauge,canonical)`
on R. They contain `certificateLevel="transfer-R"` and have no literal
bar `boundary` field. `searchComplete` refers to the native affine family
examined.

## Exact transport and worker

The sparse comparison first checks \(fg=1\) over the integral group
ring on every basis generator through cochain degree k+2. This is
stronger than checking only trivial and sign characters. The group bar
passes only for one degree-zero generator and a contraction with
\(\partial h(e_j)=0\) on the generators; otherwise the model switches to
the cell comparison and checks it the same way.

### The cell comparison

Let \(J\) index the degree-zero generators of R and let EX be the
simplicial set of ordered tuples in \(X=G\times J\), with G acting on
the first factor. Its normalized chains form a free resolution. Put
\(g(e_i^0)=(1,i)\), \(g(e_j)=\mathrm{cone}_b\,g(\partial e_j)\) with
\(b=(1,1)\), and \(f(\sigma)=K(f(\partial\sigma))\) on anchored simplices,
extended G-equivariantly, where

\[
K=h+\partial\psi-\psi\partial,\qquad
\psi_n(x)=\sum_j\pi_j(x)\,h(e^n_j).
\]

Here \(\pi=(\pi_j)\) are the coordinates against a right inverse over Z of
the matrix of the boundaries \(\partial e^n_j\) at the identity. K is a
contraction, \(K(\partial e_j)=e_j\) and \(\partial K(e_j)=0\), so
\(fg(e_j)=K(fg(\partial e_j))=e_j\) inductively. When the boundaries of a
degree do not split off over Z, each of its generators gets a private cone
vertex \(x_j\) with its own contraction
\(k_j=h+\partial\psi_j-\psi_j\partial\), \(\psi_j=\pi_j(\cdot)h(e_j)\) and
\(\pi_j(\partial e_j)=1\). This needs \(\partial e_j\) to be primitive; a
generator with zero or non-primitive boundary leaves the degree
unresolved. Local cochain values use the sign of the twist along the
contraction path from \(e_1\) to \(f(v)\) at the first vertex v.

For a finite group on the group bar the worker keeps the multiplication
table. Otherwise vertices are integer labels assigned by GAP, which also
normalizes every simplex it is asked about; no group elements are
enumerated.

Given the checked retraction, the model constructs the normalized
homotopy on chains:

\[
q=1-gf,\qquad u=qhq,\qquad h'=u\partial u.
\]

This supplies the annihilation and square-zero side conditions, rather
than assuming they hold for the HAP comparison. The proof and cochain
signs are in [the transfer note](transfer.md). All integral arithmetic
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

The Python worker imports the bundled formulas and the shared cochain kernel.
Its GAP callbacks request sparse f or normalized-homotopy chains on demand.
The top-degree formulas are projected only along g-support during native
operations.

The formulas branch on whether a cochain vanishes (the legal, pure and
complete flags). The zero test pairs the cochain with \(g(e_j)\) for every
basis element of R in its degree, modulo two for binary cochains: a
cochain is zero when it vanishes on the resolution. For a pulled-back
cochain \(\Lambda w\) this is the same as vanishing on the whole comparison
complex, because \(\Pi\Lambda=1\); for the other intermediate cochains it
is the definition used here. No simplices are enumerated.

Native curvature returns a correctly shaped four-layer tuple whose
components are exact through the first nonzero obstruction. Later
components are zero-filled and are not obstruction values. The defining
cochain solver advances only after preceding components vanish. Full
products require flat inputs; triangular division evaluates only the
layers needed at each step and checks its final equation and flatness.

## Limits and performance

The preflight and sparse homotopy have explicit term budgets, and the
worker bounds its caches. These are implementation bounds, not
mathematical claims of nonexistence. The preflight counts g-support; it does not claim to have
built the entire face/homotopy closure in advance.

Native solving reduces matrix sizes substantially: cyclic resolutions
can have rank one in each degree even when the complete bar is too large.
The calibrated per-simplex formulas remain costly. No complete-bar size limit applies, but the native engine
has its own search and transport bounds.

The [signed C4 example](../examples/c4_signed.g) solves every degree from
-1 to 6 on a supplied resolution, and the [suspension example](../examples/suspension.g)
does the same for an infinite group on a product resolution. All results
retain `certified_ko=false` and the five-row scope.
