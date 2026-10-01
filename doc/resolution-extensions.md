# Extensions on a supplied resolution

`koFull` and `koFull_batch` solve higher extensions in the same integral
HAP resolution used for the AHSS. They accept a supplied resolution
directly and retain the resolution object and the coordinates of its
twists:

```gap
R := ResolutionFiniteGroup(CyclicGroup(4),9);;
full := koFull(R,[1],0,6);;
full.invariants;                    # [ 4 ], from the relation 2B=D
full.degreeResult.modelId;          # "transferred-normalized-bar"
batch := koFull_batch(R,[1],0,6);;
batch.degreeResults[8].modelId;     # package degree 6
batch.degreeResults[5].certificateLevel;   # package degree 3: "direct-sum"
```

Package degree three of this case has a single nonzero layer and measures
no relation, so it builds no model and records no `modelId`.

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

Defining-cochain solves and integer matrices use R in package degrees 1–6,
for finite and infinite groups and
any number of degree-zero generators. Nonlinear formulas are evaluated on
simplices of a comparison complex, lazily through sparse comparison
chains: the normalized group bar when it retracts onto R, and otherwise
the cell complex described below. No complete-bar model
is constructed for extension certification or fallback. After split and
primary-operation rows, the default path evaluates
[light residues](extensions.md#light-rows), without a flat-lift search or
reflected gauge comparison. The model, its
worker and the comparison chains are built only when a relation of the
degree needs them: degrees -1 and 0, which have only D, and degrees whose
relations all split or have zero rows by their lower groups (a line with
at most one nonzero layer, for example) build none of them. In degree six,
the cutoff of the all-cochain differential, the D-layer terms \(J_6\) and
the legal \(\gamma_6\) are evaluated on the legal lower locus only, which
is the only locus the native operations reach; a request outside it leaves
the degree unresolved (see [extensions.md](extensions.md)).

## Completion and the gauge-completeness assumption

The model measurement path, selected when light rows are disabled or a
light prime part fails for a non-resource reason, assumes gauge completeness:
two flat bar states
represented by native states are bar-gauge equivalent exactly when their
native states are equivalent under the transferred gauge action. This is
an assumption of the computation, not a theorem established by the runtime
checks; the [transfer note](transfer.md) gives a counterexample to its
derivation from those checks alone. `koFull_batch` results and completed
degree results that measured a relation in the native model record
`gaugeCompletenessAssumed=true`.

Light completion records `certificateLevel="light-R"` and
`gaugeCompletenessAssumed=false`. Its checks are the defining equations on
R, normalized transport, residue integrality and closedness, and a final
audit of the live marking versions, required precisions and target layers.
The returned `heavyMeasurements` counts entries into model measurement,
including failed attempts; light worker pairings do not count. Restarted
prime parts are recorded in `lightFallbacks`.

The engine retains the following exact checks on R:

- Flatness of every flat lift it solves, through the target layer for a
  lift solved only that far, with its defining equations.
- Each measured stacking relation against the ordered lower product,
  with a native gauge satisfying `target = act(gauge,canonical)`.
- The identity \(fg(e_j)=e_j\) over the integral group ring on every
  comparison chain \(g(e_j)\) it uses.

Stacking is also assumed to be commutative and associative on gauge classes,
so the relations determine the group and no finite multiplication table is
audited. Completed native degree results record `abelianQuotientAssumed=true`.
Free quotient coordinates split in the intended abelian abutment category.

If the model setup, the verification of a comparison chain or a search
reaches its bounds, the degree stays unresolved.
There is no retry with a complete-bar engine. Failed exact identities remain
errors; they are not treated as split extensions.

| Result field | Meaning |
| --- | --- |
| `batch.gaugeCompletenessAssumed` | `true`; identifies the assumption used for native completion |
| `degree.modelId` | Identifier of the native model, when a relation of the degree needed it |
| `degree.gaugeCompletenessAssumed` | `true` when a relation of the completed degree was measured in the native model or read from a primary operation on R |
| `degree.certificateLevel` | `"light-R"` when light rows were used without model measurement; `"transfer-R"` for model measurement; otherwise `"primary-R"` for primary-operation rows, `"prime-split"` for split odd relations, or `"direct-sum"` when every row is zero |
| `degree.abelianQuotientAssumed` | `true` with `"transfer-R"` and `"primary-R"`; the group law on gauge classes is assumed commutative and associative |
| `degree.primaryFallbacks` | The relations whose primary-operation row was unavailable and which were measured in the model, when there are any |
| `degree.heavyMeasurements`, `degree.lightFallbacks` | Number of entries into model measurement; light prime parts restarted in that model, when there are any |
| `degree.primes`, `degree.primeParts` | Per prime: the relations and their models; the generators and the primary part of the invariants |
| `degree.singleLayer` | `true` when the line has at most one nonzero layer |

Native `canonicalComparison` records verify `target = act(gauge,canonical)`
on R. They contain `certificateLevel="transfer-R"` and have no literal
bar `boundary` field. `searchComplete` refers to the native affine family
examined.

## Exact transport and worker

The comparison is chosen on R alone, before any comparison chain is built
(`KOAHSS_ExtensionGroupBarRetraction`). On the normalized bar \(g(e_j)\) is
the cone at the identity over \(g(\partial e_j)\), and f of such a cone over
a normalized cycle c is \(h(fc)\), so \(fg(e_j)=h(\partial e_j)\) once
\(fg=1\) below the degree of \(e_j\); in degree zero f sends every vertex to
the first generator. The group bar therefore satisfies \(fg=1\) over the
integral group ring through cochain degree k+2, which is stronger than
checking only trivial and sign characters, exactly when R has one
degree-zero generator and \(h(\partial e_j)=e_j\), that is
\(\partial h(e_j)=0\), for every basis element \(e_j\) of degrees 1 to k+2;
the first basis element that fails is the one at which the complete audit
of the [transfer note](transfer.md#preflight) reports the retraction
failure. The model then uses the group bar, and otherwise (or with the
testing override `cells`) the cell comparison.

Each comparison chain \(g(e_j)\) is built when it is first used, checked
against the term and support budgets and for \(fg(e_j)=e_j\) over the
integral group ring (`KOAHSS_ExtensionTransportVerifier`), and kept; a
refusal is final and leaves the degree unresolved. `model.transportAudit`
records these checks as they happen (`mode="on first use"`, the comparison
`selection`, and per degree the chains checked so far).

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
\(\pi_j(\partial e_j)=1\). This needs \(\partial e_j\) to be primitive. A
zero or non-primitive boundary has no right inverse, so its degree is
private, and the cell comparison refuses up front, before it builds any
degree, when some boundary of R through degree k+2 is zero or has
coefficients with gcd other than one ("the boundary of resolution
generator j in degree n is not primitive"); a degree with a relation to
measure in the model then stays unresolved, while the primary-operation
rows of relations with adjacent targets need no comparison
([extensions.md](extensions.md), "Primary-operation rows"). Each degree of the cell comparison (the
Smith form of its boundaries, the right inverse and the private vertex
types) is built when it is first used, after every lower degree, so every
value is the one of the construction of all degrees at once. Local cochain
values use the sign of the twist along the contraction path from \(e_1\)
to \(f(v)\) at the first vertex v.

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
Its setup carries the ranks and coboundary matrices of R but no comparison
chain. Its GAP callbacks request, on demand, sparse f or
normalized-homotopy chains of a simplex and the comparison chain
\(g(e_j)\) of one basis element
(`{"operation":"transport","kind":"g","degree":n,"basis":j}`, with
0-based j); the worker keeps every chain it receives. The top-degree
formulas are projected only along g-support during native operations.

A projection to R pairs a cochain with every chain \(g(e_j)\) of its
degree. The pairing is linear and is computed once per cochain: a sum is
paired term by term, the lift \(\Lambda v\) pairs to \(v\) (\(fg=1\)) and a
normalized-homotopy image pairs to zero (\(h'g=0\)), each identity being
checked once per degree on the chains themselves; otherwise the cochain is
evaluated on the chains.

The formulas branch on whether a cochain vanishes (the legal, pure and
complete flags). The zero test pairs the cochain with \(g(e_j)\) for every
basis element of R in its degree, modulo two for binary cochains: a
cochain is zero when it vanishes on the resolution. For a pulled-back
cochain \(\Lambda w\) this is the same as vanishing on the whole comparison
complex, because \(\Pi\Lambda=1\); for the other intermediate cochains it
is the definition used here. The test uses a pairing already computed for
the cochain; otherwise it requests the chains one at a time and stops at
the first nonzero pairing, so the later chains of the degree are not
requested. No simplices are enumerated.

Native curvature returns a correctly shaped four-layer tuple whose
components are exact through the first nonzero obstruction. Later
components are zero-filled and are not obstruction values. The defining
cochain solver advances only after preceding components vanish. Full
products require flat inputs; triangular division evaluates only the
layers needed at each step and checks its final equation and flatness.

## Limits and performance

The verification of the comparison chains and the sparse homotopy have
explicit term budgets, and the worker bounds its caches. These are
implementation bounds, not mathematical claims of nonexistence. The
support budget counts the distinct g-support simplices of the chains
verified in each degree; it does not include the face/homotopy closure.

Native solving reduces matrix sizes substantially: cyclic resolutions
can have rank one in each degree even when the complete bar is too large.
The calibrated per-simplex formulas remain costly. No complete-bar size limit applies, but the native engine
has its own search and transport bounds.

The [signed C4 example](../examples/c4_signed.g) solves every degree from
-1 to 6 on a supplied resolution, and the [suspension example](../examples/suspension.g)
solves degree six of the infinite group Z/3 × Z on a product resolution.
All results retain `certified_ko=false` and the five-row scope.
