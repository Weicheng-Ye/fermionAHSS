# Extension API and implemented scope

`koAHSS` computes the five-row associated graded through E6. `koFull`
attempts to assemble its stacking extensions in every requested degree.
The abstract abelian-group assembler accepts general finite or free
layers. Its higher-degree relation oracle solves full flat cochain tuples,
measures their actual stacking powers, and verifies their reduction in a
common marked basis. The finite cochain model and its searches are bounded;
missing relations remain unresolved. Every detailed
result retains `certified_ko=false`; see [mathematical status](mathematical-status.md).

## Detailed AHSS results

```gap
koAHSS(group, s, omega, k[, n][, options]);
```

The four-argument call still returns one raw E6 table. With integer
`n` in `[1..5]`, it still returns the raw E2-first list of `n` tables.
The physical cutoff remains `-1 <= k <= 6`, and `n` counts pages.

Add `rec(details:=true)` as the fifth argument, or as the sixth argument
after `n`, to retain representatives and the exact computation context.
`rec(details:=false)` and an empty options record retain the corresponding
raw output. Unknown option names and nonboolean `details` are errors.

```gap
ahss := koAHSS(CyclicGroup(2), 0, 0, 2, rec(details:=true));;
ahss.kind;                        # "koAHSSResult"
ahss.pages.pageNumbers;           # [6]
ahss.pages.tables[1];             # the legacy E6 invariant table
```

| Field | Meaning |
| --- | --- |
| `kind` | `"koAHSSResult"` |
| `maxDegree` | Physical cutoff `k` |
| `computedThrough` | 6 without `n`, or `n+1` with `n` |
| `pages` | Tagged page labels and raw invariant tables |
| `pageData` | Cell tables parallel to `pages.tables`, retaining groups and lift/projection maps |
| `status` | `"computed"` if every requested cell is determined; otherwise `"unresolved"` |
| `modelId`, `twists`, `calibrationId` | Group-bar model, actual normalized twist vectors and calibration convention |
| `_context` | In-memory resolution, backend, hidden cells, maps and caches |
| `scope`, `certified_ko` | `"five-row-associated-graded"` and `false` |

The private context retains `getCell(r,p,q)`, `getMap(r,p,q)`, `backend`,
`resolution`, and the same twist basis. Its `maxDegree` is the physical
cutoff; `cohomologyCapacity` is the cohomological preparation requirement.
This object is an in-memory calculation handle, not a portable JSON
serialization. Requesting details does not itself evaluate extensions.

The tagged page payload is

```gap
rec(kind:="koAHSSPages", pageNumbers:=[6], tables:=[ahss.pages.tables[1]])
```

Labels are distinct integers in `[2..6]`, stored in the intended display
order. All tables have the same display window. `koAHSSDisplay` and
`koAHSSFormat` accept raw tables/lists, tagged pages, detailed AHSS
results, and full results. Their existing table text is unchanged.

```gap
koAHSSDisplay(ahss);
koAHSSDisplay(ahss.pages, 6);
```

For tagged input, an optional page number selects its actual label and
an absent page is an error. For a raw single table, the optional number
still labels that table; a raw list still begins at E2. Display never
calculates extensions or appends an extension summary.
`koAHSSDisplay(full)` and `koAHSSFormat(full)` select E6 by default;
passing `full.pages` instead displays all retained pages. An explicit
page number can select an earlier stored page from a full result.

## The main calculation

```gap
full := koFull(group, s, omega, k);;
full := koFull(ahss);;
```

The first form creates one resolution and backend, computes E6 once and
attempts every package degree `j` in `[-1..k]`. It uses the same group,
twist and cutoff conventions as `koAHSS`. There is no page-count argument.

The reuse form requires a detailed E6 result with its retained cochain
context. It uses the same resolution and representative basis, preserving
the input's tagged page selection. Thus a detailed `n=5` input continues
to contain E2 through E6. Raw tables, an earlier-page result and a saved
page payload without the context cannot be used for extension work.

| Field | Meaning |
| --- | --- |
| `kind` | `"koFullResult"` |
| `maxDegree`, `degrees` | `k` and `[-1..k]` |
| `invariants` | Parallel lists for completed groups; explicit status records otherwise |
| `degreeResults` | Detailed extension calculation at each degree |
| `ahss`, `pages` | Shared detailed AHSS result and its tagged pages |
| `status` | `"computed"` if all degrees complete, `"partial"` if only some complete, otherwise `"unresolved"` |
| `scope`, `certified_ko` | `"five-row-stacking-model"` and `false` |

Degree `j` is at index `j+2`, including degree -1 at index 1. Every degree
record also stores `degree:=j` and its associated-graded `layers`.

```gap
full := koFull(CyclicGroup(2), 0, 0, 2);;
full.invariants;
result := full.degreeResults[2+2];;
if result.status="computed" then
    Print(result.invariants, "\n", result.basis.orders, "\n");
fi;
koAHSSDisplay(full);
```

A completed degree stores `group`, `invariants`, `basis`, `relationMatrix`,
`extensionVectors`, `smith`, `filtration` and `lowerModel`. An unresolved
degree stores its reason, pending layer and (when applicable) generator,
completed stages, measured vectors and lower model. It has no fabricated
complete group or invariant list. The public `invariants` view uses `[]`
only for a proved zero group and `[0]` for Z.

In the higher-degree engine, `layers.<name>.fullLifts[i]` retains the
immutable full representative and its defining-equation witnesses for
each marked generator. Completed higher-degree results also retain
`algebraAudit`. If that audit cannot complete, the degree remains
unresolved; `candidatePresentation` preserves the measured presentation
without presenting it as a completed group.

## Layers and measured relations

For package degree `j`, extract the E6 cells

\[
 A_j=E_6^{j-3,0},\qquad B_j=E_6^{j-2,-1},\qquad
 C_j=E_6^{j-1,-2},\qquad D_j=E_6^{j+1,-4}.
\]

Negative cohomological indices are absent zero layers. The row `q=-3`
adds no layer. Assemble D, then C, then B, then A. The leading classes
are represented by actual cochains lifted through the retained page maps.

If the next quotient has generators `q_i` of orders `m_i`, choose lifts
and measure

\[
 m_i\widetilde q_i=t_i\in H.
\]

Every `t_i` is a full integer row vector in the same identified lower
presentation. The oracle's response includes that presentation ID and a
witness. Recording only the order of a lift loses the correlations
between these rows. With lower generators `h_1,h_2` of order two:

| Measured relations | Result |
| --- | --- |
| `2*c_1=2*c_2=h_1` | Z/4 + Z/2 + Z/2 |
| `2*c_1=h_1`, `2*c_2=h_2` | Z/4 + Z/4 |

Both chosen lifts have order four in both examples. Their common versus
independent lower images determine different groups.

If `R_H` presents H and T has rows `t_i`, the enlarged presentation is

\[
 R=\begin{pmatrix}R_H&0\\-T&\operatorname{diag}(m_i)\end{pmatrix}.
\]

Free quotient generators add columns but no power-relation rows and
require no torsion query. A free quotient splits because the intended
abutment category is that of abelian groups. The higher-degree engine
still solves and retains its full representative, including a free A
generator's B, C and D choices. The implementation keeps the named presentation
columns across all four stages, so later vectors can refer to any earlier
generator without losing its embedding.

## Smith coordinates and filtration maps

The returned Smith data satisfy `smith.U * relationMatrix * smith.V =
smith.S`. Coordinates are row vectors. A vector `x` in the named
presentation maps to Smith coordinates as `x * smith.V`; retain
`smith.activeIndices` to omit unit-order factors and reduce finite
coordinates modulo the corresponding orders.

`basis.expressions` contains the active rows of `smith.inverseV`,
expressing the cyclic/free Smith generators in the named columns
`basis.generatorIds`. `basis.orders` records their actual orders, with
zero for free generators. These Smith orders can differ from the public
prime-power list `invariants`: one Smith factor of order 6 has public
invariants `[2,3]`. Do not pair those public entries directly with the
Smith basis.

These expressions are formal integer combinations of the selected layer
lifts. The production witnesses contain concrete cochain data on the
adapter's supported domain. The abstract assembler does not reconstruct
flat cochain representatives for an arbitrary supplied relation oracle.

Each record in `filtration` retains its layer name, relation matrix,
group and Smith data, together with:

- `inclusionMatrix`: previous active Smith coordinates to the new active
  Smith coordinates.
- `quotientMatrix`: new active Smith coordinates to the quotient layer's
  recorded generator coordinates.
- `quotientOrders`: orders for reducing the quotient coordinates, with
  zero denoting a free coordinate.

These maps preserve, for example, a basis vector represented by
`c_2-c_1` and the inclusion of Z as `2*Z` in a nonsplit extension of Z/2
by Z. Each measured vector also retains its `lowerGeneratorIds` and the
oracle's `lowerPresentationId`.

## Abstract assembler interface

```gap
koAHSSExtensionFromLayers(layers, oracle[, rec()]);
```

`layers` is a record with any of `A`, `B`, `C`, `D`; omitted layers are
zero. A layer is a list of independent generator orders, or a record
containing `orders` and optional caller metadata. Orders must be integers
greater than one or zero for free generators. The only currently accepted
options record is empty.

For each required torsion relation the assembler calls
`oracle(layer,i,m,lower)`. A computed response must contain
`status:="computed"`, the matching `lowerPresentationId`, an integer
`lowerCoordinates` vector of length `lower.generatorCount`, and a
`witness`. The coordinates refer to `lower.generatorIds`, not a separately
chosen basis. To defer a relation, return
`rec(status:="unresolved",reason:="...")`; `fail` in place of an oracle
also leaves required relations unresolved.

```gap
sameImage := function(layer, i, m, lower)
    return rec(status:="computed",
        lowerPresentationId:=lower.presentationId,
        lowerCoordinates:=[1,0],
        witness:=rec(kind:="supplied-relation"));
end;;
extension := koAHSSExtensionFromLayers(
    rec(D:=[2,2],C:=[2,2]),sameImage);;
extension.invariants;               # [2,2,4]
```

The assembler inserts the D-layer orders directly when there is no lower
presentation. It queries all torsion generators of a later layer against
one common lower presentation before adjoining that layer's generators.
It then computes the joint Smith form. Wrong presentation IDs, malformed
vectors and inconsistent exact data are errors, not unresolved results.

## Complete flat representatives and boundary comparisons

In degrees 3–5, the production engine uses the selected four-cochain
`d` and `xtimes` on a complete normalized bar model of the finite group.
The leading E6 cochains are transported into that fixed basis once.
All independent generators of D, C, B and A receive a full flat lift,
including free generators and generators over a zero lower group.

For an A generator, it solves

\[
 \delta_s A=0,\qquad \delta B+P(A)=0,\qquad
 \delta C+\tau'(A;B)=0,\qquad \delta_sD+J(A,B,C)=0.
\]

The B and C equations are solved over F2 and the D equation over Z.
The solver searches the affine B and C solution families lazily: if a
chosen B does not admit C, or a chosen C does not admit D, it changes
that defining choice and retries the dependent equations. It evaluates
the actual full differential to verify the resulting tuple is flat.
B- and C-layer generators undergo the same lower-equation procedure.
The solver stores the chosen primitives, their adjustments, the exact
equation data and the final zero-curvature witness. No universal helper
is replaced by a locally solved formula.

Each `fullLifts[i].state` is an immutable tuple `(A,B,C,D)` in this complete
bar basis. The exact same tuple is reused in all powers, lower reductions
and basis comparisons. A generator's lower components are not reset to
zero or reconstructed from their cohomology classes in a later query.

For a torsion generator of order `m`, the oracle stacks this complete
tuple to obtain its measured power. It reduces the result using the
previously stored lower lifts and retains all integral carries and
boundary data. The resulting coordinates refer to the named D/C/B/A
columns used by the common relation matrix.

An ordered reduction alone does not justify rearranging cochain products.
The engine therefore constructs the canonical lower product in its
recorded generator order and solves the exact comparison

\[
 \text{measured power}=d(g)\mathbin{\times}\text{canonical lower product}.
\]

The gauge `g` includes all four components in the preceding degree. Its
support is searched in the order

1. `(0,0,0,D)`, using an ordinary integral D primitive;
2. `(0,0,C,D)`, allowing closed C gauges and their full D contributions;
3. `(0,B,C,D)`, solving consistent C/D choices for a B gauge;
4. `(A,B,C,D)`, solving all remaining defining choices for an A gauge.

Earlier components are fixed to literal zero at each stage. The solver
continues to larger support when a smaller stage fails or reaches its
share of the search budget. Native cohomology representatives are tried
first as candidate directions; the complete cochain kernel remains a
fallback, including ordinary coboundary directions and their nonlinear
carries. A bounded integer-kernel search is not reported as exhaustive.
The leading primitive and subsequent B/C/D choices are solved consistently;
the comparison retains `g`, `d(g)` and the product, verifies `d(d(g))=0`,
and checks the displayed equality component by component. It does not
assume strict associativity or commute factors at the cochain level.
Failure to find this exact witness within the search bounds remains
unresolved.

If ordinary layer reduction fails at D, the reducer projects the remaining
D cocycle into the retained E6 cell to propose coordinates in the marked
D basis. This accounts for incoming differentials that ordinary
coboundaries alone do not remove. It then compares the **original full
stacked tuple** with the canonical product of all recorded lower lifts,
including their earlier carries, using the staged gauge search above.
The E6 projection alone never certifies a stacking relation.

For other failed ordinary reductions, the fallback tries at most 32
marked finite lower normal forms, stage by stage, using the same stored
representatives. It does not guess coordinates for a free lower factor.
Accepted relations retain `canonicalComparison.winningStage`,
`attemptedStages` and the full gauge. Unresolved relation queries retain
`pendingRelation`, including attempted comparisons and residual cochains,
instead of claiming that a particular missing gauge necessarily exists.

Before returning a higher-degree result as computed, the finite quotient
audit constructs every marked finite normal form, checks every ordered
pair product against its Smith-coordinate target by literal equality or
an exact boundary comparison, and verifies the resulting table has an
identity, inverses, commutativity and associativity. The audit handles at
most 32 finite normal forms. Free quotient coordinates are excluded from
this finite enumeration and split in the intended abelian abutment
category; this is not a finite verification of every infinite cochain
product.

The degree-six path additionally uses the fixed pointed finite section
and exact basis conversions in
[production_g6_section.py](../python/stacking_model/production_g6_section.py). Its section
search and complete-state comparisons are separate from the calibrated
E6 page computation. The existence of a degree-six formula or adapter
does not by itself certify a completed extension calculation; the same
flat-lift, relation, boundary-comparison and quotient-audit requirements
apply. Verification records state which calculations were actually run.

The runtime bundles the selected stacking sources under
[python/stacking_model](../python/stacking_model/) with
[source provenance and checksums](../python/stacking_model/provenance.json).
Loading the package does not require the separate research workspace.
No expected classification table is consulted at runtime.

## Low-degree adapter and resource limits

The current group-bar adapter implements these relation measurements:

| Package degree | Supported torsion query | Conditions |
| --- | --- | --- |
| 1 | C to D | Order-two C generator; arbitrary valid `s,omega`; integral upper lift and retained page projection must succeed |
| 2 | C to D, B to C/D | Order-two quotient generator; `omega=0`; closed binary representatives and the implemented lower reduction |

Each queried order-two generator uses one doubling with the selected
`xtimes` formula. Lower-coordinate reduction is separate work and retains
its carries and operation counts. The degree-one witness includes the
upper integral primitive and bar-comparison homotopy correction. The
degree-two phase is the exact 256-value restriction of the selected
low-degree commutative stacking formula to closed B/C inputs at
`omega=0`, bundled with provenance in
[stacking-low-phase.json](../data/stacking-low-phase.json). No expected
classification table is consulted at runtime.

Degrees -1 and 0 require only D. A general nonzero-`omega` degree-two
gauge reducer is not supplied by this low-degree adapter; a required
query outside its domain remains unresolved. The abstract assembler's
generality does not remove a production cochain requirement.

The complete-bar implementation has explicit resource bounds:

| Work | Current bound |
| --- | --- |
| Complete bar cochains in a required degree | 8192 coordinates |
| A required coboundary matrix | 2,000,000 entries |
| Flat-lift affine search | 4096 distinct differential evaluations by default |
| Gauge-comparison affine search | 4096 equation evaluations and 64 leading choices shared across requested stages; integral kernel coefficients initially bounded by absolute value 1 |
| Fallback lower-coordinate search | 32 marked finite normal forms; 4096 equation evaluations shared across stages and candidates |
| Degree-six finite-section search | 4096 candidates |
| Finite quotient audit | 32 marked normal forms |

These are implementation limits; they add no arguments to `koFull`.
Repeated integer solves now reuse exact Smith preparations, bounded by
eight entries and two million retained matrix cells including transforms.
The nonlinear higher model still uses the complete bar. The
[resolution-transfer review](transfer.md) documents the implemented
preflight and the mathematical prerequisites remaining before a new
production model can preserve its classification results.
These bounds do not make every operation inexpensive: the bar dimensions
grow with the group order and degree, and one nonlinear evaluation can
be costly. Reaching a resource bound or failing to find a gauge within
the bounded search is not evidence that an extension splits. Exact
identity failures and inconsistent coordinate data are errors; unavailable
witnesses leave the corresponding result unresolved.

The [paper comparison record](extension-paper-comparisons.md) uses spatial
dimension `d`, related to package degree by `j=d+1`. Paper dimension `d`
therefore appears in `full.invariants[d+3]`. Compare complete invariant
lists with the same twists and phase convention; neither E6 layer-order
products nor unresolved outputs establish a full-group match.
