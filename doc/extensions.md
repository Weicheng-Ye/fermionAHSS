# Extension API and implemented scope

`koAHSS` and `koAHSS_batch` compute the five-row associated graded
through E6. `koFull` attempts to assemble its stacking extension in one
degree, and `koFull_batch` in every degree up to a cutoff.
The abstract abelian-group assembler accepts general finite or free
layers. In degrees 1–6 the relations are recorded prime by prime; for each
relation that needs a measurement, the higher-degree relation oracle solves
flat cochain tuples, measures their actual stacking powers, and verifies
their reduction in the common marked basis. The finite cochain model and
its searches are bounded; missing relations remain unresolved. Every detailed
result retains `certified_ko=false`; see [mathematical status](mathematical-status.md).

## Detailed AHSS results

```gap
koAHSS_batch(groupOrR, s, omega, k[, n][, options]);
koAHSS(groupOrR, s, omega, k[, n][, options]);
```

The four-argument call of `koAHSS_batch` returns one raw E6 table. With
integer `n` in `[1..5]`, it returns the raw E2-first list of `n` tables.
The physical cutoff is `-1 <= k <= 6`, and `n` counts pages. `koAHSS`
takes the same arguments and returns the entries of those pages on the
line `p+q=k-3` as a record of kind `"koAHSSLine"`; see the
[README](../README.md#ahss-pages).

Add `rec(details:=true)` as the fifth argument, or as the sixth argument
after `n`, to retain representatives and the exact computation context.
`rec(details:=false)` and an empty options record return the corresponding
raw output. Unknown option names and nonboolean `details` are errors.

```gap
ahss := koAHSS_batch(CyclicGroup(2), 0, 0, 2, rec(details:=true));;
ahss.kind;                        # "koAHSSResult"
ahss.pages.pageNumbers;           # [6]
ahss.pages.tables[1];             # the raw E6 invariant table
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
results, `koAHSS` lines, and the results of `koFull` and `koFull_batch`.

```gap
koAHSSDisplay(ahss);
koAHSSDisplay(ahss.pages, 6);
```

For tagged input, an optional page number selects its actual label and
an absent page is an error. For a raw single table, the optional number
labels that table; a raw list begins at E2. Display never calculates
extensions. For a `koFull_batch` result, `koAHSSDisplay(full)` and
`koAHSSFormat(full)` select E6 by default and append no extension
summary; passing `full.pages` instead displays all retained pages, and an
explicit page number can select an earlier stored page. A `koFull` result
is displayed as its E6 line followed by the group of its degree.

## The main calculation

```gap
full := koFull(group, s, omega, k);;
full := koFull(R, s, omega, k);;
full := koFull(ahss);;
batch := koFull_batch(group, s, omega, k);;
batch := koFull_batch(R, s, omega, k);;
batch := koFull_batch(ahss);;
```

The forms with a group create one resolution and backend and compute E6
once. `koFull` then attempts the extension in package degree `k`, and
`koFull_batch` in every package degree `j` in `[-1..k]`. They use the same
group, twist and cutoff conventions as `koAHSS`. There is no page-count
argument. An integral HAP resolution may replace the group argument; it is
retained without reconstruction. There is no model-selection option. The
[native-resolution model](resolution-extensions.md) is built only when a
relation of the degree needs a measurement; a two-primary relation whose
target layer lies right below its generator's layer is read from a primary
operation on R instead (see "Primary-operation rows"), and the other
relations are light rows (see "Light rows") unless
`FERMIONAHSS_LIGHT_RELATIONS=0`. The model runs
higher searches on R,
checks the exact sparse comparison identity of each comparison chain when
it first uses it, and assumes gauge completeness. Runtime completion uses
native flatness and gauge checks.
It never constructs a complete-bar model for certification or retries an
unresolved degree with one. Native gauge records use the action equation
described below.

The reuse forms require a detailed E6 result with its retained cochain
context, and `koFull(ahss)` solves degree `ahss.maxDegree`. They use the
same resolution and representative basis; `koFull_batch` also preserves
the input's tagged page selection, so a detailed `n=5` input continues to
contain E2 through E6. Raw tables, an earlier-page result and a saved page
payload without the context cannot be used for extension work.

A `koFull` result has these fields:

| Field | Meaning |
| --- | --- |
| `kind` | `"koFullDegreeResult"` |
| `k` | The package degree |
| `invariants` | The abelian invariants of the group, when it is computed |
| `reason` | The reason, when the degree is unresolved |
| `status` | `"computed"` or `"unresolved"` |
| `line` | The E6 entries on the line `p+q=k-3`, as a `koAHSS` line record |
| `degreeResult` | The detailed extension calculation of degree `k` |
| `scope`, `certified_ko` | `"five-row-stacking-model"` and `false` |

A `koFull_batch` result has these fields:

| Field | Meaning |
| --- | --- |
| `kind` | `"koFullResult"` |
| `maxDegree`, `degrees` | `k` and `[-1..k]` |
| `invariants` | Parallel lists for completed groups; explicit status records otherwise |
| `degreeResults` | Detailed extension calculation at each degree |
| `gaugeCompletenessAssumed` | `true`; native equivalence is assumed to capture bar gauge equivalence |
| `ahss`, `pages` | Shared detailed AHSS result and its tagged pages |
| `status` | `"computed"` if all degrees complete, `"partial"` if only some complete, otherwise `"unresolved"` |
| `scope`, `certified_ko` | `"five-row-stacking-model"` and `false` |

Degree `j` is at index `j+2`, including degree -1 at index 1. Every degree
record, including `full.degreeResult`, also stores `degree:=j` and its
associated-graded `layers`. When a relation of the degree needed the
native model, its `modelId` is retained, and the completed calculation has
`certificateLevel="transfer-R"`, `gaugeCompletenessAssumed=true` and
`abelianQuotientAssumed=true`. A completed degree with light rows and no
measurement in the model has `certificateLevel="light-R"`,
`gaugeCompletenessAssumed=false`, `abelianQuotientAssumed=true` and the
list `lightShortcuts` of the page forms and absorptions its rows used. A
completed degree of 1–6 that built no
model has `certificateLevel="primary-R"`, with the same two assumptions,
when some relation row is the class of a primary operation on R. Otherwise
it has `"prime-split"` when some relation split at an odd prime and
`"direct-sum"` if not: every relation row is then zero and the group is the
direct sum of its layers.

```gap
full := koFull(CyclicGroup(2), 0, 0, 2);;
result := full.degreeResult;;
if result.status="computed" then
    Print(result.invariants, "\n", result.basis.orders, "\n");
fi;
koAHSSDisplay(full);
```

A completed degree stores `group`, `invariants`, `basis`, `relationMatrix`,
`extensionVectors`, `smith`, `filtration` and `lowerModel`; in degrees 1–6
also `certificateLevel`, the per-prime summary `primes` and the prime parts
`primeParts` (see "Localization at the primes"), and `singleLayer=true`
when its line has at most one nonzero layer. A degree whose
primary-operation rows fell back to the model lists those relations in
`primaryFallbacks` (see "Primary-operation rows"), and one whose light
prime part was computed again in the model lists it in `lightFallbacks`
(see "Light rows"). An unresolved
degree stores its reason, pending layer and (when applicable) generator,
completed stages, measured vectors and lower model. It has no fabricated
complete group or invariant list. The public `invariants` view uses `[]`
only for a proved zero group and `[0]` for Z.

In the higher-degree engine, `layers.<name>.fullLifts[i]` retains the
immutable full representative and its defining-equation witnesses of each
marked generator whose complete or three-local lift a measurement needed
(see "Complete flat representatives and gauge comparisons"). Completed
results with `certificateLevel="transfer-R"` or `"primary-R"` record
`abelianQuotientAssumed=true`: the measured relations are combined under
the assumption that stacking is commutative and associative on gauge
classes, and no finite multiplication table is audited.

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

The rows `t_i` need not be measured completely. The isomorphism type of
the extension of a quotient generator `q` of order `m` by the lower group
`H` is the class of `m\widetilde q` in `Ext(Z/m,H)=H/mH`. Call a lower
generator `e` *free for `q`* when `e\notin mH`, that is, when its class in
`H/mH` is nonzero (the integer system `m\,y+r\,R_H=e` has no solution).
The **target layer** of the relation is the lowest layer of the lower
presentation containing a free generator; every generator below it lies in
`mH`, so the components of `t_i` in those layers cannot change the class,
and the oracle measures the power only through the target layer
(`KOAHSS_ExtensionTargetLayer`), which it decides from the lower
presentation alone, before it needs a model or a lift. The stacked power
and the reductions are then layer-limited (`upto` in the model's `d`,
`xtimes`, `act` and `divideLeft`), so the D-layer stacking correction is
never evaluated when the target layer is B or C; the unmeasured entries of
`t_i` are recorded as zero, and the witness carries `measuredLayers`,
`truncatedBelow` and the certificate `sufficiency`, the explicit integer
combinations exhibiting every lower generator below the target layer as an
element of `mH` modulo the relations. A lower generator of the target layer
enters the compared product with its marked cocycle only (the state with
that one layer), which agrees there with its flat lift.

When the target layer `T` lies right below the generator's layer (B for an
A generator, C for a B generator), the generator's flat lift is solved only
through `T` (`KOAHSS_ExtensionFlatLift` with `upto`, the lift's witness
recording `solvedThrough`): the order is even (for odd `m`, `H/mH=D/mD`
because the binary layers are two-groups, so the target layer is D or
none) and `T` is binary, so the power and its reduction through `T` do not
depend on the choice in `T`, and no layer below `T`, in particular no
D-layer curvature, is evaluated for the relation. The target layer C of an
A generator uses the complete lift: the choices in B change the power
through C (the product term `\beta`) and must be those of an element of the
group, which needs the extension through D. When no lower generator is
free, the class is zero and the relation `m\widetilde q=0` is recorded
without any measurement, model, lift or comparison chain; when a free
generator lies in the D layer, the measurement is complete. If the
layer-limited reduction fails to express a component through ordinary
coboundaries, the complete measurement with the complete lift takes over.
The environment variable `FERMIONAHSS_LAYERED_RELATIONS=0` disables the
shortcut.

For the Pin⁻ line of `C2` with `s` at degree six, the lower group of the A
generator is `<b,c,d\mid 2b=c,2c=d,2d=0>\cong Z/8`: only `b` is free, the
target layer is B, and in the model the lift of the A generator is solved
through B only, and the B coordinate of its relation is the B component
`Q_D(rho u)` of the boundary state of the A gauge (`delta_s u=2A`; the B
layer `alpha(A,A)` of the square is a coboundary), so `Z/16` is determined
without the universal pair source of `gamma_6`; the class is `D(rho u)`,
which the primary-operation row reads on R without the model (see
"Primary-operation rows"). With a lower group
`Z/2\oplus Z/4` (`2b=0`, `2c=d`) both `b` and `c` are free and the target
layer is C: `2a=b` gives `Z/4\oplus Z/4` while `2a=b+c` gives
`Z/8\oplus Z/2`, which only the C component decides.

A truncated row presents its own stage, but it is the relation of another
lift of its generator. If the components of `t_i` below the target layer
form `t_i''\in mH` and `mh=t_i''`, the recorded row `t_i-t_i''` is the
relation of `\widetilde q_i-h`; a zero row over a lower group with `mH=H`
is likewise the relation of `\widetilde q_i-h` with `mh=t_i`. In both cases
`h` lies in the lower group of `q_i`, the layers below it. Every later
measurement multiplies the flat lift `\widetilde q_i` itself, so a later
relation `m'\widetilde u=t'` with a coefficient `c` on `q_i` is a relation
of the recorded columns only up to `ch`. When the target layer of `u` lies
at or above the layer of `q_i`, `ch` lies in the layers below that target,
all of which lie in `m'H'`, and the row of `u` stands. Otherwise, when
`c\neq0`, the relation of `q_i` is measured through D, with the complete
lift and without the shortcut (`witness.measuredThroughD`); if that
measurement is unresolved, so is the relation of `u`. The row of `u` itself
is kept: a measurement multiplies the flat lifts and reads the lower
presentation only for its target layer, which depends on the filtration of
the lower group and not on its presentation. The target-layer certificates
of the rows above a row measured again are restated in the final lower
presentations. With every layer `Z/2` and the relations `2c=d`, `2b=d`,
`2a=b` of the flat lifts, for example, the B relation has target layer C
and is recorded as `2b=0`, the relation of `b-c`; the A relation has target
layer C and the coefficient one on `b`, and without the measurement of `b`
through D the rows would present `Z/4\oplus Z/4` instead of
`Z/2\oplus Z/8`. With prime localization on, a relation without a free
lower generator takes the zero row of the zero-local shortcut, which is
exact after localization at its prime (see "Localization at the primes"),
so only a B relation with target layer C and a C relation read from its
primary operation (see "Primary-operation rows") are ever measured again.

Free quotient generators add columns but no power-relation rows and
require no torsion query. A free quotient splits because the intended
abutment category is that of abelian groups. Free generators occur in A
and D only, and neither needs a flat lift: an A generator is never a lower
generator, and a D generator enters every measured product with its marked
cocycle. The implementation keeps the named presentation
columns across all four stages, so later vectors can refer to any earlier
generator without losing its embedding.

## Primary-operation rows

A two-primary relation `m\widetilde q=t` whose target layer lies right
below the layer of `q` (A over B, B over C, C over D) is not measured in
the model. The class of its target-layer component is a primary operation
of the marked cocycle of `q`
([extension_cup_i_formulas.md](extension_cup_i_formulas.md), Sections 3–5),
which the engine evaluates on R with the cup-i products of its diagonals
and the Bockstein of the 0/1 lift (`backend.nativePrimary`), and projects
to the target layer with its E6 cell (`layers.<name>.cell.project`). With
`D=Sq^2+s\,Sq^1+\omega`, the operation of `d_2` on the row `q=-1`:

| Relation | Class of the target-layer component | Evaluation on R |
| --- | --- | --- |
| `a` over B | `D(z)` in `E_6^{j-2,-1}`, where `\beta_s z=(m/2)[a]` | `z=\rho u`, where `\delta_s u=m\,a` is solved over the integers |
| `b` over C | `\rho\beta_s(b)=(Sq^1+s)b` in `E_6^{j-1,-2}` | `(\delta_s\widetilde b)/2` modulo two, `\widetilde b` the 0/1 lift of `b` |
| `c` over D | an integral lift of `D(c)` in `E_6^{j+1,-4}` | `x-2w`, `x` the 0/1 lift of `D(c)` and `\delta_s w=(\delta_s x)/2` |

The classes are well defined in the E6 cells: changing `u` by a cocycle
`v` changes `D(z)` by `\operatorname{Dbar}[v]`, an incoming image; the
representatives and the cup-i structure change the cochains by
coboundaries; and the integral lift of `D(c)` is determined modulo
`2H^{j+1}(BG;Z_s)`. The cell projection quotients by all incoming images,
so the rows never need the gauge search of a component that meets them.
Only the target-layer coordinates are filled in. The row of an A or B
generator is the row of the measurement through its target layer and, like
it, the relation of a lift shifted below that layer (`truncatedBelow` is C
or D). The row of a C generator agrees with the measurement through D
modulo `2D`: it is the relation
`2(\widetilde c-h)=t` of a lift shifted by some `h` in D, its class in
`D/2D` is exact, and its witness records `shiftedLift="D"`. A later
relation measured through D with a nonzero coefficient on `c` therefore
measures the relation of `c` through D in the model, as for a truncated row
(previous section).

These rows need no model, flat lift, gauge comparison, bar transport or
universal value. On the Pin⁻ line of `C2` with `s` at degree six
(`2a=b`, `2b=c`, `2c=d`) every relation has its target layer right below
its generator, and `Z/16` is determined on R alone. Each row records
`witness.model="primary-R"` and the evaluated `witness.method`, and a
completed degree that built no model and has such a row records
`certificateLevel="primary-R"`, with the assumptions of the model whose
relations it reads (`gaugeCompletenessAssumed`, `abelianQuotientAssumed`).

When an ingredient is unavailable (a backend without the operations on R,
a layer without its E6 cell) or the evaluation fails, the relation is
measured in the model, and the degree result lists it in
`primaryFallbacks`. The odd relations keep their split rows. The relations
whose target layer lies two or three layers below the generator's (A over C
or D, B over D), and the measurements through D requested by a later
relation, are light rows (next section); their classes involve secondary
and tertiary operations
([extension_cup_i_formulas.md](extension_cup_i_formulas.md), Section 6).
The environment variable `FERMIONAHSS_NATIVE_RELATIONS=0` measures every
relation in the model; in a GAP session,
`KOAHSS_EXTENSION_RELATION_OVERRIDE.native` set to `false` or `true` takes
precedence until it is unbound.

## Light rows

A relation that neither its target layer nor a primary operation settles is
read from one fixed residue cochain of transported defining data
(`gap/extension_light.gi`, `python/extension_light.py`). Each prime part of a
degree is light as a whole: the two-primary part, and the three-primary part
in degrees five and six. The bar cochains are evaluated by the worker of the
transferred model; no flat lift, reflected product, gauge search or nonzero-A
D completion is formed.

**Defining data.** Every non-closed defining cochain is the primitive
`P(z;r)=\Lambda r+Hz` of a closed source `z`, where `\delta r=\Pi z` is solved
on R, `\Lambda` lifts R-cochains to the bar, `\Pi` pairs with the comparison
chains and `H` is the comparison homotopy. Then `\delta P(z;r)=z` holds
literally, so every branch flag of the stacking formulas is known by
construction; nothing tests a cochain for zero on R. The sources are the
model's own: `Q_D(\rho A)` for B, `f^\sharp(A,B)` for C, `Q_D(b)` for an A=0
state. Defining data are selected by linear algebra on the cohomology of R
with the actual classes: a C cochain of a B atom is corrected by the
`\widetilde D` image until its curvature class vanishes (flat-admissible),
the lower system of an A relation by the D image, the A=0 systems of
`\ker D` and the `\widetilde D` image, a relation gauge by the `\bar D`
image, and a gauge class by the D image after the page Tau map has
selected its E3 part. Rational cochains whose coboundary enters a residue
are paired on the chains of their own degree and the twisted coboundary of
R is applied (as for T on the pages). This identity requires normalized
cochains. The evaluator checks omitted degenerate faces before using the
chain-map identity, and degeneracies around the simplices consumed by a
primitive's homotopy. A nonzero value is an error, not a correction to the
residue. Eager and lazy comparison chains use the same tuple vertices.

**Residues.**

| Relation | Row |
| --- | --- |
| C over D, through D | `2D_c+\gamma_C(c,c)` for the marked C state `(c,D_c)`, `\delta_sD_c=-J(0,0,c)` |
| B over D, leading C part zero, k=4,5,6 | `\rho[T_b]=[\Phi_{\mathcal P}(b,P;c,\pi)]`, the page form of the B atom with its actual gauge phase; exact for `(b,c)` with some D completion |
| B over D, otherwise, k=2,3, and an atom an A-over-D row refers to | the exact residue `2D_b+\gamma_k(\hat b,\hat b)-D_C-J_{k-1}(u,y,\pi)-\gamma_k((0,0,t),(0,0,C))` against the pure-C reference of the leading part |
| A over C, m=2, k=5,6, leading part zero | `[\Phi']+[\omega a]`, `\Phi'=H_{k-1}(U,y'')+S_{page}(A,B_A)` |
| A over C, m=2, k=4, leading part zero, `\omega=0` | `[h^D(B_A,B_A)]` |
| A over C, otherwise | the exact `R_C=X_C+H_{k-1}(U,Y)+\beta(mA,X_B+B_0;0,B_0)+C_0` |
| A over D, k=4,5,6 | `\rho_m\Pi\Psi_{rel}`, with at most five diagonal phases (four in degree six) and, in degree six, the unary primitive with the 22 constants of `data/unary-gamma6-coefficients.json` |
| A over D at the prime three, k=5,6 | `t_D=2\cdot3^{e-1}Y`, `\rho_3Y=P^1_s\rho_3A` (the cube in degree five) |

A residue that should be an integral class modulo m is lifted by the
coefficient reduction map and read in the D cell; the Smith form of the
presentation then takes the quotient `D/(D\cap mH)`.

**Markings.** A C generator denotes its marked state `(c,D_c)`, a B generator
an integer combination of A=0 atoms `(0,b,c,D_b)` and other B generators.
A light row records what it is exact for: `truncatedBelow=fail` and
`shiftedWithin=fail` for the element, `shiftedWithin="D"` for its lower part
with some D completion, `truncatedBelow="D"` for a row in the absorption
frame. KOAHSS_ExtensionPrimeRows asks for the row of a lower generator again
through D (a complete request, answered for the current marking) when a later
row measured through the layer `lastMeasured` refers to it and
`lastMeasured` reaches the index of `shiftedWithin` (never for `fail`).

**Absorption.** When every two-primary or free D generator has order two,
the B rows with independent leading C parts are `(L,0)` (light-transport-proof
(27)–(28)): a change of the C generators by a homomorphism into `D[2]`
absorbs their D parts, and a dependent leading part is a formal sum of a
kernel atom and those pivots. Otherwise, or with
`FERMIONAHSS_LIGHT_ABSORPTION=0`, the D residues are retained. A kernel atom
may still use its page form, which fixes its lower part with some D
completion; an element reference requests its exact residue.

**The A plan.** At the first A request the A-over-D relations are ordered by
exponent and adapted (`a_i\mapsto a_i+2^{e_j-e_i}a_j`) until their leading B
vectors are independent. A relation with one leading B generator refers to
its atom; otherwise the B basis is changed unitriangularly so that a fresh
atom with cocycle `\sum_iL_ib_i` is a basis element. The reference atom is
re-marked by the pure-C state of the relation's relative C class, and a
relation without leading part refers to that pure-C state (light-transport-proof
(30)). The changed B generators are listed in `witness.light.remarked` and
recorded again. References to other B generators are live: a marking change
or precision upgrade invalidates every dependent row transitively and bumps
its version. Absorption is kept only if it is compatible with these
references; otherwise the B rows become exact residues.

Before accepting a prime part, a final frame audit checks the recorded
reference versions and precisions, the absorption frame, the target layers
in the final lower presentations, and the exact atom rows used by A-over-D
references. A failed audit restarts the whole prime part in the model;
rows from the two frames are never mixed.

Each light row records `witness.model="light-R"`, its `method` and
`witness.light` (row type, precision, shortcuts). A completed degree with a
light row and no measurement in the model has `certificateLevel="light-R"`
(see "The main calculation"). A resource limit leaves the degree
unresolved. Any other failure of a light row computes its prime part again
in the model and lists it in `lightFallbacks` with the failed step.
`degreeResult.heavyMeasurements` counts entries into the model measurement
path, including unsuccessful attempts. Zero means that no relation entered
that path; a worker used only for light residue pairings does not increase
the count.
`FERMIONAHSS_LIGHT_RELATIONS=0`, or
`KOAHSS_EXTENSION_RELATION_OVERRIDE.light:=false` in a GAP session, measures
these relations in the model. `FERMIONAHSS_NATIVE_RELATIONS=0` and
`FERMIONAHSS_LAYERED_RELATIONS=0` switch the light rows off as well, and
`FERMIONAHSS_PRIME_LOCAL=0` those of the prime three.

## Localization at the primes

The E6 layer generators are the independent generators of the cells
(`IndependentGeneratorsOfAbelianGroup`), so every one of them has prime-power
or infinite order (any other order is an error) and every relation
`m\widetilde q=t` belongs to one prime. The binary rows `q=-1,-2` are
two-groups, so odd-primary generators occur in the layers A and D only, and
every odd-primary relation belongs to A.

Right after the layers are read, the generators of the resolved layers
(those below the first unresolved layer) are grouped by the prime of their
order (`KOAHSS_ExtensionPrimeParts`). A free generator belongs to no prime:
it is a lower generator of the relations of every prime and carries no
relation itself. The rows are recorded prime by prime, two first, and for
each prime one generator at a time from D upward
(`KOAHSS_ExtensionPrimeRows`); a generator over the zero lower group is not
queried (the assembler records it with `kind="zero-lower-group"`). A
relation of the prime `p` over a nonzero lower group `H` has the zero row,
without a model or a lift, when it splits at `p` (below) and, otherwise,
when no lower generator has `p`-power order and none is free
(`witness.kind="zero-local-lower-group"`): `H` is then finite of order
prime to `m`, so `mH=H` and the class in `H/mH` vanishes. Every other
relation is measured: the target-layer shortcut (see "Layers and measured
relations") first; for `p=2`, the primary-operation row when the target
layer lies right below the generator's (see "Primary-operation rows"),
otherwise the complete model; for `p=3` in degrees five and six, the
three-local model.

A measured row keeps only its coordinates on the generators of its prime
and on the free generators; when this drops a nonzero coordinate, the full
row stays in `witness.fullLowerCoordinates`, with `witness.localizedAt`.
This is exact. The lower generators of odd order of a two-primary relation
(of B, C or A) lie in D, where they are torsion of odd order and vanish
after localization at two. An odd-primary relation belongs to A, is split or
three-local and has no B or C coordinate, since the three-local model has
no binary layers, so the D generators of other primes that it drops vanish
after localization at its prime as well. The rows of all primes therefore
present the group: its localization at each prime is presented by the rows
of that prime on the generators of that prime and the free generators, and
its rank is the number of free generators. The rows are assembled into one
presentation by `koAHSSExtensionFromLayers` with an oracle that replays
them (`KOAHSS_ExtensionReplayOracle`), which gives the Smith staging and
the fields of the next sections over all generators; the lower
presentation of a measurement is built the same way from the rows below it
(`KOAHSS_ExtensionLowerPresentation`). `degreeResult.primeParts` records,
for each prime, its generators and the `p`-primary part of the invariants.

The stacking model of a relation is chosen by the prime `p` of the order of
its generator (`KOAHSS_ExtensionRelationModel`):

- `p=2`: the complete transferred model, exactly as in the previous sections,
  or the primary-operation row of a relation whose target layer lies right
  below its generator's.
- `p\ge 5`: no measurement. `ko_{(p)}` is a sum of Adams summands of period
  `2(p-1)\ge 8`, so the rows `q=0` and `q=-4` lie in different summands and
  the window carries no k-invariant and no stacking correction at `p`; the
  relation is `m\widetilde q=0` in every degree, also over free lower
  generators, and the witness records `model="split"` with the certificate.
  A split relation needs neither the stacking model nor a flat lift; a
  completed degree with a split relation that built no model records
  `certificateLevel="prime-split"`.
- `p=3`: the three-local window is the two-stage tower of the rows `q=0`
  (layer A) and `q=-4` (layer D) with k-invariant `2\beta_3P^1\rho_3`; the
  potential is the coefficient-two phase `\Omega(A)=\tfrac23\,\mathrm{lift}(P^1_s\rho_3A)`
  of the [tertiary operations](tertiary_operations.md). In input degrees zero
  and one `P^1=0`, so below package degree five the tower splits and the
  relation is `m\widetilde q=0` (`model="split"`). In package degree five,
  `P^1_s\rho_3A=\rho_3(A\cup A\cup A)` and `\Omega` differs from
  `\Omega'(A)=\tfrac23\,A\cup A\cup A` (transported cup products) by an
  integral cochain, an isomorphism of stacking models. In the model of
  `\Omega'` the curvature of a cocycle vanishes, the stacking correction is
  the integral polynomial

  \[
   \gamma(A,A')=-2\,(A\cup A\cup A'+A\cup A'\cup A'),
  \]

  whose difference from `\Omega'(A)+\Omega'(A')-\Omega'(A+A')` is
  `\tfrac23\,\delta_s\Xi` with the cup-one expression
  `\Xi=2A(A\cup_1A')+(A\cup_1A')A+2(A\cup_1A')A'+A'(A\cup_1A')`
  (the Hirsch identity `A'A-AA'=\delta_s(A\cup_1A')` word by word), and the
  boundary state of a gauge `(u,w)` is `(\delta_su,\delta_sw)`. The relation
  of a generator of order `3^a` with `3^aA=\delta_su` and local flat lift
  `(A,0,0,0)` is

  \[
   x_D=2\cdot 3^{a-1}\,A\cup A\cup A
       +\tfrac23\,\delta_s\Bigl[\sum_{j=1}^{3^a-1}\Xi(jA,A)-u\cup\delta_su\cup\delta_su\Bigr]
       \pmod{3^aH_D},
  \]

  whose last term is a Bockstein class and is computed, not dropped.
- **k = 6 (|A| = 3):** `P^1_s\rho_3A` is the nineteen-term cyclic-diagonal
  formula (`mod3_power.reduced_power_terms`, the evaluation of
  `a\otimes a\otimes a` on `D_2`), and the model uses the potential
  `\Omega(A)=\tfrac23\,L(P^1_s\rho_3A)` itself, `L` the 0,1,2 lift. The
  curvature of a cocycle is `J(A)=\delta_s\Omega(A)=2\beta_3P^1_s\rho_3A`,
  the three-primary d5 term, so the flat lift `(A,0,0,D)` solves
  `\delta_sD=-J(A)`, which is possible exactly when the class survives to
  E6. The stacking correction and the gauge boundary use the two natural
  primitives of the reduced power modulo three that the cyclic diagonals
  provide (`mod3_power.cross_effect_primitive`,
  `mod3_power.coboundary_primitive`):

  \[
   \gamma(A,A')=\tfrac23\bigl[L(P^1a)+L(P^1a')-L(P^1(a+a'))+\delta_sL(\varphi(a,a'))\bigr],
   \qquad D_u=-\tfrac23\bigl[L(P^1\rho_3\delta_su)-\delta_sL(\chi(u))\bigr],
  \]

  where `\varphi(a,a')` is the sum of the mixed words
  `a\otimes a\otimes a'+a\otimes a'\otimes a'` evaluated on
  `(\rho^2+2\rho)D_3`, and `\chi(u)` is `u\otimes\delta u\otimes\delta u` on
  `D_2` minus `u\otimes\delta u\otimes u-u\otimes u\otimes\delta u` on `D_1`
  plus `u\otimes u\otimes u` on `D_0`. The first is a primitive of the cross
  effect because the six mixed words are the norm `N=1+\rho+\rho^2` of the
  two words on `D_2`, `N=(\rho-1)(\rho^2+2\rho)` modulo three and
  `(\rho-1)D_2=dD_3+D_3d` on cocycle tensors; the second is a primitive of
  `P^1(\delta u)` because `(\delta u)^{\otimes3}` is the tensor coboundary of
  `u\otimes\delta u\otimes\delta u`, `dD_2-D_2d=ND_1`, the norm of that
  tensor is three times itself plus the tensor coboundary of
  `u\otimes\delta u\otimes u-u\otimes u\otimes\delta u`, and `(\rho-1)` of
  the latter is the tensor coboundary of `u^{\otimes3}` modulo three, whose
  value on `D_0` is a coboundary. Both brackets are divisible by three, which
  the model checks on every value; the boundary state of a gauge `(u,w)` is
  `(\delta_su,\delta_sw+D_u)`, and a gauge acts by the product with its
  boundary state. Associativity and commutativity hold up to signed
  D-coboundaries.

The model (`python/extension_three_local.py`, requests with `prime=3` on the
same worker process, `model.primeLocal(3)` in GAP) evaluates `\gamma` on
the resolution through the comparison lift and projection of the
transferred model, so associativity and commutativity hold up to signed
D-coboundaries as in the complete model; the flat lift of a three-primary
A generator is solved in the two-layer model when its relation is measured
(`(A,0,0,0)` in degree five, `(A,0,0,D)` in degree six), the power is
reduced in the A layer by coboundaries and in the D layer by the marked D
cocycles and coboundaries, the binary layers are absent, and the relation
keeps its exact gauge comparison (`witness.model="three-local"`,
`witness.prime=3`). If the three-local reduction does not express the D
residual through the marked D cocycles and ordinary coboundaries, the
relation stays unresolved (`code="three-local-reduction"`). The complete
four-layer measurement is not used for it: it writes the relation in the
lifts of the binary generators, whose three-local contribution needs their
complete, untruncated relations.

Lifts of one prime's model never enter a measurement of another model: a
three-local lift enters only its own relation, and the lower products of
both models consist of complete lifts and marked cocycles of lower
generators (in the three-local model, marked D cocycles only). The
per-prime summary of a degree is `degreeResult.primes` (prime, number of
relations over a nonzero lower group, and their models), and
`FERMIONAHSS_PRIME_LOCAL=0` restores the complete measurement for every
prime: every relation uses the complete model, a row keeps all its
coordinates, and the zero rows of a lower group without generators of the
prime are not used.

For cyclic groups the window is `F^2/F^{10}` of the skeletal filtration of
`ko^2(BZ/n)`, which is the `I`-adic filtration of the anti-invariant part
of the representation ring under complex conjugation, and this gives the
three-local results independently: `Z/9` for `Z/3` (`3a=2d` modulo 3, the
relation `3A=-16D` of the model), `Z/3\oplus Z/27` for `Z/9` (`9a=6d`
modulo 9), `Z/9\oplus Z/81` for `Z/27`, and the split groups
`Z/5\oplus Z/5`, `Z/7\oplus Z/7`, `Z/25\oplus Z/25` at the primes five and
seven. The complete four-layer formulas are not defined on every
three-torsion input (their half-lift carries fail for the generator of
`Z/9` and for the sign-twisted `S_3`), so for these groups the three-local
model is the only measurement.

Degree six is checked by suspension: the degree-six window of `G\times Z`
is the suspension of the degree-five window of `G` (`ko^3(\Sigma BG)=ko^2(BG)`,
the layers being `H^3(G\times Z)\supset H^2(G)\otimes H^1(Z)` and
`H^7\supset H^6(G)\otimes H^1(Z)`), and on the product resolution
(`ResolutionDirectProduct` of the cyclic and the `Z` resolution) the
three-local model gives `Z/9` for `Z/3\times Z` (`3A=14D`) and
`Z/3\oplus Z/27` for `Z/9\times Z` (`9A=42D`, that is `9a=6d` modulo nine),
the degree-five values of `Z/3` and `Z/9`. Finite groups rarely have
three-torsion survivors in the degree-six A layer: for elementary abelian
groups d5 does not vanish on the Tor classes of `H^3`.

## Smith coordinates and filtration maps

In degrees 1–6 `relationMatrix` is the one presentation of the rows of all
primes ("Localization at the primes"). The returned Smith data satisfy
`smith.U * relationMatrix * smith.V = smith.S`. Coordinates are row
vectors. A vector `x` in the named
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
lifts. The production witnesses contain concrete cochain data. The
abstract assembler does not reconstruct flat cochain representatives for
an arbitrary supplied relation oracle.

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
greater than one or zero for free generators. The only accepted options
record is empty.

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
`koFull` records its rows prime by prime first and then calls the assembler
with an oracle that replays them.

## Complete flat representatives and gauge comparisons

In degrees 1–6, the production engine uses the transferred four-cochain
curvature and product in the retained resolution basis. The calibrated
formulas are evaluated lazily through sparse bar transport. The leading
E6 cochains keep their native coordinates. In degrees 4–6 the D-layer
product of two lower-legal states is `closed_ab_upper.gamma` with the
half-lift carry of the C residuals and the K rephasing, without the pure-C
normalization of the production correction: that normalization is
\(\delta_s\) of a phase correction with an explicit integral primitive
(`pure_c_normalization.integer_phase_correction`), hence an integral
coboundary, and a D-gauge with that primitive carries one product to the
other. The stacking group on gauge classes is unchanged, the second phase
evaluation at A=B=0 is saved, and native D representatives differ from the
finite-section reference model's by D-coboundaries, which the D stage of
the gauge comparison absorbs.

In degrees 5 and 6 the A=0 sector uses the direct formulas of
`a0_high_gamma` and `compatible_sector` (`extension_native_upper`). For a
lower-legal state with A=0 on the resolution the D-layer curvature
`closed_ab_upper.J(0,B,C)` equals, as a cochain,
\(\delta_s[\Omega_{a0}(B,C)+\mathsf h(t\smile_{k-1}Q_D(B))]-\mathsf h(E(t))\)
with \(\Omega_{a0}\) the phase of `HigherA0Stacking`, \(t\) the binary C
residual \(\rho(\delta C+Q_D(B))\), and `pure_c_g(C)` when B=0 as well:
the source splitting of A=0 vanishes, so `production_phase(n,0,B,C)` is
\(\Omega_{a0}(B,C)\) for closed B, \(f^\sharp(0,B)=Q_D(B)\), and every
comparison gauge and universal primitive vanishes at A=0. No universal
simplex is registered for such a curvature. For two fully legal states
with A=B=0 the product correction `closed_ab_upper.gamma` equals
`pure_c_gamma(C,C')` plus \(\delta_s\) of the explicit integral cochain
\(-I\,\gamma_C(C\ell,C'\ell)+I\,\widetilde{\mathrm{pol}(\delta(C\ell),\delta(C'\ell))}-[s\smile C]\,[s\smile C']\)
(\(I\) the right prism, \(\ell\) the interval coordinate,
\(\gamma_C\) = `pure_c_gamma`, the last term the carry of the half lifts),
an integral coboundary absorbed by a D-gauge; the native product uses
`pure_c_gamma` there, keeping the half-lift carry and the K change. With a
nonzero B the difference of the two corrections is \(\delta_s\) of a
rational cochain whose integrality is only tested, so that sector keeps
`closed_ab_upper.gamma`. In degree 6 a product of which exactly one factor
has a nonzero A layer is evaluated: its universal pair primitive vanishes
identically, since every term of the pair contraction keeps a zero fiber
and the relative small basis has no word of one fiber alone.

Flat lifts are solved when a measurement first needs them. The measured
generator gets its complete lift, its lift through the target layer when
that layer lies right below the generator's ("Layers and measured
relations"), or its two-layer lift at the prime three. A lower generator
gets its complete lift only when it enters the compared product with a
nonzero coefficient in a layer above the last measured layer (the target
layer, or D in a complete or three-local measurement); in the last measured
layer it enters with its marked cocycle, the state with that one layer,
since the leading component of a flat lift is the marked cochain itself
and the products through that layer read nothing else. The rows of a layer
reduction are these marked cochains, so a zero coefficient needs no lift,
and no D generator and no free generator ever needs one. The complete and
three-local lifts that were needed are kept in `layers.<name>.fullLifts`,
a lift through the target layer in the witness of its relation
(`witness.flatLift`); a needed lift that fails leaves its relation
unresolved with the lift's reason.

For an A generator, the complete lift solves

\[
 \delta_s A=0,\qquad \delta B+P_R(A)=0,\qquad
 \delta C+\tau_R(A;B)=0,\qquad \delta_sD+J_R(A,B,C)=0.
\]

Here the subscript R denotes the transferred nonlinear terms, including
the homotopy corrections; it does not mean direct substitution into the
bar formula. The B and C equations are solved over F2 and the D equation over Z.
The solver searches the affine B and C solution families lazily: if a
chosen B does not admit C, or a chosen C does not admit D, it changes
that defining choice and retries the dependent equations. It evaluates
the actual full differential to verify the resulting tuple is flat.
B- and C-layer generators undergo the same lower-equation procedure.
A lift through the layer `T` solves the equations through `T` only, with
the layer-limited differential, and leaves the layers below `T` zero.
The solver stores the chosen primitives, their adjustments, the exact
equation data and the final zero-curvature witness. No universal helper
is replaced by a locally solved formula.

Each `fullLifts[i].state` is an immutable tuple `(A,B,C,D)` in the retained
resolution basis. The exact same tuple is reused in all powers, lower reductions
and basis comparisons. A complete lift's lower components are not reset to
zero or reconstructed from their cohomology classes in a later query; only
in the last measured layer, where nothing below it is read, does a
generator enter with its marked cocycle alone.

For a torsion generator of order `m`, the oracle stacks its lift to obtain
its measured power. It reduces the result using the lifts and marked
cocycles of the lower generators described above and retains all integral
carries and boundary data. The resulting coordinates refer to the named
D/C/B/A columns of the lower presentation, and the rows of all primes are
combined in one relation matrix.

An ordered reduction alone does not justify rearranging cochain products.
The engine therefore constructs the canonical lower product in its
recorded generator order and solves the exact comparison

\[
 \text{measured power}=\operatorname{act}_R(g,\text{canonical lower product}).
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
first as candidate directions; the search also uses the complete native
cochain kernel, including ordinary coboundary directions and their nonlinear
carries. A bounded integer-kernel search is not reported as exhaustive.
The leading primitive and subsequent B/C/D choices are solved consistently;
the comparison retains `g` and its native action, verifies flatness of
the action result, and checks the displayed equality component by component.
The action reflects the full boundary of the embedded gauge, not the
truncated native curvature of a possibly nonflat gauge. It does not
assume strict associativity or commute factors at the cochain level.
Failure to find this exact witness within the search bounds remains
unresolved.

If ordinary layer reduction fails at D, the reducer projects the remaining
D cocycle into the retained E6 cell to propose coordinates in the marked
D basis. This accounts for incoming differentials that ordinary
coboundaries alone do not remove. It then compares the **original full
stacked tuple** with the canonical product of the lower generators with
nonzero proposed coordinates (complete lifts, including their earlier
carries, and marked D cocycles), using the staged gauge search above.
The E6 projection alone never certifies a stacking relation.

For other failed ordinary reductions, the alternative search tries at most 32
marked finite lower normal forms, stage by stage, with the same
representatives. It does not guess coordinates for a free lower factor.
Accepted relations retain `canonicalComparison.winningStage`,
`attemptedStages` and the full gauge. Unresolved relation queries retain
`pendingRelation`, including attempted comparisons and residual cochains,
instead of claiming that a particular missing gauge necessarily exists.

A higher-degree result is returned as computed once every relation has
its exact comparison. Commutativity and associativity of stacking on gauge
classes are assumed rather than audited: products of the marked finite
normal forms are not enumerated. Free quotient coordinates split in the
intended abelian abutment category.

Degree six is the cutoff of the all-cochain differential of
[all_cochain_differential.md](all_cochain_differential.md). Its lower
three components are the formulas of degrees three to five with input
degree three: \(\delta_sA\), \(\delta B+P_6(A)\) and
\(\delta C+\tau'_3(A;B)\) on the legal lower locus \(L_6\), with the
prism \(H_6\) outside it. On \(L_6\) the D-layer differential is
\(J_6\), the integral cochain (A5) of the fixed degree-three phase
\(\Omega_6\) with the residual correction, for arbitrary C; on a full
defining system it is \(T_3\). The D-layer stacking correction of two
lower-legal triples is the legal production \(\gamma_6\) with the common
pure-C normalization and the rephasing of the successor, exactly as the
finite-section reference model evaluates it on legal data. The native
engine evaluates these D-layer terms only on lower-legal triples: native
curvature stops at the first nonzero layer, and products and gauge actions
receive flat states and differential images, whose lower pairs are legal.
The section retraction \(\mathcal R_6\) that the note uses outside
\(L_6\) is not evaluated on the resolution; a request for that branch
leaves the degree unresolved with that reason. The degree-six section of
the complete-bar reference model is used only by the tests.

The stacking correction \(\gamma_6\) of two states whose A layers are
not both zero on the resolution evaluates the universal pair source of
`production_gamma6`. Each universal term of that source contains three V3
source contractions on eight-vertex universal simplices, and the pair
contraction itself has thousands of terms per output simplex, so this one
correction costs orders of magnitude more than every degree-five term. The
worker therefore refuses it and the degree stays unresolved with that
reason, unless the environment variable `FERMIONAHSS_DEGREE_SIX_A_STACKING=1`
asks for the evaluation regardless of its running time. The degree-six
differential \(J_6\) with nonzero A, and every degree-six term of states
with A zero (the B, C and D layers), are evaluated.

The runtime bundles the selected stacking sources under
[python/stacking_model](../python/stacking_model/) with their
[upstream provenance](../python/stacking_model/provenance.json). Helpers that
are identical to the package kernel (`cochain_tools.py`, `phase_eval.py`) are
imported from it; the bundled hashes are not verified at runtime.
Loading the package does not require the separate research workspace.
No expected classification table is consulted at runtime.

## Low degrees and resource limits

Degrees -1 and 0 have only the D layer: the group is that layer and no
relation is measured, so, as in every degree without a relation to
measure, no stacking model is built. In degrees one and two every
two-primary relation with an adjacent target (C over D, and B over C) is
read from a primary operation on R (see "Primary-operation rows"); a
relation measured in the model there (B over D, a fallback, or a
measurement through D requested by a later relation) uses the native
engine with the differentials and products of
[low_degree_stacking.md](low_degree_stacking.md), for arbitrary valid
`s` and `omega`. A degree-one state is `(C,D)` and its gauges are
degree-zero states `(D)`; a degree-two state is `(B,C,D)` and its gauges
are degree-one states. Flat lifts, relation measurements, gauge
comparisons, the comparison complex and the zero test are those of degrees
3–6. No expected classification table is consulted at runtime. The
abstract assembler's generality does not remove a production cochain
requirement.

The native implementation has explicit resource bounds:

| Work | Bound |
| --- | --- |
| Distinct g-support simplices per degree, over the comparison chains verified in that degree | 8192 |
| Processed sparse transport terms: g and f terms of the verified comparison chains, and normalized-homotopy expansion terms | 2,000,000 each |
| Flat-lift affine search | 4096 distinct differential evaluations by default |
| Gauge-comparison affine search | 4096 equation evaluations and 64 leading choices shared across requested stages; integral kernel coefficients initially bounded by absolute value 1 |
| Alternative lower-coordinate search | 32 marked finite normal forms; 4096 equation evaluations shared across stages and candidates |
| Degree-six stacking correction with a nonzero A layer | refused, unresolved, unless `FERMIONAHSS_DEGREE_SIX_A_STACKING=1` |

These are implementation limits, not additional tuning arguments to `koFull`
or `koFull_batch`.
Repeated integer solves reuse exact Smith preparations, bounded by eight
entries and two million retained matrix cells including transforms.
States and matrices use R; complete-bar dimension and matrix limits do
not gate native completion. The fixed formulas still use sparse
transport to the comparison complex, and branch predicates are tested on
the basis of R. See [resolution extensions](resolution-extensions.md).

The native model assumes gauge completeness. The
[transfer note](transfer.md) explains why a strict retraction and sound
individual relations do not establish that claim.
The assumption is recorded in results; the relation checks do not turn
it into a proof.

A nonlinear evaluation can still be costly. Reaching a resource bound or
failing to find a gauge within the bounded search is not evidence that an
extension splits. Exact identity failures and inconsistent coordinate data
are errors; unavailable witnesses leave the corresponding result unresolved.

The [paper comparisons](extension-paper-comparisons.md) use spatial
dimension `d`, related to package degree by `j=d+1`. Paper dimension `d`
is therefore `koFull(...,d+1).invariants`, or `full.invariants[d+3]` for a
`koFull_batch` result. Compare complete invariant
lists with the same twists and phase convention; neither E6 layer-order
products nor unresolved outputs establish a full-group match.
