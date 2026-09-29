# Mathematical status and provenance

fermionAHSS uses fixed formulas on the normalized group bar resolution,
evaluated through resolution comparison maps. Its page calculation covers rows `q=-4,-3,-2,-1,0`
through E6 and physical cutoff `p+q+3 <= k <= 6`. E6 is terminal within
this strip; it is not automatically the full ko E-infinity page. Other ko
rows remain outside the implementation. The extension assembler and
native stacking engine described below do not establish a
general abutment identification. Detailed results retain
`certified_ko: false`.

## Extension implementation scope

`koFull(group,s,omega,k)` computes E6 and attempts the four-layer
extension in package degree `k`; `koFull_batch(group,s,omega,k)` reuses
one E6 calculation for every degree from -1 through `k`. Both also accept
an already computed detailed E6 result and reuse its exact context. They
do not change the row window or the calibrated page differentials. See
the [extension API](extensions.md).

A supplied integral HAP resolution can replace the group argument.
Native states and solves on that resolution are the only extension
mode in degrees 1–6; there is no model-selection option. The comparison
(the group bar, or else the cell complex on the generators of R) must pass
the strict retraction check, and the sparse homotopy is normalized.
The calculation assumes gauge completeness: native gauge equivalence
captures the bar equivalence relation. Results record
`gaugeCompletenessAssumed=true`; this is an assumption, not a runtime proof.
There is no complete-bar certification or automatic complete-bar fallback.
See [resolution extensions](resolution-extensions.md) for the exact scope.

The abstract assembler accepts arbitrary finitely generated abelian
A/B/C/D layers with a relation-vector oracle. For each torsion quotient
generator of order `m`, the oracle records the full coordinates of
`m*lift(generator)` in one identified lower presentation. The assembler
combines these rows into a single integer presentation, retains its Smith
basis transformations, and carries the filtration inclusions and
quotients to subsequent stages. Equal and independent nonzero lower
images remain distinguishable.

For degrees 1–6, the production engine retains the marked E6
representatives in the supplied resolution basis. It
solves and retains immutable full `(A,B,C,D)` lifts for every generator,
including free A generators. Binary B/C defining choices are changed when
needed to solve later equations, and D is solved integrally. The actual
nonlinear differential verifies flatness of each complete tuple. The
selected `xtimes` then measures powers of those same tuples; lower
reductions reuse the recorded B/C/D basis and all integral carries. In
degrees 4–6 the native D-layer product omits the pure-C normalization of
the production correction, an integral coboundary, so its D
representatives differ from the reference model's by D-gauges while the
stacking group on gauge classes is unchanged.

Each measured relation must admit an exact native gauge comparison
with its canonical lower product in its measured layers. A relation is
measured only through its target layer, the lowest layer of the lower
presentation with a generator outside `m*H` (`m` the order of the
generator, `H` the lower group): the components below it lie in `m*H` and
cannot change the class of the extension in `Ext(Z/m,H)=H/mH`, which the
recorded certificate exhibits, and they are left unmeasured. The witness
retains a native gauge tuple, verifies flatness of its action result
through the measured layers, and checks the action equality component by
component there. Commutativity and associativity of
stacking on gauge classes are assumed, so the measured relations determine
the group; no finite multiplication table of normal forms is audited. Free
quotients split in the intended abelian abutment category.

The measurement is localized at the prime of the generator's order. At the
primes five and above no k-invariant of `ko` links the rows `q=0` and `q=-4`
of the window (they lie in different Adams summands), so those relations are
`m*g=0` without a measurement and a degree whose relations all split needs
no stacking model; for cyclic groups this agrees with the `I`-adic
filtration of the representation ring (`Z/5+Z/5`, `Z/7+Z/7`, `Z/25+Z/25` in
degree five). Three-primary relations of degree five are measured in the
two-layer model of the rows `q=0` and `q=-4` with the coefficient-two phase
`(2/3) lift(P^1_s rho_3 A)`, in degree five the cube modulo three, whose
stacking correction is the polynomial `-2(AAA'+AA'A')` and whose gauge
boundaries are `(delta_s u, delta_s w)`; the relation keeps its exact gauge
comparison in that model. It reproduces the representation-ring values
`Z/9` for `Z/3`, `Z/3+Z/27` for `Z/9`, and `Z/9` for the sign-twisted
`S_3`, `Z/6` and `A_4`, where the complete four-layer formulas raise a
non-integrality error on the three-torsion generator. Three-primary
relations of degree six, whose reduced power is the nineteen-term formula,
stay on the complete path (refused for two nonzero A layers), and
`FERMIONAHSS_PRIME_LOCAL=0` restores the complete measurement at every
prime. See [extensions](extensions.md), "Localization at the primes".

Gauge comparison tries D-only, C/D, B/C/D and A/B/C/D support in that
order. This includes higher-layer boundaries responsible for incoming E6
differentials. If a final D residual is not generated by the surviving D
lifts and ordinary coboundaries, its E6 projection proposes lower
coordinates, which still require an exact comparison of the original
full stacked tuple. Native cohomology lifts prioritize the search but do
not replace the complete cochain kernel. Unsuccessful bounded searches
remain unresolved, with their diagnostics retained.

Degree six uses the same engine with the degree-six D-layer terms of the
[all-cochain note](all_cochain_differential.md): \(J_6\) and the legal
production \(\gamma_6\) with the pure-C rephasing, evaluated on the
legal lower locus, which is the only locus the native operations reach.
The section retraction used by the note outside that locus is not
evaluated on the resolution; such a request leaves the degree unresolved.
The stacking correction of degree-six states with a nonzero A layer is
refused for its running time unless `FERMIONAHSS_DEGREE_SIX_A_STACKING=1` is
set, again leaving the degree unresolved; see [extensions.md](extensions.md).
The degree-six section of the complete-bar model serves only as a test
reference and is not a runtime fallback. The E6 page calculation covers
degrees up to six.

Degrees one and two use the same engine with the degree-one and
degree-two formulas of [low_degree_stacking.md](low_degree_stacking.md),
for arbitrary valid `s` and `omega`. The stacking formulas are bundled
with their [upstream provenance](../python/stacking_model/provenance.json);
runtime loading does not require the separate stacking research workspace.
Helpers identical to the package kernel are imported from it, and the
bundled hashes are not verified at runtime.

Degrees -1 and 0 have only D, so no relation is measured. Incomplete page data
and exhausted resource bounds remain unresolved. Native transfer
requires a comparison passing `fg=id` over the integral group ring; the
cell comparison provides one for any resolution whose generators have
primitive boundaries, including infinite groups and several degree-zero
generators. Branch flags test vanishing on the basis of R. Sparse
transport has term bounds. Flat-lift searches default to 4096 differential
evaluations. Gauge comparisons have additional finite search bounds
documented in [extensions.md](extensions.md). Reaching a bound leaves the
degree unresolved; it is not substituted for missing cochain evidence.

The [paper comparisons](extension-paper-comparisons.md) distinguish
computed groups, unresolved cases and literature expectations. Their
spatial dimension `d` corresponds to package degree `d+1`. Finite matches
support these cases and do not establish all-degree naturality or a
complete ko/SPT identification.

## Fixed secondary family

The runtime family is `chi7_tail` in input degrees zero through seven,
with epsilon `(1,0,0)` and eta `(1,0,1)`. Its word identities, right-product
suspension, and Thom calibration are the mathematical inputs to the
implementation. The exact common integral lift satisfies
`rho M'=tau'` and `d_s M'=2 psi'`. The square-zero assertion `Dtilde[Tau]=0`
means an integral coboundary, not a literally zero cochain.

The source-workspace provenance is
`note/extra/chi_suspension_degree7/`, including `tail_words.json`, and
`note/extra/secondary_operations_psi_note.tex`, appendix equations
`Th-squarezero` and theorem `global-squarezero`. These proof files are not
bundled. The packaged [ANF data](../data/chi-calibrated-degree7-anf.g) and
[helper reference](universal_helpers.md) fully specify runtime evaluation.
Other chi families, such as the head-suspension family and constrained
degree-nine references, have different normalizations and cannot be
interchanged with `chi7_tail`. Local residual solves and diagnostics based
on injectivity of `Dtilde` do not select production maps.

## Fixed tertiary family

The universal helper R and final T are implemented only in input degrees
zero through three. The general signed R0, R1, R2, and R3 family is used;
the compact even-input family is not silently spliced into it. R0 retains
full and quarter-input carries. R2 is `R2sharp`, which includes the term
`-A^cup3/4`; the R3 coefficient prescriptions use the source `V2fin`,
without that term, in `C0` and xi. Applying the cubic correction twice or
absorbing it into the binary L6 helper changes the operation.

The source supplements establish low-degree first-b descent, additivity,
and the specified suspension comparisons. These are separate claims from
residual solvability, exact boundary, and absolute normalization. The
rank-six `BPSO(6)=BPU(4)` comparison fixes the `q(omega)A` correction to
zero for the formulas in [tertiary_operations.md](tertiary_operations.md).
The signed Euler calculation gives `mu_R=0`.
Thus final T adds `2 beta_3,s P^1_s rho_3,s`. On the page it can
contribute only in input degree three; in input degree two the reduced
power is the cube modulo three, whose Bockstein vanishes on cocycles, but
the phase itself enters the degree-five stacking correction, where its
failure of additivity is the three-primary extension (the untwisted `Z/3`
in package degree five is `Z/9`, not `Z/3 + Z/3`).

The nonzero finite data are essential: tertiary low-selector vector
\(\boldsymbol{\zeta}=(\zeta_1,\zeta_2,\zeta_3)=(0,1,0)\), V2 source
values `(0,3/4,0,3/4,1/4)`, R3 selectors `(1,0,1,1,1)`, xi `3/4`,
and suspension periods `(3/4,1/4,0,1/2,1/2)`. The R3 source calculation
has a 46-by-94 matrix and all 63 cycle-obstruction checks pass; source
denominators reach sixteen. The final rational phase may also have a
factor of three. See [tertiary formulas](tertiary_operations.md) and
[universal helpers](universal_helpers.md).

No all-degree universal R extension is asserted. Degree-six
quarter- and eighth-denominator obstructions concern a different range and
do not contradict the implemented low-degree family.

GAP solves only the permitted defining cochains b and c. R is evaluated
from the fixed universal operators and calibrated source values, never
chosen by a local residual solve. A solution using the last nullhomotopy,
such as R=-U(c), would not define the required universal R(A,b).

## Sources

The formulas come from the following research notes, which are not
bundled:

- General R formulas and suspension arguments: `R0_fixed.md`, `r3.md`,
  and `r2_cubic_normalization.md`.
- `notes/extra/tertiary_T_degree3/`: low-degree T descent, additivity,
  rank-six calibration, and normalization arguments.
- `notes/extra/tertiary_T_degree3/mu_verification/`: exact Euler calibration
  with `Xi(D)=V2(D)=13/4`, `kappa6 V2(Ss)(Tor)=3/2`, and
  `Uraw(Tor)=-1/2`; this is distinct from the production R3 computation.

The universal R3 source data in [r3_source.json](../python/r3_source.json)
record the hashes of the modules that produced them; these hashes are
not checked when the data are loaded.

Finite executable checks support the stacking engine and arithmetic. They do not
replace the universal arguments or establish an all-degree theorem.
