# Implemented formulas and calibration data

These documents describe the fixed bar formulas and resolution-based
implementation of fermionAHSS.

| Document | Contents |
| --- | --- |
| [conventions.md](conventions.md) | Signed coefficients, interval-cut words, integral signs, binary lifts, primary maps, bar comparison, and defining-cochain corrections |
| [secondary_operations.md](secondary_operations.md) | Complete Tau and Psi formulas, exact common lift, domains, and indeterminacy |
| [tertiary_operations.md](tertiary_operations.md) | Implemented R0–R3 phases, the cubic term of R2, and all 19 prime-three terms |
| [universal_helpers.md](universal_helpers.md) | Chi ANF decoding and degree table, zeta words, universal contractors, finite source tables, and every normalization selector |
| [universal_value_growth.md](universal_value_growth.md) | Growth and periodicity of the two stored universal source values: exact decomposition of the pair source, closed form of the degree-one contraction, affine growth of `V_1` on residue classes, and the finite data that determine them |
| [backends.md](backends.md) | HAP, cochain, and page interfaces; direct-operation audits and exact quotient conventions |
| [extensions.md](extensions.md) | Detailed AHSS results, `koFull` and `koFull_batch`, retained flat tuples, exact gauge comparisons, common-coordinate Smith presentations, and the abelian-quotient assumption |
| [transfer.md](transfer.md) | Resolution transfer: the retraction check, the normalized homotopy, conditional flatness and reflection, a counterexample to gauge completeness, and the preflight |
| [resolution-extensions.md](resolution-extensions.md) | Native extension API, sparse normalized transport, gauge-completeness assumption and native completion checks |
| [extension-paper-comparisons.md](extension-paper-comparisons.md) | Full-group literature fixtures, dimension and twist conventions, and the measured relations of the C2 cases |
| [mathematical-status.md](mathematical-status.md) | Implemented range, mathematical assumptions, sources, and unresolved scope |
| [dimension_indexed_differentials.md](dimension_indexed_differentials.md) | Explicit nonlinear differential through k=6 on the first-two-layer domain, the residual correction, and a literal square-zero proof |
| [low_degree_stacking.md](low_degree_stacking.md) | Stacking products in package degrees one and two (gamma1, beta2, gamma2) and their use in `koFull` |
| [all_cochain_differential.md](all_cochain_differential.md) | Piecewise extension to arbitrary cochains with first component δ_sA, complete degree table, exact square-zero proof, and the defining-system obstructions; not a natural local extension |

The large chi word lists are encoded exactly by
[data/chi-calibrated-degree7-anf.g](../data/chi-calibrated-degree7-anf.g)
and the decoding rule in the helper reference. This is the complete
runtime coefficient data; the millions of word summands are not needed
separately.

## Names and fixed coefficients

The notes and GAP use the same differential names throughout.

| Differential | Rows | Function |
| --- | --- | --- |
| d2 | 0 to -1 | `Dbar` |
| d2 | -1 to -2 | `D` |
| d3 | -2 to -4 | `Dtilde` |
| d3 | 0 to -2 | `Tau` |
| d4 | -1 to -4 | `Psi` |
| d5 | 0 to -4 | `T` |

The primary functions are defined in [conventions.md](conventions.md).
The formulas for `Tau` and `Psi` are in
[secondary_operations.md](secondary_operations.md), and those for `T`
are in [tertiary_operations.md](tertiary_operations.md).

| Datum | Value |
| --- | --- |
| Chi family | `chi7_tail` |
| Secondary epsilon, eta | `(1,0,0)`, `(1,0,1)` |
| Tertiary low-selector vector \(\boldsymbol{\zeta}=(\zeta_1,\zeta_2,\zeta_3)\) | `(0,1,0)` |
| R2 | `R2sharp`, including the term `-A^cup3/4` |
| Final T rank correction and mu_R | `0`, `0` |
| R3 selectors c4,cN,cO,cM,epsilon_c | `(1,0,1,1,1)` |
| R3 xi | `3/4` |
| V2 source values | `(0,3/4,0,3/4,1/4)` |
| R3 suspension periods | `(3/4,1/4,0,1/2,1/2)` |
| Prime-three coefficient | `2`; contributes only in input degree three here |

The entries of \(\boldsymbol{\zeta}\) select the R1 rank normalization and
the two R2 suspension terms, respectively. This vector is distinct from
the secondary epsilon and eta vectors and from the cochain helpers
\(\zeta_{i,n}\).

## Runtime source map

| Source | Role |
| --- | --- |
| [natural_bar.gi](../gap/natural_bar.gi) | Integral equivariant bar comparison and homotopy |
| [natural_words.gi](../gap/natural_words.gi) | Interval cuts, integral cup signs, and compiled chi evaluation |
| [natural_secondary.gi](../gap/natural_secondary.gi) | Complete lower formulas and matched lift |
| [natural_tertiary.gi](../gap/natural_tertiary.gi) | Defining systems, Python worker, phase projection, and exact boundary |
| [worker.py](../python/worker.py) | Exact JSON line protocol, one process per GAP session, face-equation audits, and the shared `V1` table |
| [low_phases.py](../python/low_phases.py) | Degree-zero, -one, and -two phases |
| [high_phase.py](../python/high_phase.py) | Degree-three phase |
| [cochain_tools.py](../python/cochain_tools.py) | The single interval-cut engine: cup-i words, integral signs, coboundary, interval pullback and prism, Q |
| [phase_eval.py](../python/phase_eval.py) | Shared lower phase, source splitting, chi, polarization and hD, integrality checks |
| [mod3_power.py](../python/mod3_power.py) | Prime-three formula and signed transport |
| [pages.gi](../gap/pages.gi) | Exact homology, surviving representatives, and page quotients |
| [group_api.gi](../gap/group_api.gi) | `koAHSS` and `koAHSS_batch` for finite groups and supplied resolutions |
| [extensions.gi](../gap/extensions.gi) | Marked abelian extension presentations, Smith transformations, `koFull` and `koFull_batch` |
| [extension_lifts.gi](../gap/extension_lifts.gi) | Immutable full generator lifts, affine B/C defining choices, integral D solves, and exact flatness witnesses |
| [extension_transfer.gi](../gap/extension_transfer.gi) | Native extension model, exact retraction preflight and sparse normalized transport |
| [extension_transfer.py](../python/extension_transfer.py) | Lazy transferred curvature, products and gauge actions using the fixed formulas |
| [extension_native_six.py](../python/extension_native_six.py) | Degree-six D-layer terms of the native model: `J_6` and the legal `gamma_6` with the pure-C rephasing, on the legal lower locus |
| [extension_acceleration.py](../python/extension_acceleration.py) | Exact evaluation policy of the extension worker: structural zeros, identity-memoized builders, persistent universal values, unstored coboundaries, scalar multiples and chi, bounded memo tables and chain caches for universal evaluations, and the recycled source registry |
| [universal_values.py](../python/universal_values.py) | Store of universal values shared by the extension and page workers: bundled values, the cache file and its source hash |
| [generate_universal_values.py](../python/generate_universal_values.py) | Recomputes the bundled universal values in [universal-values.json](../data/universal-values.json) for the current formula sources |
| [extension_bar.gi](../gap/extension_bar.gi) | Bounded complete finite-bar model, used only as a reference by the tests |
| [extension_relations.gi](../gap/extension_relations.gi) | Powers of fixed full lifts, common lower coordinates, and retained reduction carries |
| [extension_equivalence.gi](../gap/extension_equivalence.gi) | Exact ordered native gauge comparison, without assumed cochain associativity |
| [extension_degree_six.gi](../gap/extension_degree_six.gi) | Basis data for the degree-six section of the complete-bar reference model |
