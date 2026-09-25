# Implemented formulas and calibration data

These documents describe the fixed bar-resolution edition of fermionAHSS.
They adapt the consolidated formula references and relevant backend/status
documentation from the koAHSS research workspace. They preserve the
formulas and calibration data while updating transport and Python paths.

| Document | Contents |
| --- | --- |
| [conventions.md](conventions.md) | Signed coefficients, interval-cut words, integral signs, binary lifts, primary maps, bar comparison, and defining-cochain corrections |
| [secondary_operations.md](secondary_operations.md) | Complete Tau and Psi formulas, exact common lift, domains, and indeterminacy |
| [tertiary_operations.md](tertiary_operations.md) | Implemented R0–R3 phases, current cubic correction, and all 19 prime-three terms |
| [universal_helpers.md](universal_helpers.md) | Chi ANF decoding and degree table, zeta words, universal contractors, finite source tables, and every normalization selector |
| [backends.md](backends.md) | HAP, cochain, and page interfaces; direct-operation audits and exact quotient conventions |
| [mathematical-status.md](mathematical-status.md) | Implemented range, mathematical assumptions, historical verification, and unresolved scope |
| [dimension_indexed_differentials.md](dimension_indexed_differentials.md) | Explicit nonlinear differential through k=6 on the first-two-layer domain, the residual correction, and a literal square-zero proof |
| [all_cochain_differential.md](all_cochain_differential.md) | Piecewise extension to arbitrary cochains with first component δ_sA, complete degree table, exact square-zero proof, and unchanged defining-system obstructions; not a natural local extension |
| [g6_repair.md](g6_repair.md) | Section repair retaining full T3, with exact square-zero and degree-six stacking identities, a finite implementation, and the historical projection construction |
| [Four-cochain stacking](../../fermionAHSS_stacking/ALL_COCHAIN_STACKING.md) | Separate research implementation through k=6: explicit corrections, strict all-cochain compatibility, and legal exchange coboundaries; global choices are stated |
| [Worked (Z/4)^4 obstruction](all_cochain_differential.md#the-z4-four-factor-obstruction) | Incorporated in both differential notes: A=2 admits B,C but no D, including all lower choices, the detecting cycle, and verification scope |

The large chi word lists are encoded exactly by
[data/chi-calibrated-degree7-anf.g](../data/chi-calibrated-degree7-anf.g)
and the decoding rule in the helper reference. This is the complete
runtime coefficient data; the original millions of word summands are
not needed separately.

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
| Current R2 | `R2sharp = R2old - A^cup3/4` |
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
| [worker.py](../python/worker.py) | Exact JSON protocol and face-equation audits |
| [low_phases.py](../python/low_phases.py) | Degree-zero, -one, and -two phases |
| [high_phase.py](../python/high_phase.py) | Degree-three phase |
| [phase_eval.py](../python/phase_eval.py) | Shared lower phase and source splitting |
| [mod3_power.py](../python/mod3_power.py) | Prime-three formula and signed transport |
| [pages.gi](../gap/pages.gi) | Exact homology, surviving representatives, and page quotients |
