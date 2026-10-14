# Documentation

These notes give, explicitly and degree by degree, the two things a
rewrite of the stacking computation has to reproduce:

1. the **generalized differential** $\mathfrak d_k$ on four-cochain states of
   package degree $k=-1,\dots,6$ (with the endpoint $k=7$), and
2. the **stacking rules**: the triangular product $\times_k$ with its
   corrections $\alpha_k,\beta_k,\gamma_k$, the identity
   $\mathfrak d_k(x\times y)=\mathfrak d_kx\times\mathfrak d_ky$, the gauge
   relation, and the two-layer model at the prime three.

The formulas were taken from the notes of the previous implementation
(`../fermionAHSS/doc`) and, where those notes left a formula to code, from
its Python modules (`../fermionAHSS/python`, `.../python/stacking_model`);
each such formula names its source function. The secondary differentials
`Tau`, `Psi` and the tertiary `T` are written out on their own pages, with
their helpers, universal source primitives and calibration constants; the
stacking pages use them through the interface table of
[conventions.md](conventions.md), Section 8.

## Reading order

| File | Content |
| --- | --- |
| [conventions.md](conventions.md) | notation (lifts, $\mathsf h$, $\delta_s$, cup-$i$ words, $\zeta$ words, $Q^j,E,Q_D$, polarizations, prism and ramp), and the interface table of the secondary and tertiary data |
| [secondary_operations.md](secondary_operations.md) | `Tau` and `Psi`: the defining cochain, $F$, $q$, $\psi'$, $\tau'$, the matched lift, the A-only splitting $(g,k_0,\lambda)$, the helpers $\zeta,\chi$ with the data file, domains, page-class reductions |
| [tertiary_operations.md](tertiary_operations.md) | `T`: the phase $\mathcal O_n$, $\Theta_m$, $R_0$ (even/odd), $\widehat R_1,\widehat R_2,\widehat R_3$, the nineteen-term $P^1_s$, the source primitives $V_1,V_2,V_3$ with their models and tables, the selectors and the independence of the contraction |
| [differential.md](differential.md) | $\mathfrak d_k$ written out for $k=-1,0,1,2,3,4,5,6,7$: lower components, the phase $\Omega_k$ (explicit for $k\le3$), $J_k$ on the legal locus, the prism $G_k$ or the section branch outside; the pure-C rule; square zero; the E6 layers |
| [product_lower.md](product_lower.md) | $\alpha_k$ and $\beta_k$ in every degree, including the polynomial $\mathcal V_0$, the extension $\beta^\sharp$, the comparison gauge $\xi_k$, the arbitrary-cochain branch and the natural degree-six product |
| [product_upper_k1_k3.md](product_upper_k1_k3.md) | $\gamma_1,\gamma_2,\gamma_3$: the degree-one closed form, the degree-two prism phase, the paper phase with its repair and coordinate gauge, the even and odd degree-three phases, the calibration, the off-legal branch and $q_3$; the upper identity, gauges and the group of a degree |
| [product_upper_k4_k6.md](product_upper_k4_k6.md) | $\gamma_4,\gamma_5,\gamma_6$: the common ingredients ($f^\sharp$, $J_0$, $j_0$, $\Theta$-pair phase, $\varphi$, $R_A$, $Q_k$), the pair phases $e^{\rm prod}_k$, normalization and rephasing, the off-legal branch with $q_4,q_5$, the direct $A=0$ formulas |
| [universal_primitives.md](universal_primitives.md) | the construction behind every universal primitive ($\mathcal V_n$, $\lambda_A$, $v_n$, $\mathfrak q$'s, $\mathcal U^\Theta_m$, $\mathcal U^{A0}_k$, $\mathcal P_k$, $V^{\rm PS}$): models, sections, contractions, sources and fixed coefficients |
| [degree_six_section.md](degree_six_section.md) | $g_6$, $\gamma_6$, $\gamma_7$ outside the legal locus: section, retraction, two legal cylinders, C-shift transport |
| [prime_three.md](prime_three.md) | the two-layer model at the prime three in degrees five and six |
| [extension_assembly.md](extension_assembly.md) | how the group of a degree is assembled: relations, target layer, the adjacent rows that need no stacking, localization at the primes, transport to a resolution |
| [power_relations.md](power_relations.md) | the relation of an A generator worked out: powers by squaring telescope, only diagonal pairs meet a pair primitive, and the two-fiber primitive on the diagonal is replaced by a one-fiber primitive plus 3, 4 or 22 constants |
| [implementation.md](implementation.md) | where each formula is evaluated in the package, the shortcuts, the worker and the caches |
| [natural_evaluation.md](natural_evaluation.md) | which formulas can be evaluated on an arbitrary resolution with its own cup-$i$ structure, rewritten in that form; what a homotopy package on a resolution must provide for the secondary and tertiary operations; why the universal primitives need a simplicial section; the two transport-free routes |

## Verification status

- The explicit $\tau'_0$ of (C12)–(C13), the identity $f^\sharp=\tau'$ on
  legal $B$, and the agreement of the two sources $g$, $g^\flat$ were
  checked numerically against the old modules on random cochains.
- All other formulas are transcriptions of the notes or of the code, with
  the function named; the old implementation checks the integrality of
  every $\gamma_k$, $g_k$, $Z_k$, $L$ at evaluation time, and these notes
  record those identities as the properties a reimplementation must verify.

## What is not in these notes

The proofs behind the fixed data are not bundled in the old repository
either: the degree-seven suspension companion that fixes `chi7_tail`,
$\epsilon$ and $\eta$; the source notes for $R_0$–$R_3$ and the Euler and
rank-six calibrations; and the stacking research workspace with the
integrality and exchange proofs. The data files the formulas read,
`data/chi-calibrated-degree7-anf.g` ($\chi_n$), `python/r3_source.json`
($V_3$), the one-fiber coefficient tables and the universal-value bundle,
are shipped with the package; see [implementation.md](implementation.md).
