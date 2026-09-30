# Degree-indexed four-cochain stacking formulas

This note fixes one set of cochain representatives for the stacking
product of the four-cochain model in package degrees \(k=-1,\ldots,6\),
together with the successor terms at \(k=7\) that the degree-six identities
use. The representatives are chosen together: every correction of degree
\(k\) uses the same literal successor that degree \(k+1\) uses, so phase
gauges that are valid only separately are never mixed. Coefficients are
exact, negative cochain degrees are zero, and every integral carry is
kept.

The differential is the all-cochain differential \(\mathfrak d_k\) of
[all_cochain_differential.md](all_cochain_differential.md), with the phases
\(\Omega_k\) of
[dimension_indexed_differentials.md](dimension_indexed_differentials.md).
The products of package degrees one and two are written out in
[low_degree_stacking.md](low_degree_stacking.md). The degree-six section
construction is derived in [G6_GENERAL_SECTION.md](G6_GENERAL_SECTION.md),
and [ALL_COCHAIN_OBSTRUCTION.md](ALL_COCHAIN_OBSTRUCTION.md) shows why
\(g_6\) must be nonzero outside \(L_6\). The conventions for cup-\(i\)
products \(\smile_i\), the signed coboundary \(\delta_s\), reduction
\(\rho\), lifts \(\widetilde x\), the half lift
\(\mathsf h(x)=\widetilde x/2\), the prism \(I\) and the ramp \(\ell\) are
those of [conventions.md](conventions.md) and Section 2 of the all-cochain
note. The formulas are implemented in
[../python/stacking_model](../python/stacking_model/), with their
[provenance](../python/stacking_model/provenance.json). Sections 8 and 9
state which parts `koFull` evaluates and how the formulas are called.

## 1. States, differential and product

A state of package degree \(k\) is \(x=(u,D)\) with \(u=(A,B,C)\) and

\[
(A,B,C,D)\in C^{k-3}(X;\mathbf Z_s)\times C^{k-2}(X;\mathbf F_2)
\times C^{k-1}(X;\mathbf F_2)\times C^{k+1}(X;\mathbf Z_s),
\]

one cochain per layer A \((k-3,0)\), B \((k-2,-1)\), C \((k-1,-2)\) and
D \((k+1,-4)\). The differential (A13) and the product are

\[
\mathfrak d_k(x)=\bigl(F_ku,\ \delta_sD+g_k(u)\bigr),\qquad
F_k(u)=\bigl(\delta_sA,\ \delta B+P_k(A),\ \delta C+f_k(A,B)\bigr),
\]

\[
x\times_k y=\bigl(\mu_k(u,v),\ D+D'+\gamma_k(u,v)\bigr),\qquad
\mu_k(u,v)=\bigl(A+A',\ B+B'+\alpha_k(A,A'),\ C+C'+\beta_k(A,B;A',B')\bigr),
\]

where \(P_k(A)=Q_D(\rho A)\), and where \(f_k\) and \(g_k\) are the entries
of the table in Section 3 of the all-cochain note. As in (A4),
\(p_k(A,B)=(\delta_sA,\ \delta B+P_k(A))\) and
\(L_k=\{(A,B):p_k(A,B)=0\}\), an equality of whole cochains on \(X\). On
\(L_k\), \(\tau_k=\tau'_{k-3}(A;B)\) and \(t_k=\delta C+\tau_k\). Binary
additions are taken modulo two. Integral lifts are made separately where
they are indicated, and an integral carry is never discarded. For a
function \(Z\) of one lower triple, write

\[
\Delta Z(u,v)=Z(u)+Z(v)-Z\bigl(\mu(u,v)\bigr),
\]

with the lower product \(\mu\) stated for each branch. The integral cocycle
\(\mathcal K_k=\delta_s\Omega_k\) of (A7) equals \(T_{k-3}\) on a full
defining system. The roman \(K_k\) of Section 3 is a different cochain.

## 2. The lower corrections

### 2.1 The B-layer correction

For binary cochains \(x,y\) of degree \(r\), the polarizations are

\[
h^j_r(x,y)=x\smile_{r-j+1}y+\delta x\smile_{r-j+2}y,\qquad
\delta h^j_r(x,y)+h^j_{r+1}(\delta x,\delta y)=Q^j(x+y)+Q^j(x)+Q^j(y),
\]

(`phase_eval.polarization`), and the stacking cross term is
(`phase_eval.hD`)

\[
h^D_r(a,a')=h^2_r(a,a')+s\smile h^1_r(a,a')
=a\smile_{r-1}a'+\delta a\smile_ra'
+s\smile\bigl(a\smile_ra'+\delta a\smile_{r+1}a'\bigr).
\]

Then

\[
\alpha_k(A,A')=h^D_{k-3}(\rho A,\rho A'),
\]

with \(\alpha_k=0\) when \(A\) is absent (\(k<3\))
(`off_shell_beta.alpha`, `stacking_lower.alpha`). This is an explicit
all-cochain formula.

### 2.2 The legal C-layer correction

Let \(g_A\) be the binary component `g` of the source splitting of the
secondary representative \(\tau'\) at \(A\) (`phase_eval.source_splitting`).
At \(k=3\) it is \(g_A=t\smile\omega+a\smile s\smile s\), with \(a=\rho A\)
and \(t=\rho\bigl((A-\widetilde a)/2\bigr)\) (`n0_beta.source0`). The
splitting carry is

\[
\nu_k(A,B)=\rho\left(\frac{\widetilde{g_A}+\widetilde{s\smile B}
-\widetilde{g_A+s\smile B}}2\right).
\]

Let \(\mathcal V_n(A,A')\) be the fixed universal pair source primitive of
input degree \(n\): a binary cochain of degree \(n+2\) whose coboundary is
the A-pair source. It is evaluated by a fixed finite contraction of
universal chains and is never solved on \(X\). It is distinct from the
source primitives \(V_1,V_2,V_3\) of the tertiary phases in
[universal_helpers.md](universal_helpers.md). For both inputs in \(L_k\),
put \(A_S=A+A'\) and \(B_S=B+B'+\alpha_k(A,A')\). The calibrated legal
correction is

\[
\beta^L_k=h^D(B,B')+h^D\bigl(B+B',\alpha_k(A,A')\bigr)
+\nu_k(A,B)+\nu_k(A',B')+\nu_k(A_S,B_S)+\mathcal V_{k-3}(A,A')
+[k=5]\,\rho A\smile\rho A'.
\]

It satisfies \(\delta\beta^L_k=\tau_k(A,B)+\tau_k(A',B')+\tau_k(A_S,B_S)\)
and does not depend on \(C\). The closed term \(\rho A\smile\rho A'\) at
\(k=5\) is fixed by the periods of the upper source. It changes no
secondary boundary identity and vanishes when \(A=0\) or \(A'=0\).
Because it is part of \(\beta^L_5\), it also enters the degree-four raw
correction below, where \(\beta^L_5\) is the successor.

| k | n | Pair primitive \(\mathcal V_n\) | Implementation |
| ---: | ---: | --- | --- |
| 3 | 0 | the polynomial \(\mathcal V_0\) below | `n0_beta.beta_legal`, `n0_beta.source_primitive` |
| 4 | 1 | face-sum contraction | `v1_pair_shared.legal_beta`, `v1_pair_shared.primitive` |
| 5 | 2 | face-sum contraction, plus \(\rho A\smile\rho A'\) | `v2_pair_shared.legal_beta` |
| 6 | 3 | face-sum contraction | `v3_pair_shared.legal_beta` |
| 7 | 4 | face-sum contraction; successor of degree six and \(k=7\) endpoint | `v4_pair_shared.legal_beta` |

The dispatcher is `stacking_lower.legal_beta`. On a face \((012)\), write
\(A(0)\equiv a_0+2a_1\) and \(A'(0)\equiv b_0+2b_1\pmod 4\), \(x=s(01)\),
\(y=s(12)\) and \(w=\omega(012)\). Then

\[
\mathcal V_0(A,A')(012)=a_0b_1\,xy+w\bigl(a_1b_1+a_0b_0(x+y)+(a_1b_0+a_0b_1)\,xy\bigr)\pmod 2,
\]

so the degree-three correction depends on \(A\) only modulo four.

### 2.3 The C-layer correction on arbitrary cochains

For \(3\le k\le 6\), let \(\xi_k(A,B)\) be the fixed comparison cochain with
\(\delta\xi_k=\tau_k(A,B)+H_k(A,B)\) on legal pairs, where \(H_k\) is the
prism formula (A9) (`h_tau_primitive.comparison_gauge`; \(\xi_3=0\)). Its
local extension \(q_{\rm loc}(A,B)\) takes the value of \(\xi_k\) on a
simplex whose faces all have vanishing curvature \(p_k(A,B)\), and zero on
other simplices (`off_shell_beta.local_comparison_extension`). Put

\[
\widehat f_k=H_k+\delta q_{\rm loc},\qquad
e(A,B)=\begin{cases}0,&(A,B)\in L_k,\\ q_{\rm loc}(A,B),&(A,B)\notin L_k,\end{cases}
\]

(`off_shell_beta.natural_f`, `off_shell_beta.coordinate_gauge`). The map
\(\widehat f_k\) agrees with \(\tau_k\) on \(L_k\), and \(\phi(u)=(A,B,C+e(A,B))\)
converts the coordinates of (A13) to the natural coordinates in which
\(f_k\) is replaced by \(\widehat f_k\). Accordingly, for \(3\le k\le6\) the
natural lower differential is

\[
\widehat F_k(u)=\bigl(\delta_sA,\ \delta B+P_k(A),\ \delta C+\widehat f_k(A,B)\bigr).
\]

The binary raw correction is the explicit prism (`off_shell_beta.raw_beta`)

\[
\beta^{\rm raw}_k(u,v)=I\Bigl[\widehat f_k(A_S,B_S)+\widehat f_k(A,B)
+\widehat f_k(A',B')+\beta^L_{k+1}\bigl(p_k(A,B),p_k(A',B')\bigr)\Bigr],
\]

in which the bracket is evaluated on the ramped cochains
\(A\ell,B\ell,A'\ell,B'\ell\) over \(X\times I\), with
\(A_S=A\ell+A'\ell\) and \(B_S=B\ell+B'\ell+\alpha_k(A\ell,A'\ell)\). The
successor \(\beta^L_{k+1}\) is evaluated only on curvature pairs, which lie
in \(L_{k+1}\) because \(p_{k+1}p_k=0\). The complete correction for
\(k=3,4,5,6\) is (`off_shell_beta.beta`, `stacking_lower.all_cochain_beta`)

\[
\beta_k(u,v)=\begin{cases}
\beta^L_k(A,B;A',B'),&(A,B),(A',B')\in L_k,\\[2pt]
\beta^{\rm raw}_k(u,v)+e(A,B)+e(A',B')+e(A_S,B_S),&\text{otherwise},
\end{cases}
\]

with \(A_S=A+A'\) and \(B_S=B+B'+\alpha_k(A,A')\). The value \(e\) on the
total depends only on its A and B components, so the definition is not an
implicit equation for \(C\). No \(\alpha\) or \(\beta\) depends on \(C\) or
\(D\). In the lower degrees, \(\beta_0=\beta_1=0\) and

\[
\beta_2(B,B')=\delta B\smile B'+s\smile(B\smile B')
\]

(`a0_gamma.beta2`). Its first term is needed when \(B\) is not closed.
These formulas satisfy

\[
F_k\mu_k(u,v)=\mu_{k+1}(F_ku,F_kv)
\]

on all lower cochains.

## 3. The upper correction

### 3.1 Normalization and rephasing, \(k\ge4\)

For \(4\le k\le6\), let \(e^{\rm prod}_k(u,v)\) be the production pair phase
(`closed_ab_upper.phase`, which is `production_gamma4.phase4`,
`production_gamma5.phase` or `production_gamma6.phase`), and let
\(\Gamma^{\rm prod}_k=\Delta\Omega_k+\delta_se^{\rm prod}_k\) be the
production legal product (`production_gamma4.gamma4`,
`production_gamma5.gamma`, `production_gamma6.gamma`). The common
normalization is (`pure_c_normalization`)

\[
\chi_k=\mathsf h\bigl(h^2(C,C')\bigr)-e^{\rm prod}_k(0,0,C;0,0,C'),\qquad
e^{\rm norm}_k=e^{\rm prod}_k+\chi_k,\qquad
\Gamma^{\rm norm}_k=\Delta\Omega_k+\delta_se^{\rm norm}_k.
\]

The cochain \(\chi_k\) is closed modulo one on unrestricted \(C,C'\), so

\[
N_k=\chi_k-\delta_sI\chi_k(A,B,C\ell;A',B',C'\ell)
\]

is integral with \(\delta_sN_k=\delta_s\chi_k\)
(`pure_c_normalization.integer_phase_correction`). The normalization
therefore changes the product by the integral coboundary \(\delta_sN_k\):
\(\Gamma^{\rm norm}_k=\Gamma^{\rm prod}_k+\delta_sN_k\).

Define the integral closed cochain

\[
K_k(u)=\begin{cases}\delta_s\mathsf h(s\smile C),&k\ge4,\ A=0,\ B=0,\ \delta C=0\ \text{on }X,\\
0,&\text{otherwise},\end{cases}
\]

(`all_cochain_upper.K`), so that \(K_k=0\) for \(k\le3\). The legal product
is

\[
\Gamma^L_k=\Gamma^{\rm norm}_k-\Delta K_k.
\]

For inputs in \(L_k\) with arbitrary \(C\), put
\(\Phi_k=\Omega_k+\mathsf h(t_k\smile_{k-1}\tau_k)\), as in
[README.md](README.md) (4). The closed-AB product is
(`closed_ab_upper.gamma`)

\[
\Gamma^{\rm closed}_k(u,v)=\Phi_k(u)+\Phi_k(v)-\Phi_k\bigl(\mu_k(u,v)\bigr)
+\delta_se^{\rm prod}_k(u,v)+\mathsf h\bigl(h^2(t_k(u),t_k(v))\bigr),
\qquad
\Gamma^{\rm closed,norm}_k=\Gamma^{\rm closed}_k+\delta_s\chi_k
\]

(`pure_c_normalization.closed_gamma`). The plus sign of the last term of
\(\Gamma^{\rm closed}_k\) fixes the pure-C successor. On full defining
systems \(t_k=0\), so \(\Gamma^{\rm closed,norm}_k=\Gamma^{\rm norm}_k\).
At \(k=3\), \(\Gamma^{\rm closed}_3\) is the same expression with
\(e^{\rm prod}_3\) replaced by the componentwise even/odd phase
\(e^{\rm cl}_3\) of `closed_ab_degree3.phase`, and no normalization is
applied: \(\Gamma^{\rm closed,norm}_3=\Gamma^{\rm closed}_3\). On
components where \(A\) is even, \(e^{\rm cl}_3\) contains the term
\(\mathsf h\bigl((H+z)\smile_3t'\bigr)\) of the canonical \(\tau\) of the
first input (`even_a_gamma3.canonical_data`). On components where \(A\)
is odd, it transports the complete curvature source along an orientation
prism and a flattening prism. These curvature terms vanish on full legal
data.

### 3.2 The correction used in each degree

| k | \(\gamma_k\) | Reference API (`four_cochain_stacking`) | Native `koFull` engine |
| ---: | --- | --- | --- |
| −1, 0 | \(0\): ordinary signed integral addition | `Stacking.xtimes` | D layer only; no stacking model is built |
| 1 | \(\Delta\mathsf h(E(C))+\delta_se_1+e_2(F_1u,F_1v)\), \(e_1(v_0v_1)=\tfrac12C(v_0)\bigl[C'(v_0)+C'(v_1)\bmod2\bigr]\) | `coherent_low_commutative.DegreeOneCommutativeStacking.gamma` | the same |
| 2 | \(\Delta\Phi_2+\delta_se_2+e_3(F_2u,F_2v)\); \(e_2\) is the successor prism with its C-transport term | `coherent_low_commutative.DegreeTwoCommutativeStacking.gamma` | the same |
| 3 | the branches of Section 3.3 with \(q_3\), the phase \(e^{\rm cl}_3\) and the calibration of Section 3.4 | `all_cochain_upper.gamma` (`closed_ab_degree3`, `degree3_legal_recalibration`) | the same |
| 4, 5 | the branches of Section 3.3 | `all_cochain_upper.gamma` | the same outside the both-legal branch; on it, without \(\delta_s\chi_k\), and, at \(k=5\) with \(A=B=0\), through `pure_c_gamma` (Section 8) |
| 6 | the section expression of Section 4 | `production_g6_section.ProductionG6Section.gamma` | lower-legal pairs only, without \(\delta_s\chi_6\) and through `pure_c_gamma` when \(A=B=0\); refused with two nonzero A layers unless `FERMIONAHSS_DEGREE_SIX_A_STACKING=1`; otherwise unresolved (Section 8) |

Here \(E(C)=\omega\smile C\) in degree one, and \(\Phi_2\) is the potential
of [low_degree_stacking.md](low_degree_stacking.md). The degree-three
successor \(e_3\) is the commutative paper-calibrated production phase at
\(A=A'=0\), `paper_commutative_production.phase`. The module
`coherent_low_commutative` regenerates \(e_2\) and \(e_1\) from this
literal successor, including all integral carries.

### 3.3 Degrees three to five

For \(k=3,4,5\), let \(\widehat F_k\) be the natural lower differential of
Section 2.3, and let \(\widehat\mu_k\) be the product with
\(\beta^{\rm raw}_k\) in place of \(\beta_k\) (`natural_upper.product`).
Define (`natural_upper.ghat`, `natural_upper.gamma`, and for
\(\widehat\Phi_k\) `nonzero_degree3_comparison.raw_potential` at \(k=3\)
and `nonzero_offshell_comparison.raw_potential` at \(k=4,5\))

\[
\widehat g_k(u)=-I\,\mathcal K_{k+1}\bigl(\widehat F_k(u\ell)\bigr),\qquad
\widehat\Phi_k(u)=I\,\Omega_{k+1}\bigl(\widehat F_k(u\ell)\bigr),
\]

\[
\widehat\gamma_k(u,v)=I\Bigl[\Delta_{\widehat\mu}\widehat g_k(u,v)
+\Gamma^{\rm norm}_{k+1}\bigl(\widehat F_ku,\widehat F_kv\bigr)\Bigr]\bigl(u\ell,v\ell\bigr).
\]

The argument of \(\mathcal K_{k+1}\) is a full defining system, since
\(\delta\widehat f_k=\tau_{k+1}(p_k)\). The explicit integral comparison
\(Z_k\) (`all_cochain_upper.integer_gauge`) satisfies

\[
g_k(u)=\widehat g_k(\phi u)+\delta_sZ_k(u)+K_{k+1}(F_ku).
\]

The last term vanishes off \(L_k\), where the first two output components
are nonzero. On \(L_k\),

\[
\begin{aligned}
Z_3&=\Phi_3-\widehat\Phi_3-\mathsf h(s\smile t_3)-\delta_sq_3,\\
Z_4&=\Phi_4-\widehat\Phi_4-\mathsf h(s\smile t_4)-\delta_sq_4,\\
Z_5&=\Phi_5-\widehat\Phi_5-\mathsf h(s\smile t_5)+\tfrac83A^{\smile3}-\delta_sq_5,
\end{aligned}
\qquad t_k=\delta C+\tau_k(A,B),
\]

where \(A^{\smile3}\) is the signed ordered cube \((A\smile_sA)\smile_sA\)
of [conventions.md](conventions.md) (C3). The gauges are
`nonzero_degree3_comparison.gauge` (\(q_3\)), in which the rank-zero, even
and odd cases are selected at vertices, and
`nonzero_offshell_comparison.legal_pair_gauge4` and `legal_pair_gauge5`
(\(q_4,q_5\)). The source coefficients of \(q_4,q_5\) and both degree-five
B transgressions are fixed there. The cube is closed, and its presence
changes no value of \(g_5\). Off \(L_k\), \(Z_3=0\), and in degrees four and
five \(Z_k\) is the explicit integral C-shift prism
`nonzero_offshell_comparison.off_branch_integer_gauge`.

For \(u,v\) both in \(L_k\), the carry
\(\Delta\mathsf h(s\smile t_k)\) is an integral cochain, though not a
coboundary in general, and

\[
\gamma_k=\begin{cases}
\Gamma^{\rm closed,norm}_k-\Delta\mathsf h(s\smile t_k)-\Delta K_k,
&(A,B),(A',B')\in L_k,\\[2pt]
\widehat\gamma_k(\phi u,\phi v)+\Delta Z_k-\Delta K_k,&\text{otherwise},
\end{cases}
\]

where \(\Delta\) uses the product \(\mu_k\) (`all_cochain_upper.gamma`).
The identity of Section 5 holds exactly in both branches; the proof is in
the stacking research workspace (not bundled). In degree three the same
proof applies with \(K_3=0\), \(q_{\rm loc}=0\) and the explicit \(q_3\)
(`nonzero_degree3_comparison.gauge`), and the calibration of Section 3.4
is added. In degrees four and five the full legal branch is exactly
\(\Gamma^L_k\).

### 3.4 Degree-three legal phase calibration

This closed legal phase correction fixes the exchange calibration of the
degree-three product. Let \(e^{\rm cl}_3\) be `closed_ab_degree3.phase`,
and let \(e^{\rm pap}_3\) be the explicit production-coordinate transport
of the paper exchange phase, `paper_commutative_production.phase`. Put
\(m(A)=(A^2-1)/8\bmod 2\). The cochain \(\chi^{\rm cal}_3\) is defined
componentwise, by the values of \(A\) and \(A'\) at the first vertex
(`degree3_legal_recalibration.phase_change`):

\[
\chi^{\rm cal}_3=\begin{cases}
e^{\rm pap}_3-e^{\rm cl}_3,&A=A'=0,\\
\tfrac12\,m(A)\,(B+B')^{\smile3},&A,A'\ \text{both odd},\\
0,&\text{otherwise},
\end{cases}
\]

(`gamma3_odd_symmetry.phase_correction` for the odd case). The correction
is selected only when both input triples are full defining systems,
\(F_3u=F_3v=0\) on \(X\). Then \(\delta_s\chi^{\rm cal}_3\) is integral and
closed, and on this branch the degree-three product is

\[
\gamma_3=\Gamma^{\rm closed}_3-\Delta\mathsf h(s\smile t_3)+\delta_s\chi^{\rm cal}_3
\]

(`degree3_legal_recalibration.gamma_change`). The other branches of
Section 3.3 do not contain this term. Adding a closed integral cochain
changes neither the strict degree-three upper identity nor any \(g_k\).
Every degree-two differential image is a full defining system with \(A=0\),
so the degree-two and degree-one prisms of Section 3.2 use exactly this
successor. The branch must be global. The legal phase difference has a
nonzero period on \(B\mathbf Z/2\times K(\mathbf F_2,2)\), so this phase
change cannot be extended naturally to arbitrary \(C\) with the same
successor. The period computation (a six-term cycle) is in the stacking
research workspace (not bundled). The odd term is closed because \(B+B'\)
is a binary cocycle on components where both ranks are odd, and \(m(A)\) is
invariant under the sign of \(A\). Its exchange primitive is
`gamma3_odd_symmetry.exchange`; the proof that the odd term is necessary is
in the stacking research workspace (not bundled).

## 4. Degree six: \(g_6\) and \(\gamma_6\)

The derivation is [G6_GENERAL_SECTION.md](G6_GENERAL_SECTION.md). Its
engine is `g6_general_section` (`FiniteCurrentSection`,
`NaturalSectionStacking`, `CurrentSectionStacking`), and its operations are
bound to the fixed formulas by `production_g6_section.ProductionG6Section`.

**Section and natural product.** Choose a pointed section \(\sigma\) of
\(\widehat F_6\) onto its image: \(\sigma(0)=0\), and
\(\widehat F_6\sigma=\mathrm{id}\) on \(\operatorname{im}\widehat F_6\). The
finite implementation computes it by exact lower solves and a finite
enumeration of integral cycles modulo four, as described in Section 3 of
[all_cochain_differential.md](all_cochain_differential.md). Let \(\ast\) be
the natural triangular product with correction \(\alpha_6\) and

\[
\widehat\beta_6=\beta^{\rm raw}_6+\delta\lambda_{\rm loc}
\]

(`lower_raw_comparison.natural_beta`). Here \(\lambda_{\rm loc}\) is the
local extension (`lower_raw_comparison.local_extension`) of the explicit
comparison \(\lambda_3\) with \(\delta\lambda_3=\beta^L_6+\beta^{\rm raw}_6\)
on legal pairs (`lower_raw_comparison.primitive`). The product \(\ast\)
agrees with the calibrated product on legal inputs and makes
\(\widehat F_6\) a homomorphism. Let \(a\backslash u\) be the explicit
triangular solution of \(a\ast z=u\), and put

\[
\mathcal R(u)=\sigma(\widehat F_6u)\backslash u,\qquad \mathcal R_6(u)=\mathcal R(\phi u).
\]

Left cancellation gives \(\widehat F_6\mathcal R(u)=0\), so \(\mathcal R(u)\) is always a full
defining system. In the coordinates of (A13), the upper differential is

\[
g_6(u)=\begin{cases}J_6(u),&(A,B)\in L_6,\\ T_3\bigl(\mathcal R_6(u)\bigr),&\text{otherwise},\end{cases}
\]

where \(T_3=\delta_s\Omega_6\) is the fixed integral representative (A7),
including its coefficient-two prime-three term. Every value on \(L_6\) is
\(J_6\). Both adjacent square-zero identities hold exactly for arbitrary
\(D\). With \(g_6=0\) outside \(L_6\), no normalized unital triangular
stacking product satisfies
\(\mathfrak d_6(x\times_6y)=\mathfrak d_6(x)\times_7\mathfrak d_6(y)\) on all
cochains. On \(B(\mathbf Z/3)^2\) the full defining system
\((4\beta_{\mathbf Z,3}(u_1u_2),0,0,0)\) has a nonzero \(T_3\) class, and
stacking it with an input outside \(L_6\) would force that class to be a
coboundary ([ALL_COCHAIN_OBSTRUCTION.md](ALL_COCHAIN_OBSTRUCTION.md)). With
\(g_6=T_3\mathcal R_6\) the off-shell class is the same nonzero class, so
there is no contradiction.

**The pair correction.** In natural coordinates, write

\[
w=\widehat F_6u,\quad w'=\widehat F_6v,\quad a=\sigma(w),\quad b=\sigma(w'),\quad
z=\mathcal R(u),\quad z'=\mathcal R(v),
\]
\[
q=a\ast b,\quad q_0=\sigma(\widehat F_6q),\quad \kappa=q_0\backslash q,\quad
r=q\backslash(u\ast v).
\]

The explicitly parenthesized cylinders over \(X\times I\)

\[
r_\ell=(a\ell\ast b\ell)\backslash\bigl[(a\ell\ast\pi^*z)\ast(b\ell\ast\pi^*z')\bigr],\qquad
k_\ell=(q_0\ell)\backslash\bigl[\bigl((q_0\ell)\ast\pi^*\kappa\bigr)\ast\pi^*r\bigr]
\]

are legal, with endpoints \(i_0^*r_\ell=z\ast z'\), \(i_1^*r_\ell=r\),
\(i_0^*k_\ell=\kappa\ast r\) and \(i_1^*k_\ell=\mathcal R(u\ast v)\)
(`g6_general_section.legal_cylinder_transports`; \(\pi^*\) is the plain
pullback from \(X\)).
Let

\[
\Gamma^{\rm sec}_6=\Gamma^{\rm closed,norm}_6-\Delta\mathsf h(s\smile t_6)-\Delta K_6
\]

be the closed-AB product of the section construction
(`ProductionG6Section.closed_gamma`), and let \(\Gamma^L_6\) be
`ProductionG6Section.legal_gamma`; they are \(\gamma_6^{\rm cl}\) and
\(\gamma_6^{\rm leg}\) of [G6_GENERAL_SECTION.md](G6_GENERAL_SECTION.md),
Section 7. Define

\[
j(w)=\begin{cases}J_6(\sigma w),&\text{the first two components of }w\text{ vanish},\\0,&\text{otherwise},\end{cases}
\qquad
\mathcal Q(u)=\begin{cases}\Gamma^{\rm sec}_6\bigl(\sigma(\widehat F_6u),\mathcal R(u)\bigr),&(A,B)\in L_6,\\0,&\text{otherwise},\end{cases}
\]

\[
\mathcal G=\Gamma^L_6(z,z')+\Gamma^L_6(\kappa,r)-I\,T_3(r_\ell)-I\,T_3(k_\ell),
\qquad
\delta_s\mathcal G=T_3(z)+T_3(z')-T_3\bigl(\mathcal R(u\ast v)\bigr)+T_3(\kappa).
\]

When \(u,v\) are not both in \(L_6\), the natural correction and its
successor on differential images are

\[
\gamma_6=\mathcal G-\mathcal Q(u)-\mathcal Q(v)+\mathcal Q(u\ast v),\qquad
\gamma_7(w,w')=T_3(\kappa)+j(\widehat F_6q)-j(w)-j(w').
\]

When both \(u,v\) lie in \(L_6\), put
\(\gamma_6=\Gamma^{\rm sec}_6(u,v)\). On images \(w=(0,0,c)\) and
\(w'=(0,0,c')\), the successor is then the rephased pure-C product

\[
\gamma_7(w,w')=\gamma_C(c,c')-K_7(c)-K_7(c')+K_7(c+c'),
\]
\[
\gamma_C(C,C')=\mathsf h(E(C))+\mathsf h(E(C'))-\mathsf h(E(C+C'))
+\delta_s\mathsf h\bigl(h^2(C,C')\bigr)+\mathsf h\bigl(h^2(\delta C,\delta C')\bigr),
\]

with \(K_7(c)=\delta_s\mathsf h(s\smile c)\)
(`ProductionG6Section.pure_successor`, `compatible_sector.pure_c_gamma`).
The half-lift carry in \(\mathcal Q(u)\) is zero because \(\mathcal R(u)\) has
\(t_6=0\), but it is kept for general closed-AB pairs. The branches are
selected by the output pair \((w,w')\) itself, so \(\gamma_7\) is well
defined. Outside the image of \(F_6\), its unused endpoint values may be set
to zero; no degree-seven classification is asserted. When both inputs are
section representatives (\(z=z'=0\)), then \(r=0\), the first cylinder
vanishes and the second is the constant \(\kappa\). The implementation then
evaluates \(\gamma_6=-\mathcal Q(u)-\mathcal Q(v)+\mathcal Q(u\ast v)\)
directly.

**Return to the product of (A13).** The natural product differs from the
current one by the C shift \(\delta\lambda\), where \(\lambda=0\) when both
inputs lie in \(L_6\), and \(\lambda=\lambda_{\rm loc}\) otherwise. Add to
\(\gamma_6\) the integral transport (`g6_general_section.C_shift_transport`)

\[
\mathcal L(y,\lambda)=\mathsf h\bigl(E(C_y+\delta\lambda)\bigr)-\mathsf h\bigl(E(C_y)\bigr)
-\delta_s\mathsf h\bigl(E(\lambda)+h^2(C_y,\delta\lambda)\bigr),
\]

where \(y\) is the current total when it lies in \(L_6\), and
\(y=\mathcal R_6(\text{total})\) otherwise. A C-coboundary changes neither
\(\widehat F_6\) nor its section, so the retract undergoes the same shift.
The plus sign is fixed by
\(\delta_s\gamma_6^{\rm nat}=g_6(u)+g_6(v)-g_6(\text{natural total})+\gamma_7\):
the transport adds \(g_6(\text{natural total})-g_6(\text{current total})\)
(`CurrentSectionStacking.gamma`, which also checks the lower
representative comparison).

No associativity of the lower cochain product is used. For a general
\(X\), a fixed pointed set section suffices, and the proof never chooses a
section on \(X\times I\). The finite free algorithm makes the choice
computable; its search is finite but exponential in the integral cycle rank
and in the dimension of the binary kernel. The resulting extension depends
on the global section and need not be natural under pullback.

## 5. The two identities

With these conventions, the upper condition is

\[
\delta_s\gamma_k(u,v)=g_k(u)+g_k(v)-g_k\bigl(\mu_k(u,v)\bigr)+\gamma_{k+1}(F_ku,F_kv).
\]

Together with the lower identity \(F_k\mu_k=\mu_{k+1}(F_k,F_k)\), it
gives

\[
\mathfrak d_k(x\times_ky)=\mathfrak d_k(x)\times_{k+1}\mathfrak d_k(y)
\]

for arbitrary integral \(D,D'\). The square-zero condition (A21) is

\[
F_{k+1}F_k=0,\qquad \delta_sg_k+g_{k+1}F_k=0.
\]

For \(u\) outside \(L_6\), \(T_3(\mathcal R_6u)\) is closed, and the
degree-seven endpoint \(g_7\) of (A15) vanishes on its output, whose first
two components are nonzero. On \(L_6\) the cancellation of \(J_6\) against
the pure-C endpoint (A18) applies. Degree-five outputs always lie in
\(L_6\), so the square-zero identity \(\delta_sg_5+g_6F_5=0\) uses only the
\(J_6\) branch of \(g_6\) and holds by (A18) for \(u\) in \(L_5\) and by
(A20) otherwise.

## 6. Prime exchange

Literal symmetry of the cochain products is not required. On legal
defining systems, the exchange gives a lower cylinder, an integral D
transport \(L\) and an explicit integral \(M\) with

\[
\gamma_k(v,u)-\gamma_k(u,v)+L=\delta_sM.
\]

Thus the two orders represent the same stacked class, and the change of
the lower representative is retained. The exchange data per degree:

- \(k\le0\): literal symmetry.
- \(k=1\): on legal inputs \(\gamma_1=\omega\smile C\smile C'\), which is literally symmetric.
- \(k=2\): \(\varsigma=\mathsf h(C\smile_1C'+C\smile B')\), \(L=0\) and
  \(M=e_2^{\rm op}-e_2-\delta_s\varsigma\). The exchange gauge
  \(\varsigma\) is `paper_stacking.phase_swap_gauge_k2`. The integrality of
  \(M\) for the coherent phase \(e_2\) is established in the stacking
  research workspace (not bundled).
- \(k=3\): `degree3_legal_recalibration.exchange` returns the tuple
  `(lam, K, sigma, L, M)`, whose entry `sigma` is the degree-three exchange
  gauge, the analogue of \(\varsigma\). It selects, by the rank at each
  vertex, the paper-calibrated zero-rank formula
  (`paper_commutative_production.exchange`), the orientation formula for
  nonzero even rank (`gamma3_even_symmetry.exchange`), or the two-prism
  formula for odd rank (`gamma3_odd_symmetry.exchange`). On legal data every
  selector is constant on connected components.
- \(k=4,5,6\): the fixed source contraction and a lower cylinder,
  transported through the pure-C normalization \(\chi_k\) and the integral
  rephasing \(K_k\). Both have their own explicit transport terms. This
  construction is in the stacking research workspace (not bundled), and no
  bundled module evaluates it.

Without its curvature terms, a raw all-cochain commutator need not be a
coboundary. Neither the reference model nor `koFull` evaluates these
exchange data: `koFull` assumes that stacking is commutative and associative
on gauge classes (`abelianQuotientAssumed`).

## 7. Relation to arXiv:2310.19058

The source of the stacking law is Ren, Ning, Qi, Wang and Gu,
[arXiv:2310.19058v2](https://arxiv.org/html/2310.19058v2), here written in
the repository's calibrated Wang–Gu, `Tau`, `Psi` and `T` cochain
conventions. The paper's lower layers and phase for \(A=0\) and
\(k=1,2,3\) are transcribed in `paper_stacking.paper_corrections`. In these
cup conventions, the transcribed degree-three phase fails the cochain
consistency equations for a generic \(s\). The local degree-three
correction `paper_stacking.antiunitary_phase_repair` (solved over
\(\mathbf Z/4\) from the universal consistency equations, and zero when
\(s\), \(B\) or \(B'\) vanishes) is added in
`paper_stacking.repaired_paper_corrections`. This correction is established
on legal data only and does not cover nonzero \(A\). The present
construction establishes the stated cochain identities. It does not by
itself establish completeness of the gauge relation, naturality in all
degrees, or a ko classification.

## 8. Use in the package

**Reference model.** `four_cochain_stacking.Stacking` evaluates exactly
the formulas of this note: degrees \(-1\) to \(5\) and the endpoint
\(\mathfrak d_7\) of (A15) directly, and degree six and the endpoint
product \(\times_7\) (on images of \(\mathfrak d_6\)) through
`ProductionG6Section`. It is the engine of the complete-bar reference
model: [../python/extension_worker.py](../python/extension_worker.py)
constructs `Stacking`, [../python/extension_degree_six.py](../python/extension_degree_six.py)
attaches the degree-six adapter, and GAP drives both from
`gap/extension_bar.gi` and `gap/extension_degree_six.gi`. This model is
used only by the tests (`tst/extension_degree_six.tst` is opt-in), never by
`koFull`.

**Native engine.** `koFull` evaluates the formulas on a supplied resolution
through the transfer of [transfer.md](transfer.md) and
[resolution-extensions.md](resolution-extensions.md), in
[../python/extension_transfer.py](../python/extension_transfer.py) and
[../python/extension_native_upper.py](../python/extension_native_upper.py):

- Degrees \(-1\) and \(0\) have only the D layer, and no stacking model is
  built.
- Degrees one and two use `DegreeOneCommutativeStacking`,
  `DegreeTwoCommutativeStacking` and \(\beta_2\), as in the reference
  model.
- In degrees three to six the lower layers use \(\alpha_k\) and the
  complete \(\beta_k\) of Section 2 (`stacking_lower.all_cochain_beta`).
- Degree three uses `all_cochain_upper.g` and `all_cochain_upper.gamma`,
  including the calibration of Section 3.4.
- In degrees four and five, \(g_k\) is `all_cochain_upper.g`, and the
  product outside the both-legal branch is `all_cochain_upper.gamma`.
- On the both-legal branch of degrees four to six, the product is
  \(\Gamma^{\rm closed}_k-\Delta\mathsf h(s\smile t_k)-\Delta K_k\)
  (`extension_native_upper.legal_gamma`), without the normalization
  \(\delta_s\chi_k\). It therefore differs from \(\gamma_k\) of this note
  by the integral coboundary \(\delta_sN_k\), which a D-gauge absorbs; the
  gauge classes are unchanged.
- In degrees five and six, for a lower-legal triple with \(A=0\) on the
  resolution, \(g_k\) is evaluated by the direct formula
  `extension_native_upper.a0_curvature`, equal to \(g_k=J_k\) there as a
  cochain. For two full defining systems with \(A=B=0\), the product uses
  `compatible_sector.pure_c_gamma` with the same carry and \(K\) change,
  which again differs by an integral coboundary
  ([extensions.md](extensions.md)).
- In degree six, the engine evaluates \(g_6=J_6\) (`closed_ab_upper.J`) and
  the legal product without normalization (with the \(A=B=0\) dispatch
  above) only on the legal lower locus. Native curvature stops at the first
  nonzero layer, and products and gauge actions receive flat states and
  differential images, whose lower pairs are legal. A request outside
  \(L_6\) raises `SectionBranchRequired` and leaves the degree unresolved.
  The section \(\sigma\), the retraction \(\mathcal R_6\), the legal
  cylinders and \(\mathcal L\) are never evaluated on the resolution.
- The degree-six correction of two states whose A layers are both nonzero
  nests the V3 contraction in its universal pair source. It is refused
  (`PairSourceLimit`, the degree unresolved) unless
  `FERMIONAHSS_DEGREE_SIX_A_STACKING=1` is set. When at most one A layer is
  nonzero, the pair primitive vanishes identically and the product is
  evaluated.

The complete model of this note measures the two-primary relations. After
localization, relations at primes at least five, and three-primary
relations below degree five, split; three-primary relations of degrees
five and six use the two-layer model
[../python/extension_three_local.py](../python/extension_three_local.py),
whose potential is the three-primary term
\(\tfrac23\widetilde{P^1_s\rho_3A}\) of \(\Omega_5\) and \(\Omega_6\);
`FERMIONAHSS_PRIME_LOCAL=0` restores the complete model for every relation
([extensions.md](extensions.md), "Localization at the primes").

## 9. Calling the implemented family

`four_cochain_stacking.State(A,B,C,D)` stores cochains of degrees
\((k-3,k-2,k-1,k+1)\), with \(k\) read from \(D\). Construct
`Stacking(s,omega,is_zero,degree_six=adapter)`, then use `rule.d(x)` and
`rule.xtimes(x,y)`. Cochains are exact callable face functions.
`is_zero` must test a whole cochain on \(X\), not a single simplex value,
because it decides the global branches. Without an adapter, the API covers
\(\mathfrak d_k\) and \(\times_k\) for \(k=-1,\ldots,5\) and
\(\mathfrak d_7\); \(\mathfrak d_6\), \(\times_6\) and \(\times_7\) need the
adapter `production_g6_section.ProductionG6Section`, which is supplied only
with the linear data of a finite normalized cochain basis: coordinates, the
binary differentials and the conversions between cochains and vectors. At
\(k=7\), \(\mathfrak d_7\) is the endpoint (A15), and \(\times_7\) is
defined only on images of \(\mathfrak d_6\). No nonlinear upper operation is
left for the caller to supply. The complete-bar reference model
(`extension_degree_six.configure_degree_six`) builds the adapter from a
normalized bar basis supplied by GAP.
