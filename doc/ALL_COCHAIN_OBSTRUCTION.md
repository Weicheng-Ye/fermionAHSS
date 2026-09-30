# The degree-six obstruction outside the legal locus

This note explains why the degree-six D-layer correction \(g_6\) of the
all-cochain differential is not zero outside the legal locus \(L_6\). The
degree-six row of the table in
[all_cochain_differential.md](all_cochain_differential.md), Section 3,
sets

\[
g_6(A,B,C)=\begin{cases}
J_6(A,B,C),&(A,B)\in L_6,\\
T_3\bigl(\mathcal R_6(A,B,C)\bigr),&(A,B)\notin L_6,
\end{cases}
\tag{O1}
\]

where \(\mathcal R_6\) is the section retraction onto full legal triples.
If the second line were replaced by zero, the differential would still
square to zero and would have the same defining-system equations, but it
would admit no normalized unital triangular stacking product compatible
with it on all cochains. The proof uses an explicit input on
\(B((\mathbf Z/3)^2)\). The defining system
\((4\beta_{\mathbf Z,3}(u_1u_2),0,0,0)\) has a nonzero \(T_3\) class.
Stacking it with an input outside \(L_6\) would force that class to be a
coboundary. The stacking law that uses (O1) is
[ALL_COCHAIN_STACKING.md](ALL_COCHAIN_STACKING.md), and its degree-six
product is constructed in [G6_GENERAL_SECTION.md](G6_GENERAL_SECTION.md).
Notation is that of [all_cochain_differential.md](all_cochain_differential.md)
and of the formula sheet in [README.md](README.md).

## 1. Statement and scope

Write

\[
p_k(A,B)=(\delta_sA,\ \delta B+P_k(A)),\qquad
L_k=\{(A,B):p_k(A,B)=0\},
\]

as in (A4). Let \(\mathfrak d^0\) be the differential (A13) through package
degree \(k=6\) (spatial dimension \(d=5\)) in which the degree-six D-layer
correction is replaced by

\[
g^0_6(A,B,C)=\begin{cases}J_6(A,B,C),&(A,B)\in L_6,\\
0,&(A,B)\notin L_6.
\end{cases}
\tag{O2}
\]

It squares to zero by the argument of Section 4 of
[all_cochain_differential.md](all_cochain_differential.md). Zero is a
closed integral value, and the endpoint \(g_7\) of (A15) vanishes on the
outputs, whose first pair is nonzero. Its zero-output equations are the
defining-system equations, as in Section 5 there.

A stacking product of package degree \(k\) is triangular:

\[
x\times_k y=\bigl(\mu_k(u,u'),\ D+D'+\gamma_k(u,u')\bigr),\qquad
x=(u,D),\quad y=(u',D'),\quad u=(A,B,C),
\]

with lower product

\[
\mu_k(u,u')=\bigl(A+A',\ B+B'+\alpha_k(A,A'),\
C+C'+\beta_k(A,B;A',B')\bigr).
\]

**Proposition.** No normalized unital triangular stacking product satisfies

\[
\mathfrak d^0_k(x\times_k y)=\mathfrak d^0_k(x)\times_{k+1}\mathfrak d^0_k(y)
\]

on all cochains through \(k=6\). This holds even when neither symmetry nor
associativity is demanded.

The obstruction comes from the value zero in (O2) alone. It is not an
obstruction to stacking of flat states, to the stacking law of Ren, Ning,
Qi, Wang and Gu ([arXiv:2310.19058](https://arxiv.org/abs/2310.19058)) on
solutions, or to the classification. It shows that an all-cochain
extension admitting a normalized unital triangular product with the strict
identity must keep a nonzero \(g_6\) outside \(L_6\); (O1) does so. An
extension that keeps the value zero must give up either the strict
identity on all cochains or the normalized stacking unit.

The unit is the zero tuple. *Normalized unital* means that the triangular
correction cochains vanish when one entire input is zero. Since
\(\alpha_k\) does not depend on \(B,C,D\), this includes
\(\alpha_k(A,0)=0\) for all \(A\). This holds for the explicit
\(\alpha_k=h^D_{k-3}(\rho A,\rho A')\) of
[ALL_COCHAIN_STACKING.md](ALL_COCHAIN_STACKING.md), Section 2.1, where
\(h^D=h^2+s\smile h^1\) (`phase_eval.hD`, used by `compatible_sector.alpha`).
Since \(\gamma_k\) does
not depend on \(D,D'\), this includes \(\gamma_k(0,u)=0\) for every lower
triple \(u\), even when the first input has a nonzero D component. These
are the ordinary normalization conditions of a stacking operation.

## 2. Formulas used

On a full defining system, that is, with \(p_6(A,B)=0\) and
\(\delta C+\tau'_3(A;B)=0\), the correction \(J_6\) of (A5) is the
repository's integral tertiary representative:

\[
J_6(A,B,C)=\mathcal K_6(A,B,C)=\delta_s\Omega_6(A,B,C)=T_3(A;B,C),
\qquad
\Omega_6=\mathsf h(E_5C)+\widehat R_3(A,B)+\tfrac23\widetilde{P^1_s\rho_3A},
\]

by (A7) and formula (T) of [tertiary_operations.md](tertiary_operations.md).
Its prime-three term is

\[
2\beta_{3,s}P^1_s\rho_3A=\frac23\,\delta_s\widetilde{P^1_s\rho_3A},
\tag{O3}
\]

with the \(0,1,2\) lift and the nineteen-term reduced power (P3) of
[tertiary_operations.md](tertiary_operations.md), Section 7. At \(s=0\),
(O3) is \(2\beta_{\mathbf Z,3}P^1\rho_3A=\tfrac23\delta\widetilde{P^1\rho_3A}\).
The remaining phase \(\mathsf h(E_5C)+\widehat R_3(A,B)\) is dyadic: its
denominators are powers of two (Section 6 of
[tertiary_operations.md](tertiary_operations.md)). The coefficient two and
the reduced power are listed in [mathematical-status.md](mathematical-status.md).
The argument below uses these displayed formulas only. It does not certify
the whole implemented spectrum.

## 3. A nonzero degree-six class on \(B((\mathbf Z/3)^2)\)

Take \(G=(\mathbf Z/3)^2\), \(X=BG\), \(s=\omega=0\). The mod-three
cohomology ring is

\[
H^*(BG;\mathbf F_3)=\mathbf F_3[v_1,v_2]\otimes\Lambda(u_1,u_2),
\qquad |u_i|=1,\quad |v_i|=2.
\]

Distinguish the two Bocksteins:

- \(\beta_{\mathbf Z,3}:H^m(-;\mathbf F_3)\to H^{m+1}(-;\mathbf Z)\) is the
  connecting map of \(0\to\mathbf Z\xrightarrow{3}\mathbf Z\to\mathbf F_3\to0\),
  represented by \(x\mapsto\delta\widetilde x/3\) with the \(0,1,2\) lift.
  This is \(\beta_{3,s}\) of (O3) at \(s=0\).
- \(\beta_{\mathbf F_3}=\rho_3\beta_{\mathbf Z,3}\) is the mod-three
  Bockstein, with \(\beta_{\mathbf F_3}(u_i)=v_i\) and
  \(\beta_{\mathbf F_3}(v_i)=0\).

Let \(U_i=\beta_{\mathbf Z,3}(u_i)\), so that \(\rho_3U_i=v_i\), and let

\[
z=\beta_{\mathbf Z,3}(u_1u_2)\in H^3(BG;\mathbf Z).
\]

The derivation rule gives

\[
\rho_3z=v_1u_2-u_1v_2.
\]

Choose an integral cocycle \(Z\) representing \(z\) and take

\[
\boxed{x=(A,B,C,D)=(4Z,0,0,0)\quad\text{at }k=6.}
\tag{O4}
\]

One explicit cochain representative is
\(Z=\delta\bigl(\widetilde{u_1\smile u_2}\bigr)/3\), where the mod-three
cup product is lifted to the values \(0,1,2\). The coboundary of this lift
is divisible by three because \(u_1\smile u_2\) is a mod-three cocycle.

The integral cochain \(A\) is divisible by four. Hence \(\rho A=0\), and
the binary carry \(\rho\bigl((A-\widetilde{\rho A})/2\bigr)=\rho(2Z)=0\).
In (S11) of the formula sheet every term of \(\tau'_3(A;B)\) vanishes: the
inputs \(a=\rho A\) and \(t=\rho K\), the defining cochain \(b\) (here the
layer B) and \(B(a)=d\widetilde a/2\) vanish. Also
\(\delta_sA=4\delta Z=0\) and \(P_6(A)=Q_D(\rho A)=0\). Consequently \(x\)
is a full defining system and

\[
\mathfrak d_6x=(0,0,0,T),\qquad T=g_6(4Z,0,0)=J_6(4Z,0,0)=T_3(4Z;0,0).
\]

The same holds for \(\mathfrak d^0_6x\), since (O1) and (O2) agree on
\(L_6\).

Since \(4\equiv1\pmod3\), \(\rho_3A=\rho_3z\). At \(s=0\) the cochain
formula (P3) induces the reduced power \(P^1\) on cohomology. It is the
cyclic-diagonal operation of
[../python/mod3_power.py](../python/mod3_power.py), whose normalization
gives \(P^1(uv)=uv^3\) on \(B\mathbf Z/3\). A unit multiple would not
change the argument. By the Cartan formula and instability
(\(P^1u_i=0\), \(P^1v_i=v_i^3\)),

\[
P^1\rho_3A=v_1^3u_2-u_1v_2^3.
\]

Since \(\rho_3\beta_{\mathbf Z,3}=\beta_{\mathbf F_3}\), reducing the
three-primary term (O3) modulo three gives, with the vanishing of the
dyadic class shown below,

\[
\boxed{\rho_3[T]=2\bigl(v_1^3v_2-v_1v_2^3\bigr)\ne0.}
\tag{O5}
\]

Equivalently, the three-primary integral class is

\[
[T]_{(3)}=2\bigl(U_1^3U_2-U_1U_2^3\bigr)\ne0
\quad\text{in }H^8(BG;\mathbf Z).
\]

The two polynomial monomials are distinct and nonzero. The dyadic phase
contributes no integral class on this group. Both summands of
\(T=\delta\bigl(\mathsf h(E_5C)+\widehat R_3\bigr)+2\beta_{\mathbf Z,3}P^1\rho_3A\)
are integral, the second because \(\rho_3A\) is a mod-three cocycle. If
\(2^N\) clears the denominators of the dyadic phase, then \(2^N\) kills
the class of its coboundary. Positive-degree integral cohomology of a
finite group of odd order is killed by \(|G|=9\). These annihilators are
coprime, so that class is zero. This proves the nonvanishing without a
numerical approximation and without a cochain chosen by solving an
equation.

## 4. An input outside \(L_6\)

Let \(a=(1,0)\) and \(b=(0,1)\) in \(G\). Define a normalized
inhomogeneous group-bar 4-cochain \(B_{\rm off}\) by

\[
B_{\rm off}(g_1,g_2,g_3,g_4)=
\begin{cases}1,&(g_1,g_2,g_3,g_4)=(a,a,a,a),\\0,&\text{otherwise}.
\end{cases}
\]

Then

\[
\delta B_{\rm off}(a,a,a,a,b)=1\pmod2.
\]

The last face \((a,a,a,a)\) is the only nonzero summand. Merging adjacent
copies of \(a\) produces \(2a\), and merging the last \(a\) with \(b\)
produces \(a+b\). Put

\[
y=(0,B_{\rm off},0,0).
\]

Then \(p_6(y)=(0,\delta B_{\rm off})\ne0\), so \(y\notin L_6\).
Formula (O2) gives \(g^0_6(y)=0\), so the D component of
\(\mathfrak d^0_6y\) is zero. Formula (O1) instead gives
\(g_6(y)=T_3(\mathcal R_6(y))\).

## 5. The contradiction

The product \(x\times_6y\) has A component \(4Z\) and B component
\(B_{\rm off}+\alpha_6(4Z,0)=B_{\rm off}\). It therefore lies outside
\(L_6\) as well, whatever its C component. Hence
\(g^0_6(x\times_6y)=0\).

Write \(\gamma=\gamma_6(u_x,u_y)\) for the lower triples \(u_x,u_y\) of
\(x,y\). Both D components vanish, so the D component of
\(\mathfrak d^0_6(x\times_6y)\) is \(\delta\gamma\).

On the other hand \(\mathfrak d^0_6x=(0,0,0,T)\). The correction
\(\gamma_7\) depends only on the lower triples and is normalized at a zero
input. The D component of
\(\mathfrak d^0_6x\times_7\mathfrak d^0_6y\) is therefore
\(T+g^0_6(y)=T\). Strict compatibility would force

\[
\boxed{\delta\gamma=T.}
\]

This contradicts (O5). The degrees match: \(\gamma\) is an integral
7-cochain and \(T\) is an integral 8-cocycle. This proves the Proposition.

## 6. Consequence: \(g_6\) keeps the \(T_3\) obstruction outside \(L_6\)

The argument uses three facts: the value of \(g_6\) on the full defining
system \(x\), its values on \(y\) and on \(x\times y\) outside \(L_6\), and
the normalization. For any differential with the same lower components and
the same \(J_6\) on \(L_6\), and any normalized triangular product, the D
component of the compatibility identity at the pair \((x,y)\) reads

\[
g_6(x\times_6y)-g_6(y)-J_6(x)=-\delta_s\gamma_6(u_x,u_y),
\tag{O6}
\]

since \(\gamma_7(0,\cdot)=0\). A compatible extension therefore keeps the
nonzero \(T_3\) obstruction when a full defining system is translated by
an arbitrary lower input outside \(L_6\). In particular, adding the
corrections \(\alpha,\beta,\gamma\) to the differential with (O2) cannot
satisfy the identity.

The section branch of (O1) satisfies (O6). Use the notation of
[G6_GENERAL_SECTION.md](G6_GENERAL_SECTION.md), Sections 1–3 and 5: the
involution
\(\phi(A,B,C)=(A,B,C+e(A,B))\), with \(e=0\) on \(L_6\) and
\(e=q_{\rm loc}\) elsewhere, the natural product \(*\), the pointed
section \(\sigma\) of \(\widehat F\) with \(\sigma(0)=0\), and the natural
retraction \(\mathcal R(u)=\sigma(\widehat Fu)\backslash u\), so that
\(\mathcal R_6=\mathcal R\circ\phi\) (Section 3 of
[all_cochain_differential.md](all_cochain_differential.md)). Work in
natural coordinates with the inputs \(u=x\) and
\(v=\phi y=(0,B_{\rm off},q_{\rm loc}(0,B_{\rm off}))\). Since \(x\) lies
in \(L_6\), \(\phi x=x\); since \(\sigma(0)=0\) and the zero triple is the
unit of \(*\), \(\mathcal R(x)=x\). Also \(\widehat F(\phi y)=F_6y\) and
\(\mathcal R(\phi y)=\mathcal R_6y\). The reference product is
\(q=\sigma(0)*\sigma(F_6y)=\sigma(F_6y)\), so its comparison triple
\(\kappa=\sigma(\widehat Fq)\backslash q\) is zero and
\(\gamma_7(0,F_6y)=T_3(\kappa)=0\). The comparison terms \(j\) and
\(\mathcal Q(\cdot)\) of the mixed product vanish on this pair. The
two-cylinder identity
\(\delta_sG=T_3(z)+T_3(z')-T_3\bigl(\mathcal R(u*v)\bigr)+T_3(\kappa)\), with
\(z=\mathcal R(u)\) and \(z'=\mathcal R(v)\), becomes

\[
\delta_sG=T+T_3(\mathcal R_6y)-T_3\bigl(\mathcal R(x*\phi y)\bigr),
\]

which is the natural-coordinate form of (O6) with \(\gamma_6=G\). Section 6
("Return to the current product") of
[G6_GENERAL_SECTION.md](G6_GENERAL_SECTION.md) adds the integral
transport \(L(y',\lambda_{\rm loc})\), with
\(y'=\mathcal R_6(x\times_6y)\), and obtains (O6) for the product
\(\times_6\), where \(g_6(x\times_6y)=T_3\bigl(\mathcal R_6(x\times_6y)\bigr)\)
and \(\gamma_6=G+L(y',\lambda_{\rm loc})\). The off-\(L_6\) values of
\(g_6\) on \(y\) and on \(x\times_6y\) thus differ, up to the coboundary
\(\delta_s\gamma_6\), by \(T\), whose prime-three class \([T]\) is nonzero,
and there is no contradiction.

The explicit pure-C product
([../python/stacking_model/compatible_sector.py](../python/stacking_model/compatible_sector.py),
`compatible_sector.pure_c_g` and `compatible_sector.pure_c_gamma`) is not
affected by this obstruction. Inputs with \(A=B=0\) lie in \(L_k\) and
never reach the branch outside \(L_6\).

## 7. Compatibility identities for any differential

Let \(u=(A,B,C)\) and \(u'=(A',B',C')\), and let

\[
F_k(u)=(a,b,c)=\bigl(\delta_sA,\ \delta B+P_k(A),\ \delta C+f_k(A,B)\bigr),
\]

with \(F_k(u')=(a',b',c')\). Write
\(\mu_k(u,u')=(M_A,M_B,M_C)\) for the lower product of Section 1. The
homomorphism \(\mathfrak d_k(x\times_ky)=\mathfrak d_k(x)\times_{k+1}\mathfrak d_k(y)\)
holds for all \(D,D'\) if and only if

\[
\delta\alpha_k(A,A')+\alpha_{k+1}(a,a')
=P_k(A+A')+P_k(A)+P_k(A'),
\]

\[
\delta\beta_k(A,B;A',B')+\beta_{k+1}(a,b;a',b')
=f_k(M_A,M_B)+f_k(A,B)+f_k(A',B'),
\]

\[
\delta_s\gamma_k(u,u')-\gamma_{k+1}(F_ku,F_ku')
=g_k(u)+g_k(u')-g_k(\mu_k(u,u')).
\]

The first two equalities hold over \(\mathbf F_2\). The third is integral.
The first two are the lower homomorphism
\(F_k\mu_k=\mu_{k+1}(F_k,F_k)\), and the third is the upper condition of
Section 5 ("The two identities") of
[ALL_COCHAIN_STACKING.md](ALL_COCHAIN_STACKING.md).

## 8. The rational-phase construction and its integrality condition

Suppose the differential has the telescoping form (A25),

\[
g_k(u)=\delta_s\Phi_k(u)-\Phi_{k+1}(F_ku),\qquad
\Phi_k(u)\in C^{k+1}(X;\mathbf Q_s),\qquad F_{k+1}F_k=0,
\]

and the first two compatibility identities hold. Define

\[
\Delta_\mu\Phi_k(u,u')=\Phi_k(u)+\Phi_k(u')-\Phi_k(\mu_k(u,u')).
\]

The lower homomorphism gives
\(\delta_s\Delta_\mu\Phi_k=g_k(u)+g_k(u')-g_k(\mu_k(u,u'))
+\Delta_\mu\Phi_{k+1}(F_ku,F_ku')\). Hence
\(\gamma_k=\Delta_\mu\Phi_k+h_k\) obeys the third compatibility identity
exactly when

\[
\delta_sh_k(u,u')=h_{k+1}(F_ku,F_ku').
\]

One explicit sufficient form is

\[
\boxed{\gamma_k=\Delta_\mu\Phi_k+
\delta_se_k(u,u')+e_{k+1}(F_ku,F_ku'),}
\tag{O7}
\]

where \(e_k\in C^k(X;\mathbf Q_s)\) and \(e_j(0,0)=0\). The successor term
has a plus sign. The identity follows from \(F^2=0\) and
\(\delta_s^2=0\). One must separately prove that the correction is
integral:

\[
\Delta_\mu\Phi_k+\delta_se_k+e_{k+1}(F_ku,F_ku')
\in C^{k+1}(X;\mathbf Z_s).
\tag{O8}
\]

For the differential \(\mathfrak d^0\) of (O2), Section 5 shows that no
normalized choice satisfies all these conditions. Formal phase
cancellation cannot bypass the integral class \([T]\).

For the differential \(\mathfrak d\) of (A13), with (O1), (A28) of
[all_cochain_differential.md](all_cochain_differential.md) writes every
\(g_k\), \(0\le k\le6\), in the form (A25). Outside \(L_6\) it uses
\(\Phi_6=\Omega_6(\mathcal R_6u)\) and \(\Phi_7(F_6u)=0\), so
\(g_6=T_3(\mathcal R_6u)\). The degree-six product \(\gamma_6\) is not
constructed from (O7); it is the two-cylinder construction of
[G6_GENERAL_SECTION.md](G6_GENERAL_SECTION.md).

The pure-C sector is an instance of (O7). At \(A=B=0\) the phase is
\(\Phi_k(0,0,C)=\mathsf h(E_{k-1}C)\), and with
\(e_k(C,C')=\mathsf h\bigl(h^2(C,C')\bigr)\), where \(h^2\) is the
polarization of [ALL_COCHAIN_STACKING.md](ALL_COCHAIN_STACKING.md),
Section 2.1 (`compatible_sector.polarization`, `phase_eval.polarization`),
which satisfies, for binary cochains \(x,y\) of degree \(r\),

\[
\delta h^2_r(x,y)+h^2_{r+1}(\delta x,\delta y)=Q^2(x+y)+Q^2(x)+Q^2(y),
\]

(A25) and (O7) become

\[
g_k(0,0,C)=\tfrac12\bigl(\delta_s\widetilde{E(C)}-\widetilde{E(\delta C)}\bigr),
\]

\[
\gamma_k=\tfrac12\bigl(\widetilde{E(C)}+\widetilde{E(C')}-\widetilde{E(C+C')}
+\delta_s\widetilde{h^2(C,C')}+\widetilde{h^2(\delta C,\delta C')}\bigr).
\]

These are `compatible_sector.pure_c_g` and `compatible_sector.pure_c_gamma`.
Their integrality (O8) follows from the polarization identity modulo two.

## 9. Relation to the implementation

The vendored stacking formulas in
[../python/stacking_model](../python/stacking_model) (their origin is
recorded in
[../python/stacking_model/provenance.json](../python/stacking_model/provenance.json))
evaluate (O1):

| Formula | Implementation |
| --- | --- |
| \(g_6\): \(J_6\) on \(L_6\), \(T_3(\mathcal R_6)\) outside | [../python/stacking_model/g6_general_section.py](../python/stacking_model/g6_general_section.py): `g6_general_section.CurrentSectionStacking.correction` (conjugation by \(\phi\), then `NaturalSectionStacking.correction`), `NaturalSectionStacking.correction` and `retract`, with the section `g6_general_section.FiniteCurrentSection` |
| \(\gamma_6\) and the successor \(\gamma_7\) | [../python/stacking_model/g6_general_section.py](../python/stacking_model/g6_general_section.py): `g6_general_section.NaturalSectionStacking.gamma` and `successor_gamma`, and `g6_general_section.CurrentSectionStacking.gamma` (natural \(\gamma_6\) plus the transport \(L(y',\lambda_{\rm loc})\)) |
| \(T_3\) with the prime-three term (O3), and \(J_6\) | [../python/stacking_model/production_g6_section.py](../python/stacking_model/production_g6_section.py): `production_g6_section.ProductionG6Section.T`, the signed coboundary of `upper_phase_diagnostic.production_phase(3,…)` ([../python/stacking_model/upper_phase_diagnostic.py](../python/stacking_model/upper_phase_diagnostic.py)), which adds `mod3_power.tertiary_three_primary_phase` ([../python/mod3_power.py](../python/mod3_power.py)) through `high_phase.build_phase` ([../python/high_phase.py](../python/high_phase.py)); and `ProductionG6Section.J`, which is `closed_ab_upper.J` ([../python/stacking_model/closed_ab_upper.py](../python/stacking_model/closed_ab_upper.py)) |
| \(\mathfrak d_6\) and \(\times_6\) on four cochains | [../python/stacking_model/four_cochain_stacking.py](../python/stacking_model/four_cochain_stacking.py): `four_cochain_stacking.Stacking` with `degree_six=ProductionG6Section` |
| \(\alpha_k\), pure-C \(g_k\) and \(\gamma_k\) | [../python/stacking_model/compatible_sector.py](../python/stacking_model/compatible_sector.py): `compatible_sector.alpha`, `compatible_sector.pure_c_g`, `compatible_sector.pure_c_gamma` |

The module
[../python/stacking_model/g6_repair.py](../python/stacking_model/g6_repair.py)
supplies the finite linear algebra of the section (`IntegralCycleCoordinates`,
`f2_solve`). Its class `g6_repair.G6Repair` is another retraction onto
full legal triples, by projection to the primary and secondary
admissibility lattices, with the same branch form, \(J_6\) on \(L_6\) and
\(T_3\) of its own retraction outside, whose values off \(L_6\) differ in
general from \(T_3\circ\mathcal R_6\). The section construction does not
call it.

The complete-bar reference model evaluates this section construction on
the normalized bar resolution. It consists of
[../gap/extension_bar.gi](../gap/extension_bar.gi) and
[../gap/extension_degree_six.gi](../gap/extension_degree_six.gi), with
[../python/extension_degree_six.py](../python/extension_degree_six.py)
(`extension_degree_six.configure_degree_six` attaches
`ProductionG6Section`). It is used only by tests, never by `koFull`. Its
section branch outside \(L_6\) is checked by
`test_nonzero_off_shell_A_and_successor` in
[../python/test_extension_degree_six.py](../python/test_extension_degree_six.py).
The opt-in [../tst/extension_degree_six.tst](../tst/extension_degree_six.tst)
and the opt-in comparison in
[../python/test_extension_transfer.py](../python/test_extension_transfer.py)
compare the reference model with the native model on flat states, that is,
on the \(J_6\) branch.

`koFull` uses the native resolution model
[../python/extension_transfer.py](../python/extension_transfer.py) with
[../python/extension_native_upper.py](../python/extension_native_upper.py).
In degree six it evaluates \(J_6\) and the legal \(\gamma_6\) on the legal
lower locus only (`extension_native_upper.g`, `extension_native_upper.gamma`).
The legal \(\gamma_6\) of two states whose A layers are both nonzero is
refused (`extension_native_upper.PairSourceLimit`, degree unresolved)
unless `FERMIONAHSS_DEGREE_SIX_A_STACKING=1`; see
[extensions.md](extensions.md). Neither the zero value (O2) nor the
section branch of (O1) is evaluated on the resolution. A request outside
\(L_6\) raises `extension_native_upper.SectionBranchRequired`, and the
degree is recorded as unresolved. Such a request does not arise. The
native curvature stops at the first nonzero layer, and products and gauge
actions receive flat states and differential images. Flat states lie in
\(L_6\) because \(p_6=0\) on them. A gauge of a degree-six state is a
degree-five state \(g\), acted on through \(\mathfrak d_5\); the lower pair
of \(\mathfrak d_5g\) is \(p_5(A_g,B_g)\), which lies in \(L_6\) because
\(p_6p_5=0\) ((A3), (A16) of
[all_cochain_differential.md](all_cochain_differential.md)). The
obstruction of this note concerns the strict identity on arbitrary
cochains, such as the input \(y\) of Section 4, which is neither flat nor
an image. It does not enter the groups that `koFull` computes.

For lower-legal triples with \(A=0\) in degrees five and six, the native
engine evaluates the D-layer curvature by the exact direct formula of
`a0_high_gamma` (`a0_high_gamma.HigherA0Stacking`,
[../python/stacking_model/a0_high_gamma.py](../python/stacking_model/a0_high_gamma.py)).
When also \(B=0\) this formula is `compatible_sector.pure_c_g` of
Section 8. For two fully legal triples with \(A=B=0\), the native product
uses `compatible_sector.pure_c_gamma` in place of `closed_ab_upper.gamma`
and keeps the half-lift carry and the K change; the two base corrections
differ by an integral coboundary that a D-gauge absorbs. With a nonzero B
the base is `closed_ab_upper.gamma`, with the same carry and K terms. See
[extensions.md](extensions.md), Section 10 of
[G6_GENERAL_SECTION.md](G6_GENERAL_SECTION.md) and
[mathematical-status.md](mathematical-status.md).
