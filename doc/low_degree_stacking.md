# Stacking in package degrees one and two

This note gives the differential and the stacking product in package
degrees one and two in the form the native extension engine
([extension_transfer.py](../python/extension_transfer.py)) and the
complete-bar reference model
([four_cochain_stacking.py](../python/stacking_model/four_cochain_stacking.py))
evaluate them. Both engines use the classes `DegreeOneCommutativeStacking`
and `DegreeTwoCommutativeStacking` of
[coherent_low_commutative.py](../python/stacking_model/coherent_low_commutative.py).
The differential is (A13)–(A14) of
[all_cochain_differential.md](all_cochain_differential.md). In degree one
it is (6) of
[dimension_indexed_differentials.md](dimension_indexed_differentials.md).
In degree two it is (7) there, but only on the locus \(L_2\) where
\(\delta B=0\). The products are the rows \(k=1,2\) of the table in
Section 3.2 of [ALL_COCHAIN_STACKING.md](ALL_COCHAIN_STACKING.md). This note
writes out their phases \(e_1,e_2,e_3\).

Cup-\(i\) products, the signed coboundary \(\delta_s\), binary lifts
\(\widetilde x\) and the half lift \(\mathsf h(x)=\widetilde x/2\) follow
[conventions.md](conventions.md). Each \(\mathsf h\) lifts its whole binary
argument once, and different occurrences of \(\mathsf h\) are lifted
separately. The ramp and the prism follow Section 2 of the all-cochain note:

- \(u\ell\) is the pullback of \(u\) to \(X\times I\), multiplied by the
  interval coordinate of the last vertex;
- \(s\) and \(\omega\) are pulled back without \(\ell\);
- \(I\) is the signed right prism (A8), and
  \(\delta_sI+I\delta_s=i_1^*-i_0^*\).

## States and differentials

In package degree k the state is \((A,B,C,D)\), with cochain degrees
\((k-3,k-2,k-1,k+1)\). Layers of negative degree are absent.

| k | State | Differential |
| ---: | --- | --- |
| 1 | \((C,D)\), \(C\in C^0(\mathbf F_2)\), \(D\in C^2(\mathbf Z_s)\) | \(\mathfrak d_1\) of (A14): \((\delta C,\ \delta_sD+\tfrac12(\delta_s\widetilde{\omega C}-\widetilde{\omega\delta C}))\) |
| 2 | \((B,C,D)\), \(B\in C^0(\mathbf F_2)\), \(C\in C^1(\mathbf F_2)\), \(D\in C^3(\mathbf Z_s)\) | \(\mathfrak d_2\) of (A13): \((\delta B,\ \delta C+H_2(B),\ \delta_sD+g_2(B,C))\) |

In these degrees the operations (A2) are

\[
E_0(x)=\omega\smile x,\qquad
E_1(x)=x\smile\delta x+\omega\smile x,\qquad
E_2(x)=x\smile x+x\smile_1\delta x+\omega\smile x,
\]

and (A10) is \(H_2(B)=Q_{D,0}(B)=\omega\smile B+s\smile(B\smile\delta B)\).
Write

\[
F_1(C)=(0,\ \delta C),\qquad F_2(B,C)=\bigl(\delta B,\ \delta C+H_2(B)\bigr)
\]

for the lower components of \(\mathfrak d_1\) and \(\mathfrak d_2\), as
lower pairs of the next degree; \(F_2\) is `a0_gamma.lower_d2`. These are
the \(F_k\) of Section 1 of
[ALL_COCHAIN_STACKING.md](ALL_COCHAIN_STACKING.md) with the zero
components left out. By (A3) and (A11), every \((b,c)=F_2(B,C)\) is a full
defining system of degree three at \(A=0\), which means \(\delta b=0\) and
\(\delta c=Q_D(b)=\tau_3(0,b)\).

Both \(g_1\) and \(g_2\) have the form (A25),
\(g_k=\delta_s\Phi_k-\Phi_{k+1}(F_k\,\cdot)\), with the potentials of
(A28).

**Degree one.** Here \(g_1=J_1\) of (A5), with \(\tau_1=0\) and
\(t_1=\delta C\):

\[
g_1(C)=\delta_s\mathsf h(\omega C)-\mathsf h(\omega\,\delta C).
\]

This is (A14) (`compatible_sector.pure_c_g`). In the form (A25) it reads
\(\Phi_1=\mathsf h(\omega C)\) and
\(\Phi_2(F_1C)=\Phi_2(0,\delta C)=\mathsf h(E_1(\delta C))=\mathsf h(\omega\,\delta C)\).

**Degree two.** Here \(g_2=J_2\) of (A5) when \(B\) is closed on the whole
complex, which is the locus \(L_2\) of (A4). Otherwise \(g_2=G_2\) of (A12).
Both branches use the potential

\[
\Phi_2(B,C)=
\begin{cases}
\mathsf h(C\smile\delta C+\omega C)+\mathsf h(t_2\smile_1\omega B),
  & \delta B=0,\quad t_2=\delta C+\omega B,\\
I\,\Omega_3\bigl(F_2(B\ell,C\ell)\bigr),
  & \delta B\ne0,
\end{cases}
\qquad
g_2=\delta_s\Phi_2-\Omega_3\bigl(F_2(B,C)\bigr)
\]

(`a0_gamma.phi2`, `a0_gamma.g2`). Here \(\Omega_3(b,c)\) is the
degree-three phase at \(A=0\) from the table in Section 1 of the
all-cochain note, for closed \(b\):

\[
\Omega_3(b,c)=\mathsf h(E_2c)+R_0(0,b).
\]

Its two terms are lifted separately (`a0_gamma.omega3`). At \(A=0\) the
additional terms of (R0e) in [tertiary_operations.md](tertiary_operations.md)
vanish, so \(R_0(0,b)=R_{\partial,1}(b)=\Theta_1(b,0)\) of (B) there. The
function `a0_gamma.theta_degree_one` evaluates (B) at \(m=1\), \(v=0\),
where \(\chi_1=0\), \(\operatorname{Sq}^2b=0\) and
\(B_b=\delta\widetilde b/2=\widetilde{b\smile b}\).

The last term of \(g_2\) is always \(\Omega_3(F_2(B,C))\), which is
\(\Phi_3(F_2(B,C))\) because \(t_3=0\) on a full defining system. Only the
potential depends on the branch:

- **On \(L_2\):** \(F_2(B,C)=(0,t_2)\), because \(s\smile(B\smile\delta B)=0\).
  Since \(\delta t_2=\omega\smile\delta B=0\) and \(R_0(0,0)=0\), we get
  \(\Omega_3(0,t_2)=\mathsf h(t_2\smile t_2+\omega t_2)\). So \(g_2\) is (7).
- **Off \(L_2\):** put
  \((b_I,c_I)=F_2(B\ell,C\ell)=\bigl(\delta(B\ell),\ \delta(C\ell)+H_2(B\ell)\bigr)\),
  a full defining system on \(X\times I\). Prism Stokes gives
  \(\delta_s\Phi_2=\Omega_3(F_2(B,C))-\Omega_3(0,0)-I\,\delta_s\Omega_3(b_I,c_I)\),
  and \(\Omega_3(0,0)=0\), so

  \[
  g_2=-I\,\delta_s\Omega_3(b_I,c_I)=-I\,\mathcal K_3(b_I,c_I)=G_2.
  \]

The branch of \(g_2\), and of each potential in the product below, depends
on whether \(\delta B\) vanishes on the whole complex. The reference model
tests this on the complete bar complex. The native engine uses the zero
test on the resolution.

## Products

The product is triangular. All values are exact rationals, and
\(\gamma_1\) and \(\gamma_2\) are checked to be integral on evaluation
(`compatible_sector.integral`).

### Degree one

\[
(C,D)\times(C',D')=\bigl(C+C',\ D+D'+\gamma_1(C,C')\bigr),
\]
\[
\gamma_1=\tfrac12\bigl(\widetilde{E_0(C)}+\widetilde{E_0(C')}-\widetilde{E_0(C+C')}\bigr)
+\delta_se_1(C,C')+e_2\bigl(0,\delta C;\,0,\delta C'\bigr).
\]

The degree-one phase is

\[
e_1(C,C')=\mathsf h(C\smile\delta C'),\qquad
e_1(v_0v_1)=\tfrac12\,C(v_0)\bigl[C'(v_0)+C'(v_1)\bmod2\bigr]
\]

(`DegreeOneCommutativeStacking.phase`). The successor term is the
degree-two phase \(e_2\) below, evaluated on \(F_1C=(0,\delta C)\) and
\(F_1C'=(0,\delta C')\). A closed form of \(\gamma_1\) is given
at the end of this section.

### Degree two

\[
(B,C,D)\times(B',C',D')=
\bigl(B+B',\ C+C'+\beta_2(B,B'),\ D+D'+\gamma_2\bigr),
\qquad
\beta_2(B,B')=\delta B\smile B'+s\smile(B\smile B'),
\]
\[
\gamma_2=\Phi_2(B,C)+\Phi_2(B',C')-\Phi_2\bigl(B+B',\,C+C'+\beta_2(B,B')\bigr)
+\delta_se_2(B,C;B',C')+e_3\bigl(F_2(B,C);\,F_2(B',C')\bigr).
\]

Here \(\beta_2\) is `a0_gamma.beta2`, the lower correction of Section 2 of
[ALL_COCHAIN_STACKING.md](ALL_COCHAIN_STACKING.md). The three potentials
take their branches from the closedness of \(B\), \(B'\) and \(B+B'\),
respectively.

### The phases \(e_3\) and \(e_2\)

Let \((b,c)\) and \((b',c')\) be full defining systems of degree three at
\(A=0\). Write

\[
\beta_3(b,b')=b\smile b'+s\smile(b\smile_1b')
\]

for the C-layer correction \(\beta_3\) at \(A=A'=0\)
(`a0_gamma.beta3_closed`). The degree-three phase is
`paper_commutative_production.phase`:

\[
e_3(b,c;b',c')=e^{\rm rep}_3(b,c;b',c')+\kappa\bigl(b+b',\,c+c'+\beta_3(b,b')\bigr)-\kappa(b,c)-\kappa(b',c')
+\mathsf h\bigl(\omega\smile(b\smile_1b')\bigr),
\]
\[
\kappa(b,c)=\tfrac18\,\widetilde b\smile\widetilde b\smile\widetilde b
-\mathsf h\bigl(s\smile c+c\smile_2\delta c\bigr)
+\mathsf h\bigl((b+s)\smile c+(b\smile_1s)\smile(b\smile b)+(b\smile_1\omega)\smile b\bigr).
\]

The coordinate gauge \(\kappa\) is
`paper_commutative_production.coordinate_gauge`; the last two terms of its
final \(\mathsf h\) are `paper_commutative_production.reorder_gauge`. The
term \(e^{\rm rep}_3\) is the phase of
`paper_stacking.repaired_paper_corrections` in degree three. That is the
degree-three phase of [arXiv:2310.19058](https://arxiv.org/abs/2310.19058)
as transcribed in these cup conventions (`paper_stacking.paper_corrections`),
plus the local correction `paper_stacking.antiunitary_phase_repair`, which
is not taken from the paper; see Section 7 of
[ALL_COCHAIN_STACKING.md](ALL_COCHAIN_STACKING.md). The degree-three
product at \(A=0\) is

\[
\gamma_3=\Omega_3(b,c)+\Omega_3(b',c')-\Omega_3\bigl(b+b',\,c+c'+\beta_3(b,b')\bigr)+\delta_se_3,
\]

which is integral on full defining systems
(`paper_commutative_production.gamma`). On these data it is the degree-three
product of Section 3.4 of [ALL_COCHAIN_STACKING.md](ALL_COCHAIN_STACKING.md),
where \(e_3\) is written \(e^{\rm pap}_3\). The successor \(e_3\) of
\(\gamma_2\) is therefore the phase of the degree-three product itself.
Here \(e_3\) is evaluated only on full defining systems: \(F_2(B,C)\) on
\(X\), and \((b_I,c_I)\) on \(X\times I\).

The degree-two phase is the prism transport of \(e_3\)
(`DegreeTwoA0Stacking.phase`):

\[
e_2(B,C;B',C')=-I\Bigl[e_3\bigl(b_I,c_I;\,b'_I,c'_I\bigr)
+\mathsf h\bigl(E_1(\eta)+c_I^{\rm tot}\smile_1\delta\eta
+\delta c_I^{\rm tot}\smile_2\delta\eta\bigr)\Bigr].
\]

The ingredients are:

- \((b_I,c_I)=F_2(B\ell,C\ell)\) and \((b'_I,c'_I)=F_2(B'\ell,C'\ell)\);
- \(\eta=\beta_2(B\ell,B'\ell)+\beta_2(B,B')\ell\);
- \(c_I^{\rm tot}=c_I+c'_I+\beta_3(b_I,b'_I)=c_I+c'_I+b_I\smile b'_I+s\smile(b_I\smile_1b'_I)\);
- \(E_1(\eta)=\eta\smile\delta\eta+\omega\smile\eta\).

The \(\mathsf h\) term is one lift of its whole binary argument. It is the
C-transport `a0_gamma.phase3_C_transport`, which moves the successor's C
argument from \(c_I^{\rm tot}\) to \(c_I^{\rm tot}+\delta\eta\). To see
why, write \(u=(B,C)\) and \(v=(B',C')\) for the lower pairs, and, as in
Section 1 of [ALL_COCHAIN_STACKING.md](ALL_COCHAIN_STACKING.md), write the
lower products as

\[
\mu_2(u,v)=\bigl(B+B',\ C+C'+\beta_2(B,B')\bigr),\qquad
\mu_3\bigl((b,c),(b',c')\bigr)=\bigl(b+b',\ c+c'+\beta_3(b,b')\bigr).
\]

The product of the two ramps, \(\mu_2(u\ell,v\ell)\), and the ramp of the
product, \(\mu_2(u,v)\ell\), have the same B layer, and their C layers
differ by \(\eta\). Also \(F_2\mu_2(u,v)=\mu_3(F_2u,F_2v)\). So \(F_2\)
takes the first of these products to C layer \(c_I^{\rm tot}\) and the
second to \(c_I^{\rm tot}+\delta\eta\).

### Closed forms in degree one

The terms of \(\gamma_1\) have closed forms.

**The first term.** Since
\(E_0(C)(v_0v_1v_2)=\omega(v_0v_1v_2)\,C(v_2)\), the first term of
\(\gamma_1\) is the lift \(\widetilde{\omega\smile C\smile C'}\), where
\(C\smile C'\) is the pointwise product.

**The successor term.** Take \(b=b'=0\) and closed binary \(c,c'\) of
degree two. Then the gauges \(\kappa\) cancel, and every other term of
\(e_3\) contains \(b\), \(b'\), \(\delta c\), \(\delta c'\) or
\(\beta_3=0\). So

\[
e_3(0,c;0,c')=\mathsf h(c\smile_1c').
\]

Now take closed binary \(a,a'\) of degree one. In \(e_2(0,a;0,a')\) the
C-transport vanishes (\(\eta=0\)), and \(c_I=\delta(a\ell)=a\smile\delta\ell\)
modulo two. Of the three prism simplices over \((v_0v_1v_2)\), only
\(((v_0,0),(v_1,0),(v_2,0),(v_2,1))\) contributes, with value
\(\mathsf h(a(v_0v_1)\,a'(v_1v_2))\). Hence

\[
e_2(0,a;0,a')=-\mathsf h(a\smile a').
\]

**The resulting formula.** Therefore

\[
\gamma_1=\widetilde{\omega\smile C\smile C'}+\delta_s\mathsf h(C\smile\delta C')
-\mathsf h(\delta C\smile\delta C').
\]

The last two terms are integral together, since
\(\delta(C\smile\delta C')=\delta C\smile\delta C'\) modulo two. On closed
\(C,C'\), \(\gamma_1=\widetilde{\omega\smile C\smile C'}\).

**The phase \(e_1\).** The same evaluation over an edge shows that \(e_1\)
is the prism transport of its successor:

\[
e_1(C,C')=-I\,e_2\bigl(0,\delta(C\ell);\,0,\delta(C'\ell)\bigr).
\]

Only the first prism simplex \(((v_0,0),(v_0,1),(v_1,1))\) contributes.
This is the construction of \(e_2\) from \(e_3\), without a C-transport
term because \(\beta\) vanishes in degree one. The code evaluates the
closed form. The prism form is `DegreeOneA0Stacking.phase_from_successor`.

## Implementation

Both engines construct `DegreeTwoCommutativeStacking()` and
`DegreeOneCommutativeStacking(successor=...)` from it.

- **Degree one.** The class overrides the phase with \(e_1\) and inherits
  `gamma` and `g` from `DegreeOneA0Stacking` in
  [a0_degree1.py](../python/stacking_model/a0_degree1.py).
- **Degree two.** The class inherits `phase` (\(e_2\)), `gamma` and `g` from
  `DegreeTwoA0Stacking` in
  [a0_degree2.py](../python/stacking_model/a0_degree2.py), with the successor
  `CommutativeLegalDegreeThree`. That successor evaluates
  [paper_commutative_production.py](../python/stacking_model/paper_commutative_production.py);
  it raises unless its three closedness flags are set; `DegreeTwoA0Stacking`
  sets them, since it evaluates the successor only on \(F_2\)-images, which
  are full defining systems.
- **Helpers.** [a0_gamma.py](../python/stacking_model/a0_gamma.py) supplies
  \(\beta_2\) (`beta2`), \(F_2\) (`lower_d2`), \(\beta_3\) (`beta3_closed`),
  \(\Phi_2\) (`phi2`, `raw_phi2`), \(g_2\) (`g2`), \(\Omega_3\) (`omega3`,
  `theta_degree_one`) and the C-transport (`phase3_C_transport`). The
  function \(g_1\) is `compatible_sector.pure_c_g`.
- **Not evaluated.** `DegreeTwoCommutativeStacking` replaces the default
  successor `a0_degree3.DegreeThreeA0Stacking` of `DegreeTwoA0Stacking` by
  `CommutativeLegalDegreeThree`, and `DegreeOneCommutativeStacking`
  overrides `a0_degree1.edge_phase`. Neither is evaluated in degrees one
  and two; `a0_degree3` enters the engines only in degree three, through
  the rank-zero case of the gauge \(q_3\)
  (`nonzero_degree3_comparison.gauge`).

In the native engine:

- `nonlinear` fills the D layer of \(\mathfrak d_1\) with \(g_1\), and the
  C and D layers of \(\mathfrak d_2\) with \(H_2(B)\) and \(g_2\);
- `bar_product` fills the C layer of the degree-two product with
  \(\beta_2\), and the D layers with \(\gamma_1\) and \(\gamma_2\).

In the reference model, the same calls are made by `Stacking.d` and
`Stacking.xtimes`.

## Use in `koFull`

The native engine of [resolution-extensions.md](resolution-extensions.md)
evaluates these formulas when `koFull` measures a relation of degree one
or two in the model. They hold for arbitrary binary cocycles \(s\) and
\(\omega\) on the resolution. The states have degree one and two, and the
gauges have degree zero and one. The comparison complex, zero test, flat
lifts, relation measurements and gauge comparisons are the same as in
degrees 3–6.

A two-primary relation whose target layer lies right below its
generator's layer is read from a primary operation on R instead
([extensions.md](extensions.md), "Primary-operation rows"). This covers
every relation of degree one (C over D) and, in degree two, B over C and
C over D. Such a relation goes to the model only on a fallback (listed in
`primaryFallbacks`), when a later relation needs it measured through D, or
with `FERMIONAHSS_NATIVE_RELATIONS=0`. Odd relations split in these
degrees. So in degree two the model is needed for B over D and for those
cases, and in degree one only on a fallback or with the variable set to 0.

- A degree-zero gauge is a single \(D\in C^1(X;\mathbf Z_s)\). By (A14),
  \(\mathfrak d_0(D)=(0,\delta_sD)\) as a degree-one state.
- A degree-one gauge \((C,D)\) of a degree-two state has
  \(\mathfrak d_1(C,D)=(0,\ \delta C,\ \delta_sD+g_1(C))\).
