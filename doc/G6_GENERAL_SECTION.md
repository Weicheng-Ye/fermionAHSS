# A degree-six section construction without associativity

This note constructs the degree-six branch of the stacking model outside
the legal lower locus \(L_6\). The branch has two parts. The first is the
D-layer differential \(g_6=T_3\circ\mathcal R_6\) of (A13) in
[all_cochain_differential.md](all_cochain_differential.md), Section 3. The
second is the product correction \(\gamma_6\) with its successor
\(\gamma_7\) on differential images. Together they satisfy the strict
upper identity of [ALL_COCHAIN_STACKING.md](ALL_COCHAIN_STACKING.md),

\[
\delta_s\gamma_6(u,v)=g_6(u)+g_6(v)-g_6(u\cdot v)+\gamma_7(F_6u,F_6v).
\]

The construction combines three pieces: an explicit natural lower product,
a fixed pointed section of the lower differential, and two legal
cylinders. It does not assume associativity of any lower product and it
does not solve an upper cochain equation. Only lower defining equations
are solved, and every upper term is an evaluation of a fixed formula.
[ALL_COCHAIN_OBSTRUCTION.md](ALL_COCHAIN_OBSTRUCTION.md) shows why \(g_6\)
must be nonzero outside \(L_6\): with \(g_6=0\) there, no normalized
unital triangular product satisfies
\(\mathfrak d(x\times y)=\mathfrak d x\times\mathfrak d y\) on all cochains.

The construction is implemented in three modules of the vendored stacking
model, whose origin is recorded in
[../python/stacking_model/provenance.json](../python/stacking_model/provenance.json):

- [../python/stacking_model/g6_general_section.py](../python/stacking_model/g6_general_section.py):
  the section and the section identity;
- [../python/stacking_model/lower_raw_comparison.py](../python/stacking_model/lower_raw_comparison.py):
  the natural lower product;
- [../python/stacking_model/production_g6_section.py](../python/stacking_model/production_g6_section.py):
  the binding to the fixed production operations.

The complete-bar reference model evaluates this branch, and only the tests
use that model. The native `koFull` engine evaluates degree six on the
legal lower locus only (Section 10).

Cup products \(\smile_i\), the signed coboundary \(\delta_s\), the reduction
\(\rho\), the lift \(\widetilde x\), the half lift
\(\mathsf h(x)=\widetilde x/2\) and the operations \(E\), \(Q_D\), \(Q^j\) are
those of [conventions.md](conventions.md) and (A2) of
[all_cochain_differential.md](all_cochain_differential.md). The prism
\(I\) is (A8) there, with \(\delta_sI+I\delta_s=i_1^*-i_0^*\). On
\(X\times I\), \(u\ell\) is the pullback of \(u\) multiplied by the
interval coordinate of the last vertex (the ramp, which vanishes at
\(i_0\) and restricts to \(u\) at \(i_1\)), and \(\pi^*u\) is the plain
pullback. Backgrounds \(s,\omega\) are pulled back. The prism of a cochain
pulled back from \(X\) vanishes, because its projected simplices are
degenerate.

## 1. Degree-six data

In package degree \(k=6\) a state is \((A,B,C,D)\) with
\(A\in C^3(X;\mathbf Z_s)\), \(B\in C^4(X;\mathbf F_2)\),
\(C\in C^5(X;\mathbf F_2)\) and \(D\in C^7(X;\mathbf Z_s)\). The lower
triple \(u=(A,B,C)\) has the lower differential, the first three
components of (A13),

\[
F_6(A,B,C)=\bigl(\delta_sA,\ \delta B+P_6(A),\ \delta C+f_6(A,B)\bigr),
\qquad
f_6(A,B)=\begin{cases}\tau'_3(A;B),&(A,B)\in L_6,\\ H_6(A,B),&(A,B)\notin L_6,\end{cases}
\]

where \(P_6(A)=Q_D(\rho A)\) and \(H_6\) is the prism (A9). Its values lie
in \(C^4(\mathbf Z_s)\times C^5(\mathbf F_2)\times C^6(\mathbf F_2)\) and
are legal in degree seven, since \(F_7F_6=0\). A triple has **legal AB**
when \((A,B)\in L_6\), that is, when the first two components of \(F_6\)
vanish as cochains on all of \(X\). It is a **full legal triple** when
\(F_6u=0\). An image \(w=(a,b,c)\) is **pure** when \(a=b=0\). Thus
\(F_6u\) is pure exactly when \(u\) has legal AB. The corrections have
degrees \(g_6\in C^8\), \(\gamma_6\in C^7\) and \(\gamma_7\in C^8\).

The lower operations used below are the following. The code column names
each function with its module of the vendored stacking model; most of them
are in [../python/stacking_model/off_shell_beta.py](../python/stacking_model/off_shell_beta.py)
and [../python/stacking_model/lower_raw_comparison.py](../python/stacking_model/lower_raw_comparison.py).

| Symbol | Definition | Code |
| --- | --- | --- |
| \(q_{\rm loc}(A,B)\) | the fixed comparison gauge \(h\) with \(\delta h=\tau'_3+H_6\) on each simplex all of whose faces have vanishing curvature \((\delta_sA,\delta B+P_6A)\), and zero elsewhere | `off_shell_beta.local_comparison_extension`, `h_tau_primitive.comparison_gauge` |
| \(\widehat f(A,B)\) | \(H_6(A,B)+\delta q_{\rm loc}(A,B)\), a natural cochain equal to \(f_6=\tau'_3\) on \(L_6\), that is, a natural extension of \(\tau'_3\) from \(L_6\) | `off_shell_beta.natural_f` |
| \(e(A,B)\) | \(0\) on \(L_6\), \(q_{\rm loc}(A,B)\) outside | `off_shell_beta.coordinate_gauge` |
| \(\phi\) | \(\phi(A,B,C)=(A,B,C+e(A,B))\), an involution | `g6_general_section.FiniteCurrentSection.conjugate` |
| \(\widehat F\) | \(F_6\circ\phi=\bigl(\delta_sA,\ \delta B+P_6(A),\ \delta C+\widehat f(A,B)\bigr)\), so \(F_6=\widehat F\circ\phi\) | `g6_general_section.FiniteCurrentSection.natural_lower` |
| \(\alpha(A,A')\) | \(h^D_3(\rho A,\rho A')\), the binary cross term of \(Q_D\) ([ALL_COCHAIN_STACKING.md](ALL_COCHAIN_STACKING.md), Section 2.1) | `off_shell_beta.alpha` |
| \(\beta^{\rm leg}\) | the calibrated legal lower correction of input degree three | `stacking_lower.legal_beta` |
| \(\beta_{\rm raw}\) | the prism of the closed pair residual, below | `off_shell_beta.raw_beta` |
| \(\lambda_3\), \(\lambda_{\rm loc}\) | the comparison of \(\beta^{\rm leg}\) and \(\beta_{\rm raw}\), below | `lower_raw_comparison.primitive`, `lower_raw_comparison.local_extension` |
| \(\widehat\beta\) | \(\beta_{\rm raw}+\delta\lambda_{\rm loc}\) | `lower_raw_comparison.natural_beta` |
| \(T_3\) | \(\delta_s\Omega_6\) on full legal triples, including the coefficient-two prime-three term | `upper_phase_diagnostic.production_phase` followed by \(\delta_s\) |
| \(J_6\) | (A5) on legal AB, with \(C\) arbitrary | `closed_ab_upper.J` |

With \(A_S=A+A'\) and \(B_S=B+B'+\alpha(A,A')\), the closed pair residual
is

\[
\mathrm{res}(A,B;A',B')=\widehat f(A_S,B_S)+\widehat f(A,B)+\widehat f(A',B')
+\beta^{\rm leg}_7\bigl(p_6(A,B),p_6(A',B')\bigr),
\]

where \(p_6\) is the first two components of \(F_6\) and \(\beta^{\rm leg}_7\) is the
legal correction of input degree four. The raw correction is its prism
over the ramped inputs on \(X\times I\), with pulled-back \(s,\omega\):

\[
\beta_{\rm raw}(A,B;A',B')=I\,\mathrm{res}(A\ell,B\ell;A'\ell,B'\ell)
\]

(`off_shell_beta._closed_pair_residual`, `off_shell_beta.raw_beta`). The
current lower product of
[ALL_COCHAIN_STACKING.md](ALL_COCHAIN_STACKING.md) is

\[
u\cdot v=\bigl(A+A',\ B+B'+\alpha(A,A'),\ C+C'+\beta_6(A,B;A',B')\bigr),
\qquad
\beta_6=\begin{cases}\beta^{\rm leg},&(A,B),(A',B')\in L_6,\\
\beta_{\rm raw}+e(A,B)+e(A',B')+e(A_S,B_S),&\text{otherwise}\end{cases}
\]

(`stacking_lower.all_cochain_beta`, `off_shell_beta.beta`). The
comparison primitive is

\[
\lambda_3(A,B;A',B')=\lambda_A(A,A')+I\bigl[\beta_{\rm ext}+\beta_{\rm raw}\bigr]\bigl(\pi^*A,B\ell;\pi^*A',B'\ell\bigr),
\qquad
\delta\lambda_3=\beta^{\rm leg}+\beta_{\rm raw}\ \text{on legal pairs}.
\]

Here \(\beta_{\rm ext}\) (`lower_raw_comparison.extended_beta`) is the
lower correction for \(\widehat f\) on signed-closed A and arbitrary B:
the closed-A correction \(\beta^\sharp\) (`closed_a_upper.beta_sharp`)
plus the lower gauges `nonzero_offshell_comparison.lower_gauge` of the two
inputs and of their sum, plus \(b\smile_5b'\), with
\(b=\delta B+P_6(A)\) and \(b'=\delta B'+P_6(A')\). \(\lambda_A\)
(`lower_raw_comparison.A_primitive`) is the fixed shared-sign pair
contraction of the A-only source. The certificate of \(\lambda_A\) is
`relative_acyclicity` and `small_source_values`: the relative universal
complex in degree five has no lower cells, two middle cells and twelve
upper cells, it is acyclic modulo two, and the A-only source vanishes on
both middle cells. No coefficient is fitted on \(X\). The local form
\(\lambda_{\rm loc}\) equals \(\lambda_3\) on each simplex all of whose faces
have vanishing curvature of both pairs, and zero elsewhere.

## 2. An explicit pointed section

Write a finite lower input as \((A,B,C)\) with integral \(A\) and binary
\(B,C\). Fix a unimodular integral basis of \(C^3(X;\mathbf Z_s)\) whose
first \(r\) columns span \(\ker(\delta_s:C^3\to C^4)\). The remaining
columns span a complement on which \(\delta_s\) is injective, which proves
that the first block is the whole integral kernel
(`g6_repair.IntegralCycleCoordinates`). Let \(Z\) be the matrix of the
first block. For an output \((a,b,c)\), solve \(\delta_sA_0=a\) in the
complement by exact elimination. Every integral solution is
\(A_0+Zt\). Enumerate \(t\in\{0,1,2,3\}^r\). For each \(A=A_0+Zt\), solve

\[
\delta B=b+P_6(A),
\]

and enumerate every \(B\) in the coset of \(\ker(\delta:C^4\to C^5)\). For
each such \(B\), solve

\[
\delta C=c+f_6(A,B).
\]

The first solution in this order is \(\sigma_0(a,b,c)\). The order is:
\(t\) lexicographic; \(B=B_0+\sum_i\epsilon_ik_i\) over the binary
counting of \(\epsilon\) in a fixed kernel basis; \(B_0\) and \(C\) the
reduced row-echelon solutions with free variables zero. Hence
\(\sigma_0(0)=0\). Only the lower defining equations are solved.

The enumeration is complete on \(\operatorname{im}F_6\). The secondary
formula \(\tau'_3(A;B)\) depends on \(A\) only through \(\rho A\) and the
carry \(\rho\bigl((A-\widetilde{\rho A})/2\bigr)\), so it is four-periodic in
\(A\). The same holds for \(H_6\): its suspended integral input changes by
a multiple of four and its suspended primary input does not change. The
global branch condition \((A,B)\in L_6\) is unchanged by adding four times
an integral cycle. Therefore the finite enumeration finds a preimage
whenever one exists. In the worst case it examines
\(4^r\cdot2^{\dim\ker(\delta:C^4\to C^5)}\) candidates, each with exact
lower linear solves.

The coordinate gauge \(e\) is not assumed four-periodic. The section is
conjugated instead:

\[
\sigma=\phi\circ\sigma_0.
\]

Then \(\widehat F\sigma=F_6\phi\phi\sigma_0=F_6\sigma_0\) is the identity on
\(\operatorname{im}F_6=\operatorname{im}\widehat F\), and \(\sigma(0)=0\)
(`FiniteCurrentSection.current_section` is \(\sigma_0\),
`FiniteCurrentSection.section` is \(\sigma\)).

## 3. The natural triangular lower product

In natural coordinates the lower product is

\[
u\ast v=\bigl(A+A',\ B+B'+\alpha(A,A'),\ C+C'+\widehat\beta(A,B;A',B')\bigr),
\qquad
\widehat\beta=\beta_{\rm raw}+\delta\lambda_{\rm loc}.
\]

It has three properties. First, on a pair of inputs with legal AB,
\(\lambda_{\rm loc}=\lambda_3\), so \(\widehat\beta=\beta^{\rm leg}\) and \(\ast\)
agrees with the calibrated product. Second, \(\widehat F\) is a
homomorphism: \(\widehat F(u\ast v)\) is the degree-seven lower product of
\(\widehat Fu\) and \(\widehat Fv\), because \(\delta\beta_{\rm raw}\) equals
the closed residual at \((A,B;A',B')\) (prism Stokes; the residual is
closed and vanishes at \(i_0\)) and \(\delta\lambda_{\rm loc}\) is a
coboundary. Third, \(\alpha\) and \(\widehat\beta\) vanish when either pair
\((A,B)\) is zero, so the zero triple is a two-sided unit. The product is
triangular, so left division \(v\backslash w=z\), the unique solution of
\(v\ast z=w\), is explicit:

\[
\begin{aligned}
A_z&=A_w-A_v,\\
B_z&=B_w+B_v+\alpha(A_v,A_z),\\
C_z&=C_w+C_v+\widehat\beta(A_v,B_v;A_z,B_z).
\end{aligned}
\]

The last two lines are binary. The section retraction is

\[
\mathcal R(u)=\sigma(\widehat Fu)\backslash u\quad\text{in natural coordinates},\qquad
\mathcal R_6(u)=\mathcal R(\phi u)=\sigma(\widehat F\phi u)\backslash\phi u\quad\text{in current coordinates},
\]

which is the \(\mathcal R_6\) of [all_cochain_differential.md](all_cochain_differential.md).
With \(z=\mathcal R(u)\), \(u=\sigma(\widehat Fu)\ast z\). The homomorphism
property gives \(\widehat Fu=\widehat Fu\cdot\widehat Fz\), and left cancellation
in the triangular output product gives \(\widehat F\mathcal R(u)=0\). The
retraction fixes every full legal triple, because \(\sigma(0)=0\). Neither
this argument nor the construction below uses associativity.

## 4. Two legal cylinders

For natural lower triples \(u,v\) write

\[
\begin{aligned}
&w=\widehat Fu,\quad w'=\widehat Fv,\quad a=\sigma(w),\quad b=\sigma(w'),\quad
z=\mathcal R(u),\quad z'=\mathcal R(v),\\
&q=a\ast b,\quad q_0=\sigma(\widehat Fq),\quad \kappa=q_0\backslash q,\quad
r=q\backslash(u\ast v).
\end{aligned}
\]

The homomorphism property gives \(\widehat F(u\ast v)=\widehat Fq\), so
\(\kappa\) and \(r\) are full legal triples. On \(X\times I\), form the two
explicitly parenthesized cylinders

\[
\begin{aligned}
r_\ell&=(a\ell\ast b\ell)\backslash\bigl[(a\ell\ast\pi^*z)\ast(b\ell\ast\pi^*z')\bigr],\\
k_\ell&=(q_0\ell)\backslash\bigl[\bigl((q_0\ell)\ast\pi^*\kappa\bigr)\ast\pi^*r\bigr].
\end{aligned}
\]

They are legal by the lower homomorphism and left cancellation. Since the
ramps vanish at \(i_0\) and the product is natural, their endpoints are

\[
i_0^*r_\ell=z\ast z',\qquad i_1^*r_\ell=r,\qquad
i_0^*k_\ell=\kappa\ast r,\qquad i_1^*k_\ell=\mathcal R(u\ast v).
\]

For the last endpoint: \(q\ast r=u\ast v\) and \(q_0=\sigma(\widehat F(u\ast v))\).

Let \(\gamma_6^{\rm leg}\) be the integral pair correction of full legal
triples (Section 7), with

\[
\delta_s\gamma_6^{\rm leg}(x,y)=T_3(x)+T_3(y)-T_3(x\ast y).
\]

Both cylinders are full legal triples on \(X\times I\), so
\(T_3(r_\ell)\) and \(T_3(k_\ell)\) are integral cocycles. Signed prism
Stokes gives

\[
\delta_sI\,T_3(r_\ell)=T_3(r)-T_3(z\ast z'),\qquad
\delta_sI\,T_3(k_\ell)=T_3(\mathcal R(u\ast v))-T_3(\kappa\ast r).
\]

Consequently

\[
G=\gamma_6^{\rm leg}(z,z')+\gamma_6^{\rm leg}(\kappa,r)-I\,T_3(r_\ell)-I\,T_3(k_\ell)
\]

satisfies exactly

\[
\delta_sG=T_3(z)+T_3(z')-T_3(\mathcal R(u\ast v))+T_3(\kappa).
\]

The functions `natural_product`, `natural_left_division` and
`legal_cylinder_transports` of
[../python/stacking_model/g6_general_section.py](../python/stacking_model/g6_general_section.py)
evaluate \(\ast\), \(\backslash\) and the two prisms
\(I\,T_3(r_\ell)\) and \(I\,T_3(k_\ell)\), after checking that \(T_3\) is
integral on each cylinder.

## 5. The correction \(g_6\) and the products \(\gamma_6\) and \(\gamma_7\)

For an image \(w=(a,b,c)\) and a natural lower triple \(u\), put

\[
j(w)=\begin{cases}J_6(\sigma w),&a=b=0,\\0,&\text{otherwise},\end{cases}
\qquad
\mathcal Q(u)=\begin{cases}\gamma_6^{\rm cl}\bigl(\sigma(\widehat Fu),\mathcal R(u)\bigr),&u\text{ has legal AB},\\0,&\text{otherwise},\end{cases}
\]

where \(\gamma_6^{\rm cl}\) is the closed-AB pair correction of Section 7.
On inputs with legal AB and arbitrary \(C\) it satisfies

\[
\delta_s\gamma_6^{\rm cl}(x,y)=J_6(x)+J_6(y)-J_6(x\ast y)+\gamma_7^{\rm pure}(t_x,t_y),
\]

with \(t_x\) the C residual of \(x\). Apply this to
\(u=\sigma(\widehat Fu)\ast\mathcal R(u)\), whose second factor has
\(t=0\). This gives

\[
g(u)=j(\widehat Fu)+T_3(\mathcal Ru)-\delta_s\mathcal Q(u).
\]

Thus \(g\) equals \(J_6\) on legal AB and \(T_3(\mathcal Ru)\) elsewhere;
\(g\circ\phi\) is the \(g_6\) of (A13), since \(\phi\) is the identity on
legal AB and \(\mathcal R_6=\mathcal R\circ\phi\)
(`NaturalSectionStacking.correction`, `CurrentSectionStacking.correction`).
It satisfies the square-zero identity of
[all_cochain_differential.md](all_cochain_differential.md), Section 4.
When the inputs \(u,v\) do not both have legal AB, the product correction
and its successor on the images are

\[
\gamma_6(u,v)=G-\mathcal Q(u)-\mathcal Q(v)+\mathcal Q(u\ast v),\qquad
\gamma_7(w,w')=T_3(\kappa)+j(\widehat Fq)-j(w)-j(w'),
\]

with \(q=\sigma(w)\ast\sigma(w')\) and \(\kappa=\sigma(\widehat Fq)\backslash q\).
Here \(\widehat Fq\) is the degree-seven product of \(w\) and \(w'\).
Substituting \(\delta_sG\) and the formula for \(g\) gives the strict
identity

\[
\delta_s\gamma_6(u,v)=g(u)+g(v)-g(u\ast v)+\gamma_7(\widehat Fu,\widehat Fv).
\]

When both inputs have legal AB, \(\gamma_6=\gamma_6^{\rm cl}(u,v)\) and
\(\gamma_7=\gamma_7^{\rm pure}(c,c')\) on the pure images. The branch is
selected by the output pair \((w,w')\) itself, so \(\gamma_7\) is well
defined. Both branches obey the same identity and agree with the legal
formulas on full legal inputs.

When \(z=z'=0\), both inputs are section representatives. Then \(r=0\),
the first cylinder is zero and the second is the constant \(\pi^*\kappa\),
whose prism vanishes. The two normalized legal terms vanish as well, so
\(\gamma_6=-\mathcal Q(u)-\mathcal Q(v)+\mathcal Q(u\ast v)\). The class
`NaturalSectionStacking` uses this unit case directly.

## 6. Return to the current product

The natural product and the current product differ by an explicit C
shift. For current-coordinate triples \(u,v\),

\[
\phi(\phi u\ast\phi v)=u\cdot v+(0,0,\delta\lambda),\qquad
\lambda=\begin{cases}0,&u,v\text{ both have legal AB},\\
\lambda_{\rm loc}(A,B;A',B'),&\text{otherwise}.\end{cases}
\]

Off the both-legal branch, the three values of \(e\) cancel between the
two sides and \(\widehat\beta-\beta_{\rm raw}=\delta\lambda_{\rm loc}\). On the
both-legal branch, \(e=0\) and \(\widehat\beta=\beta^{\rm leg}=\beta_6\).

For a triple \(y=(A,B,C)\) with legal AB and arbitrary \(C\), and a binary
\(\lambda\in C^4\), define the integral cochain

\[
L(y,\lambda)=\Omega_6(A,B,C+\delta\lambda)-\Omega_6(A,B,C)
-\delta_s\mathsf h\bigl[E_4(\lambda)+E_\times(C,\delta\lambda)\bigr],
\qquad
E_\times(C,x)=C\smile_4x+\delta C\smile_5x.
\]

For closed \(x\), \(\delta E_\times(C,x)=E(C+x)+E(C)+E(x)\), and
\(E_5(\delta\lambda)=\delta E_4(\lambda)\) by (A3). This makes \(L\) integral.
Since \(\widehat R_3(A,B)\) and the prime-three term of \(\Omega_6\) do not
depend on \(C\), the first difference equals
\(\mathsf h(E_5(C+\delta\lambda))-\mathsf h(E_5C)\), which is the form the
code evaluates (`C_shift_transport`). The residual \(t_6\) is unchanged by
\(C\mapsto C+\delta\lambda\), so

\[
\delta_sL(y,\lambda)=J_6(A,B,C+\delta\lambda)-J_6(A,B,C).
\]

On a full legal triple this is the difference of the \(T_3\) values. The
current correction is \(g_6(u)=g(\phi u)\). For a total outside \(L_6\),
take \(y=\mathcal R_6(u\cdot v)\): a C coboundary changes neither
\(\widehat F\) nor its section, and \(e\) depends only on \((A,B)\), so
\(\mathcal R_6(u\cdot v+\delta\lambda)=\mathcal R_6(u\cdot v)+(0,0,\delta\lambda)\).
For a total with legal AB, take \(y=u\cdot v\). In both cases

\[
g_6(u\cdot v+\delta\lambda)-g_6(u\cdot v)=\delta_sL(y,\lambda),
\]

and therefore

\[
\gamma_6^{\rm current}(u,v)=\gamma_6^{\rm natural}(\phi u,\phi v)+L(y,\lambda).
\]

The plus sign follows from
\(\delta_s\gamma_6^{\rm natural}=g_6(u)+g_6(v)-g_6(u\cdot v+\delta\lambda)+\gamma_7\).
Adding \(\delta_sL\) replaces the natural total by the current total.
`CurrentSectionStacking` implements this step. It checks the displayed
comparison of the two lower products before adding \(L\), and it adds
nothing when \(\lambda=0\). Its successor \(\gamma_7\) is that of
`NaturalSectionStacking`.

## 7. The production operations and the common pure-C normalization

`ProductionG6Section` in
[../python/stacking_model/production_g6_section.py](../python/stacking_model/production_g6_section.py)
takes from its caller the integral cycle coordinates, the binary
differentials on \(C^4\) and \(C^5\), the twists \(s\) and \(\omega\), and
the conversions between coordinate vectors and cochains, which must
implement the normalized simplicial pullbacks. These linear finite-basis
maps are its only caller-supplied callbacks. It selects every nonlinear
operation itself:

- \(P_6\), \(f_6\), \(e\) and \(\alpha\) from `off_shell_beta`;
- \(\widehat\beta\) and \(\lambda_{\rm loc}\) from `lower_raw_comparison`;
- the current \(\beta_6\) from `stacking_lower.all_cochain_beta`;
- \(T_3\) and \(J_6\);
- the two legal-cylinder prisms and the C transport \(L\) of
  `g6_general_section`.

The pair corrections carry the common pure-C normalization and the
rephasing by the global closed-pure-C selector. For a function \(M\) of
lower triples, write \(\Delta M(u,v)=M(u)+M(v)-M(u\cdot v)\), and put

\[
K_6(A,B,C)=\begin{cases}\delta_s\mathsf h(s\smile C),&A=0,\ B=0,\ \delta C=0\text{ on }X,\\
0,&\text{otherwise},\end{cases}
\qquad
K_7(c)=\delta_s\mathsf h(s\smile c).
\]

The three corrections are

\[
\begin{aligned}
\gamma_6^{\rm leg}&=\gamma_6^{\rm norm}-\Delta K_6,\\
\gamma_6^{\rm cl}&=\gamma_6^{\rm cl,norm}-\Delta\,\mathsf h(s\smile t)-\Delta K_6,\\
\gamma_7^{\rm pure}(c,c')&=\gamma_7^{\rm can}(c,c')-\Delta K_7(c,c').
\end{aligned}
\]

Here \(\gamma_6^{\rm norm}\) is `pure_c_normalization.gamma`, the
production correction of `production_gamma6` plus the pure-C
normalization, an integral coboundary.
\(\gamma_6^{\rm cl,norm}\) is `pure_c_normalization.closed_gamma`, the
closed-AB correction `closed_ab_upper.gamma` with the same normalization.
\(\gamma_7^{\rm can}\) is `compatible_sector.pure_c_gamma`, and \(t\) is the
C residual of each argument. On legal AB,
\(t_{u\cdot v}=t+t'\), so \(\Delta\,\mathsf h(s\smile t)\) is an integral
carry, and its coboundary is exactly \(\Delta K_7(t,t')\). Subtracting the
carry in \(\gamma_6^{\rm cl}\) and \(\Delta K_7\) in \(\gamma_7^{\rm pure}\)
therefore gives the common successor. The term \(-\Delta K_6\) is integral
and closed and does not change any identity. In the section comparison
\(\mathcal Q(u)\), the second argument \(\mathcal R(u)\) has \(t=0\), so the
half-lift carry vanishes identically and Section 5 is unchanged. The
differential uses \(J_6\) (`closed_ab_upper.J`) on all of \(L_6\), without
the \(K\) rephasing.

## 8. Arbitrary X

The algebraic construction needs only a fixed pointed set section
\(\sigma:\operatorname{im}\widehat F\to C^3(X;\mathbf Z_s)\times C^4(X;\mathbf F_2)\times C^5(X;\mathbf F_2)\)
with \(\sigma(0)=0\). It does not need an integral complement or a finite
cochain basis. For an arbitrary cochain model of \(X\), such a section can
be fixed by set-theoretic choice, and the two cylinder expressions prove
the same identities. The cylinders pull back and ramp cochains already
chosen on \(X\). They never evaluate a section on \(X\times I\), and no
naturality of \(\sigma\) is assumed.

The finite algorithm of Section 2 makes this choice computable by exact
finite enumeration. It is an implementation of the general construction,
not a hypothesis of its proof. The resulting extension depends on the
global section and is not a natural local operation in \(X\). It shares
this global character with the branch selection of
[all_cochain_differential.md](all_cochain_differential.md).

## 9. A worked example on the \(C_2\) bar

Take the normalized bar complex of \(C_2\), with one cell in each degree,
\(s=0\) and \(\omega=1\). The lower triples \(u=(1,0,0)\) and
\(v=(1,1,0)\) are their own section representatives, so \(z=z'=0\). Their
images are \(F_6u=(2,0,1)\) and \(F_6v=(2,0,0)\), and their product is
\(u\cdot v=(2,0,1)\), with \(F_6(u\cdot v)=(4,0,0)\). Here \(r=0\), the first
cylinder is identically zero and the second is the constant \(\kappa\).
Both prism terms and both normalized legal terms vanish, and
\(\mathcal Q=0\) because no image is pure. The result is \(\gamma_6=0\) and

\[
\delta_s\gamma_6+g_6(u\cdot v)=1=g_6(u)+g_6(v)+\gamma_7(F_6u,F_6v),
\]

with \(g_6(u)=g_6(v)=0\), \(\mathcal R_6(u\cdot v)=\kappa\) and
\(\gamma_7=T_3(\kappa)=1\). The first two values hold because \(u\) and
\(v\) are \(\sigma_0\)-representatives, so \(\mathcal R_6u=\mathcal R_6v=0\).
For the third, binary coboundaries vanish on the one-cell \(C_2\) bar, so
\(\delta\lambda=0\), \(\phi(u\cdot v)=\phi u\ast\phi v=q\) and
\(\mathcal R_6(u\cdot v)=q_0\backslash q=\kappa\). In particular the
example has nonzero \(A\) outside \(L_6\) and a nonzero production
obstruction. The test `test_nonzero_off_shell_A_and_successor` of
[../python/test_extension_degree_six.py](../python/test_extension_degree_six.py)
takes the states \(x=(u,0)\) and \(y=(v,0)\), whose product is
\(x\times y=(u\cdot v,0)\) because \(\gamma_6=0\), and checks
\(\mathfrak d_6(x\times y)=\mathfrak d_6x\times\mathfrak d_6y=(4,0,0,1)\)
and \(\mathfrak d_7\mathfrak d_6x=0\).

The selector \(K_6=\delta_s\mathsf h(s\smile C)\) vanishes identically
when \(s=0\), so the example above does not exercise it. With
\(s=\omega=1\) on the same bar it is nonzero: for the pure closed \(C=1\)
in degree five, the interior faces of \(\delta_s\) are degenerate in the
normalized bar; the first face carries the twist sign \((-1)^{s}=-1\) and
the last face the sign \((-1)^7\), so
\(K_6(0,0,1)=-\tfrac12-\tfrac12=-1\). Any nonzero \(A\) or \(B\) turns the
selector off. On this bar with \(s=1\), \(K_6\) therefore enters
\(\gamma_6^{\rm leg}\) and \(\gamma_6^{\rm cl}\) through \(-\Delta K_6\).

## 10. The finite implementation and its use in the package

The classes of
[../python/stacking_model/g6_general_section.py](../python/stacking_model/g6_general_section.py)
evaluate the formulas above in finite coordinates. An integral vector
represents \(A\), and bit strings represent \(B\) and \(C\).

| Class or function | Evaluates |
| --- | --- |
| `FiniteCurrentSection` | \(F_6\) (`lower`), \(\widehat F\) (`natural_lower`), \(\phi\) (`conjugate`), \(\sigma_0\) (`current_section`), \(\sigma\) (`section`) |
| `natural_product`, `natural_left_division` | \(\ast\) and \(\backslash\) on cochains |
| `legal_cylinder_transports` | \(I\,T_3(r_\ell)\) and \(I\,T_3(k_\ell)\) |
| `C_shift_transport` | \(L(y,\lambda)\) |
| `NaturalSectionStacking` | \(\ast\), \(\backslash\), \(\mathcal R\), \(j\), \(\mathcal Q\), \(g\), \(\gamma_6\) and \(\gamma_7\) in natural coordinates |
| `CurrentSectionStacking` | \(u\cdot v\), \(g_6=g\circ\phi\), \(\gamma_6^{\rm current}\) and \(\gamma_7\) |

Each step checks its defining equation: the section's lower equations,
\(\widehat F\mathcal R=0\), the legality of \(\kappa\) and \(r\), the
compatibility of the lower product with \(\widehat F\), the integrality of
\(T_3\) on each cylinder and of every transport, and the comparison of
Section 6.

`extension_degree_six.configure_degree_six` in
[../python/extension_degree_six.py](../python/extension_degree_six.py)
binds `ProductionG6Section` to the normalized finite bar of a group. The
GAP side, `KOAHSS_ExtensionDegreeSixData` in
[../gap/extension_degree_six.gi](../gap/extension_degree_six.gi),
supplies the unimodular cycle basis from a Smith normal form, and the
adapter checks every differential against the bar basis. The section
search is bounded: a request whose exponent \(2r+\dim\ker(\delta:C^4\to C^5)\)
reaches the bit length of `maxSectionCandidates` (default 4096) is refused
with `DegreeSixResourceLimit`. The zero image returns the zero triple
without a search. `four_cochain_stacking.Stacking(degree_six=...)` uses the
adapter for \(\mathfrak d_6\), for the degree-six product, and for the
degree-seven endpoint product on legal images. The worker
[../python/extension_worker.py](../python/extension_worker.py) attaches
it for the complete-bar reference model of
[../gap/extension_bar.gi](../gap/extension_bar.gi) and
[../gap/extension_degree_six.gi](../gap/extension_degree_six.gi). Only
the tests use that model: `tst/extension_degree_six.tst` (opt-in) compares
it with the native model on states with A zero, and
[../python/test_extension_degree_six.py](../python/test_extension_degree_six.py)
covers the example of Section 9, the pointed zero, the search bound and
the rejection of mismatched coordinates.

The native `koFull` engine,
[../python/extension_transfer.py](../python/extension_transfer.py) with
[../python/extension_native_upper.py](../python/extension_native_upper.py),
does not evaluate this construction. In degree six it evaluates \(J_6\)
and the legal \(\gamma_6\) on the legal lower locus only. That product
uses `closed_ab_upper.gamma` with the carry and \(K\) terms of Section 7
and without the pure-C normalization, which is an integral coboundary, so
native D representatives differ from the reference model's by
D-coboundaries. In degrees five and six the A=0 sector uses the direct
formulas: the curvature of a lower-legal triple with A=0 is evaluated by
`a0_high_gamma` (and by `compatible_sector.pure_c_g` when B=0 as well),
which equals \(J_k\) (\(J_5\), respectively \(J_6\)) as a cochain. For two
fully legal triples with A=B=0, `compatible_sector.pure_c_gamma` replaces
`closed_ab_upper.gamma` as the base term (they differ by an integral
coboundary), and the half-lift carry and the \(K\) change of Section 7 are
kept. The degree-six product
of two states whose A layers are both nonzero is refused
(`PairSourceLimit`), and the degree is left unresolved, unless
`FERMIONAHSS_DEGREE_SIX_A_STACKING=1` ([extensions.md](extensions.md)).
Native curvature stops at the first nonzero layer, and products and gauge
actions receive flat states and differential images, whose lower pairs are
legal. The native operations therefore request \(g_6\) and \(\gamma_6\)
only on \(L_6\). A request outside \(L_6\) would need \(\mathcal R_6\),
which is not evaluated on the resolution; the engine raises
`SectionBranchRequired` and the degree is left unresolved
([extensions.md](extensions.md)).

The stacking research workspace (not bundled) has further checks of the
construction:

- the generic-twist legal test of the \(\lambda_3\) constructor with
  nonzero mismatch;
- the section test of signed integral image membership;
- the modulo-four enumeration, with a comparison gauge that is not
  four-periodic;
- a synthetic triangular loop with a nonzero associator, on which the
  two-cylinder algebra holds for all sixteen selected pairs;
- a synthetic pair with a nonzero local pair gauge, on which the natural
  \(\gamma_6\) is zero, the current \(\gamma_6\) is one, and the identity
  confirms the plus sign of \(L\);
- the \(C_2\) bar with \(s=\omega=1\), on which the rephased legal and
  closed-AB pair corrections both equal two.

The synthetic checks isolate the section identity. They do not test the
production formulas.
