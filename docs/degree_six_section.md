# Degree six outside the legal locus: the section construction

This page specifies the degree-six branch for inputs whose lower pair is
**not** legal: the D-layer differential $g_6=T_3\circ\mathcal R_6$ and the
product correction $\gamma_6$ with its successor $\gamma_7$ on
differential images. Together they satisfy the strict upper identity

$$
\delta_s\gamma_6(u,v)=g_6(u)+g_6(v)-g_6(u\cdot v)+\gamma_7(F_6u,F_6v)
$$

without assuming associativity of any lower product and without solving
an upper equation. Notation: [conventions.md](conventions.md); here
$A\in C^3(\mathbb Z_s)$, $B\in C^4$, $C\in C^5$, $D\in C^7$, $g_6\in C^8$,
$\gamma_6\in C^7$, $\gamma_7\in C^8$. (`g6_general_section`,
`production_g6_section`.)

## 1. Data

- $F_6(A,B,C)=(\delta_sA,\ \delta B+P(A),\ \delta C+f_6(A,B))$ ([differential.md](differential.md), degree 6). A triple has **legal AB** when $(A,B)\in L_6$; it is a **full legal triple** when $F_6u=0$; an image $w=(a,b,c)$ is **pure** when $a=b=0$.
- $q_{\rm loc}$, $\widehat f=H_6+\delta q_{\rm loc}$, $e$, $\phi(A,B,C)=(A,B,C+e(A,B))$, $\widehat F=F_6\circ\phi$ ([product_lower.md](product_lower.md) Section 5).
- $u\cdot v$: the current lower product with $\alpha_6,\beta_6$; $u\ast v$: the natural product with $\widehat\beta_6=\beta^{\rm raw}_6+\delta\lambda_{\rm loc}$ (Section 7 there). Left division $v\backslash w=z$, the unique solution of $v\ast z=w$:
  $$A_z=A_w-A_v,\qquad B_z=B_w+B_v+\alpha_6(A_v,A_z),\qquad C_z=C_w+C_v+\widehat\beta_6(A_v,B_v;A_z,B_z).$$
- $T_3=\delta_s\Omega_6$ on full legal triples (with the prime-three term); $J_6$ of (J) on legal AB with arbitrary $C$.
- The pair corrections of [product_upper_k4_k6.md](product_upper_k4_k6.md) Section 5: $\Gamma^L_6$ (full legal triples), $\Gamma^{\rm sec}_6:=\gamma_6$ of that section (legal AB, arbitrary $C$), and the rephased pure-C successor $\gamma_7^{\rm pure}(c,c')=\gamma_C(c,c')-K_7(c)-K_7(c')+K_7(c+c')$.

## 2. A pointed section

Choose $\sigma$ with $\sigma(0)=0$ and $\widehat F\sigma=\mathrm{id}$ on
$\operatorname{im}\widehat F$. For a finite free cochain model:

1. Fix a unimodular integral basis of $C^3(X;\mathbb Z_s)$ whose first $r$ columns $Z$ span $\ker\delta_s$; for an output $(a,b,c)$ solve $\delta_sA_0=a$ in the complement (exact elimination) and enumerate $A=A_0+Zt$, $t\in\{0,1,2,3\}^r$ (lexicographic).
2. For each $A$ solve $\delta B=b+P(A)$ and enumerate $B$ over the coset of $\ker(\delta:C^4\to C^5)$ (binary counting in a fixed kernel basis; particular solutions in reduced row echelon form with free variables zero).
3. For each $B$ solve $\delta C=c+f_6(A,B)$.

The first solution is $\sigma_0(a,b,c)$ (the enumeration modulo four is
complete because $\tau'_3$ and $H_6$ are four-periodic in $A$ and the
branch test is unchanged by adding $4\times$ an integral cycle); then
$\sigma=\phi\circ\sigma_0$. The **retraction**

$$
\boxed{\;\mathcal R(u)=\sigma(\widehat Fu)\backslash u\ \text{(natural coordinates)},\qquad \mathcal R_6(u)=\mathcal R(\phi u)\ \text{(current coordinates)}\;}
$$

always returns a full legal triple ($\widehat F$ is a homomorphism for
$\ast$, and left cancellation gives $\widehat F\mathcal R(u)=0$) and fixes
every full legal triple. Hence $g_6=T_3(\mathcal R_6u)$ outside $L_6$ is an
integral cocycle.

## 3. Two legal cylinders

For natural triples $u,v$ put

$$
w=\widehat Fu,\quad w'=\widehat Fv,\quad a=\sigma(w),\quad b=\sigma(w'),\quad z=\mathcal R(u),\quad z'=\mathcal R(v),\quad
q=a\ast b,\quad q_0=\sigma(\widehat Fq),\quad \kappa=q_0\backslash q,\quad r=q\backslash(u\ast v);
$$

$\kappa$ and $r$ are full legal triples. On $X\times I$ form

$$
r_\ell=(a\ell\ast b\ell)\backslash\bigl[(a\ell\ast\pi^*z)\ast(b\ell\ast\pi^*z')\bigr],\qquad
k_\ell=(q_0\ell)\backslash\bigl[\bigl((q_0\ell)\ast\pi^*\kappa\bigr)\ast\pi^*r\bigr],
$$

full legal triples on the cylinder with ends
$i_0^*r_\ell=z\ast z'$, $i_1^*r_\ell=r$, $i_0^*k_\ell=\kappa\ast r$, $i_1^*k_\ell=\mathcal R(u\ast v)$.
By prism Stokes

$$
\mathcal G=\Gamma^L_6(z,z')+\Gamma^L_6(\kappa,r)-I\,T_3(r_\ell)-I\,T_3(k_\ell)
\qquad\text{satisfies}\qquad
\delta_s\mathcal G=T_3(z)+T_3(z')-T_3\bigl(\mathcal R(u\ast v)\bigr)+T_3(\kappa).
$$

## 4. $g$, $\gamma_6$ and $\gamma_7$ in natural coordinates

$$
j(w)=\begin{cases}J_6(\sigma w),&w\text{ pure},\\0,&\text{otherwise},\end{cases}\qquad
\mathcal Q(u)=\begin{cases}\Gamma^{\rm sec}_6\bigl(\sigma(\widehat Fu),\mathcal R(u)\bigr),&u\text{ has legal AB},\\0,&\text{otherwise},\end{cases}
$$

$$
g(u)=j(\widehat Fu)+T_3(\mathcal Ru)-\delta_s\mathcal Q(u)\quad(=J_6\text{ on legal AB},\ T_3(\mathcal Ru)\text{ elsewhere}),\qquad g_6=g\circ\phi .
$$

When $u,v$ do **not** both have legal AB:

$$
\boxed{\;
\gamma_6(u,v)=\mathcal G-\mathcal Q(u)-\mathcal Q(v)+\mathcal Q(u\ast v),\qquad
\gamma_7(w,w')=T_3(\kappa)+j(\widehat Fq)-j(w)-j(w'),\;}
$$

with $q=\sigma(w)\ast\sigma(w')$ and $\kappa=\sigma(\widehat Fq)\backslash q$
($\widehat Fq$ is the degree-seven product of $w,w'$). When both have
legal AB, $\gamma_6=\Gamma^{\rm sec}_6(u,v)$ and $\gamma_7=\gamma_7^{\rm pure}(c,c')$
on the pure images. The branch is selected by $(w,w')$, so $\gamma_7$ is
well defined; its values outside $\operatorname{im}F_6$ may be set to
zero. When both inputs are section representatives ($z=z'=0$): $r=0$,
the first cylinder vanishes, the second is the constant $\kappa$, and
$\gamma_6=-\mathcal Q(u)-\mathcal Q(v)+\mathcal Q(u\ast v)$.

## 5. Return to the current product

$\phi(\phi u\ast\phi v)=u\cdot v+(0,0,\delta\lambda)$ with $\lambda=0$ when
both inputs have legal AB and $\lambda=\lambda_{\rm loc}(A,B;A',B')$
otherwise. For a triple $y$ with legal AB and a binary $\lambda\in C^4$,
the integral transport (`C_shift_transport`)

$$
L(y,\lambda)=\mathsf h\bigl(E(C_y+\delta\lambda)\bigr)-\mathsf h\bigl(E(C_y)\bigr)
-\delta_s\mathsf h\bigl(E(\lambda)+h^2(C_y,\delta\lambda)\bigr),\qquad
\delta_sL(y,\lambda)=J_6(y+(0,0,\delta\lambda))-J_6(y),
$$

and with $y=\mathcal R_6(u\cdot v)$ when the total is outside $L_6$, $y=u\cdot v$ otherwise,

$$
\boxed{\;\gamma_6^{\rm current}(u,v)=\gamma_6^{\rm natural}(\phi u,\phi v)+L(y,\lambda).\;}
$$

## 6. Scope

Only a fixed pointed set section with $\sigma(0)=0$ is needed; the finite
enumeration makes it computable (worst case
$4^r\cdot2^{\dim\ker(\delta:C^4\to C^5)}$ candidates). The result depends
on the global section and is not natural under pullback. The previous
implementation evaluated this branch only in its complete-bar reference
model, never in the production engine, which only ever meets legal lower
pairs. With $g_6=0$ outside $L_6$ no normalized product could satisfy the
upper identity: on $B(\mathbb Z/3)^2$ the full defining system
$(4\beta_{\mathbb Z,3}(u_1u_2),0,0,0)$ has nonzero $T_3$ class
$2(U_1^3U_2-U_1U_2^3)$, and stacking it with an input outside $L_6$ would
force that class to be a coboundary.
