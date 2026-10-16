# The D-layer correction $\gamma_k$ in degrees 4, 5 and 6

This page writes out the pair phases of degrees four to six, which the
notes of the previous implementation left to code. Every formula names
the function it was transcribed from (`python/stacking_model`, or
`python/` for `phase_eval`, `closed_a_upper`). Notation:
[conventions.md](conventions.md); $n=k-3\in\{1,2,3\}$, $|A|=n$,
$|B|=n+1$, $|C|=n+2$, $|P(A)|=n+2$, $|k_0(A)|=n+3$; lower products with
$\alpha_k$, $\beta^L_k$, $\beta^\sharp$: [product_lower.md](product_lower.md).
Universal primitives ($\mathcal U^\Theta_m$, $\mathcal P_k$, $\mathfrak q$,
$V^{\rm PS}$): [universal_primitives.md](universal_primitives.md).

The result has the same structure in the three degrees:

$$
\boxed{\;\Gamma^{\rm prod}_k(u,v)=\Omega_k(u)+\Omega_k(v)-\Omega_k(\mu^\sharp(u,v))+\delta_se^{\rm prod}_k(u,v)\quad\text{(integral on full defining systems)},\;}
$$

followed by a normalization and a rephasing (Section 5) that give the
$\gamma_k$ actually used. Sections 1–4 define $e^{\rm prod}_k$.

## 1. Common ingredients

**Secondary data** ((S) and (S23) of [secondary_operations.md](secondary_operations.md), Section 6): $P=P(A)$, $t=t(A)$,
$g=g(A)$, $k_0=k_0(A)$, $\lambda(A,B)=g\cdot(s\smile B)$. For $n=1$ the
degree-four phases use the *adjusted* source (`phase_eval.source`, $n=1$)

$$
e_0=s\smile t\ (\deg 2),\qquad q=(\omega+s\smile s)\smile a\ (\deg 3),\qquad
k=k_0+Q_D(e_0)+\delta e_0\smile_2P+s\smile(\delta e_0\smile_3P)\ (\deg 4);
$$

for $n=2,3$ put $q=P$, $k=k_0$.

**Extension of $\tau'$ to arbitrary $B$** (`closed_a_upper.fsharp`), for a signed cocycle $A$:

$$
f^\sharp(A,B)=Q_D(B)+k_0(A)+\delta\lambda(A,B)+h^D(\delta B,P(A))+\delta\bigl(Q^1B+s\smile B\bigr),
$$

equal to $\tau'_n(A;B)$ when $\delta B=P(A)$. The **curvature** of a
triple with closed $A$ is $(b,c)=\bigl(\delta B+P(A),\ \delta C+f^\sharp(A,B)\bigr)$
(`closed_a_upper.curvature`); the lower product with $\alpha_k,\beta^\sharp$
is written $\mu^\sharp$.

**Binary helpers** (`production_upper_binary_comparison`). For a binary $B$ of degree $n+1$,

$$
J_0(B)=Q^1(Q_DB)+s\smile Q_DB+E(s\smile B)\qquad(\deg n+4),
$$

and for binary $x,y$ of equal degree $r$, with $h=h^D(x,y)$, $t_h=h^D(\delta x,\delta y)$,
$L_x=Q_Dx$, $L_y=Q_Dy$, $z=L_x+L_y$, $w=\delta h+t_h$ and $K(p,q)=p\smile_{|p|}q$,

$$
j_0(x,y)=K(L_x,L_y)+Q^1(h)+K(z,w)+K(\delta h,t_h)+s\smile h+h^2(s\smile x,s\smile y),
\qquad \delta j_0(x,y)+j_0(\delta x,\delta y)=J_0(x+y)+J_0(x)+J_0(y).
$$

With $\alpha=\alpha_k(A,A')$: the A-only residue
$\mathfrak J(A,A')=J_0(\alpha)+j_0(P,P')+j_0(P+P',\delta\alpha)$ and the B
primitive $\mathfrak j(A,B;A',B')=j_0(B,B')+j_0(B+B',\alpha)$.

**C-transport** (`production_gamma6_comparison.C_transport`): for $|L|=|c|-1$,
$\mathcal T(c,L)=\mathsf h\bigl(E(L)+h^2(c,\delta L)\bigr)$.

**The $\Theta$-pair phase** (`theta_pair_phase.ThetaPairPhase(m)`), for
closed binary $q,q'$ of degree $m\in\{3,4,5\}$ and $v,v'$ of degree $m+1$:

$$
\boxed{\;
\Theta^{\rm pair}_m(q,v;q',v')=\mathsf h\Bigl(h^2(v,v')+h^2\bigl(v+v',h^D(q,q')\bigr)\Bigr)+\mathcal U^\Theta_m(q,q'),\;}
$$

where $\mathcal U^\Theta_m(q,q')$ (degree $m+2$) is the universal primitive of

$$
S^\Theta_m(q,q')=\theta_m(q+q',0)-\theta_m(q,0)-\theta_m(q',0)
+\mathsf h\Bigl(h^2(Q_Dq,Q_Dq')+E\bigl(h^D(q,q')\bigr)+h^2\bigl(Q_Dq+Q_Dq',\,\delta h^D(q,q')\bigr)\Bigr)
$$

(Section 7 of the universal-primitives page; $\theta_m(q,0)$ is the
boundary phase with $v=0$).

**The A-only source primitive** $V_P(A)$ (`closed_a_upper.source_primitive`; the tertiary $V_n$ of [tertiary_operations.md](tertiary_operations.md), Section 8):
the tertiary $V_n(A)$ for $n=2,3$; for $n=1$ the tertiary $V_1(A)$ minus the prism

$$
\Pi^{\rm src}(A)=I\,\theta_3\Bigl(\pi^*P+\delta e_I,\ \pi^*k_0+Q_D(e_I)+h^D(\delta e_I,\pi^*P)\Bigr),\qquad e_I=e_0\ell,
$$

so that in all three cases $\delta_sV_P(A)\equiv\theta_{n+2}(P,k_0)$.

**The auxiliary potential** (`closed_a_upper.phi`, `upper_pair_source.phi_without_A_source`):

$$
\varphi(A,B,C)=\mathsf h\Bigl(E(C)+h^2\bigl(\delta C,f^\sharp(A,B)\bigr)\Bigr)-V_P(A)-I\,\theta_{n+2}\bigl(b_I,\ f^\sharp(\pi^*A,B\ell)\bigr),
\qquad b_I=\delta(B\ell)+P(\pi^*A),\qquad \varphi_0=\varphi+V_P(A),
$$

with $\delta_s\varphi\equiv-\theta_{n+2}(b,c)$ on the curvature.

**The pair source and its ramp** (`upper_pair_source.source`, `reduced_full_source`):

$$
S^{\rm aux}(A,A')=-V_P(A)-V_P(A')-\varphi\bigl(A+A',\ \alpha,\ \beta^\sharp(A,0;A',0)\bigr)-\Theta^{\rm pair}_{n+2}(P,k_0;P',k_0'),
$$

$$
\mathrm{RFS}(u,v)=\varphi_0(u)+\varphi_0(v)-\varphi_0\bigl(\mu^\sharp(u,v)\bigr)-\Theta^{\rm pair}_{n+2}(b_u,c_u;b_v,c_v),
$$

the latter always evaluated on the cylinder data $(\pi^*A,B\ell,C\ell;\ \pi^*A',B'\ell,C'\ell)$
with pulled-back twists (there $V_P(\pi^*A)$ has zero prism, which is why
$\varphi_0$ replaces $\varphi$).

**The comparison with the production phase.** In each degree there are an
A-only phase $R^{(k)}_A(A)$ and a gauge $Q_k$ (depending on $A,B$, and on
$C$ for $k=4$) such that

$$
\Omega_k(A,B,C)\equiv\varphi(A,B,C)+R^{(k)}_A(A)+\mathsf h\bigl(J_0(B)\bigr)+\delta_sQ_k\pmod{\mathbb Z}.
$$

The common part of $R_A$ (`production_gamma6_comparison.reduced_A_phase`), with $r=Q^1(q)+s\smile q$:

$$
R^{\rm red}(A)=\mathcal U^\Theta_{n+2}(q,q)+\mathsf h\Bigl(s\smile k+Q^1(r)+\delta k+E(q)+h^2\bigl(Q^1q,\,s\smile q\bigr)\Bigr)
$$

(the binary sum inside one lift; $\delta k$ reduced modulo two). The
common part of $Q_k$ (`production_gamma6_comparison.comparison_gauge`),
with $m=n+2$ and, on $X\times I$, $u=\delta(B\ell)$, $y=Q_D(B\ell)$,
$c_0=\pi^*k+y+h^D(\pi^*q,u)$, $L=\lambda(\pi^*A,B\ell)+Q^1(B\ell)+s\smile B\ell+\pi^*q\smile_mu$:

$$
Q^{\rm gen}(A,B)=\mathsf h\bigl(\mathcal Q(B,k,\lambda(A,B))\bigr)-I\Bigl[\Theta^{\rm pair}_m(\pi^*q,\pi^*k;\,u,y)+\mathcal T(c_0,L)\Bigr],
$$

where, for a binary $B$ of degree $m-1$ and $\lambda$ of degree $m$, with
$P_B=\delta B$, $y=Q_DB$, $t_\lambda=\delta\lambda$, $r_B=Q^1(P_B)+s\smile P_B$,
$V_0=Q^1B+s\smile B$, $Z_0=V_0+P_B$ (`comparison_terms`):

$$
\mathcal Q(B,k,\lambda)=Q^1(\lambda)+s\smile\lambda+h^2(\lambda,Z_0)+t_\lambda\smile_{m+1}r_B+(y+k)\smile_{m+1}t_\lambda+h^2(V_0,P_B)+h^2\bigl(Q^1B,\,s\smile B\bigr).
$$

## 2. Degree 4 ($n=1$)

Because the degree-four tertiary phase uses the adjusted source $(q,k)$,
every object is first *shifted* by $e_0$ (`production_gamma4_comparison`).
Write $b^\circ=B+e_0$ and

$$
\lambda^{\rm sh}(A,B)=\lambda(A,B)+B\smile_1e_0+P\smile_2e_0+P\smile_3\delta e_0+s\smile\bigl(B\smile_2e_0+P\smile_3e_0+P\smile_4\delta e_0\bigr),
$$
$$
f^{\rm sh}(A,B)=Q_D(b^\circ)+k+\delta\lambda^{\rm sh}(A,B)+h^D(\delta b^\circ,q)+\delta\bigl(Q^1b^\circ+s\smile b^\circ\bigr).
$$

On $X\times I$ put $b^\circ_t=b^\circ\ell$, $B_t=b^\circ_t+\pi^*e_0$,
$u=\delta b^\circ_t$, $y=Q_D(b^\circ_t)$. Then

$$
\begin{aligned}
\Phi^{\rm sh}(A,B,C)&=\mathsf h\Bigl(E(C)+h^2\bigl(\delta C,f^{\rm sh}(A,B)\bigr)\Bigr)-\mathcal T\bigl(k,\lambda^{\rm sh}(A,e_0)\bigr)-I\,\theta_3\bigl(u+\pi^*q,\ f^{\rm sh}(\pi^*A,B_t)\bigr),\\
D^{\rm cl}(A,B,C)&=\Phi^{\rm sh}(A,B,C)+\mathcal T\bigl(c,\ b\smile_2e_0+b\smile_3\delta e_0\bigr)-\varphi_0(A,B,C)-\Pi^{\rm src}(A),\qquad (b,c)\text{ the curvature},\\
\Sigma^{\rm sh}(A)&=\mathcal U^\Theta_3(q,q)+\mathsf h\Bigl(s\smile k+Q^1(r)+\delta k+E(q)+h^2(Q^1q,s\smile q)\Bigr),\qquad r=Q^1(q)+s\smile q,\\
R^{(4)}_A(A)&=D^{\rm cl}(A,0,0)+\Sigma^{\rm sh}(A)+\mathsf h\Bigl(J_0(e_0)+j_0(P,\delta e_0)\Bigr)\quad(\text{the selector }\epsilon=0\text{ removes the term }-\tfrac\epsilon4P_\omega\smile_sA),\\
Q^{\rm sh}_4(A,B)&=\mathsf h\bigl(\mathcal Q(b^\circ,k,\lambda^{\rm sh}(A,B))\bigr)-I\Bigl[\Theta^{\rm pair}_3(\pi^*q,\pi^*k;\,u,y)+\mathcal T(c_0,L)\Bigr]+\mathsf h\bigl(Q^1(\lambda^{\rm sh}(A,B))\bigr),\\
&\qquad c_0=\pi^*k+y+h^D(\pi^*q,u),\quad L=\lambda^{\rm sh}(\pi^*A,B_t)+Q^1(b^\circ_t)+s\smile b^\circ_t+\pi^*q\smile_3u,\\
Q_4(A,B,C)&=Q^{\rm sh}_4(A,B)+I\,D^{\rm cl}(\pi^*A,B\ell,C\ell)+\mathsf h\bigl(j_0(B,e_0)\bigr).
\end{aligned}
$$

The pair source and its universal primitive (`production_gamma4.source`, `primitive`):

$$
S_4(A,A')=S^{\rm aux}(A,A')+R^{(4)}_A(A)+R^{(4)}_A(A')-R^{(4)}_A(A+A')+\mathsf h\bigl(\mathfrak J(A,A')\bigr),
$$

$\mathcal P_4(A,A')$ = the universal primitive of $S_4$ on the shared-sign
degree-one pair model with the two fixed coefficients $(\tfrac18,-\tfrac18)$
(universal-primitives page, Section 8). Finally (`production_gamma4.phase4`)

$$
\boxed{\;
e^{\rm prod}_4(u,v)=-\mathcal P_4(A,A')-I\,\mathrm{RFS}(u,v)-\mathsf h\bigl(\mathfrak j(A,B;A',B')\bigr)-\bigl[Q_4(u)+Q_4(v)-Q_4(\mu^\sharp(u,v))\bigr].\;}
$$

## 3. Degree 5 ($n=2$)

Here $q=P$, $k=k_0$, $m=4$, and (`production_gamma5_comparison`, with the selector $\alpha=1$)

$$
\begin{aligned}
R^{(5)}_A(A)&=R^{\rm red}(A)+\mathsf h\Bigl(\zeta_{1,4}(s,P)+(\omega\smile_1s)\smile P\Bigr)-\tfrac14\,\widetilde{Q_D(P)}
+\tfrac14\,P_\omega\smile_sA-\tfrac14\,A^{\smile3}+\tfrac23\,\mathrm L\bigl(\rho_3A^{\smile3}\bigr),\\
Q_5(A,B)&=Q^{\rm gen}(A,B)+\mathsf h\bigl(s\smile\lambda(A,B)\bigr)-I\,\mathsf h\bigl(E(Q^1(B\ell))\bigr)+\tfrac14\,\widetilde{Q_D(B)},\\
\mathfrak B(A,A')&=\mathsf h\Bigl(s\smile\bigl(Q_D(\alpha)+h^D(P,P')+h^D(P+P',\delta\alpha)\bigr)\Bigr),\qquad
\mathfrak B'(A,B;A',B')=\mathsf h\Bigl(s\smile\bigl(h^D(B,B')+h^D(B+B',\alpha)\bigr)\Bigr),\\
\mathfrak T^{\rm s}(A,A')&=\mathsf h\bigl(a\smile P'+a'\smile P\bigr),\qquad
\mathfrak T(A,B;A',B')=\mathsf h\bigl(a\smile B'+a'\smile B\bigr).
\end{aligned}
$$

(The term $\mathsf h(s\smile k_0)$ of the code cancels against the same
term inside $U_A$ and is omitted above; $\alpha=\alpha_5(A,A')$.) Then
(`production_gamma5.source`, `phase`)

$$
S_5(A,A')=S^{\rm aux}(A,A')+R^{(5)}_A(A)+R^{(5)}_A(A')-R^{(5)}_A(A+A')+\mathfrak B(A,A')-\mathfrak T^{\rm s}(A,A'),
$$

$\mathcal P_5(A,A')$ = the universal primitive of $S_5$ on the shared-sign
degree-two pair model with the 48 fixed coefficients of Section 8 there, and

$$
\boxed{\;
e^{\rm prod}_5(u,v)=-\mathcal P_5(A,A')-I\,\mathrm{RFS}(u,v)-\mathfrak B'(A,B;A',B')-\bigl[Q_5(A,B)+Q_5(A',B')-Q_5(A_S,B_S)\bigr]-\mathfrak T(A,B;A',B').\;}
$$

The lower product includes the closed term $a\smile a'$ of $\beta^L_5$.

## 4. Degree 6 ($n=3$)

$q=P$, $k=k_0$, $m=5$; with the selectors $(c_4,c_N,c_O,c_M,\epsilon_c)=(1,0,1,1,1)$ (`production_gamma6_comparison.A_phase`):

$$
R^{(6)}_A(A)=R^{\rm red}(A)+\tfrac14\,P_\omega\smile_sA+\mathsf h\bigl((s\smile s)\smile\operatorname{Sq}^2a\bigr)+2\,\mathsf h\bigl(\omega\smile\operatorname{Sq}^2a\bigr)+\mathsf h\bigl((s\smile s\smile\omega)\smile a\bigr)+\tfrac23\,\mathrm L\bigl(P^1_s\rho_3A\bigr),
\qquad Q_6(A,B)=Q^{\rm gen}(A,B),
$$

$$
S_6(A,A')=S^{\rm aux}(A,A')+R^{(6)}_A(A)+R^{(6)}_A(A')-R^{(6)}_A(A+A')+\mathsf h\bigl(\mathfrak J(A,A')\bigr),
$$

$\mathcal P_6(A,A')$ = the universal primitive of $S_6$ on the degree-three
pair model with the 12 fixed coefficients, and (`production_gamma6.phase`)

$$
\boxed{\;
e^{\rm prod}_6(u,v)=-\mathcal P_6(A,A')-I\,\mathrm{RFS}(u,v)-\mathsf h\bigl(\mathfrak j(A,B;A',B')\bigr)-\bigl[Q_6(A,B)+Q_6(A',B')-Q_6(A_S,B_S)\bigr].\;}
$$

When one of $A,A'$ is zero, $\mathcal P_6$ vanishes identically. When both
are nonzero its evaluation nests three $V_3$ contractions per universal
simplex, which is why the previous implementation refused it by default.
On the diagonal $A'=A$, which is all that the relations of A generators
need, the two-fiber primitive is replaced by a one-fiber one:
[power_relations.md](power_relations.md).

## 5. Normalization, rephasing and the product actually used

For $k=4,5,6$ (`pure_c_normalization`, `all_cochain_upper`, `extension_native_upper`):

$$
\chi_k(C,C')=\mathsf h\bigl(h^2(C,C')\bigr)-e^{\rm prod}_k(0,0,C;\,0,0,C'),\qquad
e^{\rm norm}_k=e^{\rm prod}_k+\chi_k,\qquad
\Gamma^{\rm norm}_k=\Delta\Omega_k+\delta_se^{\rm norm}_k .
$$

$\chi_k$ is closed modulo one, so $N_k=\chi_k-\delta_sI\chi_k(C\ell,C'\ell)$
is integral and $\Gamma^{\rm norm}_k=\Gamma^{\rm prod}_k+\delta_sN_k$.
Define the integral closed cochains

$$
K_k(A,B,C)=\begin{cases}\delta_s\mathsf h(s\smile C),&A=0,\ B=0,\ \delta C=0\ \text{on }X,\\0,&\text{otherwise},\end{cases}
\qquad K_7(c)=\delta_s\mathsf h(s\smile c).
$$

**Full defining systems:** $\Gamma^L_k=\Gamma^{\rm norm}_k-\Delta K_k$.

**Legal AB, arbitrary $C$** (`closed_ab_upper.gamma`): with
$\tau=\tau'_n(A;B)$, $t=\delta C+\tau$, $\Phi_k=\Omega_k+\mathsf h(t\smile_{k-1}\tau)$,

$$
\Gamma^{\rm cl}_k(u,v)=\Phi_k(u)+\Phi_k(v)-\Phi_k\bigl(\mu(u,v)\bigr)+\delta_se^{\rm prod}_k(u,v)+\mathsf h\bigl(h^2(t,t')\bigr),
\qquad \Gamma^{\rm cl,norm}_k=\Gamma^{\rm cl}_k+\delta_s\chi_k ,
$$

and the correction used when both inputs have legal lower pairs is

$$
\boxed{\;
\gamma_k(u,v)=\Gamma^{\rm cl,norm}_k(u,v)-\Bigl[\mathsf h(s\smile t)+\mathsf h(s\smile t')-\mathsf h(s\smile t_S)\Bigr]-K_k(u)-K_k(v)+K_k\bigl(\mu(u,v)\bigr),\qquad k=4,5,6,\;}
$$

($t_S=t+t'$; the bracket is an integral carry). The successor on pure
images $(0,0,c)$, $(0,0,c')$, with the pure-C product

$$
\gamma_C(c,c')=\mathsf h(E(c))+\mathsf h(E(c'))-\mathsf h(E(c+c'))+\delta_s\mathsf h\bigl(h^2(c,c')\bigr)+\mathsf h\bigl(h^2(\delta c,\delta c')\bigr),
$$

is $\gamma_{k+1}(c,c')=\gamma_C(c,c')-K_{k+1}(c)-K_{k+1}(c')+K_{k+1}(c+c')$.
The pure-C product $\gamma_C$, the rephasing $K$ and the carry are natural
for any cup-$i$ structure on a resolution; the pair phases $e^{\rm prod}_k$
are not ([natural_evaluation.md](natural_evaluation.md), Sections 3.2 and 5).

The previous native engine used, on this branch, $\Gamma^{\rm cl}_k$ in place of
$\Gamma^{\rm cl,norm}_k$ (dropping the integral coboundary $\delta_sN_k$, which a
D-gauge absorbs), and for two full defining systems with $A=B=0$ it used
$\gamma_C$ as the base (an integral coboundary away). With $A=0$ it
evaluated $g_k$ and $\Omega_k$ by the direct formulas of Section 7.

## 6. Inputs with a non-legal lower pair, $k=4,5$

With $\widehat F_k$, $\phi$, $q_{\rm loc}$, $\beta^{\rm raw}_k$ as in
[product_lower.md](product_lower.md) Section 5 and $\widehat\mu_k$ the
product with $\beta^{\rm raw}_k$ (`natural_upper`):

$$
\widehat g_k(u)=-I\,\mathcal K_{k+1}\bigl(\widehat F_k(u\ell)\bigr),\qquad
\widehat\Phi_k(u)=I\,\Omega_{k+1}\bigl(\widehat F_k(u\ell)\bigr),\qquad
\widehat\gamma_k(u,v)=I\Bigl[\widehat g_k+\widehat g_k'-\widehat g_k\circ\widehat\mu_k+\Gamma^{\rm norm}_{k+1}\bigl(\widehat F_k\,\cdot,\widehat F_k\,\cdot\bigr)\Bigr](u\ell,v\ell),
$$

$$
\boxed{\;\gamma_k(u,v)=\widehat\gamma_k(\phi u,\phi v)+Z_k(u)+Z_k(v)-Z_k(\mu(u,v))-K_k(u)-K_k(v)+K_k(\mu(u,v))\quad\text{otherwise},\;}
$$

with the integral comparison (`all_cochain_upper.integer_gauge`)

$$
Z_k(u)=\begin{cases}
\Phi_k(u)-\widehat\Phi_k(u)-\mathsf h(s\smile t)+[k=5]\,\tfrac83A^{\smile3}-\delta_sq_k(u), & (A,B)\in L_k,\\[3pt]
-\,I\Bigl[\mathsf h\bigl(E(c^{\rm raw}+\delta M)\bigr)-\mathsf h\bigl(E(c^{\rm raw})\bigr)-\delta_s\mathcal T(c^{\rm raw},M)\Bigr], & \text{otherwise},
\end{cases}
$$

where on $X\times I$ $c^{\rm raw}=\delta(C\ell)+H_k(A\ell,B\ell)$ and
$M=q_{\rm loc}(A,B)\ell+q_{\rm loc}(A\ell,B\ell)$ (`off_branch_integer_gauge`).

**The gauges $q_4,q_5$** (`nonzero_offshell_comparison.legal_pair_gauge4/5`).
Let $\Psi(u)=\Phi_k(u)-\widehat\Phi_k(u)-\mathsf h(s\smile t)$ and define the
output gauge and the closed comparison, for signed-closed $A$ and arbitrary $B,C$ with curvature $(b,c)$:

$$
K^{\rm out}(A,B,C)=\mathfrak q^{\rm HTC}_{m}(b)+\mathsf h\bigl(h^2(c,r_b)\bigr)+\mathcal T\bigl(c+r_b,\ L_k(A,B)\bigr)-\Theta^{\rm pair}_m(b,c;b,c),\qquad r_b=Q^1(b)+s\smile b,
$$
$$
L_k(A,B)=\xi_k(A,B)+q_{\rm loc}(A,B)+P\smile_mb,\qquad
\mathcal D(A,B,C)=\widehat\Phi_k(A,B,C)-\varphi(A,B,C)-K^{\rm out}(A,B,C),\qquad S_A(A)=\mathcal D(A,0,0),
$$

with $\mathfrak q^{\rm HTC}_m$ the higher $\Theta$ comparison primitive
(universal-primitives page, Section 9). Then, with $V^{\rm PS}$ the primary
suspension primitive (same page), $\mathfrak q^U(B)=I\bigl[\mathsf h(Q^1(Q_D(B\ell))+E(s\smile B\ell))-V^{\rm PS}(\delta(B\ell))\bigr]$ and the five-cell
A-only sources

$$
\mathcal D^{(4)}(A)=R^{(4)}_A(A)-S_A(A)+V^{\rm PS}(P)-\mathsf h(s\smile k_0),\qquad
\mathcal D^{(5)}(A)=R^{(5)}_A(A)-S_A(A)-\mathsf h(s\smile k_0)-\mathsf h(a\smile P)-\mathsf h\bigl((s\smile s)\smile P\bigr)+\tfrac83A^{\smile3},
$$

whose universal primitives are $\mathfrak q^A_4(A)$ (dihedral group model) and
$\mathfrak q^A_5(A)$ (degree-two signed model), one has

$$
\begin{aligned}
q_4(A,B,C)&=I\,\Psi(\pi^*A,\pi^*B,C\ell)+\mathfrak q^A_4(A)+\mathfrak q^U(B)-\mathsf h\bigl(s\smile\lambda(A,B)\bigr)+Q_4(A,B,0)-I\,\mathcal D(\pi^*A,B\ell,0),\\
q_5(A,B,C)&=I\,\Psi(\pi^*A,\pi^*B,C\ell)+\mathfrak q^A_5(A)+\mathsf h(a\smile B)+\mathsf h\bigl((s\smile s)\smile B\bigr)-\mathsf h\bigl(s\smile\lambda(A,B)\bigr)+Q_5(A,B)-I\,\mathcal D(\pi^*A,B\ell,0).
\end{aligned}
$$

## 7. The direct formulas at $A=0$ (`a0_degree4`, `a0_high_gamma`)

For a legal triple with $A=0$ (closed $B$) the phase and the pair phase
have a direct form without any A-source. With $m=n+1=|B|$ and the
boundary phases

$$
\mathrm{bp}_4(B)=-I\,\theta_3\bigl(\delta(B\ell),Q_D(B\ell)\bigr)+\mathsf h\bigl(E(Q^1B)+s\smile Q_DB\bigr),\quad
\mathrm{bp}_5(B)=-I\,\theta_4\bigl(\delta(B\ell),Q_D(B\ell)\bigr)+\mathsf h\bigl(E(Q^1B)+s\smile Q_DB\bigr),\quad
\mathrm{bp}_6(B)=I\,\theta_5\bigl(\delta(B\ell),Q_D(B\ell)\bigr)-\mathsf h\bigl(E(s\smile B)\bigr),
$$

$$
\Omega_k(0,B,C)=\mathsf h(E(C))+\mathrm{bp}_k(B),\qquad
e^{A0}_k(B,C;B',C')=\mathsf h\Bigl(h^2(C,C')+h^2\bigl(C+C',h^D(B,B')\bigr)\Bigr)+\mathcal U^{A0}_k(B,B'),
$$

where $\mathcal U^{A0}_k(B,B')$ is the universal primitive of

$$
S^{A0}_k(B,B')=\mathrm{bp}_k(B+B')-\mathrm{bp}_k(B)-\mathrm{bp}_k(B')+\mathsf h\Bigl(h^2(Q_DB,Q_DB')+E\bigl(h^D(B,B')\bigr)+h^2\bigl(Q_DB+Q_DB',\delta h^D(B,B')\bigr)\Bigr)
$$

(Section 7 of the universal-primitives page: $k=4$ on the diagonal
product model, $k=5,6$ on the reduced binary cocycle model). Then
$\gamma^{A0}_k=\Delta\Omega_k+\delta_se^{A0}_k$ on legal data, and for a
lower-legal triple with $A=0$ and arbitrary $C$,
$g_k(0,B,C)=\delta_s\bigl[\Omega_k(0,B,C)+\mathsf h(t\smile_{k-1}Q_DB)\bigr]-\mathsf h(E(t))$,
$t=\delta C+Q_DB$, which is the pure-C rule (P) when $B=0$.

**The production phase at $A=0$, written out.** With both A layers zero
every A-source term of $e^{\rm prod}_k$ is identically zero: $\mathcal P_k(0,0)$,
$\lambda$, $k_0$, $\alpha$, the transgression, $\mathcal U^\Theta_m(0,u)$
inside $Q_k$ and $\mathcal V_n(0,0)$ inside $\beta^\sharp$. With
$f_0(B)=Q_D(B)+\delta(Q^1B+s\smile B)$ and, on $X\times I$, $\beta=h^D(B\ell,B\ell)+\delta(B\ell)\smile_{k-1}\delta(B\ell)$,
$c_u=\delta(C\ell)+f_0(B\ell)$, $L=Q^1(B\ell)+s\smile B\ell$:

$$
\begin{aligned}
e^{\rm prod}_k(0,B,C;0,B,C)&=-I\Bigl[2\varphi_0(B\ell,C\ell)-\mathsf h(E(\beta))-\mathsf h\bigl(h^2(c_u,c_u)\bigr)-\mathcal U^\Theta_{k-1}\bigl(\delta(B\ell),\delta(B\ell)\bigr)\Bigr]-\mu_k(B)-2Q_k(0,B),\\
\varphi_0(B,C)&=\mathsf h\bigl(E(C)+h^2(\delta C,f_0(B))\bigr)-I\,\theta_{k-1}\bigl(\delta(B\ell),f_0(B\ell)\bigr),\qquad
\mu_5=\mathsf h\bigl(s\smile h^D(B,B)\bigr),\quad \mu_4=\mu_6=\mathsf h\bigl(j_0(B,B)\bigr),\\
Q_k(0,B)&=\mathsf h\bigl(h^2(Q^1B,s\smile B)\bigr)-I\,\mathsf h\bigl(E(L)+h^2(Q_D(B\ell),\delta L)\bigr)\quad\bigl[k=5:\ -I\,\mathsf h(E(Q^1(B\ell)))+\tfrac14\widetilde{Q_D(B)}\bigr],
\end{aligned}
$$

equal to the production phase on every simplex (random inputs, $k=4,5,6$;
`a0_square_phase` of `python/light_rows.py`, used by the B-over-D row).
On legal data it differs from $e^{A0}_k(B,C;B,C)$ above by quarter-integers
pointwise, so the two are distinct primitives and the direct form is not a
substitute for the production one ([power_relations.md](power_relations.md), Section 11).

## 8. Degree 6 with a non-legal lower pair

$\gamma_6$ and the successor $\gamma_7$ on non-pure images are the section
construction of [degree_six_section.md](degree_six_section.md), built from
$\Gamma^L_6$, $\Gamma^{\rm sec}_6:=\gamma_6$ of Section 5, $T_3$, $J_6$ and
two legal cylinders.
