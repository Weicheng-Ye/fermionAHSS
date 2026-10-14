# The tertiary differential `T`

`T` is $d_5:E_5^{n,0}\to E_5^{n+5,-4}$ in input degrees $0\le n\le3$ (package
degrees $3$ to $6$). It is the signed coboundary of a rational phase
$\mathcal O_n$ built from the defining system, a universal source primitive
$V_n$ and calibrated selectors. Notation: [conventions.md](conventions.md);
the secondary data $\tau',k_0,g,\lambda$: [secondary_operations.md](secondary_operations.md).
This page transcribes `low_phases.py`, `high_phase.py`, `phase_eval.theta`,
`mod3_power.py` and the universal-helper note. The selectors are
$\boldsymbol\zeta=(\zeta_1,\zeta_2,\zeta_3)=(0,1,0)$, $(c_4,c_N,c_O,c_M,\epsilon_c)=(1,0,1,1,1)$,
$\xi=\tfrac34$, rank correction $0$, $\mu_R=0$, prime-three coefficient $2$.

## 1. Defining system and the returned cocycle

Fix $s,\omega$ and a signed cocycle $A\in Z^n(X;\mathbb Z_s)$, $a=\rho A$. Choose

$$
b\in C^{n+1}(X;\mathbb F_2),\quad c\in C^{n+2}(X;\mathbb F_2),\qquad \delta b=D(a)=\operatorname{Sq}^2a+\omega a+s\rho(\delta\widetilde a/2),\quad \delta c=\tau'_n(A;b).
$$

(The solver changes $b$ by a closed cochain when needed to make $\tau'_n(A;b)$ exact.) Then

$$
\boxed{\;T_n(A;b,c)=\delta_s\mathcal O_n(A,b,c),\qquad
\mathcal O_n=\mathsf h\bigl(E(c)\bigr)+\widehat R_n(A,b)+\mathbf 1_{n\ge2}\,\tfrac23\,\mathrm L\bigl(P^1_s\rho_3A\bigr),\;}
\tag{T}
$$

an integral cocycle of degree $n+5$ (the lifts are specified below, so the
cochain, not only its class, is fixed; an integral change of a phase
changes $T_n$ by a coboundary). The page domain is the classes in
$\ker\operatorname{Dbar}$ with vanishing `Tau`; the target is
$H^{n+5}(\mathbb Z_s)/\bigl(\operatorname{Dtilde}H^{n+2}(\mathbb F_2)+\operatorname{Psi}_{n+1}\ker(D:H^{n+1}\to H^{n+3})\bigr)$.
In degrees $1,2,3$ the kernel omits a final term $\delta_s\mathcal D_n$ of
the notes whose further differential is zero.

## 2. The boundary phase $\Theta_m$, the prism and the source data

For a closed binary $q$ of degree $m$ and $v$ with $\delta v=D(q)$, with
$B_q=\tfrac12\delta\widetilde q$, $e_q=\rho B_q$, $C_q=\tfrac12(B_q+\widetilde{e_q})$,
$u_q=\operatorname{Sq}^2q$, $v_q=\omega q$, $w_q=se_q$:

$$
H_m(q)=\zeta_{2,m}(\omega,q)+\chi_m(q)+u_q\smile_{m+1}v_q+u_q\smile_{m+1}w_q+v_q\smile_{m+1}w_q+\zeta_{1,m+1}(s,e_q)+(\omega\smile_1s)e_q+su_q+s\smile s\smile\rho C_q,
$$

$$
\boxed{\;\Theta_m(q,v)=\mathsf h\bigl(E(v)+H_m(q)\bigr)+\tfrac14\bigl(\widetilde\omega\smile B_q+B_q\smile_{m-1}B_q\bigr)+\mathsf h\bigl(s\operatorname{Sq}^2q+\omega\operatorname{Sq}^1q\bigr).\;}
\tag{B}
$$

$E(v)+H_m(q)$ is lifted together and the last binary sum once; the
quarter term uses integral cup signs. $\Theta_m(q,0)=:R_{\partial,m}(q)$
(this is $\theta_1$ of [differential.md](differential.md) for $m=1$, where
$\chi_1=0$, $\operatorname{Sq}^2q=0$ and $B_q=\widetilde{q\smile q}$).
$\Theta_m$ is a rational cochain of degree $m+3$ whose coboundary is integral on its domain.

The prism $I$ is (C11); $b\ell$ is the ramp; every lift and carry inside
$\Theta$ is recomputed on each prism simplex, and the prism of a rational
phase is the exact signed sum (no reduction modulo one).

The A-only data are (S) and (S23) of the secondary page: $p=D(a)$,
$g$, $k_0$, $\lambda(A,b)=g\cdot(s\smile b)$, $y=Q_D(b)$, $\tau'=y+k+\delta\lambda$.

## 3. Input degree zero

$\mathcal O_0=\mathsf h(E(c))+R_0(A,b)$; $\rho A$ is locally constant, and the
branch is chosen by the parity of $A$ on each component.

**Even $A$.** $K=A/2$, $t=\rho K$, $q_K=(K-\widetilde t)/2$, $Z=t\omega$,
$N_b=\tfrac12(\widetilde Z+\widetilde{sb}-\widetilde{Z+sb})=Z\cdot(s\smile b)$,
$\lambda=\rho\bigl(q_K\smile\widetilde\omega+N_b\bigr)$, $H=Q_D(b)$,
$r=\rho(\delta\widetilde\omega/2)$, $z=tr$:

$$
\boxed{\;R_0(A,b)=\Theta_1(b,0)-\tfrac18\,K\smile P_\omega+\mathsf h\bigl(H\smile_2z+E(\lambda)+(H+z)\smile_2\delta\lambda\bigr).\;}
\tag{R0e}
$$

**Odd $A$.** $N=|A|$, $\epsilon_A=[A<0]$, $t=\rho((N-1)/2)$,
$\lambda_+=b^2+tQ^1b$, $R_+(N,b)=-\mathsf h(E(\lambda_+))+\tfrac{3+2\widetilde t}8\,\delta(\widetilde b^{\,3})$
(ordinary integral cube). On the cylinder with the gauge $\gamma=\epsilon_A\ell$:
$s_I=s+\delta\gamma$, $A_I=(-1)^\gamma A$, $b_I=\pi^*b$, and $\tau_I=\tau'_0(A_I;b_I)$ computed with $(s_I,\pi^*\omega)$:

$$
\boxed{\;R_0(A,b)=(-1)^{\epsilon_A(v_0)}R_+(N,b)+\mathsf h\bigl(\rho\,I\,E(\tau_I)\bigr).\;}
\tag{R0o}
$$

(The parity of the prism sum is taken before the half lift; $s=\delta\epsilon_A$ on an odd component, so the top end is the positive oriented input.)

**Page classes with trivial sign** (native shortcuts, $s=0$): for $A=2K$
with $K$ odd and a $c$ with $\delta c=-\rho(\delta W/2)$, $W=\widetilde\omega$,
the class is $[T_0(A)]=-K\,\beta_8P_4(W+2\widetilde c)$ with
$P_4(z)=z\smile z+z\smile_1\delta z$ modulo eight; for $A=4m$,
$[T_0(A)]=-m\,\beta_4P_2(\omega)$; for positive odd $A$ the class is zero.

## 4. Input degree one

With the adjusted source of Section 6 of the secondary page ($p_0=p$):

$$
\begin{gathered}
e_0=s\smile t,\quad b^+=b+e_0,\quad q=(\omega+s\smile s)a=\delta b^+,\quad
k=k_0+Q_D(e_0)+\delta e_0\smile_2p_0+s\smile(\delta e_0\smile_3p_0),\quad y=Q_D(b^+),\\
\lambda=\lambda(A,b)+b\smile_1e_0+p_0\smile_2e_0+p_0\smile_3\delta e_0+s\smile\bigl(b\smile_2e_0+p_0\smile_3e_0+p_0\smile_4\delta e_0\bigr),\quad
\tau_1=y+k+\delta\lambda,\quad \operatorname{Pol}_1(y,k)=y\smile_3k+k\smile_4\delta y+Q^1k,
\end{gathered}
$$

$$
\boxed{\;\widehat R_1=-I\,\Theta_3\bigl(\delta(b^+\ell),Q_D(b^+\ell)\bigr)-V_1
+\mathsf h\bigl(\operatorname{Pol}_1(y,k)+E(\lambda)+\tau_1\smile_3\delta\lambda+E(Q^1b^+)+s\smile\tau_1\bigr),\qquad \mathcal O_1=\mathsf h(E(c))+\widehat R_1 .\;}
\tag{R1}
$$

$V_1=V_1(A,s,\omega)$ is the degree-five source primitive of
$\Theta_3(q,k)$ (Section 8). The rank term $-\zeta_1\tfrac14P_\omega\smile_sA$ is absent since $\zeta_1=0$.

## 5. Input degree two

With $k=k_0$, $y=Q_D(b)$, $\lambda=\lambda(A,b)$, $\tau_2=y+k+\delta\lambda$,
$\operatorname{Pol}_2(y,k)=y\smile_4k+k\smile_5\delta y+Q^1k$,
$U_A=\zeta_{1,4}(s,p)+s\smile k+(\omega\smile_1s)\smile p$,
$\Pi_2=I\,\Theta_4(\delta(b\ell),Q_D(b\ell))$, $L_2=\tfrac14P_\omega\smile_sA$:

$$
\boxed{\;\widehat R_2=\mathsf h(\operatorname{Pol}_2(y,k))+\Pi_2-V_2+\mathsf h\bigl(E(\lambda)+(y+k)\smile_4\delta\lambda\bigr)-\mathsf h\bigl(E(s\smile b)\bigr)+\mathsf h(U_A)+L_2-\tfrac14A^{\smile3},
\qquad \mathcal O_2=\mathsf h(E(c))+\widehat R_2+\tfrac23\,\mathrm L\bigl(\rho_3A^{\smile3}\bigr).\;}
\tag{R2}
$$

Each $\mathsf h$ is a separate half lift. The cubic term is present
($R_2=R_2^\sharp$); $\zeta_2=1$ puts $L_2$ in, $\zeta_3=0$ removes a term
$\tfrac12\widetilde{s^2\omega a}$. On a six-simplex
$A^{\smile3}(v_0..v_6)=(-1)^{s_{02}+s_{04}}A_{012}A_{234}A_{456}$. The three-primary
term is $\tfrac23$ times the $\{0,1,2\}$ lift of the cube modulo three; its
coboundary is exact on cocycles, so it does not change $T_2$ on the page,
but it is kept for the stacking correction.

## 6. Input degree three

With $\operatorname{Pol}_3(y,k)=y\smile_5k+k\smile_6\delta y+Q^1k$,
$\Pi_3^+=I\,\Theta_5(\delta(b\ell),Q_D(b\ell))$, $L_3=\tfrac14P_\omega\smile_sA$,
$N_3=\mathsf h(s\smile s\smile\operatorname{Sq}^2a)$, $O_3=\mathsf h(\omega\smile\operatorname{Sq}^2a)$,
$M_3=\mathsf h(s\smile s\smile\omega\smile a)$:

$$
\boxed{\;\widehat R_3=\mathsf h\bigl(\operatorname{Pol}_3(y,k)+E(\lambda)+(y+k)\smile_5\delta\lambda\bigr)-\Pi_3^+-V_3+\mathsf h\bigl(E(Q^1b)+s\smile\tau_3\bigr)+L_3+N_3+2O_3+M_3,
\qquad \mathcal O_3=\mathsf h(E(c))+\widehat R_3+\tfrac23\,\mathrm L\bigl(P^1_s\rho_3A\bigr).\;}
\tag{R3}
$$

The coefficients are $c_4=1$ on $L_3$, $c_N+\epsilon_c=1$ on $N_3$,
$c_O+\epsilon_c=2$ on $O_3$, $c_M=1$ on $M_3$; $2O_3$ is integral but kept.
$V_3=V_3(A,s,\omega)$ is the degree-seven source primitive (Section 8).

## 7. The reduced power $P^1_s$ modulo three

For $a_3=\rho_3A$ of degree three with the $\{0,1,2\}$ lift, on an ordered
face $F=(i,j,k,l)$ of a seven-simplex write $a^{(0)}_F=(-1)^{s(v_0,v_i)}a_3(v_i,v_j,v_k,v_l)$ (sign $1$ for $i=0$). Then

$$
P^1_s(a_3)(v_0,\dots,v_7)=\rho_3\sum_{(F_1,F_2,F_3;\nu)}\nu\,a^{(0)}_{F_1}a^{(0)}_{F_2}a^{(0)}_{F_3}
\tag{P3}
$$

over the nineteen terms ($\nu$; $F_1$; $F_2$; $F_3$, faces as vertex strings):

| | | | | | | | |
| --- | --- | --- | --- | --- | --- | --- | --- |
| $-1$;0123;3456;0167 | $-1$;0123;3456;1267 | $-1$;0123;3456;2367 | $-1$;0123;3456;3467 | $-1$;0123;3456;4567 | $-1$;0123;3467;4567 | $+1$;0123;3567;3457 | $-1$;0124;4567;2347 |
| $+1$;0125;2567;2345 | $+1$;0126;2367;3456 | $+1$;0127;2347;4567 | $+1$;0134;4567;1237 | $+1$;0145;1567;1234 | $-1$;0156;1267;2345 | $+1$;0167;1237;3456 | $-1$;0234;4567;0127 |
| $+1$;0345;0567;0123 | $+1$;0456;0167;1234 | $+1$;0567;0127;2345 | | | | | |

This is the tridegree-$(3,3,3)$ part of the cyclic diagonal $D_2$, built by
the recursion $D_i(\sigma)=h\bigl(\mathcal R_iD_{i-1}(\sigma)+(-1)^iD_i(\partial\sigma)\bigr)$,
$\mathcal R_1=\varrho-1$, $\mathcal R_2=1+\varrho+\varrho^2$, from the
threefold Alexander–Whitney $D_0$ ($h$ the tensor contraction to the first
vertex, $\varrho$ the Koszul-signed rotation). In degree two $P^1_s\rho_3A=\rho_3(A^{\smile3})$;
in degrees $0,1$ it vanishes. The Bockstein is positive:
$2\beta_{3,s}P^1_s\rho_3A=\tfrac23\delta_s\mathrm L(P^1_s\rho_3A)$, the
$k$-invariant of $ko$ at the prime three (the factor $2$ from
$c(\alpha_4)=2u^2$).

## 8. The source primitives $V_1,V_2,V_3$

Each $V_n$ is a rational cochain of degree $n+4$ in $(A,s,\omega)$ with
$\delta_sV_n\equiv\Theta_{n+2}(q_n,k_n)\pmod{\mathbb Z}$ ($q_1,k_1$ the
adjusted source, $q_n=p$, $k_n=k_0$ for $n=2,3$), computed by a fixed
contraction of a universal model, never by solving on $X$
([universal_primitives.md](universal_primitives.md), Section 1, is the
general pattern). The universal source $\phi_n=\Theta_{n+2}(q_n,k_n)$ has
degree $n+5$.

### 8.1 Chain comparison conventions

Alexander–Whitney is the sum of front/back cuts; the shuffle map sums
monotone lattice paths with the inversion sign; degenerate terms are
discarded. The product homotopy is the cone recursion
$\widehat H_0=0$, $\widehat H_n=C_{(0,0)}(Q-\widehat H_{n-1}\partial)$, $u=Q\widehat HQ$,
$H=u\partial u$ with $Q=1-\mathrm{EZ}\circ\mathrm{AW}$ (the dihedral model
uses the raw $\widehat H$, the matrix models the normalized $H$). Tensor
and composite homotopies: $H_{X\times Y}=H_{\rm prod}+G_{\rm prod}(H_X\otimes1+G_XF_X\otimes H_Y)F_{\rm prod}$
(Koszul sign on the second term), $H_{21}=H_1+G_1H_2F_1$.

For the signed two-cocycle matrix model the small target has generators
$u_{p,k}$ of degree $p+2k$ with $\partial u_{p,k}=-2u_{p-1,k}$ when $p>0$ and
$p+k$ odd, else $0$; the background target has bar words $[x_{i_1}|\dots|x_{i_l}]$
of degree $\sum(i_j+1)$ of the normalized bar of the $C_2$ chain algebra.
The signed comparison first contracts unit-edge chains to the crossed
exterior algebra $\mathbb Z\{1,\sigma,e,e\sigma\}$ ($\sigma^2=1$, $e^2=0$,
$\sigma e=-e\sigma$; reduced letters $t=\sigma-1$, $e$, $e\sigma$, $t^2=-2t$),
with the recursive comparison $f=h_Pf\partial$, $g=h_Qg\partial$,
$K=h_Q(1-gf-K\partial)$, free boundary
$\partial u_{p,k}=(-1)^peu_{p,k-1}+(\sigma+(-1)^{p+k})u_{p-1,k}$, and contractions
$h_P(\sigma u_{p,0})=u_{p+1,0}$, $h_P(eu_{p,k})=(-1)^pu_{p,k+1}$,
$h_P(e\sigma u_{p,k})=(-1)^{p+1}\sigma u_{p,k+1}$, $h_Q(a[w])=-[a-\epsilon(a)|w]$;
finally the sign module sends $1,\sigma,e,e\sigma$ to $1,-1,0,0$ (`exterior_bar.py`).

### 8.2 $V_1$ (dihedral model)

For signed integral $A$ of degree one, the homogeneous vertices of a
simplex are $g_i=(A_{0i},s_{0i})\in\mathbb Z\rtimes C_2$, $g_0=(0,0)$,
identified with the Cayley line by $(k,\epsilon)\mapsto2k-\epsilon$
(reconstruction $s_{ij}=\epsilon_i+\epsilon_j$, $A_{ij}=(-1)^{\epsilon_i}(k_j-k_i)$);
the background $\omega$ uses the matrix section of 8.3. The retraction and
homotopy on group simplices are

$$
R[g]=K_{g_0}(R\partial[g]),\qquad h[g]=C_{g_0}\bigl([g]-R[g]-h\partial[g]\bigr),\qquad R[g_0]=[g_0],\ h[g_0]=0,
$$

with $K_{\rm root}$ the oriented unit-edge path from the root to a vertex,
and, on an alternating simplex in one unit edge, the prepending of its
endpoint nearest to the root. With $P$ the product retraction ($R$ in the
group factor, identity in the background) and
$K=\widehat H_{\rm prod}+\mathrm{shuffle}(h\otimes1)\mathrm{AW}$,

$$
\boxed{\;V_1=\phi_1K+U_{\rm wedge}P,\;}
$$

each evaluation multiplied by $(-1)^{\epsilon_0}$ when the first vertex is
not in the zero sign component. $U_{\rm wedge}=0$ on the branch with
increment $(A,s)=(0,1)$; on the branch $(1,1)$ (the literal odd branch
$A=\widetilde s$) it is $U_t$: append a prism coordinate $\tau$,
$s_I=\rho(s+\delta\tau)$, $A_I=\widetilde{s_I}$, form the source $\Phi_I$,
and with $I$ the raw signed right prism

$$
\kappa=\rho\,I(\delta_{s_I}\Phi_I),\qquad U_{\rm raw}=-\tfrac12I\Phi_I+\tfrac12\widetilde{\operatorname{Fill}(\kappa)},\qquad
\operatorname{Fill}(\gamma)(z)=\gamma\bigl((H+Gh_{\rm sm}F)j(z)\bigr)\bmod2,
$$

where $F,G,H$ is the product comparison of $BC_2^{\,s}$ with the background
matrix model and the only nonzero $h_{\rm sm}$ rules are
$[3]\mapsto[1|2]$, $[5]\mapsto[1|4]$, $[1|3]\mapsto[1|1|2]$, $[3|1]\mapsto[2|1|1]$.
All eight odd-source detector bits are zero (no quarter-valued Bockstein
corrections). The literal base normalization is essential:

$$
U_t=U_{\rm raw}-\delta_sE_0,\qquad E_0=(U_{\rm raw}|_{s=0})H_\omega+e_0F_\omega,
$$

with source values $(B,C)=(7/8,0)$, base period $0$ (so $\eta=0$), and the
sparse table $e_0([x_3])=-7/16$, $e_0([x_1|x_1])=-7/32$, all other values
zero. $V_1$ grows affinely in the integer labels on residue classes modulo
four and is not periodic.

### 8.3 $V_2$ (signed matrix model)

Two diagonal matrix models: signed integral rows $(\sigma_i,M_i)$ and
binary rows $W_i$, with reconstruction on a face $(r,l,t)$
$s_{rt}=\sum_{i=r}^{t-1}\sigma_i$, $A_{rlt}=\sum_{i=r}^{l-1}\sum_{j=l}^{t-1}(-1)^{s_{ri}}M_{ij}$,
$\omega_{rlt}=\sum\sum W_{ij}$, and the section $j$: $\sigma_i=s_{i,i+1}$,
$M_{ij}=A_{i,i+1,j+1}-A_{i,i+1,j}$, $M_{ji}=-(-1)^{s_{ij}}M_{ij}$,
$W_{ij}=\omega_{i,i+1,j+1}+\omega_{i,i+1,j}$ (diagonals zero; $rj=1$). With
the product contraction $(F_{\rm tot},G_{\rm tot},H_{\rm tot})$ of 8.1,

$$
\boxed{\;V_2=j^*\bigl((r^*\phi_2)H_{\rm tot}+e_fF_{\rm tot}\bigr),\;}
$$

where in the small basis $u_{p,k}\otimes[x_i|\dots]$ the five source values are
$C=\phi(u_{2,1}\otimes[x_2])=0$, $G=\phi(u_{0,1}\otimes[x_1|x_2])=\tfrac34$,
$H=\phi(u_{0,1}\otimes[x_2|x_1])=0$, $A_{\rm src}=\phi(u_{3,2}\otimes[\,])=\tfrac34$,
$B_{\rm src}=\phi(u_{1,2}\otimes[x_1])=\tfrac14$ (modulo one), and the complete degree-six $e_f$ is

| basis | formula | value |
| --- | --- | ---: |
| $u_{2,1}\otimes[x_1]$ | $-C/2$ | $0$ |
| $u_{0,1}\otimes[x_3]$ | $(H-G)/2$ | $-3/8$ |
| $u_{0,1}\otimes[x_1|x_1]$ | $-(G+H)/4$ | $-3/16$ |
| $u_{2,2}\otimes[\,]$ | $-A_{\rm src}/2$ | $-3/8$ |
| $u_{0,2}\otimes[x_1]$ | $-B_{\rm src}/2$ | $-1/8$ |

all other coefficients zero. The suspension normalization uses
$V_{2,\rm fin}=V_2-L_2$ without the cubic term; the assembled $R_2$ has the cubic term separately.

### 8.4 $V_3$ (crossed $U_2$ model)

The signed model is $\operatorname{diag}B(U_2\rtimes C_2)$: a row is an
integral two-cocycle $z_i$ (as a skew matrix, strict upper-triangle section)
with sign $\sigma_i$, multiplication $(z,\epsilon)(z',\epsilon')=(z+(-1)^\epsilon z',\epsilon+\epsilon')$;
reconstruction and section

$$
rA(v_0,v_1,v_2,v_3)=\sum_{i=v_0}^{v_1-1}(-1)^{s(v_0,i)}z_i(v_1,v_2,v_3),\qquad
\bar A(v_0..v_3)=(-1)^{s(0,v_0)}A(v_0..v_3),\qquad
z_i=(-1)^{s(0,i)}\bigl(\iota_i\bar A-\iota_{i+1}\bar A\bigr),\ \sigma_i=s(i,i+1),
$$

($\bar A$ extended alternatingly, zero on repeated vertices). The $U_2$
contraction projects to divided powers $u_k$ (degree $2k$) and keeps the
$C_2$ component ($t=\sigma-1$, $u_k=(0,k)$, $v_k=(1,k)$, suspended degrees
$1,2k+1,2k+1$); the actual differential is the **full** transfer with all
crossed terms: $L=\sum(-\delta h_T)^m$, $R=\sum(-h_T\delta)^m$, $F=f_TL$,
$G=Rg_T$, $H=h_TL$, $d_{\rm small}=f_T\delta Rg_T$ (`r3_chain.py`), tensored
with the background word complex and taken relative to $A=0$. Basis
elements $[\ell_1|\dots|\ell_r]\otimes[x_{j_1}|\dots]$ are ordered by signed
degree, length, lexicographic $(\text{component},k)$, then background degree,
length, word: 46 in degree seven, 94 in degree eight. With $D$ the $46\times94$
boundary matrix and $f=(r^*\phi_3)G_{\rm tot}$ the degree-eight source row,
the deterministic Smith recipe (least nonzero pivot, nonnegative remainders
promoted immediately) gives $UDV=\operatorname{diag}(1^{27},2,2,2,12,0,\dots)$,
rank $31$, all 63 cycle values of $fV$ zero modulo one; for $d_i=2^{v_i}o_i$,
$n_i=4\operatorname{rep}((fV)_i)$, $m_i=o_i^{-1}n_i\bmod4$, $e'_i=m_i/2^{v_i+2}$
(a **dyadic** division, also at the pivot $12$), $e_f=e'U$, and

$$
\boxed{\;V_3=j^*\bigl((r^*\phi_3)H_{\rm tot}+e_fF_{\rm tot}\bigr).\;}
$$

The 34 nonzero source values ($w_k$ stands for $u_k$ or $v_k$, one $w$ per row):
$[w_1]\!\otimes\![x_4]:\tfrac34$, $[w_1]\!\otimes\![x_1|x_2]:\tfrac34$, $[w_2]\!\otimes\![x_2]:\tfrac14$,
$[t|t|w_1]\!\otimes\![x_2]:\tfrac12$, $[u_1|w_1]\!\otimes\![x_1]:\tfrac14$, $[t|w_2]\!\otimes\![x_1]:\tfrac14$,
$[v_1|w_1]\!\otimes\![x_1]:\tfrac34$, $[u_1|w_2]:\tfrac34$, $[u_2|w_1]:\tfrac34$, $[v_1|w_2]:\tfrac34$,
$[v_2|w_1]:\tfrac14$, $[u_1|t|t|w_1]:\tfrac12$, $[t|u_1|t|w_1]:\tfrac34$, $[t|t|t|w_2]:\tfrac14$,
$[t|t|v_1|w_1]:\tfrac12$, $[t|v_1|t|w_1]:\tfrac34$, $[v_1|t|t|w_1]:\tfrac12$.
The 16 nonzero primitive coefficients:
$[w_1]\!\otimes\![x_3]:\tfrac{19}8$, $[w_1]\!\otimes\![x_1|x_1]:-\tfrac{13}{16}$, $[t|w_1]\!\otimes\![x_2]:-\tfrac{11}4$,
$[w_2]\!\otimes\![x_1]:-\tfrac38$, $[w_3]:-\tfrac74$, $[u_1|t|w_1]:-\tfrac12$, $[t|u_1|w_1]:-\tfrac14$, $[t|t|w_2]:\tfrac18$.
The matrices are in `python/r3_source.json`.

### 8.5 Selectors and independence of the contraction

The selectors are evaluated on fixed small cycles: the R1 rank period
(period $2$, $\zeta_1=0$); the R2 suspension periods ($\tfrac14$, orientation
$j=1$, $\zeta_2=1$; $\kappa_6V_2(Ss)=\tfrac32$, $U_{\rm raw}(\mathrm{Tor})=-\tfrac12$, $\zeta_3=0$);
the odd base period $0$; the rank-six $BPSO(6)=BPU(4)$ comparison (rank
correction $0$); the Euler value $\Xi(D)=V_2(D)=\tfrac{13}4$ ($\mu_R=0$); and
the R3 row on the degree-six cycles $F=u_{0,3}$, $Y=u_{0,1}[x_1|x_1]$, $W=u_{2,2}$,
$V=u_{0,2}[x_1]$, $Z+T=u_{1,1}[x_2]+u_{2,1}[x_1]$ with
$C_0=\kappa_7V_3(SA)-V_{2,\rm fin}$ ($\kappa_7=+I$, $\kappa_6=-I$), periods
$(\tfrac34,\tfrac14,0,\tfrac12,\tfrac12)$, and
$\xi=a_FC_0(F)$, $c_4=j_L\,4C_0(Y)\bmod4$, $(c_N,c_O,c_M)=2(C_0(W),C_0(V),C_0(Z+T))\bmod2$,
$\epsilon_c=2(\xi-\tfrac14)\bmod2$, with $a_F=j_L=1$.

Any contraction yields a primitive ($\delta_sV\equiv\phi$ because $\phi$ is
closed modulo one), and two primitives differ by a $\mathbb Q/\mathbb Z$-cocycle
on the universal model whose class lies in a **finite** group: for $V_3$,
$H^7(U,U_0;\mathbb Q/\mathbb Z)=(\mathbb Z/2)^3\oplus\mathbb Z/12$. The dyadic
convention fixes the $\mathbb Z/3$ part, and the universal cocycles
$L_3,N_3,O_3,M_3$ (values on the four torsion cycles $(0,0,0,\tfrac34)$,
$(\tfrac12,0,0,0)$, $(0,\tfrac12,0,0)$, $(0,0,\tfrac12,\tfrac12)$) realize every
two-primary class exactly once, so the selector family compensates the
choice of contraction: the operation is the same for every admissible
contraction once its selectors are evaluated by the suspension periods.

## 9. Scope

$R$ is implemented in input degrees $0$–$3$ only; no all-degree family is
asserted. GAP solves only the defining cochains $b,c$ and transports the
whole phase to the bar; a failed identity is an error. The three native
degree-zero reductions of Section 3 bypass the transport where their
hypotheses hold.
