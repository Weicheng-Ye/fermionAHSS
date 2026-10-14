# The generalized differential, degree by degree

A state of package degree $k$ is $x=(A,B,C,D)$ with

$$
A\in C^{k-3}(X;\mathbb Z_s),\quad B\in C^{k-2}(X;\mathbb F_2),\quad
C\in C^{k-1}(X;\mathbb F_2),\quad D\in C^{k+1}(X;\mathbb Z_s),
$$

one cochain per layer of the line $p+q=k-3$ of the spectral sequence
(A on row $q=0$, B on $q=-1$, C on $q=-2$, D on $q=-4$; row $-3$ is zero).
Layers of negative degree are absent. Write $u=(A,B,C)$ and $n=k-3$.
Notation: [conventions.md](conventions.md).

The differential has the same shape in every degree,

$$
\boxed{\;
\mathfrak d_k(A,B,C,D)=\bigl(\delta_sA,\ \ \delta B+P(A),\ \ \delta C+f_k(A,B),\ \ \delta_sD+g_k(A,B,C)\bigr),
\qquad P(A)=Q_D(\rho A),\;}
$$

and this page writes $f_k$ and $g_k$ out in each degree. Three general
definitions first.

**Legal locus.** $p_k(A,B)=(\delta_sA,\ \delta B+P(A))$ and
$L_k=\{(A,B):p_k(A,B)=0\text{ on all of }X\}$. Membership is a global
test, and the formulas branch on it. By (C7), $p_{k+1}\circ p_k=0$.

**Flat states.** $\mathfrak d_kx=0$ means
$\delta_sA=0$, $\delta B=P(A)$, $\delta C=\tau_k(A,B)$, $\delta_sD=-\mathcal K_k(A,B,C)$:
the defining-system equations of the spectral sequence, with
$\mathcal K_k$ the tertiary cocycle of Section 2.

**The residual and the legal D component.** For $(A,B)\in L_k$ and
arbitrary $C$, with $\tau_k$ and $\Omega_k$ the cocycle and phase of each
degree below,

$$
\boxed{\;
t_k=\delta C+\tau_k(A,B),\qquad
\Phi_k=\Omega_k+\mathsf h\bigl(t_k\smile_{k-1}\tau_k\bigr),\qquad
J_k=\delta_s\Phi_k-\mathsf h\bigl(E(t_k)\bigr).\;}
\tag{J}
$$

$J_k$ is integral, $\delta_sJ_k=-\delta_s\mathsf h(E(t_k))$, and when
$t_k=0$, $J_k=\delta_s\Omega_k=:\mathcal K_k$. The only property of
$\Omega_k$ that this uses is the residual identity
$\delta_s\Omega_k\equiv\mathsf h(E(\delta C))+\mathsf h(E(\tau_k))\pmod{\mathbb Z}$
on $L_k$.

The phases $R_0$ (even and odd input), $\widehat R_1,\widehat R_2,\widehat R_3$
and the reduced power $P^1_s$ are written out in
[tertiary_operations.md](tertiary_operations.md), Sections 3–7 ((R0e), (R0o),
(R1), (R2), (R3), (P3)); the secondary representative $\tau'_n(A;B)$ is (S11)
and (S23) of [secondary_operations.md](secondary_operations.md).

## Degrees $-1$ and $0$

Only $D$ is present ($D\in C^0$, resp. $C^1$):

$$
\boxed{\;\mathfrak d_{-1}(D)=\delta_sD,\qquad \mathfrak d_0(D)=\delta_sD.\;}
$$

## Degree 1

State $(C,D)$, $C\in C^0(\mathbb F_2)$, $D\in C^2(\mathbb Z_s)$. Here
$\tau_1=0$, $\Omega_1=\mathsf h(\omega\smile C)$, $E(C)=\omega\smile C$:

$$
\boxed{\;
\mathfrak d_1(C,D)=\Bigl(\delta C,\ \ \delta_sD+\tfrac12\bigl(\delta_s\widetilde{\omega\smile C}-\widetilde{\omega\smile\delta C}\bigr)\Bigr).\;}
$$

## Degree 2

State $(B,C,D)$, $B\in C^0(\mathbb F_2)$, $C\in C^1(\mathbb F_2)$,
$D\in C^3(\mathbb Z_s)$. The C component is the same formula whether or
not $B$ is closed:

$$
f_2(B)=H_2(B)=Q_D(B)=\omega\smile B+s\smile(B\smile\delta B).
$$

With $p=\omega\smile B$ and $t_2=\delta C+p$,

$$
\boxed{\;
\mathfrak d_2(B,C,D)=\bigl(\delta B,\ \ \delta C+H_2(B),\ \ \delta_sD+g_2(B,C)\bigr),\;}
$$

$$
g_2(B,C)=\begin{cases}
\delta_s\Bigl[\mathsf h(C\smile\delta C+\omega\smile C)+\mathsf h(t_2\smile_1p)\Bigr]
-\mathsf h\bigl(t_2\smile t_2+\omega\smile t_2\bigr), & \delta B=0\text{ on }X,\\[6pt]
-\,I\,\mathcal K_3\bigl(0,\ \delta(B\ell),\ \delta(C\ell)+H_2(B\ell)\bigr), & \delta B\ne0,
\end{cases}
$$

where in the second line everything lives on $X\times I$ and
$\mathcal K_3(0,b,c)=\delta_s\Omega_3(0,b,c)$ with
$\Omega_3(0,b,c)=\mathsf h(E(c))+\theta_1(b)$ (Section "Degree 3" below,
with $A=0$; the two summands are lifted separately). The argument of
$\mathcal K_3$ is a full defining system on the cylinder, so the integrand
is an integral cocycle and $g_2$ is integral.

## Degree 3

State $A\in C^0(\mathbb Z_s)$, $B\in C^1$, $C\in C^2$, $D\in C^4$.
Write $a=\rho A$, $t=t(A)$.

**C component.**

$$
f_3(A,B)=\begin{cases}\tau'_0(A;B)\ \text{of (C12)}, & (A,B)\in L_3,\\[2pt]
H_3(A,B)=I\,\tau'_1\bigl(\delta_s(A\ell);\ \delta(B\ell)+P(A\ell)\bigr), & \text{otherwise.}\end{cases}
$$

**Phase.** $\Omega_3(A,B,C)=\mathsf h(E(C))+R_0(A,B)$ with
$E(C)=C\smile C+C\smile_1\delta C+\omega\smile C$. The degree-zero tertiary
phase $R_0$ is evaluated componentwise by the parity of $A$ (the sign
cocycle condition forces $|A|$ to be locally constant). With
$P_\omega$ of (C3), $r=\rho(\delta\widetilde\omega/2)$ and the lifts taken
separately:

*Even $A$:* $K=A/2$, $t=\rho K$, $\kappa=(K-\widetilde t)/2$,
$Z=t\smile\omega$, $\lambda^{\rm ev}=\kappa\smile\omega+Z\cdot(s\smile B)$
(the first term integral, then reduced), $H=Q_D(B)$, $z=t\smile r$:

$$
\Omega_3^{\rm ev}=\mathsf h(E(C))+\theta_1(B)-\tfrac18\,K\smile P_\omega
+\mathsf h\Bigl(H\smile_2z+E(\lambda^{\rm ev})+(H+z)\smile_2\delta\lambda^{\rm ev}\Bigr).
$$

*Odd $A$:* with $\ell_A=[A<0]$, $t^{\rm abs}=\rho\bigl((|A|-\widetilde{\rho|A|})/2\bigr)$,
$\lambda^{\rm od}=B\smile B+t^{\rm abs}\smile Q^1(B)$, $c_B=\widetilde B\smile\widetilde B\smile\widetilde B$
(integral), and on $X\times I$ the orientation gauge $\gamma=\ell_A\ell$,
$s_I=\pi^*s+\delta\gamma$, $A_I=(-1)^{\gamma}\pi^*A$, $\psi_I=\tau'_0(A_I;\pi^*B)$
computed with the twists $(s_I,\pi^*\omega)$:

$$
\Omega_3^{\rm od}=\mathsf h(E(C))+(-1)^{\ell_A}\Bigl[-\mathsf h\bigl(E(\lambda^{\rm od})\bigr)+\tfrac{3+2t^{\rm abs}}8\,\delta c_B\Bigr]
+\tfrac12\,\rho\bigl(I\,E(\psi_I)\bigr).
$$

Here $\theta_1(b)$ is the degree-one boundary phase, fully explicit for a
closed binary $b$ of degree one, with $b^2=b\smile b$:

$$
\boxed{\;
\theta_1(b)=\tfrac12\Bigl(\zeta_{2,1}(\omega,b)+(\omega b)\smile_2(s\,b^2)+\zeta_{1,2}(s,b^2)+(\omega\smile_1s)\,b^2+s\,s\,b^2\Bigr)^{\sim}
+\tfrac34\,\widetilde{\omega\,b^2}+\tfrac14\,\widetilde{b^2\smile b^2}.\;}
$$

**D component.** On $L_3$, $J_3$ of (J) with $\tau_3=\tau'_0(A;B)$:

$$
g_3=\delta_s\Bigl[\Omega_3+\mathsf h(t_3\smile_2\tau'_0)\Bigr]-\mathsf h\bigl(t_3\smile_1t_3+\omega\smile t_3\bigr),\qquad t_3=\delta C+\tau'_0(A;B).
$$

Outside $L_3$:

$$
g_3=G_3=-\,I\,\mathcal K_4(a_I,b_I,c_I),\qquad
(a_I,b_I)=p_3(A\ell,B\ell),\quad c_I=\delta(C\ell)+H_3(A\ell,B\ell),
$$

with $\mathcal K_4=T_1=\delta_s\Omega_4$ on the full defining system
$(a_I,b_I,c_I)$ of the cylinder ($H_3$ uses a second interval).

## Degree 4

State $A\in C^1(\mathbb Z_s)$, $B\in C^2$, $C\in C^3$, $D\in C^5$.

$$
f_4(A,B)=\begin{cases}\tau'_1(A;B), & (A,B)\in L_4,\\
H_4(A,B)=I\,\tau'_2\bigl(\delta_s(A\ell);\ \delta(B\ell)+P(A\ell)\bigr), & \text{otherwise,}\end{cases}
\qquad
\Omega_4=\mathsf h(E(C))+\widehat R_1(A,B),\quad E(C)=C\smile_1C+C\smile_2\delta C+\omega\smile C,
$$

$$
g_4=\begin{cases}
\delta_s\Bigl[\Omega_4+\mathsf h(t_4\smile_3\tau'_1)\Bigr]-\mathsf h\bigl(t_4\smile_2t_4+\omega\smile t_4\bigr),\quad t_4=\delta C+\tau'_1(A;B), & (A,B)\in L_4,\\[4pt]
-\,I\,\mathcal K_5\bigl(p_4(A\ell,B\ell),\ \delta(C\ell)+H_4(A\ell,B\ell)\bigr), & \text{otherwise.}
\end{cases}
$$

## Degree 5

State $A\in C^2(\mathbb Z_s)$, $B\in C^3$, $C\in C^4$, $D\in C^6$.

$$
f_5(A,B)=\begin{cases}\tau'_2(A;B), & (A,B)\in L_5,\\
H_5(A,B)=I\,\tau'_3\bigl(\delta_s(A\ell);\ \delta(B\ell)+P(A\ell)\bigr), & \text{otherwise,}\end{cases}
$$

$$
\Omega_5=\mathsf h(E(C))+\widehat R_2(A,B)+\tfrac23\,\mathrm L\bigl(\rho_3(A\smile_sA\smile_sA)\bigr),\qquad
E(C)=C\smile_2C+C\smile_3\delta C+\omega\smile C .
$$

$\widehat R_2$ contains the term $-\tfrac14A^{\smile3}$ exactly once; the
three-primary term is $\tfrac23$ times the $\{0,1,2\}$ lift of the signed
cube modulo three. Its coboundary is integral and exact on cocycles (so it
does not change $d_5$), but the term is kept because its failure of
additivity is the three-primary part of $\gamma_5$.

$$
g_5=\begin{cases}
\delta_s\Bigl[\Omega_5+\mathsf h(t_5\smile_4\tau'_2)\Bigr]-\mathsf h\bigl(t_5\smile_3t_5+\omega\smile t_5\bigr),\quad t_5=\delta C+\tau'_2(A;B), & (A,B)\in L_5,\\[4pt]
-\,I\,\mathcal K_6\bigl(p_5(A\ell,B\ell),\ \delta(C\ell)+H_5(A\ell,B\ell)\bigr), & \text{otherwise.}
\end{cases}
$$

## Degree 6

State $A\in C^3(\mathbb Z_s)$, $B\in C^4$, $C\in C^5$, $D\in C^7$.

$$
f_6(A,B)=\begin{cases}\tau'_3(A;B), & (A,B)\in L_6,\\
H_6(A,B)=I\,\tau'_4\bigl(\delta_s(A\ell);\ \delta(B\ell)+P(A\ell)\bigr), & \text{otherwise,}\end{cases}
$$

$$
\Omega_6=\mathsf h(E(C))+\widehat R_3(A,B)+\tfrac23\,\mathrm L\bigl(P^1_s\rho_3A\bigr),\qquad
E(C)=C\smile_3C+C\smile_4\delta C+\omega\smile C,
$$

with $P^1_s$ the nineteen-term reduced power; here
$\tfrac23\delta_s\mathrm L(P^1_s\rho_3A)=2\beta_{3,s}P^1_s\rho_3A$ is the
three-primary $d_5$ term.

$$
g_6=\begin{cases}
\delta_s\Bigl[\Omega_6+\mathsf h(t_6\smile_5\tau'_3)\Bigr]-\mathsf h\bigl(t_6\smile_4t_6+\omega\smile t_6\bigr),\quad t_6=\delta C+\tau'_3(A;B), & (A,B)\in L_6,\\[4pt]
T_3\bigl(\mathcal R_6(A,B,C)\bigr), & \text{otherwise,}
\end{cases}
$$

where $T_3=\mathcal K_6=\delta_s\Omega_6$ on full defining systems and
$\mathcal R_6$ is the section retraction of
[degree_six_section.md](degree_six_section.md), which sends any triple to
a full defining system. There is no $T_4$ to integrate, and $g_6$ must be
nonzero outside $L_6$ (with $g_6=0$ there no normalized product satisfies
$\mathfrak d(x\times y)=\mathfrak dx\times\mathfrak dy$).

## Degree 7 (endpoint only)

To state $\mathfrak d_7\mathfrak d_6=0$ use the same shape with

$$
f_7=\begin{cases}\tau'_4(A;B),&(A,B)\in L_7,\\H_7(A,B)=I\,\tau'_5(p_7(A\ell,B\ell)),&\text{otherwise},\end{cases}
\qquad
g_7=\begin{cases}\tfrac12\bigl(\delta_s\widetilde{E(C)}-\widetilde{E(\delta C)}\bigr),&A=B=0,\\0,&\text{otherwise}.\end{cases}
$$

## The pure-C rule and the telescoping form

For $A=B=0$ in every degree the D correction is, with $E=E$ of the degree,

$$
\boxed{\;
\mathfrak d_k(0,0,C,D)=\Bigl(0,\ 0,\ \delta C,\ \delta_sD+\tfrac12\bigl(\delta_s\widetilde{E(C)}-\widetilde{E(\delta C)}\bigr)\Bigr).\;}
\tag{P}
$$

This rule is natural: it can be evaluated with any Steenrod cup-$i$
structure on a resolution ([natural_evaluation.md](natural_evaluation.md), Section 3.2).

Every $g_k$ above is $g_k(u)=\delta_s\Phi_k(u)-\Phi_{k+1}(F_ku)$ for the
pointed potentials ($F_k$ the first three components)

$$
\Phi_k=\begin{cases}
\Omega_k+\mathsf h(t_k\smile_{k-1}\tau_k),&(A,B)\in L_k,\\
I\,\Omega_{k+1}\bigl(p_k(A\ell,B\ell),\ \delta(C\ell)+H_k(A\ell,B\ell)\bigr),&(A,B)\notin L_k,\ 2\le k\le5,\\
\Omega_6(\mathcal R_6(A,B,C)),&k=6,\ (A,B)\notin L_6,
\end{cases}
\qquad \Phi_7=\begin{cases}\mathsf h(E(C)),&A=B=0,\\0,&\text{else},\end{cases}
$$

($\Omega_k$ normalized by $\Omega_k(0)=0$). Hence
$\mathfrak d_{k+1}\mathfrak d_k=0$ on all cochains for $0\le k\le6$: the
B component by (C7), the C component by $\delta\tau_k=0$ on $L_k$ and
$\delta H_k=\tau_{k+1}(p_k)$ off it, the D component by
$\delta_sJ_k=-\delta_s\mathsf h(E(t_k))=-g_{k+1}(0,0,t_k)$ on $L_k$ and by
$\delta_sG_k=-\mathcal K_{k+1}(F_ku)$ off it. The integrality of each
$g_k$, not the cancellation, is the content: an integral $\Phi_k$ would
kill the tertiary obstruction class.

## Relation to the E6 page

For $3\le k\le6$ and a full defining system, $\mathcal K_k=T_{k-3}(A;B,C)$
represents $d_5$; for $k=2$, $\mathcal K_2=\psi'_0$ represents $d_4$ in
input degree zero. The E6 layers are

$$
\begin{aligned}
A_k&=\{a\in H^{k-3}(\mathbb Z_s):\ \operatorname{Dbar}a=0,\ \operatorname{Tau}a=0,\ Ta=0\},\\
B_k&=\bigl(\ker D\cap\ker\operatorname{Psi}\subset H^{k-2}(\mathbb F_2)\bigr)/\operatorname{Dbar}H^{k-4}(\mathbb Z_s),\\
C_k&=\ker\bigl(\operatorname{Dtilde}:H^{k-1}(\mathbb F_2)\to H^{k+2}(\mathbb Z_s)\bigr)/\bigl(DH^{k-3}+\operatorname{Tau}\ker\operatorname{Dbar}\bigr),\\
D_k&=H^{k+1}(\mathbb Z_s)/\bigl(\operatorname{Dtilde}H^{k-2}+\operatorname{Psi}_{k-3}(\ker D)+T_{k-4}(E_5^{k-4,0})\bigr).
\end{aligned}
$$
