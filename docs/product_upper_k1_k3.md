# The D-layer correction $\gamma_k$ in degrees $\le3$

$$
x\times_ky=\bigl(\mu_k(u,v),\ D+D'+\gamma_k(u,v)\bigr),\qquad \gamma_k\in C^{k+1}(X;\mathbb Z_s).
$$

The mechanism in every degree: with $\Phi_k$ the potential of the
differential and $e_k(u,v)\in C^k(X;\mathbb Q_s)$ a *pair phase* with
$e_k(0,0)=0$,

$$
\gamma_k=\Delta\Phi_k+\delta_se_k(u,v)+e_{k+1}(F_ku,F_kv)
\tag{$\gamma$}
$$

satisfies the upper identity of Section 6 by cancellation; the content of
each degree is a pair phase that makes $\gamma_k$ **integral**, and the
rule that the successor $e_{k+1}$ used in degree $k$ is literally the phase
of degree $k+1$. Changing $e_k$ by a $\mathbb Q/\mathbb Z$-cocycle changes
$\gamma_k$ by the Bockstein of that cocycle, which is not a coboundary in
general: the specific phases below are part of the formula.
Notation: [conventions.md](conventions.md); lower products:
[product_lower.md](product_lower.md).

## 1. Degrees $-1$ and $0$

$\gamma=0$: $D$ layers add.

## 2. Degree 1

States $(C,D)$ with $C\in C^0(\mathbb F_2)$. $\Phi_1=\mathsf h(\omega\smile C)$,
$E(C)=\omega\smile C$, $F_1C=(0,\delta C)$,

$$
e_1(C,C')=\mathsf h(C\smile\delta C'),\qquad e_1(v_0v_1)=\tfrac12\,C(v_0)\bigl[(C'(v_0)+C'(v_1))\bmod2\bigr],
$$

$$
\gamma_1=\tfrac12\bigl(\widetilde{E(C)}+\widetilde{E(C')}-\widetilde{E(C+C')}\bigr)+\delta_se_1(C,C')+e_2\bigl(0,\delta C;\,0,\delta C'\bigr),
$$

with $e_2$ the degree-two phase of Section 3. Everything collapses to

$$
\boxed{\;
\gamma_1(C,C')=\widetilde{\omega\smile(C\cdot C')}+\delta_s\mathsf h(C\smile\delta C')-\mathsf h(\delta C\smile\delta C'),\;}
$$

($C\cdot C'$ the pointwise product of the two $0$-cochains; the last two
terms are integral together). On closed $C,C'$, $\gamma_1=\widetilde{\omega\smile(C\cdot C')}$.

## 3. Degree 2

States $(B,C,D)$, $B\in C^0$, $C\in C^1$. $F_2(B,C)=(\delta B,\ \delta C+Q_D(B))$; every
such image is a full defining system of degree three with $A=0$.

$$
\boxed{\;
\gamma_2=\Phi_2(B,C)+\Phi_2(B',C')-\Phi_2\bigl(B+B',\,C+C'+\beta_2(B,B')\bigr)
+\delta_se_2(B,C;B',C')+e_3\bigl(F_2(B,C);\,F_2(B',C')\bigr),\;}
$$

$$
\Phi_2(B,C)=\begin{cases}
\mathsf h(C\smile\delta C+\omega\smile C)+\mathsf h\bigl(t_2\smile_1(\omega\smile B)\bigr),\quad t_2=\delta C+\omega\smile B, & \delta B=0\text{ on }X,\\[3pt]
I\,\Omega_3\bigl(0,\ \delta(B\ell),\ \delta(C\ell)+Q_D(B\ell)\bigr), & \text{otherwise},
\end{cases}
$$

(each of the three potentials takes its own branch) with
$\Omega_3(0,b,c)=\mathsf h(E(c))+\theta_1(b)$ ([differential.md](differential.md), degree 3).
The phase $e_2$ is the prism of the degree-three phase with a C-transport:

$$
\boxed{\;
e_2(B,C;B',C')=-I\Bigl[e_3\bigl(b_I,c_I;\,b'_I,c'_I\bigr)
+\mathsf h\bigl(E(\eta)+c^{\rm tot}_I\smile_1\delta\eta+\delta c^{\rm tot}_I\smile_2\delta\eta\bigr)\Bigr],\;}
$$

on $X\times I$ with $(b_I,c_I)=F_2(B\ell,C\ell)$, $(b'_I,c'_I)=F_2(B'\ell,C'\ell)$,
$\eta=\beta_2(B\ell,B'\ell)+\beta_2(B,B')\ell$,
$c^{\rm tot}_I=c_I+c'_I+\beta_3(b_I,b'_I)$, $\beta_3(b,b')=b\smile b'+s\smile(b\smile_1b')$
(the legal $\beta^L_3$ at $A=A'=0$), and $E(\eta)=\eta\smile\delta\eta+\omega\smile\eta$.

## 4. Degree 3 with $A=A'=0$: the paper-calibrated phase $e_3$

For full defining systems $(0,b,c)$, $(0,b',c')$ ($b$ closed, $\delta c=Q_D(b)$),
write $h=b\smile_1b'$, $b_S=b+b'$, $m=\beta_3(b,b')=b\smile b'+s\smile h$,
$c_S=c+c'+m$. The transcription of the degree-three exchange phase of
arXiv:2310.19058 into these cup conventions is $e^{\rm paper}_3=\tfrac12\widetilde\epsilon+\vartheta_0+\vartheta_\omega+\vartheta_s+\vartheta_{\rm mix}$ with

$$
\begin{aligned}
\epsilon&=c\smile_1c'+\delta c\smile_2c+c'\smile_2\delta c+\delta c'\smile_2c'+(\delta c+\delta c')\smile_2c_S+m\smile_1c_S+m\smile_2\delta c_S,\\
\vartheta_0&=\tfrac12\,\widetilde{(b\smile h)\smile b'}+\tfrac14\,\widetilde{(b\smile b')\smile b'},\qquad
\vartheta_\omega=\tfrac12\,\widetilde{(\omega\smile b_S)\smile_2(b\smile b')}+\tfrac34\,\widetilde{\omega\smile h},\\
\vartheta_s&=\tfrac12\Bigl((s\smile s)\smile h+\bigl(s\smile((s\smile_1b)\smile_1b')\bigr)\smile h+(b\smile(s\smile_1b'))\smile h+((s\smile_1b)\smile b_S)\smile h+((s\smile_1b)\smile h)\smile b+(s\smile h)\smile b'\Bigr)^{\sim}\\
&\quad+\tfrac14\,\widetilde{(s\smile b')\smile b}+\tfrac14\,\widetilde{(s\smile h)\smile b_S}+\tfrac58\,\widetilde{(s\smile b)\smile b'},\\
\vartheta_{\rm mix}&=\tfrac12\Bigl((\omega\smile_1s)\smile h+(\omega\smile_2(s\smile h))\smile b_S+(\omega\smile_2(s\smile b'))\smile b+(\omega\smile_2(s\smile b))\smile(b'+h)\Bigr)^{\sim}.
\end{aligned}
$$

For a generic $s$ this transcription fails the cochain consistency
equations; the repair, solved once over $\mathbb Z/4$ from the universal
consistency equations and zero when $s$, $b$ or $b'$ vanishes, is the
degree-three cochain, on $(v_0v_1v_2v_3)$ with $x_i=b(v_{i-1}v_i)$,
$y_i=b'(v_{i-1}v_i)$, $t_1=s(v_0v_1)$, $t_2=s(v_1v_2)$ (integer sums of
$0/1$ values, not reduced):

$$
\chi^{\rm rep}(b,b')=\tfrac14\Bigl(t_2x_1y_3+t_1(x_2+x_3)y_2y_3\Bigr)
+\tfrac12\Bigl(t_2y_3\bigl(y_1(x_1+x_2+x_3)+x_1(x_2+x_3+y_2+t_1)\bigr)+t_1x_3y_2(y_1+x_1)\Bigr).
$$

With the coordinate gauge (lifts separate; $\widetilde b^{\,3}=\widetilde b\smile\widetilde b\smile\widetilde b$ integral)

$$
\kappa(b,c)=\tfrac18\,\widetilde b^{\,3}-\mathsf h\bigl(s\smile c+c\smile_2\delta c\bigr)
+\mathsf h\Bigl((b+s)\smile c+(b\smile_1s)\smile(b\smile b)+(b\smile_1\omega)\smile b\Bigr),
$$

the selected legal phase and product at $A=A'=0$ are

$$
\boxed{\;
e^{\rm pap}_3(b,c;b',c')=e^{\rm paper}_3+\chi^{\rm rep}(b,b')+\kappa(b_S,c_S)-\kappa(b,c)-\kappa(b',c')+\mathsf h(\omega\smile h),
\qquad
\gamma_3=\Omega_3(0,b,c)+\Omega_3(0,b',c')-\Omega_3(0,b_S,c_S)+\delta_se^{\rm pap}_3 .\;}
$$

This $e^{\rm pap}_3$ is the successor used by $e_2$ and $e_1$. (Closed
forms: $e^{\rm pap}_3(0,c;0,c')=\mathsf h(c\smile_1c')$ for closed $c,c'$;
$e_2(0,a;0,a')=-\mathsf h(a\smile a')$ for closed $a,a'$ of degree one.)

The older "production" phase at $A=A'=0$, still used inside the even-$A$ phase of Section 5, is

$$
e^{\rm prod}_3(b,c;b',c')=e^{\rm paper}_3+\chi^{\rm rep}-\kappa_c(c_S)+\kappa_c(c)+\kappa_c(c')+\tfrac18\bigl(\widetilde b_S^{\,3}-\widetilde b^{\,3}-\widetilde b'^{\,3}\bigr)+\sigma(b,b'),
$$

$$
\kappa_c(c)=\mathsf h(s\smile c+c\smile_2\delta c),\qquad
\sigma(b,b')=\tfrac14\,\widetilde s\smile(\widetilde b\smile\widetilde b')+\mathsf h\Bigl(s\smile\bigl(h\smile b+b'\smile h\bigr)+(s\smile s)\smile h\Bigr).
$$

## 5. Degree 3 with legal AB and arbitrary $C$

Let $u=(A,B,C)$, $v=(A',B',C')$ with $(A,B),(A',B')\in L_3$ (so $|A|$ and
$|A'|$ are locally constant), $\tau=\tau'_0(A;B)$, $t=\delta C+\tau$,
$\Phi_3(u)=\Omega_3(u)+\mathsf h(t\smile_2\tau)$, lower product $u\cdot v$ with
$\alpha_3,\beta^L_3$. The **closed-AB product** is

$$
\boxed{\;
\Gamma^{\rm cl}_3(u,v)=\Phi_3(u)+\Phi_3(v)-\Phi_3(u\cdot v)+\mathsf h\bigl(t\smile_2t'\bigr)+\delta_se^{\rm cl}_3(u,v),
\qquad
\gamma_3=\Gamma^{\rm cl}_3(u,v)-\Bigl[\mathsf h(s\smile t)+\mathsf h(s\smile t')-\mathsf h(s\smile t_S)\Bigr]+[\,u,v\text{ full legal}\,]\,\delta_s\chi^{\rm cal}_3,\;}
$$

where the bracket is an integral carry ($t_S=t+t'$ on legal AB), the
pure-C successor term is $\mathsf h(h^2(t,t'))=\mathsf h(t\smile_2t')$, and
$e^{\rm cl}_3$ is chosen **componentwise by the parity of $A,A'$ at the
first vertex**: the even phase where both are even, the odd phase where
at least one is odd.

**Calibration.** When both inputs are full defining systems ($F_3u=F_3v=0$ on $X$),
with $m(A)=\bigl((A^2-1)/8\bigr)\bmod2$ at the first vertex,

$$
\chi^{\rm cal}_3=\begin{cases}
e^{\rm pap}_3(B,C;B',C')-e^{\rm cl}_3(u,v), & A=A'=0,\\
\tfrac12\,m(A)\,\widetilde{(B+B')\smile(B+B')\smile(B+B')}, & A,A'\text{ both odd},\\
0, & \text{otherwise}.
\end{cases}
$$

### 5.1 The even phase

For even $A$ define the canonical data ($K=A/2$, $t_A=\rho K$,
$\kappa_A=(K-\widetilde{t_A})/2$, $r=\rho(\delta\widetilde\omega/2)$, $r_c=\rho\bigl((\delta\widetilde\omega/2-\widetilde r)/2\bigr)$):

$$
Z=t_A\smile\omega,\quad \lambda=\rho\bigl(\kappa_A\smile\widetilde\omega\bigr)+Z\cdot(s\smile B),\quad H=Q_D(B),\quad z=t_A\smile r,\quad
\bar C=C+\lambda,\quad \mathsf g=\mathsf h\bigl(h^2(C,\lambda)+Q^1(\lambda)\bigr),
$$

and the same with primes for $v$ and with subscript $S$ for the total
$(A+A',\ B+B',\ C+C'+\beta^L_3)$. With $m_B=\alpha_4(B,B')=h^D(B,B')$ and
$e^{\rm prod}_3$ of Section 4,

$$
\begin{aligned}
e^{\rm ev}_3(u,v)={}&\mathsf h\bigl(h^2(\bar C,\bar C')+h^2(\bar C+\bar C',m_B)\bigr)
+e^{\rm prod}_3(B,0;B',0)+\mathsf h\bigl((H+H')\smile_2m_B+H\smile_3H'\bigr)\\
&+\mathsf h\bigl(z\smile_3H'+(z+z')\smile_3\delta m_B+(t_A\smile t_{A'})\smile r_c\bigr)
+\mathsf g_S-\mathsf g-\mathsf g'
+\mathsf h\bigl((H+z)\smile_3t'\bigr),
\end{aligned}
$$

where the last term carries the C residual of the second input and
$e^{\rm cl}_3=e^{\rm ev}_3$ on even components. (Here $\Omega_3=\Omega_3^{\rm ev}$ of
[differential.md](differential.md).)

### 5.2 The odd phase

On a component where $A$ or $A'$ is odd, let the *anchor* be $A$ if it is
odd, else $A'$; $\ell=[\text{anchor}<0]$; $U$ the degree-one cochain equal
to $B$ where $A$ is odd and to $B'$ elsewhere (on legal data $\delta U=\omega$
there). Two legal cylinders transport the data to a pure-C situation:

- *Orientation.* On $X\times I$ with gauge $\gamma=\ell\,\ell_I$ (ramp of $\ell$):
  $s_I=\pi^*s+\delta\gamma$, $A_I=(-1)^{\gamma}\pi^*A$, $B_I=\pi^*B$,
  $C_I=\pi^*C+\rho\bigl(\mathcal J\,\tau'_0(A_I;B_I)\bigr)$, where $\mathcal J$
  is the contraction of $X\times I$ onto $X\times\{0\}$: $\mathcal Jc$ is the
  prism of the cochain $((v,\epsilon),t)\mapsto c((v,\epsilon t))$ over a
  second interval. $u_I=(A_I,B_I,C_I)$ is a full defining system for the
  twists $(s_I,\pi^*\omega)$, equal to $u$ at $t=0$ and to the oriented
  triple $u^{(1)}$ (with $A^{(1)}=|A|$ and twist $s+\delta\ell=0$) at $t=1$;
  the same for $v$.
- *Flattening.* For a triple $x=(A,B,C)$ with twists $(0,\omega)$ and a
  $U$ with $\delta U=\omega$: $\bar b=B+a\smile U$,
  $\bar c=C+\Lambda(A,B,U)$, with
  $\Lambda(A,B,U)=U\smile B+(a\smile\bar b)\smile U+t(A)\cdot\rho\bigl((\delta\widetilde U-\widetilde\omega)/2\bigr)$;
  on $X\times I$: $u_I=\pi^*U+U\ell$, $w_I=\delta u_I$, $A_I=\pi^*A$,
  $B_I=\pi^*\bar b+a\smile u_I$, $C_I=\pi^*\bar c+\Lambda(A_I,B_I,u_I)$ (twists $(0,w_I)$).
  At $t=1$ the twist $\omega$ has been removed.

With $\mathrm{src}(x,y)=\Phi_3(x)+\Phi_3(y)-\Phi_3(x\cdot y)+\mathsf h(t_x\smile_2t_y)$
(evaluated with the twists of the cylinder in question), $u^{(1)},v^{(1)}$
the $t=1$ ends of the orientation cylinders, $u^{\flat},v^{\flat}$ the
$t=1$ ends of their flattening cylinders, and $u^{\flat}_0=(0,B^\flat,C^\flat)$
(the odd rank replaced by zero, twists $(0,0)$):

$$
\boxed{\;
e^{\rm od}_3(u,v)=I\,\mathrm{src}(u_I,v_I)+(-1)^{\ell}\Bigl[I\,\mathrm{src}(u^{\flat}_I,v^{\flat}_I)+e^{\rm ev}_3\bigl(u^\flat_0,v^\flat_0\bigr)\Bigr],\;}
$$

where the first prism is over the orientation cylinders and the second
over the flattening cylinders.

## 6. Degree 3 with a non-legal lower pair

Let $\widehat F_3$, $\phi$ be as in [product_lower.md](product_lower.md) Section 5
(with $q_{\rm loc}=0$, so $\phi=\mathrm{id}$ in degree three) and $\widehat\mu_3$
the lower product with $\beta^{\rm raw}_3$. Natural off-legal data by prisms
of degree four:

$$
\widehat g_3(u)=-I\,\mathcal K_4\bigl(\widehat F_3(u\ell)\bigr),\qquad
\widehat\gamma_3(u,v)=I\Bigl[\widehat g_3(u\ell)+\widehat g_3(v\ell)-\widehat g_3\bigl(\widehat\mu_3(u\ell,v\ell)\bigr)
+\Gamma^{\rm norm}_4\bigl(\widehat F_3u\ell,\widehat F_3v\ell\bigr)\Bigr],
$$

with $\Gamma^{\rm norm}_4$ the normalized legal product of degree four
([product_upper_k4_k6.md](product_upper_k4_k6.md)). The comparison with the
prescribed $g_3$ is the integral cochain

$$
Z_3(u)=\begin{cases}\Phi_3(u)-\widehat\Phi_3(u)-\mathsf h(s\smile t)-\delta_sq_3(u),&(A,B)\in L_3,\\0,&\text{otherwise},\end{cases}
\qquad \widehat\Phi_3(u)=I\,\Omega_4\bigl(\widehat F_3(u\ell)\bigr),
$$

and

$$
\boxed{\;\gamma_3(u,v)=\widehat\gamma_3(u,v)+Z_3(u)+Z_3(v)-Z_3(u\cdot v)\qquad\text{when $u,v$ are not both lower-legal.}\;}
$$

**The gauge $q_3$.** Put $\Psi(u)=\Phi_3(u)-\widehat\Phi_3(u)-\mathsf h(s\smile t)$,
a closed phase on legal AB. $q_3$ is selected at each vertex by the rank
$A(v_0)$, with $\ell_A=[A<0]$, $\gamma=\ell_A\ell$, $A_I=(-1)^\gamma\pi^*A$,
$s_I=\pi^*s+\delta\gamma$ and the orientation prism
$O(u)=I\,\Psi(A_I,\pi^*B,\pi^*C)$ (twists $s_I,\pi^*\omega$):

| rank at $v_0$ | $q_3(u)$ |
| --- | --- |
| $0$ | $I\,\Psi(\pi^*B,C\ell)\;+\;\mathfrak q_3^{(0)}(B)$, the universal primitive of the source $\Psi(b,0)$ on the model $(b,(s,\omega))$ ([universal_primitives.md](universal_primitives.md), Section 6) |
| even, nonzero | $(-1)^{\ell_A}\Bigl[I\,\Psi\bigl(\pi^*|A|,\pi^*B,C\ell\bigr)_{s=0}+\mathfrak q_3^{\rm ev}(|A|,B;\omega)\Bigr]-O(u)$ |
| odd | $(-1)^{\ell_A}\Bigl[-I\,\Psi\bigl(\pi^*|A|,\ \pi^*B+B\ell,\ \pi^*C+C\ell\bigr)_{s=0,\ \omega=\delta(\pi^*B+B\ell)}\Bigr]-O(u)$ |

where $\mathfrak q_3^{\rm ev}$ is the universal primitive of the rank-$r$ source
$\Psi(r,b,0)_{s=0}$ on the model $(b,\omega)$ with the fixed small
coefficients of the even certificate (same section). On legal AB each
selector is constant on connected components.

## 7. The upper identity (all degrees)

$$
\boxed{\;\delta_s\gamma_k(u,v)=g_k(u)+g_k(v)-g_k\bigl(\mu_k(u,v)\bigr)+\gamma_{k+1}\bigl(F_ku,F_kv\bigr),\;}
$$

which, with the lower identity, gives $\mathfrak d_k(x\times_ky)=\mathfrak d_k(x)\times_{k+1}\mathfrak d_k(y)$
for arbitrary integral $D,D'$; products of flat states are flat.

**Gauges and the group of a degree.** A gauge is a state $g$ of degree
$k-1$ acting by $\operatorname{act}(g,x)=\mathfrak d_{k-1}(g)\times_kx$; two flat
states are equivalent when a gauge carries one to the other; the group of
degree $k$ is flat states modulo gauge. Commutativity and associativity on
gauge classes are assumed (on legal data the exchange of the two factors
is a coboundary: $k\le1$ literally symmetric; $k=2$:
$\gamma_2(v,u)-\gamma_2(u,v)=\delta_sM$ with $M=e^{\rm op}_2-e_2-\delta_s\mathsf h(C\smile_1C'+C\smile B')$).
