# The stacking product: shape and the lower corrections $\alpha_k,\beta_k$

The product of two states of package degree $k$ is triangular:

$$
\boxed{\;
(A,B,C,D)\times_k(A',B',C',D')=\bigl(A+A',\ \ B+B'+\alpha_k(A,A'),\ \ C+C'+\beta_k(A,B;A',B'),\ \ D+D'+\gamma_k(A,B,C;A',B',C')\bigr).\;}
$$

Binary layers add modulo two. $\alpha_k$ depends only on $A,A'$; $\beta_k$
only on $A,B,A',B'$; $\gamma_k$ only on the lower triples. The zero state is
a two-sided unit and every correction vanishes when one entire input is
zero. This page gives $\alpha_k$ and $\beta_k$ in every degree; $\gamma_k$ is
in [product_upper_k1_k3.md](product_upper_k1_k3.md) and
[product_upper_k4_k6.md](product_upper_k4_k6.md). Write $n=k-3$,
$a=\rho A$, $a'=\rho A'$, and for the total $A_S=A+A'$,
$B_S=B+B'+\alpha_k(A,A')$. Notation: [conventions.md](conventions.md).

## 1. The B-layer correction $\alpha_k$

$$
\boxed{\;\alpha_k(A,A')=h^D(a,a')=a\smile_{n-1}a'+\delta a\smile_na'+s\smile\bigl(a\smile_na'+\delta a\smile_{n+1}a'\bigr),\qquad \alpha_k=0\ (k\le2).\;}
$$

| $k$ | $n$ | $\alpha_k(A,A')$ for arbitrary $A,A'$ | for signed cocycles ($\delta a=0$) |
| ---: | ---: | --- | --- |
| 3 | 0 | $(\delta a)\,a'+s\smile\bigl(a\,a'+\delta a\smile_1a'\bigr)$ | $s\smile(a\,a')$ |
| 4 | 1 | $a\smile a'+\delta a\smile_1a'+s\smile(a\smile_1a'+\delta a\smile_2a')$ | $a\smile a'+s\smile(a\smile_1a')$ |
| 5 | 2 | $a\smile_1a'+\delta a\smile_2a'+s\smile(a\smile_2a'+\delta a\smile_3a')$ | $a\smile_1a'+s\smile(a\smile_2a')$ |
| 6 | 3 | $a\smile_2a'+\delta a\smile_3a'+s\smile(a\smile_3a'+\delta a\smile_4a')$ | $a\smile_2a'+s\smile(a\smile_3a')$ |
| 7 | 4 | $a\smile_3a'+\delta a\smile_4a'+s\smile(a\smile_4a'+\delta a\smile_5a')$ | $a\smile_3a'+s\smile(a\smile_4a')$ |

Identity (C10): $P(A+A')=P(A)+P(A')+\delta\alpha_k(A,A')+\alpha_{k+1}(\delta_sA,\delta_sA')$.
$\alpha_k$ and this identity are natural for any cup-$i$ structure on a
resolution ([natural_evaluation.md](natural_evaluation.md), Section 3.3);
$\beta_k$ for $k\ge3$ is not (Section 4 there).

## 2. The C-layer correction $\beta_k$, degrees $\le2$

$$
\beta_{-1}=\beta_0=\beta_1=0,\qquad
\boxed{\;\beta_2(B,B')=\delta B\smile B'+s\smile(B\smile B').\;}
$$

## 3. The legal C-layer correction $\beta^L_k$, degrees $3$ to $7$

For $(A,B)$ and $(A',B')$ both in $L_k$ ($\delta_sA=0$, $\delta B=P(A)$, and the same for the primed pair):

$$
\boxed{\;
\beta^L_k=h^D(B,B')+h^D\bigl(B+B',\alpha_k(A,A')\bigr)
+\lambda(A,B)+\lambda(A',B')+\lambda(A_S,B_S)+\mathcal V_n(A,A')+[k=5]\,a\smile a' .\;}
$$

- $h^D$ is (C10) in degree $n+1$: $h^D(B,B')=B\smile_nB'+\delta B\smile_{n+1}B'+s\smile(B\smile_{n+1}B'+\delta B\smile_{n+2}B')$.
- $\lambda(A,B)=g(A)\cdot(s\smile B)$ is the splitting carry (S23) of [secondary_operations.md](secondary_operations.md), Section 6, with the degree-zero value (C13) (pointwise product of the two binary cochains of degree $n+2$); at $k=3$, $g=g^\flat(A)=\omega\smile t(A)+a\,(s\smile s)$.
- $\mathcal V_n(A,A')$ is the **universal pair primitive**: a binary cochain of degree $n+2$ with $\mathcal V_n(A,0)=\mathcal V_n(0,A')=0$ and
  $$\delta\mathcal V_n(A,A')=\mathcal S^\beta_n(A,A'):=k_0(A_S)+k_0(A)+k_0(A')+Q_D(\alpha_k)+h^D\bigl(P(A),P(A')\bigr)+h^D\bigl(P(A)+P(A'),\delta\alpha_k\bigr).$$
  For $n=0$ it is the explicit polynomial below; for $n=1,2,3,4$ it is produced by the fixed universal contraction described in [universal_primitives.md](universal_primitives.md), Section 3 (no equation is solved on $X$).
- The closed term $a\smile a'$ at $k=5$ is fixed by the periods of the upper source; it vanishes when $A=0$ or $A'=0$.

**The polynomial $\mathcal V_0$.** On a face $(v_0v_1v_2)$ write
$A(v_0)\equiv a_0+2a_1$ and $A'(v_0)\equiv b_0+2b_1\pmod4$ (so $a_0=a$, $a_1=t(A)$),
$x=s(v_0v_1)$, $y=s(v_1v_2)$, $w=\omega(v_0v_1v_2)$:

$$
\mathcal V_0(A,A')(v_0v_1v_2)=a_0b_1\,xy+w\bigl(a_1b_1+a_0b_0(x+y)+(a_1b_0+a_0b_1)\,xy\bigr)\pmod2 .
$$

Per degree:

| $k$ | $n$ | $\beta^L_k$ |
| ---: | ---: | --- |
| 3 | 0 | $h^D(B,B')+h^D(B+B',\alpha_3)+\lambda(A,B)+\lambda(A',B')+\lambda(A_S,B_S)+\mathcal V_0(A,A')$, with $g^\flat$ and the polynomial above |
| 4 | 1 | the same shape with $g(A)$, $\mathcal V_1$ |
| 5 | 2 | the same shape with $\mathcal V_2$, plus $a\smile a'$ |
| 6 | 3 | the same shape with $\mathcal V_3$ |
| 7 | 4 | the same shape with $\mathcal V_4$ (successor only) |

Property: $\delta\beta^L_k=\tau_k(A,B)+\tau_k(A',B')+\tau_k(A_S,B_S)$, and
$\beta^L_k$ does not depend on $C$.

## 4. Extension to closed $A$ and arbitrary $B$ (used inside the upper corrections)

For signed cocycles $A,A'$ and arbitrary binary $B,B'$, with
$b=\delta B+P(A)$, $b'=\delta B'+P(A')$ (degree $n+2$):

$$
\beta^\sharp(A,B;A',B')=\beta^L_k(A,B;A',B')+b\smile_{n+2}\delta B'+\delta\alpha_k\smile_{n+2}(b+b'),
$$

where $\beta^L_k$ is the formula of Section 3 evaluated on the given data.
On legal pairs $b=b'=0$ and $\beta^\sharp=\beta^L_k$. It is the lower
correction for the extension $f^\sharp$ of $\tau'$ to arbitrary $B$
([product_upper_k4_k6.md](product_upper_k4_k6.md), Section 1).

## 5. Arbitrary inputs, degrees $3$ to $6$

**The comparison gauge $\xi_k$.** On legal pairs, $\tau_k(A,B)$ and the
prism $H_k(A,B)$ differ by a coboundary, $\tau_k+H_k=\delta\xi_k$. Fix

$$
\xi_k(A,B)=\lambda(A,B)+I\Bigl[Q_D(B\ell)+h^D\bigl(\delta(B\ell),P(A\ell)\bigr)+\lambda\bigl(\delta_s(A\ell),\ \delta(B\ell)+P(A\ell)\bigr)\Bigr]+v_n(A)\qquad(n\ge1),\qquad \xi_3=0,
$$

where the bracket is on $X\times I$ (both inputs ramped) and $v_n(A)$ is
the universal primitive of the closed binary source

$$
Z_n(A)=k_0(A)+I\Bigl[k_0\bigl(\delta_s(A\ell)\bigr)+Q_D\bigl(P(A\ell)\bigr)\Bigr]+h^D\bigl(P(A),P(A)\bigr)
$$

([universal_primitives.md](universal_primitives.md), Section 4; its small
restriction vanishes identically, so there are no coefficients).

**Local extension and natural coordinates.** Let $q_{\rm loc}(A,B)$ be the
cochain equal to $\xi_k(A,B)$ on every simplex all of whose faces have
vanishing curvature $p_k(A,B)$, and zero on the other simplices. Put

$$
\widehat f_k=H_k+\delta q_{\rm loc},\qquad
e(A,B)=\begin{cases}0,&(A,B)\in L_k,\\ q_{\rm loc}(A,B),&\text{otherwise},\end{cases}\qquad
\phi(A,B,C)=(A,B,\ C+e(A,B)),
$$

so that $\widehat f_k=\tau_k$ on $L_k$ and $\widehat F_k=(\delta_sA,\ \delta B+P(A),\ \delta C+\widehat f_k)$
satisfies $F_k=\widehat F_k\circ\phi$ ($F_k$ the lower three components of $\mathfrak d_k$).

**The raw correction** is the prism of the closed pair residual, on the
ramped inputs $A\ell,B\ell,A'\ell,B'\ell$ (with $A_S=A\ell+A'\ell$,
$B_S=B\ell+B'\ell+\alpha_k(A\ell,A'\ell)$ on the cylinder):

$$
\beta^{\rm raw}_k(u,v)=I\Bigl[\widehat f_k(A_S,B_S)+\widehat f_k(A,B)+\widehat f_k(A',B')
+\beta^L_{k+1}\bigl(p_k(A,B),\,p_k(A',B')\bigr)\Bigr].
$$

The successor $\beta^L_{k+1}$ is evaluated on curvature pairs, which are
legal because $p_{k+1}p_k=0$ (for $k=6$ this is $\beta^L_7$).

**The complete correction:**

$$
\boxed{\;
\beta_k(u,v)=\begin{cases}
\beta^L_k(A,B;A',B'), & (A,B),(A',B')\in L_k,\\[3pt]
\beta^{\rm raw}_k(u,v)+e(A,B)+e(A',B')+e(A_S,B_S), & \text{otherwise},
\end{cases}\qquad k=3,4,5,6 .\;}
$$

The value of $e$ on the total depends only on its A and B layers, so this
is not an implicit equation for $C$.

## 6. The lower identity

With these $\alpha_k,\beta_k$,

$$
F_k\bigl(\mu_k(u,v)\bigr)=\mu_{k+1}\bigl(F_ku,\ F_kv\bigr)
$$

on all lower cochains, i.e. $\delta\alpha_k(A,A')+\alpha_{k+1}(\delta_sA,\delta_sA')=P(A+A')+P(A)+P(A')$
and $\delta\beta_k+\beta_{k+1}(p_ku,p_kv)=f_k(A_S,B_S)+f_k(A,B)+f_k(A',B')$.

## 7. The natural degree-six product $\ast$

Degree six additionally uses a *natural* lower product (no global branch
in its C correction), needed by the section construction:

$$
u\ast v=\bigl(A+A',\ B+B'+\alpha_6,\ C+C'+\widehat\beta_6\bigr),\qquad
\widehat\beta_6=\beta^{\rm raw}_6+\delta\lambda_{\rm loc},
$$

where $\lambda_{\rm loc}$ is the local extension (as for $q_{\rm loc}$,
testing the curvature of both pairs on every face) of the comparison

$$
\lambda_3(A,B;A',B')=\lambda_A(A,A')+I\Bigl[\beta_{\rm ext}+\beta^{\rm raw}_6\Bigr]\bigl(\pi^*A,B\ell;\ \pi^*A',B'\ell\bigr),
\qquad \delta\lambda_3=\beta^L_6+\beta^{\rm raw}_6\ \text{on legal pairs},
$$

with the extended correction for closed $A$ and arbitrary $B$

$$
\beta_{\rm ext}=\beta^\sharp(A,B;A',B')+L_6(A,B)+L_6(A',B')+L_6(A_S,B_S)+b\smile_5b',\qquad
L_6(A,B)=\xi_6(A,B)+q_{\rm loc}(A,B)+P(A)\smile_5b,
$$

($b=\delta B+P(A)$, $b'=\delta B'+P(A')$), and $\lambda_A(A,A')$ the
A-only pair primitive of the residual
$\beta_{\rm ext}(A,0;A',0)+\beta^{\rm raw}_6(A,0;A',0)$
([universal_primitives.md](universal_primitives.md), Section 5). On a pair
with legal AB, $\widehat\beta_6=\beta^L_6$; $\widehat F_6$ is a homomorphism for $\ast$.
