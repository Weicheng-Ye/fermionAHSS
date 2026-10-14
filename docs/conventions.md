# Conventions and the interface data

Everything in these notes is written on normalized ordered simplicial
cochains of a space $X$ (in practice the normalized bar resolution of a
group, or a free resolution compared with it). This page fixes the
notation and lists, once and for all, the pieces of the secondary and
tertiary cohomology operations that the stacking formulas consume as
black boxes (Section 8). Everything else in the notes is written out.

## 1. Symbols

| Symbol | Meaning |
| --- | --- |
| $s\in Z^1(X;\mathbb F_2)$, $\omega\in Z^2(X;\mathbb F_2)$ | the two fixed twist cocycles |
| $\mathbb Z_s$, $\mathbb Q_s$ | integers (rationals) with monodromy $(-1)^s$, the sign local system |
| $\rho$ | reduction modulo two of an integral cochain; $\rho_3$ reduction modulo three |
| $\widetilde x$ | the pointwise $\{0,1\}$-valued integral lift of a binary cochain $x$; $\mathrm L(y)$ the pointwise $\{0,1,2\}$ lift of a mod-three cochain |
| $\mathsf h(x)=\tfrac12\widetilde x$ | the *half lift*, a specific rational cochain (not a class in $\mathbb Q/\mathbb Z$) |
| $\delta$ | ordinary coboundary, binary or integral |
| $\delta_s$ | the signed coboundary, (C1) |
| $\smile_i$ | the interval-cut cup-$i$ product; $\smile=\smile_0$ |
| $x\cdot y$ | pointwise product of two binary cochains of the same degree (value $x(\sigma)y(\sigma)$) |
| $k$ | package degree; $n=k-3$ is the degree of the A layer; the four layers have degrees $n,n+1,n+2,n+4$ |

Rules that every formula depends on:

- A binary equation holds modulo two. A rational *phase* identity holds
  modulo integral cochains unless it is called a literal identity.
- A binary expression is reduced modulo two **as a whole** and then
  lifted; **different occurrences of $\mathsf h$ or of a lift are lifted
  separately**. Lifting a sum is not the sum of the lifts; the difference is
  the pointwise product:
  $$\tfrac12\bigl(\widetilde x+\widetilde y-\widetilde{x+y}\bigr)=\widetilde{x\cdot y}.\tag{C0}$$
- A division by $2$, $4$, $8$ or $3$ in an integral formula is exact
  integer division of the fully assembled numerator (the code raises an
  error otherwise). Integral carries are never discarded.
- Negative cochain degrees and negative cup indices are zero. A word in
  which a degree-zero input occurs twice is zero.
- Cochains are normalized: they vanish on a simplex with a repeated
  adjacent vertex.
- Products keep their written order; the coefficient system of a product
  is the tensor product of those of its factors.

## 2. Signed coboundary and transport

For $c$ of degree $r$, in the frame of the first vertex,

$$
(\delta_sc)(v_0,\dots,v_{r+1})=(-1)^{s(v_0,v_1)}\,c(v_1,\dots,v_{r+1})
+\sum_{j=1}^{r+1}(-1)^j\,c(v_0,\dots,\widehat{v_j},\dots,v_{r+1}),
\qquad\text{i.e.}\quad \delta_sc=\delta c-2\,\widetilde s\smile c .
\tag{C1}
$$

The **transported product** of $x$ (degree $r$) with a sign-valued right
factor $y$ is

$$
(x\smile_sy)(v_0,\dots,v_{r+t})=x(v_0,\dots,v_r)\,(-1)^{s(v_0,v_r)}\,y(v_r,\dots,v_{r+t}),
\tag{C2}
$$

and the **signed cube** of an integral $A$ is $A^{\smile3}=(A\smile_sA)\smile_sA$.
The integral **Pontryagin cochain** of the twist is

$$
P_\omega=\widetilde\omega\smile\widetilde\omega+\widetilde\omega\smile_1\delta\widetilde\omega\in C^4(X;\mathbb Z).
\tag{C3}
$$

## 3. Cup-$i$ products and words

Let $w=(w_1,\dots,w_\ell)$ be a nondegenerate surjection word onto
$\{1,\dots,r\}$ and let the inputs $x_j$ have degrees $d_j$. The output
degree is $N=\sum_jd_j-\ell+r$. On an ordered simplex $[v_0,\dots,v_N]$ sum
over cuts $0=t_0\le t_1\le\dots\le t_\ell=N$, assign the closed interval
$[t_{a-1},t_a]$ of vertices to the input $x_{w_a}$, keep only the cuts for
which the vertex list of every input (the concatenation of its intervals)
is strictly increasing with exactly $d_j+1$ vertices, multiply the input
values on those faces and sum modulo two. This is $w(x_1,\dots,x_r)$.

- $x\smile_iy$ is the alternating word $1,2,1,2,\dots$ of length $i+2$.
- The **Cartan helpers**, for $|x|=1$ resp. $|x|=2$ and $|a|=n\ge1$ (zero for $n=0$):
  $$\zeta_{1,n}(x,a)=w_{1232\,\mathrm{alt}_n(4,3)}(x,x,a,a),\qquad
    \zeta_{2,n}(x,a)=w_{1231\,\mathrm{alt}_{n+1}(3,4)}(x,x,a,a),\tag{C4}$$
  where $\mathrm{alt}_m(u,v)$ is $m$ letters alternating $u,v,u,\dots$. For
  example $\zeta_{2,1}(\omega,b)=w_{123134}(\omega,\omega,b,b)$ and
  $\zeta_{1,2}(s,y)=w_{123243}(s,s,y,y)$. Output degrees $n+2$ and $n+3$.
- The integral cup-$i$ product uses the same cuts with the sign
  $$(-1)^{\,i(p+q)+\binom i2+\sum_a\epsilon_at_a+\sum_{a<b,\ w_a>w_b}m_am_b},\quad
    \epsilon_a=[\,w_a\text{ occurs again later}\,],\ m_a=t_a-t_{a-1}+\epsilon_a .\tag{C5}$$
- Binary cup-$i$ products satisfy
  $\delta(x\smile_iy)=\delta x\smile_iy+x\smile_i\delta y+x\smile_{i-1}y+y\smile_{i-1}x$.

## 4. Squares, the primary operation, Bocksteins

For a binary cochain $x$ of degree $r$, closed or not,

$$
Q^j(x)=x\smile_{r-j}x+x\smile_{r-j+1}\delta x,\qquad
E(x)=Q^2(x)+\omega\smile x,\qquad
Q_D(x)=E(x)+s\smile Q^1(x).
\tag{C6}
$$

(Degrees $r+j$, $r+2$, $r+2$.) For a cocycle $a$ of degree $n$,
$\operatorname{Sq}^ja=a\smile_{n-j}a$ and $Q_D(a)=D(a)$ with
$D=\operatorname{Sq}^2+s\operatorname{Sq}^1+\omega$, the primary
operation of the window. Identities:

$$
\delta Q_D(x)=Q_D(\delta x),\qquad \delta E(x)=E(\delta x),
\tag{C7}
$$
$$
E(u+v)+E(u)+E(v)=\delta(u\smile_{r-1}v)\quad\text{for closed }u,v\text{ of degree }r .
\tag{C8}
$$

Bocksteins: $\beta_s[z]=[\delta_s\widetilde z/2]$ for a binary cocycle
$z$; $e(a)=\rho(\delta\widetilde a/2)$, and $e(a)=a\smile_{n-1}a$ literally
for a binary cocycle. The names of the other primary operations are
$\operatorname{Dbar}=D\rho$ and $\operatorname{Dtilde}=\beta_s(\operatorname{Sq}^2+\omega)$.

## 5. Polarizations

For binary $x,y$ of degree $r$:

$$
h^j(x,y)=x\smile_{r-j+1}y+\delta x\smile_{r-j+2}y,\qquad
\delta h^j(x,y)+h^j(\delta x,\delta y)=Q^j(x+y)+Q^j(x)+Q^j(y),
\tag{C9}
$$

$$
h^D(x,y)=h^2(x,y)+s\smile h^1(x,y)
=x\smile_{r-1}y+\delta x\smile_ry+s\smile\bigl(x\smile_ry+\delta x\smile_{r+1}y\bigr),
\tag{C10}
$$

with $Q_D(x+y)+Q_D(x)+Q_D(y)=\delta h^D(x,y)+h^D(\delta x,\delta y)$. For
a closed $x$, $h^D(x,x)$ represents $(\operatorname{Sq}^1+s)[x]=\rho\beta_s[x]$.

## 6. The prism and the ramp

For a cochain $v$ of degree $r$ on the cylinder $X\times I$,

$$
(Iv)(v_0,\dots,v_{r-1})=\sum_{j=0}^{r-1}(-1)^j\,
v\bigl((v_0,0),\dots,(v_j,0),(v_j,1),\dots,(v_{r-1},1)\bigr),
\qquad \delta_sI+I\delta_s=i_1^*-i_0^* .
\tag{C11}
$$

The twists are pulled back to $X\times I$ unchanged. For a cochain $u$ on
$X$: $\pi^*u$ is the plain pullback (its prism is zero), and the **ramp**
$u\ell$ is the pullback multiplied by the interval coordinate of the last
vertex (zero at $i_0$, equal to $u$ at $i_1$). A nested prism uses a second
independent interval.

## 7. Phases, carries and three conventions used everywhere

- A *phase* is a rational cochain written with a specified lift; a
  statement "$\equiv$ mod $\mathbb Z$" is an equality of phases.
- $\frac14\,\widetilde x$ and $\frac18\,\widetilde x$ denote the stated
  lift divided by $4$ or $8$; $\frac14 X$ for an integral $X$ (such as
  $P_\omega$ or a signed cube) is the rational cochain with those values.
- For a function $Z$ of one lower triple and a lower product $\mu$,
  $\Delta Z(u,v)=Z(u)+Z(v)-Z(\mu(u,v))$ is its cross effect.

## 8. The secondary and tertiary data used by the stacking model

The secondary differential `Tau` and the tertiary `T` are written out in
[secondary_operations.md](secondary_operations.md) and
[tertiary_operations.md](tertiary_operations.md). The stacking formulas use
them only through the following named objects, for a signed cocycle $A$ of
degree $n$ ($\delta_sA=0$), with $a=\rho A$ (the third column says where
each is defined):

| Object | Degree | Definition / property used |
| --- | --- | --- |
| $P(A)=Q_D(a)$ | $n+2$ | explicit, (C6); closed |
| $t(A)=\rho\bigl((A-\widetilde a)/2\bigr)$ | $n$ | the parity carry of $A$; explicit |
| $g(A)$ | $n+2$ | the binary **source** of the secondary splitting, (S) of the secondary page, Section 6 ($n=0$ below) |
| $k_0(A)$ | $n+3$ | the binary **A-only part** of the secondary representative, (S) there ($n=0$ below); $\delta k_0=D(P(A))$ |
| $\lambda(A,B)=g(A)\cdot(s\smile B)$ | $n+2$ | the **splitting carry**, i.e. $\rho\bigl(\tfrac12(\widetilde g+\widetilde{s\smile B}-\widetilde{g+s\smile B})\bigr)$ by (C0); explicit given $g$ |
| $\tau'_n(A;B)=Q_D(B)+k_0(A)+\delta\lambda(A,B)$ | $n+3$ | the secondary representative with defining cochain $B$, valid when $\delta B=P(A)$ (**verified numerically for $n=0$**); it depends on $A$ only through $a$ and $t(A)$, so it is four-periodic in $A$; $\tau'_n(0;b)=Q_D(b)$ for closed $b$ |
| $\theta_m(q,v)$ | $m+3$ | the boundary phase (B) of the tertiary page, Section 2, for a binary cocycle $q$ of degree $m$ and $v$ with $\delta v=D(q)$; $\theta_m(q,0)$ is also used |
| $V_P(A)$ | $n+4$ | the universal A-only source primitive, $\delta_sV_P(A)\equiv\theta_{n+2}(P(A),k_0(A))\pmod{\mathbb Z}$, $n=1,2,3$: the tertiary $V_n$ of Section 8 there, with the prism adjustment at $n=1$ of [product_upper_k4_k6.md](product_upper_k4_k6.md), Section 1 |
| $\Omega_k(A,B,C)=\mathsf h(E(C))+\widehat R_{k-3}(A,B)\,[+\tfrac23\mathrm L(P^1_s\rho_3A)]$ | $k+1$ | the tertiary phase: $R_0$ and $\widehat R_1,\widehat R_2,\widehat R_3$ of the tertiary page, Sections 3–6 |
| $P^1_s\rho_3A$ | $n+4$ | the signed reduced power modulo three: $\rho_3(A\smile_sA\smile_sA)$ for $n=2$, the nineteen-term formula (P3) of the tertiary page, Section 7, for $n=3$ |

For $n=0$ ($A$ a signed $0$-cocycle, so $a$ is locally constant and
$t=t(A)$ satisfies $\delta t=s\smile a$):

$$
g(A)=\omega\smile t+a\,(s\smile s)+a\,(\omega\smile_1s),\qquad
\tau'_0(A;B)=B\smile\delta B+\omega\smile B+\rho\Bigl(\tfrac12\,\delta\,\widetilde{g'}\Bigr)+s\smile g'+s\smile s\smile s\smile a,
\quad g'=g(A)+s\smile B,
\tag{C12}
$$

(the numerator $\delta\widetilde{g'}$ is even), and equivalently
$\tau'_0(A;B)=Q_D(B)+k_0(A)+\delta\lambda(A,B)$ with the shorter source

$$
g^\flat(A)=\omega\smile t+a\,(s\smile s),\qquad
k_0(A)=\rho\Bigl(\tfrac12\bigl(-\delta_s\widetilde{g^\flat}+\widetilde s\smile\widetilde{\omega a}\bigr)\Bigr)+s\smile s\smile s\smile a,\qquad
\lambda=g^\flat\cdot(s\smile B).
\tag{C13}
$$

Both forms of (C12)–(C13) agree with the old implementation on random
data (the degree-three legal $\beta$ uses $g^\flat$).

Finally, the **universal primitives** $\mathcal V_n$, $\mathcal U^\Theta_m$,
$\mathcal P_k$, $v_n$, $q_k$ that appear in the product are *not* black
boxes: they are specified in [universal_primitives.md](universal_primitives.md)
by their source formula, their universal model and their fixed coefficients.
