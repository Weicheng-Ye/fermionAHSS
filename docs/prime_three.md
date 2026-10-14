# The two-layer model at the prime three

Away from the prime two the binary layers B and C vanish, and the window
of $ko$ at the prime three is the two-stage tower of the rows $q=0$ (layer
A) and $q=-4$ (layer D) with $k$-invariant $2\beta_3P^1\rho_3$. Relations of
a generator of order $3^a$ (which lives in A) are measured in this model
instead of the four-layer one. States are $(A,D)$ with
$A\in C^{k-3}(X;\mathbb Z_s)$, $D\in C^{k+1}(X;\mathbb Z_s)$. Below package
degree five $P^1=0$ on the A degrees and the tower splits
($3^a\widetilde q=0$). Notation: [conventions.md](conventions.md); the
integral triple products are transported ones, $x\smile y\smile z=(x\smile_sy)\smile_sz$.

## Package degree five ($|A|=2$)

$P^1_s\rho_3A$ is the cube $\rho_3(A\smile A\smile A)$, and the potential
$\tfrac23\mathrm L(P^1_s\rho_3A)$ differs from $\Omega'(A)=\tfrac23A\smile A\smile A$
by an integral cochain, an isomorphism of stacking models. In the model of
$\Omega'$:

$$
\boxed{\;
\mathfrak d_5(A,D)=(\delta_sA,\ \delta_sD),\qquad
(A,D)\times(A',D')=\bigl(A+A',\ D+D'+\gamma(A,A')\bigr),\qquad
\gamma(A,A')=-2\,\bigl(A\smile A\smile A'+A\smile A'\smile A'\bigr),\;}
$$

(the curvature of a cocycle vanishes). $\gamma$ differs from the cross
effect $\Omega'(A)+\Omega'(A')-\Omega'(A+A')$ by $\tfrac23\delta_s\Xi$ with

$$
\Xi(A,A')=2A(A\smile_1A')+(A\smile_1A')A+2(A\smile_1A')A'+A'(A\smile_1A')
$$

(from the Hirsch identity $A'A-AA'=\delta_s(A\smile_1A')$, word by word). A
gauge $(u,w)$ of degree four has boundary state $(\delta_su,\delta_sw)$ and
acts by the product with it. The relation of a generator of order $3^a$
with $3^aA=\delta_su$ and flat lift $(A,0)$ is

$$
\boxed{\;
x_D=2\cdot3^{a-1}\,A\smile A\smile A
+\tfrac23\,\delta_s\Bigl[\sum_{j=1}^{3^a-1}\Xi(jA,A)-u\smile\delta_su\smile\delta_su\Bigr]\pmod{3^aH_D},\;}
$$

whose last term is a Bockstein class and is computed, not dropped.

## Package degree six ($|A|=3$)

$P^1_s\rho_3A$ is the nineteen-term cyclic-diagonal formula and the model
uses $\Omega(A)=\tfrac23\mathrm L(P^1_s\rho_3A)$ itself. With $a=\rho_3A$,
$a'=\rho_3A'$, and the two natural primitives of the reduced power built
from the cyclic diagonals $D_0,D_1,D_2,D_3$ (boundaries $\rho-1$ and
$N=1+\rho+\rho^2$ alternately; $\rho$ rotates the tensor factors):

- $\varphi(a,a')$: the mixed words $a\otimes a\otimes a'+a\otimes a'\otimes a'$ evaluated on $(\rho^2+2\rho)D_3$, with $\delta_s\varphi=P^1(a+a')-P^1(a)-P^1(a')$ modulo three;
- $\chi(u)$ for a degree-two $u$: $u\otimes\delta u\otimes\delta u$ on $D_2$, minus $(u\otimes\delta u\otimes u-u\otimes u\otimes\delta u)$ on $D_1$, plus $u\otimes u\otimes u$ on $D_0$, with $\delta_s\chi=P^1(\delta_su)$ modulo three

(each factor transported to the first vertex of the simplex, products
reduced modulo three):

$$
\boxed{\;
\begin{aligned}
\mathfrak d_6(A,D)&=\bigl(\delta_sA,\ \delta_sD+J(A)\bigr),\qquad J(A)=\delta_s\Omega(A)=\tfrac23\,\delta_s\mathrm L(P^1_s\rho_3A)=2\beta_3P^1_s\rho_3A,\\
\gamma(A,A')&=\tfrac23\Bigl[\mathrm L(P^1a)+\mathrm L(P^1a')-\mathrm L\bigl(P^1(a+a')\bigr)+\delta_s\mathrm L\bigl(\varphi(a,a')\bigr)\Bigr],\\
\text{gauge }(u,w)&\ \text{has boundary state }\bigl(\delta_su,\ \delta_sw+D_u\bigr),\qquad
D_u=-\tfrac23\Bigl[\mathrm L\bigl(P^1\rho_3\delta_su\bigr)-\delta_s\mathrm L\bigl(\chi(u)\bigr)\Bigr].
\end{aligned}\;}
$$

Both brackets are divisible by three (checked on every value). The flat
lift $(A,0,0,D)$ solves $\delta_sD=-J(A)$, possible exactly when the class
survives to E6; associativity and commutativity hold up to signed
D-coboundaries. Both degrees are natural on a resolution with transported
integral cups and cyclic diagonals through stage three
([natural_evaluation.md](natural_evaluation.md), Section 3.4).

## Results

Degree five: $\mathbb Z/9$ for $\mathbb Z/3$ (the relation $3A=-16D$),
$\mathbb Z/3\oplus\mathbb Z/27$ for $\mathbb Z/9$ ($9a=6d$ modulo nine),
$\mathbb Z/9\oplus\mathbb Z/81$ for $\mathbb Z/27$: the $I$-adic filtration
of the representation ring. Degree six by suspension
($ko^3(\Sigma BG)=ko^2(BG)$): $\mathbb Z/9$ for $\mathbb Z/3\times\mathbb Z$
($3A=14D$), $\mathbb Z/3\oplus\mathbb Z/27$ for $\mathbb Z/9\times\mathbb Z$.
At the primes $p\ge5$ the rows $q=0$ and $q=-4$ lie in different Adams
summands and every relation splits.
