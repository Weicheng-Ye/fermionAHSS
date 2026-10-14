# Relations of an A generator: squaring, the diagonal, and one-fiber primitives

The expensive object in a production run is the pair primitive
$\mathcal P_k(A,A')$ inside the D-layer pair phase of
[product_upper_k4_k6.md](product_upper_k4_k6.md): a contraction over a
two-fiber universal chain, each term evaluating a source that itself
contains three one-fiber primitives $V_n$ and two $\Theta$-pair primitives.
This page works out what the relation of an A generator actually needs
and shows that every two-fiber rational primitive in it can be replaced by
a one-fiber primitive of an explicit **diagonal source**, plus finitely
many universal constants. The previous implementation did this in degree
six only (its "unary table" of 22 constants); here it is made explicit in
every degree. Notation: the stacking pages; $k\in\{4,5,6\}$, $n=k-3$.

## 1. Setting

Let $q$ be an A generator of order $m=2^e$ (odd orders are handled by the
prime-three rows and $p\ge5$ needs no measurement), with marked cocycle
$A$ and flat lift $u=(A,B,C,D)$: $\delta B=P(A)$, $\delta C=\tau'_n(A;B)$,
$\delta_sD=-J_k(u)=-\delta_s\Omega_k(u)$. Hence

$$
\Phi:=D+\Omega_k(u)\ \text{ is a rational cocycle of degree }k+1 .
$$

On closed states the D correction of the product is

$$
\gamma_k(x,y)=\Omega_k(x)+\Omega_k(y)-\Omega_k(x\cdot y)+\delta_se_k(x,y),
\qquad e_k=e^{\rm prod}_k+\chi_k ,
$$

with $e^{\rm prod}_k$ the pair phase of the product page and $\chi_k(C,C')$
the explicit pure-C normalization (no A, hence no pair primitive).

## 2. Powers by squaring telescope

Put $u_0=u$ and $u_{i+1}=u_i\cdot u_i$, so $u^m=u_e$. The lower layers are
explicit:

$$
A_i=2^iA,\qquad B_0=B,\ B_1=h^D(a,a),\ B_i=0\ (i\ge2),\qquad
C_0=C,\ C_{i+1}=\beta^L_k(A_i,B_i;A_i,B_i),
$$

because $B_{i+1}=2B_i+\alpha_k(A_i,A_i)$ and $\alpha_k(A_i,A_i)=h^D(\rho A_i,\rho A_i)$
vanishes for $i\ge1$. For $\Phi_i:=D_i+\Omega_k(u_i)$ the product rule gives
$\Phi_{i+1}=2\Phi_i+\delta_se_k(u_i,u_i)$, hence

$$
\boxed{\;(u^m)_D=m\,\Phi-\Omega_k(u^m)+\delta_s\sum_{i=0}^{e-1}2^{\,e-1-i}\,e_k(u_i,u_i).\;}
\tag{P1}
$$

No phase of a product $u_i\cdot u_j$ with $i\ne j$ ever appears. The relation
is read $2$-locally, and the $2$-adic denominators of $e_k$ divide $2^{t_k}$
with $t_k=5,4,4$ for $k=4,5,6$ (the implementation's tail bound), so only
the last $t_k$ squares contribute: $i\ge e-t_k$.

## 3. The gauge step has a zero fiber

The A layer of $u^m$ is the coboundary $mA$. The measurement writes

$$
u^m=\mathfrak d_{k-1}(g)\times x_0,\qquad g=(U,Y,W),\ \delta_sU=mA,\qquad x_0=(0,B_0,C_0,D_0),
$$

with $Y,W$ solved by linear algebra on the resolution so that the lower
layers agree; $x_0$ carries the relation
([extension_assembly.md](extension_assembly.md), Section 3). Since
$\mathfrak d_{k-1}(g)$ is closed, the product rule and (P1) give

$$
\boxed{\;D_0=m\,\Phi-\Omega_k\bigl(\mathfrak d_{k-1}g\bigr)-\Omega_k(x_0)-G_{k-1}(g)
+\delta_s\Bigl[\sum_{i=e-t_k}^{e-1}2^{\,e-1-i}e_k(u_i,u_i)-e_k\bigl(\mathfrak d_{k-1}g,\ x_0\bigr)\Bigr],\;}
\tag{P2}
$$

read modulo $\delta_sC^k(\mathbb Z_s)$; the term $\Omega_k(u^m)$ cancels.
Here $G_{k-1}(g)=-I\,\mathcal K_k(\text{cylinder of }g)$ is the off-legal
prism of the differential (one-fiber: $V_n$ on the cylinder), and the pair
phase $e_k(\mathfrak d_{k-1}g,x_0)$ has A layers $(mA,0)$: every pair
primitive ($\mathcal P_k$, $\mathcal V_n$, $\mathcal U^\Theta$) vanishes
identically on a pair with a zero fiber, so this phase is explicit.
Therefore

$$
\text{the only two-fiber primitives in the relation are }\ \mathcal P_k(A_i,A_i)\ \text{and}\ \mathcal V_n(A_i,A_i),\qquad A_i=2^iA,\ e-t_k\le i<e .
$$

## 4. The diagonal source, explicitly

Put $A'=A$ in the pair source $S_k(A,A')$ of the product page
($S^{\rm aux}+\Delta R^{(k)}_A+\ldots$) and use $\varphi=\varphi_0-V_P$:

$$
\boxed{\;
S^\Delta_k(A):=S_k(A,A)=-2V_P(A)+V_P(2A)+2R^{(k)}_A(A)-R^{(k)}_A(2A)
-\varphi_0\bigl(2A,\ \alpha,\ \beta^\sharp(A,0;A,0)\bigr)-\Theta^{\rm pair}_{n+2}(P,k_0;P,k_0)+X_k(A),\;}
\tag{P3}
$$

$$
\alpha=h^D(a,a),\qquad
X_4(A)=X_6(A)=\mathsf h\bigl(\mathfrak J(A,A)\bigr)=\mathsf h\bigl(J_0(\alpha)+j_0(P,P)+j_0(0,\delta\alpha)\bigr),\qquad
X_5(A)=\mathfrak B(A,A)=\mathsf h\Bigl(s\smile\bigl(Q_D(\alpha)+h^D(P,P)\bigr)\Bigr)
$$

($\mathfrak T^{\rm s}(A,A)=\mathsf h(a\smile P+a\smile P)=0$). Since
$\Theta^{\rm pair}$ contains $\mathcal U^\Theta_{n+2}(P,P)$ once and
$2R^{(k)}_A(A)-R^{(k)}_A(2A)$ contains $\mathcal U^\Theta_{n+2}(q,q)$ twice
(with $q=P$ for $k=5,6$ and the adjusted $q=(\omega+s^2)a$ for $k=4$,
while $q(2A)=0$), the diagonal source holds the $\Theta$-pair primitive
**at one argument $(P,P)$ for $k=5,6$**, the value the curvature already
uses, and at the two arguments $(q,q)$, $(P,P)$ for $k=4$.

The doubled argument is simpler than $A$: $\rho(2A)=0$, so
$P(2A)=0$, $\alpha(2A)=0$, $t(2A)=a$, $g^\flat(2A)=\omega\smile a$, and

$$
k_0(2A)=\rho\Bigl(-\tfrac12\,\delta_s\widetilde{\omega\smile a}\Bigr),\qquad
R^{\rm red}(2A)=\mathsf h\bigl(s\smile k_0(2A)+\delta k_0(2A)\bigr),
$$

while the integral terms of $R^{(k)}_A$ scale: $\tfrac14P_\omega\smile_s2A=\tfrac12P_\omega\smile_sA$,
$\tfrac14(2A)^{\smile3}=2A^{\smile3}$, and the three-primary terms use
$\rho_3(2A)=-\rho_3A$. For the higher squares $A_i=2^iA$, $i\ge1$, the
input itself is even: then $P=0$, $\alpha=0$,
$\Theta^{\rm pair}_{n+2}(0,k_0;0,k_0)=\mathsf h\bigl(h^2(k_0,k_0)\bigr)$ and
$\mathcal U^\Theta=0$, so (P3) contains only $V_P$ at two arguments,
$R^{\rm red}$ at two arguments, explicit binary words in $k_0$, $t$, and the
integral terms.

## 5. The one-fiber primitive and its constants

Let $(\mathcal M_n;F,G,H)$ be the one-fiber universal model of a closed
signed degree-$n$ cocycle with the $\omega$ factor, with its contraction
onto the small complex $\mathfrak S_n$ (the model of $V_n$:
[tertiary_operations.md](tertiary_operations.md), Section 8; the dihedral
wedge for $n=1$, the signed degree-two model for $n=2$, the degree-three
model for $n=3$), $j$ its closed-cocycle section, and
$\Delta:\mathcal M_n\to\mathcal M^{\rm pair}_n$ the diagonal, so that the
pair section of $(A,A)$ is $\Delta j(A)$. Define, for a simplex $\sigma$,

$$
\boxed{\;
\widehat{\mathcal P}_k(A)(\sigma)=\sum_{\tau}c_\tau\,S^\Delta_k(A_\tau)+\sum_{b}c_b\,\langle Fj_\sigma(A),\,b\rangle,
\qquad Hj_\sigma(A)=\sum_\tau c_\tau\,\tau,\qquad
c_b:=\mathcal P_k\bigl(\Delta\,G b\bigr)\in\mathbb Q,\;}
\tag{P4}
$$

where $\tau$ runs over the universal simplices of the chain $Hj_\sigma(A)$,
$A_\tau$ is the cocycle labelled by $\tau$, and $b$ runs over the cells of
$\mathfrak S_n$ in degree $k$ with nonzero A fiber. The constants $c_b$
are the values of the two-fiber primitive on the finitely many universal
chains $\Delta Gb$: they are computed once per $k$ and stored.

**Identity.** For every $k$-chain $x$ of $\mathcal M_n$,
$x=\partial Hx+H\partial x+GFx$, so, modulo integers,

$$
\Delta^*\mathcal P_k(x)=(\delta\Delta^*\mathcal P_k)(Hx)+\Delta^*\mathcal P_k(H\partial x)+\sum_b\langle Fx,b\rangle\,c_b
\equiv S^\Delta_k(Hx)+(\delta W)(x)+\sum_b\langle Fx,b\rangle\,c_b,\qquad W:=\Delta^*\mathcal P_k\circ H .
$$

Hence, as natural cochains in $A$,

$$
\boxed{\;\mathcal P_k(A,A)=\widehat{\mathcal P}_k(A)+\delta_sW(A)+N(A),\qquad W\ \text{rational of degree }k-1,\ N\ \text{integral},\;}
\tag{P5}
$$

and $\delta_s\mathcal P_k(A,A)=\delta_s\widehat{\mathcal P}_k(A)+\delta_sN(A)$:
using $\widehat{\mathcal P}_k$ in (P2) changes $D_0$ by the coboundary of a
natural integral cochain, a pure-D gauge, which the relation does not see.
Nothing else in (P2) changes, because the lower layers $B_i,C_i$ are kept
exact. (This is the degree-six "unary primitive" of the previous code,
whose 22 constants are exactly the $c_b$ on `r3_chain.total_basis(6)`;
its proof is the three lines above.)

The sizes of the small complexes fix the number of constants:

| $k$ | one-fiber small cells in degree $k$ | constants $c_b$ | two-fiber small cells (old coefficients) |
| --- | ---: | ---: | ---: |
| 4 | dihedral wedge, degree 4 | 3 | 2 |
| 5 | `basis2(5)` of the signed degree-two model | 4 | 48 |
| 6 | `r3_chain.total_basis(6)` | 22 | 12 |

The binary pair primitive $\mathcal V_n(A,A)$ inside $C_{i+1}$ admits the
same identity with $\mathbb F_2$ constants, but there the coboundary
$\delta W'$ is not harmless: it moves the C layer by a natural
coboundary, i.e. by a pure-C gauge whose D companion needs $W'$ itself.
Computing $W'=\Delta^*\mathcal V_n\circ H$ needs the pair chains again, so
$\mathcal V_n(A_i,A_i)$ is kept as a binary pair contraction. Its source
is a cup word (no nested primitive) and its chain is one degree lower than
that of $\mathcal P_k$, so it is cheap by comparison.

## 6. The relation, written out

With (P2), (P4) and the explicit parts of $e^{\rm prod}_k$ on the diagonal,

$$
\boxed{\;
D_0\equiv m\,\Phi-\Omega_k(\mathfrak d_{k-1}g)-\Omega_k(x_0)-G_{k-1}(g)
+\delta_s\Bigl[\sum_{i=e-t_k}^{e-1}2^{\,e-1-i}\,\widehat e_k(u_i,u_i)-e_k(\mathfrak d_{k-1}g,x_0)\Bigr]
\pmod{\delta_sC^k(\mathbb Z_s)},\;}
$$

$$
\widehat e_k(u_i,u_i)=-\widehat{\mathcal P}_k(A_i)-I\,\mathrm{RFS}(u_i,u_i)-\Xi_k(u_i,u_i)-\bigl[2Q_k(u_i)-Q_k(\mu^\sharp(u_i,u_i))\bigr]+\chi_k(C_i,C_i),
$$

where $\Xi_k$ is the explicit middle term of $e^{\rm prod}_k$ on the
diagonal: $\Xi_4=\mathsf h(\mathfrak j(A,B;A,B))$,
$\Xi_5=\mathfrak B'(A,B;A,B)+\mathfrak T(A,B;A,B)=\mathsf h\bigl(s\smile(h^D(B,B)+h^D(0,\alpha))\bigr)$,
$\Xi_6=\mathsf h(\mathfrak j(A,B;A,B))$. The relation $m\widetilde q=t$ is
then read from $x_0$ layer by layer as on the assembly page. Everything
on the right is a one-fiber quantity: $\Omega_k$ at three states, $V_n$ at
the arguments $A_i$, $2A_i$, the gauge cylinder, and $\mathcal U^\Theta_{n+2}(P,P)$
once.

Since $H^{k+1}(X;\mathbb Q)=0$ for a finite group, the rational cocycle
$m\Phi-\ldots$ above has zero class, and the integral class of $D_0$ is the
Bockstein of the $\mathbb Q/\mathbb Z$ class of the bracket: this is why the
universal primitives cannot be dropped, and why only their values modulo
integers, on the small cycles, matter.

## 7. What is left to evaluate

| Object | Before | After |
| --- | --- | --- |
| $\mathcal P_k(A_i,A_i)$ | one two-fiber chain per simplex and per square, each term with $3V_n+2\,\mathcal U^\Theta$ inside | none |
| $\widehat{\mathcal P}_k(A_i)$ | — | one one-fiber chain per simplex and per square, each term with $V_n$ at two labels and $\mathcal U^\Theta$ at one |
| constants | — | $c_b$ on 3, 4, 22 cells, computed once with the two-fiber contraction |
| $\mathcal V_n(A_i,A_i)$ | binary pair chain | unchanged (cup-word source) |
| $V_n$, $\Theta_m$ leaves | per two-fiber term | per one-fiber term, cached by label |

The remaining nesting is two levels ($\widehat{\mathcal P}\to V_n\to\theta_m$)
instead of three with a pair chain on top. The measured chain sizes are in
Section 9.

## 8. What has to be verified before use

1. The squaring recursion (P1) is exact for the products as implemented,
   including the pure-C normalization $\chi_k$ (closed states, zero residual).
2. The constants $c_b$: computed by the two-fiber contraction on the chains
   $\Delta Gb$, once per $k$, and checked on a group with a known answer
   (C4 untwisted in degree five returns $\mathbb Z/4$ with $2\widetilde q=D$;
   the degree-six table of the previous code is the $k=6$ instance).
3. The tail bound $t_k$: the $2$-adic denominators of
   $\widehat e_k$ on the squares (they are those of $e_k$, by (P5)).

## 9. Measured chain sizes

Random closed signed cocycles $A=\delta_sx$ (labels of $x$ in $\{-1,0,1\}$)
with random $s,\omega$, built with the previous implementation's universal
models: the one-fiber chain $Hj(A)$ of the model of $V_n$ against the
two-fiber chain of $(A,A)$ over which $\mathcal P_k(A,A)$ is summed.
Degree five ($n=2$, a 5-simplex; `source_primitive.Htot` against
`r2_pair_chain.Htot`):

| input | max label | one-fiber terms (build) | two-fiber terms (build) | ratio |
| --- | ---: | ---: | ---: | ---: |
| seed 0, $A$ | 3 | 2,146 (0.2 s) | 15,996 (5.0 s) | 7.5 |
| seed 0, $2A$ | 6 | 3,635 (0.3 s) | 34,844 (10.9 s) | 9.6 |
| seed 0, $4A$ | 12 | 7,693 (0.6 s) | — | |
| seed 1, $A$ | 3 | 4,197 (0.3 s) | 41,077 (11.8 s) | 9.8 |
| seed 1, $2A$ | 6 | 7,671 (0.5 s) | 90,787 (26.0 s) | 11.8 |
| seed 1, $4A$ | 12 | 16,563 (1.2 s) | — | |
| seed 2, $A$ | 3 | 832 (0.2 s) | 3,727 (1.6 s) | 4.5 |
| seed 2, $2A$ | 6 | 1,138 (0.1 s) | 6,169 (1.2 s) | 5.4 |

Degree six ($n=3$, a 6-simplex; `r3_chain.Htot` against `r3_pair_chain.Htot`):

| input | max label | one-fiber terms (build) | two-fiber terms (build) | ratio |
| --- | ---: | ---: | ---: | ---: |
| seed 0, $A$ | 4 | 4,124 (6.1 s) | 23,003 (70.1 s) | 5.6 |
| seed 0, $2A$ | 8 | 7,841 (11.2 s) | 34,707 (99.0 s) | 4.4 |
| seed 1, $A$ | 4 | 17,871 (22.9 s) | 59,213 (154.2 s) | 3.3 |
| seed 1, $2A$ | 8 | 36,729 (47.9 s) | 110,194 (286.9 s) | 3.0 |
| seed 2, $A$ | 3 | 3,459 (4.5 s) | 14,988 (40.9 s) | 4.3 |
| seed 2, $2A$ | 6 | 7,771 (9.6 s) | 28,757 (71.2 s) | 3.7 |

Per term, the two-fiber source needs $V_n$ at three labels ($A$, $A'$,
$A+A'$) and $\mathcal U^\Theta$ at two pairs; the diagonal source (P3) needs
$V_n$ at two labels ($A$, $2A$) and $\mathcal U^\Theta$ at one. So per
square and per simplex the number of leaf evaluations drops by a factor of
about $1.5\times$ the term ratio: $7$ to $18$ in degree five, $4.5$ to $8$ in
degree six, before any caching; the chain construction itself drops by
$10$ to $50$. Both chains grow with the labels, and the squares
$(2^iA,2^iA)$ of the tail are exactly the large-label inputs.

**Where the remaining time is.** Evaluating the diagonal source (P3) on
every term of the seed-0 degree-five one-fiber chain, on a cold cache,
took 775 s for the 2,146 terms (all with distinct labels): 0.36 s per
term, almost all of it in the two $V_2$ values and their $\theta_4$ words.
On a smaller input (labels at most 2, a 6-simplex, Section 10) the
one-fiber value of a 5-face took 110 to 128 s against 207 to 415 s for
the two-fiber value, both cold. So after the reduction the cost per
simplex is the number of one-fiber terms times the cost of a $V_n$ leaf;
the leaf is the next target, by caching $V_n$ by label as the previous
code did and, beyond that, by tabulating it as the quasi-polynomial in
the labels that it is.

## 10. Measured saving on the C4 benchmark

The reduction was implemented in the previous code base exactly as its
degree-six unary table, as a timing patch on a copy of the frozen
benchmark source (`runs/unary5-benchmark-20261013/` of the previous
repository: `python/unary_gamma5.py`, an environment switch in
`extension_light.diagonal_phase`, the source registered in the persistent
universal-value store, and the table `data/unary-gamma5-coefficients.json`).
The four constants (P4) on `basis2(5)`, each computed in 12 to 249 s from
the two-fiber primitive on the chains $\Delta Gb$ (60 to 180 terms, labels
at most 2):

| cell $b$ | $((0,1),(2))$ | $((1,1),(1))$ | $((3,1),())$ | $((1,2),())$ |
| --- | ---: | ---: | ---: | ---: |
| $c_b$ (exact) | $1/4$ | $13/8$ | $-63/8$ | $19/4$ |
| $c_b \bmod 1$ | $1/4$ | $5/8$ | $1/8$ | $3/4$ |

The benchmark is the one of the previous repository's
`runs/c4-benchmark-20261003`: C4, $s=\omega=0$, package degree 5, rank-one
resolution, one CPU, fresh processes, cold = empty persistent store with
the bundled values, warm = the cold run's store. Same node for all four
runs (AMD EPYC 7543; the original benchmark's node was faster, 806 s cold):

| run | wall | extension stage | peak RSS | store written |
| --- | ---: | ---: | ---: | ---: |
| two-fiber (old), cold | 1231 s | 1229 s | 3.7 GB | 12 MB |
| one-fiber, cold | 331 s | 329 s | 2.2 GB | 1.5 MB |
| two-fiber (old), warm | 29 s | 27 s | 0.54 GB | |
| one-fiber, warm | 24 s | 21 s | 0.54 GB | |

All four runs return the same result: invariants $\mathbb Z/4$, relation
matrix $[[2,0],[-1,2]]$ ($2D=0$, $2A=D$), Smith diagonal $(1,4)$,
certificate level light-R, no heavy measurements. On a second, faster
node the one-fiber cold run took 267 s and 276 s and the warm run 20 s.

**Direct check of (P5).** On a random closed signed degree-two cocycle on
a 6-simplex (labels at most 2, random $s,\omega$), both primitives were
evaluated on all seven 5-faces (one-fiber 70 to 284 s per face, two-fiber
186 to 896 s) and their coboundaries compared at the 6-simplex:

$$
\delta_s\mathcal P_5(A,A)-\delta_s\widehat{\mathcal P}_5(A)=36,\qquad
\delta_s\mathcal P_5(A,A)-S_5(A,A)=34,\qquad
\delta_s\widehat{\mathcal P}_5(A)-S_5(A,A)=-2,
$$

all integers, as (P5) requires; the face values themselves differ by
non-integers ($\tfrac12$, $6$, $\tfrac{15}2$, $-\tfrac{17}2$, $\tfrac52$, $\tfrac{21}4$, $\tfrac{113}4$),
which is the coboundary term $\delta_sW$.

So the cold relation is 3.7 times faster, and the warm run is unchanged
within noise: once the universal values are stored, the remaining time
is the group-independent scaffolding and the binary chains, which the
previous profile already identified. The cold factor is bounded by what
is left per simplex, the $V_2$ leaves of the diagonal source; it is not
the chain-size ratio of Section 9 because the old run's pair source
values were shared between the two squares through the store while the
one-fiber run recomputes its leaves.

## 11. The other slow cases

The same question was asked of every case that the previous repository's
run reports found slow. Each is a different formula.

**Twisted C4, package degree 5** (1812 s in the transport-reduction
report, "universal-phase dominated"). The same A-over-D row as the
benchmark, with $\omega\ne0$: the relation of the order-four A generator
through its squares. The one-fiber primitive applies unchanged. On the
benchmark node it returns the documented $\mathbb Z/2\oplus\mathbb Z/32$
with the old relation matrix. Same node, one CPU, cold store:

| run | wall | extension stage | peak RSS | store written |
| --- | ---: | ---: | ---: | ---: |
| two-fiber (old), cold | 1614 s | 1612 s | 3.3 GB | 14 MB |
| one-fiber, cold | 518 s | 516 s | 2.3 GB | 2.2 MB |

Both return $\mathbb Z/2\oplus\mathbb Z/32$ with the relation matrix
$[[4,0,0,0],[-1,2,0,0],[-3,0,2,0],[0,0,-1,4]]$ of the previous
repository's saved result: a factor 3.1, the same order as the untwisted
case.

**Degree six, C4$\times\mathbb Z$** (750 s warm in the same report). The
profile of the extension worker (sampled every 10 ms of CPU, four worker
processes, the frozen benchmark source) shows no pair primitive at all;
the degree-six A phase already used the unary table. The time is the
one-fiber $V_3$ chain itself: 71 to 79 percent of all samples lie inside
`r3_chain.Htot` (the contraction of the degree-three model:
`product_shuffle`, `_rect_pull`, `_u2_pull`, `_v_pull`), reached from the
A-step of the lower product and from the potential $\Omega_6$; the self
time is 48 percent in the chain operators and 24 to 28 percent in the
degree-three model, and only 7 to 8 percent in the interval-cut leaves.
So in degree six the nesting that costs is $\Omega_6\to V_3\to$ chain of
$4\,000$ to $37\,000$ terms per simplex (Section 9), and the leaf
$\theta_5$ is already cached by label. The simplification here is not
another diagonal: it is to stop building the chain. The growth note of the
previous repository proves that $V_n$ on a simplex is a quasi-polynomial
in the labels (for the one-parameter family it computed,
$V_3(m)=\tfrac54m+(\tfrac14,-\tfrac{17}8,\tfrac14,-\tfrac98)_{m\bmod4}$),
so $V_n$ is a finite table of residue classes and linear coefficients,
evaluated without any chain. The same table removes the $V_n$ leaves of
Section 10, which are what bounds the degree-five factor at 3.7.

**Space groups, package degree 4** (the bottleneck profile of four
spin-half space groups, 240 s samples each). The sampled GAP stacks sit
in the normalized bar homotopy applied to degree-three defining cochains
(70 of 80 samples, task `atom_curvature` of a B-over-D row), and the SG 5
case was dominated by the construction of the primary comparison tensors
(330 s of GAP time before the primary-comparison optimization, 21 s
after). Neither is a universal primitive: both are the transport
$\Lambda,\Pi,H$ between the resolution and the bar. In the natural
formulation ([natural_evaluation.md](natural_evaluation.md)) the primary
rows and the pure-C curvature are evaluated with the resolution's own
cup-$i$ structure and the transport does not exist, so these costs
disappear rather than shrink.

**B generators** (B-over-D rows; both states have $A=0$). The pair phase
is $e^{A0}_k(B,C;B',C')$ of the product page, Section 7, whose only
universal primitive is $\mathcal U^{A0}_k(B,B')$ on the two-fiber reduced
binary cocycle model. On the square $B'=B$ the source collapses exactly
as (P3) did, since $B+B=0$:

$$
S^{A0}_k(B,B)=-2\,\mathrm{bp}_k(B)+\mathsf h\Bigl(h^2(Q_DB,Q_DB)+E\bigl(h^D(B,B)\bigr)+h^2\bigl(0,\delta h^D(B,B)\bigr)\Bigr),
$$

a one-fiber source on the model of a single binary cocycle with
$\theta_{k-1}$ leaves, so $\mathcal U^{A0}_k(B,B)$ is a one-fiber primitive
plus constants on that model's small cells, by the same three-line
identity (P5). The $\Theta$-pair primitive $\mathcal U^\Theta_m(q,q)$ inside
$R^{\rm red}$ collapses the same way ($q+q=0$). Neither was measured as a
bottleneck in the reports, so this is a formula statement, not a timing.

**What does not collapse.** The binary pair primitive $\mathcal V_n(A,A)$
in the C layer of a square (Section 5), the curvature prisms on the gauge
cylinder, and the leaves $\theta_m$. The first is a chain with a cup-word
source, the second and third are one-fiber words; all three are
candidates for label tables, not for a further algebraic reduction.
