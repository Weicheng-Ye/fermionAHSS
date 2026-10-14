# Evaluating the formulas on an arbitrary resolution

The formulas of these notes are written for normalized ordered simplicial
cochains. The previous implementation evaluated every nonlinear formula on
the normalized bar resolution and moved states to and from the supplied
resolution $R$ with comparison maps ("bar transport"). This page goes
through every site where that transport is used, says which formulas can be
evaluated directly on $R$ with its own structure ("natural"), rewrites those
in that form, and states exactly what the others need. Sections 7 and 8 do
the same for the secondary and tertiary operations and complete the
picture for the extension formulas.

The short version:

- everything **primary** is natural for any cup-$i$ structure on $R$ (Section 3);
- the secondary and tertiary operations and the stacking phases are natural
  for **any ordered simplicial model** of $BG$, and on a non-simplicial $R$
  need a *homotopy package* (Section 7);
- the **universal primitives** ($V_1,V_2,V_3$ of `T`, and $\mathcal V_n$,
  $\mathcal U$, $\mathcal P_k$, $\mathfrak q$ of the product) are evaluated
  through a simplicial section and cannot be replaced by a local solve on
  $R$ (Section 5); they are the extension data.

## 1. What structure a resolution can carry

Let $R$ be a free $\mathbb ZG$-resolution of finite type with a contracting
homotopy. Cochains $C^n(R;\mathbb F_2)=\mathbb F_2^{r_n}$, $C^n(R;\mathbb Z_s)=\mathbb Z^{r_n}$
are coordinate vectors on the basis of $R_n$; the twisted coboundary is the
matrix $M_s$ built with the sign character. The following can be
constructed on $R$ by acyclic carriers, inductively along the contracting
homotopy:

| Structure | Definition | Status in the old repo |
| --- | --- | --- |
| **Steenrod cup-$i$ structure** | equivariant $D_i:R\to R\otimes R$, $i\ge0$, $\partial D_i+D_i\partial=(1+T)D_{i-1}$, $D_{-1}=0$, $D_0$ a diagonal approximation; $x\smile^R_iy=(x\otimes y)\circ D_i$ | `cupMod2` (HAP), `native_coherence.gi` |
| integral cups with signs and transport | $D_0$ with the sign (C5) and the transport (C2) in the coefficient action | `nativePrimary`, `transported` |
| **cyclic diagonals** $D^{(3)}_0..D^{(3)}_3$ | $R\to R^{\otimes3}$, $\partial D^{(3)}_i-(-1)^iD^{(3)}_i\partial=\mathcal R_iD^{(3)}_{i-1}$ | stages $0,1,2$ built; stage 3 not |
| a **cylinder** $R^I$ | $i_0,i_1:R\to R^I$, $\pi:R^I\to R$, ramp $u\mapsto u\ell$, prism $I$ with $\delta_sI+I\delta_s=i_1^*-i_0^*$, $I\pi^*=0$, $(u\ell)|_{i_0}=0$, $(u\ell)|_{i_1}=u$, cup-$i$ structures restricting to $R$'s at both ends | not built |
| **Cartan and Adem homotopies** | four-input $\zeta^R_{1,n},\zeta^R_{2,n}$ and quadratic $\chi^R_n$ with the identities of Section 7 of [secondary_operations.md](secondary_operations.md) | not built |
| **universal sections** | for each universal model $M$, a map from the simplices of the model of $X$ to $M$ with the reconstruction identities ([universal_primitives.md](universal_primitives.md), Section 4) | exist only for simplicial models |

A formula is **natural on $R$** if it is built from these by identities that
hold for every such structure. Three kinds of identity occur:

- **generic identities**, valid for any Steenrod structure (Section 2);
- **interval-cut identities**, true for the cups on simplicial cochains and
  false on $R$ except up to coboundaries: $x\smile_nx=x$ for $|x|=n$,
  $e(a)=a\smile_{n-1}a$ literally, the vanishing of a word with a repeated
  degree-zero input, the "two top vertices" lemma of `closed_prism`;
- **universal identities**: $\delta V=S$ for a universal primitive, proved on
  the universal model and transported along a simplicial section.

## 2. Identities that hold for every Steenrod structure

Modulo two, with $\smile_i=\smile^R_i$ and $s,\omega$ closed:

- (G1) $\delta(x\smile_iy)=\delta x\smile_iy+x\smile_i\delta y+x\smile_{i-1}y+y\smile_{i-1}x$ (the cochain form of $\partial D_i+D_i\partial=(1+T)D_{i-1}$).
- (G2) bilinearity of every $\smile_i$.
- (G3) $\delta E(x)=E(\delta x)$, $\delta Q_D(x)=Q_D(\delta x)$. *Proof:* by (G1), $\delta(x\smile_{r-2}x)=\delta x\smile_{r-2}x+x\smile_{r-2}\delta x$ and $\delta(x\smile_{r-1}\delta x)=\delta x\smile_{r-1}\delta x+x\smile_{r-2}\delta x+\delta x\smile_{r-2}x$; the sum is $Q^2(\delta x)$; $\omega,s$ closed.
- (G4) $\delta h^j(x,y)+h^j(\delta x,\delta y)=Q^j(x+y)+Q^j(x)+Q^j(y)$ for $h^j(x,y)=x\smile_{r-j+1}y+\delta x\smile_{r-j+2}y$. *Proof:* (G1), (G2); the cross terms of $Q^j(x+y)$ are $x\smile_{r-j}y+y\smile_{r-j}x+x\smile_{r-j+1}\delta y+y\smile_{r-j+1}\delta x$.
- (G5) $E(u+v)+E(u)+E(v)=\delta(u\smile_{r-1}v)$ for closed $u,v$; $Q_D(x+y)+Q_D(x)+Q_D(y)=\delta h^D(x,y)+h^D(\delta x,\delta y)$.
- (G6) the carry identity (C0) is coordinatewise on any basis.
- (G7) the signed Leibniz rule for $\smile^R$ with transport; the Hirsch identity $A'\smile A-A\smile A'=\pm\delta_s(A\smile_1A')$ for even-degree cocycles.
- (G8) for closed $x$: $h^D(x,x)$, $x\smile_{r-1}x+s\smile x$ and $\rho(\delta_s\widetilde x/2)$ represent the same class $(\operatorname{Sq}^1+s)[x]$; on $R$ they differ by coboundaries.

## 3. The natural formulas

Throughout, $\smile_i$ means $\smile^R_i$, lifts and carries are
coordinatewise on the basis of $R$, and classes are projected to E6 cells
with `cell.project`.

### 3.1 Adjacent relation rows (no stacking at all)

$$
\boxed{\;
\begin{aligned}
&\text{A over B } (k=4,5,6,\ m=2^e):\ \ \delta_su=m\,A\ \text{over }\mathbb Z,\ z=\rho u,\qquad x_B=\bigl[z\smile_{k-6}z+s\smile(z\smile_{k-5}z)+\omega\smile z\bigr]\in B_k,\\
&\text{B over C } (k=2,\dots,6):\qquad x_C=\bigl[\rho(\widetilde BM_s/2)\bigr]=\bigl[B\smile_{k-3}B+s\smile B\bigr]\in C_k,\\
&\text{C over D } (k=1,\dots,6):\qquad E_R=C\smile_{k-3}C+\omega\smile C,\ Y_R=C\smile_{k-2}C,\ \delta_sD_R=-\tfrac12\delta_s\widetilde{E_R},\qquad x_D=\bigl[2D_R+\widetilde{E_R}+\tfrac12\delta_s\widetilde{Y_R}\bigr]\in D_k/2D_k .
\end{aligned}\;}
$$

These are the classes $D(z)$, $(\operatorname{Sq}^1+s)b$, $\rho^{-1}(\operatorname{Sq}^2c+s\operatorname{Sq}^1c+\omega c)$,
independent of the cup-$i$ structure, derived from the product formulas
(B layer of the $m$-th power $\binom m2\alpha(A,A)$ plus the gauge boundary
$Q_D(\rho u)$; C layer $h^D(b,b)$; D layer $2D_c+\gamma_C(c,c)$) and (G8).

### 3.2 The pure-C sector

For $A=B=0$ in every degree, with $E(C)=C\smile_{r-2}C+C\smile_{r-1}\delta C+\omega\smile C$, $h^2(x,y)=x\smile_{r-1}y+\delta x\smile_ry$:

$$
\boxed{\;
g(C)=\tfrac12\bigl(\delta_s\widetilde{E(C)}-\widetilde{E(\delta C)}\bigr),\qquad
\gamma_C(C,C')=\tfrac12\Bigl(\widetilde{E(C)}+\widetilde{E(C')}-\widetilde{E(C+C')}+\delta_s\widetilde{h^2(C,C')}+\widetilde{h^2(\delta C,\delta C')}\Bigr),\;}
$$

integral by (G3)–(G5), with $\delta_sg(C)+g(\delta C)=0$ and
$\delta_s\gamma_C(C,C')=g(C)+g(C')-g(C+C')+\gamma_C(\delta C,\delta C')$. The
rephasing $K(C)=\delta_s\mathsf h(s\smile C)$ and the C-shift carry
$I_\eta=\tfrac12(\widetilde{E(C+\delta\eta)}-\widetilde{E(C)}-\delta_s\widetilde v_\eta)$,
$v_\eta=E(\eta)+h^2(C,\delta\eta)$, are natural by the same identities. In the
native marking of the old audit a C generator with $\delta_sD_R=-\beta^R_sE_R(c)$
has $[2\widetilde c]=[2D_R+\widetilde{E_R(c)}+\beta^R_s\operatorname{Sq}^1c]$ in $D_k$.

### 3.3 The B-layer correction and the lower identity

$$
\boxed{\;\alpha_k(A,A')=h^D(a,a'),\qquad P(A+A')=P(A)+P(A')+\delta\alpha_k+\alpha_{k+1}(\delta_sA,\delta_sA')\;}
$$

by (G5); the B layer of every product and power is natural, as are
$\gamma_1$, $\beta_2$, $g_2$ and $\Phi_2$ on closed $B$ (cup-$0$, cup-$1$, lifts).

### 3.4 The model at the prime three

With the integral transported cups and the cyclic diagonals of $R$:

$$
\boxed{\;
\begin{aligned}
k=5:&\quad \gamma(A,A')=-2(A\smile A\smile A'+A\smile A'\smile A'),\qquad
x_D=2\cdot3^{a-1}A\smile A\smile A+\tfrac23\delta_s\Bigl[\sum_{j=1}^{3^a-1}\Xi^R(jA,A)-u\smile\delta_su\smile\delta_su\Bigr],\\
k=6:&\quad P^1_s\rho_3A=\langle a\otimes a\otimes a,D^{(3)}_2\rangle,\quad
\gamma=\tfrac23\bigl[\mathrm LP^1a+\mathrm LP^1a'-\mathrm LP^1(a+a')+\delta_s\mathrm L\varphi\bigr],\quad
D_u=-\tfrac23\bigl[\mathrm LP^1\rho_3\delta_su-\delta_s\mathrm L\chi(u)\bigr],
\end{aligned}\;}
$$

$\Xi^R$ the Hirsch primitive of (G7), $\varphi,\chi$ the words of
[prime_three.md](prime_three.md) on $R$'s $D^{(3)}_3,\dots,D^{(3)}_0$. Natural
once $D^{(3)}_3$ is built (the old engine stopped at stage 2).

### 3.5 Small natural facts

In degree four the torsion of $H^1(X;\mathbb Z_s)$ is $\{0,[\widetilde s]\}$, on
$R$ represented by $-\delta_s\widetilde1/2$ with relation primitive $U=-1$ and
B row $[\omega]$; for an A gauge $u$ with $\delta_su=mA$ the bar primitive is
$\Lambda u$ literally and its B boundary $Q_D(\rho u)$ is evaluated on $R$;
the degree $-1,0,1$ differentials and products use cup-$0$ and lifts only.

## 4. The audit

Status: **cup** = natural with a Steenrod structure only; **pkg** = natural
on a non-simplicial $R$ once the homotopy package and cylinder of Section 7
exist; **univ** = needs a universal section (simplicial model or transport).

| Site | Used for | Needs | Status |
| --- | --- | --- | --- |
| adjacent rows A–B, B–C, C–D | relations with adjacent target | (G1)–(G8) | **cup** (3.1) |
| pure-C $g,\gamma_C,K,I_\eta$ | C generators, pure-C references | (G3)–(G5) | **cup** (3.2) |
| $\alpha_k$, lower identity, B layer of powers | every product | (G5) | **cup** (3.3) |
| $\beta_2,\gamma_1$, legal $\Phi_2,g_2$ | degrees 1, 2 | cup-$0$, cup-$1$ | **cup** |
| three-local $\gamma,x_D,J,D_u$ | three-primary A over D | cups, $D^{(3)}_{\le3}$ | **cup** (3.4) |
| $g,k_0,\lambda$, $\tau'=Q_D(B)+k_0+\delta\lambda$, $f^\sharp$, $\beta^\sharp$ | C layer of lifts and products | $\zeta,\chi$ words, carries | **pkg** (7.1) |
| $\Psi$'s cochain, the matched lift $M'$ | `Psi`, B atoms | the same | **pkg** |
| $H_k$, $G_k$, $\beta^{\rm raw}_k$, $\widehat g_k$, $\widehat\gamma_k$ | off-legal branches | prisms of the degree above | **pkg** + the degree above |
| $\Theta_m$, $\mathrm{bp}_k(B)$, $\Omega_k(0,B,C)$ | B atoms, B-over-D rows | $\zeta,\chi$, prism | **pkg** |
| $R_0$ (even/odd), $\Omega_3$, $e^{\rm ev}_3,e^{\rm od}_3,e^{\rm pap}_3$ | degree 3 | words, prisms, and the face polynomial $\chi^{\rm rep}$ | **pkg** (the repair needs a simplicial reading or a cup-word rewrite) |
| $\widehat R_1,\widehat R_2,\widehat R_3$ without $V_n$ | $\Omega_k$, $J_k$, $k\ge4$ | words, prisms, cubes | **pkg** |
| $V_1,V_2,V_3$ | inside $\widehat R_n$ | universal sections | **univ** |
| $\beta^L_k$ ($k\ge4$): $\mathcal V_n$ | C layer with nonzero A | universal sections | **univ** |
| $\mathcal V_0$ | degree-3 C layer | face polynomial | **pkg** (cup-word rewrite of a face polynomial) |
| $e^{\rm prod}_k$, $e^{A0}_k$: $\mathcal P_k$, $\mathcal U^\Theta_m$, $\mathcal U^{A0}_k$ | D layer of products | universal sections | **univ** |
| $\xi_k$, $q_{\rm loc}$, $q_3,q_4,q_5$, $\lambda_3$ | branch comparisons | $v_n$, $\mathfrak q$'s | **univ** |
| section branch of degree six | strict identity off $L_6$ | a global section | not needed for the groups |

## 5. Why a universal primitive cannot be replaced by a local solve

Two primitives of the same source differ by a cocycle $z$. For the binary
primitives ($\mathcal V_n$, $v_n$, $\lambda_3$) the C layer of a product
changes by $z(A,A')$, so the C coordinate of a relation changes by the class
$\sum_jz(jA,A)$; for the rational ones ($V_n$, $\mathcal P_k$, $\mathcal U$,
$\mathfrak q$) the phase changes by a $\mathbb Q/\mathbb Z$-cocycle and the
integral output by its Bockstein. These classes are the extension data and
the $d_5$ values the computation exists to find.

On the **universal model** the ambiguity is harmless: it is a finite group
(for $V_3$, $(\mathbb Z/2)^3\oplus\mathbb Z/12$), and the calibration selectors
compensate it exactly (Section 8.5 of [tertiary_operations.md](tertiary_operations.md)).
On a **specific** $X$ the ambiguity of a solution of $\delta_sV\equiv\phi$ is
all of $H^{n+4}(X;\mathbb Q/\mathbb Z_s)$, and nothing on $X$ singles out the
pullback of the universal choice. This is why the primitives are evaluated
through a universal section, and why no natural operation on $R$ can
replace them: $R$ has no map to the universal model (a free resolution is
not a simplicial set), and the comparison with a simplicial model *is* the
bar transport. Two further facts:

- **Pushing a primitive up a degree does not remove it.** The prism
  constructions express degree $k$ through degree $k+1$ on the cylinder, but
  the ramp of a flat state is not flat, so $\Gamma^{\rm norm}_{k+1}$ is needed
  on nonzero curvature with its own $\mathcal P_{k+1}$; $\mathcal P_6$ is irreducible.
- **The calibration lives on the universal model.** `chi7_tail`, $\epsilon,\eta$,
  the tertiary selectors and the small coefficients were fixed by suspension
  and Thom comparisons there; a package on $R$ computes the same operations
  only after it has been compared with the universal one.

## 6. What this changes for production

The natural sites of Section 3 are those the old engine already evaluated
on $R$ (primary rows, light rows, three-primary rows). The remaining
transport in a production run is the B-atom curvature, the exact B-over-D
residue, A over C and two-primary A over D, i.e. exactly the **pkg** and
**univ** rows. A transport-free production path therefore has two routes:

1. **Compute on a simplicial model of $BG$.** Every formula in these notes
   is natural for any ordered simplicial set (normalized cochains with
   interval cuts); the normalized bar is one, the cell comparison complex of
   the old repo another, and any $G$-triangulation of a contractible free
   $G$-complex with a $G$-invariant vertex order gives one with cells
   $C_*(E)$, cochains $\operatorname{Hom}_G(C_*(E),\cdot)$ and universal sections
   from the face values. Nothing is deferred and nothing is transported; the
   price is the number of simplices ($(|G|-1)^n$ nondegenerate $n$-simplices
   for the bar).
2. **Compute on a HAP resolution** with the natural part of Section 3
   native and the **pkg** part made native by the package of Section 7;
   only the **univ** rows then need transport, and only on the phases that
   contain a universal primitive.

## 7. The secondary and tertiary operations on a resolution

### 7.1 What `Tau` and `Psi` use

Besides (G1)–(G6), the formulas (S1)–(S14) of
[secondary_operations.md](secondary_operations.md) use:

- (W1) the Cartan identities $\delta\zeta_{1,n}(x,a)=\operatorname{Sq}^2(xa)+x\operatorname{Sq}^2a+x^2e_a$, $\delta\zeta_{2,n}(x,a)=\operatorname{Sq}^2(xa)+x\operatorname{Sq}^2a+x^2a+e_xe_a$;
- (W2) the Adem identity $\delta\chi_n(a)=\operatorname{Sq}^2\operatorname{Sq}^2a+e_a\smile_{n-2}e_a$, $\chi_n(0)=0$;
- (L1) the literal identity $e(a)=\rho(\delta\widetilde a/2)=a\smile_{n-1}a$ for a binary cocycle, which makes $D(a)=u+v+w$ and ties the carries $B,C_B$ to the cup words;
- (L2) the vanishing of a word with a repeated degree-zero input (needed only for $n=0$: $\zeta_{i,0}=0$);
- (X) the exactness of the divisions in $q$, $L$, $\psi'$, $L_A$, proved in the unbundled notes from (W1)–(W2) and (L1).

On a non-simplicial $R$: (W1)–(W2) become **requirements** on the package.
Equivariant homotopies $\zeta^R_{1,n},\zeta^R_{2,n},\chi^R_n$ with these
coboundaries exist and can be built by acyclic carriers (the right-hand
sides are null-homotopic because the Cartan and Adem relations hold in
cohomology, and $R^{\otimes4}$ is acyclic in positive degrees). (L1) fails;
the consistent replacement is the device of the old light-transport note:
keep $e(a)=\rho(\delta\widetilde a/2)$ wherever an integral carry needs it,
use the cup-defined $\operatorname{Sq}^1a=a\smile_{n-1}a$ wherever a word needs
it, and insert a chosen $t_R$ with $\delta t_R=a\smile_{n-1}a+e(a)$ (a
primary comparison) into $w$, $F$ and $G$ so that (S2) and the closedness
of $F$ hold on $R$; (X) must then be re-verified for the modified formulas.
With those four steps `Tau` and `Psi` are evaluable on $R$ as cochain
formulas. What is **not** automatic is that they are the calibrated
operations: a package built by acyclic carriers is one system of
null-homotopies among many, and on $R$ the difference between two systems
is not a natural operation. Matching it to `chi7_tail` is the "calibrated
native comparison" that the old repo lists as its missing step; it has to
be done once per package construction (by a universal comparison, or by
agreement on test groups where the answer is known). The secondary
operations are therefore **pkg**: natural as formulas on $R$, with a
one-time calibration.

### 7.2 What `T` uses in addition

$\mathcal O_n$ uses, beyond the secondary data: prisms of $\Theta$ on the
ramped defining cochain (a cylinder), the explicit cup and carry terms
$\operatorname{Pol}_n$, $E(\lambda)$, $L_n,N_3,O_3,M_3$, the signed cube and
$P^1_s$ (cyclic diagonals), and the source primitives $V_1,V_2,V_3$. All but
the last are **pkg** (the cylinder has to carry a cup-$i$ structure
compatible with both ends; the prism of a word of ramps is then a word on
$R^I$). $V_n$ is **univ**: it is the value at the universal section of a
contraction on the universal model, with the small coefficients of Section
8 of [tertiary_operations.md](tertiary_operations.md). On a simplicial model
it is natural and explicit; on $R$ it is transported.

### 7.3 Consequences

- On a simplicial model of $BG$ the entire system (secondary, tertiary,
  stacking) is natural and explicit: no transport, no package, no
  calibration beyond the fixed constants of the notes.
- On a HAP resolution, `Tau`, `Psi`, the explicit parts of `T` and all
  **pkg** rows of the stacking become native once the package exists and
  is calibrated; the **univ** rows ($V_n$ and the stacking primitives) stay
  transported.

## 8. The extension formulas, completed

With [secondary_operations.md](secondary_operations.md) and
[tertiary_operations.md](tertiary_operations.md) in place, every object in
the stacking pages is written out: the C layer through $\tau'=Q_D(B)+k_0+\delta\lambda$
and $f^\sharp$, $\beta^\sharp$, $H_k$; the phases $\Omega_k$ through
$\widehat R_n$ and $V_n$; the pair phases through $R^{(k)}_A$, $Q_k$ and the
primitives of [universal_primitives.md](universal_primitives.md). The only
inputs that are not closed cochain formulas are the universal primitives,
each specified by its source, model and coefficients.

For the group computation this gives one explicit rule: **the relation of a
generator $q$ of order $m$ is the $m$-th power of its flat lift under
$\times_k$, reduced in the lower group through the target layer**
([extension_assembly.md](extension_assembly.md), Section 1), with the
shortcuts of Section 3.1 for adjacent targets and the two-layer model at
the prime three. The old engine's further "light-row" residues for the
non-adjacent relations are optimizations of this rule that evaluate the same
product on transported defining data; their derivations are not bundled, and
the rule above does not depend on them.

Finally, the status of each extension formula on a resolution is the status
of its ingredients in the audit of Section 4: natural (**cup**) for the
adjacent rows, the pure-C sector, $\alpha_k$ and the prime-three model;
natural given the package (**pkg**) for every secondary-level cochain and
every prism; universal (**univ**) for the primitives. There is no site whose
status is unknown.
