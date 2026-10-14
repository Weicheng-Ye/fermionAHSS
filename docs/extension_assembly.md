# Using the stacking model: assembling the group of a degree

This page records how the differential and product of
[differential.md](differential.md) and
[product_upper_k1_k3.md](product_upper_k1_k3.md) are used to turn the four E6
layers of a degree into one abelian group, and which relations never need
the stacking formulas at all. Notation: [conventions.md](conventions.md).

## 1. Layers and relations

For package degree $j$ the E6 cells on the line $p+q=j-3$ are

$$
A_j=E_6^{\,j-3,0},\qquad B_j=E_6^{\,j-2,-1},\qquad C_j=E_6^{\,j-1,-2},\qquad D_j=E_6^{\,j+1,-4},
$$

assembled from D upward (D, then C, then B, then A). Each layer generator
is marked by a cocycle lifted through the retained page maps. The group is
filtered, $D\subset\langle C,D\rangle\subset\langle B,C,D\rangle\subset G$,
and each step is an extension problem: for a generator $q$ of order $m$ of
the next layer and a flat lift $\widetilde q$ (a flat state of the model
whose leading component is the marked cocycle of $q$), the **relation**

$$
m\,\widetilde q=t\in H
$$

is measured by forming the $m$-th stacking power of $\widetilde q$ and
reducing it in the lower group $H$ (the layers below $q$). With $R_H$
presenting $H$ and $T$ the matrix of the rows $t_i$, the group is the
cokernel of

$$
R=\begin{pmatrix}R_H&0\\-T&\operatorname{diag}(m_i)\end{pmatrix}.
$$

Rows are recorded as full coordinate vectors in one identified
presentation of $H$ (two relations with the same versus independent lower
images give different groups), and Smith normal form gives the invariants.
Free generators (only in A and D) add columns but no rows and split off.

**Target layer.** Only the class of $t$ in $\operatorname{Ext}(\mathbb Z/m,H)=H/mH$
matters. A lower generator $e$ is *free for $q$* when $e\notin mH$. The
target layer is the lowest layer containing a free generator; the power is
measured only through it (products and reductions are layer-limited, so no
D-layer correction is evaluated when the target is B or C), the components
below are recorded as zero with a certificate, and a relation without a
free lower generator is $t=0$ with no measurement at all. Such a truncated
row is the relation of a lift shifted within the lower group; when a later
relation with a nonzero coefficient on that generator is measured through a
layer below, the earlier relation is re-measured through D.

**Reduction.** The power is reduced layer by layer from the top: at each
layer the component is written as a combination of the marked lower
cocycles plus a coboundary (over $\mathbb Z$ with the signed coboundary in A
and D, modulo two in B and C), the state is divided on the left by the
product of the lower lifts with those coefficients and by the boundary
state of the gauge made from the coboundary primitive, and the quotient
continues downward. The exact comparison that is finally verified is

$$
\text{measured power}=\operatorname{act}(g,\ \text{canonical lower product}),
\qquad \operatorname{act}(g,x)=\mathfrak d(g)\times x,
$$

with a gauge $g$ of the preceding degree whose support is searched in the
order $(0,0,0,D)$, $(0,0,C,D)$, $(0,B,C,D)$, $(A,B,C,D)$. The lower product
is formed in the recorded generator order; the engine never commutes
factors at the cochain level.

**Lifts.** A measured A generator's flat lift solves, over $\mathbb F_2$ in B
and C and over $\mathbb Z$ in D,

$$
\delta_sA=0,\qquad \delta B+P(A)=0,\qquad \delta C+\tau(A;B)=0,\qquad \delta_sD+J(A,B,C)=0,
$$

changing the binary defining choices when a later equation fails. A lower
generator needs a complete lift only when it enters the compared product
with a nonzero coefficient above the last measured layer; in that layer it
enters with its marked cocycle alone.

## 2. Relations that need no stacking: the adjacent pairs

A two-primary relation whose target layer lies right below its generator's
layer is the class of a **primary operation** of the marked cocycle,
evaluated with the cup-$i$ products of the resolution and projected into
the E6 cell of the target layer (which quotients by all incoming images).
With $D=\operatorname{Sq}^2+s\operatorname{Sq}^1+\omega$:

| Relation | Class of the target-layer component | Evaluation on the resolution |
| --- | --- | --- |
| $a$ over B ($k=4,5,6$, $m=2^e$) | $D(z)\in B_k$, where $\beta_sz=\tfrac m2[a]$ | solve $\delta_su=m\,a$ over $\mathbb Z$, put $z=\rho u$; $z\smile_{k-6}z+s\smile(z\smile_{k-5}z)+\omega\smile z$ |
| $b$ over C ($k=2,\dots,6$) | $(\operatorname{Sq}^1+s)\,b=\rho\beta_s(b)\in C_k$ | $\rho\bigl(\delta_s\widetilde b/2\bigr)$, or $b\smile_{k-3}b+s\smile b$ |
| $c$ over D ($k=1,\dots,6$) | $\rho^{-1}\bigl(\operatorname{Sq}^2c+s\operatorname{Sq}^1c+\omega c\bigr)\in D_k/2D_k$ | $x-2w$ with $x$ the $0/1$ lift of $D(c)$ and $\delta_sw=\delta_sx/2$; equivalently $2D_R+\widetilde{E_R}+\delta_s\widetilde Y_R/2$ with $E_R=c\smile_{k-3}c+\omega c$, $Y_R=c\smile_{k-2}c$, $\delta_sD_R=-\delta_s\widetilde{E_R}/2$ |

These classes are consequences of the product formulas: the $m$-th power
of a lift has B component $\binom m2\alpha(a,a)$ (a coboundary) plus the B
boundary $Q_D(\rho u)$ of the A gauge; the square of a B lift has C
component $\beta(0,b;0,b)=h^D(b,b)$; the square of a C lift has D component
$2D_c+\gamma_C(c,c)$ with $\gamma_C(c,c)=\widetilde{E(c)}+\tfrac12\delta_s\widetilde{\operatorname{Sq}^1c}$
up to an integral coboundary. A C-over-D row is exact only modulo $2D$; it
is the relation of a lift shifted within D. These three rows, the
prime-split and zero-local rows, and the three-primary rows are the
relations that need no bar transport; the others do
([natural_evaluation.md](natural_evaluation.md)).

Non-adjacent pairs involve the secondary operation attached to the relation
$D\circ D=\rho\operatorname{Dtilde}\circ(\operatorname{Sq}^1+s)$ (A over C: C
part; B over D) or the tertiary operation (A over D at the prime two). Their
cochains need the natural primitives of the stacking model (pair
primitives, prisms, phases); the notes give no cup-$i$ formula on the
resolution for them, and the previous implementation measured them in the
model or read them from fixed residues of transported defining data.

## 3. Localization at the primes

Layer generators have prime-power or infinite order; the binary layers B
and C are two-groups, so every odd-primary relation belongs to A. Rows are
recorded prime by prime, two first, each prime from D upward, and a row
keeps only its coordinates on generators of its prime and on free
generators (the dropped coordinates lie on D generators of order prime to
$m$, which vanish after localization). A relation whose lower group has no
generator of its prime and no free generator is $t=0$ without measurement.
The model of a relation is chosen by the prime $p$ of its generator:

- $p=2$: the complete four-layer model, or the primary-operation row of
  Section 2 for an adjacent target.
- $p\ge5$: no measurement, $m\widetilde q=0$. The rows $q=0$ and $q=-4$
  lie in different Adams summands of $ko_{(p)}$, so the window carries no
  $k$-invariant at $p$.
- $p=3$: the two-layer model of rows $q=0$ and $q=-4$ of
  [prime_three.md](prime_three.md). Below package
  degree five it splits; in degrees five and six the light row is
  $t_D=2\cdot3^{e-1}Y$ with $\rho_3Y=P^1_s\rho_3A$ ($Y=A\smile_sA\smile_sA$ in
  degree five; a cyclic diagonal for $P^1$ in degree six).

## 4. Evaluating the formulas on a supplied resolution

The previous implementation evaluated the cochain formulas on the
normalized bar resolution (or a cell comparison complex) and transported
to the resolution $R$ with chain maps $f:B\to R$, $g:R\to B$ and a homotopy
$h$ with $1-gf=\partial h+h\partial$ and, after normalization
($q=1-gf$, $u=qhq$, $h'=u\partial u$), the side conditions
$fh'=h'g=0$, $(h')^2=0$. On cochains, with $\Lambda=f^*$, $\Pi=g^*$,
$H=(h')^*$:

$$
\Pi\Lambda=1,\qquad 1-\Lambda\Pi=\delta H+H\delta,\qquad H\Lambda=\Pi H=0 .
$$

A native state $w$ is embedded triangularly, with $N_\ell$ the nonlinear
part of the differential in layer $\ell$,

$$
\Phi(w)_A=\Lambda w_A,\qquad
\Phi(w)_\ell=\Lambda w_\ell-H\,N_\ell\bigl(\Phi(w)_{<\ell}\bigr),\qquad
\kappa_R(w)_\ell=\delta w_\ell+\Pi\,N_\ell\bigl(\Phi(w)_{<\ell}\bigr),
$$

so that $\kappa_R(w)=0$ layer by layer makes $\Phi(w)$ flat. Products are
formed on the bar and reflected back ($\Pi$ of the result, with the gauge
$H$ of the residual), gauge actions reflect $\mathfrak d(\Phi(e))\times\Phi(\text{canonical})$,
and the global branch flags are tested by pairing with the comparison chains
$g(e_j)$ of every basis element (a cochain is zero when it vanishes on $R$).
Defining cochains of a closed source $z$ were built as
$P(z;r)=\Lambda r+Hz$ with $\delta_Rr=\Pi z$, so that $\delta P=z$ literally
and every branch flag is known by construction; for a primary source the
homotopy term was replaced by a binary diagonal comparison
$C=\Lambda c_R+K_D(b_R)$ with $\delta K_D(b_R)=Q_D(\Lambda b_R)+\Lambda D_Rb_R$.

Two assumptions are built into the completed groups and are not proved
from the cochain identities: **gauge completeness** (native gauge
equivalence captures bar gauge equivalence; there is an algebraic
counterexample to deriving it from the retraction identities alone) and
the **abelian quotient** (stacking is commutative and associative on gauge
classes).
