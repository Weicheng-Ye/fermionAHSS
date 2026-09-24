# Fixed universal helpers used by the implemented operations

This reference records the formulas and finite coefficient data used by the
current production evaluators. Read [conventions.md](conventions.md) for the
coefficient systems and [secondary_operations.md](secondary_operations.md)
and [tertiary_operations.md](tertiary_operations.md) for the assembled
differentials. All statements about the current implementation refer to the
`chi7_tail`, epsilon `(1,0,0)`, eta `(1,0,1)` convention.

Binary expressions are reduced modulo two before their indicated binary
lift. Rational expressions are evaluated in **Q**, using the displayed
representatives, and only then interpreted as phases in Q/Z. In particular,
an equation such as `d_s V = phi` below is an equation of phases unless
explicitly described as an equality of rational cochains. Replacing a
displayed source value by another lift before dividing can change the
normalization.

## 1. Interval-cut words and Cartan helpers

The executable definitions are
[natural_words.gi](../gap/natural_words.gi) and
[cochain_tools.py](../python/cochain_tools.py).
For a nondegenerate surjection word `w=(w_1,...,w_L)` on `1,...,r`, and
cochains of degrees `d_1,...,d_r`, its output degree is

\[
N=\sum_jd_j-L+r.
\]

On an ordered `N`-simplex sum over cuts
`0=i_0 <= i_1 <= ... <= i_L=N`. For label `j`, concatenate the intervals
`[i_{ell-1},i_ell]` with `w_ell=j`, in their order. Retain a cut only when
this concatenation is strictly increasing and contains exactly `d_j+1`
vertices for every label. Its contribution is the product of the input
cochains on those faces. Sum modulo two. Negative output degrees and
infeasible cuts contribute zero. This specifies `word_op` without any
choice of a cochain primitive.

The alternating word of length `i+2`, beginning with `1`, implements
`cup_i`. Integral `cup_i` uses the same cuts, with sign

\[
(-1)^{i(d_1+d_2)+\binom i2+
 \sum_{\ell\text{ nonfinal}}i_\ell+
 \sum_{\ell<t,\,w_\ell>w_t}\mu_\ell\mu_t},
\qquad
\mu_\ell=i_\ell-i_{\ell-1}
 +\mathbf1_{\ell\text{ nonfinal}}.
\]

Here “nonfinal” means that the same label occurs later in the word.
Transport for an actual sign-valued factor is separate from this ordinary
integral interval-cut sign; see the conventions reference.

For `n>=1`, the exact Cartan helpers are

\[
\begin{aligned}
\zeta_{1,n}(x,a)
 &=\operatorname{word}_{1232\,\mathrm{alt}_n(4,3)}(x,x,a,a),
 &&|x|=1,\ |a|=n,\\
\zeta_{2,n}(x,a)
 &=\operatorname{word}_{1231\,\mathrm{alt}_{n+1}(3,4)}(x,x,a,a),
 &&|x|=2,\ |a|=n.
\end{aligned}
\]

`alt_l(u,v)` means `l` letters alternating `u,v`, beginning with `u`.
Their output degrees are `n+2` and `n+3`, respectively. Both helpers are
defined to be zero when `n=0`. For example, their first words are
`12324` and `123134`.

For binary cocycles, set `e_y=rho_2(d tilde y/2)` and similarly `e_x`.
Their actual cochain boundary identities are

\[
\begin{aligned}
d\zeta_{1,n}(x,a)&=Sq^2(xa)+xSq^2a+x^2e_a,\\
d\zeta_{2,n}(x,a)&=Sq^2(xa)+xSq^2a+x^2a+e_xe_a.
\end{aligned}
\]

The degree of `x` is respectively one or two. These identities use the
ordinary binary cocycle differential and the displayed interval-cut
convention.

## 2. The exact finite chi formula

Production GAP and Python both read the **same** compiled data file:
[chi-calibrated-degree7-anf.g](../data/chi-calibrated-degree7-anf.g).
It specifies a finite algebraic normal form, rather than calling a solver.
For a binary degree-`n` cochain `a`, let `F_1,...,F_M` be all ordered
`(n+1)`-vertex faces of the `(n+3)`-simplex, in lexicographic order, and set
`x_j=a(F_j)`. The implemented formula is

\[
\boxed{\quad
\chi_n(a)=\sum_{c\in\mathcal C_n}\prod_{j\in J(c)}x_j\pmod2.
\quad}
\]

The complete explicit coefficient list `C_n` is the `monomials` array of
the degree-`n` record. To decode an integer `c`, expand it in radix
`2^bitsPerFace`; its nonzero digits are the **one-based face indices** in
`J(c)`. GAP stores vertex positions starting at one. Python subtracts one
from those vertex positions, but not from the conceptual ordering of the
faces. Products are square-free because `x_j^2=x_j` over F2. Thus the data
file and this decoding rule specify every coefficient; the large word
lists are not additional runtime inputs.

| Input n | Output degree | Faces | Bits per face | Original words | Nonzero ANF monomials |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 3 | 4 | 7 | 0 | 0 |
| 1 | 4 | 10 | 7 | 0 | 0 |
| 2 | 5 | 20 | 7 | 24 | 33 |
| 3 | 6 | 35 | 7 | 267 | 496 |
| 4 | 7 | 56 | 7 | 2,896 | 5,680 |
| 5 | 8 | 84 | 7 | 22,147 | 32,720 |
| 6 | 9 | 120 | 7 | 160,350 | 149,704 |
| 7 | 10 | 165 | 8 | 1,248,069 | 516,329 |

The current family comes from
tail_words.json (source-workspace provenance: `note/extra/chi_suspension_degree7/tail_words.json`; not bundled),
with its fixed suspension comparison. The optional `head` and `tail6`
branches still visible in `koAHSSNaturalChiValue` are separate diagnostic
families. They are not the default or interchangeable with this table.
For binary cocycles the boundary identity is

\[
d\chi_n(a)=Sq^2Sq^2a+e_a\cup_{n-2}e_a,
\qquad e_a=\rho_2(d\widetilde a/2),\qquad\chi_n(0)=0.
\]

## 3. The prescribed boundary phase Theta

The rational-phase helper \(\Theta_m\) is separate from the differential `Psi`.
The exact helper is `theta` in
[phase_eval.py](../python/phase_eval.py). Write

\[
Q^i(x)=\rho_2\bigl(x\cup_{|x|-i}x+
             x\cup_{|x|-i+1}\widetilde{dx}\bigr),\qquad
E(x)=Q^2(x)+\omega x,\qquad Q_D(x)=E(x)+sQ^1(x).
\]

The argument to `d` inside `Q^i` is reduced modulo two; `x` is binary.
For a cocycle, `Q^i(x)=Sq^i(x)` in the cochain convention used here.
Let `q` be a binary cocycle of degree `m`, and let `v` have degree `m+1`
with `dv=Dq`. On a binary cocycle the chosen interval-cut convention has
the exact identity `Sq^1q=q cup_(m-1)q=rho_2(d tilde q/2)`. Thus the literal
primary value is `Dq=Sq^2q+omega q+s e_q`. Define the integral and binary
carries

\[
B=\frac{d\widetilde q}{2},\qquad
e=\rho_2B,\qquad C_B=\rho_2\frac{B+\widetilde e}{2},
\qquad u=Sq^2q,\quad w=\omega q,\quad z=se.
\]

The **plus** in `B+e` is part of the convention. Put

\[
\begin{aligned}
H_m(q)={}&\zeta_{2,m}(\omega,q)+\chi_m(q)
 +u\cup_{m+1}w+u\cup_{m+1}z+w\cup_{m+1}z\\
 &+\zeta_{1,m+1}(s,e)+(\omega\cup_1s)e
 +su+s^2C_B,\\
Z_m(q)={}&sSq^2q+\omega Sq^1q.
\end{aligned}
\]

All terms in `H_m` and `Z_m` are binary. Then

\[
\boxed{\displaystyle
\Theta_m(q,v)=
\frac12\widetilde{E(v)+H_m(q)}
 +\frac14\left(\widetilde\omega\cup B+
                         B\cup_{m-1}B\right)
 +\frac12\widetilde{Z_m(q)}.}
\]

Its degree is `m+3`. The final term is exactly eta `(1,0,1)`.
The products involving `B` in the quarter-valued term are integral
interval-cut products. `Theta_m` is a fixed rational lift of the phase;
it is not the result of division in Q/Z.
In particular, the last term is one binary lift of the **sum**
`s Sq^2q+omega Sq^1q`. The secondary GAP assembler's separately lifted
correction summands can give a different rational lift by an integral
cochain; the Python `theta` above is the lift actually used by the
tertiary source contractors.

The universal A-only source used for `V_n` is
`phi_n=Theta_{n+2}(q_n,k_n)`, of degree `n+5`. The source splitting
`q_n,k_n` is given with the tertiary assemblers. For `n=2,3`, `q_n=p=Da`;
for `n=1`, the adjusted source is `q_1=(omega+s^2)a`, not the unadjusted
`p`. The general finite primitive `V_n` has degree `n+4`.

## 4. Fixed chain comparison conventions

The contractors below are prescribed operators on universal simplicial
models. They never solve a residual equation on the supplied group or
space. Their chain maps use exact integer coefficients.

For a product, the Alexander–Whitney map is the sum of front/back cuts;
the shuffle map sums monotone lattice paths with sign equal to the
parity of the inversions between horizontal and vertical steps. Degenerate
terms are discarded. The product homotopy in
[chain_models.py](../python/chain_models.py) is fixed
recursively, rather than selected from its homotopy class. On universal
vertices, with `Q=1-shuffle*AW`, set

\[
\widehat H_0=0,\qquad
\widehat H_n=C_{(0,0)}(Q-\widehat H_{n-1}\partial),\qquad
u=Q\widehat H Q,\qquad H=u\partial u.
\]

`C_(0,0)` prepends the initial vertex in both factors. The homogeneous
dihedral construction in the next section uses the specified raw product
homotopy `Hhat`; the matrix-model contractions use the normalized `H`.
For two contractions, the tensor homotopy and composite homotopy are

\[
H_{X\times Y}=H_{\rm prod}+
G_{\rm prod}(H_X\otimes1+G_XF_X\otimes H_Y)F_{\rm prod},
\qquad
H_{21}=H_1+G_1H_2F_1,
\]

with the Koszul sign `(-1)^|x|` on the second tensor-homotopy term.

For the signed two-cocycle matrix model, the small target has generators
`u_(p,k)` of degree `p+2k`, and

\[
\partial u_{p,k}=
\begin{cases}
-2u_{p-1,k},&p>0\text{ and }p+k\text{ odd},\\
0,&\text{otherwise}.
\end{cases}
\]

The background target has bar words `[x_(i1)|...|x_(il)]`, of degree
`sum(i_j+1)`. Its products and differentials are those of the normalized
bar of the C2 chain algebra, as implemented in `chain_models.py`.
The signed comparison first contracts the unit-edge chains to the crossed
exterior algebra

\[
\mathbb Z\{1,\sigma,e,e\sigma\},\qquad
\sigma^2=1,\quad e^2=0,\quad\sigma e=-e\sigma.
\]

Its reduced letters are `t=sigma-1,e,e sigma`, with `t^2=-2t`.
[exterior_bar.py](../python/exterior_bar.py) fixes the next
comparison by `f=h_P f partial`, `g=h_Q g partial` and
`K=h_Q(1-gf-K partial)`, with identities in degree zero. The free small
boundary and contractions are

\[
\begin{aligned}
\partial u_{p,k}&=(-1)^p e u_{p,k-1}
 +(\sigma+(-1)^{p+k})u_{p-1,k},\\
h_P(\sigma u_{p,0})&=u_{p+1,0},&
h_P(eu_{p,k})&=(-1)^p u_{p,k+1},&
h_P(e\sigma u_{p,k})&=(-1)^{p+1}\sigma u_{p,k+1},\\
h_Q(a[w])&=-[a-\epsilon(a)|w].
\end{aligned}
\]

Terms with a negative index and all unlisted `h_P` values are zero.
Finally the sign module sends `1,sigma,e,e sigma` to `1,-1,0,0`.
The unit-edge and diagonal-to-bar stages retain their finite perturbation
series; they are not replaced by a multiplicative guess.

## 5. V1 and the strictly normalized odd primitive

The executable definitions are `V1_pair`, `odd_base_normalization`, and
their helpers in [low_phases.py](../python/low_phases.py), with
[odd_primitive.py](../python/odd_primitive.py) and
[odd_comparison.py](../python/odd_comparison.py).

For signed integral `A` of degree one, the homogeneous dihedral vertices
attached to a simplex are

\[
g_i=(A_{0i},s_{0i})\in\mathbb Z\rtimes C_2,
\qquad g_0=(0,0),\qquad (k,\epsilon)\longmapsto 2k-\epsilon\in\mathbb Z.
\]

The last map identifies its Cayley line. The reconstruction is
`s_ij=epsilon_i+epsilon_j` and
`A_ij=(-1)^epsilon_i(k_j-k_i)`. The background uses the matrix section
specified in the next section.

On homogeneous group simplices, the retraction and homotopy are

\[
R[g]=K_{g_0}(R\partial[g]),\qquad
h[g]=C_{g_0}([g]-R[g]-h\partial[g]),
\]

starting with `R[g_0]=[g_0]`, `h[g_0]=0`. On a vertex, `K_root` is the
oriented unit-edge path in the integer Cayley line from `root` to that
vertex. On an alternating simplex in one unit edge, prepend its endpoint
nearest the root. Actual reverse bar edges remain different generators;
they are not collapsed by an unoriented-edge identification.

Let `P` be the product retraction using `R` in the group factor and the
identity in the background, and let

\[
K=\widehat H_{\rm prod}+
\operatorname{shuffle}(h\otimes1)\operatorname{AW}.
\]

The fixed source primitive is

\[
\boxed{V_1=\phi_1K+U_{\rm wedge}P.}
\]

On the branch with increment `(A,s)=(0,1)`, `U_wedge=0`; on the branch
with increment `(1,1)`, it is the following `U_t`. Every evaluation is
multiplied by `(-1)^epsilon_0` when the homogeneous first vertex is not in
the zero sign component.

For the literal odd branch `A=tilde s`, append a right-prism coordinate
`tau`, set `s_I=rho_2(s+d tau)`, `A_I=tilde s_I`, and form the corresponding
source `Phi_I`. Write `I` for the raw signed right-prism sum. Then

\[
\kappa=\rho_2 I(d_{s_I}\Phi_I),\qquad
U_{\rm raw}=-\tfrac12 I\Phi_I+\tfrac12\widetilde{\operatorname{Fill}(\kappa)}.
\]

All eight odd-source detector bits are zero; therefore the three possible
additional quarter-valued Bockstein corrections have zero coefficients.
`Fill` is fixed by

\[
\operatorname{Fill}(\gamma)(z)=
\gamma\bigl((H+Gh_{\rm sm}F)j(z)\bigr)\pmod2.
\]

Here `F,G,H` are the explicit product comparison for `BC2_s` times the
background matrix model. On the background word, the only nonzero
`h_sm` rules are

\[
[3]\mapsto[1|2],\quad [5]\mapsto[1|4],\quad
[1|3]\mapsto[1|1|2],\quad [3|1]\mapsto[2|1|1].
\]

The group degree is retained. The further **literal base normalization**
is essential:

\[
U_t=U_{\rm raw}-d_sE_0,
\qquad
E_0=(U_{\rm raw}|_{s=0})H_\omega+e_0F_\omega.
\]

The base period is zero, so its possible eta correction vanishes. The two
source values are `(B,C)=(7/8,0)`, and the complete sparse `e_0` table is

| Background basis | Literal rational value |
| --- | ---: |
| `[x_3]` | `-7/16` |
| `[x_1|x_1]` | `-7/32` |

All other values are zero. This makes the chosen `U_t` the zero phase
cochain on `s=0`, not merely a phase with zero cohomology class. Its
prescribed rational lift can still be integral there. Dropping `E_0`
would change the fixed source convention used by the general `V1`
evaluator.

## 6. V2: section, contraction, and complete coefficient table

Use two diagonal matrix models: signed integral rows `(sigma_i,M_i)` and
binary background rows `W_i`. Their reconstruction on a face `(r,l,t)` is

\[
s_{rt}=\sum_{i=r}^{t-1}\sigma_i\pmod2,\qquad
A_{rlt}=\sum_{i=r}^{l-1}\sum_{j=l}^{t-1}(-1)^{s_{ri}}M_{ij},\qquad
\omega_{rlt}=\sum_{i=r}^{l-1}\sum_{j=l}^{t-1}W_{ij}\pmod2.
\]

The section `j` takes `sigma_i=s_(i,i+1)` and, for `i<j`,

\[
\begin{aligned}
M_{ij}&=A_{i,i+1,j+1}-A_{i,i+1,j},
&M_{ji}&=-(-1)^{s_{ij}}M_{ij},\\
W_{ij}&=\omega_{i,i+1,j+1}+\omega_{i,i+1,j}\pmod2,
&W_{ji}&=W_{ij}.
\end{aligned}
\]

Diagonals are zero; repeated-vertex cochain values are zero. These formulas
give the strict retraction `rj=1`. With the fixed product contraction in
Section 4,

\[
\boxed{V_2=j^*\bigl((r^*\phi_2)H_{\rm tot}+e_fF_{\rm tot}\bigr).}
\]

In the small basis `u_(p,k) tensor [x_i|...]`, the five source evaluations
used by [source_primitive.py](../python/source_primitive.py) are

| Name | Degree-seven source basis | Value modulo one |
| --- | --- | ---: |
| C | `u_(2,1) tensor [x_2]` | `0` |
| G | `u_(0,1) tensor [x_1|x_2]` | `3/4` |
| H | `u_(0,1) tensor [x_2|x_1]` | `0` |
| A_src | `u_(3,2) tensor []` | `3/4` |
| B_src | `u_(1,2) tensor [x_1]` | `1/4` |

The **complete** degree-six `e_f` prescription is

| Basis | Formula | Literal rational value |
| --- | --- | ---: |
| `u_(2,1) tensor [x_1]` | `-C/2` | `0` |
| `u_(0,1) tensor [x_3]` | `(H-G)/2` | `-3/8` |
| `u_(0,1) tensor [x_1|x_1]` | `-(G+H)/4` | `-3/16` |
| `u_(2,2) tensor []` | `-A_src/2` | `-3/8` |
| `u_(0,2) tensor [x_1]` | `-B_src/2` | `-1/8` |

Every unlisted coefficient is zero. The current suspension normalization
uses `V2fin=V2-L2` where `L2=P(omega) cup_s A/4`. The R2 assembler separately
contains the adopted `-A^3/4` term. Do not put that cubic term into `V2fin`
when forming the R3 calibration source.

## 7. V3: full crossed transfer and fixed dyadic source

The implemented source contractor is
[r3_chain.py](../python/r3_chain.py) together with
[r3_source.py](../python/r3_source.py).
Its full fixed matrix, transformations, source evaluations and coefficients
are [r3_source.json](../python/r3_source.json).

The signed model is `diag B(U2 semidirect C2)`. A row is an ordinary
integral two-cocycle `z_i` with sign component `sigma_i`; multiplication is

\[
(z,\epsilon)(z',\epsilon')=(z+(-1)^\epsilon z',\epsilon+\epsilon').
\]

Its reconstruction and section are

\[
\begin{aligned}
rA(v_0,v_1,v_2,v_3)
 &=\sum_{i=v_0}^{v_1-1}(-1)^{s(v_0,i)}z_i(v_1,v_2,v_3),\\
\bar A(v_0,v_1,v_2,v_3)&=(-1)^{s(0,v_0)}A(v_0,v_1,v_2,v_3),\\
z_i&=(-1)^{s(0,i)}(\iota_i\bar A-\iota_{i+1}\bar A),
\qquad\sigma_i=s(i,i+1).
\end{aligned}
\]

Extend `Abar` alternatingly, and by zero on repeated vertices. Each `z_i`
is represented by its ordinary skew matrix, using the same strict
upper-triangle section as for `V2`. Together with the unchanged background
matrix section these define `rj=1` literally.

The ordinary `U2` contraction projects to divided-power generators `u_k`
of degree `2k`. Retain the C2 component: the augmented alphabet is

\[
t=(1,0)=\sigma-1,\qquad
u_k=(0,k),\quad v_k=(1,k)\quad(k\geq1).
\]

After bar suspension their degrees are `1,2k+1,2k+1`. This is only a
description of the graded module. The actual differential is obtained by
the **full** transfer, including all crossed multiplication and initial
sign-module action terms. For the suspended tensor contraction
`f_T,g_T,h_T` and all bar multiplication terms `delta`, it is

\[
\begin{aligned}
L&=\sum_{m\geq0}(-\delta h_T)^m,&
R&=\sum_{m\geq0}(-h_T\delta)^m,\\
F&=f_TL,&G&=Rg_T,&H&=h_TL,&
d_{\rm small}&=f_T\delta Rg_T.
\end{aligned}
\]

Each series terminates because the perturbation lowers bar length. The
preceding diagonal-to-bar stage has its own finite series
`sum(-d_h h_0)^m`; it is also retained. The local normalized homotopies
use `Q=1-ip`, `u=Q H_raw Q`, `h=u d u`. Tensor homotopies have the actual
suspended Koszul signs. This prescription avoids assuming that the
ordinary U2 contraction is multiplicative or C2-equivariant.

Tensor with the background word complex, and quotient by words containing
no `u_k` or `v_k`: that is the relative `A=0` quotient. A basis element is
written `[ell_1|...|ell_r] tensor [x_(j1)|...|x_(jt)]`. Order first by the
total signed-word degree, then signed-word length, then lexicographically
by `(component,k)`, followed by background degree, length and word. There
are 46 relative basis elements in degree seven and 94 in degree eight.

Let `D` be the resulting `46 by 94` boundary matrix and
`f=(r^*phi_3)G_tot` its degree-eight phase row. The deterministic Smith
recipe chooses the least nonzero absolute value in the remaining block,
breaking ties by row then column; it uses nonnegative Euclidean
remainders, promotes a nonzero remainder immediately, and repeats if the
current pivot does not divide a remaining entry. This gives

\[
UDV=\operatorname{diag}(1^{\times27},2,2,2,12,0,\ldots),
\qquad \operatorname{rank}D=31.
\]

All 63 transformed cycle values of `fV` are zero modulo one. For the
nonzero diagonal entries write `d_i=2^(v_i)o_i`, `o_i` odd, and define

\[
n_i=4\operatorname{rep}((fV)_i),\qquad
m_i=o_i^{-1}n_i\pmod4\quad(0\leq m_i<4),\qquad
e'_i=\frac{m_i}{2^{v_i+2}},\qquad e_f=e'U.
\]

Other `e'_i` are zero. This is a **dyadic** division prescription even
for the pivot `12`; arbitrary division by `12` would select a different
phase. Finally

\[
\boxed{V_3=j^*\bigl((r^*\phi_3)H_{\rm tot}+e_fF_{\rm tot}\bigr).}
\]

The following tables completely specify the sparse phase row and sparse
primitive row. In them `w_k` means either `u_k` or `v_k`, independently
where stated. Unlisted entries are zero.

| Degree-eight basis | f, modulo one |
| --- | ---: |
| `[w_1] tensor [x_4]` | `3/4` |
| `[w_1] tensor [x_1|x_2]` | `3/4` |
| `[w_2] tensor [x_2]` | `1/4` |
| `[t|t|w_1] tensor [x_2]` | `1/2` |
| `[u_1|w_1] tensor [x_1]` | `1/4` |
| `[t|w_2] tensor [x_1]` | `1/4` |
| `[v_1|w_1] tensor [x_1]` | `3/4` |
| `[u_1|w_2] tensor []` | `3/4` |
| `[u_2|w_1] tensor []` | `3/4` |
| `[v_1|w_2] tensor []` | `3/4` |
| `[v_2|w_1] tensor []` | `1/4` |
| `[u_1|t|t|w_1] tensor []` | `1/2` |
| `[t|u_1|t|w_1] tensor []` | `3/4` |
| `[t|t|t|w_2] tensor []` | `1/4` |
| `[t|t|v_1|w_1] tensor []` | `1/2` |
| `[t|v_1|t|w_1] tensor []` | `3/4` |
| `[v_1|t|t|w_1] tensor []` | `1/2` |

Each row has exactly one `w_k`, so this displays all 34 nonzero values.

| Degree-seven basis | e_f, literal rational lift |
| --- | ---: |
| `[w_1] tensor [x_3]` | `19/8` |
| `[w_1] tensor [x_1|x_1]` | `-13/16` |
| `[t|w_1] tensor [x_2]` | `-11/4` |
| `[w_2] tensor [x_1]` | `-3/8` |
| `[w_3] tensor []` | `-7/4` |
| `[u_1|t|w_1] tensor []` | `-1/2` |
| `[t|u_1|w_1] tensor []` | `-1/4` |
| `[t|t|w_2] tensor []` | `1/8` |

This displays all 16 nonzero primitive coefficients. The largest
denominator is 16. The matrices and all zero entries remain available in
the JSON, and the runtime checks their source hashes and basis ordering.

## 8. Evaluated normalization selectors

The tertiary low selectors form the vector
\(\boldsymbol{\zeta}=(\zeta_1,\zeta_2,\zeta_3)=(0,1,0)\).
Its entries select the R1 rank normalization and the two R2 suspension
terms, respectively. This vector is distinct from the lower-operation
calibration vectors epsilon `(1,0,0)` and eta `(1,0,1)`, and from the
cochain helpers \(\zeta_{i,n}\).
The packaged selector data are
[low_calibration.json](../python/low_calibration.json) and
[high_calibration.json](../python/high_calibration.json).

| Selector | Actual value | Fixed period or convention |
| --- | ---: | --- |
| R1 rank selector \(\zeta_1\) | `0` | rank period `2` |
| R2 suspension selector \(\zeta_2\) | `1` | order-four period `1/4`, orientation `j=1` |
| R2 suspension selector \(\zeta_3\) | `0` | `kappa_6 V2(Ss)=3/2`, `Uraw(Tor)=-1/2` |
| Odd base selector `eta` | `0` | base period `0`, source `(B,C)=(7/8,0)` |
| New final-T rank ambiguity coefficient | `0` | Danus reference |
| New final-T twist ambiguity `mu_R` | `0` | Euler value `Xi(D)=V2(D)=13/4`, with full J-plus-Psi indeterminacy |
| R3 `xi` | `3/4` | cubic period |
| R3 `(c4,cN,cO,cM,epsilon_c)` | `(1,0,1,1,1)` | current `R2sharp` family |
| Final-T prime-three coefficient | `2` | only input degree three in this range |

For clarity, the R3 suspension row is evaluated on the degree-six small
cycles

\[
F=u_{0,3},\quad Y=u_{0,1}[x_1|x_1],\quad
W=u_{2,2},\quad V=u_{0,2}[x_1],\quad
Z+T=u_{1,1}[x_2]+u_{2,1}[x_1].
\]

For `C0=kappa_7 V3(SA)-V2fin`, where `kappa_7=+I`, their periods are
`(3/4,1/4,0,1/2,1/2)`. The orientations are `aF=jL=1`, and the selector
rules are

\[
\xi=a_FC_0(F),\quad c_4=j_L\,4C_0(Y)\pmod4,\quad
(c_N,c_O,c_M)=2(C_0(W),C_0(V),C_0(Z+T))\pmod2,\quad
\epsilon_c=2(\xi-\tfrac14)\pmod2.
\]

The preceding comparison uses `kappa_6=-I`. None of these fixed source
normalizations is an additional local choice in a GAP calculation. The
old `TReference` coefficient one belongs to a different reference and is
not applied to the new final natural T.

## 9. Provenance and scope

This is a summary of the actual packaged formulas and finite data, checked
against their source on 2026-09-23. It does not replace the longer source
proofs in
tertiary_R_Danus (source-workspace provenance: `notes/extra/tertiary_R_Danus/README.md`; not bundled),
tertiary_T_degree3 (source-workspace provenance: `notes/extra/tertiary_T_degree3/README.md`; not bundled), or the
saved executable
verification record (source-workspace provenance: `notes/extra/tertiary_T_degree3/gap_verification_20260923/README.md`; not bundled).
Some copied module header comments still describe verification prototypes;
their functions are now the kernels imported by the production worker.
The production domain of the tertiary assembler remains inputs `0,1,2,3`.
No all-degree R extension or abutment-extension computation is asserted.
