# Growth and periodicity of the universal source values

The extension worker keeps two universal source functions in a persistent
store (`UNIVERSAL` in
[extension_acceleration.py](../python/extension_acceleration.py)): the
degree-five pair source `production_gamma4.source_value` of the D-layer
stacking correction, and the degree-one primitive `low_phases.V1_pair` of
[universal_helpers.md](universal_helpers.md) Section 5. Both are exact
rational functions of finitely many integer labels, and a computation can
request them at labels of any size. This note answers whether they have a
finite closed form. The pair source is **not** periodic in its labels
modulo any integer, because it contains the non-additive part of `V_1`;
`V_1` grows linearly in the labels, exactly affinely on residue classes
modulo four; and both functions are determined by finite data through an
explicit closed form of the degree-one contraction. No formula,
calibration table or stored value is changed here. Notation follows
[conventions.md](conventions.md) and the helper reference; the statements
marked as proved follow from the implemented definitions, the others are
finite exact checks recorded in Section 8.

## 1. The two functions

**Pair source.** A key is `(Borel(sigma,a,a'),W)`: `sigma` in `{0,1}^5`
are the sign edges and `a,a'` in `Z^5` the two integral edge labels of the
shared-sign Borel 5-simplex (`r1_pair_chain_signed.Borel`), and `W` is the
binary background diagonal. On the model simplex `(0,...,5)`
(`v1_pair_shared.from_diag`)

\[
s_{ij}=\sum_{i\le t<j}\sigma_t,\qquad
A_{ij}=(-1)^{\sigma_0+\cdots+\sigma_{i-1}}\sum_{i\le t<j}a_t,\qquad
A'_{ij}=(-1)^{\sigma_0+\cdots+\sigma_{i-1}}\sum_{i\le t<j}a'_t,
\]

and `omega` is reconstructed from `W` as in Section 6 of the helper
reference. The value (`production_gamma4.source`,
`upper_pair_source.source`) is the degree-five rational cochain

\[
S(A,A')=\Phi(A,0,0)+\Phi(A',0,0)-\Phi(A+A',\alpha,\beta)
-\Theta^{\rm pair}(P,k,P',k')
+R_A(A)+R_A(A')-R_A(A+A')+\tfrac12\widetilde{J_0}
\]

on `(0,...,5)`. Here `Phi` is the natural potential `closed_a_upper.phi`
with `d_s Phi = -Theta(curvature)` modulo integers,
`alpha=hD(rho_2A,rho_2A',s)`, `beta` is `closed_a_upper.beta_sharp`,
`P,k` and `P',k'` are the binary primary and secondary data of `A` and
`A'`, `Theta^pair` is the pair phase of `theta_pair_phase.py`, `R_A` is
the production rephasing `production_gamma4_comparison.A_phase`, and
`J_0` is the binary residue `production_upper_binary_comparison.affine_J0_source`.
By definition the stored value is `0` when `a=0` or `a'=0`.

**Degree-one primitive.** A key is `((g_0,...,g_5),W)` with
`g_i=(k_i,epsilon_i)` in `Z x| C_2`, `g_0=(0,0)`, and the value is

\[
V_1=\phi_1K+U_{\rm wedge}P,\qquad \phi_1=\Theta_3(q_1,k_1),
\]

exactly as in Section 5 of the helper reference. The pair source calls it
three times through `closed_a_upper.source_primitive`
(`low_phases.V1`, `to_pair1`): with

\[
g_i=(A_{0i},s_{0i}),\qquad g_i=(A'_{0i},s_{0i}),\qquad g_i=((A+A')_{0i},s_{0i}),
\qquad A_{0i}=a_0+\cdots+a_{i-1},
\]

and the symmetrised upper triangle of `W` as background. Write
`c(g)=2k-epsilon` for the Cayley coordinate of `g`, and `V_1(m)` for the
value after replacing `k_j` by `k_j+m` at a fixed vertex `j>=1`.

## 2. Structure of the pair source

**Theorem 1 (decomposition).** As exact rational numbers,

\[
\boxed{\;S(\sigma,a,a',W)=S_{\rm rest}(\sigma,a,a',W)
+V_1(A+A')-V_1(A)-V_1(A'),\;}
\]

where `S_rest` is the same formula with its three `V_1` calls replaced by
zero. Moreover

1. `S_rest` depends on `a,a'` only through their residues modulo four,
   except through the binary cochain `beta`;
2. `beta` is a mod-two sum, over the unit-edge expansion of the labels, of
   binary source values that each depend on the labels modulo four.

*Proof.* `production_gamma4.source` is the displayed sum, and
`closed_a_upper.phi(A,B,C) = cross - source_primitive(A) - prism Theta(b_I,f_I)`.
In degree one, `source_primitive(A) = V - prism Theta(shifted P, shifted k)`
with `V(sigma)=V_1(to_pair1(sigma))`; on the model simplex `to_pair1` gives
the three keys displayed in Section 1. The three `Phi` terms therefore
contribute exactly `-V_1(A)-V_1(A')+V_1(A+A')`, and every other term of
`Phi` and of `S` is built from `A` in one of the following ways:
through `phase_eval.source_splitting`, which uses only `a=rho_2A` and the
carry `t=rho_2((A-a)/2)`, that is `rho_4A`; through `rho_2A` directly
(`closed_a_upper.primary`, the `hD` in `alpha`, `affine_J0_source`); and
then through binary operations only (interval-cut words, `zeta_1`,
`zeta_2`, the chi ANF, `Q^i`, `E`, `Q_D`, `hD`, `hE`, and the binary `KB`
models of the pair phase in `a0_high_gamma.HigherA0Stacking`), or through
integer operations on integer cochains that are themselves built from
binary data (`B=d a/2`, the integral cup terms of `Theta` divided by
four, the carry `L_A`). The interval lift and the right prism only
re-index vertices. The single term linear in the integral `A`, the
Pontryagin term `-(epsilon/4) P(omega) smile_s A` of `A_phase`, is linear
in `A` and cancels in `R_A(A)+R_A(A')-R_A(A+A')`; besides,
`epsilon=0` in the packaged selectors. This proves 1 except for `beta`.

For 2, `beta_sharp(A,0,A',0)` reduces to `v1_pair_shared.legal_beta`,
whose only non-binary ingredient is `v1_pair_shared.primitive`: a mod-two
sum of the binary pair source over the chain `r1_pair_chain.Htot`, whose
vertical homotopy `h0` is the unit-edge contraction `r2_pair_chain.unit_h`.
Its raw step `chain_models._unit_raw` replaces an edge of label `L` by the
`|L|` unit edges between its endpoints, so the chain has a number of terms
growing with the labels, each term carrying labels that are partial sums
of the original ones; the source value of a term depends on its labels
modulo four by the first part. `beta` enters `S` only through
`cross = (1/2) rho_2(E(beta)+hE(d beta,f))`. QED.

The mod-two sum in 2 is a sum of a period-four sequence over a number of
terms growing with the label, so its period in a label is a power of two
that need not be four. In the finite checks of Section 8 the binary table
of `beta` on the fifteen faces had period 1, 2 or 4 in 54 of 60
(key,label) series and period exactly 8 in the other 6, and `S_rest` had
period 8 (in 16 of 400 series) or 4 (384 of 400).

The short-circuit `S=0` for `a=0` or `a'=0` agrees with the formula, whose
value is also `0` at every zero-`a` key tested; but the first step away
from the zero vector can already start a drift (Section 5).

## 3. The degree-one contraction in closed form

Fix `g=(g_0,...,g_i)` with coordinates `c_t=c(g_t)`. Let `seg(x,y)` be
the set of unit edges of the Cayley line between the coordinates `x` and
`y`, and

\[
E_{l,i}=\bigcap_{t=l}^{i-1}\operatorname{seg}(c_t,c_{t+1}),
\]

the edges crossed by **every** step of the path `c_l -> ... -> c_i`. For a
unit edge `e={v,v'}` write `alt_d(e;v)=(v,v',v,v',...)` for the alternating
`d`-simplex on `e` starting at `v`, and `near_l(e)` for the endpoint of `e`
nearest to `c_l` (the two endpoints are never equidistant). Recall from the
helper reference that `K_root` sends a vertex to the oriented unit-edge
path from `root`, and an alternating simplex on one unit edge to the same
simplex with its endpoint nearest to `root` prepended, degenerate results
being dropped (`low_phases._K`, `_cone`, `_gbasis`).

**Lemma 1.** For `i>=1`, `R[g_0..g_i]=sum_{e in E_{0,i}} alt_i(e;near_0(e))`.

*Proof.* For `i=1`, `R partial[g_0,g_1]=[g_1]-[g_0]` and `K_{g_0}` returns
the path from `c_0` to `c_1`: the edges of `seg(c_0,c_1)=E_{0,1}`, each
oriented away from `c_0`, which is `alt_1(e;near_0(e))`. For the induction
step, every face `R[g_0..hat g_t..g_i]` with `t>=1` consists, by
hypothesis, of alternating simplices starting at `near_0(e)`; prepending
`near_0(e)` makes them degenerate, so `K_{g_0}` kills them. The face `t=0`,
with sign `+1`, is `sum_{e in E_{1,i}} alt_{i-1}(e;near_1(e))`, and
prepending `near_0(e)` is non-degenerate exactly when
`near_0(e) != near_1(e)`, i.e. when `c_0` and `c_1` lie on opposite sides
of `e`, i.e. when `e` is in `seg(c_0,c_1)`. Hence
`R[g_0..g_i]=sum_{e in E_{1,i} cap seg(c_0,c_1)} alt_i(e;near_0(e))`, and
`E_{1,i} cap seg(c_0,c_1)=E_{0,i}`. In particular every simplex reaching
`K_root` lies on one unit edge, so the guard in `_K` never fires. QED.

**Lemma 2.** `h[g_0..g_i]=sum_{l=0}^{i-1}(-1)^{l+1}(g_0..g_l)*R[g_l..g_i]`,
where `(g_0..g_l)*sigma` prepends the vertices `g_0,...,g_l` to `sigma`
(degenerate results dropped).

*Proof.* In `h[g]=C_{g_0}([g]-R[g]-sum_t(-1)^t h[face_t])` the cone over
`[g]` and over every `h[face_t]` with `t>=1` is degenerate, since those
chains already start with `g_0`. So `h[g]=-C_{g_0}R[g]-C_{g_0}h[g_1..g_i]`,
and unrolling the recursion gives the formula. QED.

**Proposition 3 (closed form of `V_1`).** With the Alexander-Whitney
fronts `(g_0..g_i)` and backs `W|[i..5]` of `group_product_AW`
(degenerate fronts or backs omitted),

\[
\begin{aligned}
V_1(g,W)={}&\phi_1\bigl(\widehat H_{\rm prod}(g,W)\bigr)
+\sum_{i=1}^{5}\sum_{l=0}^{i-1}(-1)^{l+1}\sum_{e\in E_{l,i}}
 \phi_1\Bigl(\operatorname{shuffle}\bigl((g_0..g_l)*\operatorname{alt}_{i-l}(e;\operatorname{near}_l(e)),\,W|_{[i..5]}\bigr)\Bigr)\\
&+\sum_{i=1}^{5}\sum_{e\in E_{0,i}}
 U_{\rm wedge}\Bigl(\operatorname{shuffle}\bigl(\operatorname{alt}_i(e;\operatorname{near}_0(e)),\,W|_{[i..5]}\bigr)\Bigr).
\end{aligned}
\]

This is `V_1=phi_1 K+U_wedge P` with Lemmas 1 and 2 substituted into
`K=Hhat_prod+shuffle(h tensor 1)AW` and `P=shuffle(R tensor 1)AW`
(`low_phases.V1_pair`). Two properties of the summands follow from the
definitions:

* `phi_1` sees the group vertices only modulo `Z/4 x Z/2`: it is evaluated
  on the registered source with `A_{ij}=(-1)^{epsilon_i}(k_j-k_i)` and
  `s_{ij}=epsilon_i+epsilon_j` (`low_phases.from_pair1`), and
  `Theta_3(q_1,k_1)` depends on `A` only through `rho_4A`, by the argument
  of Theorem 1.
* `U_wedge` is translation invariant: on an alternating unit-edge simplex
  all edge labels are `(0,1)` or `(1,1)`, so the registered source is the
  same for every `(1,1)` edge, and the only position dependence is the sign
  `(-1)^{epsilon_0}` of the first vertex (`low_phases._wedge_primitive`,
  `_evaluate_pair`), which is determined by the parity of the edge and by
  the side on which the root lies.

## 4. Growth of `V_1`

**Theorem 2 (exact affine law on residue classes).** Fix `(g,W)` and a
vertex `j>=1`, let `M=max_{l!=j}c_l`, `M'=min_{l!=j}c_l`, and

\[
m_+=\min\{m:\ c_j+2m>M\},\qquad m_-=\max\{m:\ c_j+2m<M'\}.
\]

Then for all `m>=m_+`, `V_1(m+4)-V_1(m)=P_+(m mod 4)` depends on `m` only
through `m mod 4`; likewise `V_1(m+4)-V_1(m)=P_-(m mod 4)` for all
`m+4<=m_-`. Consequently

\[
\boxed{\;V_1(m)=a_\pm(m\bmod4)\,m+b_\pm(m\bmod4)\quad
\text{for }m\ge m_+\ \text{resp. }m\le m_-,\qquad a_\pm=P_\pm/4.\;}
\]

The class-dependence of `a_+` (and of `a_-`) can only come from the family
`(l,i)=(j,j+1)` of Proposition 3, so for the last vertex `j=5` the four
classes have one common slope. The slopes `a_+` and `a_-` are unrelated in
general.

*Proof.* For `m>=m_+` the moving coordinate `c_j(m)=c_j+2m` lies to the
right of all other vertices. Steps of the path not involving `j` have
fixed segments; the steps `j-1 -> j` and `j -> j+1` have
`seg=[c_{j-1},M) cup B(m)` and `[c_{j+1},M) cup B(m)`, where
`B(m)=[M,c_j(m))` is the *bridge*. Hence every `E_{l,i}` is a fixed set of
edges below `M`, together with the bridge exactly when
`{l,...,i-1}` is contained in `{j-1,j}`, that is for the three families
`(j-1,j)`, `(j,j+1)`, `(j-1,j+1)` (the last two only if `j<5`); and on
the bridge `near_l(e)` is the left endpoint unless `l=j`, when it is the
right one. A fixed term depends on `m` only through the values
`A(g_t,g_j)=±(k_j+m-k_t)`, hence through `m mod 4`. A bridge term at the
edge `e` depends on the parity of `e` (edge type and near endpoint), on
`(k_e-k_t) mod 4` for the prefix vertices `t<=l`, and, for the family
`(j,j+1)` only, on `(k_j+m-k_t) mod 4`; as a function of the coordinate of
`e` it has period 8. Passing from `m` to `m+4` appends eight edges, one
full period, to every bridge and leaves `m mod 4` unchanged, so the
difference is the sum over one period, a constant for fixed `m mod 4`.
The `Hhat_prod` term has no bridge. The left side is the same argument
with the bridge `[c_j(m),M')`, on which the near endpoints and the
orientations of the unit edges are the other ones (reverse bar edges are
different generators), which is why `a_-` is unrelated to `a_+`. QED.

Between the two thresholds, `m_-<m<m_+`, the ordering of the Cayley
vertices changes and no simple law holds; this transient region has at
most `(max_l c_l-min_l c_l)/2+1` values of `m`, which for the stored keys
(all `|k|<=2`) means `|m|<=3`. Beyond the thresholds the law holds at
every checked `m` (often it already holds one to three units before the
threshold; never later).

**Example.** For the stored key `g=((0,0),(1,1),(0,0),(1,1),(0,0),(1,1))`
with `W_ij=1` for `i!=j` and the vertex `j=1` (`V_1(0)=-65/16`),

\[
V_1(m)=\frac{15}{16}m+\Bigl(-\frac{65}{16},-\frac{49}{16},-\frac{57}{16},-\frac{41}{16}\Bigr)_{m\bmod4}
\quad(m\ge0),
\]

so `V_1(1,2,3,4,8,32)=-17/8,-27/16,1/4,-5/16,55/16,415/16`, while for
`m<=-1`

\[
V_1(m)=\frac{39}{16}m+\Bigl(-\frac{49}{16},\cdot,-\frac{41}{16},\cdot\Bigr)_{m\bmod4}
\ (m\text{ even}),\qquad
V_1(m)=\frac{31}{16}m+\Bigl(\cdot,-\frac{33}{16},\cdot,-\frac{25}{16}\Bigr)_{m\bmod4}
\ (m\text{ odd}),
\]

e.g. `V_1(-1,-2,-3,-4,-5,-6)=-7/2,-119/16,-63/8,-205/16,-45/4,-275/16`:
the mean slope on the left is `35/16`, not `15/16`. For the vertex `j=4`
of the same key the right slopes are `(0,5/8,1/4,3/8)` on the classes
`m=0,1,2,3 mod 4`: `V_1(0,4,8)=-65/16`, `V_1(1,5,9)=-9/16,31/16,71/16`,
`V_1(2,6,10)=-53/16,-37/16,-21/16`, `V_1(3,7)=-9/16,15/16`.

**Slopes of the checked keys.** Eight keys (six stored nonzero keys and two
synthetic ones) times five vertices; `m` in `[-24,24]` and
`±32,±40,±48,±64`. A four-tuple lists the classes `m=0,1,2,3 mod 4`; the
22 rows below are those with a nonzero slope on at least one side,
all other `(key,j)` have `a_+=a_-=0`.

| key | `g` (background) | `j` | `m_+` | `a_+` | `m_-` | `a_-` |
|---|---|---|---|---|---|---|
| 0 | `((0,0),(1,1),(0,0),(1,1),(0,0),(1,1))`, `W_ij=1` for `i!=j` | 1 | 1 | `15/16` | -1 | `(39/16, 31/16, 39/16, 31/16)` |
| 0 | `((0,0),(1,1),(0,0),(1,1),(0,0),(1,1))`, `W_ij=1` for `i!=j` | 2 | 1 | `(3/4, 11/8, 1, 13/8)` | -1 | `(3/4, 3/8, 7/8, 1/4)` |
| 0 | `((0,0),(1,1),(0,0),(1,1),(0,0),(1,1))`, `W_ij=1` for `i!=j` | 3 | 1 | `-5/4` | -1 | `3` |
| 0 | `((0,0),(1,1),(0,0),(1,1),(0,0),(1,1))`, `W_ij=1` for `i!=j` | 4 | 1 | `(0, 5/8, 1/4, 3/8)` | -1 | `(1/4, -1/8, 1/4, -1/8)` |
| 0 | `((0,0),(1,1),(0,0),(1,1),(0,0),(1,1))`, `W_ij=1` for `i!=j` | 5 | 1 | `-5/8` | -1 | `11/8` |
| 1 | `((0,0),(2,1),(0,0),(2,1),(0,0),(2,1))`, `W_ij=1` for `i!=j` | 1 | 1 | `15/16` | -2 | `(31/16, 39/16, 31/16, 39/16)` |
| 1 | `((0,0),(2,1),(0,0),(2,1),(0,0),(2,1))`, `W_ij=1` for `i!=j` | 2 | 2 | `(7/4, 3/4, 11/8, 1)` | -1 | `(3/4, 5/8, 5/4, 1)` |
| 1 | `((0,0),(2,1),(0,0),(2,1),(0,0),(2,1))`, `W_ij=1` for `i!=j` | 3 | 1 | `-3/4` | -2 | `7/4` |
| 1 | `((0,0),(2,1),(0,0),(2,1),(0,0),(2,1))`, `W_ij=1` for `i!=j` | 4 | 2 | `(1/4, 0, 0, 1/4)` | -1 | `0` |
| 1 | `((0,0),(2,1),(0,0),(2,1),(0,0),(2,1))`, `W_ij=1` for `i!=j` | 5 | 1 | `-1/4` | -2 | `3/4` |
| 2 | `((0,0),(2,0),(1,0),(0,0),(0,0),(0,0))`, `W` rows (00000),(00000),(00000),(00001),(00010) | 2 | 2 | `(-1/8, 1/4, 0, -1/8)` | -2 | `(-3/8, -3/4, -3/4, -3/8)` |
| 2 | `((0,0),(2,0),(1,0),(0,0),(0,0),(0,0))`, `W` rows (00000),(00000),(00000),(00001),(00010) | 3 | 3 | `-3/8` | -1 | `3/8` |
| 2 | `((0,0),(2,0),(1,0),(0,0),(0,0),(0,0))`, `W` rows (00000),(00000),(00000),(00001),(00010) | 4 | 3 | `(0, -1/4, 0, 0)` | -1 | `(0, 1/4, 1/4, 1/4)` |
| 3 | `((0,0),(2,0),(2,1),(0,1),(2,1),(0,1))`, `W=0` | 4 | 1 | `1/4` | -3 | `0` |
| 4 | `((0,0),(2,0),(0,0),(0,1),(2,1),(0,1))`, `W=0` | 4 | 1 | `(-1/4, -1/4, 0, 0)` | -3 | `(0, 1/4, 0, 1/4)` |
| 4 | `((0,0),(2,0),(0,0),(0,1),(2,1),(0,1))`, `W=0` | 5 | 3 | `0` | -1 | `1/4` |
| 5 | `((0,0),(2,1),(0,1),(0,0),(2,0),(0,0))`, `W=0` | 4 | 0 | `(-1/4, 1/8, 1/8, -3/8)` | -3 | `(1/8, -1/8, -1/8, -1/8)` |
| 5 | `((0,0),(2,1),(0,1),(0,0),(2,0),(0,0))`, `W=0` | 5 | 3 | `-1/4` | -1 | `3/8` |
| 6 | `((0,0),(2,1),(-1,0),(1,1),(3,0),(-2,1))`, `W` rows (01010),(10001),(00000),(10000),(01000) | 4 | -1 | `(1/8, 1/4, 0, -1/8)` | -6 | `(1/8, 0, 1/8, 0)` |
| 6 | `((0,0),(2,1),(-1,0),(1,1),(3,0),(-2,1))`, `W` rows (01010),(10001),(00000),(10000),(01000) | 5 | 6 | `-3/8` | 1 | `7/8` |
| 7 | `((0,0),(1,0),(1,1),(-2,1),(2,0),(0,1))`, `W` rows (00110),(00001),(10000),(10000),(01000) | 4 | 0 | `(0, -1/8, 0, 0)` | -5 | `(1/8, 0, 0, 0)` |
| 7 | `((0,0),(1,0),(1,1),(-2,1),(2,0),(0,1))`, `W` rows (00110),(00001),(10000),(10000),(01000) | 5 | 3 | `-1/8` | -3 | `0` |

All slopes have denominators dividing 16. Of the 80 (key,vertex,side)
cases, 59 have a class-independent slope; `a_+=a_-` in 19 of 40 keys and
vertices (mostly both zero) and `a_+=-a_-` in 21.

## 5. Consequences for the pair source

**Theorem 3.** `S` is not periodic in its labels modulo any integer. More
precisely, along a shift `d` of one label `a_i` (which adds `d` to
`A_{0l}` for `l>i`, `d` to nothing in `A'`, and `d` to `(A+A')_{0l}`), the
three `V_1` terms of Theorem 1 are eventually affine on residue classes
with slopes `a_A, a_{A'}, a_{A+A'}` that depend on the fixed vertices and
on `W`; the combination has slope `a_{A+A'}-a_A-a_{A'}`, nonzero in
general, and even when it vanishes the two-sided patterns `b_±` of the
three terms differ, so `S(d)=S(d+16)` fails across the transient region.

Exact instances (all values are literal rationals):

1. *Pure drift.* `sigma=(0,0,0,0,1)`, `a'=(2,-2,2,-2,0)`, `W=0`,
   `a=(0,0,0,0,L)`. Here `S_rest=0` and `V_1(A)=V_1(A')=0`, so `S` is the
   single family `phi(4,5)=-(g_0..g_4)*R[g_4,g_5]` of `V_1(A+A')` summed
   over the `|2L-1|` unit edges between `c_4=0` and `c_5=2L-1`:

   \[
   S(L)=-\tfrac12\lfloor L/2\rfloor\ (L\ge0),\qquad
   S(L)=-\tfrac12\lceil |L|/2\rceil\ (L\le0),
   \]

   for every `L` in `[-24,40]`; `S(1)=0`, `S(17)=-4`, `S(33)=-8`.
2. *Stored key, drift in both directions.* `sigma=0`, `a=(0,1,0,0,0)`,
   `a'=(1,0,0,0,0)`, `W` rows `(0,0,0,0,0),(0,0,0,0,0),(0,0,0,1,1),(0,0,1,0,1),(0,0,1,1,0)`,
   label `a_2+=d`: `S(0)=0`, `S(8)=-4`, `S(16)=-8`, `S(32)=-16`, `S(40)=-20`,
   `S(-8)=-2`, `S(-16)=-4`; `S_rest=0` throughout.
3. *Transition without drift.* `sigma=(1,0,1,0,0)`, `a=a'=(1,0,-1,0,0)`,
   `W` rows `(0,1,0,1,1),(1,0,0,1,1),(0,0,0,0,0),(1,1,0,0,1),(1,1,0,1,0)`,
   label `a_1+=d`: `S(0)=-5/4`, `S(8)=S(16)=S(32)=3/4`, `S(-8)=S(-16)=-5/4`.
   The `V_1` combination is `(0,-3/2,0,-1)_{d mod 4}` for `d<=0` and
   `(2,7/4,5/2,9/4)_{d mod 4}` for `d>=2` (zero slope on both sides), while
   `S_rest` has period 8 throughout:
   `S_rest(0..7)=-5/4,-3/2,-3/4,-3/2,-3/4,-2,-5/4,-2`.

In the finite checks (Section 8) the value was unchanged under `+16` in
3180 of 3370 (stored key, label) pairs; the 190 failures lie in 84 keys,
are the same for `+8`, `+16` and `+32`, and in every one of them
`S(16)-S(0)` equals the change of the `V_1` combination, so `S_rest` was
16-periodic in all 3370 pairs. On dense random keys the failure rate is
much higher (46 of 120 single-label shifts, 19 of 36 joint shifts).

## 6. Mathematical intuition

*Why almost everything is four-periodic.* The secondary and tertiary
formulas are cochain representatives of operations that see an integral
class only through its reduction modulo four: the binary lift `rho_2A`,
the integral carry `(A-rho_2A)/2` reduced modulo two, the Bockstein
`d rho_2A/2`, and the quarter-valued Pontryagin-square terms
`(omega smile B+B smile_{m-1}B)/4` of `Theta`. Everything downstream is
a function of these binary and integer cochains, so the labels enter
through `Z/4`. The universal models are therefore, for all these terms,
pulled back from a finite quotient: the Cayley line of `Z x| C_2` modulo
translation by four in `k` is a cycle of eight unit edges (two edge types
`(0,1)` and `(1,1)`, alternating).

*Why `V_1` grows linearly.* `V_1` is a potential: `d_sV_1=phi_1` modulo
integers. Contracting the dihedral group to its Cayley line turns the
potential into a path sum over unit edges from the base vertex, and the
integrand, by the previous paragraph, is periodic along the line with
period eight. A discrete integral of a periodic sequence is an affine
function plus a periodic one; on each residue class the discrete integral
is exactly affine, with slope the average of the integrand over a
period. The finer structure of Theorem 2 reflects two features of the
contraction: the cone `C_{g_0}` attaches the moving unit edges to the
fixed vertices, so the integrand contains the cup products of the fixed
labels with the moving ones (through `A(g_t,v) mod 4`), and when the
moving vertex itself belongs to the cone prefix, the phase of the
integrand shifts with `m mod 4`, which is the only source of a
class-dependent slope; and reverse bar edges are different generators, so
the left and right bridges are different chains with different averages.

*Why the pair source is not periodic.* `phi_1=Theta_3(q_1,k_1)` contains
terms quadratic in the label data (products such as `u smile w` and the
Pontryagin-square term), so its potential behaves like the potential of a
quadratic form, and `V_1(A+A')-V_1(A)-V_1(A')` is its polarization. For a
genuinely quadratic function the polarization is bilinear, hence linear
in each label. The pair source is exactly this polarization plus a
four-periodic part: its linear growth is the bilinear part, whose
coefficient `a_{A+A'}-a_A-a_{A'}` is a period average that depends on the
fixed vertices, and its periodic part is the torsion part seen modulo four.
Crossing another Cayley vertex changes which unit edges are covered by
every step of the folded path, i.e. which alternating simplices occur in
`R` and `h`; this wall-crossing produces the one-time transitions of
Theorem 3(3), where two different periodic branches meet.

*Why a finite description still exists.* The whole infinite family of
values is generated by finitely many `phi_1` and `U_wedge` evaluations on
reduced simplices (vertices in `Z/4 x Z/2` up to translation, edges by
parity), combined by the interval arithmetic of the bridges. The growth
is in the combinatorics of the path, not in the transcendence of the
values.

## 7. Finite data and the universal store

* Reducing the labels of a pair-source key modulo 16 (or any modulus)
  before a lookup is **not** exact: Theorem 3 gives stored keys whose
  value changes under `+16` and is unbounded in the label. A store keyed
  by the raw `Borel` key can therefore not be finite. What is exact is
  Theorem 1: `S_rest` is a function of `(sigma, a mod 4, a' mod 4, W)` up
  to the `beta` term, and of `(sigma, a mod 8, a' mod 8, W)` in every
  checked case (proved only modulo four for the `beta`-free part), plus
  three closed-form `V_1` values.
* `V_1` has the closed form "finite table plus explicit linear term":
  every family `E_{l,i}` of Proposition 3 is an interval of unit edges,
  and its summands depend on the edge only through its coordinate modulo
  8, so

  \[
  \sum_{e\in E_{l,i}}F(e)=q\cdot\sum_{\text{one period}}F+\sum_{\text{partial}}F,
  \qquad q=\lfloor |E_{l,i}|/8\rfloor,
  \]

  with at most eight distinct values of `F` per family. The finite data
  that determine `V_1` completely are the values of `phi_1` on the product
  simplices of Proposition 3 with group vertices reduced modulo
  `Z/4 x Z/2` and the values of `U_wedge` on alternating `(1,1)`-edge
  simplices shuffled with background backs; the raw `(g,W)`-keyed store is
  unbounded and its evaluation time grows linearly with the labels.
* The pair source is determined by the finite `S_rest` data together with
  the closed form of `V_1`.

## 8. Verification record (2026-09-27)

Base revision `f9eb8b6` with the uncommitted
[shared-kernel changes](verification/shared-kernel-20260927.md), which
that record shows return identical universal values. The checks ran in
fresh processes with `FERMIONAHSS_CACHE_DIR` empty (no persistent store
written), with the exact evaluation policy of the extension worker
installed, on the 1129 stored `source_value` keys (337 with both `a` and
`a'` nonzero) and the 714 stored `V1_pair` keys of a degree-four
`koFull(CyclicGroup(2),[1],[1],4)` run, plus synthetic keys. About two
hours of wall time on a loaded 32-thread machine; the scripts were external
scratch tools and are not bundled. Observed:

| Check | Result |
| --- | --- |
| The user's test: first 40 nontrivial stored keys x 10 labels, shifts `+4,+8,+16` | Unchanged in 374, 388, 388 of 400; the 12 exceptions to `+8` are the 12 exceptions to `+16`, scanned for every `d` in `[-24,40]` (6 drifts, 6 transitions) |
| Theorem 1, `S = S[V_1:=0] + (V_1(A+A')-V_1(A)-V_1(A'))` | 400 of 400 evaluations (40 keys x 10 labels x shifts `0,4,8,-4,16`); the three `V_1` calls recorded inside `S` equal the constructed keys |
| Period of `S[V_1:=0, beta:=0]` under `±4,+8` | 400 of 400 unchanged |
| Period of `S[V_1:=0]` | unchanged under `+8` and `+16` in 400 of 400, under `±4` in 384 of 400 |
| Binary table of `beta` on the 15 faces, 6 stored keys x 10 labels, `d` in `[-20,20]` | period 1, 2 or 4 in 54 series, exactly 8 in 6; those 6 re-checked for `d` in `[-40,40]`: `beta(d)=beta(d+8)` throughout |
| All 337 nontrivial stored keys x 10 labels, shifts `+8,+16,+32,-16` | unchanged in 3180, 3180, 3180, 3215 of 3370; identical set of 190 failing pairs in 84 keys for the three positive shifts; failures by label index `a: (0,7,16,11,64)`, `a': (1,8,13,12,58)`; in all 190, `S(16)-S(0)` equals the change of the `V_1` combination |
| 12 dense random keys (labels in `[-2,2]`), shifts `±16` of each label and three joint shifts | unchanged in 74 and 74 of 120, and 17 of 36 |
| Formula value at 12 stored keys with `a=0`, without the short-circuit | `0` in all 12; `S(a=16e_4)=-4` for 3 of them |
| Proposition 3 against `V1_pair` (independent evaluator of the closed form) | equal in 317 of 317 keys, 160 of them transients with `-3<=m<=3`; the front-by-front split of the code's own chains agrees in 96 of 96 |
| Theorem 2: 8 keys x 5 vertices, `m` in `[-24,24]` and `±32,±40,±48,±64` (2280 values) | affine law on every residue class at every `m` beyond the threshold in 80 of 80 sides; slope denominators 1, 4, 8, 16 |
| Family attribution of the drifts (values of the closed form per family) | drift of instance 1 entirely in `phi(4,5)`; class-dependent increments only in the family `(j,j+1)` |
| Re-evaluation on the final working tree (fresh process, 90 s): instance 1 at `L=1,17,33`; instance 3 at `d=0,8,16`; instance 2 at `d=0,8,16,32`; the example key at `m=0,1,2,3,4,8,-1,-2,-5,-6`; the closed form at two shifted keys; Theorem 1 at one shifted key | `0,-4,-8`; `-5/4,3/4,3/4`; `0,-4,-8,-16`; `-65/16,-17/8,-27/16,1/4,-5/16,55/16,-7/2,-119/16,-45/4,-275/16`; equal; equal, all as stated above |

The numerical rows support but do not replace the proofs of Theorems 1
and 2 and Proposition 3; the period of `beta` beyond the checked ranges
is not established. `certified_ko` and the scope statements of
[mathematical-status.md](mathematical-status.md) are unchanged.
