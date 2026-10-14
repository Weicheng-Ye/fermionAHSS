# Universal primitives

Several corrections of the product are *primitives*: binary or rational
cochains $V$ with a prescribed coboundary $\delta V=S$ (modulo $2$, or
modulo integers) for an explicitly given closed *source* $S$. They are
never obtained by solving $\delta V=S$ on $X$ (two solutions differ by a
cocycle, and a change of a rational primitive by a $\mathbb Q/\mathbb Z$
cocycle changes the integral product by the Bockstein of that cocycle).
Instead each is the value of one fixed chain contraction of a *universal
model* of the input data. This page specifies that construction once
(Sections 1–3), the models and sections (Section 4), and every primitive
used by the stacking formulas (Sections 5–9). Function names refer to
`python/stacking_model` and `python/chain_models.py` of the previous implementation.

## 1. The construction

Let $M$ be a simplicial set of universal input data and
$j:\{\text{simplices of }X\}\to M$ a *section*: it sends an ordered
$q$-simplex $\sigma=(v_0,\dots,v_q)$ of $X$ to a $q$-simplex $j(\sigma)$ of
$M$ built from the values of the input cochains on the faces of $\sigma$,
in such a way that the canonical cochains of $M$ pull back along $j$ to
the input cochains (the *reconstruction* identities of Section 4). Let
$(\mathfrak S,\partial)$ be a small chain complex with chain maps
$F:C_*(M)\to\mathfrak S$, $G:\mathfrak S\to C_*(M)$ and a homotopy
$H:C_*(M)\to C_{*+1}(M)$ with

$$
\partial H+H\partial=1-GF .
$$

Let $S$ be the source, an explicit cochain formula in the input cochains,
evaluated on the canonical cochains of $M$; its pullback $S\circ j$ is the
source on $X$. For fixed coefficients $e=(e_w)_{w}$ on the basis of
$\mathfrak S$ in degree $|S|-1$, define

$$
\boxed{\;V(\sigma)=S\bigl(H\,j(\sigma)\bigr)+\sum_we_w\,\langle F\,j(\sigma),w\rangle .\;}
$$

If $S$ is closed on $M$ then $\delta V=S\circ j-(S\circ G-e\circ\partial)\circ F\circ j$,
so $\delta V=S\circ j$ exactly when the **small system** $e\circ\partial=S\circ G$
holds on the basis of $\mathfrak S$ in degree $|S|$. The *periods* of $S$,
its values on the small cycles, must vanish for this to be solvable; the
coefficients are determined once from the values $S(G(w))$ and are part
of the formula (binary primitives: a solution of the $\mathbb F_2$ system;
rational primitives: a specified Smith-form division). When the relevant
part of $\mathfrak S$ is zero there are no coefficients.

## 2. Product contractions

Every model is a product of factor models, contracted factor by factor and
then as a product. For two contractions $(F_X,G_X,H_X)$, $(F_Y,G_Y,H_Y)$:

$$
F_{X\times Y}=(F_X\otimes F_Y)\circ\mathrm{AW},\qquad
G_{X\times Y}=\mathrm{EZ}\circ(G_X\otimes G_Y),\qquad
H_{X\times Y}=H_{\rm prod}+\mathrm{EZ}\circ\bigl(H_X\otimes1+(G_XF_X)\otimes H_Y\bigr)\circ\mathrm{AW}
$$

(Koszul sign $(-1)^{|x|}$ on the second term; `Htot`, `product_homotopy`).
Here, on the product of simplicial sets (diagonal simplices $(x,y)$ of
equal dimension $n$): $\mathrm{AW}(x,y)=\sum_{p}(x|_{0..p},y|_{p..n})$;
$\mathrm{EZ}(x,y)$ for $|x|=p$, $|y|=q$ is the sum over $(p,q)$-shuffles
$(h,v)$ of $\pm(x\circ h,y\circ v)$ with the parity of the inversions; and
$H_{\rm prod}$ is the canonical product homotopy, built once on the
universal pair of nondecreasing vertex maps by the cone recursion

$$
\widehat H_0=0,\qquad \widehat H_n=C_0\Bigl(Q_n-\sum_i(-1)^id_i\widehat H_{n-1}\Bigr),\qquad Q=1-\mathrm{EZ}\circ\mathrm{AW},\qquad
u=Q\widehat HQ,\qquad H_{\rm prod}=u\,\partial\,u,
$$

$C_0$ prepending the initial vertex in both factors (`_raw_h_terms`,
`_h_terms`); it is applied to a pair by pulling both factors along the
resulting vertex maps. Degenerate simplices are discarded throughout.

## 3. Bar-model contractions (`chain_models`)

A **diagonal bar simplex** of dimension $q$ is a list of $q$ rows
$(\sigma_i,e_i)$ with $\sigma_i\in\{0,1\}$ (sign edge) and $e_i\in\mathbb Z^q$
(mode `zsign`) or $\mathbb F_2^q$ (mode `c2`); a **rectangle** has $p$ rows
of $q$-simplices. Faces: horizontal faces multiply adjacent rows
($(\sigma,a)(\tau,b)=(\sigma+\tau,\ a+(-1)^\sigma b)$, with the frame
transport $(-1)^{\sum_{i<v_0}\sigma_i}$ in mode `zsign`), vertical faces are
simplicial faces of each row; $\partial_{\rm rect}=\partial_h+(-1)^p\partial_v$.
The contraction of the diagonal simplices onto *bar words* is

$$
F=F_V\circ F_D,\qquad G=G_D\circ G_V,\qquad H=H_D+G_D\,H_V\,F_D ,
$$

- *diagonal to rectangle:* $F_D(x)=\sum_p x|_{(0..p)\times(p..q)}$, $G_D$ = shuffles of a rectangle into diagonal simplices, $H_D$ = the product homotopy $H_{\rm prod}$ applied to the diagonal (`_F_D`, `_G_D`, `_H_D`);
- *vertical (bar) stage*, the perturbation series in the horizontal direction with
  $f_0$ = iterated Alexander–Whitney of the rows into a tensor word (words containing a unit atom dropped), $g_0$ = iterated shuffle, $h_0=(-1)^p\times$ the iterated product homotopy, and the bar sign $\varepsilon(w)=(-1)^{p+\sum_j(p-j-1)|a_j|}$:
  $$F_V=\sum_{m\ge0}\varepsilon\,f_0\,(-\partial_hh_0)^m,\qquad G_V=\sum_{m\ge0}(-h_0\partial_h)^m\,\varepsilon\,g_0,\qquad H_V=\sum_{m\ge0}h_0\,(-h_0\partial_h)^m$$
  (`_F_V`, `_G_V`, `_H_V`; the series terminate because each step shortens the horizontal length).

In mode `c2` the small complex is the bar complex of words $(n_1,\dots,n_r)$
of positive edge-lengths with the bar differential (`omega_words`,
`omega_boundary`); this is the model of the twist $\omega$ in every
primitive below ("the $\omega$-factor"). A pure-sign factor
(rows $(\sigma_i,0)$) is never contracted. In mode `zsign` with nonzero
edges the words are further transferred to the crossed exterior algebra of
the tertiary note; this transfer enters the stacking only through the
tertiary $V_2,V_3$ and through the degree-three pair model (Section 8.1 of
[tertiary_operations.md](tertiary_operations.md)).

A single row (an **atom** $(\sigma,e)$, $e\in\mathbb Z^q$: a $K(\mathbb Z,1)$-simplex
with sign) is contracted onto the unit edge by the *circle contraction*
(`_unit_h`): $p(\text{atom})$ keeps only the empty atom ($\mapsto t$ if
$\sigma=1$) and the single unit edge ($\mapsto e$ or $es$), $i$ is the
inclusion of these letters, $\widehat h$ replaces the last edge by the chain
of unit edges from its start to its end, $Q=1-ip$, $u=Q\widehat hQ$,
$h=u\partial u$.

## 4. The universal models and their sections

For a simplex $\sigma=(v_0,\dots,v_q)$ of $X$ with $s_{ij}=s(v_i,v_j)$:

| Model | Simplex | Section $j(\sigma)$ | Reconstruction |
| --- | --- | --- | --- |
| $\omega$-factor (`Diag`, `c2`) | symmetric binary $q\times q$ matrix $W$, rows $(0,W_i)$ | $W_{ij}=\omega(v_i,v_{i+1},v_{j+1})+\omega(v_i,v_{i+1},v_j)$ for $i<j$ (second term absent if $j=i+1$) | $\omega_{rlt}=\sum_{i=r}^{l-1}\sum_{j=l}^{t-1}W_{ij}$ |
| degree-two signed cocycle (`Diag`, `zsign`) | rows $(\sigma_i,M_i)$, $M$ integral, skew up to sign | $\sigma_i=s_{i,i+1}$; $M_{ij}=A(v_i,v_{i+1},v_{j+1})-A(v_i,v_{i+1},v_j)$ ($i<j$), $M_{ji}=-(-1)^{s_{ij}}M_{ij}$ | $s_{rt}=\sum_{i=r}^{t-1}\sigma_i$, $A_{rlt}=\sum_{i=r}^{l-1}\sum_{j=l}^{t-1}(-1)^{s_{ri}}M_{ij}$ |
| degree-one binary cocycle (`Atom`, `c2`) | edge bits $(b_1,\dots,b_q)$ | $b_i=b(v_i,v_{i+1})$ | $b_{rt}=\sum_{i=r}^{t-1}b_i$ |
| shared-sign degree-one pair (`r1_pair_chain.Borel`) | $(\sigma,\ell,r)$, $\ell,r\in\mathbb Z^q$ | $\sigma_i=s_{i,i+1}$, $\ell_i=(-1)^{\sum_{j<i}\sigma_j}A(v_i,v_{i+1})$, likewise $r_i$ from $A'$ | $A_{ij}=(-1)^{s_{0i}}\sum_{i\le m<j}\ell_m$ |
| shared-sign degree-two pair (`r2_pair_chain.Diag2`) | rows $(\sigma_i,M_i,M'_i)$ | as the signed degree-two model for $A$ and $A'$ with one common sign | as above for both |
| shared-sign degree-three pair (`r3_pair_chain.Diag3`) | rows $(\sigma_r,Z^{(r)},Z'^{(r)})$, skew integral matrices | $\bar A(i_0..i_3)=\operatorname{sgn}\cdot(-1)^{s_{0,\min}}A(\text{sorted})$; $z_r(i,j,k)=(-1)^{s_{0r}}\bigl(\bar A(r,i,j,k)-\bar A(r+1,i,j,k)\bigr)$; $Z^{(r)}_{ij}=z_r(i,i+1,j+1)-z_r(i,i+1,j)$ | $A_{ijkl}=\sum_{r=i}^{j-1}(-1)^{s_{ir}}\sum_{t=j}^{k-1}\sum_{v=k}^{l-1}Z^{(r)}_{tv}$ |
| shared-sign degree-four pair (`r4_pair_ez`) | $(\sigma,\kappa,\kappa')$ with $\kappa,\kappa'$ integral $K(\mathbb Z,4)$ cocycle simplices (values on all $5$-faces) | $\kappa_f=(-1)^{\sum_{j<f_0}\sigma_j}A(v_f)$ | direct |
| reduced binary cocycle (`a0_high_gamma.KB`) | the table of a closed binary cochain of degree $m$ on all $(m+1)$-faces of $\Delta^q$ | $c(v_f)$ on each face $f$ | direct |
| iterated binary bar (`binary_iterated_bar.Simplex`) | level-$m$ diagonal bar simplex of $K(\mathbb F_2,m)$ | the alternating section `coordinate`: level 1 is the edge section; level $m$ has rows built from $c(i,\cdot)+c(i{+}1,\cdot)$ recursively | rectangular sums |
| dihedral group model (`low_phases.to_pair1`) | homogeneous vertices $g_i=(A(v_0,v_i),s(v_0,v_i))\in\mathbb Z\rtimes C_2$, $g_0=(0,0)$ | direct | $s_{ij}=\epsilon_i+\epsilon_j$, $A_{ij}=(-1)^{\epsilon_i}(k_j-k_i)$ |

Every pair model is taken in the **relative** form: words with a zero
left or right fiber are excluded from the small basis, and the primitives
vanish identically when one fiber is zero (every term of the contraction
keeps a zero fiber). The pair models are contracted by the unit
contraction of each row as a pair of fibers sharing the sign (the tensor
product, Section 2, of the circle contractions for degree one and two;
`unit_p`, `unit_i`, `unit_h` in `r2_pair_chain`), followed by the bar
transfer of Section 3 with the crossed bar differential (the leading
letter $t=\sigma-1$ acts by $2$: `_bar_delta`), and finally by the product
with the $\omega$-factor (Section 2). Small letters are $(\sigma,k,l)$: the
sign bit and the fiber degrees ($k,l\in\{0,1\}$ for degree-one and
degree-two pairs, divided-power indices for degree three); $(1,0,0)$ is
$t$. The small basis in total degree $d$ is the set of pairs
$(\text{word},\ \omega\text{-word})$ of total degree $d$ (`total_basis`);
its boundary is `total_boundary`.

## 5. The pair primitives $\mathcal V_n$ and $\lambda_A$ (lower corrections)

Binary. Source (`v3_pair_shared.source_pair`), with $\alpha=h^D(a,a')$:

$$
\mathcal S^\beta_n(A,A')=k_0(A+A')+k_0(A)+k_0(A')+Q_D(\alpha)+h^D\bigl(P,P'\bigr)+h^D\bigl(P+P',\delta\alpha\bigr),\qquad \deg n+3 .
$$

| $n$ | Model | Coefficients |
| ---: | --- | --- |
| 1 | shared-sign degree-one pair $\times$ $\omega$-factor (`r1_pair_chain`, trivial coefficient chains) | the $\mathbb F_2$ solution of the small system on `total_basis(3)` (`v1_pair_shared.fixed_certificate`) |
| 2 | shared-sign degree-two pair $\times$ $\omega$-factor (`r2_pair_chain`) | the $\mathbb F_2$ solution on `total_basis(4)` (`v2_pair_shared.fixed_certificate`) |
| 3 | shared-sign degree-three pair $\times$ $\omega$-factor (`r3_pair_chain`) | the $\mathbb F_2$ solution of a $12\times2$ system on `total_basis(5)` (`v3_pair_shared.fixed_certificate`) |
| 4 | shared-sign degree-four pair $\times$ $\omega$-factor (`r4_pair_ez`) | none: the relative complex starts in degree eight |

$\mathcal V_n(A,A')=\mathcal S^\beta_n(H\,j(\sigma))+\sum e_w\langle F\,j(\sigma),w\rangle$ with these data; degree $n+2$.

$\lambda_A(A,A')$ (degree 4, $n=3$; `lower_raw_comparison.A_primitive`): the
same model as $\mathcal V_3$, source
$\beta_{\rm ext}(A,0;A',0)+\beta^{\rm raw}_6(A,0;A',0)$
([product_lower.md](product_lower.md) Section 7), no coefficients (the
relative complex in degree five has two middle cells, is acyclic modulo
two, and the source vanishes on both).

## 6. The comparison primitives $v_n$ and the degree-three gauges

$v_n(A)$ (degree $n+2$, binary; `h_tau_primitive.source_primitive`): the
primitive of $Z_n(A)=k_0(A)+I[k_0(\delta_s(A\ell))+Q_D(P(A\ell))]+h^D(P,P)$
on the one-fiber model of $A$ with the $\omega$-factor: $n=1$ the dihedral
group model (contraction `group_homotopy`: the Cayley-line retraction of
the tertiary note tensored with the $\omega$-factor), $n=2$ the signed
degree-two model with the tertiary $V_2$ contraction, $n=3$ the degree-three
model with the tertiary $V_3$ contraction. Its small restriction vanishes
identically, so there are no coefficients.

$\mathfrak q^{(0)}_3(B)$ (degree 3, rational; `a0_degree3.closed_pair_gauge`):
model $(b,(s,\omega))$ with $b$ a degree-one binary cocycle atom (not
contracted; its bar words are the inclusion), $s$ the pure-sign factor,
$\omega$ the $\omega$-factor; source $\Psi(0,b,0)$ of
[product_upper_k1_k3.md](product_upper_k1_k3.md) Section 6 on the
$4$-simplex; the three calibration values are integral, so there are no coefficients.

$\mathfrak q^{\rm ev}_3(r,B;\omega)$ (degree 3, rational; `even_oriented_gauge`):
model $(b,\omega)$, source $\Psi(r,b,0)$ at $s=0$ for the constant even rank
$r$; small bases in degree four $\{(0,(3)),(0,(1,1)),(1,(2)),(2,(1)),(4,())\}$
and degree three $\{(0,(2)),(1,(1)),(3,())\}$ (first entry the length of the
$b$-word, second the $\omega$-word); coefficients $(0,\ \tfrac12\{c_3\},\ \tfrac12\{c_5\})$
where $c_i$ are the source values on the degree-four basis and $\{\cdot\}$ the
fractional part (the periods $c_1,c_2,c_3-c_4$ must be integral).

## 7. The $\Theta$-pair and $A=0$ primitives

$\mathcal U^\Theta_m(q,q')$ (degree $m+2$, rational, $m=3,4,5$;
`theta_pair_phase.ThetaPairPhase.primitive`): the primitive of
$S^\Theta_m(q,q')$ of [product_upper_k4_k6.md](product_upper_k4_k6.md)
Section 1 on the model $((q,q'),(s,\omega))$ with $q,q'$ reduced binary
cocycle simplices (uncontracted leaves), $s$ pure sign, $\omega$ the
$\omega$-factor; the contraction is the product contraction of Section 2
applied recursively over the tree (`forward`, `inclusion`, `homotopy` of
`a0_high_gamma`); no coefficients (for $m=3$ the minimal period is checked
to be integral).

$\mathcal U^{A0}_k(B,B')$ (degree $k$, rational): the primitive of
$S^{A0}_k$ ([product_upper_k4_k6.md](product_upper_k4_k6.md) Section 7);
for $k=5,6$ the same model and contraction as $\mathcal U^\Theta_{k-2}$
(`a0_high_gamma.HigherA0Stacking.primitive`); for $k=4$ the model
$((b,b'),(s,\omega))$ with $b,b'$ degree-two binary cocycles as `c2`
diagonal bar simplices contracted by Section 3 (`a0_degree4`), with the
three calibration values integral, hence no coefficients.

## 8. The pair-phase primitives $\mathcal P_k$

Rational. Sources $S_4,S_5,S_6$ of [product_upper_k4_k6.md](product_upper_k4_k6.md);
models: the shared-sign pair of degree $n$ with the $\omega$-factor
(signed coefficient chains, `r1_pair_chain_signed`, `r2_pair_chain`, `r3_pair_chain`).

- $\mathcal P_4$: lower basis `total_basis(4)` $=\{((0,1,1),(1)),\ ((2,1,1),())\}$, coefficients $(\tfrac18,-\tfrac18)$ (`production_gamma4.primitive`; the auxiliary primitive of `upper_pair_source` uses $(\tfrac18,\tfrac18)$ and is not part of the product).
- $\mathcal P_5$: lower basis `r2_pair_chain.total_basis(5)` (48 elements, ordered by word degree, length, word, then $\omega$-degree, length, word); coefficients, in that order:
  $\tfrac78,\tfrac38,\tfrac{11}{12},\tfrac{11}{12},\tfrac56,\tfrac56,\tfrac56,\tfrac5{12},\tfrac56,\tfrac5{12},\tfrac5{12},\tfrac5{12},\tfrac56,\tfrac56,\tfrac{11}{12},\tfrac13,\tfrac{11}{12},\tfrac13,0,\tfrac34,\tfrac34,\tfrac38,0,\tfrac14,\tfrac14,\tfrac38,0,0,0,0,\tfrac14,\tfrac78,0,\tfrac14,\tfrac34,\tfrac34,\tfrac34,0,\tfrac18,\tfrac58,0,0,0,\tfrac58,0,0,0,0$
  (`production_gamma5.FIXED_COEFFICIENTS`; the thirds come from the three-primary phase inside $R^{(5)}_A$).
- $\mathcal P_6$: lower basis `r3_pair_chain.total_basis(6)` (12 elements), coefficients
  $\tfrac{19}{24},\tfrac7{24},\tfrac{11}{12},\tfrac5{12},0,\tfrac{11}{12},\tfrac{11}{12},\tfrac38,\tfrac38,0,0,0$ (`production_gamma6.FIXED_COEFFICIENTS`).

The coefficient vectors were obtained by the Smith-form recipe of
`universal_pair6/7` (periods on the universal cycles vanish, then
$e=(\text{coordinates})\cdot V$) and are frozen.

## 9. The primitives inside $q_4,q_5$

- $\mathfrak q^{\rm HTC}_m(b)$, $m=3,4$ (degree $m+2$; `higher_theta_comparison`): model $(b,(s,\omega))$ with $b$ an iterated binary bar simplex of level $m$, contraction of Section 3 for the iterated bar plus the transfer to small words; source $\mathrm{bp}_{m+2}(b)-\theta_m(b,0)$; lower basis of 5 elements, coefficients $(0,0,0,0,0)$ for $m=3$ and $(\tfrac14,0,0,\tfrac14,0)$ for $m=4$.
- $V^{\rm PS}(b)$ (degree 5; `primary_suspension_comparison`): the same model with $m=3$, source $\mathsf h\bigl(Q^1(Q_Db)+E(s\smile b)\bigr)$ for closed $b$, coefficients $(\tfrac14,0,0,\tfrac14,0)$.
- $\mathfrak q^A_4(A)$ (degree 4; `degree_four_calibration_primitive`): dihedral group model with the $\omega$-factor, source $\mathcal D^{(4)}(A)$ on the $5$-simplex (frame sign $(-1)^{\epsilon_0}$), retracted part with the small coefficients $(1,(2))\mapsto0$, $(2,(1))\mapsto\tfrac14$, $(4,())\mapsto\tfrac14$ on the dihedral wedge.
- $\mathfrak q^A_5(A)$ (degree 5; `degree_five_calibration_primitive`): the signed degree-two model with the tertiary $V_2$ contraction, source $\mathcal D^{(5)}(A)$ on the $6$-simplex, coefficients $(0,\tfrac38,\tfrac18,0)$ on `basis2(5)` of `h_tau_primitive`.

## 10. The tertiary source primitives $V_1,V_2,V_3$

The same construction, with the rational sources $\phi_n=\Theta_{n+2}(q_n,k_n)$
of degree $n+5$ and $\delta_sV_n\equiv\phi_n$ modulo integers: $V_1$ on the
dihedral group model (Section 8.2 of [tertiary_operations.md](tertiary_operations.md),
with the odd-branch normalization and the two-entry table $e_0$), $V_2$ on
the signed degree-two matrix model with the five source values and the
five-entry $e_f$ table (Section 8.3 there), $V_3$ on the crossed $U_2$ model
with the $46\times94$ Smith reduction, 34 source values and 16 coefficients
(Section 8.4 there). Their contractions are the ones reused by $v_2,v_3$,
$\mathfrak q^A_5$ and the degree-two and degree-three pair models above. Any
admissible contraction gives a primitive, and the finite ambiguity on the
universal model is compensated by the calibration selectors (Section 8.5
there).
