# Documentation guide and formula sheet

fermionAHSS computes the five rows `q = -4,-3,-2,-1,0` of the twisted
connective real K-theory Atiyah–Hirzebruch spectral sequence of a group
resolution through E6, and assembles the E6 layers of a degree into one
group by solving the stacking extension problem. These notes record the
mathematics behind the two steps: the fixed cochain formulas of the page
differentials, the stacking model of the extension problem, the
assumptions, and the comparisons with the literature. This page gives the
layout of the notes and a sheet of the formulas that the package
evaluates; each formula points to the note that derives it.

## Layout

```
doc/
├── README.md                          this guide and the formula sheet
│
│   Cochains and the page differentials
├── conventions.md                     twists, coefficient systems, cup-i words, lifts, primary maps
├── secondary_operations.md            d3 = Tau (row 0 to -2) and d4 = Psi (row -1 to -4)
├── tertiary_operations.md             d5 = T (row 0 to -4) in input degrees 0-3, prime-three term
├── universal_helpers.md               the fixed helpers: chi words, Theta, V1, V2, V3, selectors
├── universal_value_growth.md          growth and periodicity of the stored universal values
│
│   The stacking extension problem
├── extensions.md                      koFull: light residues, markings, target layers, prime localization, Smith bases
├── low_degree_stacking.md             the products of package degrees one and two
├── dimension_indexed_differentials.md the nonlinear differential of the stacking model, k = 0..6
├── all_cochain_differential.md        its extension to arbitrary cochains and the square-zero proof
├── ALL_COCHAIN_STACKING.md            the stacking products alpha, beta, gamma_k of every degree, k = -1..6
├── G6_GENERAL_SECTION.md              the degree-six section construction of g6 and gamma6 outside L6
├── ALL_COCHAIN_OBSTRUCTION.md         why g6 must be nonzero outside the legal locus L6
├── extension_cup_i_formulas.md        the relation classes of each layer pair as cup-i formulas on R
├── transfer.md                        evaluating the bar formulas on a supplied resolution
├── resolution-extensions.md           the native model on a resolution: transport, checks, limits
│
│   Interfaces, comparisons and status
├── backends.md                        HAP, cochain and page backends, direct operation calls
├── extension-paper-comparisons.md     the literature fixtures and what the package finds
└── mathematical-status.md             implemented range, assumptions, sources, open scope
```

Read [conventions.md](conventions.md) first: every other note uses its
cochain conventions. The page differentials are in
[secondary_operations.md](secondary_operations.md) and
[tertiary_operations.md](tertiary_operations.md), with their helpers in
[universal_helpers.md](universal_helpers.md). The extension problem is in
[extensions.md](extensions.md), whose formulas are derived in
[dimension_indexed_differentials.md](dimension_indexed_differentials.md),
[low_degree_stacking.md](low_degree_stacking.md) and
[all_cochain_differential.md](all_cochain_differential.md), with the products of
every degree in [ALL_COCHAIN_STACKING.md](ALL_COCHAIN_STACKING.md), the
degree-six section construction in
[G6_GENERAL_SECTION.md](G6_GENERAL_SECTION.md) and its necessity in
[ALL_COCHAIN_OBSTRUCTION.md](ALL_COCHAIN_OBSTRUCTION.md), with the relation
classes of each pair of layers as cup-i formulas on the resolution in
[extension_cup_i_formulas.md](extension_cup_i_formulas.md), and evaluated on
a resolution as described in [transfer.md](transfer.md) and
[resolution-extensions.md](resolution-extensions.md). The assumptions and
the limits are collected in [mathematical-status.md](mathematical-status.md).

## Names and fixed data

The notes and GAP use the same names for the differentials.

| Differential | Rows | Name | Cohomological expression |
| --- | --- | --- | --- |
| d2 | 0 to -1 | `Dbar` | \(D\rho\) |
| d2 | -1 to -2 | `D` | \(\operatorname{Sq}^2+s\operatorname{Sq}^1+\omega\) |
| d3 | -2 to -4 | `Dtilde` | \(\beta_s(\operatorname{Sq}^2+\omega)\) |
| d3 | 0 to -2 | `Tau` | secondary, formula (S11) below |
| d4 | -1 to -4 | `Psi` | secondary, formula (S7) below |
| d5 | 0 to -4 | `T` | tertiary, formula (T) below |

The binary representative of `Tau` is \(\tau'\), the integral
representative of `Psi` is \(\psi'\). The helper \(\Theta_m\) is a rational
phase, not a differential.

| Datum | Value |
| --- | --- |
| Chi family | `chi7_tail` |
| Secondary epsilon, eta | `(1,0,0)`, `(1,0,1)` |
| Tertiary low-selector vector \(\boldsymbol{\zeta}=(\zeta_1,\zeta_2,\zeta_3)\) | `(0,1,0)` |
| R2 | `R2sharp`, including the term `-A^cup3/4` |
| Final T rank correction and mu_R | `0`, `0` |
| R3 selectors c4,cN,cO,cM,epsilon_c | `(1,0,1,1,1)` |
| R3 xi | `3/4` |
| V2 source values | `(0,3/4,0,3/4,1/4)` |
| R3 suspension periods | `(3/4,1/4,0,1/2,1/2)` |
| Prime-three coefficient | `2`; on the page it contributes only in input degree three, in degree two it enters the degree-five stacking correction |

The entries of \(\boldsymbol{\zeta}\) select the R1 rank normalization and
the two R2 suspension terms; the vector is distinct from the secondary
epsilon and eta vectors and from the cochain helpers \(\zeta_{i,n}\). The
chi word lists are encoded exactly by
[data/chi-calibrated-degree7-anf.g](../data/chi-calibrated-degree7-anf.g)
and the decoding rule of [universal_helpers.md](universal_helpers.md).

## Formula sheet

### Twists, coefficients and operations

The twists are fixed cocycles \(s\in Z^1(X;\mathbf F_2)\) and
\(\omega\in Z^2(X;\mathbf F_2)\); the rows `q=0,-4` carry the sign system
\(\mathbf Z_s\) with monodromy \((-1)^s\), the rows `q=-1,-2` carry
\(\mathbf F_2\), and the row `q=-3` is zero. \(\rho\) is reduction modulo
two, \(\widetilde z\) the pointwise \(0,1\)-valued lift of a binary
cochain, \(\mathsf h(z)=\widetilde z/2\), \(\smile_i\) the interval-cut
cup-\(i\) product and \(d_s\) the coboundary of \(\mathbf Z_s\)-valued
cochains, transported from the first vertex:

\[
(d_sc)(v_0,\ldots,v_{r+1})=(-1)^{s(v_0,v_1)}c(v_1,\ldots,v_{r+1})
+\sum_{j\ge1}(-1)^jc(v_0,\ldots,\widehat{v_j},\ldots,v_{r+1}).
\]

For a binary cochain \(x\) of degree \(r\), possibly not closed,

\[
Q^j(x)=x\smile_{r-j}x+x\smile_{r-j+1}dx,\qquad
E(x)=Q^2(x)+\omega\smile x,\qquad
Q_D(x)=E(x)+s\smile Q^1(x),
\]

and for a cocycle \(a\), \(\operatorname{Sq}^j(a)=a\smile_{n-j}a\),
\(B(a)=d\widetilde a/2\), \(e(a)=\rho B(a)\). The primary maps are

\[
D(a)=\operatorname{Sq}^2(a)+s\smile e(a)+\omega\smile a,\qquad
\operatorname{Dbar}(A)=D(\rho A),\qquad
\operatorname{Dtilde}(a)=\beta_s\bigl(\operatorname{Sq}^2[a]+\omega\smile[a]\bigr),
\]

with \(\beta_s[z]=[d_s\widetilde z/2]\)
([conventions.md](conventions.md), (C5)–(C7)).

### The secondary differentials Tau and Psi

For an input \(a\) of degree \(n\) (binary for `Psi`, \(a=\rho A\) for
`Tau`) choose a defining cochain \(b\) with \(db=Da\), and set
\(u=\operatorname{Sq}^2a\), \(v=\omega a\), \(w=se\), \(e=e(a)\),
\(B=B(a)\), \(C_B=(B+\widetilde e)/2\). The common binary cochain and
integer cochain of degree \(n+3\) are

\[
\begin{aligned}
F={}&Q^2(b)+\omega\smile b+\zeta_{2,n}(\omega,a)+\chi_n(a)
 +u\smile_{n+1}v+u\smile_{n+1}w+v\smile_{n+1}w\\
 &+\zeta_{1,n+1}(s,e)+(\omega\smile_1s)\smile e+s\smile u+s^2\smile\rho C_B,
\qquad q=\widetilde\omega\smile B+B\smile_{n-1}B.
\end{aligned}
\]

With \(Z=\widetilde{s\smile u}+\widetilde{\omega\smile e}\) (an integer
sum of two lifts), the integral representative of `Psi` is

\[
\psi'=\frac{d_s\bigl(2\widetilde F+q+2Z\bigr)}4,
\qquad
\operatorname{Psi}_n(a)=\operatorname{Psi}_{{\rm ref},n}(a)
+\beta_s\bigl(s\operatorname{Sq}^2[a]\bigr)+\beta_s\bigl(\omega\operatorname{Sq}^1[a]\bigr).
\tag{S7}
\]

For `Tau`, with \(d_sA=0\), \(K=(A-\widetilde a)/2\), \(t=\rho K\),
\(x=dt\), \(y=s\smile a\), the binary
\(G=Q^2(t)+\omega\smile t+x\smile_ny+\zeta_{1,n}(s,a)+(\omega\smile_1s)\smile a+s\smile b\)
and \(L=(q-d_s\widetilde G)/2\), the binary representative is

\[
\tau'=F+\rho\Bigl(\frac{q-d_s\widetilde G}{2}\Bigr)+s^3\smile a,
\qquad
\operatorname{Tau}_n(A)=\operatorname{Tau}_{{\rm ref},n}(A)+s^3\rho A.
\tag{S11}
\]

The helpers \(\zeta_{1,n},\zeta_{2,n},\chi_n\) are the fixed words of
[universal_helpers.md](universal_helpers.md); the divisions are exact for a
valid defining system, and the code checks the projected numerators
([secondary_operations.md](secondary_operations.md)).

### The tertiary differential T

For a signed cocycle \(A\) of degree \(n\le3\) with defining cochains
\(b,c\) (\(db=Da\), \(dc=\tau_n(A,b)\)), the integral differential is the
signed coboundary of a rational phase,

\[
T_n(A;b,c)=d_s\mathcal O_n(A,b,c),\qquad
\mathcal O_n=\mathsf h(Ec)+\widehat R_n(A,b)+\mathbf 1_{n\ge2}\,\tfrac23\,\widetilde{P^1_s\rho_3A}.
\tag{T}
\]

\(\widehat R_n\) is the calibrated boundary-normalized phase of input
degree \(n\) (Sections 3–6 of
[tertiary_operations.md](tertiary_operations.md)), built from the boundary
phase \(\Theta_m\), the prism and the universal sources \(V_1,V_2,V_3\) of
[universal_helpers.md](universal_helpers.md); \(P^1_s\) is the signed
reduced power modulo three, the cube of \(\rho_3A\) in degree two and the
nineteen-term cyclic-diagonal formula in degree three (Section 7 there),
and its lift takes the values `0,1,2`. The three-primary term
\(\tfrac23\widetilde{P^1_s\rho_3A}\) has coboundary
\(2\beta_3P^1_s\rho_3A\), the k-invariant of `ko` at the prime three.

### The stacking model

A state of package degree \(k\) is \((A,B,C,D)\) with
\(A\in C^{k-3}(\mathbf Z_s)\), \(B\in C^{k-2}(\mathbf F_2)\),
\(C\in C^{k-1}(\mathbf F_2)\), \(D\in C^{k+1}(\mathbf Z_s)\), one cochain
per layer of the line `p+q=k-3`. With \(\tau_k=\tau'_{k-3}(A;B)\)
(\(\omega B\) at \(k=2\), zero below) and the phase \(\Omega_k\) of the
table, its differential is

\[
t_k=\delta C+\tau_k,\qquad
\Phi_k=\Omega_k+\mathsf h\bigl(t_k\smile_{k-1}\tau_k\bigr),\qquad
J_k=\delta_s\Phi_k-\mathsf h\bigl(E_k(t_k)\bigr),\qquad
\mathfrak d_k(A,B,C,D)=\bigl(\delta_sA,\ \delta B+Q_D(\rho A),\ t_k,\ \delta_sD+J_k\bigr),
\tag{4}
\]

| k | \(\Omega_k\) |
| ---: | --- |
| 1 | \(\mathsf h(\omega C)\) |
| 2 | \(\mathsf h(C\smile\delta C+\omega C)\) |
| 3, 4 | \(\mathsf h(E\,C)+\widehat R_{k-3}(A,B)\) |
| 5, 6 | \(\mathsf h(E\,C)+\widehat R_{k-3}(A,B)+\tfrac23\widetilde{P^1_s\rho_3A}\) |

so that a flat state, \(\mathfrak d_k=0\), has \(\delta_sD=-J_k\) and
\(J_k=T_{k-3}(A;B,C)\) when \(t_k=0\)
([dimension_indexed_differentials.md](dimension_indexed_differentials.md),
which also proves \(\mathfrak d_k^2=0\) literally, and
[all_cochain_differential.md](all_cochain_differential.md) for arbitrary
cochains). Products are triangular,

\[
(A,B,C,D)\cdot(A',B',C',D')=\bigl(A+A',\ B+B'+\alpha(A,A'),\ C+C'+\beta(\ldots),\ D+D'+\gamma_k(\ldots)\bigr),
\]

with the D-layer correction the polarization of the phase corrected by a
universal primitive,
\(\gamma_k=\Omega_k+\Omega_k'-\Omega_k(\text{sum})+\delta_s(\text{primitive})\),
so that it is integral; the degree-one and degree-two products
\(\gamma_1,\beta_2,\gamma_2\) are written out in
[low_degree_stacking.md](low_degree_stacking.md). Two flat states are gauge
equivalent when a gauge \((u,b,c,w)\) of degree \(k-1\) carries one to the
other through the product with its boundary state; the group of a degree is
the group of flat states modulo gauge, and stacking is assumed to be
commutative and associative on gauge classes.

### The extension problem

The E6 layers of degree \(j\) are

\[
A_j=E_6^{j-3,0},\qquad B_j=E_6^{j-2,-1},\qquad C_j=E_6^{j-1,-2},\qquad D_j=E_6^{j+1,-4},
\]

assembled from D upward. A generator \(q\) of order \(m\) of the next layer
has the relation \(m\widetilde q=t\) in the lower group \(H\), measured on
the power of a flat lift \(\widetilde q\), and the group is the cokernel of
the presentation

\[
R=\begin{pmatrix}R_H&0\\-T&\operatorname{diag}(m_i)\end{pmatrix}.
\]

Only the class of \(t\) in \(\operatorname{Ext}(\mathbf Z/m,H)=H/mH\)
matters: a lower generator \(e\notin mH\) is free for \(q\), and the
relation is measured only through its **target layer**, the lowest layer
with a free generator; the components below it are recorded as zero with
the certificate that exhibits them as multiples, and a relation without a
free lower generator is \(t=0\), with no lift and no measurement. Such a
row is the relation of a lift shifted within the lower group, so it is
measured through D when a later relation with a nonzero coefficient on its
generator is measured through a layer below that generator's
([extensions.md](extensions.md)).

### Localization at the primes

The layer generators have prime-power or infinite order, and the binary
layers B and C are two-groups, so a relation of an odd prime belongs to A.
The generators are grouped by prime when the layers are read, and the
relations of each prime are recorded from D upward; a free generator is a
lower generator of every prime. A row keeps only its coordinates on the
generators of its prime and the free generators: the ones it drops lie on
D generators of order prime to \(m\), which vanish after localization at
the prime (a relation of an odd prime has no B or C coordinate), so the
rows of all primes present the group. A relation whose lower group has no
generator of its prime and no free generator is \(t=0\) without a
measurement, since then \(mH=H\). Each other relation follows the model of
its prime:

- at \(p\ge5\) the rows `q=0` and `q=-4` lie in different Adams summands of
  `ko`; the relation is \(m\widetilde q=0\) without a measurement;
- at \(p=2\) the complete model above;
- at \(p=3\) the two-layer model of the rows `q=0` and `q=-4` with the
  potential \(\Omega(A)=\tfrac23\widetilde{P^1_s\rho_3A}\). In package
  degree five, where \(P^1\) is the cube, the potential is
  \(\tfrac23A\smile A\smile A\) up to an integral cochain, the correction is
  the integral polynomial

  \[
  \gamma(A,A')=-2\,(A\smile A\smile A'+A\smile A'\smile A'),
  \]

  the gauge boundary is \((\delta_su,\delta_sw)\), and the relation of a
  generator of order \(3^a\) with \(3^aA=\delta_su\) is

  \[
  x_D=2\cdot3^{a-1}A\smile A\smile A
  +\tfrac23\,\delta_s\Bigl[\sum_{j=1}^{3^a-1}\Xi(jA,A)-u\smile\delta_su\smile\delta_su\Bigr]
  \pmod{3^aH_D},
  \]

  with \(\Xi=2A(A\smile_1A')+(A\smile_1A')A+2(A\smile_1A')A'+A'(A\smile_1A')\).
  In package degree six the correction and the gauge boundary are

  \[
  \gamma(A,A')=\tfrac23\bigl[L(P^1a)+L(P^1a')-L(P^1(a+a'))+\delta_sL(\varphi(a,a'))\bigr],
  \qquad
  D_u=-\tfrac23\bigl[L(P^1\rho_3\delta_su)-\delta_sL(\chi(u))\bigr],
  \]

  with \(L\) the `0,1,2` lift and \(\varphi,\chi\) the natural cross-effect
  and coboundary primitives of the reduced power modulo three
  ([extensions.md](extensions.md), "Localization at the primes").

These give Z/9 for Z/3, Z/3 ⊕ Z/27 for Z/9 and the split groups at the
primes five and seven in degree five, the values of the representation
ring of the cyclic group, and by suspension Z/9 for Z/3 × Z in degree six.

## Reading the results

An exact invariant list uses `[]` for zero, `[0]` for Z and `[2]` for
Z/2; the display writes `0`, `.` and `?` for zero, outside the window and
unresolved. Unresolved entries are never zero. Every result keeps
`certified_ko=false`: the window is the five rows through E6, the
identification with `ko` and with the classification of fermionic phases
rests on the assumptions listed in [mathematical-status.md](mathematical-status.md).
