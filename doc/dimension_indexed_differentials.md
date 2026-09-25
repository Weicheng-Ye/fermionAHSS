# The nonlinear differential through dimension six

For the later requirement with first component \(\delta_sA\) and arbitrary
input cochains, see [all_cochain_differential.md](all_cochain_differential.md).
This note remains the explicit defining-system part of that construction.

Revised 2026-09-25 on the requested domain: the first two layer equations
hold, while the last two cochains are arbitrary. The formulas below give
an integer-valued fourth component and a literal identity
\(\mathfrak d_{k+1}\mathfrak d_k=0\). No definition of addition or stacking
is needed for this composition calculation.

## 1. Domain and conventions

Fix \(s\in Z^1(X;\mathbf F_2)\) and
\(\omega\in Z^2(X;\mathbf F_2)\). Write \(\delta\) for binary or ordinary
integral coboundary, and \(\delta_s\) for signed coboundary. In the
first-vertex frame, \(\delta_s=\delta-2\widetilde s\smile-\).
The four inputs have the requested degrees:

\[
(A_k,B_k,C_k,D_k)\in
C^{k-3}(X;\mathbf Z_s)\times C^{k-2}(X;\mathbf F_2)
\times C^{k-1}(X;\mathbf F_2)\times C^{k+1}(X;\mathbf Z_s).
\tag{1}
\]

Negative-degree cochains are zero. The domain is

\[
\boxed{\delta_sA_k=0,\qquad \delta B_k=D(\rho A_k).}
\tag{2}
\]

There is **no condition on \(C_k\) or \(D_k\)**. In particular, the formula
must retain the possibly nonzero residual \(\delta C_k+\tau\).
When \(A_k\) is absent, (2) says \(\delta B_k=0\); when \(B_k\) is also
absent, it imposes no further condition.

Here \(\rho\) is reduction modulo two, \(\widetilde x\) is the binary
integer lift, and \(\mathsf h(x)=\widetilde x/2\) is a specified rational
cochain, not a division operation in \(\mathbf Q/\mathbf Z\). Define

\[
\begin{aligned}
Q^j_r(x)&=x\smile_{r-j}x+x\smile_{r-j+1}\delta x,\\
E_r(x)&=Q^2_r(x)+\omega\smile x,\\
Q_{D,r}(x)&=E_r(x)+s\smile Q^1_r(x).
\end{aligned}
\tag{3}
\]

All three expressions are binary. For closed \(x\),
\(D(x)=Q_D(x)=\operatorname{Sq}^2x+s\operatorname{Sq}^1x+\omega x\).
Each argument of \(\mathsf h\) is reduced as a whole before lifting;
**different occurrences of \(\mathsf h\) are lifted separately**.
All products, integral signs and sign transports are those of
[conventions.md](conventions.md). Negative cup indices give zero.

## 2. One formula for every k

Let \(\tau_k\in Z^k(X;\mathbf F_2)\) and the rational
\(\Omega_k\in C^{k+1}(X;\mathbf Q_s)\) be specified in the table below.
The table uses the input variables in dimension \(k\), and
\(\tau'_n(A;B)\), \(\widehat R_n(A,B)\) retain their current calibrated
meaning. Their definitions are given in Sections 4–5.

| k | Input degrees (A,B,C,D) | \(\tau_k\) | \(\Omega_k\) |
| ---: | --- | --- | --- |
| 0 | absent, absent, absent, 1 | 0 | 0 |
| 1 | absent, absent, 0, 2 | 0 | \(\mathsf h(\omega C_1)\) |
| 2 | absent, 0, 1, 3 | \(\omega B_2\) | \(\mathsf h(C_2\smile\delta C_2+\omega C_2)\) |
| 3 | 0, 1, 2, 4 | \(\tau'_0(A_3;B_3)\) | \(\mathsf h(E_2C_3)+R_0(A_3,B_3)\) |
| 4 | 1, 2, 3, 5 | \(\tau'_1(A_4;B_4)\) | \(\mathsf h(E_3C_4)+\widehat R_1(A_4,B_4)\) |
| 5 | 2, 3, 4, 6 | \(\tau'_2(A_5;B_5)\) | \(\mathsf h(E_4C_5)+\widehat R_2(A_5,B_5)\) |
| 6 | 3, 4, 5, 7 | \(\tau'_3(A_6;B_6)\) | \(\mathsf h(E_5C_6)+\widehat R_3(A_6,B_6)+\frac23\widetilde{P^1_s\rho_3A_6}\) |

The last lift has values `0,1,2`. The degree-four signed reduced power
\(P^1_s\) is the 19-term operation in
[tertiary_operations.md](tertiary_operations.md), Section 7.

For \(1\le k\le6\), put

\[
\boxed{\begin{aligned}
t_k&=\delta C_k+\tau_k\in Z^k(X;\mathbf F_2),\\
\Phi_k&=\Omega_k+\mathsf h\bigl(t_k\smile_{k-1}\tau_k\bigr),\\
J_k&=\delta_s\Phi_k-\mathsf h\bigl(E_k(t_k)\bigr),\\
\mathfrak d_k(A_k,B_k,C_k,D_k)
  &=\bigl(0,0,t_k,\ \delta_sD_k+J_k\bigr).
\end{aligned}}
\tag{4}
\]

Both extra terms in (4) are needed when \(t_k\ne0\). Section 6 proves
that \(J_k\) is an **integer** cochain of degree \(k+2\), and that this
map squares to zero as a literal cochain identity. All four outputs have
the required degrees at \(k+1\).

At \(k=0\), the formula is simply

\[
\mathfrak d_0(0,0,0,D_0)=(0,0,0,\delta_sD_0).
\tag{5}
\]

The optional package degree \(k=-1\) uses the same expression with
\(D_{-1}\in C^0(X;\mathbf Z_s)\).

When \(t_k=0\), the two corrections in (4) vanish literally. Thus
\(J_2=\psi'_0(B_2;C_2)\), and
\(J_k=T_{k-3}(A_k;B_k,C_k)\) for \(3\le k\le6\), with the same integral
representatives as the current formulas. The choice of plus sign in
\(\delta_sD_k+J_k\) makes a zero of the fourth component satisfy
\(\delta_sD_k=-J_k\).

## 3. Expanded low-degree and high-degree answers

For \(k=1\), set \(t_1=\delta C_1\). Formula (4) reads

\[
\boxed{\mathfrak d_1(0,0,C_1,D_1)=
\left(0,0,\delta C_1,
\delta_sD_1+\frac{
\delta_s\widetilde{\omega C_1}-\widetilde{\omega\delta C_1}}2\right).}
\tag{6}
\]

For \(k=2\), \(\delta B_2=0\). Set
\(p=\omega B_2\), \(t_2=\delta C_2+p\). Then

\[
\boxed{\begin{aligned}
\mathfrak d_2(0,B_2,C_2,D_2)
  =\bigl(0,0,t_2,\ \delta_sD_2
  &+\delta_s\{\mathsf h(C_2\smile\delta C_2+\omega C_2)
                    +\mathsf h(t_2\smile_1p)\}\\
  &-\mathsf h(t_2\smile t_2+\omega t_2)\bigr).
\end{aligned}}
\tag{7}
\]

The four remaining cases can be written without leaving any cup index
implicit:

\[
\boxed{\begin{aligned}
\mathfrak d_3&=(0,0,t_3,\ \delta_sD_3+
 \delta_s[\Omega_3+\mathsf h(t_3\smile_2\tau'_0)]
 -\mathsf h(t_3\smile_1t_3+\omega t_3)),\\
\mathfrak d_4&=(0,0,t_4,\ \delta_sD_4+
 \delta_s[\Omega_4+\mathsf h(t_4\smile_3\tau'_1)]
 -\mathsf h(t_4\smile_2t_4+\omega t_4)),\\
\mathfrak d_5&=(0,0,t_5,\ \delta_sD_5+
 \delta_s[\Omega_5+\mathsf h(t_5\smile_4\tau'_2)]
 -\mathsf h(t_5\smile_3t_5+\omega t_5)),\\
\mathfrak d_6&=(0,0,t_6,\ \delta_sD_6+
 \delta_s[\Omega_6+\mathsf h(t_6\smile_5\tau'_3)]
 -\mathsf h(t_6\smile_4t_6+\omega t_6)).
\end{aligned}}
\tag{8}
\]

Here each \(\tau'_n\) has the arguments \((A_{n+3};B_{n+3})\), and
\(t_{n+3}=\delta C_{n+3}+\tau'_n\). In the table,

\[
\begin{aligned}
E_2C_3&=C_3\smile C_3+C_3\smile_1\delta C_3+\omega C_3,\\
E_3C_4&=C_4\smile_1C_4+C_4\smile_2\delta C_4+\omega C_4,\\
E_4C_5&=C_5\smile_2C_5+C_5\smile_3\delta C_5+\omega C_5,\\
E_5C_6&=C_6\smile_3C_6+C_6\smile_4\delta C_6+\omega C_6.
\end{aligned}
\tag{9}
\]

## 4. The explicit Tau expressions used in all four cases

Use \(n=k-3\), \(A=A_k\), \(b=B_k\), \(a=\rho A\). The following
single formula is evaluated at \(n=0,1,2,3\). To keep the letter \(B_k\)
free for the second layer, write \(\mathcal B\) for its auxiliary
Bockstein cochain:

\[
\begin{gathered}
\mathcal B=\frac{\delta\widetilde a}{2},\quad e=\rho\mathcal B,
\quad \mathcal C_B=\frac{\mathcal B+\widetilde e}{2},\\
u=a\smile_{n-2}a,\quad v=\omega a,\quad w=se,
\qquad p=u+v+w=\delta b,\\
t=\rho\frac{A-\widetilde a}{2},\quad x=\delta t=e+sa,\quad y=sa.
\end{gathered}
\tag{10}
\]

The auxiliary \(t\) in (10) is the parity of the integral carry, distinct
from the third-layer residual \(t_k\) in (4). Define

\[
\begin{aligned}
F={}&E_{n+1}(b)+\zeta_{2,n}(\omega,a)+\chi_n(a)\\
 &+u\smile_{n+1}v+u\smile_{n+1}w+v\smile_{n+1}w\\
 &+\zeta_{1,n+1}(s,e)+(\omega\smile_1s)e+su+s^2\rho\mathcal C_B,\\
q={}&\widetilde\omega\smile\mathcal B+
                         \mathcal B\smile_{n-1}\mathcal B,\\
G={}&E_n(t)+x\smile_n y+\zeta_{1,n}(s,a)
          +(\omega\smile_1s)a+sb.
\end{aligned}
\tag{11}
\]

Then the required expression is

\[
\boxed{\tau'_n(A;b)=F+
\rho\left(\frac{q-\delta_s\widetilde G}{2}\right)+s^3a.}
\tag{12}
\]

All binary sums are reduced before lifting. The products containing
\(\mathcal B\) in (11) have the prescribed integral cup signs. The plus
in \(\mathcal C_B\) is intentional. For \(n=0\),
\(\zeta_{1,0}=\zeta_{2,0}=\chi_0=0\). The full fixed word and ANF
coefficients are specified in
[universal_helpers.md](universal_helpers.md), Sections 1–2.
Equation (12) is exactly (S11) of
[secondary_operations.md](secondary_operations.md), with `chi7_tail` and
secondary epsilon \((1,0,0)\). In particular \(\delta\tau'_n=0\)
under (2).

## 5. The phases, with k=5 and k=6 written out

The known \(k=3,4\) phases use \(R_0\) and \(\widehat R_1\) from
[tertiary_operations.md](tertiary_operations.md), Sections 3–4, equations
(R0e), (R0o), and (R1), including the even-input carry and odd-input
signed gauge. These are explicit fixed formulas, not variables to solve.
The following gives the remaining \(k=5,6\) phases in the same convention.

For \(n=2,3\), continue using (10). Split off the \(b\)-independent
part of \(F\) and the binary \(g\) by

\[
\begin{gathered}
F_A=F+E_{n+1}(b),\qquad g=G+sb,\\
L_A=\frac{q-\delta_s\widetilde g+\widetilde s\smile\widetilde p}{2},
\qquad \kappa_n=F_A+\rho L_A+s^3a,\\
y_n=Q_D(b),\qquad
\lambda_n=\rho\frac{\widetilde g+\widetilde s\smile\widetilde b
                                  -\widetilde{g+sb}}2,\\
\tau'_n=y_n+\kappa_n+\delta\lambda_n,\qquad
\delta\kappa_n=D(p).
\end{gathered}
\tag{13}
\]

Here \(F_A,g,\kappa_n,y_n,\lambda_n\) are binary; \(L_A\) is integral.
Define

\[
\begin{gathered}
\operatorname{Pol}_n(y,\kappa)=
 y\smile_{n+2}\kappa+\kappa\smile_{n+3}\delta y+Q^1\kappa,\\
\Pi_n=I\Theta_{n+2}\bigl(\delta(b\ell),Q_D(b\ell)\bigr),\qquad
P_\omega=\widetilde\omega\smile\widetilde\omega+
           \widetilde\omega\smile_1\delta\widetilde\omega.
\end{gathered}
\tag{14}
\]

The right interval coordinate is \(\ell\), with
\((b\ell)(z)=b(z)\ell(z_{\rm last})\). On a degree-\(r\) cochain, the
raw right prism is the fixed finite sum

\[
(Iu)(v_0,\ldots,v_{r-1})=
\sum_{j=0}^{r-1}(-1)^j
u((v_0,0),\ldots,(v_j,0),(v_j,1),\ldots,(v_{r-1},1)).
\tag{15}
\]

The phase \(\Theta_m\) in (14) is explicitly

\[
\Theta_m(z,v)=
\mathsf h(Ev+H_m(z))+
\frac14\bigl(\widetilde\omega\smile\mathcal B_z+
                     \mathcal B_z\smile_{m-1}\mathcal B_z\bigr)
+\mathsf h(s\operatorname{Sq}^2z+\omega\operatorname{Sq}^1z),
\tag{16}
\]

where \(z\) is closed, \(\delta v=Dz\),
\(\mathcal B_z=\delta\widetilde z/2\), and \(H_m(z)\) is the expression
\(F_A\) in (11),(13) with input \(z\) in degree \(m\).
\(V_n\) denotes the fixed calibrated source primitive of
\(\Theta_{n+2}(p,\kappa_n)\), evaluated by the contractors and full finite
source tables in [universal_helpers.md](universal_helpers.md), Sections
6–7. These data fix every coefficient; no primitive on \(X\) is selected.

For **k=5**, take \(n=2\), \(A=A_5\), \(b=B_5\), \(C=C_5\), and put
\(U_A=\zeta_{1,4}(s,p)+s\kappa_2+(\omega\smile_1s)p\). The complete
phase to use in (8) is

\[
\boxed{\begin{aligned}
\Omega_5={}&\mathsf h(C\smile_2C+C\smile_3\delta C+\omega C)\\
 &+\mathsf h(\operatorname{Pol}_2(y_2,\kappa_2))+\Pi_2-V_2\\
 &+\mathsf h\bigl(E\lambda_2+(y_2+\kappa_2)\smile_4\delta\lambda_2\bigr)\\
 &-\mathsf h(E(sb))+\mathsf h(U_A)
   +\frac14 P_\omega\smile_s A-\frac14 A^{\smile3}.
\end{aligned}}
\tag{17}
\]

The last product is the signed ordered cube, with its coefficient-system
transports. The cubic correction is included exactly once.

For **k=6**, take \(n=3\), \(A=A_6\), \(b=B_6\), \(C=C_6\). Define
\(N_3=\mathsf h(s^2\operatorname{Sq}^2a)\),
\(O_3=\mathsf h(\omega\operatorname{Sq}^2a)\), and
\(M_3=\mathsf h(s^2\omega a)\). Then

\[
\boxed{\begin{aligned}
\Omega_6={}&\mathsf h(C\smile_3C+C\smile_4\delta C+\omega C)\\
 &+\mathsf h\bigl(\operatorname{Pol}_3(y_3,\kappa_3)+E\lambda_3
                  +(y_3+\kappa_3)\smile_5\delta\lambda_3\bigr)
       -\Pi_3-V_3\\
 &+\mathsf h(EQ^1b+s\tau'_3)
   +\frac14 P_\omega\smile_s A+N_3+2O_3+M_3
   +\frac23\widetilde{P^1_s\rho_3A}.
\end{aligned}}
\tag{18}
\]

Keep the integer term \(2O_3\) when computing literal representatives.
These formulas use tertiary low selectors
\(\boldsymbol\zeta=(\zeta_1,\zeta_2,\zeta_3)=(0,1,0)\), secondary eta
\((1,0,1)\), zero rank correction, \(\mu_R=0\), and the coefficient-two
prime-three term. The symbols \(O_3\) in (18), \(\Omega_k\), and the
phase helper \(\Theta_m\) are distinct; `Psi` continues to mean the
integral degree-four secondary operation.

## 6. Explicit proof that d squared is zero

The proof consists of integrality, invariance of the domain, and literal
cancellation. It applies to every \(X,s,\omega\) on the domain (2).

The two cochain identities needed are

\[
\delta E(x)=E(\delta x),\qquad
E(u+v)+E(u)+E(v)=\delta(u\smile_{k-1}v)
\quad(u,v\in Z^k(X;\mathbf F_2)).
\tag{19}
\]

Both follow by expanding the cup-i boundary identity; the \(\omega\)
terms cancel in the second identity. In phase coefficients,
\(\delta_s\mathsf h(z)=\mathsf h(\delta z)\), since their difference
is integral. Signed transport causes no extra phase term on order-two
values.

The established lower/tertiary identities give

\[
\delta\tau_k=0,\qquad
\delta_s\Omega_k\equiv
\mathsf h(E(\delta C_k))+\mathsf h(E\tau_k)
\pmod{C^{k+2}(X;\mathbf Z_s)}.
\tag{20}
\]

For \(k\ge3\), (20) is the residual identity
\(\delta_s\widehat R_{k-3}\equiv-\mathsf h(E\tau_k)\); a half-valued
phase equals its negative. The prime-three phase in (18) has integral
coboundary because \(\rho_3A\) is a signed cocycle. The residual identity
also follows locally from the existing formula's integrality with
\(\delta c=\tau_k\): such a \(c\) exists on each simplex, so no global
primitive of \(\tau_k\) is being assumed.

At \(k=1\), \(\tau_1=0\), so (20) follows directly from (19). At
\(k=2\), \(\tau_2=\omega B_2\) is closed and

\[
E(\omega B_2)=(\omega B_2)^2+\omega^2B_2=0
\tag{21}
\]

literally, because \(B_2\) is locally constant and \(B_2^2=B_2\).
Thus (20) holds in the low degrees too.

Since \(t_k=\delta C_k+\tau_k\), both \(t_k\) and \(\tau_k\) are
closed. Substitute \(u=t_k,v=\tau_k\) in (19):

\[
\delta(t_k\smile_{k-1}\tau_k)
 =E(t_k)+E(\delta C_k)+E(\tau_k).
\tag{22}
\]

Combining (20) and (22) cancels the last two terms twice, proving

\[
\delta_s\Phi_k\equiv\mathsf h(E t_k)\pmod{\mathbf Z_s}.
\]

Consequently \(J_k=\delta_s\Phi_k-\mathsf h(E t_k)\) is integral.
Now take its signed coboundary **as an exact rational equality**:

\[
\boxed{\delta_sJ_k=-\delta_s\mathsf h(E t_k).}
\tag{23}
\]

Write the first output as \((0,0,t_k,D')\), where
\(D'=\delta_sD_k+J_k\). It satisfies (2) automatically. The pointed
normalizations of the helpers give the pure-C rule

\[
\mathfrak d_{k+1}(0,0,t_k,D')
 =\bigl(0,0,0,\delta_sD'+\delta_s\mathsf h(E t_k)\bigr),
\qquad \delta t_k=0.
\tag{24}
\]

Indeed, \(\tau(0,0)=0\), and the signed boundary of the corresponding
\(R(0,0)\) is zero; all lower source terms are reduced at the zero input.
Thus the second application has fourth component

\[
\begin{aligned}
\delta_sD'+\delta_s\mathsf h(E t_k)
 &=\delta_s^2D_k+\delta_sJ_k+\delta_s\mathsf h(E t_k)\\
 &=0-\delta_s\mathsf h(E t_k)+\delta_s\mathsf h(E t_k)=0.
\end{aligned}
\tag{25}
\]

Its first three components are also zero. This proves

\[
\boxed{\mathfrak d_{k+1}\mathfrak d_k(A_k,B_k,C_k,D_k)=(0,0,0,0).}
\tag{26}
\]

For the endpoint \(k=6\), only (24) on the image in degree seven is
needed. It is the ordinary primary
\(\operatorname{Dtilde}(t_6)=\delta_s\widetilde{E(t_6)}/2\), not a
new \(R_4\) or \(T_4\). More generally the pure-C rule on arbitrary
cochains is available in every degree as

\[
\mathfrak d(0,0,C,D)=
\left(0,0,\delta C,\ \delta_sD+
\frac{\delta_s\widetilde{E(C)}-\widetilde{E(\delta C)}}2\right).
\tag{27}
\]

At \(k=0\), the second iterate vanishes directly by \(\delta_s^2=0\).

## 7. Explicit verification and source record

The symbolic proof (19)–(27) establishes the identity in general.
Exact-arithmetic diagnostics additionally checked the representative
conventions, correction indices, integrality, and signs:

- 100 generic-twist simplex inputs each for \(k=1,2\): all integrality
  and second-iterate checks passed.
- 12 cases for \(k=3,4,5,6\): 96 integer boundary-face evaluations and
  12 literal second iterates passed. Six cases had nonzero \(B\),
  \(\tau\), and the new cup correction, with nonclosed \(C\).
- The degree-zero tertiary cases include signed odd \(A=\pm3\), even
  \(A=\pm2\), and even \(A=\pm6\) with a nonzero quarter-input carry.
- The full fixed production phase assemblers were used. For \(k=3,4,5\),
  the second application also evaluated the next production assembler at
  \((0,0,t_k)\); for \(k=6\), it used the explicit rule (24).

The higher-degree samples deliberately use sparse signed \(A\). Their
cubic and prime-three terms can vanish, so these are checks of the new
square-zero formula, not an independent audit of every calibration term.
No calibration data or runtime operation was changed.

Reproducible scripts, full inputs, rational values, and results are under
the [local verification archive](/Users/victor/Documents/miscellaneous/fermionAHSS/20260925-nonlinear-square-zero/README.md).
The archive has been moved outside the repository; the results and their
scope are retained here. This consolidation did not rerun calculations.

### 7.1. Separate evidence for arbitrary lower inputs

The same archive also contains checks of the later
[all-cochain extension](all_cochain_differential.md), whose first output
is \(\delta_sA\). Those checks do **not** enlarge the domain (2) of
the defining-system formula (4):

- The first all-cochain report has five signed, nonzero-\(\omega\)
  inputs in \(k=2,3,4,5,6\). The \(k=2\) input has nonclosed \(B\);
  the others have nonclosed \(A\). In \(k=2\) through \(5\), all 200
  cylinder T values and 30 output \(g\) face values were integral,
  and all four second iterates were exactly zero. There were 31 nonzero
  cylinder values. Two next-T values were 1 and \(-2\), so the checks
  include cancellation of nonzero summands.
- The dense-input audit independently populated integer \(A,D\) and
  binary \(B,C\), including mixed signs and parity in \(A\), nonclosed
  \(D\), and a family with signed-closed \(A\) but invalid \(B\).
  Six full compositions passed, two each in \(k=2,3,6\), with 144
  integral cylinder values, of which 70 were nonzero. Four further
  \(k=4,5\) cases checked only the lower three components and
  \(\delta_s^2D\); their top components remain **uncomputed**, not
  passing full checks. Three recorded top-evaluation attempts timed
  out or were interrupted after a preceding timeout.
- The arbitrary-input \(k=6\) entries in these two reports used the
  **historical zero branch outside the defining domain** (including
  ten zero output face values in the first report). They are not tests
  of the repaired \(g_6=T_3(\mathcal R_6u)\) branch. The current repair
  and its separate verification are recorded in
  [all_cochain_differential.md](all_cochain_differential.md). The
  defining-domain \(k=6\) tests above still apply to the unchanged
  \(J_6\), including its pure-C successor rule (24).

The archived scripts were code-reviewed; the low-degree report and the
\((\mathbf Z/4)^4\) report below were independently reproduced, and
the all-cochain \(k=2\) case was rerun after its H2 expression was
aligned with (A10) of the all-cochain note. Finite checks neither
establish naturality nor replace the algebraic square-zero proofs. The
sparse high-degree samples do not independently test the cubic or
prime-three calibration terms.

The current representatives are those in
[secondary_operations.md](secondary_operations.md),
[tertiary_operations.md](tertiary_operations.md), and
[universal_helpers.md](universal_helpers.md). The source comparison also
read [Wang–Gu, arXiv:1811.00536v3](https://arxiv.org/html/1811.00536v3)
and the supplied companion PDFs:
[generic parity](../../fermionAHSS_math/O5_generic_symmetry_multilayers.pdf),
[degree-four/five progress](../../fermionAHSS_math/O4_O5_progress_report.pdf),
[untwisted O6](../../fermionAHSS_math/O6_detailed_note.pdf), and
[two-layer p+ip](../../fermionAHSS_math/double_pip_O3_O4_obstruction_note.pdf).
The new step here is (4), which extends the last two components to arbitrary
\(C,D\) on the specified lower-stage domain and verifies their literal
square-zero composition.

## 8. Nonzero obstruction for A=2 on B(Z/4)^4

This example supplies the classification check: the first three equations
can be solved, but the fourth cannot, for **every** allowed lower choice.
Take \(G=(\mathbf Z/4)^4\), \(s=0\), and
\([\omega]=x_1x_2+x_3x_4\). Write the coordinates of group elements in
\(\{0,1,2,3\}\), with addition modulo four, and set

\[
\begin{gathered}
F(g,h)=g_1h_2+g_3h_4,\qquad w=F\bmod2\in\{0,1\},
\qquad \omega=\rho_2w,\\
\boxed{A_3=2,\qquad B_3=0,\qquad
C_3=C_0=\rho_2\bigl((F-w)/2\bigr).}
\end{gathered}
\tag{28}
\]

Here \(F\) is an integer cochain, distinct from the helper \(F\) in
(11). Its reduction modulo four is bilinear, so \(\delta F\in4\mathbf Z\).
The calibrated Tau formula (12) therefore gives

\[
\tau'_0(2;0)=\rho_2(-\delta w/2)
 =\omega\smile_1\omega=\delta C_0.
\tag{29}
\]

Thus (28) solves all three lower equations. The right side of (29) is
not identically zero: it is one on \([e_1|e_1|e_2]\), so setting
\(C_3=0\) would lose the required defining cochain. Here and below,
\(e_i\) denotes the generator of the ith cyclic factor.

Using the actual even-input production phase (R0e), with its full carry
equal to one and quarter-input carry zero, gives

\[
\begin{gathered}
P_\omega=w\smile w+w\smile_1\delta w,\qquad
E_2(C)=C\smile C+C\smile_1\delta C+\omega\smile C,\\
\boxed{\Omega_3(C)=\mathsf h(E_2(C))-\frac18P_\omega,
\qquad J_3(2,0,C)=g_3(2,0,C)=\delta\Omega_3(C)}
\quad(\delta C=\tau'_0(2;0)).
\end{gathered}
\tag{30}
\]

The third-layer residual is zero here, so (4) contributes no residual
correction. In particular this example is independent of the later
all-cochain extension away from the defining domain.

The oriented bar cycle

\[
z=\sum_{\sigma\in S_4}\operatorname{sgn}(\sigma)
[e_{\sigma(1)}|e_{\sigma(2)}|e_{\sigma(3)}|e_{\sigma(4)}]
\tag{31}
\]

has 24 terms and zero boundary. In the repository's cup conventions,
\(\langle\widetilde{E_2(C_0)},z\rangle=0\) and
\(\langle P_\omega,z\rangle=2\): the ordinary cup square contributes
only the two positive permutations \((1,2,3,4)\), \((3,4,1,2)\), and
the cup-one correction contributes zero. Hence

\[
\boxed{\langle\Omega_3(C_0),z\rangle=-\frac14,
\qquad \exp(2\pi i\langle\Omega_3(C_0),z\rangle)=-i\ne1.}
\tag{32}
\]

The supplied two-layer p+ip PDF uses a positive quadratic-phase
convention and gives \(+i\) for its canonical choice. Both conclusions
are nonzero; no equality of their complete cochain representatives is
asserted. Since \(G\) is finite, its positive-degree rational cohomology
vanishes. The coefficient connecting map

\[
\beta:H^4(BG;\mathbf Q/\mathbf Z)\xrightarrow{\cong}H^5(BG;\mathbf Z)
\tag{33}
\]

therefore takes \([\Omega_3]\) to the nonzero \([g_3]\). The cycle
detects the degree-four phase, not the integral degree-five cocycle
directly; its pullback to \(T^4\) in degree five is zero. Equation (32)
proves that no integral \(D_3\) solves \(\delta D_3+g_3=0\).

This conclusion survives all lower-cochain choices. Every possible
\(B\) is closed and
\(\tau'_0(2;B)=\omega\smile B+\rho_2(-\delta w/2)\). In

\[
H^*(BG;\mathbf F_2)=\Lambda(x_1,x_2,x_3,x_4)
 \otimes\mathbf F_2[y_1,y_2,y_3,y_4],\qquad |x_i|=1,\quad|y_i|=2,
\tag{34}
\]

multiplication by \(\omega\) sends the four \(x_i\) to the four
independent exterior triples. Solvability for \(C\) thus forces
\([B]=0\). Normalized group-bar degree-one coboundaries with trivial
coefficients are zero, so \(B=0\) literally. This also removes the Psi
indeterminacy, because \(\ker(D:H^1\to H^3)=0\).

Every remaining choice has \(C=C_0+b\) for a binary two-cocycle \(b\).
The first term of (30) always pairs with \(z\) to a half-integer,
whereas its second term always pairs to \(-1/4\). Consequently

\[
\langle\Omega_3(C_0+b),z\rangle\in-\tfrac14+\tfrac12\mathbf Z,
\qquad O_4(z)\in\{i,-i\}.
\tag{35}
\]

This includes arbitrary exact shifts of \(C\), so the Dtilde
indeterminacy cannot kill the class either. In the actual page quotient,
\(\boxed{T_0(2)=d_5(2)\ne0}\): the physical dimension is \(k=3\),
the page differential is \(d_5\), and the PDF calls the same obstruction
\(O_4\) by the degree of its phase.

The archived production-formula check verified all 24 detecting
simplices, the zero integer bar boundary, and the exact periods in
(32). It also checked 128 sampled lower equations and 32 integral top
values. All 1024 shifts by the ten basis classes \(x_ix_j\) and \(y_i\)
of \(H^2(G;\mathbf F_2)\) were evaluated: 512 had phase residue
\(1/4\), and 512 had residue \(3/4\); none was zero. These recorded
checks support (28)–(35); the argument above covers all cochain choices.
They are not a full GAP AHSS table calculation for this 256-element
group. No runtime source or calibration data were changed.
