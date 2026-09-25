# An explicit all-cochain square-zero extension through k=6

2026-09-25. This answers the revised requirement with first component
\(\delta_s A\), arbitrary four input cochains, and a third-component
correction independent of \(C\). It extends the calibrated formulas in
[dimension_indexed_differentials.md](dimension_indexed_differentials.md).

**Scope of the construction.** The formulas below are piecewise maps on
cochains. The branch is selected by whether the first two equations hold on
all of \(X\), not separately on each simplex. They square to zero on every
input and retain exactly the repository's obstruction equations on
defining systems. This global branch need not commute with restriction or
pullback. Thus this is an explicit algebraic extension with the requested
obstructions, not a claim to have derived a natural local all-cochain
formula. The revised degree-six branch additionally uses a fixed pointed
section of the lower differential, as specified in
[the degree-six repair](g6_repair.md). The algebraic construction works
for any such choice on X; the implementation computes the section on a
finite free cochain model. A separate degree-seven endpoint
convention is stated below; it is not a degree-seven ko operation.

## 1. Cochains and fixed operations

Fix binary cocycles \(s\in Z^1(X;\mathbf F_2)\) and
\(\omega\in Z^2(X;\mathbf F_2)\), with the ordered simplicial cochain
and signed-lift conventions of [conventions.md](conventions.md). Set

\[
(A,B,C,D)\in C^{k-3}(X;\mathbf Z_s)\times
C^{k-2}(X;\mathbf F_2)\times C^{k-1}(X;\mathbf F_2)
\times C^{k+1}(X;\mathbf Z_s).
\tag{A1}
\]

Negative degrees are zero. The following operations apply to arbitrary
binary cochains \(x\) of degree \(r\):

\[
\begin{aligned}
Q^j_r(x)&=x\smile_{r-j}x+x\smile_{r-j+1}\delta x,\\
E_r(x)&=Q^2_r(x)+\omega\smile x,\\
Q_{D,r}(x)&=E_r(x)+s\smile Q^1_r(x),\\
P_k(A)&=Q_{D,k-3}(\rho A),\qquad
\mathsf h(x)=\widetilde x/2.
\end{aligned}
\tag{A2}
\]

Here \(\rho\) is reduction modulo two; each \(\mathsf h\) lifts its
whole binary argument separately. Negative cup indices give zero.
The cup identities give

\[
\delta Q_D(x)=Q_D(\delta x),\qquad
\delta E(x)=E(\delta x).
\tag{A3}
\]

Define the first two output components and their vanishing locus by

\[
p_k(A,B)=(\delta_s A,\ \delta B+P_k(A)),
\qquad L_k=\{(A,B):p_k(A,B)=0\}.
\tag{A4}
\]

The symbol \(p_k(A,B)=0\) means equality of entire cochains on \(X\).
Equation (A3) implies \(p_{k+1}p_k=0\) exactly. In particular every
output pair belongs to \(L_{k+1}\), even when its input is arbitrary.

On \(L_k\), let \(\tau_k\) and the rational phase \(\Omega_k\) be
the following fixed expressions:

| k | \(\tau_k(A,B)\) | \(\Omega_k(A,B,C)\) |
| ---: | --- | --- |
| 0 | 0 | 0 |
| 1 | 0 | \(\mathsf h(\omega C)\) |
| 2 | \(\omega B\) | \(\mathsf h(C\smile\delta C+\omega C)\) |
| 3 | \(\tau'_0(A;B)\) | \(\mathsf h(E_2C)+R_0(A,B)\) |
| 4 | \(\tau'_1(A;B)\) | \(\mathsf h(E_3C)+\widehat R_1(A,B)\) |
| 5 | \(\tau'_2(A;B)\) | \(\mathsf h(E_4C)+\widehat R_2(A,B)\) |
| 6 | \(\tau'_3(A;B)\) | \(\mathsf h(E_5C)+\widehat R_3(A,B)+\frac23\widetilde{P^1_s\rho_3A}\) |

The full finite \(\tau'_n\) expression is (10)–(12) of the
[defining-system formula note](dimension_indexed_differentials.md).
Its \(\Omega_5,\Omega_6\) are explicitly expanded in (17)–(18) there.
In particular \(\widehat R_2\) includes \(-A^{\smile3}/4\) exactly
once, and \(\Omega_6\) includes the coefficient-two prime-three term.
These are the calibrated `chi7_tail` operations, with their specified
integral carries and universal source tables; no new primitive on
\(X\) is chosen.

We also use \(\tau_7=\tau'_4\), and, only to define the endpoint map
on arbitrary degree-seven inputs, \(\tau_8=\tau'_5\). Both are in the
documented secondary range. On their defining domains,
\(\delta\tau_k=0\) and \(\tau_k(0,0)=0\).

On \(L_k\), with \(C\) still arbitrary, put

\[
\begin{aligned}
t_k&=\delta C+\tau_k(A,B),\\
J_k(A,B,C)&=
\delta_s\left[\Omega_k(A,B,C)
 +\mathsf h(t_k\smile_{k-1}\tau_k(A,B))\right]
-\mathsf h(E_k(t_k)).
\end{aligned}
\tag{A5}
\]

The defining-system note proves that \(J_k\) is integral for every
\((A,B)\in L_k\), without requiring \(t_k=0\). Its exact boundary is

\[
\delta_sJ_k=-\delta_s\mathsf h(E_k(t_k)).
\tag{A6}
\]

On a full defining system, where also \(\delta C=\tau_k\), write

\[
\mathcal K_k(A,B,C)=\delta_s\Omega_k(A,B,C)=J_k(A,B,C).
\tag{A7}
\]

This is an integral cocycle: \(\mathcal K_2=\psi'_0\), and
\(\mathcal K_k=T_{k-3}\) for \(3\le k\le6\), at the level of
the repository's specified cochain representatives.

## 2. Two finite prism formulas

Let \(\ell\) be the right interval coordinate on \(X\times I\).
All backgrounds are pulled back from \(X\). For a cochain \(u\),
\(u\ell\) means its pullback multiplied by the interval coordinate
of the last vertex. Pullbacks use normalized cochains, so projected
degenerate simplices contribute zero. The prism operator on a cochain
\(v\) of degree \(r\) is the explicit finite sum

\[
(Iv)(v_0,\ldots,v_{r-1})=
\sum_{j=0}^{r-1}(-1)^j
v((v_0,0),\ldots,(v_j,0),(v_j,1),\ldots,(v_{r-1},1)).
\tag{A8}
\]

Use the same formula with pulled-back local coefficients. Vertical sign
transport is trivial. Its identity is
\(\delta_s I+I\delta_s=i_1^*-i_0^*\); modulo two use \(\delta\).

For \(3\le k\le7\), define the binary cochain

\[
\boxed{H_k(A,B)=
 I\,\tau_{k+1}\bigl(p_k(A\ell,B\ell)\bigr).}
\tag{A9}
\]

The argument of \(\tau_{k+1}\) is always legal since \(p^2=0\).
Thus (A9) evaluates an existing fixed formula on a valid pair; it is
not a request to solve a secondary equation. For \(k=2\), use the
simpler expression

\[
\boxed{H_2(B)=Q_{D,0}(B)
=\omega\smile B+s\smile(B\smile\delta B).}
\tag{A10}
\]

Both definitions satisfy

\[
\delta H_k(A,B)=\tau_{k+1}(p_k(A,B)),\qquad H_k(0,0)=0.
\tag{A11}
\]

For (A9), this follows directly from prism Stokes and
\(\delta\tau_{k+1}=0\). For (A10), use (A3) and the literal identity
\(\tau'_0(0;b)=Q_D(b)\) for closed \(b\).

For \(2\le k\le5\), define on \(X\times I\)

\[
\begin{aligned}
(a_I,b_I)&=p_k(A\ell,B\ell),\\
c_I&=\delta(C\ell)+H_k(A\ell,B\ell),\\
\boxed{G_k(A,B,C)&=-I\,\mathcal K_{k+1}(a_I,b_I,c_I).}
\end{aligned}
\tag{A12}
\]

When evaluating \(H_k\) in the second line, (A9) uses a **second,
independent interval**. Thus (A12) is an explicit double-prism formula
for \(k\ge3\); for \(k=2\), (A10) applies directly.
By (A11), \(\delta c_I=\tau_{k+1}(a_I,b_I)\). Consequently the
argument of \(\mathcal K_{k+1}\) is always a full defining system.
The integrand is an integral cocycle, so \(G_k\) is an integral cochain.
The four integrands used in (A12) are exactly
\(T_0,T_1,T_2,T_3\) for \(k=2,3,4,5\), respectively.

## 3. The answer in every dimension

For all input cochains define

\[
\boxed{\mathfrak d_k(A,B,C,D)=
\bigl(\delta_sA,\ \delta B+P_k(A),\
      \delta C+f_k(A,B),\ \delta_sD+g_k(A,B,C)\bigr).}
\tag{A13}
\]

Here are the complete choices of \(f_k\) and \(g_k\):

| k | \(f_k\) on \(L_k\) | \(f_k\) outside \(L_k\) | \(g_k\) on \(L_k\) | \(g_k\) outside \(L_k\) |
| ---: | --- | --- | --- | --- |
| 0 | 0 | no such inputs | 0 | no such inputs |
| 1 | 0 | no such inputs | \(J_1\) | no such inputs |
| 2 | \(\omega B\) | \(H_2(B)\) | \(J_2\) | \(G_2\) |
| 3 | \(\tau'_0(A;B)\) | \(H_3(A,B)\) | \(J_3\) | \(G_3\) |
| 4 | \(\tau'_1(A;B)\) | \(H_4(A,B)\) | \(J_4\) | \(G_4\) |
| 5 | \(\tau'_2(A;B)\) | \(H_5(A,B)\) | \(J_5\) | \(G_5\) |
| 6 | \(\tau'_3(A;B)\) | \(H_6(A,B)\) | \(J_6\) | \(T_3(\mathcal R_6(A,B,C))\) |

In degree two, the two entries for \(f_2\) combine into the single
formula \(f_2=H_2\), since \(\delta B=0\) on \(L_2\).
In degrees zero and one, \(A,B\) are absent. In particular

\[
\begin{aligned}
\mathfrak d_0(0,0,0,D)&=(0,0,0,\delta_sD),\\
\mathfrak d_1(0,0,C,D)&=
\left(0,0,\delta C,
\delta_sD+\frac{\delta_s\widetilde{\omega C}
                       -\widetilde{\omega\delta C}}2\right).
\end{aligned}
\tag{A14}
\]

The optional degree \(k=-1\) also has only \(D\), and uses
\(\delta_sD\). In every row of the table, \(f_k\) depends only on
\(A,B,s,\omega\), never on \(C\).

Here \(\mathcal R_6\) is the section retraction onto full legal triples
specified in [g6_repair.md](g6_repair.md). In the natural auxiliary lower
coordinates, put \(\phi(A,B,C)=(A,B,C+e(A,B))\), with \(e=0\) on
\(L_6\) and \(e=q_{\mathrm{loc}}\) elsewhere. Let \(\sigma\) be the
explicit pointed section of \(\widehat F\) onto its image, and let
\(\backslash\) denote triangular left division for its natural product.
The current choice is
\[
\mathcal R_6(u)=\sigma(\widehat F\phi u)\backslash\phi u.
\]
Here \(F_6\) is the first three components of (A13), and
\(\widehat F\phi=F_6\). The auxiliary natural triangular product
\(\widehat\mu\) and its explicit corrections \(\alpha,\widehat\beta\)
are specified in the
[section construction](../../fermionAHSS_stacking/G6_GENERAL_SECTION.md).
For \(v\backslash w=z\), left division means
\[
\begin{aligned}
A_z&=A_w-A_v,\\
B_z&=B_w+B_v+\alpha(A_v,A_z),\\
C_z&=C_w+C_v+\widehat\beta(A_v,B_v;A_z,B_z).
\end{aligned}
\]
The last two lines are binary. Thus no implicit upper equation is being
solved in the definition of \(\mathcal R_6\).

For a finite free cochain model, construct the section as follows. Given
an output \((a,b,c)\), solve \(\delta_s A_0=a\) in a fixed integral
complement to the cycle lattice. If \(Z\) is the matrix of an integral
cycle basis, enumerate \(A=A_0+Zt\) with \(t\in\{0,1,2,3\}^r\).
For each candidate solve
\[
\delta B=b+P_6(A),\qquad \delta C=c+f_6(A,B),
\]
enumerating the binary kernel of \(\delta\) in the B equation before
trying the C equation. Choose the first solution, with zero chosen at
zero output, to obtain \(\sigma_{\rm current}\), and set
\(\sigma=\phi\sigma_{\rm current}\). The search is complete on
\(\operatorname{im}F_6\) because the original lower differential is
four-periodic under integral cycle shifts; no periodicity of
\(q_{\mathrm{loc}}\) is assumed. Only lower defining equations are solved.
For general X the same construction uses a chosen pointed set section
onto the image; the finite basis is required only for this search algorithm.

The lower homomorphism and triangular cancellation imply
\(\widehat F\mathcal R_6=0\). On this legal locus the natural and current
coordinates agree. Hence its output
\((a_*,b_*,c_*)\) obeys
\(\delta_s a_*=0\), \(\delta b_*=P_6(a_*)\), and
\(\delta c_*=\tau'_3(a_*;b_*)\). Thus the new value is an integral
cocycle. This replaces the earlier zero outside \(L_6\), without using
or defining a degree-four-input operation \(T_4\). The retraction fixes
every full legal triple; every value on \(L_6\) remains \(J_6\).
In particular \(T_3\) here is the entire fixed integral representative
\(\delta_s\Omega_6\), including the coefficient-two prime-three term.
The earlier projection-based repair is historical; this section
retraction is the one used by the current stacking implementation.

To make the endpoint assertion \(\mathfrak d_7\mathfrak d_6=0\)
unambiguous, use (A13) at \(k=7\), with

\[
\begin{aligned}
f_7(A,B)&=
\begin{cases}\tau'_4(A;B),&(A,B)\in L_7,\\
H_7(A,B),&(A,B)\notin L_7,
\end{cases}\\
g_7(A,B,C)&=
\begin{cases}
\displaystyle\frac{\delta_s\widetilde{E_6(C)}
                         -\widetilde{E_7(\delta C)}}2,&A=B=0,\\[4pt]
0,&(A,B)\ne(0,0).
\end{cases}
\end{aligned}
\tag{A15}
\]

The first case for \(g_7\) is integral by (A3). This degree-seven map
is only an endpoint convention for the finite construction. No ko
classification in degree seven, or \(\mathfrak d_8\mathfrak d_7=0\),
is asserted.

**Historical stacking audit and repair (2026-09-25).** The earlier choice
\(g_6=0\) outside \(L_6\) prevents this particular square-zero extension
from admitting a normalized unital triangular stacking product with
\(\mathfrak d(x\times y)=\mathfrak d(x)\times\mathfrak d(y)\) on all
cochains. On \(B((\mathbf Z/3)^2)\), the defining system
\((4\beta_{\mathbf Z,3}(u_1u_2),0,0,0)\) has a nonzero \(T_3\) class;
stacking it with an input outside \(L_6\) would force that class to be a
coboundary. See the [explicit proof](../../fermionAHSS_stacking/ALL_COCHAIN_OBSTRUCTION.md).
That cap has now been replaced by the displayed \(T_3\mathcal R_6\)
formula. On the counterexample its off-shell class is the same nonzero
prime-three obstruction, so the former contradiction no longer applies.
The current section branch is square-zero and admits the explicit
[degree-six stacking construction](../../fermionAHSS_stacking/G6_GENERAL_SECTION.md).
Two legal cylinders prove its strict upper identity without assuming
associativity of the lower cochain product. The older linear-projection
repair is retained as a historical construction, not as the retraction
used by this stacking law.

## 4. Explicit verification of the second iterate

Let \((A',B',C',D')=\mathfrak d_k(A,B,C,D)\).

**First and second components.** For every input,

\[
A''=\delta_s^2A=0,\qquad
B''=\delta P_k(A)+P_{k+1}(\delta_sA)=0.
\tag{A16}
\]

Thus \((A',B')\in L_{k+1}\).

**Third component.** If \((A,B)\in L_k\), then
\((A',B')=(0,0)\), \(f_k=\tau_k\), and \(\delta f_k=0\).
Otherwise \(f_k=H_k\), and (A11) gives
\(\delta f_k=\tau_{k+1}(A',B')=f_{k+1}(A',B')\).
In both cases,

\[
C''=\delta f_k(A,B)+f_{k+1}(A',B')=0.
\tag{A17}
\]

In particular \((A',B',C')\) is a full defining system.

**Fourth component, input in \(L_k\).** Here \(C'=t_k\) is closed
and \(D'=\delta_sD+J_k\). The next correction, including the
degree-seven endpoint, is
\(g_{k+1}(0,0,t_k)=\delta_s\mathsf h(E_k(t_k))\).
Consequently (A6) gives the literal integral equality

\[
D''=\delta_sJ_k+\delta_s\mathsf h(E_k(t_k))=0.
\tag{A18}
\]

**Fourth component, input outside \(L_k\), \(2\le k\le5\).**
The raw \(H_k\), unlike the piecewise \(f_k\), commutes with
restriction to either end of the interval. The endpoints of the full
defining system in (A12) are therefore zero and
\((A',B',\delta C+H_k(A,B))=(A',B',C')\). Its integral obstruction
is a cocycle. Signed prism Stokes gives

\[
\delta_sG_k=-\mathcal K_{k+1}(A',B',C').
\tag{A19}
\]

The next input is legal and has zero third residual by (A17), so
\(g_{k+1}(A',B',C')=J_{k+1}(A',B',C')
=\mathcal K_{k+1}(A',B',C')\). Hence

\[
D''=-\mathcal K_{k+1}(A',B',C')
       +\mathcal K_{k+1}(A',B',C')=0.
\tag{A20}
\]

**Fourth component, input outside \(L_6\).** The repaired value
\(g_6=T_3(\mathcal R_6(A,B,C))\) is an integral cocycle, and
\((A',B')=p_6(A,B)\ne(0,0)\). Formula (A15) therefore gives
\(g_7(A',B',C')=0\), so
\(D''=\delta_sg_6+g_7(A',B',C')=0\).
For an input in degree five, its output first pair lies in \(L_6\),
so the unchanged \(J_6\) branch proves \(\mathfrak d_6\mathfrak d_5=0\)
exactly as before.
Together with the elementary degree-zero case, this proves

\[
\boxed{\mathfrak d_{k+1}\mathfrak d_k=0
\quad\text{on all four cochains, for every }0\le k\le6.}
\tag{A21}
\]

This proof uses cochain identities, not componentwise additivity.
It uses the same fixed helper identities as the repository; it does not
constitute a new proof of those calibrated mathematical inputs.

## 5. Agreement with the requested obstruction problem

The first two zero-output equations in (A13) force \((A,B)\in L_k\).
The third then forces \(\delta C=\tau_k(A,B)\). On precisely those
inputs the fourth equation is

\[
\delta_sD=-\mathcal K_k(A,B,C).
\tag{A22}
\]

Thus all zero-output equations through degree six are exactly the
repository's defining-system equations. The values chosen outside
\(L_k\), including the revised degree-six branch, cannot create an extra
solution: the first or second component is already nonzero there.
For a fixed legal \(A,B,C\), failure to solve (A22) is exactly the
nonvanishing of the integral class \([\mathcal K_k(A,B,C)]\).

For \(3\le k\le6\), put \(n=k-3\). The actual AHSS \(T=d_5\) obstruction accounts for
all allowed choices of \(B,C\). Its target is the documented quotient

\[
\frac{H^{n+5}(X;\mathbf Z_s)}
{\operatorname{Dtilde}H^{n+2}(X;\mathbf F_2)
 +\operatorname{Psi}_{n+1}
   \ker(D:H^{n+1}(X;\mathbf F_2)\to H^{n+3}(X;\mathbf F_2))},
\tag{A23}
\]

where the second image is taken after the first quotient. A nonzero
ordinary class for one defining system need not remain nonzero after
allowing the other defining systems. This is the same qualification
used by [tertiary_operations.md](tertiary_operations.md), Section 1.
The earlier rows retain \(\operatorname{Dbar},D,\operatorname{Dtilde},
\operatorname{Tau},\operatorname{Psi}\) in the same way.

The comparison with Wang–Gu concerns their legal decoration data and
obstructions. Their signed first cocycle equation is explicit in
[Eq. (250)](https://arxiv.org/html/1811.00536v3#S6.SS6.SSS1).
This construction leaves the existing calibrated low-dimensional
formulas on such data intact. The global piecewise extension and its
endpoint convention are constructions here, not formulas attributed
to Wang–Gu.

The compatible stacking law is defined in the separate
[degree-indexed stacking note](../../fermionAHSS_stacking/ALL_COCHAIN_STACKING.md).
Its degree-six formula uses this same section retraction and two legal
cylinders; it does not assume associativity of the lower product.
The exact classification agreement established here is the successive
obstruction/solvability test requested in the question, with the
repository's page indeterminacies. It does not upgrade the repository's
five-row associated graded to a complete ko classification.

### The Z4 four-factor obstruction

Take \(G=(\mathbf Z/4)^4\), physical degree \(k=3\), \(s=0\), and
\([\omega]=x_1x_2+x_3x_4\). For coordinates in \(\{0,1,2,3\}\), put
\[
F(g,h)=g_1h_2+g_3h_4,\qquad w=F\bmod2\in\{0,1\},
\qquad \omega=\rho_2w.
\]
Because F is bilinear modulo four, \(\delta F\in4\mathbf Z\). The
explicit lower solution is
\[
\boxed{A_3=2,\quad B_3=0,\quad
 C_3=C_0=\rho_2\frac{F-w}{2},\qquad
 \delta C_0=\rho_2(\delta w/2)=\tau'_0(2;0).}
\]
The lower obstruction is not the zero cochain: on
\([e_1|e_1|e_2]\), both \(\delta C_0\) and \(\tau'_0(2;0)\)
equal one. This is a full legal defining system, so neither the
arbitrary-input prism branch nor the repaired g6 enters this example.

Define the integral Pontryagin cochain
\(P_\omega=w\smile w+w\smile_1\delta w\). The actual production
phase and integral obstruction reduce to
\[
\boxed{\Omega_3(C)=\tfrac12\widetilde{E_2(C)}-\tfrac18P_\omega,
\qquad g_3(2,0,C)=\delta\Omega_3(C)}
\quad\bigl(\delta C=\tau'_0(2;0)\bigr),
\]
where \(E_2(C)=C\smile C+C\smile_1\delta C+\omega\smile C\) is
reduced as a whole modulo two before lifting. This is R0 at A=2 in the
repository's convention. The closed 24-term normalized bar cycle
\[
z=\sum_{\sigma\in S_4}\operatorname{sgn}(\sigma)
[e_{\sigma(1)}|e_{\sigma(2)}|e_{\sigma(3)}|e_{\sigma(4)}]
\]
has zero boundary by adjacent-transposition and endpoint cancellation,
and satisfies
\[
\langle\widetilde{E_2(C_0)},z\rangle=0,\qquad
\langle P_\omega,z\rangle=2,\qquad
\boxed{\langle\Omega_3(C_0),z\rangle=-\tfrac14.}
\]
Its phase is \(-i\). The supplied double-pip PDF uses the opposite
quadratic-phase sign and gives \(+i\) for its canonical choice; both
are nontrivial. Since G is finite, positive-degree rational cohomology
vanishes and the connecting map
\[
H^4(BG;\mathbf Q/\mathbf Z)\xrightarrow{\ \beta\ }H^5(BG;\mathbf Z)
\]
is an isomorphism. Thus \([g_3]=\beta[\Omega_3]\ne0\), and no
integral D3 solves \(\delta D_3+g_3=0\). The four-cycle detects the
rational phase, not the integral five-cocycle directly.

This conclusion survives all lower choices. Indeed
\[
\tau'_0(2;B)=\omega\smile B+\rho_2(-\delta w/2),\qquad
H^*(BG;\mathbf F_2)=\Lambda(x_1,x_2,x_3,x_4)
 \otimes\mathbf F_2[y_1,y_2,y_3,y_4],
\]
with \(|x_i|=1\), \(|y_i|=2\). Multiplication by \(\omega\) sends
the four x_i to four independent exterior triples, so solving for C
forces \([B]=0\). In the normalized group-bar complex, degree-one
coboundaries with trivial coefficients are zero; hence B=0 literally,
and the Psi indeterminacy vanishes. Every remaining C is \(C_0+b\)
for a binary two-cocycle b. For every such cochain, including exact
shifts,
\[
\langle\Omega_3(C_0+b),z\rangle\in-\tfrac14+\tfrac12\mathbf Z.
\]
It can never be integral, so the Dtilde indeterminacy cannot remove the
obstruction either. Consequently
\[
\boxed{T_0(2)=d_5(2)\ne0\quad\text{in the actual AHSS page quotient}.}
\]
Here k=3 is the physical degree, d5 is the page differential, and O4
is the name for its degree-four phase obstruction in the supplied PDF.
The mathematical conclusion of the separately archived double-pip note
is therefore part of this note, rather than dependent on that archive.

## 6. Verification record

The identities (A16)–(A21) were checked algebraically in all branches,
including the independent second interval in (A12), signed prism Stokes,
the global branch test, and the degree-six endpoint. Two independent
mathematical audits agreed with this calculation. The following ledger
retains the necessary information from the dated exact-arithmetic checks:

| Recorded checks | Completed results and scope |
| --- | --- |
| Legal first two layers, low degrees; seed 813 | 100 inputs each in k=1,2; all integrality and literal second-iterate checks passed. |
| Legal first two layers, arbitrary C; base seed 20260925 | 12 cases in k=3,4,5,6; 96 integral boundary-face evaluations and 12 exact zero second iterates. Six mixed cases had nonzero B, Tau and the residual cup correction. Signed A of degree zero included amplitudes 2,3,6. |
| Arbitrary first two layers, original prism construction | Four full checks in k=2,3,4,5: 200 integral cylinder T evaluations, 31 nonzero; 30 integral g face values; all four second iterates zero. Next T values 1 and -2 exercised nontrivial cancellation. The independent interval coordinates were retained. |
| Dense arbitrary A,B,C,D, or closed A with incorrect B; seed 925260 | Four full checks in k=2,3, two in each degree, with 144 integral cylinder T evaluations, 70 nonzero. Four further k=4,5 cases checked the lower three components and delta_s²D only; their full fourth components were uncomputed. All ten recorded inputs, including the historical k=6 pair below, had nonclosed D. |
| Historical degree-six zero cap | One sparse and two dense k=6 checks used the former g6=0 outside L6 and its endpoint. The dense report therefore records six full passes in total, of which two concern that superseded branch. The sparse report's ten zero g face values do not test the current section repair. |

The high-degree samples used the actual production R0–R3 assemblers.
Their sparse A inputs make the cubic and prime-three phases vanish, so
these samples do not independently exercise those calibration terms.
For k=3,4,5 the next production assembler was evaluated at (0,0,t);
for k=6 only the stated pure-C endpoint was used, with no R4 or T4.
Three recorded dense k=4 top-evaluation attempts timed out or were
interrupted after a preceding timeout. They are neither successful full
checks nor counterexamples. The four k=4,5 lower-only checks and the
earlier sparse full checks remain separate evidence.

The current section repair has the symbolic proof in Section 4 and the
separate finite-model verification in [g6_repair.md](g6_repair.md).
Its production adapter passed a normalized BC2 check with nonzero omega
and nonzero A outside L6: gamma6=0, g6(u times v)=1 and gamma7(Fu,Fv)=1,
giving the exact upper identity 1=1. The completed 35.750-second check and
the earlier incomplete 240-second attempt are distinguished in
[the stacking verification record](../../fermionAHSS_stacking/VERIFICATION.md).
The historical zero-cap reports must not be cited as tests of this repair.

The double-pip example in Section 5 additionally has an exact 24-term
cycle check, all 1024 H2 class shifts (512 residues 1/4 and 512 residues
3/4), 128 sampled lower equations and 32 sampled integral top values.
The analytic proof there covers every lower choice, including exact C
shifts. This is a production-formula check, not a full GAP page run for
the 256-element group.

The original scripts, full inputs, rational outputs and completion
statuses are preserved outside the repository in the
[dated evidence archive](/Users/victor/Documents/miscellaneous/fermionAHSS/20260925-nonlinear-square-zero/README.md).
The [archive index](/Users/victor/Documents/miscellaneous/fermionAHSS/README.md)
gives current paths and reproduction instructions; the original records
retain their historical path names and are unchanged. Finite checks do
not replace the cochain proofs or establish naturality. No GAP runtime
implementation or calibration payload was changed.

## 7. Rechecking completely arbitrary A, B, C and D

The input restrictions in the secondary and tertiary evaluators apply to
their **constructed arguments**, not to the four initial cochains. In
particular, arbitrary integer \(A\) is allowed, with no parity or signed
cocycle condition. Arbitrary binary \(B,C\) and arbitrary integer \(D\)
are allowed as well. The background \(s,\omega\) remains a pair of
cocycles, as required to define the twist.

Here is a compact audit of the entire composition. Write

\[
(a,b,c)=F_k(A,B,C)=
\bigl(\delta_sA,\delta B+P_k(A),\delta C+f_k(A,B)\bigr).
\]

Expanding the second application before imposing any input conditions
gives

\[
\boxed{\mathfrak d_{k+1}\mathfrak d_k(A,B,C,D)=
\left(0,\
\delta P_k(A)+P_{k+1}(\delta_sA),\
\delta f_k(A,B)+f_{k+1}(a,b),\
\delta_sg_k(A,B,C)+g_{k+1}(a,b,c)\right).}
\tag{A24}
\]

The disappearance of \(C\) from the third component uses
\(\delta^2C=0\). The disappearance of \(D\) from the fourth component
uses \(\delta_s^2D=0\); it does **not** assume \(D=0\) or
\(\delta_sD=0\). Equations (A16)–(A20), with the stated endpoint,
prove that all three remaining expressions vanish for every input.

The two possible first-two-layer outcomes exhaust all cases. If
\(p_k(A,B)=0\), only the initial \(A,B\) satisfy the first two
equations and \(C,D\) remain arbitrary. If \(p_k(A,B)\ne0\), the
original \(A\) may itself be nonclosed, or it may be closed while \(B\)
fails its equation. Both possibilities use the same branch proof:
the prism construction for \(k\le5\) and the section retraction for
\(k=6\). The equations obeyed by \((a,b,c)\) are conclusions of the
calculation, not additional hypotheses on \((A,B,C)\).

### A simple telescoping form of the correction

There is a useful general identity. Suppose \(F_{k+1}F_k=0\), and let
\(\Phi_k(x)\in C^{k+1}(X;\mathbf Q_s)\) be pointed phase lifts,
\(\Phi_k(0)=0\), where \(x=(A,B,C)\). Set

\[
\boxed{g_k(x)=\delta_s\Phi_k(x)-\Phi_{k+1}(F_kx).}
\tag{A25}
\]

Then the fourth component cancels by the exact rational identity

\[
\begin{aligned}
\delta_sg_k(x)+g_{k+1}(F_kx)
 &=-\delta_s\Phi_{k+1}(F_kx)
   +\delta_s\Phi_{k+1}(F_kx)
   -\Phi_{k+2}(F_{k+1}F_kx)\\
 &=0.
\end{aligned}
\tag{A26}
\]

For this to define the requested integer fourth component, one must also
have

\[
\delta_s\Phi_k(x)\equiv\Phi_{k+1}(F_kx)
\pmod{C^{k+2}(X;\mathbf Z_s)}.
\tag{A27}
\]

Thus subtracting the phase evaluated on the lower output is the simple
square-zero correction. Arbitrary rational phase formulas satisfy the
telescoping calculation but need not satisfy (A27).

The rational lift is also important for retaining the obstruction: on
a full defining system \(F_kx=0\), (A25) becomes
\(g_k(x)=\delta_s\Phi_k(x)\). This integral cocycle can represent a
nonzero integral cohomology class when \(\Phi_k\) is rational. Replacing
\(\Phi_k\) by an integer cochain would make that class a coboundary
and eliminate the required \(T\) obstruction.

The construction in this note provides the required integral lifts.
Normalize \(\Omega_k\) to have zero value on zero input, if necessary,
by subtracting the closed cochain \(\Omega_k(0,0,0)\). This leaves
\(\delta_s\Omega_k\), \(J_k\), and \(G_k\) unchanged. For
\(0\le k\le6\), take

\[
\Phi_k(A,B,C)=
\begin{cases}
\Omega_k+\mathsf h(t_k\smile_{k-1}\tau_k),&(A,B)\in L_k,\\
I\Omega_{k+1}(a_I,b_I,c_I),&(A,B)\notin L_k,\ 2\le k\le5,\\
\Omega_6\bigl(\mathcal R_6(A,B,C)\bigr),&k=6,\ (A,B)\notin L_6.
\end{cases}
\tag{A28}
\]

At \(k=0\) the first line means \(\Phi_0=0\). For the endpoint use
\(\Phi_7(A,B,C)=\mathsf h(E_6C)\) when \(A=B=0\), and zero
otherwise. Formula (A25) reproduces exactly every \(g_k\) in the
degree table for \(0\le k\le6\): on \(L_k\) it is (A5), and
outside \(L_k\) for \(2\le k\le5\), signed Stokes gives
\(\delta_sI\Omega_{k+1}-\Omega_{k+1}(F_kx)
=-I\delta_s\Omega_{k+1}=G_k\). Integrality follows respectively
from (A5) and from the integral cocycle inside (A12). At the cutoff,
an input outside \(L_6\) has nonzero first-two-layer output, so
\(\Phi_7(F_6u)=0\), while
\[
\delta_s\Phi_6(u)
=\delta_s\Omega_6(\mathcal R_6u)
=T_3(\mathcal R_6u).
\]
This is integral and closed because \(\mathcal R_6u\) is a full legal
defining system. Thus (A25) reproduces the repaired, generally nonzero
\(g_6\) branch as well. On \(\operatorname{im}F_6\),
the already specified \(g_7\) agrees with the needed next correction.
For the endpoint use of (A26), set \(\Phi_8(0)=0\); no other value of
\(\Phi_8\) is used or asserted. There is still no assertion about a ko
operation in degree seven.

This reformulation does not remove the global branch choice. It makes
explicit why unrestricted \(A,B,C,D\) do not require further changes
to the stated piecewise extension, and why integer lifts matter.

Additional bounded checks with independently populated integer \(A,D\)
and binary \(B,C\), including nonclosed \(D\) and closed \(A\) with
an incorrect defining \(B\), are recorded in the
[archived arbitrary-input records](/Users/victor/Documents/miscellaneous/fermionAHSS/20260925-nonlinear-square-zero/README.md#arbitrary-input-recheck).
Section 6 retains the completion counts, uncomputed components and
superseded degree-six scope. The proof above is independent of their sampling.
