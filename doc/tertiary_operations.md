# The production tertiary operation: \(d_5\) in input degrees 0–3

This is a transcription of the formulas evaluated by the current GAP
`T` callback and its bundled Python kernel, as of 2026-09-23. The
cochain notation and signs are fixed in [conventions.md](conventions.md);
the lower operation \(\psi_n\) is GAP `Tau`, described in
[secondary_operations.md](secondary_operations.md). The universal source
operators below are specified in [universal_helpers.md](universal_helpers.md).

The production entry points are
[`koAHSSNaturalTertiary`](../gap/natural_tertiary.gi) and
[`koAHSSNaturalTCallback`](../gap/natural_tertiary.gi).
Degrees outside \(0\leq n\leq3\) are rejected. The old local-R chooser
is not part of this formula.

## 1. Defining system and returned integral differential

Fix binary cocycles \(s\in Z^1(X;\mathbf F_2)\),
\(\omega\in Z^2(X;\mathbf F_2)\) and a signed integral cocycle
\(A\in Z^n(X;\mathbf Z_s)\). Let \(\rho\) denote reduction modulo
two and let \(\widetilde z\) be the standard binary integral lift.
Choose the defining cochains

\[
a=\rho A,\qquad b\in C^{n+1}(X;\mathbf F_2),\qquad
c\in C^{n+2}(X;\mathbf F_2),
\qquad db=Da,\quad dc=\psi_n(A,b),
\]

where the literal primary representative is
\(Da=Sq^2a+\omega a+s\rho(d\widetilde a/2)\), representing
\(Sq^2a+sSq^1a+\omega a\) in cohomology. The same calibrated
lower representative \(\psi_n\) is used here and in the page callback:
`chi7_tail`, secondary epsilon \((1,0,0)\), and eta \((1,0,1)\).

For a binary degree-\(r\) cochain \(z\), including noncocycles, write

\[
Q^jz=z\cup_{r-j}z+z\cup_{r-j+1}dz,\qquad
Ez=Q^2z+\omega z,\qquad Q_Dz=Ez+sQ^1z.
\]

All expressions inside these three operations are binary. Define

\[
\mathsf h(z)=\tfrac12\widetilde z,\qquad
P_\omega=\widetilde\omega\smile\widetilde\omega
 +\widetilde\omega\cup_1d\widetilde\omega.
\]

The argument of \(\mathsf h\) is always reduced modulo two **as a
whole** before lifting. It is useful to distinguish the actual rational
phase \(\mathcal O_n\) returned by the kernel from the fully boundary-normalized
helper \(R_n\) in the notes:

\[
\boxed{\quad T_n(A;b,c)=d_s\mathcal O_n(A,b,c),\qquad
\mathcal O_n=\mathsf h(Ec)+\widehat R_n(A,b)
 +\mathbf1_{n=3}\tfrac23\widetilde{P^1_s\rho_3 A}.\quad}
\tag{T}
\]

The source-primitive and residual identities are identities of phases
in \(\mathbf Q_s/\mathbf Z_s\). The formulas in this document specify
their actual rational lifts, so that taking the further differential
also specifies an integral cocycle representative.

Here \(d_s=d-2\widetilde s\smile-\); the final lift in (T) has values
\(0,1,2\), not \(0,1\). In degrees \(1,2,3\), the kernel does not
evaluate the final signed-exact boundary-normalization term
\(d_s\mathcal D_n\) displayed in the notes. Its further signed
differential is exactly zero, so omitting this term has no effect on
\(T_n\) for the same choice of lifts. The lifts of the remaining
expressions are specified below; the notes often state equalities only
modulo integers. In degree zero the full \(R_0\) below is used.

By contrast, changing a rational phase by an **integral cochain** changes
\(T_n\) by an integral coboundary. This preserves its cohomology class but
not necessarily the literal GAP output. The separate half-lifts and
integer coefficients below retain the kernel's representative convention.

The calibrated ko correction relative to this Danus reference has rank
coefficient zero and \(\mu_R=0\). The coefficient-two three-primary
term in (T) remains present. The legacy `TReference` coefficient-one
correction must not be applied to this final `T`.

On the page, the domain consists of classes in \(\ker\overline D\)
whose \(\operatorname{Tau}_n\) class vanishes. The target is

\[
\frac{H^{n+5}(X;\mathbf Z_s)}
 {JH^{n+2}(X;\mathbf F_2)
 +\operatorname{Psi}_{n+1}\!\left(\ker\bigl(D:H^{n+1}(X;\mathbf F_2)
                 \longrightarrow H^{n+3}(X;\mathbf F_2)\bigr)\right)}.
\]

The second image is understood in the quotient by the first image.
The defining-system solver adjusts \(b\) by a closed cochain when needed
to make \(\psi_n(A,b)\) exact before solving for \(c\). Changes of the
allowed defining system are handled by this target quotient.

## 2. Shared boundary phase, prism, and A-only splitting

### The fixed boundary phase

For a closed binary degree-\(m\) cochain \(q\), put

\[
B_q=\frac{d\widetilde q}{2},\quad e_q=\rho B_q,\quad
C_q=\frac{B_q+\widetilde e_q}{2},\quad
u_q=Sq^2q,\quad v_q=\omega q,\quad w_q=se_q,
\]

\[
\begin{aligned}
H_m(q)={}&\zeta_{2,m}(\omega,q)+\chi_m(q)
+u_q\cup_{m+1}v_q+u_q\cup_{m+1}w_q+v_q\cup_{m+1}w_q\\
&+\zeta_{1,m+1}(s,e_q)+(\omega\cup_1s)e_q
+su_q+s^2\rho C_q,\\
R_{\partial,m}(q)={}&\mathsf h(H_m(q))
+\tfrac14\bigl(\widetilde\omega B_q+B_q\cup_{m-1}B_q\bigr)
+\mathsf h\bigl(sSq^2q+\omega Sq^1q\bigr).
\end{aligned}
\]

The helper \(\Theta_m(q,v)\) used by the tertiary source and prism
formulas is the **phase**, not GAP's integral `Psi` output:

\[
\boxed{\Theta_m(q,v)=\mathsf h\bigl(Ev+H_m(q)\bigr)
+\tfrac14\bigl(\widetilde\omega B_q+B_q\cup_{m-1}B_q\bigr)
+\mathsf h\bigl(sSq^2q+\omega Sq^1q\bigr).}
\tag{B}
\]

The binary \(Ev+H_m(q)\) is lifted together, as in
[`phase_eval.theta`](../python/phase_eval.py). In particular,
(B) agrees with \(\mathsf h(Ev)+R_{\partial,m}(q)\) modulo integral
cochains; the latter split is useful mathematically but is not always the
same real lift. Its source input satisfies \(dv=Dq\).
The final term in (B) also lifts its binary sum **once**, whereas GAP's
secondary formula (S7) uses an integer sum of the two separate binary
lifts. This can change the rational representative by an integral cochain.
For a cocycle, \(Sq^1q=e_q\) holds literally in the chosen convention;
the difference here is in the lifting of the sum.

### The raw right prism

Let \(\ell\) be the interval vertex coordinate \(0,1\), and pull
background cochains back to \(X\times\Delta^1\). For a degree-\(r\)
cochain \(z\), define

\[
(Iz)(v_0,\ldots,v_{r-1})=
\sum_{j=0}^{r-1}(-1)^j
z((v_0,0),\ldots,(v_j,0),(v_j,1),\ldots,(v_{r-1},1)).
\tag{I}
\]

For \(b_I=b\ell\), the interval factor is evaluated at the last vertex.
Every binary lift and integral carry in \(\Theta\) is recomputed on
each prism simplex. A prism of a rational phase is the exact signed sum
in (I), with no reduction modulo one.

### Source data for n=1,2,3

The following formulas are exactly
[`phase_eval.source`](../python/phase_eval.py). All letters
except the explicitly integral \(B,C_B,L_A,q_{\rm int}\) and \(A\) are
binary:

\[
\begin{gathered}
a=\rho A,\qquad t=\rho\frac{A-\widetilde a}{2},\qquad
B=\frac{d\widetilde a}{2},\qquad e=\rho B,\qquad
C_B=\frac{B+\widetilde e}{2},\\
u=Sq^2a,\quad v=\omega a,\quad w=se,\quad p=u+v+w=Da,
\end{gathered}
\]

\[
\begin{aligned}
F_A={}&\zeta_{2,n}(\omega,a)+\chi_n(a)
+u\cup_{n+1}v+u\cup_{n+1}w+v\cup_{n+1}w\\
&+\zeta_{1,n+1}(s,e)+(\omega\cup_1s)e+su+s^2\rho C_B,\\
g={}&Q^2t+\omega t+(dt)\cup_n(sa)
+\zeta_{1,n}(s,a)+(\omega\cup_1s)a,\\
q_{\rm int}={}&\widetilde\omega B+B\cup_{n-1}B,\\
L_A={}&\frac{q_{\rm int}-d_s\widetilde g
                  +\widetilde s\smile\widetilde p}{2},\\
k_0={}&F_A+\rho L_A+s^3a.
\end{aligned}
\tag{S}
\]

The plus carry \(C_B=(B+\widetilde e)/2\) must not be changed to a minus
carry. Both divisions in (S) are exact integer divisions, including on
negative inputs. For \(n=2,3\), set \(k=k_0\), and define

\[
y=Q_Db,\quad
\lambda=\rho\frac{\widetilde g+\widetilde s\smile\widetilde b
                         -\widetilde{g+sb}}2,
\quad \psi_n=y+k+d\lambda.
\tag{S23}
\]

Then \(dk=Dp\). Both \(p,k\) depend only on \((A,s,\omega)\); the carry
\(\lambda\) also depends on \(b\). Degree one has the additional
gauge specified in section 4.

## 3. Input degree zero

Here \(A\in Z^0(X;\mathbf Z_s)\), \(b\in C^1\), \(c\in C^2\), and
\(\mathcal O_0=\mathsf h(Ec)+R_0(A,b)\). The branch choice is
componentwise, because \(\rho A\) is locally constant.

### Even A: arbitrary sign representative

Set

\[
\begin{gathered}
K=A/2,\qquad t=\rho K,\qquad q=(K-\widetilde t)/2,\qquad
Z=t\omega,\\
N_b=\frac{\widetilde Z+\widetilde{sb}-\widetilde{Z+sb}}2,
\quad\lambda=\rho(q\smile\widetilde\omega+N_b),\\
H=Q_Db,\qquad r=\rho(d\widetilde\omega/2),\qquad z=tr.
\end{gathered}
\]

Since \(A\) is even, \(db=0\), so \(R_{\partial,1}(b)=\Theta_1(b,0)\)
is defined. The literal formula is

\[
\boxed{R_0(A,b)=R_{\partial,1}(b)-\frac{K\smile P_\omega}{8}
+\mathsf h\bigl(H\cup_2z+E\lambda+(H+z)\cup_2d\lambda\bigr).}
\tag{R0e}
\]

The signed integer \(K\) and the separate quarter-input carry \(q\) are
both retained. Reducing \(K\) to its parity loses part of the formula.
At \(A=0\), all additional terms vanish literally.

### Odd A: canonical sign gauge and its prism

Let \(N=|A|\), \(\epsilon_A(v)=\mathbf1_{A(v)<0}\), and
\(t=\rho((N-1)/2)\). In the positive oriented frame set

\[
\lambda_+=b^2+tQ^1b,\qquad
R_+(N,b)=-\mathsf h(E\lambda_+)
+\frac{3+2\widetilde t}{8}\,d(\widetilde b^{\,3}).
\]

The integer cube and its differential in this expression use the
ordinary ordered integral cup product. Define on the cylinder

\[
\gamma(v,e)=e\epsilon_A(v),\quad
s_I=s+d\gamma,\quad A_I=(-1)^\gamma A,
\quad b_I=\pi^*b,\quad\omega_I=\pi^*\omega,
\]

with \(s_I\) reduced modulo two. Evaluate the **current lower formula**
\(\psi_I=\psi_0^{s_I,\omega}(A_I,b_I)\), and put

\[
\boxed{R_0(A,b)=(-1)^{\epsilon_A(v_0)}R_+(N,b)
            +\mathsf h\bigl(\rho I(E\psi_I)\bigr).}
\tag{R0o}
\]

The parity of the signed prism sum is taken before its binary half-lift;
there are five prism terms on a four-simplex. The sign in the first
summand refers to its first vertex \(v_0\). Since \(d_sA=0\), the odd
branch has \(s=d\epsilon_A\), so the top endpoint is the positive
oriented input. This is a canonical gauge determined by \(A\), not a
chosen primitive of \(s\).

Source: [`low_phases.R0`](../python/low_phases.py), with
derivation in R0_fixed.md (source-workspace provenance: `notes/extra/tertiary_R_Danus/R0_fixed.md`; not bundled).

## 4. Input degree one

Use (S) with \(n=1\), but write its original primary cochain as \(p_0=p\).
Define

\[
\begin{gathered}
e_0=st,\qquad b^+=b+e_0,\qquad q=(\omega+s^2)a=db^+,\\
k=k_0+Q_De_0+(de_0)\cup_2p_0+s\bigl((de_0)\cup_3p_0\bigr),
\qquad y=Q_Db^+,\\
N=\rho\frac{\widetilde g+\widetilde{sb}-\widetilde{g+sb}}2,\\
\lambda=N+b\cup_1e_0+p_0\cup_2e_0+p_0\cup_3de_0\\
\hphantom{\lambda={}}+s\bigl(b\cup_2e_0+p_0\cup_3e_0
                              +p_0\cup_4de_0\bigr),\\
\psi_1=y+k+d\lambda,\qquad
\operatorname{Pol}_1(y,k)=y\cup_3k+k\cup_4dy+Q^1k.
\end{gathered}
\]

Let \(V_1(A,s,\omega)\) be the fixed degree-five A-only source primitive
of \(\Theta_3(q,k)\), with the strict background normalization
specified in [universal_helpers.md](universal_helpers.md). Then the
literal production assembler is

\[
\boxed{\begin{aligned}
\widehat R_1={}&-I\Theta_3\bigl(d(b^+\ell),Q_D(b^+\ell)\bigr)-V_1\\
&+\mathsf h\bigl(\operatorname{Pol}_1(y,k)+E\lambda
+\psi_1\cup_3d\lambda+EQ^1b^++s\psi_1\bigr),\\
\mathcal O_1={}&\mathsf h(Ec)+\widehat R_1.
\end{aligned}}
\tag{R1}
\]

The rank-normalization term would be
\(-\varepsilon_R(P_\omega\cup_s A)/4\); its evaluated coefficient is
\(\varepsilon_R=0\). This scalar is distinct from the lower
secondary epsilon vector \((1,0,0)\). The negative prism sign in (R1)
is part of the fixed right-prism convention.

Source: [`low_phases.R1_raw` and `build_phase`](../python/low_phases.py).

## 5. Input degree two

Use (S23), and put

\[
\begin{gathered}
\operatorname{Pol}_2(y,k)=y\cup_4k+k\cup_5dy+Q^1k,\\
U_A=\zeta_{1,4}(s,p)+sk+(\omega\cup_1s)p,\\
\Pi_2=I\Theta_4\bigl(d(b\ell),Q_D(b\ell)\bigr),\qquad
L_2=\tfrac14P_\omega\cup_s A.
\end{gathered}
\]

The evaluated selectors are \(\alpha=1,\ \beta=0\). With the
fixed A-only source \(V_2\) of \(\Theta_4(p,k)\), the actual phase is

\[
\boxed{\begin{aligned}
\widehat R_2={}&\mathsf h(\operatorname{Pol}_2(y,k))+\Pi_2-V_2
+\mathsf h\bigl(E\lambda+(y+k)\cup_4d\lambda\bigr)\\
&-\mathsf h(E(sb))+\mathsf h(U_A)+L_2-\tfrac14 A^{\cup3},\\
\mathcal O_2={}&\mathsf h(Ec)+\widehat R_2.
\end{aligned}}
\tag{R2}
\]

Each \(\mathsf h\) in (R2) is a separate half-lift, exactly as in
[`phase_eval.phase2`](../python/phase_eval.py). Combining
these brackets changes the real phase by an integral cochain in general.
The potential \((\beta/2)\widetilde{s^2\omega a}\) term is zero.
The quarter-cubic term is **present**, and this is the current
\(R_2=R_2^{\rm sharp}\) family.

The local-system cup product \(P_\omega\cup_s A\) transports the signed
right factor from the join vertex to the first vertex. The ordered cubic
\(A^{\cup3}\) has all three sign coefficients multiplied with
their transports. On a six-simplex, its value is

\[
(-1)^{s_{02}+s_{04}}A_{012}A_{234}A_{456}.
\]

The historical source
\(V_{2,\mathrm{fin}}=V_2-L_2\), before the cubic change, is the
one used in the degree-three selector calibration. The current source in
the assembled R2 is
\(V_{2,\mathrm{cur}}=V_2-L_2+A^{\cup3}/4\). They must not be
interchanged in that calibration.

Source: [`low_phases.build_phase`](../python/low_phases.py).

## 6. Input degree three

Use (S23), and define

\[
\begin{gathered}
\operatorname{Pol}_3(y,k)=y\cup_5k+k\cup_6dy+Q^1k,\\
\Pi_3^+=I\Theta_5\bigl(d(b\ell),Q_D(b\ell)\bigr),\\
L_3=\tfrac14P_\omega\cup_s A,\qquad
N_3=\mathsf h(s^2Sq^2a),\qquad
O_3^{\mathrm{corr}}=\mathsf h(\omega Sq^2a),\qquad
M_3=\mathsf h(s^2\omega a).
\end{gathered}
\]

The universal degree-seven source primitive \(V_3\) is evaluated from
the fixed finite source table; it depends on \(A,s,\omega\), not \(b,c\).
The selectors from the current R2 suspension comparison are

\[
(c_4,c_N,c_O,c_M,\epsilon_c)=(1,0,1,1,1),\qquad\xi=3/4.
\]

The literal production formula is therefore

\[
\boxed{\begin{aligned}
\widehat R_3={}&\mathsf h\bigl(\operatorname{Pol}_3(y,k)+E\lambda
                              +(y+k)\cup_5d\lambda\bigr)
                   -\Pi_3^+-V_3\\
&+\mathsf h(EQ^1b+s\psi_3)
  +L_3+N_3+2O_3^{\mathrm{corr}}+M_3,\\
\mathcal O_3={}&\mathsf h(Ec)+\widehat R_3
                         +\tfrac23\widetilde{P^1_s\rho_3 A}.
\end{aligned}}
\tag{R3}
\]

In particular, the sign in front of \(V_3\) is **negative** in the final
assembler. The coefficients in (R3) were substituted as ordinary integers:
\(c_N+\epsilon_c=1\), \(c_O+\epsilon_c=2\). Although
\(2O_3^{\mathrm{corr}}=\widetilde{\omega Sq^2a}\) is an integral
cochain and can be suppressed in a phase class modulo integers, it is
retained by [`high_phase.build_phase`](../python/high_phase.py).
Its removal would generally change the returned integral representative
by a coboundary.

The calibration periods on \((F,Y,W,V,Z+T)\) are
\((3/4,1/4,0,1/2,1/2)\); see the exact
[`high_calibration.json`](../python/high_calibration.json).
The source table uses dyadic denominators up to 16. The final phase can
also contain a factor of three from the following term.

## 7. Explicit degree-three prime-three term

Let \(a_3=\rho_3 A\), and write its standard lift with values \(0,1,2\).
For an ordered face \(F=(i,j,k,l)\) of the seven-simplex, write

\[
a_F^{(0)}=(-1)^{s(v_0,v_i)}a_3(v_i,v_j,v_k,v_l),
\]

where the sign is \(1\) for \(i=0\). The actual reduced power is

\[
P^1_s(a_3)(v_0,\ldots,v_7)
=\rho_3\sum_{(F_1,F_2,F_3;\nu)\in\mathcal T}
          \nu\,a_{F_1}^{(0)}a_{F_2}^{(0)}a_{F_3}^{(0)},
\tag{P3}
\]

with the following complete 19-term list. A string such as `0123`
denotes the face \((v_0,v_1,v_2,v_3)\).

| \(\nu\) | \(F_1\) | \(F_2\) | \(F_3\) |
| ---: | --- | --- | --- |
| -1 | 0123 | 3456 | 0167 |
| -1 | 0123 | 3456 | 1267 |
| -1 | 0123 | 3456 | 2367 |
| -1 | 0123 | 3456 | 3467 |
| -1 | 0123 | 3456 | 4567 |
| -1 | 0123 | 3467 | 4567 |
| +1 | 0123 | 3567 | 3457 |
| -1 | 0124 | 4567 | 2347 |
| +1 | 0125 | 2567 | 2345 |
| +1 | 0126 | 2367 | 3456 |
| +1 | 0127 | 2347 | 4567 |
| +1 | 0134 | 4567 | 1237 |
| +1 | 0145 | 1567 | 1234 |
| -1 | 0156 | 1267 | 2345 |
| +1 | 0167 | 1237 | 3456 |
| -1 | 0234 | 4567 | 0127 |
| +1 | 0345 | 0567 | 0123 |
| +1 | 0456 | 0167 | 1234 |
| +1 | 0567 | 0127 | 2345 |

This table is the tridegree-\((3,3,3)\) part of the finite integral
cyclic diagonal \(D_2\). With \(D_0\) the threefold Alexander–Whitney
diagonal, \(\varrho\) the Koszul-signed rotation moving the last
factor first, and \(h\) the tensor contraction to the first vertex, the
constructor uses

\[
D_i(\sigma)=h\bigl(R_iD_{i-1}(\sigma)+(-1)^iD_i(\partial\sigma)\bigr),
\quad R_1=\varrho-1,\quad R_2=1+\varrho+\varrho^2.
\]

Thus the sign and normalization in (P3) are fixed by the implementation,
not by fitting examples. The Bockstein convention is positive:

\[
\boxed{2\beta_{3,s}P^1_s\rho_3 A
       =\frac23d_s\widetilde{P^1_s\rho_3 A}.}
\]

In input degrees zero and one this operation vanishes by instability.
For an integral input in degree two, its cube has an integral lift, so
the Bockstein class vanishes; the implementation chooses its phase to be
zero. It only adds (P3) for \(n=3\).
Source: [`mod3_power.py`](../python/mod3_power.py).

## 8. Source primitives, transport, and scope

The three source primitives in the formulas above are fixed functions of
\((A,s,\omega)\). Schematically, for the fixed comparison \(j,r\) and
finite chain contraction \((F,G,H)\), the degree-two and degree-three
implementations use

\[
V_n(\sigma)=\Theta_{n+2}(p,k)\bigl(rHj\sigma\bigr)
                      +e_n(Fj\sigma),\qquad n=2,3.
\]

Degree one uses its dihedral Cayley-line contraction and the normalized
odd-branch phase. The functions \(e_n\), normalizations, finite operators,
and exact data are specified in [universal_helpers.md](universal_helpers.md).
They are not placeholders for solving \(d_sV_n=\Theta_{n+2}(p,k)\) on
the supplied space. In particular, the nonzero V2 source values
\((0,3/4,0,3/4,1/4)\) and the R3 universal source table remain part of
the operation.

GAP solves only the allowed defining-system equations for \(b,c\). It
then transports the **full phase** to and from the normalized bar model,
including the homotopy corrections to both defining cochains, and takes
the signed boundary in the resolution model. The kernel checks
\(ds=d\omega=0\), \(d_sA=0\), \(db=Da\), \(dc=\psi_n(A,b)\) on all
relevant faces; GAP checks exact divisibility and closure of the resulting
integral cochain. A failed identity is an error, not a request for a
different local R.

The universal choice and normalization are mathematical inputs to this
implementation. The saved
verification record (source-workspace provenance: `notes/extra/tertiary_T_degree3/gap_verification_20260923/README.md`; not bundled)
distinguishes direct degree-three coverage from page integration cases.
This document changes no formulas or certificates. The five-row AHSS
output remains an associated-graded calculation with `certified_ko:false`;
it does not resolve all ko rows or abutment extensions.
