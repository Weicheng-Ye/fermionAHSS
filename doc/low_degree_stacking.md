# Stacking in package degrees one and two

This note records the stacking product in package degrees one and two,
as it is evaluated by the native extension engine and by the complete-bar
reference model. The differential in these degrees is (A13)–(A14) of
[all_cochain_differential.md](all_cochain_differential.md), equivalently
(6)–(7) of [dimension_indexed_differentials.md](dimension_indexed_differentials.md).
Conventions for cup products, binary lifts \(\widetilde x\), the half lift
\(\mathsf h(x)=\widetilde x/2\) and the prism \(I\) are those of these notes;
\(\ell\) is the right interval coordinate and \(u\ell\) the pullback of u
times the interval coordinate of the last vertex.

## States and differentials

In package degree k the state is \((A,B,C,D)\) with cochain degrees
\((k-3,k-2,k-1,k+1)\); negative degrees are absent.

| k | State | Differential |
| ---: | --- | --- |
| 1 | \((C,D)\), \(C\in C^0(\mathbf F_2)\), \(D\in C^2(\mathbf Z_s)\) | \(\mathfrak d_1\) of (A14): \((\delta C,\ \delta_sD+\tfrac12(\delta_s\widetilde{\omega C}-\widetilde{\omega\delta C}))\) |
| 2 | \((B,C,D)\), \(B\in C^0\), \(C\in C^1\) binary, \(D\in C^3(\mathbf Z_s)\) | \(\mathfrak d_2\) of (A13): \((\delta B,\ \delta C+H_2(B),\ \delta_sD+g_2(B,C))\) |

Here \(H_2(B)=Q_{D,0}(B)=\omega B+s(B\,\delta B)\) is (A10), and \(g_2=J_2\)
of (A5) when B is closed on the whole complex (the locus \(L_2\)), and
\(g_2=G_2\) of (A12) otherwise. Both are written with the potential

\[
\Phi_2(B,C)=
\begin{cases}
\mathsf h(C\,\delta C+\omega C)+\mathsf h(t_2\smile_1\omega B),
  & \delta B=0,\quad t_2=\delta C+\omega B,\\
I\,\Omega_3\bigl(\delta(B\ell),\ \delta(C\ell)+H_2(B\ell)\bigr),
  & \delta B\ne0,
\end{cases}
\qquad
g_2=\delta_s\Phi_2-\Omega_3\bigl(\mathfrak d_2(B,C)_{B,C}\bigr),
\]

where \(\Omega_3(b,c)=\mathsf h(E_2c)+R_0(0,b)\) is the degree-three phase
at \(A=0\), \(R_0(0,b)=\Theta_1(b,0)\) of
[tertiary_operations.md](tertiary_operations.md). On \(L_2\) the second
argument of \(\Omega_3\) is \((0,t_2)\) and this is (7); off \(L_2\) prism
Stokes turns it into \(G_2=-I\,\delta_s\Omega_3\). In degree one the
same shape reads \(g_1=\delta_s\mathsf h(\omega C)-\mathsf h(\omega\,\delta C)\),
which is (A14).

## Products

The product is triangular. In degree one

\[
(C,D)\times(C',D')=\bigl(C+C',\ D+D'+\gamma_1(C,C')\bigr),
\]
\[
\gamma_1=\tfrac12\bigl(\widetilde{E(C)}+\widetilde{E(C')}-\widetilde{E(C+C')}\bigr)
+\delta_sP_1(C,C')+P_2\bigl(0,\delta C;0,\delta C'\bigr),
\]

with \(E(C)=\omega C\) in degree zero and the edge phase

\[
P_1(v_0v_1)=\tfrac32xu+\tfrac12yu-xyu+xv-xyv-2xuv-yuv+2xyuv+zyv,
\]

where \(x,y=C(v_0),C(v_1)\), \(u,v=C'(v_0),C'(v_1)\) and \(z=s(v_0v_1)\),
all as integers 0 or 1. \(P_2\) is the degree-two phase below.

In degree two

\[
(B,C,D)\times(B',C',D')=
\bigl(B+B',\ C+C'+\beta_2(B,B'),\ D+D'+\gamma_2\bigr),
\qquad
\beta_2(B,B')=\delta B\smile B'+s\smile(B\smile B'),
\]
\[
\gamma_2=\Phi_2(B,C)+\Phi_2(B',C')-\Phi_2\bigl(B+B',C+C'+\beta_2\bigr)
+\delta_sP_2(B,C;B',C')+P_3\bigl(\mathfrak d_2(B,C);\mathfrak d_2(B',C')\bigr).
\]

The three potentials use the closedness of \(B\), \(B'\) and \(B+B'\),
respectively. \(P_3\) is the legal degree-three phase at \(A=0\), the
repaired paper exchange phase with its coordinate gauges
(`paper_commutative_production.phase`), for which
\(\gamma_3=\Omega_3+\Omega_3'-\Omega_3(\text{sum})+\delta_sP_3\) on legal
lower data. The degree-two phase is its prism transport,

\[
P_2=-I\Bigl[P_3\bigl(b_I,c_I;b'_I,c'_I\bigr)
+\mathsf h\bigl(E\eta+c_I^{\rm tot}\smile_1\delta\eta
+\delta c_I^{\rm tot}\smile_2\delta\eta\bigr)\Bigr],
\]

with \((b_I,c_I)=\mathfrak d_2(B\ell,C\ell)_{B,C}\) and similarly for the
primed input, \(\eta=\beta_2(B\ell,B'\ell)+\beta_2(B,B')\ell\), and
\(c_I^{\rm tot}=c_I+c'_I+b_I\smile b'_I+s\smile(b_I\smile_1b'_I)\).
The last term transports C along the difference \(\eta\) between stacking
the two ramps and ramping their stack.

All values are exact rationals; \(\gamma_1,\gamma_2\) are checked to be
integral on evaluation. The implementation is `a0_degree1`,
`a0_degree2`, `a0_gamma` and `coherent_low_commutative` in
[python/stacking_model](../python/stacking_model/).

## Use in `koFull`

Degree-two relations are first measured by the low-degree adapter of
[extensions.md](extensions.md), whose bundled phase covers closed B and C
at \(\omega=0\). When it leaves a relation unresolved, for example for
nonzero \(\omega\), the degree is solved by the native engine of
[resolution-extensions.md](resolution-extensions.md) with these formulas:
states of degree two, gauges of degree one, the same comparison, zero test
and gauge search as in degrees 3–5. The adapter's result is kept as
`lowDegreeAttempt`.
