# The two secondary differentials used by GAP

This is a transcription of the production
[natural_secondary.gi](../gap/natural_secondary.gi), with the
`chi7_tail` helper family. See [conventions.md](conventions.md) for lifts,
integral cup signs, \(d_s\), and whole-formula transfer. The explicit helper
definitions and data are in [universal_helpers.md](universal_helpers.md).

The notes use the same names as the GAP callbacks:

| Differential | Name | Cochain output |
| --- | --- | --- |
| \(d_3:E_3^{n,0}\to E_3^{n+3,-2}\) | `Tau` | Binary \(\tau'\), degree \(n+3\) |
| \(d_4:E_4^{n,-1}\to E_4^{n+4,-4}\) | `Psi` | Integral \(\psi'\), degree \(n+4\) |

The other \(d_3\), from row \(-2\) to row \(-4\), is the primary
\(J=\beta_s(\operatorname{Sq}^2+\omega)\), called `Dtilde`.

## 1. Common input and the first defining cochain

Start with \(a\in Z^n(X;\mathbf F_2)\). For integral input take
\(a=\rho A\), where \(A\in Z^n(X;\mathbf Z_s)\). Set

\[
\begin{aligned}
B&=\frac{d\widetilde a}{2},& e&=\rho B,&
C_B&=\frac{B+\widetilde e}{2},\\
u&=a\smile_{n-2}a,&v&=\omega\smile a,&w&=s\smile e.
\end{aligned}
\tag{S1}
\]

Thus \(|B|=|e|=n+1\) and \(|u|=|v|=|w|=n+2\). The **plus** sign
in \(C_B\) is intentional. Products retain their written order.

Choose a binary \(b\in C^{n+1}(X;\mathbf F_2)\) with

\[
db=u+v+w=D(a).
\tag{S2}
\]

On a HAP resolution this equation is solved for the transferred primary,
then the bar defining cochain is corrected by the fixed chain homotopy;
see (C10) in [conventions.md](conventions.md). If the primary obstruction
does not vanish, this secondary operation is not defined on that input.

## 2. The common binary formula and integer cochain

The binary cochain \(F\) of degree \(n+3\) is

\[
\begin{aligned}
F={}&Q^2(b)+\omega\smile b
 +\zeta_{2,n}(\omega,a)+\chi_n(a)\\
 &+u\smile_{n+1}v+u\smile_{n+1}w+v\smile_{n+1}w\\
 &+\zeta_{1,n+1}(s,e)+(\omega\smile_1s)\smile e
 +s\smile u+s^2\smile\rho C_B.
\end{aligned}
\tag{S3}
\]

The entire right side is reduced modulo two. In particular
\(Q^2(b)=b\smile_{n-1}b+b\smile_n db\), because \(|b|=n+1\).
The integer cochain of the same degree is

\[
q=\widetilde\omega\smile B+B\smile_{n-1}B.
\tag{S4}
\]

Both products in (S4) use the specified **integral** cup signs. This \(q\)
is not the degree-five characteristic class also denoted \(q(\omega)\).
The uncorrected integral secondary expression is

\[
\psi_{\rm ref}
=\frac12\left(d_s\widetilde F+\frac{d_sq}{2}\right)
=\frac{d_s(2\widetilde F+q)}4.
\tag{S5}
\]

All the displayed divisions are exact for a valid defining system and
the prescribed helpers. The implementation checks the final projected
numerator rather than rounding or replacing a failed identity.

## 3. Final degree-four operation: `Psi`

Define the **integer sum**

\[
Z=\widetilde{s\smile u}+\widetilde{\omega\smile e}.
\tag{S6}
\]

The calibrated cochain is

\[
\boxed{\psi'
=\frac{d_s(2\widetilde F+q+2Z)}4
=\psi_{\rm ref}+\frac{d_sZ}{2}.}
\tag{S7}
\]

Its cohomological normalization is

\[
\operatorname{Psi}_n(a)
=\operatorname{Psi}_{{\rm ref},n}(a)
 +\beta_s\bigl(s\operatorname{Sq}^2[a]\bigr)
 +\beta_s\bigl(\omega\operatorname{Sq}^1[a]\bigr).
\tag{S8}
\]

Relative to `chi7_tail`, the coefficients of
\(\beta_s(s\operatorname{Sq}^2a),\beta_s(s^3a),
\beta_s(\omega\operatorname{Sq}^1a)\) are
\((\eta_2,\eta_3,\eta_4)=(1,0,1)\).

GAP projects \(2\widetilde F+q+2Z\) using the integral sign-coefficient
chain map, differentiates, checks divisibility by four, and returns the
quotient. It also checks \(d_s\psi'=0\).

## 4. Final degree-three operation: `Tau`

Now require \(d_sA=0\), with \(a=\rho A\). Define

\[
K=\frac{A-\widetilde a}{2},\qquad
t=\rho K,\qquad x=dt,\qquad y=s\smile a.
\tag{S9}
\]

The exact lower relation is \(dt=e+s\smile a\) modulo two. Put

\[
\begin{aligned}
G={}&Q^2(t)+\omega\smile t+x\smile_n y
 +\zeta_{1,n}(s,a)+(\omega\smile_1s)\smile a+s\smile b,\\
L={}&\frac{q-d_s\widetilde G}{2},\qquad
M=\widetilde F+L,\qquad \kappa=s^3\smile a.
\end{aligned}
\tag{S10}
\]

Here \(G\) is binary of degree \(n+2\), and \(L,M\) are integer
cochains of degree \(n+3\). The final binary result is

\[
\boxed{\tau'=\rho M+\kappa
=F+\rho\!\left(\frac{q-d_s\widetilde G}{2}\right)+s^3\smile a.}
\tag{S11}
\]

Thus, relative to this helper family,

\[
\operatorname{Tau}_n(A)=\operatorname{Tau}_{{\rm ref},n}(A)+s^3\rho A,
\qquad(\epsilon_3,\epsilon_4,\epsilon_5)=(1,0,0).
\tag{S12}
\]

The tuple refers to the basis
\(s^3\rho A,s\omega\rho A,\operatorname{Sq}^1(\omega)\rho A\).
It must not be transferred unchanged to a different chi family.

## 5. The exact matched integral lift

The two calibrated operations retain a common lift. Define the binary
degree-\(n+2\) cochain

\[
h=s\smile b+s^2\smile t+\omega\smile t
 +(\omega\smile_1s)\smile a,
\qquad
M'=M+Z-d_s\widetilde h.
\tag{S13}
\]

The cup-one term is part of the formula. The identities are

\[
dh=\rho Z+\kappa,\qquad
\rho M'=\tau',\qquad d_sM'=2\psi',\qquad
d\tau'=0,\qquad d_s\psi'=0.
\tag{S14}
\]

For integral input the evaluator stores the projected \(M'\) and computes
the matched \(\psi'\) as \(d_sM'/2\). For mod-two input it uses (S7).
These are compatible representatives. The literal bar \(\tau'\), along
with \(A,b,s,\omega\), is passed on to the tertiary defining system.

The square-zero assertion is the cohomological statement

\[
J[\tau']=0\in H^{n+6}(X;\mathbf Z_s),
\qquad J=\beta_s(\operatorname{Sq}^2+\omega).
\tag{S15}
\]

It need not be a literally zero cochain vector. The mathematical
identification and suspension argument for this calibrated family are
recorded in the degree-seven companion (source-workspace provenance: `note/extra/chi_suspension_degree7/README.md`; not bundled).
Finite implementation tests do not replace that universal argument.

## 6. Domains, indeterminacy and the actual page groups

Writing cohomology groups with their coefficients, the maps are

\[
\operatorname{Tau}_n:
\ker\!\left(\overline D:H^n(X;\mathbf Z_s)\to H^{n+2}(X;\mathbf F_2)\right)
\longrightarrow
\frac{\ker\!\left(J:H^{n+3}(X;\mathbf F_2)\to H^{n+6}(X;\mathbf Z_s)\right)}
{D H^{n+1}(X;\mathbf F_2)},
\tag{S16}
\]

\[
\operatorname{Psi}_n:
\frac{\ker\!\left(D:H^n(X;\mathbf F_2)\to H^{n+2}(X;\mathbf F_2)\right)}
{\overline D H^{n-2}(X;\mathbf Z_s)}
\longrightarrow
\frac{H^{n+4}(X;\mathbf Z_s)}{JH^{n+1}(X;\mathbf F_2)}.
\tag{S17}
\]

Negative-degree groups are zero in this range. A change of the allowed
nullhomotopy \(b\) is handled by these image quotients. The page engine
uses actual abelian-group kernels, images, lifts and quotient maps; it
does not replace them by dimension or rank counts.

## 7. Runtime range and entry points

```gap
operations := koAHSSNaturalOperations();;
space := koAHSSHAPSpace(resolution, operations);;
backend := space.koAHSS(s, omega, resolutionLength);;

tau := koAHSSNaturalSecondary(backend, n, A);;
psi := koAHSSNaturalSecondary(backend, n, a, rec(inputType := "mod2"));;
# An optional b in the options record supplies a checked resolution cochain.
```

Both secondary evaluators have the current finite helper family for
\(0\le n\le7\). Larger degrees report it unavailable. The allowed
physical cutoff \(k\le6\) needs at most `Tau_3` and `Psi_4`, including
the degree-four secondary map in the denominator of `T_3`.
`koAHSSNaturalOperations()` now also installs final `T`; its formulas are
in [tertiary_operations.md](tertiary_operations.md).
