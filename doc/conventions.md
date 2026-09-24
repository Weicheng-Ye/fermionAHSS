# Cochain conventions and transfer to GAP

This reference describes the formulas used by `koAHSSNaturalOperations()`.
The operation names and reading order are in [README.md](README.md).
All formulas below are on normalized simplicial cochains before transfer
to a GAP resolution. Coefficients, written order of products, and lifts are
part of the formulas.

## 1. Twists, coefficients, lifts and degrees

Let `s` and `omega` be fixed cocycles

\[
s\in Z^1(X;\mathbf F_2),\qquad
\omega\in Z^2(X;\mathbf F_2).
\]

The integral coefficient system is \(\mathbf Z_s\), with monodromy
\((-1)^s\). The phase coefficient system is \((\mathbf Q/\mathbf Z)_s\).
Actual phase computations use specified rational lifts, not floating point.

| Notation | Meaning |
| --- | --- |
| \(n\) | Input cochain degree; distinct from the optional page-count argument of `koAHSSpages` |
| \(p,q\) in an AHSS bidegree | Cohomological and coefficient-row indices |
| \(d\) | Ordinary integral coboundary, or its mod-two reduction where appropriate |
| \(d_s\) | Integral or rational coboundary with the sign coefficient system |
| \(\rho=\rho_2\) | Coefficient reduction modulo two |
| \([z]_2\) | Pointwise reduction of an integral-valued expression to its binary representative |
| \(\widetilde z\) | The pointwise \(0,1\)-valued integral lift of a mod-two cochain |
| \(\{r\}\) | Representative of a rational phase in \([0,1)\) |
| \(\smile_i\) | The specified interval-cut cup product; degree \(|x|+|y|-i\) |
| \(\smile=\smile_0\) | Ordered Alexander–Whitney product |

An equality of binary cochains is an equality modulo two. A rational phase
identity is modulo integral cochains unless explicitly called a literal
rational identity. An integral cocycle representative is a third object:
adding an integral cochain to a phase changes its integral boundary by a
coboundary, even though its cohomology class is unchanged.

Lifting a sum and adding its lifts are different operations:

\[
\widetilde{x+y}\ne\widetilde x+\widetilde y
\quad\text{in general.}
\]

The code distinguishes them. In particular the secondary correction
\(z=\widetilde{s\smile u}+\widetilde{\omega\smile e}\) is an **integer
sum**, while the nonlinear binary expressions \(F,G,\tau'\) are reduced
before being lifted. Division by two in an integral formula means exact
integer division after the indicated numerator has been assembled.

## 2. Coboundary and sign transport

In the frame at the first vertex, for \(c\) of degree \(r\),

\[
(d_sc)(v_0,\ldots,v_{r+1})
=(-1)^{s(v_0,v_1)}c(v_1,\ldots,v_{r+1})
 +\sum_{j=1}^{r+1}(-1)^j c(v_0,\ldots,\widehat v_j,\ldots,v_{r+1}).
\tag{C1}
\]

Equivalently, for the local numerical values used by the kernel,

\[
d_sc=dc-2\widetilde s\smile c.
\tag{C2}
\]

When the **right** factor takes values in the sign system, its ordered
product requires transport to the first frame:

\[
(x\smile_s y)(v_0,\ldots,v_{r+t})
=x(v_0,\ldots,v_r)(-1)^{s(v_0,v_r)}y(v_r,\ldots,v_{r+t}).
\tag{C3}
\]

This is `low_phases.transported(x,y,s)`. The coefficient system of the
product is the tensor product of those of its factors. For example,
the integral signed cube of a degree-two \(A\) is
\((A\smile_s A)\smile_s A\): the first two sign systems cancel, and
the third leaves a sign system. A product denoted simply \(\smile_i\)
in a displayed integer-lift formula uses its stated local numerical cup
formula; do not insert extra transports into it.

Sources: [natural_secondary.gi](../gap/natural_secondary.gi),
[phase_eval.py](../python/phase_eval.py),
[low_phases.py](../python/low_phases.py).

## 3. Word evaluation and integral cup signs

Let \(w=(w_1,\ldots,w_\ell)\) be a nondegenerate surjection onto
\(\{1,\ldots,r\}\), and let \(|x_j|=d_j\). Its output degree is
\(N=\sum_jd_j-\ell+r\). On \([v_0,\ldots,v_N]\), sum over cuts

\[
0=t_0\le t_1\le\cdots\le t_\ell=N.
\]

Assign the interval \([t_{a-1},t_a]\) to input \(x_{w_a}\). Retain
only cuts for which the concatenated vertices of each input are strictly
increasing and have exactly \(d_j+1\) vertices. Multiply the input values
on those faces and sum modulo two. This defines the word operation
\(w(x_1,\ldots,x_r)\).

For \(x\smile_i y\), use the alternating binary word of length \(i+2\),
starting with 1. Negative cup indices give zero. Integral cup evaluation
uses the same cuts with the following sign. For each interval set

\[
\epsilon_a=\begin{cases}1&\text{if }w_a\text{ occurs again later},\\0&\text{otherwise},\end{cases}
\qquad m_a=t_a-t_{a-1}+\epsilon_a.
\]

For input degrees \(p,q\), the cut sign is

\[
(-1)^{i(p+q)+\binom i2+
 \sum_a\epsilon_at_a+
 \sum_{a<b,\,w_a>w_b}m_am_b}.
\tag{C4}
\]

Thus the integral signs are fixed, including those in the secondary
\(B\smile_{n-1}B\) and Pontryagin-square formulas. This is the convention
in [natural_words.gi](../gap/natural_words.gi) and
[cochain_tools.py](../python/cochain_tools.py).

## 4. Cochain squares, Bocksteins and primary maps

For a binary cochain \(x\) of degree \(r\), define

\[
Q^j(x)=x\smile_{r-j}x+x\smile_{r-j+1}dx,
\qquad
E(x)=Q^2(x)+\omega\smile x,
\qquad
Q_D(x)=E(x)+s\smile Q^1(x).
\tag{C5}
\]

These expressions are binary, and are defined for nonclosed cochains.
For a cocycle \(a\), write \(\operatorname{Sq}^j(a)=a\smile_{n-j}a\)
for this cochain representative. The secondary formulas use a separate
ordinary Bockstein representative

\[
B(a)=\frac{d\widetilde a}{2},\qquad e(a)=\rho B(a).
\tag{C6}
\]

For a binary cocycle in this interval-cut convention, the equality
\(e(a)=Q^1(a)=a\smile_{n-1}a\) holds literally as binary cochains.
The integral cochains \(B(a)\) and \(C_B\) below still retain essential
carry information. For a noncocycle, use the full expression (C5);
\(d\widetilde x/2\) need not be integral. The literal cocycle equality is
also checked in the chi verification (source-workspace provenance: `note/extra/chi_suspension_degree6/verify.py`; not bundled).

The primary maps adjoining the higher differentials are

\[
\begin{aligned}
D(a)&=\operatorname{Sq}^2(a)+s\smile e(a)+\omega\smile a,\\
\operatorname{Dbar}(A)&=D(\rho A),\qquad A\in Z^n(X;\mathbf Z_s),\\
\operatorname{Dtilde}(a)&=\beta_s\bigl(\operatorname{Sq}^2[a]+\omega\smile[a]\bigr),
\quad
\beta_s[z]=\left[\frac{d_s\widetilde z}{2}\right].
\end{aligned}
\tag{C7}
\]

Here \(D,\operatorname{Dbar}\) have degree two and
\(\operatorname{Dtilde}\) has degree three. These names agree with the
GAP functions `D`, `Dbar`, and `Dtilde`. Their cohomological expressions
are \(D=\operatorname{Sq}^2+s\operatorname{Sq}^1+\omega\),
\(\operatorname{Dbar}=D\rho\), and
\(\operatorname{Dtilde}=\beta_s(\operatorname{Sq}^2+\omega)\).
They are evaluated on normalized group-bar cochains in the same convention
as the secondary formulas, then transferred to the supplied resolution.

The integral Pontryagin-square cochain used in the tertiary normalization is

\[
P_\omega=\widetilde\omega\smile\widetilde\omega
 +\widetilde\omega\smile_1d\widetilde\omega\in C^4(X;\mathbf Z).
\tag{C8}
\]

It is an ordinary integral lift of the mod-four Pontryagin square, so
\(dP_\omega\) is divisible by four. The legacy correction uses
\(q(\omega)=[dP_\omega/4]\). This \(q(\omega)\) is distinct from the
secondary integer cochain \(q(A)\) and the degree-three source variable
called `q` in the kernel.

## 5. Whole-formula transfer and defining systems

Let \(E_*=B_*(G)\) be the normalized homogeneous bar resolution of
the group underlying the actual HAP resolution \(R_*\). The package
constructs integral equivariant comparison maps \(f,g\) and a chain homotopy \(H\)

\[
f:E_*\longrightarrow R_*,\qquad g:R_*\longrightarrow E_*,
\qquad \partial H+H\partial=1-gf.
\tag{C9}
\]

Lift \(A_R,s_R,\omega_R\) by \(f^*\). If \(p_E=D(a_E)\), first solve
\(db_R=g^*p_E\) modulo two and set

\[
b_E=f^*b_R+H^*p_E\pmod2.
\tag{C10}
\]

Then \(db_E=p_E\). Let \(\tau_E\) denote the representative \(\tau'\)
of `Tau` on the bar resolution. For the tertiary operation, after arranging
\(dc_R=g^*\tau_E\), set

\[
c_E=f^*c_R+H^*\tau_E\pmod2.
\tag{C11}
\]

An input can survive the secondary quotient even when the initially chosen
\(b_R\) has nonzero \([g^*\tau_E]\). `koAHSSDefiningSystem` then changes it
by a closed \(h\) with \([Dh]=[g^*\tau_E]\), reevaluates the same natural
formula, and solves for \(c_R\). This chooses allowed defining cochains;
it does not choose \(\chi,\zeta\), or a residual \(R_n\) on the space.

Evaluate the **entire** secondary cocycle or tertiary rational phase on
\(E_*\), and project it with \(g^*\), including sign characters on the
group-ring coefficients. For a tertiary phase \(O_E\), the returned vector
is computed exactly as

\[
m=\operatorname{lcm}\{\text{denominators of }g^*O_E\},\qquad
T_R=\frac{d_s\bigl(mg^*O_E\bigr)}m.
\tag{C12}
\]

The comparison is built from the supplied HAP contraction and the bar cone.
A bar simplex is a homogeneous vertex list; consecutive equal vertices are
degenerate. Evaluation uses its first-vertex frame, and projection applies
the sign character to group-ring coefficients. The raw rational phase is
retained before its boundary is taken. Its audit separates the fractional
remainder, the projection carry, and the source-lift carry; the two carries
are integral and their sum with the remainder is the literal phase.
The bar comparison is fixed throughout this package and requires no model
selection option. No strict identity \(fg=1\) or normalized side conditions
on \(H\) are assumed for this comparison.

The numerator must be divisible by \(m\), and the result must be closed.
Evaluating separate universal words using an arbitrary HAP diagonal is not
this procedure. Sources:
[natural_bar.gi](../gap/natural_bar.gi),
[defining_systems.gi](../gap/defining_systems.gi), and
[natural_tertiary.gi](../gap/natural_tertiary.gi).

## 6. The interval convention

The tertiary free prisms and calibrations use the right interval factor.
For a degree-\(r\) cochain \(u\) on \(X\times I\), its prism is

\[
(Iu)(v_0,\ldots,v_{r-1})
=\sum_{j=0}^{r-1}(-1)^j
u((v_0,0),\ldots,(v_j,0),(v_j,1),\ldots,(v_{r-1},1)).
\tag{C13}
\]

The suspension-normalized slant is \(\kappa_r=(-1)^{r-1}I\), so
\(\kappa_6=-I\) and \(\kappa_7=+I\). A cochain written \(b\ell\)
in the free defining-cochain prism is evaluated as
\((b\ell)(z)=b(z)\ell(z_{\mathrm{last}})\), where \(\ell\) is the
interval vertex coordinate. This differs from suspension by the relative
degree-one interval cocycle. Preserve that distinction when reading the
R formulas and their calibration.
