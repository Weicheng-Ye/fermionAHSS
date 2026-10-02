# Reducing bar transport across the package

The reductions in [light transport](light-transport-reduction.md) are not
specific to B-over-D. This note records their scope across the page engine,
all eight light-worker tasks, the transferred-model fallback and the
universal contractors. It separates a mathematical replacement from its
implementation. The [implemented defining systems](light-transport-reduction.md#10-implemented-defining-systems)
use the primary comparison for B atoms and A systems, the A gauge's
primary-plus-Bockstein source, exact lower-degree B carry homotopies,
the integral A-gauge identity, the pure-C marked class, and native
three-primary rows when the necessary native degrees are available.
The target-character, additional absorption and calibrated higher
comparison constructions below remain opportunities, not runtime paths.

There are two different costs. The input-resolution comparison constructs
or evaluates \(f,g,H\). A formula evaluated on one of its simplices can
itself invoke a universal contraction, such as \(V_2,V_3\) or a pair
primitive. Removing an input comparison homotopy does not remove that
second cost.

## 1. Precision and the scope of the light-transport arguments

Use \(\Lambda=f^*\), \(\Pi=g^*\) and \(H=(h')^*\) for the
**normalized extension comparison**. It satisfies

\[
 1-\Lambda\Pi=\delta H+H\delta,\qquad
 \Pi\Lambda=1,\qquad H\Lambda=\Pi H=0.
 \tag{A1}
\]

These side conditions do not hold by assumption for the page engine's
raw comparison. See [the page convention](conventions.md#5-whole-formula-transfer-and-defining-systems)
and [the extension normalization](transfer.md).

Three outputs require different amounts of information:

| Required output | Permitted replacement |
| --- | --- |
| The same marked cochain | A literal identity, including every lift and carry |
| The same marked cohomology class | An integral coboundary may be omitted |
| An extension row in \(H_{\rm low}/mH_{\rm low}\) | Changes by \(mH_{\rm low}\), and coherent changes of lower marking, may also be allowed |

In particular, the adjacent primary rows already avoid transport, but
their shifted representatives can require a more precise row when an
upper relation refers to them. A reduction of that later computation must
respect its marking; recovering just its parity again is insufficient.

The arguments of the light-transport note have the following scope.

| Argument | Scope beyond B-over-D | Qualification |
| --- | --- | --- |
| Bockstein carry, (2)–(6) | Every closed binary input; all supported degrees | Its old-marking witness is a lower-degree homotopy, not literal equality of the primitives |
| \(H\Lambda=0\), (7) | Every A-relation primitive of \(mA\), in degrees 4–6 | Extension comparison only; also covers a subsequently added native cocycle |
| Binary diagonal comparison, (8)–(13) | Primary defining cochains for B atoms, A systems and closed D-gauge inputs | Match the native Sq1 representative as in (11a); does not yet compare the calibrated secondary source |
| Integer C-change carry, (13a)–(13b) | Legal A,B in degrees 3–6 | Conjugate products and gauges too if changing the model's coordinates |
| Pure-C marked class, (14) | Degrees 1–6 | Keeping \(e_R=\Pi\widetilde{E(C)}\) preserves the specified completion |
| Single torsion A and \(A=s\), (15)–(16) | Package degree 4 | The torsion statement about \(H^1\) does not extend to \(H^2,H^3\) |
| Absorption, (17)–(18) | B rows, A rows and quotients truncated at C | Requires an appropriate character and an annihilated bottom subgroup; several rows require one simultaneous change |
| Target pairings and finite precision, (19)–(21) | Every projection, including page outputs and model fallbacks | Fewer pairing targets is not a bound on their support; rational divisions retain their carries |

The Sq1 qualification is substantive. The native primary operation in
[cochains.gi](../gap/cochains.gi) uses the ordinary integral Bockstein
\(e_R=\rho(\delta\widetilde b_R/2)\). A chosen higher diagonal instead
gives \(S_R^1b_R=b_R\smile_{n-1}^Rb_R\). Their classes agree. To
substitute `backend.nativePrimary` into the primary comparison, choose
\(t_R\) with \(\delta t_R=S_R^1b_R+e_R\) and add
\(\Lambda(s_R\smile^Rt_R)\) to \(K_D\). Omitting that term need not
give the prescribed defining equation as a cochain.
The light constructor uses the cup representative consistently in both
\(D_R\) and \(K_D\); it does not substitute `nativePrimary` here.

An absorption with an unevaluated tail proves the existence of a filtered
isomorphism. It does not supply the coordinates of that isomorphism in
the old marking. Such a result can determine invariant factors while
still failing an interface that promises an explicit old-marking witness.

## 2. A common replacement for defining-system homotopies

Suppose a closed source has a **proved comparison**

\[
 z=\Lambda z_R+\delta K.
 \tag{A2}
\]

Choose \(\delta r_R=z_R\). Then

\[
 P_{\rm new}=\Lambda r_R+K,\qquad \delta P_{\rm new}=z.
 \tag{A3}
\]

There is no \(Hz\) or \(\Pi z\) in this construction. To compare to
\(P_{\rm old}=\Lambda r_{\rm old}+Hz\), where
\(\delta r_{\rm old}=\Pi z\), use

\[
 r_R=r_{\rm old}-\Pi K,\qquad
 P_{\rm old}=P_{\rm new}-\delta HK.
 \tag{A4}
\]

Indeed \(\Pi z=z_R+\delta\Pi K\) and
\(H\delta K=K-\Lambda\Pi K-\delta HK\). Over \(\mathbf F_2\)
all minus signs become plus signs. The old-marking witness uses H one
degree below the old source.

This is the common mechanism behind the Bockstein carry and \(K_D\).
It also applies to sums of sources with known comparisons. For example,
the B-relation source

\[
 x+\Lambda c_R+P,\qquad x=h^D(b,b),
\]

already has an explicit comparison for its first two terms. Only the
uncompared part P still needs H. If P is a primary D source of a closed
native input, \(K_D\) compares it too. A Tau source needs a calibrated
secondary comparison; primary diagonals alone do not provide it.

For A relations \(U=P(m\Lambda A_R;u_R)=\Lambda u_R\) literally.
The gauge's binary input \(v=\rho U\) is closed because m is even.
Although U is not an integral cocycle, v is therefore a valid input for
the primary comparison. The same applies to the closed-input `QDv` and
`ypp` branches of `a_step`. The remaining terms of the Y source,
including the lower reference, must be compared as well before removing
its entire homotopy.

For the A system itself this removes the first primitive
\(B=P(Q_D(\rho A);B_R)\). It does not by itself remove
\(C=P(f^\sharp(A,B);C_R)\), the relative-C primitive W, or their
secondary/tertiary phases. Establishing (A2) for those sources is exactly
the missing higher-comparison problem, not a license to choose a new
primitive on R.

### Matching the raw page comparison

Construction (A3) itself only needs the chain-map identity for
\(\Lambda\), so it also constructs valid page defining cochains.
Equation (A4), however, uses the extension side conditions. To obtain
its page analogue, put \(P_R=\Pi\Lambda\) and retain comparison
homotopies on R and from R to the bar cochains:

\[
 1-P_R=\delta_RS_R+S_R\delta_R,\qquad
 H\Lambda-\Lambda S_R=\delta L-L\delta_R.
 \tag{A4p}
\]

Here \(|S_R|=-1\) and \(|L|=-2\). The second operator compares the
two homotopies of \(\Lambda-\Lambda P_R\). Such chain data can be
constructed degree by degree using the free resolutions and their
contractions; they are not implied by setting \(S_R=L=0\).
For a closed \(z_R\), use

\[
 r_R=r_{\rm old}-\Pi K+S_Rz_R,\qquad
 P_{\rm old}=P_{\rm new}+\delta(Lz_R-HK).
\]

The first identity gives \(\delta_Rr_R=z_R\); the second follows by
substituting (A4p) in the calculation preceding (A4). This specifies the
additional data needed to match the old page representative. If only the
page operation is required, one may instead construct a valid complete
defining system using (A3) and evaluate the fixed bar formula with it,
respecting its stated indeterminacy. A literal cochain or matched-lift
comparison requires the extra homotopies.

## 3. Pure-C rows with a completely native new marking

The standalone C-over-D extension residue in \(D_k/2D_k\) is already
native, both in [the existing cup-i formula](extension_cup_i_formulas.md#5-c-over-d)
and in `primaryRow`. Formula (A7) below is that same native formula;
the additional argument identifies its marking by an explicit comparison
carry. It is not a new bypass for the ordinary C-over-D row.

The one-projection formula (14) of the light-transport note preserves the
old marked C class. There is also a fully native formula in an explicitly
specified new marking.

Let \(c_R\in Z^{k-1}(R;\mathbf F_2)\), \(C=\Lambda_2c_R\), and
write \(E_Rc_R=\operatorname{Sq}^2_Rc_R+\omega_R\smile^Rc_R\).
The binary primary comparison gives

\[
 \delta K_E=E(C)+\Lambda_2E_Rc_R,\qquad
 K_E=K_{k-3}(c_R,c_R)+K_0(\omega_R,c_R).
 \tag{A5}
\]

Negative-index terms vanish. Define the integral carry

\[
 L_E=\frac{\widetilde{E(C)}-\Lambda_s\widetilde{E_Rc_R}
                    -\delta_s\widetilde K_E}{2}.
 \tag{A6}
\]

Its numerator is even by (A5), and

\[
 \delta_sL_E=\beta_sE(C)-\Lambda_s\beta_s^RE_Rc_R.
\]

Solve \(\delta_s^RD_R=-\beta_s^RE_Rc_R\) natively. The bar state
with C coordinate C and D coordinate \(\Lambda_sD_R-L_E\) is flat.
For this marking its double has native cohomology class

\[
 \boxed{2D_R+\widetilde{E_Rc_R}
                  +\beta_s^R\operatorname{Sq}^1[c_R].}
 \tag{A7}
\]

To prove this, substitute (A6) into
\(2(\Lambda D_R-L_E)+\widetilde{E(C)}
 +\beta_s\operatorname{Sq}^1C\); the terms involving \(L_E\)
leave \(\Lambda(2D_R+\widetilde{E_Rc_R})
 +\delta_s\widetilde K_E\). The last term is an integral coboundary,
and the Bockstein is natural. Thus (A7) requires neither g nor H for
the pure-C row. The comparison carries define why this new marking is
valid; they need not be evaluated to obtain that isolated row.

An upper row referring to this C generator must use this same D
coordinate, or a proved symbolic cancellation of its carries. Using
(A7) as if it were the old marked C row silently loses \(\Pi L_E\).
This is why the existing native primary shortcut is sufficient for a
standalone C-over-D residue but does not automatically replace every
`c_mark` request.

The integer carry for \(C\mapsto C+\delta\eta\) in (13a)–(13b)
also extends through degree six. On legal A,B,

\[
 \Omega(A,B,C)=\tfrac12\widetilde{E(C)}+R(A,B)
       +\tfrac12\widetilde{t\smile_{k-1}\tau},
 \quad t=\delta C+\tau,
\]

where the last term is the legal-A,B correction and \(\tau=f^\sharp(A,B)\).
Only \(E(C)\) changes under \(C\mapsto C+\delta\eta\). Hence the
same carry changes D so that the curvature agrees, including at k=5,6;
no new evaluation of \(R_2\) or \(R_3\) is needed for this difference.

## 4. Absorption for several rows and for other target layers

Let H be an already computed lower group, D its bottom subgroup,
\(K=H/D\), and suppose \(pD=0\) for a prime p. Consider rows

\[
 m_i a_i=t_i+d_i,\qquad d_i\in D.
\]

If the known images \(\bar t_i\) are linearly independent in
\(K/pK\), choose characters \(\ell_i:K\to\mathbf F_p\) with
\(\ell_i(\bar t_j)=\delta_{ij}\). Then

\[
 \varphi(h)=h-\sum_i\ell_i(\bar h)d_i
 \tag{A8}
\]

is an automorphism of H, with inverse given by the plus sign, fixing D
and K. It sends every \(t_i+d_i\) to \(t_i\). This proves simultaneous
absorption of all those tails, for any orders \(m_i\). Other rows must
be transformed by the **same** automorphism. If
\(\bar t_j=\sum_i a_i\bar t_i\) modulo p, its tail becomes

\[
 d_j-\sum_i a_i d_i.
 \tag{A9}
\]

That residual tail cannot in general be dropped. This is the reason to
adapt rows to pivots and kernel rows before evaluating expensive phases.

Applications include:

- The existing B-over-D absorption when \(2D=0\), with independent
  leading C rows. It is already implemented in the light frame.
- A-over-D in degrees 5 and 6, when the leading rows are independent in
  \((H/D)/2(H/D)\). Nonzero B coordinates are sufficient for individual
  detection, but dependence between several rows still matters.
- An A-over-C computation truncated below C: use the lower group modulo
  D and its bottom C subgroup. C is killed by two; independent leading B
  rows can absorb C tails in this truncated presentation. This alone
  proves nothing about the later D tails.
- A three-primary upper row only when there is an appropriate nonzero
  lower quotient K. In the ordinary three-local A/D tower the lower
  group is D itself, so K is zero and this argument gives no reduction.
  The \(\mathbf Z/9\) extension cannot be erased this way.

At k=4 the A torsion has at most one generator, so the single-row result
needs no simultaneous top-row analysis. Higher-degree A torsion does not
have that restriction. Free D generators and D elements of order four
or larger also prevent applying the elementary-two statement without
additional hypotheses.

## 5. Three-primary operations can be evaluated natively

For a three-primary A generator of order \(m=3^e\), the light row is

\[
 2\cdot3^{e-1}Y,\qquad
 \rho_3[Y]=P_s^1\rho_3[A].
 \tag{A10}
\]

At k=5, \(|A|=2\), and take

\[
 Y=(A\smile_R A)\smile_R A.
 \tag{A11}
\]

The first cup has ordinary coefficients and the second has sign
coefficients. Naturality of the cup product identifies its integral
class with the projected bar cube. Thus this replaces `p3_power` without
an input lift or projection, at the precision of the marked cohomology
class. Different parenthesizations differ by a coboundary on cocycles.

At k=6, \(|A|=3\), use the cyclic diagonal \(\Delta_2^{(3),R}\):

\[
 P_{s,R}^1(a)=
 \langle a\otimes a\otimes a,\Delta_2^{(3),R}\rangle\pmod3.
 \tag{A12}
\]

The engine [native_coherence.gi](../gap/native_coherence.gi) already
constructs `cyclic(0)`, `cyclic(1)` and `cyclic(2)` with boundaries
\(\rho-1\) and \(1+\rho+\rho^2\), the same cyclic convention as
[mod3_power.py](../python/mod3_power.py). The tensor evaluation with three
sign characters supplies the local system. The standard comparison of
cyclic diagonals identifies this primary class with the nineteen-term
bar formula; no secondary or tertiary calibration is being selected.
Its present constructor requires total tensor degree nine and a spare
contraction degree, hence a resolution through degree ten for this call.

Lift the resulting mod-three class to an integral cocycle using the
same coefficient-sequence solve as the light frame. Survival supplies
the required integral lift. Two such lifts differ in cohomology by
\(3z\); multiplication by \(2\cdot3^{e-1}\) changes the row by
\(2m z\), which vanishes in the required extension quotient. Therefore
(A12) removes the remaining transport of the three-primary light row.

For the three-local **fallback model**, not just its row, use coherent
native cyclic diagonals in the existing formulas for curvature,
cross-effect and coboundary primitives. The cross-effect needs
\(\Delta_3^{(3),R}\), one stage beyond the currently exposed API; its
boundary is \((\rho-1)\Delta_2^{(3),R}\). The tensor filler can
construct it. The model's divisions by three and the coordinate
comparison still have to be retained. Simply substituting an unrelated
cup into the degree-five polynomial while keeping all its old gauge
formulas would not establish the same marked model.

For page T, the degree-two three-primary contribution is already zero
as a cohomology operation: the mod-three cube has an integral cocycle
lift. It remains essential as a stacking phase. The degree-three
three-primary page contribution can use
\(2\beta_3^RP_{s,R}^1\rho_3A\), with the same qualification about
cochain versus cohomology output.

## 6. Low-degree secondary and tertiary opportunities

### Even integral input to Tau

If \(A=2U\) is a cocycle, choose the primary defining cochain b=0.
In (S1)–(S11) of [the secondary formulas](secondary_operations.md),
\(a=0\), \(F=q=0\) and \(G=E(\rho U)\). Consequently

\[
 \boxed{\operatorname{Tau}(2U)
       =\rho\,\operatorname{Dtilde}(\rho U)}
 \tag{A13}
\]

as a page class, modulo the D indeterminacy. The minus sign in (S11)
disappears modulo two. This uses only a native primary operation in every
supported input degree. For a tertiary defining system, changing b or
the literal Tau representative requires its comparison and matched
integral lift; the page-class identity alone does not define T.

### Psi in input degree zero

For BG a closed binary degree-zero input is constant. The zero input
has zero value; for input one the primary condition requires
\(\delta b=\omega\). Here the secondary formula reduces to

\[
 F=b\smile\delta b+\delta b\smile b
   =\delta(b\smile b)\pmod2,\qquad \psi'=\beta_sF.
\]

The Bockstein of this coboundary is an integral coboundary. Thus
\(\operatorname{Psi}_0=0\) as a cohomology operation whenever defined,
for either sign system. Its literal cochain need not be zero. In positive
degrees Psi remains a secondary operation; this degree-zero simplification
does not justify replacing it by zero.

### Degree-four A torsion and degree-zero T

At k=4 the torsion A representative can be chosen as A=s, with U=-1.
This eliminates the relation primitive and fixes its leading B row as
\([\omega]\). When that class is nonzero and \(2D=0\), absorption
can avoid the remaining D residue for the isolated top relation. When
the detecting character does not exist, the C/D data remain essential.
With \(\omega=0\), the known A-over-C formula still contains
\(h^D(B,B)\) with \(\delta B=s^3\). B is not closed, so the
closed-input Bockstein shortcut does not apply to it.

T in input degree zero has explicit even and odd R0 formulas;
it does not call V1–V3. With trivial sign its page callback uses the
[native Pontryagin-square reductions](tertiary_operations.md#native-page-classes-with-trivial-sign)
for multiples of four and for twice-odd inputs whose omega has a mod-four
lift. Defined positive odd inputs give zero. The other sectors retain
their finite cup/carry formula and defining-cochain comparisons.
Vanishing negative-index cups and \(\chi_0=\zeta_{1,0}
=\zeta_{2,0}=0\) make this the smallest higher-comparison problem.
It is not a consequence of the primary comparison alone.

For T1–T3 and two-primary A-over-D, a fully native replacement still
requires a comparison of the calibrated higher formulas. Tensor fillers
can supply associators, Adem homotopies and further cells; their existence
does not identify the prescribed `chi7_tail`, matched lift, R selectors
and universal periods. Those choices must be matched before the runtime
can use them.

## 7. Pair only the quotient that is required

For any rational-potential expression

\[
 r_R=\Pi_sF+\delta_s^R\Pi_sP
\]

and any integral column \(\ell\),

\[
 r_R\ell=\langle F,g\ell\rangle+
                  \langle P,g(M\ell)\rangle.
 \tag{A14}
\]

This is valid for pages as well as extensions: it needs a chain map,
not the strict side conditions (A1). Select characters of the current
page quotient or of \(D/(D\cap mH)\) first. The kernel and incoming
images are computed on R. An elementary-p quotient needs only its
independent character pairings, not an ambient cochain vector.

If \(M\ell=0\) integrally, the entire P term disappears for this
target. If \(M\ell=m v\), it becomes \(\langle mP,gv\rangle\).
It must not be dropped just because \(M\ell=0\pmod m\): the rational
denominator can turn that multiple of m into the nonzero Bockstein being
measured. On a rational cohomology target, in contrast, a closed
rational-phase boundary has zero class; only torsion targets can detect
these operations.

For a cyclic target of order \(p^e\), use coefficient-sequence or Smith
presentation functionals of the appropriate precision. Binary characters
alone need not distinguish that target. Nor can arbitrary characters on
an extension quotient be obtained by reading the displayed coordinates
before accounting for incoming images and the lower presentation.

The chain builder should combine the target before expanding it. Keep
group-ring actions until after any non-equivariant cone, and collect
coefficients in the required modulus between linear stages. Computing
every g-basis column first loses the main potential saving. A scalar
adjoint evaluator can also compose \(q,h,\partial,q\) without storing
the whole normalized-homotopy chain. This is an exact linear
reorganization, not a guarantee that the scalar recursion is short.

If the assembled answer is \(N/q\in\mathbf Z\) and is needed modulo
m, keep N modulo mq. For nested divisions propagate this precision
backward through each division. Do not reduce the integer labels of a
universal contractor on this basis: the contractor is not a linear
function of those labels and may not be periodic in them.

A whole phase whose denominators are prime to p is an ordinary cochain
over \(\mathbf Z_{(p)}\); its coboundary has zero p-local cohomology
class. This can simplify a p-primary class or a coherently rephased
model. It does not allow dropping one summand inside a half-lift, or
dropping its contribution to a fixed integral marked cochain.

## 8. Inventory of page transport

Implementations: [natural_secondary.gi](../gap/natural_secondary.gi),
[natural_tertiary.gi](../gap/natural_tertiary.gi),
[natural_bar.gi](../gap/natural_bar.gi) and [worker.py](../python/worker.py).

| Formula or site | Replacement or reduction | What remains |
| --- | --- | --- |
| `primary("D")`, `primary("Dbar")`: \(\Pi Q_D(\Lambda a)\) | `nativePrimary` gives the same cohomology class using native cups and Bockstein | If this representative enters a higher defining system, include \(K_D^{\rm native}\) |
| `primary("Dtilde")`: \(\delta_s^R\Pi\widetilde E/2\) | Native \(\beta_s E_R\) gives the same primary class | A fixed cochain comparison requires the carry (A6) |
| Common secondary defining cochain \(b=\Lambda b_R+H D(a)\) | Use the primary comparison (A3); native source solve plus \(K_D\) | Page comparison lacks strict retraction; preserve its defining choices with a comparison, rather than using (A4) unchanged |
| Tau: projections of F, q, G, the matched lift and \(s^3a\) | (A13) eliminates transport for even input at page precision; otherwise quotient characters, binary precision and a calibrated native secondary comparison | General Tau is not determined by primary cup products alone |
| Psi: \(\delta_s^R\Pi(2\widetilde F+q+2Z)/4\) | Zero class in input degree zero; pair the potential with \(g(M\ell)\) in other degrees | Retain modulus \(4m\) for a residue modulo m; preserve the matched lift when used elsewhere |
| T defining cochain \(c=\Lambda c_R+H\tau'\) | Apply (A3) once a calibrated comparison for Tau is known; even-input Tau is a useful first sector | A representative of the same Tau class is insufficient without its defining-system change |
| T0 phase \(\mathsf h(Ec)+R_0\) | With trivial sign: native Pontryagin Bocksteins for multiples of four and twice-odd inputs with a mod-four omega lift; zero for defined positive odd inputs | The other sectors retain their defining-cochain comparisons; direct phase audits retain the bar formula |
| T1 phase \(\mathsf h(Ec)+R_1\) | Quotient pairings; specialize the sign A torsion sector where applicable; aggregate the V1 unit-edge sums | The prescribed odd normalization and R1 universal periods |
| T2 phase, including \(-A^3/4\) and \(\tfrac23P^1\rho_3A\) | Native cube for the three-primary class; quotient pairings and direct finite formulas for the dyadic part | The full dyadic V2 and prism phases; the three-primary stacking carry even though its page map vanishes |
| T3 phase, including V3 and \(\tfrac23P^1\rho_3A\) | Native cyclic diagonal for the three-primary class; functional projection and fixed-source contraction optimizations for the dyadic part | Calibrated R3 and V3; no general native dyadic replacement established |
| Twist/input lifts and face sampling for all these formulas | Eliminate them when the entire operation is native; otherwise memoize only requested faces and reuse common input evaluations | Nonlinear operations cannot be commuted through \(\Pi\) factor by factor |

The page engine already projects the Psi/T potential before taking its
native coboundary. That saving must not be counted as a new bypass of
the remaining projection.

## 9. Inventory of the light-worker tasks

Every task in `LightEvaluator.TASKS` is covered here. Definitions are in
[extension_light.py](../python/extension_light.py), selection and markings
in [extension_light.gi](../gap/extension_light.gi).

| Task and transported terms | Applicable reduction | Boundary of the result |
| --- | --- | --- |
| `qd`: native primary right-hand side | Implemented: D from the matched native cups; the generic \(\Pi Q_D(b)\) path is optional | The right-hand side is paired with \(\Lambda c_R+K_D\), not with an unchanged generic `atom` primitive |
| `atom_curvature`: \(C=\Lambda c_R+K_D\), then \(\Pi J_k(0,b,C)\) | Implemented: no primary-source homotopy; normalized rational potential projected in degree k+1 | Flat-admissibility and the potential still need evaluation; target-character projection is a further opportunity |
| `c_mark`: projection of \(\widetilde E(C)\) | Implemented: (14), with native Bockstein, unchanged curvature and D completion | (A7) would remove the remaining projection with a coherent new marking |
| `gauge`: D/Tau gauge sources and \(x=h^D(b,b)\) | Implemented: \(\Pi x=\rho q_R+\delta_R\Pi\rho\kappa_b\), with the projection one degree lower; closed primary defining sources use \(K_D\) | Only a calibrated secondary comparison removes a Tau-source homotopy |
| `bd_page`: secondary page residue of a kernel B atom | Carry primitive for \(\pi\); primary comparison for C and primary gauge pieces; quotient characters; absorb pivot rows before evaluating kernel atoms | The remaining secondary page form and its actual gauge phase; it is a row for some D completion |
| `bd_exact`: \(\gamma(\hat b,\hat b)-J(u,y,\pi)-\gamma(t,C)\), pure-C reference carry | The same carries, pure-C formulas and target pairings; direct A=0 formulas in high degrees | An exact marked residue can need information the page form discards |
| `a_step`, `QDa`, `fsharp`, `pot` | Remove first primitive using \(K_D\); compare C changes by the integral carry; evaluate a curvature potential before its boundary | `fsharp` is calibrated secondary data, and the potential includes R1–R3 |
| `a_step`, `Ysrc`, `QDv`, `ypp` | \(U=\Lambda u_R\) literally; primary comparisons for the closed binary inputs; combine comparisons for the whole source | Non-primary reference pieces cannot be omitted from Y's equation |
| `a_step`, `RC`, `page`, `thmC` | Native primary leading B row; special k=4 representative; simultaneous absorption in a truncated lower presentation; functional evaluation of the remaining secondary class | A-over-C is secondary in general; the nonclosed B in `thmC` is not a Bockstein input |
| `a_step`, `residue`: \(\Pi\Psi_{\rm rel}\), cylinder, cross phase and bounded diagonal tail | (A8) can remove independent tails; otherwise use (A14), carries for U/Y/W where proved and the existing unary degree-six phase | A two-primary tertiary residue remains when absorption fails; no replacement by an arbitrary native primitive |
| `p3_power`: generic projected cube or nineteen-term \(P^1\) | Implemented in GAP before the worker: native (A11) or (A12), followed by the integral lift solve | k=6 needs a supplied resolution through degree ten; shorter resolutions retain the bar path |
| Shared `combination`, `reference`, `pure_c` | Replace diagonal \(h^D(b,b)\) by its carry expression; use native coherent polarizations and (A6) for new pure-C markings | Cross terms \(h^D(b_i,b_j)\) for distinct inputs and \(N_C/2\) have real carries; do not infer them from diagonal classes |
| Shared `prim`, `bin`, `integer`, `rational`, `rational_ds` | (A3), quotient-target pairings and coefficient reduction with proved denominator bounds | Generic externally supplied cochains still require their normalization checks |

For low-degree cross terms one can evaluate the requisite fixed binary
cup comparisons and their polarizations directly, rather than normalize
an entire bar homotopy. Bilinear comparisons alone do not specify the
coordinate gauge of a rational upper product; its integral carry also
has to be transported or symbolically cancelled.

## 10. Model fallback and universal contractors

The following sites remain reachable when a light part is unavailable,
when light or native-relation options are disabled, or through explicit
model calls. They must not be omitted from a transport inventory.

| Site | Reduction strategy | Required restriction |
| --- | --- | --- |
| `TransferredModel.phi`: \(\Phi_\ell=\Lambda w_\ell-HN_\ell\) | Substitute comparisons (A2) for primary and Bockstein sources; retain triangular coordinate carries for D | Remove a higher H only after comparing its whole calibrated source |
| `kappa` / `nonlinear`: \(\delta_Rw+\Pi N(\Phi w)\) | Native primary maps, rational-potential boundary on R, target characters; propagate constructed legality certificates | A zero projection is not proof a bar cochain is identically zero |
| `bar_product`: \(\alpha,\beta,\gamma\) | Use existing A=0 and pure-C branches, finite low-degree formulas, native primary rows or light residues when only a relation is needed | General nonzero-A products require pair coherence and the fixed upper phase |
| `reflect`, `act`, `divide_left`: \(\Pi r,Hr\) and embedded gauge boundaries | Replace relation-only requests by light residues; use source comparisons and evaluate only requested layers | A nonflat gauge's curvature is not its full reflected boundary action |
| Off-legal gauge formulas and right-interval prisms | Expand the fixed low-degree prism into its interval-cut words once; cancel normalized or repeated-vertex terms before requesting any input values | Preserve the specified off-legal source and the phase's carries; a legal-locus formula cannot replace it |
| Degree-six J and gamma | Legal-locus formulas; direct A=0 formulas; unary diagonal when only powers are needed | The guarded product of two nonzero A states is not bypassed by a unary formula |
| Three-local `gamma`, `curvature_term`, `boundary` | Native cups/cyclic diagonals and their explicit cross-effect and coboundary primitives | Preserve divisions by three and the model's coordinate change; cyclic stage 3 is needed for the full k=6 product |
| Comparison setup, g construction, normalization \(h'=qhq\,\partial\,qhq\), and cochain zero tests | Avoid setup for fully native rows; combine requested chains first; dual scalar evaluation; modular collection; known-by-construction branch flags | Keep integral retraction verification, group actions and resource limits; no zero inferred from a budget refusal |
| Complete-bar reference engines | Restrict them to the existing bounded reference tests | They are outside the ordinary `koFull` path and are deliberately independent comparison oracles |

The sources for these sites are [extension_transfer.py](../python/extension_transfer.py),
[extension_transfer.gi](../gap/extension_transfer.gi),
[extension_native_upper.py](../python/extension_native_upper.py), and
[extension_three_local.py](../python/extension_three_local.py).

The universal costs need a separate treatment:

| Universal component | Exact reduction to pursue | What would be invalid |
| --- | --- | --- |
| Theta on a fixed simplex; chi/zeta/interval-cut words | Cache repeated binary face patterns; compile the lowest-degree words and known prism degeneracies | Treating Theta as an arbitrary primitive on R, or dropping its separate half-lifts |
| V1 and degree-one pair sources | Use the explicit interval-intersection formulas in [universal value growth](universal_value_growth.md); aggregate period-four source sums by residue counts instead of enumerating unit edges | Reducing the complete V1 value or pair-source labels modulo four; the values have affine growth |
| V2, its pair source and R2 prisms | Use the fixed reduced theta keys, share contractions and evaluate only the requested functional; derive residue-class summation formulas for the remaining label sums | Assuming the primitive is periodic merely because its theta source is periodic |
| V3, R3 and the degree-six pair source | The same fixed-source optimizations; reuse unary sources for diagonal powers; retain the zero-fiber short circuits | Dropping two nested V3 evaluations from an identity proved only before the actual universal section is applied |
| Binary legal-beta pair primitives | Zero when a relevant A fiber is zero; on a diagonal derive the actual specialized source, then contract or compile that source | Inferring off-diagonal cross terms from their diagonal values |
| A=0 higher pair contractions and prisms | Existing direct A=0/pure-C formulas; the explicit `closed_prism` reduction; compile further finite KB/interval patterns | Using legal closed-B formulas on arbitrary gauge inputs whose B is not closed |

The existing universal-value store and bundled tables already save repeated
values. They do not bound the number of new keys. Summation or symbolic
source simplification is the route to eliminating growth on a cold
computation; cache reuse alone does not prove that reduction.

## 11. Remaining implementation opportunities

The primary defining systems, integral A gauge, exact B carry reductions,
pure-C marked class, native three-primary rows and binary homotopy
collection are implemented. Further reductions have these dependencies:

1. Implement target-character projection and carry-aware coefficient
   reduction beyond binary H, for reuse by pages and every light task.
2. Use the binary diagonal comparisons for coherent new pure-C markings
   as in (A7), retaining their integral carries and subsequent references.
3. Extend pivot absorption to A rows under (A8), tracking all dependent
   rows and whether an explicit old-marking isomorphism is required.
4. Address the remaining calibrated secondary comparison first in the
   degree-zero/even-input sectors. Optimize the fixed universal sources
   independently of that comparison work. General two-primary A-over-D
   still needs the tertiary formulas when the presentation does not
   absorb its tail.

These priorities describe mathematical dependencies, not measured speedups.
Native tensors and combined comparison chains can themselves be large;
their cost must be measured on the resolutions and target quotients that
actually occur.
