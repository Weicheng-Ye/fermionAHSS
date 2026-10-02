# Reducing transport in the light extension relations

This note derives reductions from the package's own light residues and
comparison identities. It does not replace the calibrated secondary or
tertiary operations by a primitive chosen on the target resolution.
The implementation selects the primary comparison, the integral A-gauge
identity, the Bockstein carry with its exact marking correction, and the
pure-C class formula in
[`extension_light.py`](../python/extension_light.py). Its primary tensors
are constructed in
[`extension_primary_transport.gi`](../gap/extension_primary_transport.gi).
The additional absorption and target-character constructions below specify
further possible reductions; they are not general runtime shortcuts.
Finite executable checks are in
[`test_light_transport_reduction.py`](../python/test_light_transport_reduction.py).
The proofs, rather than those checks, specify their general scope.

The first application is spacetime dimension 3+1, package degree `k=4`:

\[
 (A,B,C,D)\in C^1(\mathbf Z_s)\times C^2(\mathbf F_2)
 \times C^3(\mathbf F_2)\times C^5(\mathbf Z_s).
\]

The Bockstein, diagonal-comparison and pairing arguments have no restriction
to this degree. The special description of the torsion A layer does.
The [transport reduction audit](transport-reduction-audit.md) extends these
arguments to the other light tasks, page operations and model fallbacks,
and distinguishes complete eliminations from reductions of their cost.

## 1. Conventions and what must be preserved

Let \(\mathcal B\) denote the normalized comparison complex, either the
group bar or the cell comparison. Write

\[
 \Lambda=f^*,\qquad \Pi=g^*,\qquad H=(h')^*,\qquad
 1-\Lambda\Pi=\delta H+H\delta,
 \quad\Pi\Lambda=1,\quad H\Lambda=0,\quad\Pi H=0.
 \tag{1}
\]

The maps are integral and used with the appropriate local system; their
reductions give the binary maps. All canonical lifts and signed
coboundaries use the same first-vertex frame. A tilde is a pointwise
integer lift, and \(\rho\) is reduction modulo two. Native cochains are
row vectors and \(\delta_s^R v=vM_n\) in degree \(n\).

These are the **extension** comparison identities. The page engine's raw
group-bar comparison does not assume \(\Pi\Lambda=1\) or \(H\Lambda=0\)
([conventions](conventions.md#5-whole-formula-transfer-and-defining-systems)).
In particular, (7) and the changes of marking below cannot be applied to
that comparison using the side conditions in (1).

There are three different required outputs:

1. an exact integer cochain in a fixed marking;
2. its cohomology class in that marking;
3. its extension residue in \(H_{\rm low}/mH_{\rm low}\), possibly after an
   explicitly recorded change of lower marking.

A simplification proved for the third output does not automatically supply
the first. In particular, a change of a lower generator by a D element has
to be propagated to every relation that refers to that generator. The
current [light frame](extensions.md#light-rows) already distinguishes
these precisions.

## 2. Bockstein primitives from carries instead of homotopies

Let \(b_R\in Z^n(R;\mathbf F_2)\). Define

\[
 b=\rho\Lambda_s\widetilde b_R,\qquad
 q_R=\frac{\delta_s^R\widetilde b_R}{2},\qquad
 \kappa_b=\frac{\Lambda_s\widetilde b_R-\widetilde b}{2}.
 \tag{2}
\]

Both divisions are integral. The first is integral because \(b_R\) is a
binary cocycle; the second because the two lifts have the same reduction.
Applying the signed coboundary gives the **cochain identity**

\[
 \frac{\delta_s\widetilde b}{2}
 =\Lambda_s q_R-\delta_s\kappa_b.
 \tag{3}
\]

On normalized ordered cochains the diagonal primary polarization is

\[
 x=h^D(b,b)=(\operatorname{Sq}^1+s)b
   =\rho\frac{\delta_s\widetilde b}{2}.
\]

Consequently, with \(\bar\kappa_b=\rho\kappa_b\),

\[
 x=\Lambda_2\rho q_R+\delta\bar\kappa_b.
 \tag{4}
\]

Suppose a B-relation gauge needs
\(\delta\pi=x+\Lambda_2 c_R+P\), where \(c_R\) is a closed reference
C cochain and \(P\) is a closed D- or Tau-gauge source. Solve only

\[
 \delta_R r=\rho q_R+c_R+\Pi_2P.
\]

Then a defining cochain is

\[
 \boxed{\pi=\Lambda_2r+\bar\kappa_b+H_2P.}
 \tag{5}
\]

Indeed, (4) and \(\delta H_2P=P+\Lambda_2\Pi_2P\) give the requested
equation, with both copies of \(\Lambda_2\Pi_2P\) cancelling. This
requires no projection of \(x\), no homotopy of \(x\), and no homotopy
of the lifted C reference. For the pure gauge \(P=0\), it requires
**no comparison homotopy at all**.

At `k=4`, \(n=2\): the carry is evaluated on triangles, rather than
building homotopy chains for the degree-three source \(x\). Its integer
lift needs \(\Lambda_s\widetilde b_R\bmod4\), not merely its parity.

### Relation to the current marking

Let the current primitive be
\(\pi_{\rm old}=\Lambda_2r_{\rm old}+H_2(x+\Lambda_2c_R+P)\).
Put \(r=r_{\rm old}+\Pi_2\bar\kappa_b\). Equation (1) gives

\[
 \pi_{\rm old}=\pi+\delta H_2\bar\kappa_b.
 \tag{6}
\]

Thus the change is an explicit coboundary, not an arbitrary closed
cochain. A class calculation can use the corresponding gauge equivalence;
an exact marked calculation must carry its induced D-coordinate change.
If that witness is needed, its homotopy acts on a degree-\(n\) carry,
one degree below the source used by the old primitive.

There is a still simpler identity for an A-relation gauge. If
\(A=\Lambda_s A_R\) and \(\delta_s^R u=mA_R\), then

\[
 \boxed{P(mA;u)=\Lambda_su+H_s(m\Lambda_sA_R)=\Lambda_su.}
 \tag{7}
\]

This is literal equality by \(H_s\Lambda_s=0\). It removes that
homotopy without changing any marking, in every supported degree.

## 3. Replace a primary-source homotopy by binary diagonal comparisons

Ordinary native cup-i products alone do not justify substituting native
defining cochains into a secondary formula. The missing data can be
specified explicitly, at arity two, for the primary source \(Q_D\).

Let \(\smile_i^R\) be native binary higher cups and
\(\smile_i\) the fixed cups on \(\mathcal B\). Introduce bilinear
comparison operators \(K_i\) of cochain degree \(-i-1\), with
\(K_{-1}=0\), satisfying

\[
\begin{aligned}
 \delta K_i(a,b)+K_i(\delta a,b)+K_i(a,\delta b)
 ={}&\Lambda_2(a\smile_i^R b)
      +\Lambda_2a\smile_i\Lambda_2b\\
    &+K_{i-1}(a,b)+K_{i-1}(b,a).
\end{aligned}
\tag{8}
\]

All terms in this section are binary. For a closed degree-\(n\) cochain
\(b_R\), put \(S_R^1b_R=b_R\smile_{n-1}^Rb_R\) and

\[
 D_Rb_R=b_R\smile_{n-2}^Rb_R+
          \omega_R\smile^Rb_R+s_R\smile^RS_R^1b_R.
\]

The comparison cochain is

\[
\boxed{\begin{aligned}
 K_D(b_R)={}&K_{n-2}(b_R,b_R)+K_0(\omega_R,b_R)\\
           &+K_0(s_R,S_R^1b_R)
             +\Lambda_2s_R\smile K_{n-1}(b_R,b_R).
\end{aligned}}
\tag{9}
\]

For `k=4` B atoms, only \(K_0,K_1\) occur. Negative-index terms are
zero. The twists on \(\mathcal B\) here are \(\Lambda_2s_R\) and
\(\Lambda_2\omega_R\), as in the light worker. Equation (8), closedness,
and the Leibniz rule give

\[
 \boxed{\delta K_D(b_R)=Q_D(\Lambda_2b_R)+\Lambda_2D_Rb_R.}
 \tag{10}
\]

To see the cancellation, the two \(K_{i-1}(b_R,b_R)\) terms cancel
in the square comparison. The coboundary of the last term of (9)
compares the inner Sq1; the preceding term compares the outer cup with s.
Their common \(\Lambda_2s_R\smile\Lambda_2S_R^1b_R\) cancels.

One may therefore solve \(\delta_Rc_R=D_Rb_R\) natively and use

\[
 \boxed{C=\Lambda_2c_R+K_D(b_R).}
 \tag{11}
\]

It satisfies \(\delta C=Q_D(\Lambda_2b_R)\) literally. This removes
both \(\Pi Q_D(\Lambda_2b_R)\) and \(H Q_D(\Lambda_2b_R)\) from
the construction of a B atom. It does not remove its secondary
flat-admissibility condition or change the calibrated phase used to test it.

Here \(D_R\) denotes the displayed **cup-defined** cochain. The runtime's
`backend.nativePrimary("D",...)` instead uses
\(e_R=\rho(\delta^R\widetilde b_R/2)\) for Sq1. On R,
\(e_R\) and \(S_R^1b_R\) have the same class, but need not be equal as
cochains. To use that runtime cochain in (11), solve the known primary
coboundary

\[
 \delta_R t_R=S_R^1b_R+e_R,\qquad
 K_D^{\rm native}=K_D+\Lambda_2(s_R\smile^Rt_R).
 \tag{11a}
\]

Then \(\delta K_D^{\rm native}=Q_D(\Lambda_2b_R)
+\Lambda_2D_R^{\rm native}b_R\). Use this corrected comparison in
(11) and (13). This finite linear solve compares two representatives of a
known primary class; it does not choose a secondary operation. Its choice
is part of the defining system and any subsequent change of marking.

### Constructing the comparison without the normalized bar homotopy

Equation (8) is not an instruction to solve for an arbitrary primitive on
R. Let \(\Delta_i^R,\Delta_i^{\mathcal B}\) be the chain diagonals and
\(T\) swap the two tensor factors. Write

\[
 L_i=\Delta_i^Rf+(f\otimes f)\Delta_i^{\mathcal B}
           +(1+T)\mathcal K_{i-1}.
\]

Construct \(\mathcal K_i:\mathcal B\to(R\otimes R)[i+1]\) on free
basis simplices, in increasing dimension, by

\[
 \mathcal K_i(\sigma)=S\bigl(L_i(\sigma)+
                         \mathcal K_i(\partial\sigma)\bigr),
 \tag{12}
\]

where \(S\) is the tensor contraction of R, and extend equivariantly.
The argument of S is a cycle: the higher-diagonal identities and
\((1+T)^2=0\) show that \(L_i\) is a cycle in the Hom complex, and
the induction cancels its boundary. At degree zero its augmentation
vanishes because both diagonals preserve augmentation. Hence (12)
has boundary \(L_i+\mathcal K_i\partial\), as required. Pairing with
\(a\otimes b\) gives (8).

These are binary tensors on R. They use f, the standard simplex diagonals
and R's contraction; they require neither g nor the normalized
\(H=(qhq\,\partial\,qhq)^*\). They can be cached by the required
bidegrees and reduced modulo two during their construction. Their cost
still needs measurement; their existence is not a bound on tensor size.

For an exact comparison with the old C marking, put
\(k_R=\Pi_2K_D(b_R)\) and \(c_R=c_{\rm old}+k_R\). Then

\[
 \delta_Rc_R=D_Rb_R,\qquad
 C_{\rm old}=C+\delta H_2K_D(b_R).
 \tag{13}
\]

This follows by applying (1) to \(K_D\). The change of native defining
vector is part of the formula. Keeping \(c_{\rm old}\) unchanged in
(11) would generally violate its defining equation.

For an A generator the same construction applies to its closed binary
reduction in degree \(k-3\), and removes the first defining-system
homotopy. Eliminating the next source \(f^\sharp(A,B)\) also requires
the calibrated higher comparison data; (8) alone does not supply them.

### The integer D carry of a C coboundary change

In degrees `k=3,4,5,6`, fix legal A and B and let
\(C'=C+\delta\eta\) in binary cochains. Define

\[
 v_\eta=E(\eta)+h^2(C,\delta\eta),\qquad
 I_\eta=\frac{\widetilde{E(C')}-\widetilde{E(C)}
                     -\delta_s\widetilde v_\eta}{2}.
 \tag{13a}
\]

Here the sum defining \(v_\eta\) is binary, and \(h^2\) is the
polarization of E, with no sign term. The identities
\(\delta E(\eta)=E(\delta\eta)\) and
\(\delta h^2(C,\delta\eta)=E(C')+E(C)+E(\delta\eta)\)
show that the numerator is even, so \(I_\eta\) is integral.
The [legal A,B potential](../python/stacking_model/closed_ab_upper.py)
in degrees four through six depends on C only through \(\widetilde{E(C)}/2\)
and through \(t=\delta C+f^\sharp(A,B)\), which is unchanged.
The degree-three potential has the same dependence. Consequently,

\[
 J_k(A,B,C')-J_k(A,B,C)=\delta_s I_\eta,\qquad
 D'=D-I_\eta.
 \tag{13b}
\]

Thus the new D coordinate preserves flatness exactly. This supplies an
explicit integer carry for changes such as (13), including nonclosed C.
To make this a change of the entire stacking model, conjugate its product
and gauge maps by the same triangular coordinate change. Flatness alone
does not identify two marked products, and simply discarding \(I_\eta\)
is not justified.

## 4. Exact pure-C relation classes with one projected source

The ordinary C-over-D extension residue in \(D_k/2D_k\) is already
computed natively by the primary-operation shortcut. This section concerns
the more precise class in \(D_k\) for a fixed C lift, which an upper
relation can require; it does not introduce native evaluation of the
standalone C-over-D relation.

Let \(C=\Lambda_2c_R\) be a closed C representative, in package degree
\(1\le k\le6\), and set

\[
 e_R=\Pi_s\widetilde{E(C)},\qquad
 \delta_s^R D_c=-\tfrac12\delta_s^Re_R.
\]

The right side is the original pure-C curvature, so this uses the same
native D completion if the same integer solver is used. Let
\(Y_R\) be any binary cocycle representing \(\operatorname{Sq}^1[c_R]\),
for example the reduction of the ordinary integral Bockstein on R.
The [pure-C product identity](extension_cup_i_formulas.md#5-c-over-d) gives
the exact **cohomology class** of its marked double:

\[
 \boxed{[2\widetilde c]=
  \left[2D_c+e_R+\frac{\delta_s^R\widetilde Y_R}{2}\right]
  \quad\text{in }D_k.}
 \tag{14}
\]

This is more than a parity formula once \(e_R,D_c\) use the fixed
comparison. Indeed the bar correction is
\(\widetilde{E(C)}+\delta_s\widetilde{Q^1C}/2\) plus an integral
coboundary. Its second summand represents
\(\beta_s\operatorname{Sq}^1[c_R]\), which can be computed on R.
The complete expression is closed by the equation for \(D_c\).

For \(k=1\), the closed degree-zero input has Sq1 zero and its square
correction is exactly \(\widetilde{E(C)}\). For \(k=2\), the closed
pure-C phase is \(-\widetilde{C\smile C}/2\), so its square correction
is \(\widetilde{E(C)}-\delta_s\widetilde{C\smile C}/2\).
Changing that minus to the plus in (14) is an integral coboundary.
Thus the same class formula applies in both low degrees; see
[low-degree stacking](low_degree_stacking.md).

Only \(e_R\) needs bar projection, in degree \(k+1\). No degree-\(k+2\)
comparison chain and no bar evaluation of the full pure-C product are
needed. At `k=4` the maximum comparison degree for this calculation drops
from six to five. Replacing \(e_R\) as well by a native representative
without its comparison carry gives only the residue modulo \(2D_k\),
as in the existing primary shortcut; that is a different precision.

## 5. The A layer in spacetime dimension 3+1

For a group and an integral sign system, the torsion in
\(H^1(G;\mathbf Z_s)\) is zero when s is trivial and is exactly one
\(\mathbf Z/2\) when s is nontrivial. Here is a direct proof. If a crossed
homomorphism a represents a torsion class, then
\(ma(g)=(\chi(g)-1)c\), with \(\chi=(-1)^s\). Thus a vanishes on
\(\ker\chi\). It has one integer value on the other coset, and
principal crossed homomorphisms change that value by an even integer.
Conversely the sign cocycle supplies its nonzero order-two class.

Passing to the surviving A subgroup preserves this conclusion. Its free
part has no extension relation. Thus at `k=4` there is **at most one
torsion A relation**, its order is two, and its cochain can be chosen as

\[
 A=s,\qquad U=-1,\qquad\delta_sU=2s.
 \tag{15}
\]

These equalities refer to the group-bar coordinate; changing a supplied
native representative to it must be recorded. A native representative of
the same torsion class is \(-\delta_s^R\widetilde1/2\).

The leading B part of the square is derived without an upper phase:
\(h^D(s,s)=s^2+s^2=0\), whereas the gauge removing \(2s\) contributes
\(Q_D(\rho U)=Q_D(1)=\omega\). Therefore

\[
 \boxed{(2\widetilde A)_B=[\omega]\in B_4.}
 \tag{16}
\]

This does not determine its C or D carry. In the \(\omega=0\) case the
existing A-over-C identity is \([h^D(B_A,B_A)]\), with
\(\delta B_A=s^3\). The following independent presentation argument
specifies when the remaining D carry need not be evaluated.

## 6. Absorb an upper D carry only when a lower character detects the row

Let H be the already computed lower stacking group, D its bottom subgroup,
and \(K=H/D\). Suppose **\(2D=0\)**. A torsion A row has the form

\[
 2a=t_0+d,\qquad d\in D,
\]

where the class \(\bar t_0\in K\) is known. If
\(\bar t_0\notin2K\), there is a character
\(\ell:K\to\mathbf F_2\) with \(\ell(\bar t_0)=1\).
For every d the map

\[
 \varphi_d(h)=h+\ell(\bar h)d
 \tag{17}
\]

is a homomorphism, fixes D pointwise, induces the identity on K, and
satisfies \(\varphi_d^2=1\). It sends \(t_0+d\) to \(t_0\). Hence

\[
 \boxed{(H\oplus\mathbf Za)/(2a-t_0-d)
       \cong(H\oplus\mathbf Za)/(2a-t_0)}
 \tag{18}
\]

as filtered groups, with the identity on each associated-graded layer.
No candidate formula for the unknown carry is needed in this proof.

This is a **change of lower marking**, not equality of two extensions
with the identity on the entire marked H. The cheaper row is exact in the
new marking and must be labelled accordingly. If d is not evaluated,
(17) proves that such an isomorphism exists; it does not compute the
isomorphism to the old marked generators or a gauge witness for them.
At `k=4` there is no second
torsion A row to which an unrecorded change could propagate. In higher
degrees every other row changes by the same \(\varphi_d\); independent
absorption of several rows is not justified.

The character is found from H's presentation after quotienting by D and
reducing modulo two. In particular, a nonzero leading B class \([\omega]\)
always makes the test pass. With zero leading B part, the C part must be
tested in \(K/2K\), not just for being a nonzero C vector.

Both hypotheses are essential:

- If \(H=D=\mathbf Z/2\), the rows \(2a=0\) and \(2a=d\) give
  \((\mathbf Z/2)^2\) and \(\mathbf Z/4\). There is no detecting character
  on \(K=0\).
- If \(H=\mathbf Z/2\{b\}\oplus\mathbf Z/4\{d\}\), the rows
  \(2a=b\) and \(2a=b+d\) give \((\mathbf Z/4)^2\) and
  \(\mathbf Z/2\oplus\mathbf Z/8\). Here \(2D\ne0\).

For a relation of order m a sufficient generalization is a homomorphism
\(\ell:K\to\mathbf Z/m\) detecting the known row by a unit and a bottom
subgroup killed by m. Rescale the character so that its value on the row
is one. The same construction is invertible, with inverse
\(h\mapsto h-\ell(\bar h)d\). This remains a statement about a coherent
change of all markings.

## 7. Project only the characters of the required residue

A light D residue has the form

\[
 r_R=\Pi_sF+\delta_s^R\Pi_sP,
 \qquad |F|=k+1,\quad |P|=k,
 \tag{19}
\]

with known native completion terms added when required. After assembly it
is integral, or has a specified integral coefficient reduction. The
presentation needs its image in \(D/(D\cap mH)\), not every coordinate
of \(C^{k+1}(R)\).

For `k=4`, m is two. Choose binary linear functionals on the space of
reduced integral cocycles that annihilate incoming classes and
\(D\cap2H\), using the existing marked presentation. Extend them to the
ambient binary cochain space while keeping them zero on coboundaries.
This is finite-dimensional linear algebra. Integral cohomology modulo
two injects into binary cohomology by the coefficient exact sequence, so
these functionals distinguish precisely the required quotient. If
\(\ell\) is one such column, then \(M_k\ell=0\pmod2\).

Lift \(\ell\) to an integer column. The exact evaluation identity is

\[
 \boxed{r_R\ell=
    \langle F,g_{k+1}\ell\rangle+
    \langle P,g_k(M_k\ell)\rangle.}
 \tag{20}
\]

It is just adjointness and the chain-map equation, with the signed
coefficient chain complex understood. **Combine the native chains before
applying g.** Computing every \(g(e_j)\) and cancelling afterwards loses
the benefit. A zero boundary \(M_k\ell\) eliminates the P term exactly.
An even boundary can be divided by two with the factor applied to P.

The recursive construction of g must retain group-ring coefficients and
actions until it has applied the bar cone. In particular,
\(M_k\ell=0\) in the signed coefficient complex does not imply
\(g_{k+1}\ell=0\). The bar cone is not equivariant, so coning a boundary
after signed augmentation is not an implementation of (20).

Thus a rank-r elementary-two target needs at most r tests in each of the
two degrees in (20), rather than full projections in those degrees.
This bounds the number of pairing targets, not their support sizes.
Exact marked integer rows may need additional information; (20) does not
justify returning them from parity data.

## 8. Carry-aware finite coefficient transport

Suppose the scalar formula has a proved denominator bound q and its
assembled value \(N/q\) is integral. To retain its residue modulo m,

\[
 \boxed{N/q\pmod m=((N\bmod mq)/q)\pmod m.}
 \tag{21}
\]

The representative \(N\bmod mq\) is divisible by q. Reducing N modulo
q instead would erase the integral carry being measured. For a rational
phase known only modulo one, modulus q is sufficient; that is a different
request.

Every *linear* stage of f, g and the normalized homotopy commutes with
coefficient reduction. Binary source homotopies may consequently collect
and cancel coefficients modulo two before each recursive expansion.
Equation (2)'s binary carry needs modulus four in its integer lift.
For (20), use the modulus dictated by (21) until the exact divisions have
been performed. Choices such as the cell comparison's integral right
inverse must remain the fixed integral choices before reduction; solving
an unrelated modular system would choose a different comparison.

This allows smaller chains without changing the calibrated scalar
formulas. It does not discard the integral retraction proof, local-system
signs, normalization requirements or marked-coordinate conversions.

## 9. A degree-four evaluation order

The reductions suggest the following order for the two-primary part:

1. Keep the existing target-layer and primary-operation shortcuts.
2. Build B defining data using native \(D_R\) and (9)–(11); use (5) for
   the relation gauge. Pure gauges need no H; other gauges retain only
   the H of their D or Tau source.
3. Use (14) when an exact marked C relation class is needed. Select the
   final D-residue characters before evaluating (19), then use (20).
4. For the possible single torsion A generator, use (7), or the recorded
   change to (15). Determine its lower relation first. Apply (18) only
   when its hypotheses hold and the change of marking is recorded.
5. Evaluate the fixed upper residue in the remaining A cases. An upper
   row with zero class in \(K/2K\) can contain essential tertiary data.

This removes identifiable transports rather than assuming all higher
operations are primary. A fully native secondary or tertiary calculation
still needs calibrated higher coherence and its comparison with the fixed
formula family. Neither an arbitrary primitive on R nor an uncorrected
substitution of native cup-i products supplies that comparison.

## 10. Implemented defining systems

Primary sources use the higher-diagonal convention in (9), including
\(b_R\smile_{n-1}^Rb_R\) for Sq1. Thus (11a)'s conversion to the
different `nativePrimary` representative is not needed. The worker returns
\(D_Rb_R\) as the defining right-hand side; GAP solves
\(\delta c_R=D_Rb_R\), and every subsequent task constructs
\(C=\Lambda c_R+K_D(b_R)\). This removes both \(\Pi Q_D(b)\)
and \(HQ_D(b)\). It applies to B atoms, the first B cochain of an A
system, closed integral inputs to D/Tau gauges, and the A-page `ypp` input.
The native tensor contraction checks its cycle and filler equations. The
simplex cuts come from the existing interval-cut engine.

These are chosen defining systems, not literal copies of the generic
\(P\) primitive. Each atom's curvature is evaluated on the chosen C,
its D completion is solved from that curvature, and every product,
reference, re-marking and later row uses the same atom. Corrections by
closed native cochains still add their \(\Lambda\)-lifts. The light
witness records `primaryDefiningSystem="Lambda r + K_D"`. The calibrated
stacking formulas and the comparison maps themselves are unchanged.

For an A relation, write \(U=\Lambda_su_R\) using (7). Its binary
reduction is closed even though U need not be an integral cocycle.
For order two the Y source is
\(Q_D(\rho U)+h^D(a,a)+B_0\), where \(a=\rho\Lambda_s A_R\)
and \(B_0=\Lambda_2b_{0,R}\). Apply (9) to \(\rho u_R\) and
(2)–(4) to \(\rho A_R\). If their comparison terms are \(K_U\)
and \(\rho\kappa_a\), solve

\[
 \delta y_R=D_R(\rho u_R)+\rho\beta_s^R(\rho A_R)+b_{0,R},
 \qquad Y=\Lambda_2y_R+K_U+\rho\kappa_a.
\]

For order at least four the polarization term is zero and its two terms
are omitted. This also avoids normalized bar H. The later secondary
source \(R_C\) retains its generic primitive. Tau corrections use their
matched defining systems as before.

The B-relation gauge retains the original primitive exactly, using (6):
its source projection uses
\(\Pi_2 h^D(b,b)=\rho q_R+\delta_R\Pi_2\rho\kappa_b\),
so the carry is projected one degree below the source. Native pure-C
references add their vectors directly.
Its carry homotopy is \(H\rho\kappa_b\), one degree lower. A remaining
D-gauge source \(P=Q_D(\Lambda y_0)\) is evaluated exactly by
\(HP=K_D-\Lambda\Pi K_D-\delta HK_D\), again with H one degree lower;
the Tau source retains its calibrated homotopy.
The integral A primitive uses (7) literally. Exact C rows use (14), keeping
the same D completion and changing their square cochain only by an
integral coboundary. Binary homotopy chains are collected modulo two at
each linear stage of the same normalized contraction. Native-lift and
homotopy-image pairings use the extension comparison's proved side
identities; requested g chains still receive their integral retraction
and resource checks.

At the prime three, the degree-five row uses the native signed cube.
The degree-six row uses the native cyclic diagonal when the supplied
resolution reaches degree ten, including the spare contraction degree;
otherwise it retains the bar reduced-power formula. Set
`FERMIONAHSS_LIGHT_TRANSPORT_REDUCTION=0` to restore the generic defining
systems and transport evaluation paths for comparisons.
