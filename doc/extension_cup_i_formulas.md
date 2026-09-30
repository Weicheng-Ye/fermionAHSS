# Cup-i formulas for the extension relations

`koFull` measures the relation \(m\widetilde q=t\) of a layer generator in
the stacking model on the normalized group bar. The relation engine
(`KOAHSS_ExtensionRelationEngine` in
[gap/extension_relations.gi](../gap/extension_relations.gi)) evaluates
every product, gauge boundary and reduction of the model through the
transport of [python/extension_transfer.py](../python/extension_transfer.py):
a state on the resolution R is embedded into bar cochains, the bar product
is formed there, and the result is reflected back to R
([resolution-extensions.md](resolution-extensions.md),
[transfer.md](transfer.md)). Below, "the engine" means this measurement in
the model.

This note takes each pair of layers X over Y and expresses the class of the
Y component of \(m\widetilde q\) on the cochains of R itself. That class is
what the engine records as Y coordinates when it measures the relation
through Y. The formulas use cup-\(i\) products and standard operations: the
reduction \(\rho\), integral lifts, the twisted Bockstein \(\beta_s\),
Steenrod squares written \(x\smile_ix\), and the twists \(s\) and
\(\omega\). The note also says when such an expression exists and why. The
formulas describe the classes the engine computes. For a two-primary
relation whose target layer lies right below the generator's, `koFull`
evaluates the formula of Sections 3–5 on R instead of measuring the
relation in the model ([extensions.md](extensions.md), "Primary-operation
rows"); the other relations are measured in the model.

The adjacent pairs A–B, B–C and C–D have primary answers, valid in every
degree where the pair occurs and for every twist (Sections 3–5). The
non-adjacent pairs A–C, B–D and A–D involve secondary or tertiary
operations. Section 6 records which parts of them are cup-\(i\) expressions
on R and which are not. Section 7 is a summary table.

## 1. Relations, target layers and recorded coordinates

The layers of package degree \(k\) have cochains in
\(C^{k-3}(R;\mathbf Z_s)\), \(C^{k-2}(R;\mathbf F_2)\),
\(C^{k-1}(R;\mathbf F_2)\) and \(C^{k+1}(R;\mathbf Z_s)\), for A, B, C and
D. Write \(X\) for the space modelled by R (BG for a group resolution).
Rows \(q>0\) are absent and row \(-3\) is zero. The only differentials
entering the four layers in the window are therefore \(\operatorname{Dbar}\),
\(D\), \(\operatorname{Dtilde}\), `Tau`, `Psi` and `T`, and the E6 layers are

\[
\begin{aligned}
B_k&=\frac{\ker D\cap\ker\operatorname{Psi}\ \subset H^{k-2}(X;\mathbf F_2)}
{\operatorname{Dbar}H^{k-4}(X;\mathbf Z_s)},\\
C_k&=\frac{\ker\bigl(\operatorname{Dtilde}:H^{k-1}(X;\mathbf F_2)\to H^{k+2}(X;\mathbf Z_s)\bigr)}
{D\,H^{k-3}(X;\mathbf F_2)+\operatorname{Tau}E_3^{k-4,0}},
\qquad E_3^{k-4,0}=\ker\operatorname{Dbar}\subset H^{k-4}(X;\mathbf Z_s),\\
D_k&=\frac{H^{k+1}(X;\mathbf Z_s)}{I_k},\qquad
I_k=\operatorname{Dtilde}H^{k-2}(X;\mathbf F_2)
+\operatorname{Psi}_{k-3}(\ker D)+T_{k-4}\bigl(E_5^{k-4,0}\bigr).
\end{aligned}
\]

\(A_k\subset H^{k-3}(X;\mathbf Z_s)\) is the subgroup on which
\(\operatorname{Dbar}\), `Tau` and `T` vanish. See the
[formula sheet](README.md), (S16), (S17) in
[secondary_operations.md](secondary_operations.md), and (A23) in
[all_cochain_differential.md](all_cochain_differential.md).

**The relation.** Let \(q\) be a generator of order \(m\) in a layer X. Let
\(\widetilde q\) be a flat lift of \(q\) whose leading component is the
marked E6 cocycle. The relation reads \(m\widetilde q=t\) in the lower group
\(H\), which is formed by the layers below X. Only the class of \(t\) in
\(\operatorname{Ext}(\mathbf Z/m,H)=H/mH\) is invariant. The relation is
measured through its target layer, the lowest layer that holds a lower
generator outside \(mH\) ([extensions.md](extensions.md), "Layers and
measured relations").

**How the engine measures it.** The engine forms the power by binary
exponentiation (`KOAHSS_ExtensionPower`) and then reduces it layer by
layer from X downward (`reduce`). At each layer it does three things:

- It solves the component as a combination of the marked lower cocycles
  plus a coboundary. The solve is over \(\mathbf Z\) with R's signed
  coboundary matrix in A and D, and modulo two in B and C.
- It divides on the left by the product of the lower lifts with those
  coefficients, and by the boundary state of the gauge made from the
  coboundary primitive.
- It continues with the quotient.

The rows of each solve are the marked cocycles and the coboundaries. They
do not contain the images of the incoming differentials. When a solve
fails:

- A layer-limited measurement falls back to the complete measurement.
- In the complete measurement, only the D layer proposes coordinates, from
  its E6 cell projection. `KOAHSS_ExtensionGaugeReduce` then searches staged
  gauges. It returns unresolved when a lower generator is free (order zero)
  or when there are more than 32 lower normal forms.

**The component class.** For a pair X over Y, let \(x_Y\) be the Y
component of the residual when the reduction reaches Y, taken as a class
in \(Y_k\). The recorded Y coordinates are the coordinates of \(x_Y\) in
the marked basis.

- For an **adjacent** pair (Y right below X), \(x_Y\) is the leading
  component once the layer X has been reduced. The generator \(q\) alone
  determines it modulo \(mY_k\): that is zero for a binary Y, and
  \(2D_k\) for \(Y=D\).
- For a **non-adjacent** pair, \(x_Y\) also depends on the coordinates
  already recorded in the layers between X and Y, and on the lifts of the
  generators of those layers.

## 2. What is needed on R

**Transport notation.** Write \(\Lambda=f^*\), \(\Pi=g^*\) and
\(H=(h')^*\) for the comparison of the normalized bar with R. They satisfy

\[
\Pi\Lambda=1,\qquad 1-\Lambda\Pi=\delta H+H\delta,\qquad
\Pi H=0,\quad H\Lambda=0,\quad H^2=0,
\]

and \(\rho\Lambda=\Lambda\rho\) ([transfer.md](transfer.md)). On cohomology
\(\Pi^*\) is the canonical isomorphism. The engine embeds a state by the
triangular \(\Phi(w)_\ell=\Lambda w_\ell-HN_\ell(\Phi(w)_{<\ell})\) and
reflects bar products back to R.

**Tools on R.** A formula on R in this note uses only the following:

1. **R's coboundary matrices.** \(M_s\) is the twisted matrix
   (`model.matrix(n,true)`, built from the sign character in
   [gap/hap.gi](../gap/hap.gi)) and \(M\) is the untwisted one. Cochains are
   row vectors, so \(\delta_sx=xM_s\).
2. **A mod-two cup-\(i\) structure on R.** This means equivariant maps
   \(D_i:R\to R\otimes R\) with
   \(\partial D_i+D_i\partial=(1+T)D_{i-1}\) modulo two, where \(D_0\) is a
   diagonal approximation. On cochains
   \[
   \delta(x\smile_iy)=\delta x\smile_iy+x\smile_i\delta y
   +x\smile_{i-1}y+y\smile_{i-1}x\pmod 2,
   \]
   and negative indices give zero. Two such structures are available:
   - HAP's recursively built diagonals, exposed as `data.cupMod2` in
     [gap/hap.gi](../gap/hap.gi) and as the backend `cup` in
     [gap/cochains.gi](../gap/cochains.gi);
   - the transported structure
     \(x\smile^\Pi_iy=\Pi(\Lambda x\smile_i\Lambda y)\).
3. **Bocksteins.** \(\beta_s[z]=[\delta_s\widetilde z/2]\), computed on R as
   \(\rho\bigl((\widetilde zM_s)/2\bigr)\). For a cocycle \(z\),
   \(\operatorname{Sq}^1z\) may be taken as \(z\smile_{n-1}z\) or as
   \(\rho(\widetilde zM/2)\); the two have the same class. The non-natural
   branch of `backend.primary("D")` in [gap/cochains.gi](../gap/cochains.gi)
   evaluates
   \(D(z)=z\smile_{n-2}z+\omega\smile z+s\smile\rho(\widetilde zM/2)\) this
   way. Under `koAHSSNaturalOperations` the backend routes \(D\) to the
   bar, so an evaluation on R calls the cup products directly.
4. **Linear solves and E6 coordinates.** Solves over \(\mathbf Z\) and
   \(\mathbf F_2\) use these matrices. The E6 cell projection
   `layers.<L>.cell.project` ([gap/pages.gi](../gap/pages.gi)) and
   `KOAHSS_ExtensionGroupCoordinates` turn a class into coordinates in the
   marked basis. Both work modulo exactly the incoming images.

**What does not depend on the cup-\(i\) structure.** The classes of
\(\operatorname{Sq}^jz=z\smile_{n-j}z\) of cocycles, of cup products of
cocycles, and of Bocksteins do not depend on the structure. The same holds
for the classes of \(D\), \(\operatorname{Dbar}\) and
\(\operatorname{Dtilde}\). Equivariant cup-\(i\) families on a free
resolution are unique up to equivariant homotopy, and \(\Pi^*\) is
canonical.

**What does depend on it, or is not available on R.**

- Cochain values themselves.
- Nested words of non-closed cochains. On R they differ from the bar
  composite by the \(H\) terms of \(1-\Lambda\Pi\).
- The interval-cut identities \(e(x)=x\smile_{n-1}x\) and \(x\smile_nx=x\)
  of (C6). On R these hold only up to coboundaries.
- The arity-four surjection words \(\zeta_{i,n}\) and the \(\chi\) words of
  [universal_helpers.md](universal_helpers.md).
- The universal sources \(V_n\), \(\Theta_m\) and the universal pair
  primitives.
- The prisms on \(X\times I\).
- The pointwise digits and carries (\(t=\rho((A-\widetilde a)/2)\), the
  splitting carries \(\lambda\), and the half lifts \(\mathsf h\) of
  non-closed words). These are functions of cochain values. They do not
  commute with chain maps.

**When a cup-\(i\) formula exists.** A component class can be written with
cup-\(i\) products on R exactly when it is a primary operation. The
operation is applied to the marked cocycle, or to a solution of a linear
equation on R such as a Bockstein preimage. In that case every non-primary
term of the engine's cochain is a coboundary, and the engine's reduction
ignores coboundaries.

When the class is the value of a secondary or tertiary operation, its
cochain needs a *natural* primitive of a relation among primary
operations. A primitive solved on R by linear algebra is ambiguous by every
cocycle of that degree. A natural one needs coherent structure beyond
binary cup-\(i\) products. The value is also fixed only relative to the
model's calibrated universal choices.

Throughout, write

\[
h^D(x,y)=x\smile_{r-1}y+\delta x\smile_ry
+s\smile\bigl(x\smile_ry+\delta x\smile_{r+1}y\bigr)
\qquad(|x|=|y|=r).
\]

This is the binary cross term of the stacking products
(`phase_eval.hD`; [ALL_COCHAIN_STACKING.md](ALL_COCHAIN_STACKING.md), Section 2). \(\alpha(A,A')=h^D(\rho A,\rho A')\), and
the legal \(\beta\) restricts to \(h^D(B,B')\) when
\(A=A'=0\). For a closed \(x\), \(h^D(x,x)\) represents
\((\operatorname{Sq}^1+s)[x]=\rho\beta_s[x]\).

## 3. A over B

**Setting.** An A generator \(a\) has marked cocycle
\(A\in Z^{k-3}(R;\mathbf Z_s)\) and order \(m=2^e\). The pair occurs in
package degrees \(k=4,5,6\). At \(k=3\) the group \(H^0(\mathbf Z_s)\) is
torsion-free, and below that there is no A layer.

- Solve \(\delta_su=mA\) over \(\mathbf Z\) on R. This is the A-layer solve
  of `reduce`, which has no marked rows.
- Put \(z=\rho u\). Because \(m\) is even, \(z\) is a mod-two cocycle, and
  \(\beta_s[z]=\tfrac m2[A]\).

Equivalently, \(z\) is any class with \(\beta_sz=\tfrac m2[A]\), or
\(z=r_{m\to2}\beta_{s,m}^{-1}[A]\).

**Cohomology formula.**

\[
x_B=D(z)=\operatorname{Sq}^2z+s\operatorname{Sq}^1z+\omega z\ \in B_k .
\tag{X1}
\]

For \(m=2\), (X1) is \(D\) applied to a Bockstein preimage of \([A]\).

| k | \(x_B\) |
| ---: | --- |
| 4 | \([\omega]\) (\(m=2\) and \(z=1\) are forced) |
| 5 | \(s z^2+\omega z\) |
| 6 | \(z^2+s\operatorname{Sq}^1z+\omega z\) |

At \(k=4\) the order two is forced. On connected \(X\) the torsion of
\(H^1(\mathbf Z_s)\) is \(\beta_sH^0(\mathbf F_2)=\{0,[\widetilde s]\}\),
which requires \([s]\neq0\). Then \(H^0(\mathbf Z_s)=0\), so there is no
indeterminacy.

When the target layer is B, every C and D generator lies in \(mH\) and
\(mB_k=0\). So \(H/mH\cong B_k\), and \(D(z)\) is the whole extension class.
When the target layer is C or D, (X1) is the B-graded part of the class.

**Cochain formula on R.**

\[
x_B^R=z\smile_{k-6}z+s\smile(z\smile_{k-5}z)+\omega\smile z
\ \in Z^{k-2}(R;\mathbf F_2),
\]

for any mod-two cup-\(i\) structure on R. Here
\(s\smile(z\smile_{k-5}z)\) may be replaced by
\(s\smile\rho(\widetilde zM/2)\). The coordinates are those of
`layers.B.cell.project` applied to the class of \(x_B^R\), in the marked
generators. Equivalently they are the unique \(c_j\) with

\[
x_B^R=\sum_jc_je_{B,j}+\sum_id_i\,D^R(\rho v_i)+\delta\pi\pmod 2,
\]

where \(v_i\) runs over representatives of \(H^{k-4}(R;\mathbf Z_s)\).

Nothing else is needed: no flat lift of \(A\), no power, no C or D layer, no
transport and no universal data.

**Relation to the engine's cochain.** The engine's cochain is

\[
x_B^{\rm eng}=\binom m2\,\Pih^D(\Lambda a,\Lambda a)+\Pi Q_D(\Lambda z)
=x_B^R(\smile^\Pi)+\delta\Bigl[\binom m2\Pi t_\Lambda
+\Pi\bigl(\Lambda s\smile H(\Lambda z\smile_{k-5}\Lambda z)\bigr)\Bigr],
\]

with \(a=\rho A\) and \(t_\Lambda=\rho\bigl((\Lambda A-(\rho\Lambda A)^\sim)/2\bigr)\).
Any other cup-\(i\) structure changes \(x_B^R\) by a further coboundary.

**Derivation.**

1. *The lift.* The lower presentation of an A relation holds D, C and B
   only. With target B the lift is \((A,B_A)\), solved through B, with
   \(\delta B_A=\Pi Q_D(\rho\Lambda A)\).
2. *Products through B.* The reflection gauge has A component
   \(H\Lambda(\cdot)=0\), and its B correction cancels the \(HQ_D\) terms of
   \(\Phi\). So through B the native product is
   \[
   (xy)_A=x_A+y_A,\qquad
   (xy)_B=x_B+y_B+\alpha_R(x_A,y_A),\qquad
   \alpha_R(p,q)=\Pih^D(\Lambda\rho p,\Lambda\rho q).
   \]
   \(\alpha_R\) is \(\mathbf F_2\)-bilinear, so every bracketing of the
   power has B component
   \(mB_A+\binom m2\alpha_R(A,A)\equiv\binom m2\alpha_R(A,A)\).
3. *The A step.* The A step divides by the boundary state of the gauge
   \((u,0,0,0)\). By (A13) and \(\delta Q_D=Q_D\delta\) (A3), its B
   component is \(\Pi Q_D(\Lambda z)\). The B layer is additive when the
   right factor has \(A=0\). So
   \(x_B^{\rm eng}=\binom m2\alpha_R(A,A)+\Pi Q_D(\Lambda z)\), and this is
   closed.
4. *The \(\alpha\) term is exact.* \(\Lambda a\) is closed on the bar,
   \(\Lambda a\smile_n\Lambda a=\Lambda a\), and
   \(e(\Lambda a)=\Lambda a\smile_{n-1}\Lambda a\). Hence
   \(h^D(\Lambda a,\Lambda a)=e(\Lambda a)+s\smile\Lambda a=\delta t_\Lambda\)
   (the lower relation \(dt=e+s\smile a\) of
   [secondary_operations.md](secondary_operations.md)). Under any
   cup-\(i\) structure its class is \((\operatorname{Sq}^1+s)\rho[A]=\rho\beta_s\rho[A]=0\).
5. *The class.* \(Q_D(\Lambda z)\) of the closed \(\Lambda z\) represents
   \(D[z]\), and \(\Pi^*\) is canonical.
6. *The Bockstein.* Writing \(u=\widetilde z+2w\) gives
   \(\delta_s\widetilde z/2=\tfrac m2A-\delta_sw\), so
   \(\beta_sz=\tfrac m2[A]\).

**Well-definedness.**

- *Changing \(u\) by a cocycle \(v\).* The cochain \(x_B^R\) changes by
  \(D^R(\rho v)+\deltah^D(z,\rho v)\), that is by
  \(\operatorname{Dbar}[v]\). Every Bockstein preimage of \(\tfrac m2[A]\)
  occurs. So the indeterminacy is exactly
  \(D(\ker\beta_s)=\operatorname{Dbar}H^{k-4}(X;\mathbf Z_s)\), which is the
  denominator of \(B_k\).
- *Changing the representative of \(A\).* Replacing \(A\) by
  \(A+\delta_sy\) changes \(u\) by \(my\) and leaves \(z\) unchanged.
- *The lift.* \(B_A\) enters only as \(mB_A\equiv0\). Shifting the lift by
  \(h\in H\) changes \(t\) by \(mh\).
- *The comparison.* The group bar or cell comparison, the choice of
  \(f,g,h'\) and the cup-\(i\) structure all leave the class unchanged.
- *Consistency.* \(D(D(z))=0\): for \(m\ge4\) because
  \(\operatorname{Sq}^1z=sz\), and for \(m=2\) by the Cartan formula and
  \(\operatorname{Sq}^1\operatorname{Sq}^2=\operatorname{Sq}^3\). So \(x_B\)
  lies in \(\ker D\), as a B class must.

**Conditions.**

- *Twists.* \(s\) and \(\omega\) are arbitrary.
- *Legality.* Every multiplied state is flat. The gauge \((u,0,0,0)\) lies
  off the legal locus, but its B component \(Q_D(\rho u)\) does not depend
  on the branch.
- *Degree six.* With target B or C, no \(\gamma_6\) and no universal pair
  source enters the B component. With target D, the engine refuses
  \(\gamma_6\) of two nonzero A layers unless
  `FERMIONAHSS_DEGREE_SIX_A_STACKING=1`, and records no row. The B part of
  the class is still (X1).
- *Odd \(m\).* The relation is split or three-local and has no B
  coordinate (Section 6.3).
- *The engine's B solve.* The rows are the marked cocycles and the
  coboundaries only. For the solver's particular \(u\), \(D[z]\) may
  therefore differ from \(\sum c_j[e_{B,j}]\) by a nonzero
  \(\operatorname{Dbar}\) class. The ordinary B solve then fails, and the
  complete measurement and the gauge search take over. The cell projection
  above, which the primary-operation row uses, works modulo exactly
  \(\operatorname{Dbar}H^{k-4}\) and has no such failure.

**Examples.**

- *Z4^Tf and Cs spinless.* C2, \(s=\omega=[1]\), \(k=4\). On the rank-one
  C2 resolution \(\delta_s:C^0\to C^1\) is \(\pm2\) and \(A=\beta_s(1)\), so
  \(z=1\) and \(x_B=\omega=x^2\), the B generator. This gives \(2A=B\) and
  Z/16 ([extension-paper-comparisons.md](extension-paper-comparisons.md)).
- *Pin⁻.* C2, \(s=[1]\), \(\omega=0\), \(k=6\). \(A=\beta_s(x^2)\), so
  \(z=x^2\) and \(x_B=x^4+x\operatorname{Sq}^1x^2=x^4=b\). This gives
  \(2A=B\) and Z/16 ([examples/pin_minus.g](../examples/pin_minus.g)).
  - On this resolution every mod-two coboundary vanishes. The B component
    \(\alpha_R(A,A)\) of the square is a coboundary, so it is zero there.
  - The generator \(b\) comes from the boundary of the A gauge,
    \(\Pi Q_D(\rho\Lambda u)\).

## 4. B over C

**Setting.** A B generator \(b\) has marked cocycle
\(B\in Z^{k-2}(R;\mathbf F_2)\) and order \(2\), in package degrees
\(k=2,\ldots,6\). The relation is measured through C, either as target
layer or inside a complete measurement. The recorded C coordinates are the
same in every case.

**Cohomology formula.**

\[
x_C=(\operatorname{Sq}^1+s)\,b=\rho\beta_s(b)\ \in C_k ,
\tag{X2}
\]

- At \(k=2\) the formula reduces to \(x_C=s\cdot b\).
- The map is \(\mathbf F_2\)-linear.
- \(\omega\) does not enter the value. It enters only through the domain
  (\(Db=0\) and \(\operatorname{Psi}b=0\)) and through the quotient
  (\(D\) and `Tau`).

**Cochain formulas on R.** The engine's cochain is

\[
r_C=\Pih^D(\Lambda B,\Lambda B)
=\Pi(\Lambda B\smile_{k-3}\Lambda B)+\Pi(\Lambda s\smile\Lambda B)\pmod 2 .
\]

At \(k=2\) only the second term is present. Three cocycles on R represent
its class:

- **(R1)** With the transported structure, literally:
  \(r_C=B\smile^\Pi_{k-3}B+s\smile^\Pi B\).
- **(R2)** With any cup-\(i\) structure on R:
  \(r'_C=B\smile_{k-3}B+s\smile B\). With HAP's diagonals this is
  `cupMod2(k-3,k-2,B,k-2,B)+cupMod2(0,1,s,k-2,B)`.
- **(R3)** With no cup product at all:
  \(r''_C=\rho\bigl((\widetilde BM_s)/2\bigr)\), the reduction of the
  twisted Bockstein.

On the bar the three expressions agree literally:
\(h^D(x,x)=e(x)+s\smile x=\rho(\delta_s\widetilde x/2)\) by
(C2) and (C6). On R they differ by coboundaries.

| k | \(r_C\) |
| ---: | --- |
| 2 | \(s\smile B\) |
| 3 | \(B\smile B+s\smile B\) |
| 4 | \(B\smile_1B+s\smile B\) |
| 5 | \(B\smile_2B+s\smile B\) |
| 6 | \(B\smile_3B+s\smile B\) |

**The coordinates.** Solve \(r=\sum_jc_je_{C,j}+\delta\pi\) modulo two
against the marked C cocycles. If that fails, the class meets
\(DH^{k-3}+\operatorname{Tau}E_3^{k-4,0}\). The coordinates are then those
of `layers.C.cell.project` applied to the class, which is the same
projection the engine already applies to D.

None of the following is needed: the lower flat-lift components, the
stacking model, the transport or the gauge search.

**Derivation.**

1. *The lift and the square.* The flat lift has \(A=0\) and leading
   component \(B\). Its embedding has \(\Phi_B=\Lambda B\), since
   \(Q_D(\rho0)=0\). The square has A component 0, B component
   \(2\Lambda B+\alpha(0,0)=0\), and C component \(2\Phi_C+\beta\equiv\beta\).
2. *The C correction \(\beta\) at \(A=A'=0\).* At \(k=2\),
   \(\beta_2=\delta B\smile B+s\smile(B\smile B)=s\smile B\)
   ([low_degree_stacking.md](low_degree_stacking.md)). For \(k=3,\ldots,6\)
   the legal \(\beta\) reduces to \(h^D(B,B)\), because its
   other terms vanish identically at \(A=A'=0\):
   - \(h^D(B+B',\alpha)=0\);
   - the three splitting carries vanish, since they are pointwise products
     with \(g(0)=0\);
   - \(V_0(0,0)=0\), because every monomial carries an A-derived factor;
   - \(V_1\) and \(V_2\) vanish, as relative contractions on zero fibers;
   - \(V_3\) vanishes by its explicit short circuit;
   - the degree-five term \(\rho A\smile\rho A'\) is zero.
3. *The reflection.* The native A and B components are zero, and
   \(\tau'(0;0)=0\). So the reflection adds nothing, and native
   \(C=\Pi\beta\).
4. *The reduction.* The A and B layers have zero right-hand sides. The
   solvers return the zero particular solution, so no A or B gauge acts.
   The C layer is then solved as above.
5. *The class.* \(\Pi^*\) is canonical, and \(\operatorname{Sq}^1\) and cup
   products are natural. Therefore \([r_C]=\operatorname{Sq}^1b+sb\).

**Well-definedness.**

- *Lower lift components.* \(C_b\) enters as \(2\Phi_C\equiv0\), and
  \(D_b\) is not read.
- *Changing the representative.* Replacing \(B\) by \(B+\delta y\) changes
  \(r_C\) by
  \(\delta[x\smile_r\delta y'+Q^1(y')+s\smile y']\), with \(x=\Lambda B\),
  \(y'=\Lambda y\) and \(r=k-2\).
- *Descent to \(B_k\).* For \(u\in H^{k-4}(\mathbf Z_s)\),
  \((\operatorname{Sq}^1+s)\operatorname{Dbar}(u)=\rho\operatorname{Dtilde}(\rho u)=\operatorname{Tau}(2u)\)
  modulo \(DH^{k-3}\). This is (10)–(12) of
  [dimension_indexed_differentials.md](dimension_indexed_differentials.md)
  evaluated at \(A=2U\): every term of \(F\) vanishes at \(a=0\), and
  \(G=\operatorname{Sq}^2\rho U+\omega\rho U\). Since \(2u\in\ker\operatorname{Dbar}\),
  the map (X2) descends to \(B_k\) exactly modulo \(\operatorname{im}\operatorname{Tau}\).
- *The value lies in \(\ker\operatorname{Dtilde}\).* By (S4), (S5) and (S7),
  \(\operatorname{Dtilde}((\operatorname{Sq}^1+s)b)=2\operatorname{Psi}(b)\)
  when \(Db=0\). This is zero on E6, since the image of
  \(\operatorname{Dtilde}\) is 2-torsion.
- *Gauges.* In the fallback, a B-stage gauge (\(y\) closed) shifts C by
  \(D[y]\), and an A-stage gauge \((u,y)\) shifts it by
  \(\operatorname{Tau}[u]\). The recorded coordinates are therefore those
  of \([r_C]\) in \(C_k\).

**Conditions.**

- *Legality* is automatic: \((0,B)\) is legal, and no off-legal prism
  branch is reached.
- *Twists.* \(s\) and \(\omega\) are arbitrary.
- *Literal agreement.* The literal equality (R1) uses the comparison
  actually built, with \(fg=1\). The class does not need it.
- *When \([r_C]\) meets \(\operatorname{im}D+\operatorname{im}\operatorname{Tau}\)*,
  the engine's ordinary C solve fails and
  `KOAHSS_ExtensionGaugeReduce` runs without proposed coordinates. It then
  refuses when a lower generator is free or when there are more than 32
  normal forms. The cell projection, which the primary-operation row uses,
  has no such limit.

**Examples.** Here \(H^*(C2;\mathbf F_2)=\mathbf F_2[x]\).

| Case | (X2) | Recorded |
| --- | --- | --- |
| C2, \(s=[1]\), \(\omega=0\), \(k=2\) | \(s\cdot1=x=C\) | \(2B=C+2D\) |
| C2, \(s=\omega=0\), \(k=3\) | \(\operatorname{Sq}^1x=x^2=C\) | \(2B=C+D\) |
| C2, \(s=\omega=[1]\), \(k=4\) | \((\operatorname{Sq}^1+x)x^2=x^3=C\) | \(2B=C\) |
| Pin⁻, \(k=6\) | \((\operatorname{Sq}^1+x)x^4=x^5=C\) | \(2B=C\) |
| Z2^f × Z4, \(s=\omega=0\), \(k=3\) | \(\operatorname{Sq}^1y=0\) | \(2B=6D\) |
| Z2², \(s=[1,1]\), \(k=2\) | \(s=x_1+x_2\) | C coordinates \([1,1]\) |

The sources are [extension-paper-comparisons.md](extension-paper-comparisons.md),
[examples/pin_minus.g](../examples/pin_minus.g) and
[tst/stacking_extensions.tst](../tst/stacking_extensions.tst).

For C4 with \(s=[1]\) at \(k=6\), \(b=z^2\) and
\((\operatorname{Sq}^1+y)z^2=yz^2=D(yz)\). The value lies in
\(\operatorname{im}D\), consistent with the zero C layer there
([examples/c4_signed.g](../examples/c4_signed.g)).

## 5. C over D

**Setting.** A C generator \(c\) has marked cocycle
\(C\in Z^{k-1}(R;\mathbf F_2)\) and order \(2\), in package degrees
\(k=1,\ldots,6\). Its lower group is \(D_k\) alone. The target layer is D
whenever some D generator has even order or is free; otherwise the row is
zero without a measurement.

**Cohomology formula.** The operation \(D\) of the \(d_2\) from row \(-1\)
to row \(-2\), applied to \(c\), lifts to \(\mathbf Z_s\):

\[
x_D=\rho^{-1}\bigl(\operatorname{Sq}^2c+s\operatorname{Sq}^1c+\omega c\bigr)
=W(\operatorname{Sq}^2c+\omega c)+\beta_s(\operatorname{Sq}^1c)
\ \in D_k/2D_k .
\tag{X3}
\]

- \(W\) is an integral lift of \(\operatorname{Sq}^2c+\omega c\). It exists
  because \(\beta_s(\operatorname{Sq}^2c+\omega c)=\operatorname{Dtilde}(c)=0\).
- \(\rho^{-1}\) is well defined for two reasons:
  - \(\rho:H^{k+1}(\mathbf Z_s)/2\to H^{k+1}(\mathbf F_2)\) is injective
    with image \(\ker\beta_s\);
  - \(\beta_sD=\operatorname{Dtilde}\), because
    \(\beta_s(s\operatorname{Sq}^1c)=\beta_s\rho\beta_s\operatorname{Sq}^1c=0\).
- The map \(H/(2H+I_k)\to H^{k+1}(\mathbf F_2)/\rho(I_k)\) is injective.
  So the recorded row, taken modulo two on the D generators of even or
  infinite order, is the unique \(t\) with
  \(D(c)\equiv\sum_jt_j\rho(d_j)\) modulo \(\rho(I_k)\).
- Equivalently,
  \(\beta_s(\operatorname{Sq}^1c)\equiv[\widetilde s]\smile\beta(c)\)
  modulo \(2H\), where \(\beta\) is the untwisted integral Bockstein.

| k | \(x_D\) |
| ---: | --- |
| 1 | \(\rho^{-1}(\omega c)\) |
| 2 | \(\rho^{-1}(sc^2+\omega c)\) |
| 3–6 | (X3) |

The term \(s\operatorname{Sq}^1c\) is essential when \(s\ne0\). The class
uses no `Tau`, `Psi` or `T` data, no universal source, no calibration
constant, no prism and no rational phase.

**Cochain formula on R.** Use any mod-two cup-\(i\) structure on R.

1. \(E_R=C\smile_{k-3}C+\omega\smile C\in Z^{k+1}(R;\mathbf F_2)\). This is
   \(\omega C\) for \(k\le2\).
2. \(Y_R=C\smile_{k-2}C\), or \(Y_R=\rho(\widetilde CM/2)\), a
   representative of \(\operatorname{Sq}^1c\). It is zero for \(k=1\).
3. Solve \(\delta_sD_R=-(\delta_s\widetilde E_R)/2\) over \(\mathbf Z\).
   The right side is the \(\operatorname{Dtilde}\) cocycle of the
   non-natural backend branch, and it is a coboundary because
   \(\operatorname{Dtilde}(c)=0\).
4. Set
   \(x_D^R=2D_R+\widetilde E_R+(\delta_s\widetilde Y_R)/2\in Z^{k+1}(R;\mathbf Z_s)\).
   Both halvings are exact.
5. Reduce \(x_D^R\) as the engine reduces a D residual: against the marked
   D cocycles plus \(\delta_sC^k\), or with the cell projection
   `layers.D.cell.project` for the incoming images. Keep the coordinates
   modulo two on the D generators of even or infinite order.

A change of \(D_R\) by a cocycle changes \(x_D^R\) by an element of
\(2H\).

*Mod-two route.* Since \(2D_R\equiv0\) and
\(\rho(\delta_s\widetilde Y/2)=e(Y)+s\smile Y\) up to a coboundary, where
\(e(Y)\) represents \(\operatorname{Sq}^1\operatorname{Sq}^1c=0\), one may
solve instead over \(\mathbf F_2\)

\[
C\smile_{k-3}C+s\smile(C\smile_{k-2}C)+\omega\smile C
=\sum_jt_j\rho(d_j)+\rho(i)+\delta z .
\]

Here \(i\) runs over integral representatives of \(I_k\). The left side is
the non-natural `backend.primary("D")` of
[gap/cochains.gi](../gap/cochains.gi).

**Relation to the engine.** The engine's cochain is
\(x_D^{\rm nat}=2D'_R+\Pi\gamma_k(\Lambda C,\Lambda C)-2\Pi Hg_k(0,0,\Lambda C)\)
with \(\delta_sD'_R=-\Pi g_k(0,0,\Lambda C)\). It agrees with (step 4)
evaluated with \(\smile^\Pi\) modulo \(2Z^{k+1}(R;\mathbf Z_s)+\delta_sC^k\).
It is not literally equal: \(\Pi\) of a 0/1 lift is an integral lift of the
transported cochain, not its 0/1 lift. Only the class modulo \(2D_k\) is
reproduced, not the engine's integer row.

**Derivation.**

1. *The lift.* It is \(\widetilde c=(0,0,C,D_c)\), with
   \(\delta_sD_c=-g_k(0,0,C)=-\tfrac12\delta_s\widetilde{E(C)}\). This is
   the pure-C rule (27) of
   [dimension_indexed_differentials.md](dimension_indexed_differentials.md),
   and it holds in every native branch:
   - \(k=1\): the degree-one product;
   - \(k=2\): \(g_2\);
   - \(k=3\): the phase \(\mathsf h(EC)+\Theta_1(0,0)\);
   - \(k=4\): `closed_ab_upper.J`;
   - \(k=5,6\): `pure_c_g`.
2. *The square.* It has A, B and C components zero, since \(\alpha(0,0)=0\)
   and \(\beta(0,0;0,0)=0\). Its D component is
   \(2D_c+\gamma_k(\widetilde c,\widetilde c)\), with
   \[
   \gamma_k(\widetilde c,\widetilde c)=\widetilde{E(C)}
   +\tfrac12\delta_s\bigl(C\smile_{k-2}C\bigr)^\sim+\delta_s\varepsilon_k ,
   \]
   where \(\varepsilon_k\) is an explicit integral cochain in each degree
   (the pure-C phase of the code, and for \(k\ge4\) the term
   \(-(s\smile C)^\sim\) of the K change of a pure-closed state). Its
   coboundary \(\delta_s\varepsilon_k\) is an integral coboundary and does
   not change the class.
3. *The class.* \(x_D=W+\tfrac12\delta_s(\operatorname{Sq}^1C)^\sim+\delta_s\varepsilon_k\),
   with \(W=2D_c+\widetilde{E(C)}\) and \(\rho W=\operatorname{Sq}^2c+\omega c\).
   Since \(\rho(\delta_s\widetilde y/2)=e(y)+s\smile y\),
   \[
   \rho[x_D]=\operatorname{Sq}^2c+\omega c+\operatorname{Sq}^1\operatorname{Sq}^1c+s\operatorname{Sq}^1c=D(c).
   \]
4. *The reduction.* In layers A, B and C the right-hand sides are zero. In
   D, the boundary states that keep a state pure-D add only four kinds of
   terms: \(\delta_sw\); \(\operatorname{Dtilde}\) cocycles; `Psi` images
   \(J_{k-1}(0,y,z)\); and `T` images \(J_{k-1}(u,y,z)\). So the recorded
   coordinates are those of \([x_D]\) in \(D_k\).

**Well-definedness.** The class lives in
\(D_k/2D_k=H^{k+1}(\mathbf Z_s)/(2H^{k+1}+I_k)\).

- *The choice of \(D_c\).* Changing \(D_c\) by a cocycle changes \(x_D\) by
  an element of \(2H\).
- *The cocycle representative and the cup-\(i\) structure.*
  \(\rho[x_D]=D(c)\) depends only on \([C]\).
- *Changing \(c\) by \(\operatorname{im}D\).* By the Cartan formula and
  \(\operatorname{Sq}^2\operatorname{Sq}^2=\operatorname{Sq}^3\operatorname{Sq}^1\),
  \[
  D\circ D=\rho\operatorname{Dtilde}\circ(\operatorname{Sq}^1+s).
  \tag{X4}
  \]
  So \(\rho^{-1}D(Dy)\equiv\operatorname{Dtilde}(\operatorname{Sq}^1y+sy)\in I_k\).
- *Changing \(c\) by \(\operatorname{im}\operatorname{Tau}\).* The relation
  class is a group invariant under `koFull`'s assumptions of gauge
  completeness and abelian gauge classes. That is the argument used here;
  no cochain identity for this case is claimed.

**Conditions.**

- *Legality* is automatic, since \(A=B=0\). The D equation of the lift is
  solvable exactly when \(\operatorname{Dtilde}(c)=0\), which holds for
  every E4 survivor.
- *Twists and comparisons.* \(s\) and \(\omega\) are arbitrary. Both the
  group-bar and the cell comparison give the class.
- *Degree six.* The formula holds at \(k=6\), with no pair-source refusal.
- *Consistency with later relations.* A row computed from (X3) is the
  relation of a lift shifted by an element of \(D_k\). Suppose a later
  relation is measured through D, has a nonzero coefficient on \(c\), and
  multiplies the engine's lift of \(c\). Then the two computations must use
  the same lift of \(c\) ([extensions.md](extensions.md), "Layers and
  measured relations").

**Examples.**

- *Table VII of Wang and Gu*, C2
  ([extension-paper-comparisons.md](extension-paper-comparisons.md)). Each
  of the following gives \(2C=D\):
  - \((0,[1])\), \(k=1\): \(D(1)=x^2\);
  - \(([1],0)\), \(k=2\): \(D(x)=x^3\);
  - \((0,0)\), \(k=3\): \(D(x^2)=x^4\);
  - \(([1],[1])\), \(k=4\): \(D(x^3)=x^5\). This is also Cs spinless.
- *Pin⁻*, \(k=6\). \(\operatorname{Sq}^2x^5=0\) but
  \(s\operatorname{Sq}^1x^5=x^7\), so \(2C=D\). Without the
  \(\beta_s(\operatorname{Sq}^1c)\) term the relations would present
  Z/8 ⊕ Z/2 instead of Z/16.
- *Z2^f × Z4*, \(k=3\). \(D(y)=y^2=\rho(u^2)\), where \(u^2\) generates
  \(H^4(\mathbf Z/4;\mathbf Z)\). So \(2C=D\).
- *Z4^f × Z2*, \(k=3\), \(\omega=a^2\). \(C_3=\langle b^2\rangle\),
  \(D_3\) is generated by the lift of \(b^4\), and
  \(D(b^2)=b^4+a^2b^2\equiv b^4\). So \(2C=D\).
- *Z2², \(s=a+b\), \(k=2\).* \(D(a)=a^3+a^2b\) and \(D(b)=ab^2+b^3\) are
  independent, and \(I_2=0\). This gives the rows \([1,0]\) and \([0,1]\)
  of [tst/stacking_extensions.tst](../tst/stacking_extensions.tst).
- *Z2^f × Z2 × Z4*, \(k=4\). The recorded row \(2C=3D_1+4D_2\) has parity
  \((1,0)\), consistent with a class defined modulo two.

## 6. The non-adjacent pairs

The binary non-adjacent components come from the relation (X4) among
primary operations. Both sides of (X4) equal

\[
\operatorname{Sq}^3\operatorname{Sq}^1+s\operatorname{Sq}^3+s\operatorname{Sq}^2\operatorname{Sq}^1
+s^3\operatorname{Sq}^1+(\operatorname{Sq}^1\omega)\operatorname{Sq}^1+s(\operatorname{Sq}^1\omega),
\]

and the middle term \(s\operatorname{Sq}^1\) of \(D\) contributes nothing to
\((\operatorname{Sq}^1+s)D(\operatorname{Sq}^1+s)\): by the Cartan formula
\((\operatorname{Sq}^1+s)(s\operatorname{Sq}^1)(\operatorname{Sq}^1+s)=0\).
Untwisted, (X4) is \(\operatorname{Sq}^2\operatorname{Sq}^2=\operatorname{Sq}^3\operatorname{Sq}^1\).

A component measured two layers below its generator is the value of the
secondary operation attached to (X4); a component three layers below is
tertiary. Their cochains need a *natural* primitive of the relation. In the
stacking model that primitive is built from the universal data: the pair
primitives \(\mathcal V_n\) of [ALL_COCHAIN_STACKING.md](ALL_COCHAIN_STACKING.md),
the words \(\chi_n\) and \(\zeta_{i,n}\) and the phases \(\Theta\) and
\(V_1,V_2,V_3\) of [universal_helpers.md](universal_helpers.md), and prisms
on \(X\times I\). A primitive solved on R by linear algebra is ambiguous by
every cocycle of its degree, so these components have no closed cup-\(i\)
formula on R in general. The parts that are primary are listed below.

### 6.1 A over C

An A generator of order \(m=2^e\) has target layer C, in degrees
\(k=4,5,6\): every D generator lies in \(mH\). With \(\delta_su=mA\) and
\(\bar u=\rho u\) as in Section 3:

- **The B coordinate** is \(c_B=D(\bar u)\) modulo
  \(\operatorname{Dbar}H^{k-4}\), exactly as in Section 3.
- **The C coordinate** is not primary. It is the value of the secondary
  operation of (X4) on the defining system formed by the B gauge
  \(y'\), with \(\delta y'=D(\bar u)+b\), where \(b\) is the combination of
  the marked B cocycles with coefficients \(c_B\), together with the B
  component of the lift for \(m=2\) or the second null-homotopy of
  \(\rho u\) for \(m\ge4\). To this value the model adds a natural primary
  class \(\theta_m\) of \([u\bmod m]\), \(s\) and \(\omega\), which its
  calibrated universal data fix and these notes do not determine. The
  value lies in
  \(H^{k-1}(\mathbf F_2)/\bigl(DH^{k-3}+\operatorname{Tau}E_3^{k-4,0}
  +[m=2]\,(\operatorname{Sq}^1+s)B_k\bigr)\).

Given the defining cochains, the terms \(D(\bar u)\), \(Q_D(y')\),
\(h^D(B_A,B_A)\), \(h^D(Q_D\bar u,b)\) and the C parts of the lower lifts are
cup-\(i\) words on R. The witness of the relation (the pair primitives on
the diagonal pairs \((jA,A)\), the prism \(H_{k-1}(u,0)\) of (A9) through the
calibrated \(\tau'_n\), and the pointwise digits and carries), the class
\(\theta_m\), and complete lifts that extend through D (which need `Tau`,
`Psi` and `T`) are not. The recorded relation \(2A=C\) of S4 spin-1/2
at \(k=4\) ([extension-paper-comparisons.md](extension-paper-comparisons.md))
is an instance, with \(m=2\), \(\bar u=1\) and \(c_B=0\), which this note
does not evaluate.

### 6.2 B over D

A B generator \(b\) has target layer D, in degrees \(k=2,\ldots,6\); the
measurement is complete, with the lift \((0,b,c,D_b)\), \(\delta c=Q_D(b)\).
Its C component \(x=h^D(b,b)\) represents \((\operatorname{Sq}^1+s)b\)
(\(x=s\cdot b\) at \(k=2\)). Since the target layer depends on the lower
presentation alone, the engine first writes

\[
x=\sum_jc_je_{C,j}+f_{k-1}(u,y)+\delta\pi
\]

with a legal gauge pair \((u,y)\) of degree \(k-1\): \(u=0\) with \(y\)
closed contributes \(D[y]\), and \(u\ne0\) contributes
\(\operatorname{Tau}[u]\).

- When \((u,y)=0\), the D coordinate, read in \(H/2H\) as a class in
  \(H^{k+1}(\mathbf F_2)/\rho(I_k)\), is the model-normalized value on
  \(b\) of the twisted secondary operation of (X4), with indeterminacy
  \(D(H^{k-1})+\rho\operatorname{Dtilde}(H^{k-2})\). Untwisted it is exactly
  that operation. With twists it is fixed only up to primary terms that the
  universal data fix and these notes do not determine.
- When \((u,y)\ne0\), a `Psi`-type (\(u=0\)) or `T`-type (\(u\ne0\)) phase of
  the gauge pair is added.

\(Q_D(c)\), the pure-C gauge boundary
\(\delta_s\mathsf h(E\pi)-\mathsf h(E\delta\pi)\), the pure-C products with
their K changes, and the C lifts of Section 5 are cup-\(i\) words and exact
halvings on R. The model's natural primitive of (X4) on \(b\), the
difference \(\rho\gamma_k(\widetilde b,\widetilde b)-Q_D(c)\), is not: it
contains the \(\zeta_{i,n}\) and \(\chi_n\) words inside \(\Theta\) and the
prisms, quarter lifts, prisms on \(X\times I\), and in degrees four to six
the Eilenberg–Zilber primitive of the \(\Theta\) pair. A primitive solved on
R instead shifts the value by an arbitrary cocycle. The relation \(2B=D\) of
C4 with \(s=[1]\) at \(k=6\) ([examples/c4_signed.g](../examples/c4_signed.g))
is an instance with \([x]=yz^2=D(yz)\), that is, with a nonzero B gauge,
which this note does not evaluate.

### 6.3 A over D

At the odd primes the answer is primary.

- \(p\ge5\), every \(k\), and \(p=3\) with \(k\le4\): the relation splits,
  \(x_D=0\) ([extensions.md](extensions.md), "Localization at the
  primes").
- \(p=3\), \(k=5,6\), \(m=3^e\): with \(Y\in H^{k+1}(X;\mathbf Z_s)\) any
  integral class with \(\rho_3Y=P^1_s\rho_3a\),
  \[
  x_D\equiv2\cdot3^{e-1}\,Y\pmod{3^eH_D}.
  \]
  At \(k=5\), \(Y=[A\smile_sA\smile_sA]\), one cup product on R. At
  \(k=6\), \(Y\) exists because the flat lift exists; evaluating
  \(P^1_s\) on R needs a cyclic diagonal, since \(P^1\) is not a
  cup-\(i\) word, and [gap/native_coherence.gi](../gap/native_coherence.gi)
  supplies the cyclic diagonals of degrees zero to two. These values
  reproduce the documented relations: \(3A=-16D\) for Z/3,
  \(9a=6d\) modulo nine for Z/9 and the group Z/9 ⊕ Z/81 for Z/27 in
  degree five, and \(3A=14D\) for Z/3 × Z and \(9A=42D\) for Z/9 × Z in
  degree six.

At the prime two (\(k=4,5,6\)) the component is tertiary: its cochain in the
model contains the sources \(V_1,V_2,V_3\) inside \(\Omega_k\), the pair
primitives on the diagonal pairs \((jA,A)\), the Eilenberg–Zilber primitive
of the \(\Theta\) pair, the pointwise digits and carries, the prisms of the
gauges, and the D component of the flat lift, whose equation is the
transfer of \(\delta_s\mathcal O_{k-3}\) of
[tertiary_operations.md](tertiary_operations.md). Multiplying by \(m\) does
not clear the dyadic phases. No cup-\(i\) formula on R is established, and
none is proved impossible. In degree six the engine refuses the product of
two states with nonzero A layers unless
`FERMIONAHSS_DEGREE_SIX_A_STACKING=1`.

## 7. Summary

| Pair | Degrees k | Class of the component | Type | Cup-\(i\) on R |
| --- | --- | --- | --- | --- |
| A–B | 4, 5, 6 | \(D(z)\), \(\beta_sz=\tfrac m2[A]\), \(z=\rho u\), \(\delta_su=mA\) | primary | yes: \(z\smile_{k-6}z+s\smile(z\smile_{k-5}z)+\omega\smile z\) |
| B–C | 2–6 | \((\operatorname{Sq}^1+s)b=\rho\beta_sb\) | primary | yes: \(B\smile_{k-3}B+s\smile B\), or \(\rho(\widetilde BM_s/2)\) |
| C–D | 1–6 | \(\rho^{-1}(\operatorname{Sq}^2c+s\operatorname{Sq}^1c+\omega c)\) modulo \(2D_k\) | primary | yes: one integral solve, or the mod-two route |
| A–C | 4, 5, 6 | B part \(D(\bar u)\); C part the secondary operation of (X4) plus \(\theta_m\) | secondary | partial: the defining cochains only |
| B–D | 2–6 | the secondary operation of (X4) on \(b\), plus a gauge phase when \([x]\) meets \(\operatorname{im}D+\operatorname{im}\operatorname{Tau}\) | secondary | partial: the pure-C and gauge terms only |
| A–D, \(p\ge5\) | all | \(0\) | split | trivially |
| A–D, \(p=3\) | 5, 6 | \(2\cdot3^{e-1}Y\), \(\rho_3Y=P^1_s\rho_3a\) | primary | \(k=5\): one cup product; \(k=6\): a cyclic diagonal for \(P^1\) |
| A–D, \(p=2\) | 4, 5, 6 | tertiary | tertiary | no formula established |

For the adjacent pairs, reading a relation on R needs only R's coboundary
matrices, a mod-two cup-\(i\) structure on R, 0/1 lifts with exact halving,
linear solves over \(\mathbf Z\) and \(\mathbf F_2\), and the E6 cell
projection of the target layer. No transported model, no flat lift beyond
the leading cocycle and no universal datum enters these classes. The
primary-operation rows of `koFull` evaluate them this way: \(D\) with HAP's
diagonals (`backend.nativePrimary`), the integral solve
\(\delta_su=mA\) for A–B, the Bockstein route (R3) for B–C, the integral
lift \(x-2w\) with \(\delta_sw=\delta_sx/2\) of the 0/1 lift \(x\) of
\(D(c)\) for C–D, and `layers.<name>.cell.project` for the coordinates.
