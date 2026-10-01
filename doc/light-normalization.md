# Normalization of light transport cochains

The internal sources and rational potentials of
[`LightEvaluator`](../python/extension_light.py) vanish exactly on
degenerate simplices. This lets the production tasks use the private
`_prim_normalized` and `_rational_ds_normalized` paths without repeating
normalization evaluations. The public `prim` and `rational_ds` helpers
retain their on-demand checks for arbitrary cochain callables.
`FERMIONAHSS_LIGHT_NORMALIZATION_CHECKS=1` also enables those checks in
the production paths; the worker reads this option when it creates its
light evaluator.

## Required property

A positive-degree cochain is normalized when

\[
c(v_0,\ldots,v_i,v_i,\ldots,v_n)=0.
\]

For rational potentials this equality is in the rationals, not merely
modulo integers. The normalized comparison identities are

\[
1-\Lambda\Pi=\delta H+H\delta,\qquad
\Pi\delta_s=\delta_s^R\Pi.
\]

Thus a normalized closed source \(z\), with
\(\delta_Rr=\Pi z\), gives
\(\delta(\Lambda r+Hz)=z\). Normalization does not replace closedness,
the defining equation on R, or the checked comparison identities.

## Closure under the formula constructors

The resolution lifts and homotopy outputs factor through normalized
chains: an adjacent repeated vertex gives an empty chain and hence zero.
This applies to the twists as well as all transported defining data.
Sums, differences, scalar multiples, pointwise products, reduction modulo
two or three, exact divisions, and integral carries preserve zero.

Coboundaries preserve normalization. On a degenerate simplex the two
faces obtained by removing one of the repeated vertices agree and cancel;
all other faces remain degenerate. The signed first-face term has the
same cancellation because the sign cochain is zero on an identity edge.

For an interval-cut operation, every consecutive edge of the output
simplex is contained in an input interval. A repeated output edge
therefore makes at least one factor in each product zero. This covers
the cup-i operations and their sums in Q, E, Q_D, hD and polarization.
Sign transport changes only the multiplier of such a product.

Interval pullbacks preserve degeneracies. Every right-prism simplex over
a degenerate base simplex has a repeated adjacent product vertex, so the
prism of a normalized cochain is normalized term by term. These statements
also cover the interval and prism operations nested in the secondary and
tertiary phases.

The packaged chi ANF and reduced-power terms admit finite symbolic
certificates, implemented in
[`test_extension_light_normalization.py`](../python/test_extension_light_normalization.py).
For chi in degrees 0–7, every elementary degeneracy substitutes zero for
degenerate face variables and identifies the other equal face variables;
the resulting Boolean polynomial is identically zero. For the
degree-two cube and nineteen-term degree-three reduced power, every
product term contains each output edge in at least one factor. Their
standard 0,1,2 lift and rational multiples therefore vanish literally.

## Universal primitives

The remaining primitives have the form

\[
V(c)(\sigma)=\operatorname{source}(H\,s(c,\sigma))
             +\operatorname{coefficients}(F\,s(c,\sigma)),
\]

where \(s\) is the universal section and \(F,H\) are its normalized
chain operators. Each section preserves an adjacent degeneracy:

- In degree one, the repeated edge has sign and A value zero, giving
  equal vertices and a zero row and column of the omega matrix.
- In degree two, the finite-difference matrices have a zero row and
  column at the repeated edge: either a face is degenerate or two
  identical face values are subtracted.
- In degree three, the alternating face extension respects the monotone
  identification of the repeated positions. Successive differences give
  a zero outer row and an inner zero row and column; the sign transport
  across the repeated edge is one.
- For binary cochain tables, faces containing both copies are zero and
  the other entries agree after one copy is removed. The table is the
  corresponding simplicial degeneracy.

These statements hold simultaneously in every fiber and background of a
pair section. The normalized homotopy kills such a universal simplex.
Every Alexander–Whitney splitting also has a degenerate factor, so its
forward chain is empty. Consequently the primitive is exactly zero
without evaluating its source or calibration coefficients.

This covers V1, V2, V3, the binary primitives in legal beta, the rational
pair phases, the Theta pair primitive, and the degree-six unary primitive.
The section tests exercise the corresponding chain operators with signed
backgrounds and normalized inputs that need not be closed. The argument
above supplies the general normalization property; finite samples alone
would not do so.

## Production scope

Every primitive source in the fixed light tasks is a sum or composition
of Q_D, fsharp, hD, legal beta, raw_H, normalized reference data, or a
signed integral multiple of A. The rational potentials use the same
constructors together with the universal phases, carries and prisms just
described. This covers the supported package degrees and both sign twists.

An arbitrary `Cochain(degree, callable)` has no such constructor contract.
For example, a constant rational degree-one cochain takes a nonzero value
on the omitted face `(v,v)`, invalidating the normalized chain-map
calculation. Such objects must use the checked public helpers. A new
formula constructor may use the private paths only after its exact
normalization has been established.

This argument establishes normalization only. The legal-input conditions,
source and primitive identities, exact carries, calibration provenance,
residue integrality, comparison-chain checks, resource limits, and final
frame audit retain their separate roles.
