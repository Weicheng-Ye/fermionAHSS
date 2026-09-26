# Resolution transfer: review and implemented prerequisites

Review date: 2026-09-26. The input was version 3 of the proposed
four-cochain transfer plan, dated 2026-09-26 and based on commit `760baf1`.
The plan was supplied as `~/Downloads/plan.md`.

**Status at the initial review:** the production transfer was not implemented.
The subsequent [native-resolution implementation](resolution-extensions.md)
adds an opt-in transfer with independent reference certification and fallback;
it does not assume the completeness theorem discussed below. The review
and prerequisite measurements in this document remain the initial findings.

At that review, `koFull` still used its complete-bar higher extension model.
The initial revision implemented exact Smith-preparation reuse, a sparse
comparison preflight, and bounded tests of several required model identities. These were
prerequisites, not an implementation of `xtimes_R` or a claim of identical
classifications on arbitrary resolutions.

The plan correctly identifies the dense bar equations and evaluations as
costs worth removing. Its proposed production certificate needs additional
mathematics. Rank-one degree zero does not imply a strict retraction, and
the stated restricted product identity does not imply completeness of
transferred gauges. Both issues have concrete counterexamples below.

## The comparison is not an arbitrary-resolution retraction

Write \(B_*\) for the normalized group-bar resolution and \(R_*\) for
the supplied integral HAP resolution. The existing comparison has chain
maps \(f:B_*\to R_*\), \(g:R_*\to B_*\), and a homotopy \(h\) with

\[
1-gf=\partial h+h\partial.
\]

On cochains put \(\Lambda=f^*\), \(\Pi=g^*\), and \(H=h^*\).
Then \(1-\Lambda\Pi=\delta H+H\delta\). Integral A and D use the
sign local system; B and C reduce the integral maps modulo two. Binary
results are reduced after each operation before they enter subsequent
integral lifts. This homotopy identity does not assert \(fg=1\).

For a concrete counterexample, add a contractible free \(\mathbf ZG\)
summand to a valid C2 resolution:

\[
R'_1=R_1\oplus\mathbf ZG u_1,\qquad
R'_2=R_2\oplus\mathbf ZG u_2,\qquad
\partial u_2=u_1,\quad\partial u_1=0.
\]

Extend the contraction by \(h_R(u_1)=u_2\), \(h_R(u_2)=0\).
This is still a resolution and still has one degree-zero generator.
The code constructs \(g(u_1)\) by coning \(g(\partial u_1)=0\), so
\(g(u_1)=0\). Consequently \(fg(u_1)=0\ne u_1\). The regression in
[extension_transfer.tst](../tst/extension_transfer.tst) constructs this
resolution explicitly and verifies \(\Pi\Lambda(0,1)=(0,0)\).

A check of \(\Pi\Lambda=1\) with just trivial and sign coefficients
is weaker than a check over \(\mathbf ZG\). The preflight checks
\(fg=1\) on every group-ring basis generator in the requested degrees.
A failed check requires a different comparison or a transfer for general
homotopy equivalences; normalizing H alone cannot repair it.

If \(fg=1\) does hold, the usual side-condition construction is available
algebraically. Set \(q=1-gf\), \(u=qhq\), and \(h'=u\partial u\).
Then \(\partial u+u\partial=q\), \(qu=uq=u\), \(fu=ug=0\), and

\[
\partial h'+h'\partial=q,\qquad fh'=h'g=0,\qquad (h')^2=0.
\]

For example, \(\partial u\partial=\partial q\) proves the first
identity, while \(\partial u^2=u^2\partial\) gives
\(u\partial u^2\partial u=0\). This is a construction to implement
and budget, not a claim that the existing h already satisfies the side
conditions. The preflight reports `sideConditionsVerified=false`.

## Flatness and reflection are conditional, not completeness

For physical degree k, write \(w=(A,B,C,D)\) with cochain degrees
\((k-3,k-2,k-1,k+1)\). Write \(N_\ell\) for the nonlinear part of
the differential in layer \(\ell\); it depends on earlier layers.
Under the retraction and normalized homotopy identities, define

\[
\Phi(w)_A=\Lambda w_A,\qquad
\Phi(w)_\ell=\Lambda w_\ell-HN_\ell(\Phi(w)_{<\ell}),
\qquad
\kappa_R(w)_\ell=\delta w_\ell+\Pi N_\ell(\Phi(w)_{<\ell}).
\]

If earlier layers are flat and \(\delta N_\ell=0\) there, then

\[
\delta\Phi(w)_\ell+N_\ell
=\Lambda\kappa_R(w)_\ell+H\delta N_\ell
=\Lambda\kappa_R(w)_\ell.
\]

This proves the proposed layerwise flatness identity (T1) under its
hypotheses. In binary layers this equation is in \(\mathbf F_2\).
In the D layer the full integral nonlinear term, including all carries,
must be used.

The plan's reflection (T2) similarly works conditionally: if a flat bar
state z has already been matched in earlier layers, subtract their known
contributions to form \(r_\ell\), and set
\(w_\ell=\Pi r_\ell\), \(g_\ell=Hr_\ell\). Square-zero, zero
normalization, and the required product identity give
\(\delta r_\ell=-\Lambda\Pi N_\ell\). Retraction gives
\(\kappa_R(w)_\ell=0\), and \(H\Lambda=0\) gives
\(\delta g_\ell+\Lambda w_\ell=r_\ell\). This establishes the
literal reflection equation \(z=d(g)\times\Phi(w)\) under those
hypotheses. It does not transfer arbitrary nonflat gauges.

In particular, `kappa_R(e)` on a gauge is not a replacement for
`Psi(d(Phi(e)))`. The existing consumers would need a separate gauge
action in all zero-unit, synthetic-equation, verification, and ordered
reduction paths. During reduction it must be `act(e,chosen)` directly;
`act(e,zero) xtimes chosen` is a different expression without a justified
coherence identity.

## A counterexample to the stated completeness assumptions

This algebraic complex tests the logical implication in the plan; it is
**not a finite-group bar complex or the calibrated package formulas**.
Use trivial sign and an integral cochain
complex R with \(R^0=\mathbf Z\) and generators e, a, c in degrees
2, 3, 4, respectively. Its only nonzero differential is \(\delta e=a\).
Add a contractible complement K with h in degree 1, b in degree 2, and
\(\delta h=b\). Let \(\mathcal C=R\oplus K\), with the split
inclusion, projection and homotopy \(H(b)=h\). They satisfy all the
strong deformation retraction identities. Use binary B/C coefficients.

Take every product to be componentwise addition. Take every nonlinear
differential term to be zero except the C output in physical degree 4:

\[
N_C(B)=ut\,c\pmod2\quad\text{for }B=ue+tb.
\]

These operations satisfy the plan's square-zero, normalization, Bianchi,
and **layer-restricted** multiplicativity assumptions:

- The image of an ordinary differential has no e coefficient. Thus the
  nonlinear term vanishes on differential images, and \(d^2=0\).
- Flatness below C forces \(\delta B=ua=0\), hence \(u=0\). The
  nonlinear term vanishes on the inputs where the plan demands
  the C-layer product identity. The other layers are linear.
- \(\delta c=0\) gives the required Bianchi identity; zero is a strict
  unit and every nonlinear term vanishes at zero.

Here \(\Phi\) is just inclusion. The physical-degree-five states
\(x=0\) and \(x'=(0,a,c,0)\) are flat. The physical-degree-four
\(\mathcal C\)-gauge with \(B=e+b\) satisfies \(d(g)=x'\).
Every R-gauge has \(B=ue\), however, and its transferred action on zero
has the form \((0,ua,0,0)\). Reflection is projection on these images.
No sequence of those R-actions can produce the c component of x'.
Thus T4 does not follow from the stated assumptions.

Full multiplicativity on nonflat gauges would exclude this example:
\(d(e+b)\ne d(e)+d(b)\) in its C layer. A proof using that stronger
identity and appropriate gauge coherence might establish completeness
for the actual formulas. That proof is still required. Alternatively,
a completed presentation could be justified by a filtered isomorphism
argument, but this needs an independently established bar quotient with
the asserted E6 filtration and an isomorphism on every marked graded
piece. A finite product audit does not establish those hypotheses.

Generic homotopy transfer results cannot be substituted without checking
their structures and coefficient assumptions. For context, transfer for
general homotopy equivalences requires additional homotopy data in
[Hogancamp's perturbation construction](https://arxiv.org/abs/1912.03843).
No such general result is asserted here for the mixed integral/binary,
globally branched stacking model.

## Implemented preflight

The internal diagnostic can be run without any nonlinear stacking call:

```gap
R := ResolutionFiniteGroup(CyclicGroup(4),8);;
backend := koAHSSHAPSpace(R,koAHSSNaturalOperations()).koAHSS([1],0,8);;
audit := KOAHSS_ExtensionTransferPreflight(backend,5);;
audit.status;                    # "checked" for this example
audit.chainRetractionVerified;   # true
audit.transferReady;             # false
audit.degrees;                   # ranks, g terms, distinct g-support counts
```

The audit accepts package degrees 3–5 and checks through k+2. It requires
resolution length at least k+3 because HAP can leave the final
contraction degree unavailable. Its optional limits record accepts
`maxSupport` (8192 distinct normalized g simplices per degree by default)
and `maxTerms` (2,000,000 processed g/f expansion terms by default).
Refusals return `status="unresolved"`, a reason, and available failure
coordinates. Support counts cover g only: they do not include the future
face/H closure or global flag tests. The limits are accounting bounds;
the underlying transport builds a chain before its size can be checked.

The initial prerequisite revision installed no new `koFull` option.
The later opt-in API is documented in [resolution extensions](resolution-extensions.md).
Degree six, low-degree adapters, calibrated formulas and numerical payloads
retain their original behavior.

## Exact Smith-preparation reuse and verification

`koAHSSSolveIntegerSystem` now reuses Smith transformations for repeated
matrices and different right-hand sides. Immutable value snapshots prevent
caller mutation from changing a cache key or a transformation. It retains
at most eight preparations and two million matrix cells, counting the
input, normal form, and square row/column transforms. This cell bound is
not a byte limit for arbitrarily large integers. Oversized preparations
are used for the current solve without being retained.

The solution formula and basis convention are unchanged. Every successful
nonempty solve still checks \(xM=b\) and every homogeneous generator
against M exactly. This optimization therefore preserves the bar
extension witnesses as well as their abstract presentations.

For the signed C4 bar coboundary with shape 243 by 729, one run measured
3558 ms of GAP CPU time for the first solve and 4 ms for a different
right-hand side using the cached preparation. Matrix construction took
2375 ms. The rank was 182 and the homogeneous rank 61. This measures
integer solving only; it is not a full `koFull` speedup measurement.

The dated [verification record](verification/extension-transfer-20260926.md)
lists the executed checks and the remaining acceptance work. In particular,
no transferred Z4T classification or arbitrary-resolution equality is
claimed. All package outputs retain `certified_ko=false`.
