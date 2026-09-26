# Backend and operation interfaces

## HAP resolutions

The high-level `koAHSS(group,s,omega,k[,n])` wrapper constructs resolutions
for finite groups. Its first argument can also be an explicit integral
HAP resolution: `koAHSS(R,s,omega,k[,n][,options])` preserves R and its
twist basis. `koFull` accepts the same explicit-resolution first argument;
see [native extension searches](resolution-extensions.md).
For an explicitly constructed integral HAP resolution, use:

```gap
operations := koAHSSNaturalOperations();;
space := koAHSSHAPSpace(R, operations);;
backend := space.koAHSS(s, omega, maxDegree);;
```

`R` must have group action, boundary, and contracting homotopy data through
the requested degrees. Both integral coefficient rows use Z_s. Multiple
degree-zero generators are supported when each has augmentation one and
the contraction is anchored at the first generator at the identity.
An explicit resolution may model an infinite group if HAP supplies the
required integral data; the finite-group convenience wrapper does not
construct such a resolution.

The factory installs matching primary, secondary, and final tertiary
conventions. It always compares with the normalized homogeneous group-bar
resolution. `backend.naturalTransport()` and `backend.naturalBar()` access
that lazy comparison. Whole formulas are lifted, evaluated, and projected,
with homotopy corrections for both defining cochains; see
[conventions](conventions.md). A local primitive is not a replacement for
this fixed operation. The separate `backend.nativeCoherence()` tensor
interface is a research utility and does not select the higher callbacks.

## Cochains and exact arithmetic

`koAHSSCochainSpace(data)` accepts a record with:

```gap
data := rec(
    dimension := function(degree) ... end,
    differential := function(degree,s) ... end,
    cupMod2 := function(i,degreeA,a,degreeB,b) ... end,
    cupIntegral := function(i,degreeA,a,degreeB,b,sFirst,sSecond) ... end,
    operations := rec(Tau := ..., Psi := ..., T := ...)
);
```

The functions above are an interface schema, not executable example code.
`dimension(n)` is the integral cochain rank. Matrices act on row vectors:
`a * differential(n,s)` is `d_s a`. Their shape is `dimension(n)` by
`dimension(n+1)`, including zero-dimensional modules. A known zero module
must not be confused with unavailable resolution data.

`cupMod2` returns `a cup_i b`; negative i gives zero. `cupIntegral` uses the
specified signed higher-cup convention, with sign cocycles `sFirst` and
`sSecond` for the two factors. The output local system is their tensor
product. For backend calls, `true` and `false` abbreviate the configured
`s` and the zero sign twist. A bare matrix complex does not determine
coherent higher operations; missing capabilities remain unavailable.

Common backend methods are:

```gap
backend.cohomology(n,q)
backend.cohomologyData(n,q)     # group, represent(class), class(cocycle)
backend.data(n,q)              # alias of cohomologyData
backend.coboundary(n,a,sign)   # true: Z_s; false: ordinary integers
backend.cupMod2(i,p,a,q,b)
backend.cupIntegral(i,p,a,q,b,sFirst,sSecond)
backend.primary(name,n,a)      # "Dbar", "D", or "Dtilde"
backend.twists                # s and omega
```

Actual integral representatives survive page transitions. If d2 removes
odd multiples of a Z generator, the next operation is evaluated on an even
representative; abstract ranks do not supply that information.

## Higher callbacks and direct audits

`operations.Tau(ctx)`, `operations.Psi(ctx)`, and `operations.T(ctx)` return
target cocycle vectors. `ctx` contains `degree`, `cochain`, `s`, `omega`,
`backend`, and current-page `source` and `target` cells. It also records the
tertiary reference and correction coefficient. Final `T` has correction
coefficient zero. The separate legacy `TReference` callback receives its
coefficient-one correction exactly once; an explicit final `T` takes precedence.

A callback may return `fail` or an unresolved record when unavailable.
Bad callback types, invalid cocycles, failed identities, and inadequate
resolution lengths are errors. A callback is responsible for its defining
cochains, indeterminacy, and normalization; the page engine does not prove
that an arbitrary user callback is a natural cohomology operation.

`koAHSSNaturalTertiary(backend,n,A[,options])` accepts optional checked
`b,c` cochains. A computed record includes its integral `cochain`, exact
rational `phase`, `phaseNumerator`, `modulus`, `definingSystem`,
`referenceConvention`, and `kernelAudit`. `carryAudit` verifies that
`remainder + projectionCarry + sourceLiftCarry = phase`, with integral
carries. Bar calls report `transportModel="bar"`,
`inputGroupBarUsed=true`, and the actual `barSamples`/`modelSamples` count.
They also report `isLocalChoice=false`, `oldCorrectionApplied=false`,
`muR=0`, and `oddPrimaryCoefficient=2`.

To keep a callback audit:

```gap
evaluations := [];;
operations := koAHSSNaturalOperations();;
operations.T := koAHSSNaturalTCallback(rec(evaluations := evaluations));;
space := koAHSSHAPSpace(R, operations);;
```

This preserves the factory's reference metadata. The backend also retains
`tertiaryEvaluations`. Completed records retain defining vectors and equation
audits, while releasing the internal secondary evaluation graph. A supplied
valid `c` can skip a redundant full solution-family calculation; in that
case `solutionFamilyComputed=false` does not mean zero ambiguity.

## Current-page backend

A custom space may instead supply
`space.koAHSS(s,omega,maxDegree)` returning a record with
`cohomology(p,q)` and `differential(r,p,q,source,target)` methods.
`maxDegree` is a capacity requirement; the engine requests dependencies
lazily. Known differentials must be checked group homomorphisms between
those exact current source and target group objects. Missing operations
return `fail` or a record with `status:="unresolved"` and reasons.

Cells retain `page`, `degree`, `q`, `group`/`current`, `base` (E2),
`lift`, and `project`. Subsequent cells include `previous`, `cycles`,
`boundaries`, `quotientMap`, `incoming`, and `outgoing`. Each transition is
simultaneous exact homology `ker(d_r)/im(d_r)`; incoming images must be
killed by the outgoing map. Structurally zero arrows may be records with
`status:="zero"` rather than homomorphisms.

Unresolved cells have coordinates, operations, reasons, and missing
`dependencies` entries `[r,p,q,name]`. They can retain `lastKnownInvariants`
and `lastKnownPage`, but never a fabricated current group. An exact
last-known page is not a classification of an unresolved later page.
