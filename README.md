# fermionAHSS

An exact GAP package for the five rows `q = -4,-3,-2,-1,0` of the twisted
connective real K-theory Atiyah–Hirzebruch spectral sequence, through E6.
It computes actual kernels, images, torsion, and representative maps. The
calibrated secondary and degree-zero-through-three tertiary operations use
a fixed normalized group-bar comparison with the supplied HAP resolution.
No transport-model setting is needed.

The output is the five-row associated-graded calculation, not the total
ko group. Other coefficient rows and abutment extensions are outside its
scope; saved calculations retain `certified_ko: false`. Unavailable
operations give unresolved entries, while invalid cochains and identities
remain errors. See [mathematical status](doc/mathematical-status.md).

## Installation and loading

Requirements: GAP 4.12 or newer, Polycyclic 2.16 or newer, HAP, GAP JSON,
and `python3` version 3.10 or newer on `PATH`. The bundled Python kernel uses
only the standard library; its modules and calibration JSON files live
directly in [python/](python/).

Place this directory, or a symlink to it, in a GAP `pkg` directory:

```sh
mkdir -p ~/.gap/pkg
ln -s ~/Workspace/fermionAHSS ~/.gap/pkg/fermionAHSS
```

Then use a fresh GAP session:

```gap
LoadPackage("fermionAHSS");
```

GAP's `ReadPackage` takes a package name **and a filename**. The supported
explicit reading form is:

```gap
ReadPackage("fermionAHSS", "load.g");
```

For an uninstalled checkout, use:

```gap
Read("/absolute/path/to/fermionAHSS/load.g");
```

Loading again in the same session is harmless, but does not reload edited
source. This package retains the `koAHSS...` API names and must not share
a GAP session with the separate koAHSS package.

## Basic use

```gap
G := CyclicGroup(2);;
table := koAHSS(G, 0, 0, 1);;       # one E6 table
pages := koAHSS(G, 0, 0, 1, 5);;    # E2, E3, E4, E5, E6
```

`koAHSS(group, s, omega, k[, n])` constructs an integral HAP resolution for
a finite group and installs the fixed calibrated operations. Its arguments
are:

- `s` and `omega`: binary cocycle vectors in degrees one and two of that
  resolution. Scalar `0` denotes the zero cocycle. Coordinates depend on the
  resolution basis; a named cohomology class is not a coordinate vector.
- `k`: largest displayed physical dimension `p+q+3`, with `-1 <= k <= 6`.
- Optional `n`: number of pages, starting at E2, with `1 <= n <= 5`.
  Omit it to return one E6 table.

Each table has rows in order `[-4,-3,-2,-1,0]`. Entry
`table[q+5][p+1]` is `E_r^(p,q)`, for `0 <= p <= k-q-3`.
Rows have lengths `max(0,k-q-2)`; an absent position lies outside the
display. At `k=6` the lengths are `[8,7,6,5,4]`.

Exact entries are GAP abelian invariant lists: `[]` is zero, `[0]` is Z,
`[2]` is Z/2, and `[0,2,4]` is Z + Z/2 + Z/4. An unresolved record is
neither zero nor a classification; inspect its reasons and last-known page.
Both integral rows use the sign local system Z_s. Hidden outgoing targets
are retained even when they are outside the displayed table.

Run the simple examples after installation:

```gap
ReadPackage("fermionAHSS", "examples/c2.g");
ReadPackage("fermionAHSS", "examples/twisted_c2.g");
```

See [c2.g](examples/c2.g) and [twisted_c2.g](examples/twisted_c2.g).
For twists with a known resolution basis, explicitly retain that resolution:

```gap
R := ResolutionFiniteGroup(CyclicGroup(2), 4);;
space := koAHSSHAPSpace(R, koAHSSNaturalOperations());;
pages := koAHSSpages(space, [1], [1], 1, 5);;
```

For general groups, establish the basis before choosing nonzero twist vectors.
A group resolution models BG, not an arbitrary space with that fundamental
group or `B^2 Z2`.

## Displaying the AHSS

Use the separate display function to draw the page with **q=0 at the top
and q=-4 at the bottom**, and p increasing from left to right:

```gap
pages := koAHSS(CyclicGroup(2), 0, 0, 1, 5);;
koAHSSDisplay(pages);       # print E2 through E6
koAHSSDisplay(pages, 6);    # print only E6

final := koAHSS(CyclicGroup(2), 0, 0, 1);;
koAHSSDisplay(final);       # a single table is labeled E6 by default
koAHSSDisplay(pages[1], 2); # label an extracted table E2
```

Each displayed cell uses `Z`, `Z/2`, `Z/4`, etc.; repeated factors use
powers and ` + ` denotes direct sum. `0` is the zero group, `.` is outside
the requested display window, and `?` is unresolved. These last two symbols
must not be read as zero. Raw data keep their original q=-4 through q=0
row order; formatting never changes or recomputes them.

`koAHSSFormat` accepts the same arguments and returns the formatted string
without printing it, for example:

```gap
out := OutputTextFile("C2-E6.txt", false);;
SetPrintFormattingStatus(out, false);
PrintTo(out, koAHSSFormat(pages, 6));
CloseStream(out);
```

A list of pages must start at E2, as returned by `koAHSS(...,n)`.
For an extracted single page supply its page number explicitly if it is
not E6. The display shows the page groups; the raw invariant lists do not
contain differential maps, so no arrows are inferred.

## Advanced evaluation

`koAHSSpages(space,s,omega,k[,n])` accepts an explicit HAP space or other
supported backend. `koAHSSPageData` has the same arguments and returns cell
records with groups and representative lift/projection maps. See
[backend interfaces](doc/backends.md).

Sufficient HAP resolution lengths are `max(3,k+2)` through E3 and
`max(3,k+3)` through E4–E6. The group wrapper selects these lengths.
A short resolution is an error. A bare resolution passed to `koAHSSpages`
does not automatically install the calibrated higher callbacks; use the
explicit `koAHSSHAPSpace` factory call above.

Direct operations use cohomological input degree `degree`, independent of
the page-count argument:

- `koAHSSNaturalSecondary(backend, degree, A)` evaluates integral-input Tau.
- `koAHSSNaturalSecondary(backend, degree, a, rec(inputType := "mod2"))`
  evaluates mod-two-input Psi.
- `koAHSSNaturalTertiary(backend, degree, A[, options])` evaluates final T.

These calls require sufficient resolution depth for their full identities:
through `degree+5` for the secondary call and through `degree+6` for the
tertiary defining system. Here `A` is a signed integral cocycle and `a` is
binary. Construct `backend := space.koAHSS(s,omega,maxDegree)` with that
capacity. Optional `b` and `c` fields supply defining cochains, which are
checked by exact coboundary equations.

## Formulas and conventions

The [formula reference](doc/README.md) records the complete implemented
secondary and tertiary expressions, chi-word encoding, lifts, signs,
calibration coefficients, and source maps. The lower helper is `chi7_tail`,
with secondary coefficients epsilon `(1,0,0)` and eta `(1,0,1)`.
Direct secondary inputs are supported through degree seven; T inputs are
supported only in degrees zero through three. The page window needs at
most `Tau_3`, `Psi_4`, and `T_3`.

Final T uses the Danus reference with zero rank correction and `mu_R=0`,
plus `2 beta_3,s P^1_s rho_3,s` in input degree three. Its universal helper
is fixed; no local residual solution or injective-J shortcut selects it.
The coefficient-one correction belongs only to the separate legacy
`TReference` interface and is not applied to this final T.

The Python worker bounds per-cochain and pure chain-operator memo tables
to 256 entries by default. Nonnegative environment variables
`KOAHSS_COCHAIN_CACHE_ENTRIES` and `KOAHSS_CHAIN_CACHE_ENTRIES` override
them; zero disables the corresponding memoization. These are entry limits,
not process-memory guarantees. Bar comparison chains and higher-degree
universal contractors can be expensive.

## Verification and development status

Run the bundled smoke check in a fresh GAP process:

```gap
TestPackage("fermionAHSS");
```

The two scripts in `examples/` also run standalone with
`gap -q --quitonbreak examples/c2.g` and
`gap -q --quitonbreak examples/twisted_c2.g` from this directory.
The package is a local development distribution; its reserved
`example.invalid` metadata URLs are placeholders, not published endpoints.
The original MIT license and attribution are preserved in [LICENSE](LICENSE).

The [verification record](doc/verification.json) and [logs](doc/verification/)
record the checks run on this distribution.
