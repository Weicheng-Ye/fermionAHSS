# fermionAHSS

An exact GAP package for the five rows `q = -4,-3,-2,-1,0` of the twisted connective real K-theory Atiyah–Hirzebruch spectral sequence (AHSS), through E6. It aims to calculate the classification of fermionic symmetry-protected topological (SPT) phases with various different fermionic symmetry groups up to (5+1)-dimension.

`koAHSS` returns the AHSS pages in these five rows. `koFull` also assembles
the supported stacking extensions and explicitly marks unsupported degrees
unresolved. The abstract extension assembler accepts arbitrary finitely
generated abelian layers and exact relation vectors; the production cochain
engine solves complete flat representatives and measures higher-layer
relations on a bounded finite group-bar model. Other coefficient rows
remain outside the scope. See the [formula reference](doc/README.md),
[extension API and limits](doc/extensions.md), and
[mathematical status](doc/mathematical-status.md).

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

`koAHSS(group, s, omega, k[, n][, options])` constructs an integral HAP resolution for
a finite group and installs the fixed calibrated operations. Its arguments
are:

- `s` and `omega`: binary cocycle vectors in degrees one and two of that
  resolution. Scalar `0` denotes the zero cocycle. Coordinates depend on the
  resolution basis; a named cohomology class is not a coordinate vector.
- `k`: largest displayed physical dimension `p+q+3`, with `-1 <= k <= 6`.
- Optional `n`: number of pages, starting at E2, with `1 <= n <= 5`.
  Omit it to return one E6 table.
- Optional `options`: `rec(details:=true)` returns a detailed result with
  labeled pages, representative maps and the retained computation context.
  The default, or `rec(details:=false)`, preserves the existing raw output.
  Unknown options and nonboolean `details` values are errors.

The dimension cutoff `k` and the page number are independent: `k=4, n=5`
computes E2 through E6 at cutoff 4. E6 does not mean dimension six.

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
A supplied resolution can also be passed directly as the first argument to
`koAHSS(R,s,omega,k[,n][,options])`; the same object and twist basis are retained.
A group resolution models BG, not an arbitrary space with that fundamental
group or `B^2 Z2`.

## Stacking extensions

The main extension entry point attempts every package degree from -1
through the requested cutoff, sharing one AHSS calculation:

```gap
full := koFull(CyclicGroup(2), 0, 0, 2);;
full.degrees;                       # [-1,0,1,2]
full.invariants;                    # one result per degree
full.degreeResults[2+2];            # detailed degree-two result
koAHSSDisplay(full);                # the same E6 table display
```

The result at degree `j` is indexed by `j+2`. A completed entry is an
abelian invariant list; an unresolved entry is a status record, never an
assumed zero. Detailed completed results retain the group, measured
relation vectors, one joint integer presentation, Smith transformations,
the cyclic/free basis, and filtration maps. Relations with the same
nonzero lower image are distinguished from relations with independent
lower images.

An existing detailed E6 calculation can be reused directly:

```gap
ahss := koAHSS(CyclicGroup(2), 0, 0, 2, 5, rec(details:=true));;
full := koFull(ahss);;
koAHSSDisplay(full.pages, 6);
```

`koFull(ahss)` retains its labeled pages and exact backend. It requires
E6 and the in-memory context; raw tables and earlier-page results cannot
supply extension representatives. Requesting detailed AHSS output alone
does not calculate extensions. `koAHSSDisplay(full)` selects E6 by default;
`koAHSSDisplay(full.pages)` displays all retained pages.

To solve the higher extension equations in a supplied resolution, select
the opt-in transferred model:

```gap
R := ResolutionFiniteGroup(CyclicGroup(4),6);;
full := koFull(R,[1],0,3,rec(extensionModel:="transfer"));;
full.degreeResults[5].modelSelection;
```

The option also works with a group or a detailed E6 result. The default
is `extensionModel:="bar"`. Transfer moves states and linear solves to R
in degrees 3–5 and evaluates the unchanged formulas on bar simplices
lazily. It checks the comparison's exact retraction identity. General
native presentations require independent complete-bar certification;
unresolved attempts fall back to the reference engine with a recorded
reason. See [resolution extensions](doc/resolution-extensions.md) for
eligibility, certificates and remaining performance limits.

For higher layers, the default engine first solves the full differential equation
for each chosen generator. An A-layer representative therefore retains
its actual B, C and D defining cochains. It changes B or C choices when a
later equation requires it, stores immutable full lifts in one common
bar basis, and uses those same lifts for every subsequent relation.
Stacking powers are identified by explicit boundary comparisons with the
ordered product of the recorded lower generators.

The complete-state runtime covers degrees 3–5. The degree-six path
additionally requires the fixed finite section described in
[the degree-six implementation](python/stacking_model/production_g6_section.py). A higher-degree result
is completed only after its bounded normal-form audit verifies all finite
products and the resulting abelian multiplication table; free quotients
split in the intended abelian abutment category. The measured C2
presentations recover `[0,8]` in unitary degree 3 and `[16]` in degree 4
with `s=omega=[1]`. Verification outcomes, including the additional audit,
are recorded in [paper comparisons](doc/extension-paper-comparisons.md).

The separate low-degree adapter retains its degree-one C-layer support
for arbitrary valid `omega`, and degree-two C/B support when `omega=0`.
Incomplete page data, missing witnesses and resource limits remain
explicitly unresolved. The implementation retains `certified_ko=false`;
see [extensions.md](doc/extensions.md) for the exact scope and limits.
The papers' spatial dimension `d` corresponds to package degree `d+1`.

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

ahss := koAHSS(CyclicGroup(2), 0, 0, 1, rec(details:=true));;
koAHSSDisplay(ahss);        # detailed result: same E6 text
koAHSSDisplay(ahss.pages);  # tagged page payload: same E6 text
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

A raw list of pages must start at E2, as returned by `koAHSS(...,n)`.
For an extracted single page supply its page number explicitly if it is
not E6. The display shows the page groups; the raw invariant lists do not
contain differential maps, so no arrows are inferred.

Detailed results use a tagged payload
`rec(kind:="koAHSSPages",pageNumbers:=[...],tables:=[...])`.
For these inputs an optional page number selects an actual stored page;
an absent page is an error. `koAHSSDisplay` and `koAHSSFormat` accept this
payload, a detailed AHSS result, or a `koFull` result. Their table text is
unchanged, and no extension summary is appended. A full result selects E6
by default; its `.pages` payload displays every stored page.

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
The notes and callbacks use the same differential names:

| Differential | Source row | Target row | Name |
| --- | ---: | ---: | --- |
| d3 | 0 | -2 | `Tau` |
| d4 | -1 | -4 | `Psi` |
| d5 | 0 | -4 | `T` |

The binary representative of `Tau` is written \(\tau'\), and the integral
representative of `Psi` is written \(\psi'\). The tertiary low selectors
form the vector \(\boldsymbol{\zeta}=(\zeta_1,\zeta_2,\zeta_3)=(0,1,0)\):
the R1 rank selector followed by the two R2 suspension selectors.
This vector is distinct from the secondary epsilon and eta vectors and
the cochain helpers \(\zeta_{i,n}\). The rational-phase helper \(\Theta_m\)
is a separate object from the differential `Psi`.

Direct secondary inputs are supported through degree seven; T inputs are
supported only in degrees zero through three. The page window needs at
most `Tau_3`, `Psi_4`, and `T_3`.

Final T has zero rank correction and `mu_R=0`, and includes
`2 beta_3,s P^1_s rho_3,s` in input degree three. Its universal helper
is fixed; no local residual solution or shortcut based on injectivity of
`Dtilde` selects it.

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

The [extension sample record](doc/extension-paper-comparisons.md) separates
literature expectations from actual computations and preserves dated
verification evidence. Its initial 2026-09-25 low-degree snapshot is
historical: the three unresolved entries there predate the complete-state
extension engine. The higher-layer follow-up supersedes that capability
assessment; consult the recorded run outcomes for completed comparisons
and any remaining limits.

## References

1. Robert E. Mosher and Martin C. Tangora,
   [*Cohomology Operations and Applications in Homotopy Theory*](https://store.doverpublications.com/products/9780486466644).
   Harper & Row (1968); Dover reprint (2008).
   Background on cohomology operations and their homotopy-theoretic applications.
2. Qing-Rui Wang and Zheng-Cheng Gu,
   [*Construction and classification of symmetry protected topological phases in interacting fermion systems*](https://arxiv.org/abs/1811.00536),
   *Physical Review X* **10**, 031055 (2020),
   [doi:10.1103/PhysRevX.10.031055](https://doi.org/10.1103/PhysRevX.10.031055).
   Table III supplies the historical finite-group page comparisons;
   Table VII supplies the full invertible-phase extension samples.
3. Shang-Qiang Ning, Xing-Yu Ren, Qing-Rui Wang, Yang Qi, and Zheng-Cheng Gu,
   [*Classification of Interacting Topological Crystalline Superconductors in Three Dimensions and Beyond*](https://arxiv.org/abs/2512.25069),
   arXiv:2512.25069 (2025).
   Source of the 230-space-group comparison tables.
4. Jian-Hao Zhang, Shang-Qiang Ning, Yang Qi, and Zheng-Cheng Gu,
   [*Construction and classification of crystalline topological superconductor and insulators in three-dimensional interacting fermion systems*](https://arxiv.org/abs/2204.13558).
   Table I supplies crystalline-superconductor extension samples, with the
   crystalline-to-internal symmetry mapping described in Section V.
