# fermionAHSS

An exact GAP package for the five rows `q = -4,-3,-2,-1,0` of the twisted connective real K-theory Atiyah–Hirzebruch spectral sequence (AHSS), through E6, including the solution of the extension problem. It aims to calculate the classification of fermionic symmetry-protected topological (SPT) phases with various different fermionic symmetry groups up to (5+1)-dimension.

`koAHSS` returns the entries of these five rows that contribute to one
degree, and `koFull` assembles them on the final page, E6, into one group
by solving the extension problem. `koAHSS_batch` and `koFull_batch` do the
same for every degree up to a cutoff. See the [formula reference](doc/README.md),
[extension API and limits](doc/extensions.md), and
[mathematical status](doc/mathematical-status.md) for relevant mathematics.

## Installation and loading

Requirements: GAP 4.12 or newer, Polycyclic 2.16 or newer, HAP, GAP JSON,
and `python3` version 3.10 or newer on `PATH`. The bundled Python kernel uses
only the standard library; its modules and calibration JSON files live
directly in [python/](python/).

Place this directory, or a symlink to it, in a GAP `pkg` directory:

```sh
mkdir -p ~/.gap/pkg
ln -s . ~/.gap/pkg/fermionAHSS
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

## Basic use

The main function computes ko-homology, or the classification of
fermionic SPT phases, in degree `k`:

```gap
full := koFull(CyclicGroup(2), 0, 0, 3);;
full.invariants;                    # [ 0, 8 ], that is Z + Z/8
koAHSSDisplay(full);                # E6 on the line p+q=k-3, and the group
```

- `s` and `omega`: binary cocycle vectors in degrees one and two of the
  resolution. Scalar `0` denotes the zero cocycle. Coordinates depend on the
  resolution basis; a named cohomology class is not a coordinate vector.
- `k`: the degree `p+q+3`, with `-1 <= k <= 6`. The papers' spatial
  dimension is `d=k-1`.

The display shows the E6 layers of degree `k` and the group they assemble to:

```text
Line p+q=0 (degree 3)
+----+---+-----+
| q  | p | E6  |
+----+---+-----+
| 0  | 0 |  Z  |
| -1 | 1 | Z/2 |
| -2 | 2 | Z/2 |
| -3 | 3 |  0  |
| -4 | 4 | Z/2 |
+----+---+-----+
0 = zero group.
Degree 3: Z + Z/8
```

A completed result has `invariants`, an abelian invariant list; an
unresolved one has `status="unresolved"` and a `reason`, never an assumed
zero. `full.degreeResult` retains the group, measured relation vectors, one
joint integer presentation, Smith transformations, the cyclic/free basis,
and filtration maps. Relations with the same nonzero lower image are
distinguished from relations with independent lower images.

`koFull_batch` solves every degree from -1 to `k` with one E6 calculation:

```gap
batch := koFull_batch(group, s, omega, k);;
batch.degrees;                      # [-1,0,...,k]
batch.invariants;                   # one result per degree
batch.degreeResults[2+2];           # detailed degree-two result
koAHSSDisplay(batch);               # the E6 page
```

The result at degree `j` is indexed by `j+2`.

### AHSS pages

The functions

```gap
koAHSS(group, s, omega, k[, n][, options])
koAHSS_batch(group, s, omega, k[, n][, options])
```

calculate the pages of the AHSS. `koAHSS_batch` returns every entry of
degree `p+q+3 <= k`; `koAHSS` returns only the entries on the line
`p+q=k-3`, which contribute to degree `k`. Both run the same page
calculation, since the differentials into that line start on the line
below it.

- `s` and `omega`: binary cocycle vectors, as above.
- `k`: the degree `p+q+3`, with `-1 <= k <= 6`; for `koAHSS_batch` the
  largest displayed degree.
- Optional `n`: number of pages, starting at E2, with `1 <= n <= 5`.
  Omit it to compute E6 only.
- Optional `options`: `rec(details:=true)` returns a detailed result with
  labeled pages, representative maps and the retained computation context.
  The default, or `rec(details:=false)`, returns raw invariant data.
  Unknown options and nonboolean `details` values are errors.

The degree `k` and the page number are independent: `k=4, n=5`
computes E2 through E6 in degrees up to 4. E6 does not mean degree six.

`koAHSS` returns a record:

```gap
line := koAHSS(CyclicGroup(2), 0, 0, 1);;
line.pageNumbers;                   # [ 6 ]
line.lines;                         # [ [ [ 2 ], [ ], [ 2 ], fail, fail ] ]
koAHSSDisplay(line);
```

For each page in `pageNumbers`, `lines` has the five entries
`E_r^(k-q-3,q)` for `q=-4,...,0`; `fail` marks `k-q-3<0`. The display has
one column per page:

```text
Line p+q=-2 (degree 1)
+----+---+-----+
| q  | p | E6  |
+----+---+-----+
| -2 | 0 | Z/2 |
| -3 | 1 |  0  |
| -4 | 2 | Z/2 |
+----+---+-----+
0 = zero group.
```

With `n`, `pageNumbers` is `[2..n+1]`. With `rec(details:=true)` the
record also has `status`, the cell records `cells` with their
representative maps, `scope` and `certified_ko`.

`koAHSS_batch` returns one E6 table without `n`, the list of tables
E2,...,E(n+1) with `n`, and a detailed result with
`rec(details:=true)`; see [extensions.md](doc/extensions.md). Each table
has rows in order `[-4,-3,-2,-1,0]`. Entry `table[q+5][p+1]` is
`E_r^(p,q)`, for `0 <= p <= k-q-3`. Rows have lengths `max(0,k-q-2)`; an
absent position lies outside the display. At `k=6` the lengths are
`[8,7,6,5,4]`.

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
`koAHSS`, `koAHSS_batch`, `koFull` and `koFull_batch`; the same object and
twist basis are retained.
A group resolution models BG, not an arbitrary space with that fundamental
group or `B^2 Z2`.

### Stacking extensions

An existing detailed E6 calculation can be reused directly for solving its extension problem:

```gap
ahss := koAHSS_batch(CyclicGroup(2), 0, 0, 2, 5, rec(details:=true));;
full := koFull_batch(ahss);;
koAHSSDisplay(full.pages, 6);
```

`koFull_batch(ahss)` retains its labeled pages and exact backend, and
`koFull(ahss)` solves degree `ahss.maxDegree` only. Both require E6 and the
in-memory context; raw tables and earlier-page results cannot supply
extension representatives. Requesting detailed AHSS output alone does not
calculate extensions. `koAHSSDisplay(full)` selects E6 by default;
`koAHSSDisplay(full.pages)` displays all retained pages.

Higher extension equations are solved in the supplied resolution:

```gap
R := ResolutionFiniteGroup(CyclicGroup(4),6);;
full := koFull(R,[1],0,3);;
full.degreeResult.modelId;
```

The same native-resolution engine is used with a group or a detailed E6
result; no model-selection option is needed or accepted. States and linear
solves use R in degrees 3–5, and the fixed formulas are evaluated on
simplices of a comparison complex lazily: the group bar when it retracts
onto R, and otherwise a cell complex built from the generators of R, which
also covers infinite groups and several degree-zero generators. The
comparison's exact retraction identity is checked.
Gauge completeness is assumed: native gauge equivalence is taken to agree
with bar gauge equivalence. Runtime completion uses native checks and does
not construct a complete bar model for certification or fallback. See
[resolution extensions](doc/resolution-extensions.md) for eligibility,
retained checks and resource limits.

For higher layers, the engine first solves the full differential equation
for each chosen generator. An A-layer representative therefore retains
its actual B, C and D defining cochains. It changes B or C choices when a
later equation requires it, stores immutable full lifts in the resolution
basis, and uses those same lifts for every subsequent relation.
Stacking powers are identified by exact native gauge comparisons with the
ordered product of the recorded lower generators.

The native complete-state runtime covers degrees 3–5. Degree-six extensions
remain explicitly unresolved because their native formula is not implemented.
A higher-degree result is completed once every stacking relation has its
exact native comparison. Stacking is assumed to be commutative and
associative on gauge classes, so no finite multiplication table is audited
(`abelianQuotientAssumed=true`); free quotients split in the intended abelian
abutment category. For C2 the unitary degree-3 group is `[0,8]`, and the
degree-4 group with `s=omega=[1]` is `[16]`; these and the other paper
fixtures are described in [paper comparisons](doc/extension-paper-comparisons.md).

Degrees -1 to 2 use a separate low-degree adapter, which supports
degree-one C-layer relations for arbitrary valid `omega`, and degree-two
C/B relations when `omega=0`.
Incomplete page data, missing witnesses and resource limits remain
explicitly unresolved. The implementation retains `certified_ko=false`;
see [extensions.md](doc/extensions.md) for the exact scope and limits.
The papers' spatial dimension `d` corresponds to package degree `d+1`.

### Displaying the AHSS

Use the separate display function to draw the page with **q=0 at the top
and q=-4 at the bottom**, and p increasing from left to right:

```gap
pages := koAHSS_batch(CyclicGroup(2), 0, 0, 1, 5);;
koAHSSDisplay(pages);       # print E2 through E6
koAHSSDisplay(pages, 6);    # print only E6

final := koAHSS_batch(CyclicGroup(2), 0, 0, 1);;
koAHSSDisplay(final);       # a single table is labeled E6 by default
koAHSSDisplay(pages[1], 2); # label an extracted table E2

ahss := koAHSS_batch(CyclicGroup(2), 0, 0, 1, rec(details:=true));;
koAHSSDisplay(ahss);        # detailed result: same E6 text
koAHSSDisplay(ahss.pages);  # tagged page payload: same E6 text

line := koAHSS(CyclicGroup(2), 0, 0, 1, 5);;
koAHSSDisplay(line);        # the line p+q=-2 on E2 through E6
koAHSSDisplay(line, 6);     # the line on E6 only
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

A raw list of pages must start at E2, as returned by `koAHSS_batch(...,n)`.
For an extracted single page supply its page number explicitly if it is
not E6. The display shows the page groups; the raw invariant lists do not
contain differential maps, so no arrows are inferred.

Detailed results use a tagged payload
`rec(kind:="koAHSSPages",pageNumbers:=[...],tables:=[...])`.
For these inputs an optional page number selects an actual stored page;
an absent page is an error. `koAHSSDisplay` and `koAHSSFormat` also accept
this payload, a detailed AHSS result, a `koAHSS` line and the results of
`koFull` and `koFull_batch`. A `koFull_batch` result selects E6 by default
and shows no extension summary; its `.pages` payload displays every stored
page. A `koFull` result shows its E6 line followed by the group of its
degree.

### Advanced evaluation

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

All Python workers evaluate cochains through one interval-cut engine,
[cochain_tools.py](python/cochain_tools.py), and one set of shared helpers in
[phase_eval.py](python/phase_eval.py); the bundled stacking formulas import
them instead of keeping copies. The extension worker adds an exact execution
policy ([extension_acceleration.py](python/extension_acceleration.py)):
canonical structural zeros, identity memoization of pure cochain builders,
and universal source values kept across processes. The last are stored
under `$XDG_CACHE_HOME/fermionAHSS` (default `~/.cache/fermionAHSS`) in a file
keyed by the hashes of all formula sources; set `FERMIONAHSS_CACHE_DIR` to
choose another directory, or to an empty string to disable the store. A
stale or unreadable store is ignored. The values that the C2 and Z4 examples
need ship in [data/universal-values.json](data/universal-values.json) and are
loaded before the store when their recorded hashes match the sources; set
`FERMIONAHSS_BUNDLED_VALUES=0` to ignore them. After a change to the formula
sources, `python3 python/generate_universal_values.py` recomputes that file.

GAP starts one page worker ([worker.py](python/worker.py)) per session and
sends it every T batch as one JSON line, so its universal values and set-up
are computed once per session. The page worker bounds per-cochain and pure
chain-operator memo tables to 256 entries by default. Nonnegative environment
variables `KOAHSS_COCHAIN_CACHE_ENTRIES` and `KOAHSS_CHAIN_CACHE_ENTRIES`
override them; zero disables the corresponding memoization. These are entry
limits, not process-memory guarantees. Bar comparison chains and
higher-degree universal contractors can be expensive. The Python workers run
without automatic garbage collection: their memory is long-lived memo tables
with almost no reference cycles, which collection only rescanned.

## Testing

Run the package tests in a fresh GAP process:

```gap
TestPackage("fermionAHSS");
```

The Python tests run with
`python3 -m unittest discover -s python -p 'test_*.py'`. The scripts in
`examples/` run standalone from this directory, for example
`gap -q --quitonbreak examples/c2.g`; `examples/extension_papers.g`
compares `koFull_batch` with the [paper fixtures](doc/extension-paper-comparisons.md).
The package is a local development distribution; its reserved
`example.invalid` metadata URLs are placeholders, not published endpoints.
The original MIT license and attribution are preserved in [LICENSE](LICENSE).

## References

1. Robert E. Mosher and Martin C. Tangora,
   [*Cohomology Operations and Applications in Homotopy Theory*](https://store.doverpublications.com/products/9780486466644).
   Harper & Row (1968); Dover reprint (2008).
   Background on cohomology operations and their homotopy-theoretic applications.
2. Qing-Rui Wang and Zheng-Cheng Gu,
   [*Construction and classification of symmetry protected topological phases in interacting fermion systems*](https://arxiv.org/abs/1811.00536),
   *Physical Review X* **10**, 031055 (2020),
   [doi:10.1103/PhysRevX.10.031055](https://doi.org/10.1103/PhysRevX.10.031055).
   Table III supplies finite-group page comparisons;
   Table VII supplies the full invertible-phase extension samples.
3. Xing-Yu Ren, Shang-Qiang Ning, Yang Qi, Qing-Rui Wang, and Zheng-Cheng Gu,
   [*Stacking group structure of fermionic symmetry-protected topological phases*](https://arxiv.org/abs/2310.19058),
   *Physical Review B* **110**, 235117 (2024),
   [doi:10.1103/PhysRevB.110.235117](https://doi.org/10.1103/PhysRevB.110.235117).
   Background on the stacking group structure of fermionic SPT phases, the
   extension problem that `koFull` solves.
