# fermionAHSS

An exact GAP package for the five rows `q = -4,-3,-2,-1,0` of the twisted connective real K-theory Atiyah–Hirzebruch spectral sequence (AHSS), through E6, including the solution of the extension problem. It aims to calculate the classification of fermionic symmetry-protected topological (SPT) phases with various different fermionic symmetry groups up to (5+1)-dimension.

`koAHSS` returns the entries of these five rows that contribute to one
degree, and `koFull` assembles the final page, E6, into one group
by solving the extension problem. `koAHSS_batch` and `koFull_batch` do the
same for every degree up to a cutoff. The mathematics is documented in
[doc/](doc/README.md): a guide to the notes with a sheet of the formulas
that the package evaluates.

## Contents

- [Installation and loading](#installation-and-loading)
- [Basic use](#basic-use): groups, pages, display, advanced evaluation
- [Examples](#examples)
- [What is computed](#what-is-computed): the window, the extension problem, the assumptions
- [Configuration](#configuration)
- [Documentation](#documentation)
- [Testing](#testing)
- [References](#references)

## Installation and loading

Requirements: GAP 4.12 or newer, Polycyclic 2.16 or newer, HAP, GAP JSON,
and `python3` version 3.10 or newer on `PATH`. The bundled Python kernel uses
only the standard library; its modules and calibration data live in
[python/](python/).

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
- `k`: the degree `p+q+3`, with `-1 <= k <= 6`, or the physics spacetime dimension. The physics spatial dimension is `d=k-1`.

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
zero. `full.degreeResult` retains the group, the measured relation vectors
with their witnesses, one joint integer presentation, Smith transformations,
the cyclic/free basis, the filtration maps, the per-prime summary `primes`
and, per prime, its generators and primary invariants `primeParts`.
Relations with the same nonzero lower image are distinguished
from relations with independent lower images.

`koFull_batch` solves every degree from -1 to `k` with one E6 calculation:

```gap
batch := koFull_batch(group, s, omega, k);;
batch.degrees;                      # [-1,0,...,k]
batch.invariants;                   # one result per degree
batch.degreeResults[2+2];           # detailed degree-two result
koAHSSDisplay(batch);               # the E6 page
```

The result at degree `j` is indexed by `j+2`.

A supplied integral HAP resolution can replace the group argument of every
function; the same object and twist basis are then retained, which also
covers infinite groups:

```gap
R := ResolutionFiniteGroup(CyclicGroup(4), 9);;
full := koFull(R, [1], 0, 6);;        # Z/4
R := ResolutionDirectProduct(ResolutionFiniteGroup(CyclicGroup(3), 9),
                             ResolutionAbelianGroup([0], 9));;
full := koFull(R, 0, 0, 6);;          # Z/3 x Z: Z/9
```

An existing detailed E6 calculation can be reused for its extension
problem, and a twist vector can be read from the cohomology data of the
resolution when a class rather than a vector is given (see
[examples/twisted_s3.g](examples/twisted_s3.g)):

```gap
ahss := koAHSS_batch(CyclicGroup(2), 0, 0, 2, 5, rec(details:=true));;
full := koFull_batch(ahss);;
koAHSSDisplay(full.pages, 6);
```

`koFull_batch(ahss)` retains its labeled pages and exact backend, and
`koFull(ahss)` solves degree `ahss.maxDegree` only. Both require E6 and the
in-memory context; raw tables and earlier-page results cannot supply
extension representatives. A group resolution models BG, not an arbitrary
space with that fundamental group or `B^2 Z2`.

### AHSS pages

```gap
koAHSS(group, s, omega, k[, n][, options])
koAHSS_batch(group, s, omega, k[, n][, options])
```

`koAHSS_batch` returns every entry of degree `p+q+3 <= k`; `koAHSS`
returns only the entries on the line `p+q=k-3`, which contribute to degree
`k`. Both run the same page calculation, since the differentials into that
line start on the line below it.

- Optional `n`: number of pages, starting at E2, with `1 <= n <= 5`. Omit
  it to compute E6 only. The degree `k` and the page number are
  independent: `k=4, n=5` computes E2 through E6 in degrees up to 4.
- Optional `options`: `rec(details:=true)` returns a detailed result with
  labeled pages, representative maps and the retained computation context.

`koAHSS` returns a record:

```gap
line := koAHSS(CyclicGroup(2), 0, 0, 1);;
line.pageNumbers;                   # [ 6 ]
line.lines;                         # [ [ [ 2 ], [ ], [ 2 ], fail, fail ] ]
koAHSSDisplay(line);
```

For each page in `pageNumbers`, `lines` has the five entries
`E_r^(k-q-3,q)` for `q=-4,...,0`; `fail` marks `k-q-3<0`. `koAHSS_batch`
returns one E6 table without `n`, the list of tables E2,...,E(n+1) with
`n`, and a detailed result with `rec(details:=true)`. Each table has rows
in order `[-4,-3,-2,-1,0]`; entry `table[q+5][p+1]` is `E_r^(p,q)`, for
`0 <= p <= k-q-3`, and rows have lengths `max(0,k-q-2)`.

Exact entries are GAP abelian invariant lists: `[]` is zero, `[0]` is Z,
`[2]` is Z/2, and `[0,2,4]` is Z + Z/2 + Z/4. An unresolved record is
neither zero nor a classification; inspect its reasons and last-known page.
Both integral rows use the sign local system Z_s.

### Displaying the AHSS

`koAHSSDisplay` draws a page with q=0 at the top, q=-4 at the bottom and p
increasing from left to right; `koAHSSFormat` returns the same text:

```gap
pages := koAHSS_batch(CyclicGroup(2), 0, 0, 1, 5);;
koAHSSDisplay(pages);       # print E2 through E6
koAHSSDisplay(pages, 6);    # print only E6
koAHSSDisplay(pages[1], 2); # label an extracted table E2
koAHSSDisplay(koAHSS(CyclicGroup(2), 0, 0, 1, 5));   # one line, E2 through E6
```

Each displayed cell uses `Z`, `Z/2`, `Z/4`, etc.; repeated factors use
powers and ` + ` denotes direct sum. `0` is the zero group, `.` is outside
the requested display window, and `?` is unresolved. These last two symbols
must not be read as zero. Both functions accept raw tables and lists (a raw
list starts at E2), the tagged page payload
`rec(kind:="koAHSSPages",pageNumbers:=[...],tables:=[...])` of a detailed
result, a detailed AHSS result, a `koAHSS` line and the results of `koFull`
and `koFull_batch`; a `koFull_batch` result selects E6 by default, its
`.pages` payload displays every stored page, and a `koFull` result shows
its E6 line followed by the group of its degree.

### Advanced evaluation

`koAHSSpages(space,s,omega,k[,n])` accepts an explicit HAP space
constructed by `koAHSSHAPSpace(R, koAHSSNaturalOperations())`, and
`koAHSSPageData` returns the cell records with groups and representative
maps; see [backends.md](doc/backends.md). Sufficient resolution lengths are
`max(3,k+2)` through E3 and `max(3,k+3)` through E4–E6; the group wrapper
selects them. The operations can be evaluated directly on cochains, with
`backend := space.koAHSS(s,omega,maxDegree)`:

- `koAHSSNaturalSecondary(backend, degree, A)` evaluates integral-input Tau;
- `koAHSSNaturalSecondary(backend, degree, a, rec(inputType := "mod2"))`
  evaluates mod-two-input Psi;
- `koAHSSNaturalTertiary(backend, degree, A[, options])` evaluates final T.

Here `A` is a signed integral cocycle and `a` is binary; the resolution
must reach `degree+5` for the secondary call and `degree+6` for the
tertiary one. Secondary inputs are supported through degree seven, T inputs
in degrees zero through three; the page window needs at most `Tau_3`,
`Psi_4` and `T_3`.

## Examples

The scripts in [examples/](examples/) run standalone from the package root,
`gap -q --quitonbreak examples/NAME.g`, or after installation with
`ReadPackage("fermionAHSS", "examples/NAME.g")`. Each states what it shows
and its running time.

| Script | Content | Time |
| --- | --- | --- |
| [c2.g](examples/c2.g) | Untwisted C2: the pages E2 through E6 in degrees up to 1 | seconds |
| [twisted_c2.g](examples/twisted_c2.g) | C2 with `s=omega=[1]`; reusing a chosen resolution for the twists | seconds |
| [c4_signed.g](examples/c4_signed.g) | C4 with `s=[1]`, every degree from -1 to 6 from one E6 calculation on a supplied resolution; Z/4 at degree 6 | under a minute |
| [odd_torsion.g](examples/odd_torsion.g) | Degree 5 of Z/3, Z/9, Z/3×Z/3, Z/5 and Z/6: the relations at the prime three and their splitting at five | under a minute |
| [suspension.g](examples/suspension.g) | The infinite group Z/3 × Z on a product resolution; degree six at the prime three gives Z/9 by suspension | a minute |
| [pin_minus.g](examples/pin_minus.g) | Pin⁻ bordism in five spatial dimensions: C2 with `s=[1]` at degree 6 gives Z/16 through the target layers | three minutes |
| [twisted_s3.g](examples/twisted_s3.g) | S3 with the sign twist read from H¹ of its resolution; degree 5 gives Z/9 | three minutes |
| [c4_omega.g](examples/c4_omega.g) | C4 with `omega=[1]` at degree 5: the complete degree-five stacking correction gives Z/2 + Z/32 | 35 minutes |
| [extension_papers.g](examples/extension_papers.g) | Comparison with the literature fixtures of [data/extension-paper-samples.json](data/extension-paper-samples.json) | five minutes by default |

## What is computed

The page calculation covers the rows `q=-4,-3,-2,-1,0` through E6 for
`-1 <= k <= 6`. The differentials `Tau` (d3 from row 0 to -2), `Psi` (d4
from row -1 to -4) and `T` (d5 from row 0 to -4) are fixed cochain formulas
with the helper family `chi7_tail`, evaluated on the group bar resolution
and transferred to the supplied resolution. Primary page arrows and the
proven low-degree page-class reductions are evaluated directly on that
resolution. Their explicit forms are in the
[formula sheet](doc/README.md#formula-sheet) and [backend interface](doc/backends.md).

`koFull` assembles the E6 layers `A=(k-3,0)`, `B=(k-2,-1)`, `C=(k-1,-2)`
and `D=(k+1,-4)` of a degree in the native stacking model on the supplied
resolution. It groups the layer generators by the prime of their order and
records, prime by prime and from D upward, the relation `m*g=t` of each
generator `g` of finite order `m` in the lower group `H`. The row `t` is
zero without a measurement when the relation splits at its prime (the
primes five and above, and three below degree five), when `H` has no
generator of that prime and no free generator, or when no lower generator
lies outside `m*H`. Otherwise the relation is determined through its
target layer, the lowest layer of the lower presentation with a generator
outside `m*H` (through D when a later relation that refers to `g` is
measured below the layer of `g`). When that layer lies right below the
layer of `g` and the relation is two-primary, its row is the class of a
primary operation of the cocycle of `g`, evaluated with the cup-i products
of the resolution (`D(z)` with `beta_s z=(m/2)[a]` for A over B,
`(Sq^1+s)b` for B over C, an integral lift of `D(c)` for C over D).
The remaining relations use [light residues](doc/extensions.md#light-rows)
of transported defining data, with no flat-lift search or reflected gauge
comparison. Versioned markings and a final frame audit keep their lower
references consistent. The degree-six A-over-D residue uses a checked table
of 22 unary constants instead of the upper two-A pair contractor. The rows
of all primes form one integer presentation, reduced to Smith form.
Workers and comparison chains are built on demand. Degrees -1 and 0 have
only the D layer.

`FERMIONAHSS_LIGHT_RELATIONS=0` restores model measurement: a flat lift is
stacked and compared with the lower product by an exact native gauge.
A non-resource failure of a light part also restarts that part in the
model, recorded in `lightFallbacks`. `heavyMeasurements` counts entries
into this measurement path. Its degree-six correction of two nonzero A
states requires `FERMIONAHSS_DEGREE_SIX_A_STACKING=1`.

Light completion records `certificateLevel="light-R"` and
`gaugeCompletenessAssumed=false`; native model measurement assumes gauge
completeness. Both assume that stacking is commutative and associative on
gauge classes (`abelianQuotientAssumed`); unresolved
searches, missing witnesses and resource limits stay unresolved and are
never reported as zero, and every result keeps `certified_ko=false`. The
exact scope and limits are in [extensions.md](doc/extensions.md),
[resolution-extensions.md](doc/resolution-extensions.md) and
[mathematical-status.md](doc/mathematical-status.md); the values found for
the literature are in
[extension-paper-comparisons.md](doc/extension-paper-comparisons.md).

## Configuration

GAP starts one Python worker for the page differentials and one for the
extension problem; both keep their universal source values in a store and
load the values bundled in
[data/universal-values.json](data/universal-values.json) when its recorded
hash matches the formula sources.
Each table's stored keys are decoded on first use, so low-degree requests
do not initialize unused higher-phase tables. The stored values and their
precedence are unchanged.

| Variable | Effect |
| --- | --- |
| `FERMIONAHSS_CACHE_DIR` | Directory of the universal-value store (default `$XDG_CACHE_HOME/fermionAHSS` or `~/.cache/fermionAHSS`); empty disables it |
| `FERMIONAHSS_BUNDLED_VALUES=0` | Ignore the bundled universal values |
| `FERMIONAHSS_NATIVE_PAGES=0` | Evaluate all page arrows through the fixed bar comparison instead of using native primary and proven low-degree page-class reductions |
| `FERMIONAHSS_PRIMARY_TENSOR_REFERENCE=1` | Use the integral tensor implementation reduced modulo two for primary-comparison checks |
| `FERMIONAHSS_LAYERED_RELATIONS=0` | Measure relations through the D layer instead of their target layer |
| `FERMIONAHSS_NATIVE_RELATIONS=0` | Measure the relations whose target layer lies right below the generator's in the model instead of reading them from a primary operation on R |
| `FERMIONAHSS_LIGHT_RELATIONS=0` | Measure the other relations in the transferred model instead of reading them from light rows (doc/extensions.md, "Light rows") |
| `FERMIONAHSS_LIGHT_ABSORPTION=0` | Disable absorption of B-row D components into the C markings |
| `FERMIONAHSS_PRIME_LOCAL=0` | Do not localize at the primes: relations of every prime use the complete model and keep every coordinate of their rows |
| `FERMIONAHSS_DEGREE_SIX_A_STACKING=1` | Evaluate the degree-six correction of two nonzero A layers (hours to days) |
| `KOAHSS_COCHAIN_CACHE_ENTRIES`, `KOAHSS_CHAIN_CACHE_ENTRIES` | Memo-table bounds of the page worker (default 256; zero disables) |

After a change to the formula sources, `python3 python/generate_universal_values.py`
recomputes the bundled values.

## Documentation

The notes in [doc/](doc/README.md) are organized as follows; the guide
there has a formula sheet with the explicit differentials, the stacking
model and the extension relations.

| Notes | Content |
| --- | --- |
| [conventions.md](doc/conventions.md) | Twists, coefficient systems, cup-i words, lifts and the primary maps |
| [secondary_operations.md](doc/secondary_operations.md), [tertiary_operations.md](doc/tertiary_operations.md) | The differentials `Tau`, `Psi` and `T` |
| [universal_helpers.md](doc/universal_helpers.md), [universal_value_growth.md](doc/universal_value_growth.md) | The fixed helpers and the stored universal values |
| [extensions.md](doc/extensions.md) | `koFull`: layers, relations, target layers, localization at the primes, Smith bases |
| [low_degree_stacking.md](doc/low_degree_stacking.md), [dimension_indexed_differentials.md](doc/dimension_indexed_differentials.md), [all_cochain_differential.md](doc/all_cochain_differential.md) | The stacking model: products and the differential in every degree |
| [transfer.md](doc/transfer.md), [resolution-extensions.md](doc/resolution-extensions.md) | Evaluating the model on a supplied resolution |
| [backends.md](doc/backends.md) | The HAP, cochain and page interfaces |
| [extension-paper-comparisons.md](doc/extension-paper-comparisons.md) | The literature fixtures and the package's results |
| [mathematical-status.md](doc/mathematical-status.md) | Implemented range, assumptions, sources and open scope |

## Testing

Run the package tests in a fresh GAP process:

```gap
TestPackage("fermionAHSS");
```

The Python tests run with
`python3 -m unittest discover -s python -p 'test_*.py'`. The degree-six
comparison with the complete-bar section model and the suspension check,
`tst/extension_degree_six.tst`, is not part of `TestPackage`; run it with
`Test` on that file (about six minutes), and set
`FERMIONAHSS_SLOW_TRANSFER_TESTS=1` for the corresponding opt-in Python
tests. The [examples](#examples) run standalone from this directory.
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
   Table III supplies the finite-group comparisons and Table VII the full
   invertible-phase groups of the C2 cases.
3. Xing-Yu Ren, Shang-Qiang Ning, Yang Qi, Qing-Rui Wang, and Zheng-Cheng Gu,
   [*Stacking group structure of fermionic symmetry-protected topological phases*](https://arxiv.org/abs/2310.19058),
   *Physical Review B* **110**, 235117 (2024),
   [doi:10.1103/PhysRevB.110.235117](https://doi.org/10.1103/PhysRevB.110.235117).
   The stacking group structure of fermionic SPT phases, the extension
   problem that `koFull` solves.
4. Zhang, Ning, Qi, and Gu, [arXiv:2204.13558](https://arxiv.org/abs/2204.13558).
   Three-dimensional crystalline topological superconductors; Table I
   supplies the point-group comparisons.
5. [arXiv:2512.25069](https://arxiv.org/abs/2512.25069). Tables I and II
   supply the E6 layers of the 230 space groups.
