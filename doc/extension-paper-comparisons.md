# Comparisons with the literature

The fixture [extension-paper-samples.json](../data/extension-paper-samples.json)
records a few printed classifications from three papers, together with the
inputs that reproduce their symmetry data here, and the package's results
for them are discussed below. The sources are Wang–Gu Tables VII and III
([arXiv:1811.00536](https://arxiv.org/abs/1811.00536)), Zhang–Ning–Qi–Gu
Table I ([arXiv:2204.13558](https://arxiv.org/abs/2204.13558)) and Tables
I–II of [arXiv:2512.25069](https://arxiv.org/abs/2512.25069). The expected
values are literature fixtures, not inputs to the solver, and every result
keeps `certified_ko=false`. The wallpaper-group tables of Ren–Ning–Qi–Wang–Gu
([arXiv:2310.19058](https://arxiv.org/abs/2310.19058)) are not compared:
`SGC_ResolutionSpaceGroup` does not support two-dimensional groups.

## Conventions

The papers' spatial dimension is `d=p+q+2` and the package degree is
`k=d+1`: paper dimension `d` is `koFull(group,s,omega,d+1).invariants`. The
group argument is the bosonic quotient `G_b=G_f/Z2^f`; `omega` is the
extension class of `G_f` and `s` its antiunitary or orientation character,
both as cocycle vectors in the basis of the resolution that
`group_expression` constructs (`ResolutionFiniteGroup` of length `cutoff+3`
for a group). Every case has the same fields:

```json
{
  "id": "wang-gu-C2-antiunitary-nonsplit",
  "source": "wang_gu_2020",
  "fermionic_symmetry": "Z4^Tf",
  "group_expression": "CyclicGroup(2)",
  "s": [1], "omega": [1], "cutoff": 4,
  "quantity": "full group",
  "expected_by_package_degree": {"1": [], "2": [2], "3": [2], "4": [16]}
}
```

`quantity` says what the expected values are. The **full group** is the
stacking group of the whole E6 line. The **group of the E6 layers with
p ≥ 1** leaves out the layer at `p=0`, the intrinsic Kitaev layer `(0,-1)`
at `k=2` and the intrinsic `p+ip` layer `(0,0)` at `k=3`, which Table III
does not count; in `koFull` this is the filtration
stage below that layer (`degreeResult.filtration`), and for `k\ge4` it is the
full group. The **E6 layers** are the cells of the line themselves, for a
paper that tabulates layers without solving the extensions. Invariant lists
are compared as multisets, in GAP's order (`[]` zero, `[0]` Z, then torsion
orders).

```sh
gap -q --quitonbreak examples/extension_papers.g
FERMIONAHSS_PAPER_CASES=wang-gu-iii gap -q --quitonbreak examples/extension_papers.g
```

The runner computes `koFull_batch` for every case, or for the cases whose
id starts with one of the listed prefixes, and only then reads the expected
values; it prints a JSON report and a summary of matches, mismatches,
unresolved and skipped cases, and asserts only that the four Table VII
groups are reproduced. The space groups are resolved by
`SGC_ResolutionSpaceGroup` of the SpaceGroupCohomology package, whose
resolutions carry the contracting homotopy that `koFull` needs (HAP's own
space-group resolutions have none); the runner loads the package when it is
installed and skips the space groups without it. The finite groups take
about five minutes altogether.

## Wang and Gu, Table VII: the four C2 cases

Table VII (printed page 67) gives the full invertible-phase groups of the
fermionic symmetry groups over `G_b=Z2`, with the twists of equations
(E23)–(E24):

| Fermionic symmetry | `s` | `omega` | Degree 1 (`d=0`) | Degree 2 (`d=1`) | Degree 3 (`d=2`) | Degree 4 (`d=3`) |
| --- | --- | --- | --- | --- | --- | --- |
| Z2^f × Z2 | `0` | `0` | `[2,2]` | `[2,2]` | `[0,8]` | `[]` |
| Z4^f | `0` | `[1]` | `[4]` | `[]` | `[0]` | `[]` |
| Z2^f × Z2^T | `[1]` | `0` | `[2]` | `[8]` | `[]` | `[]` |
| Z4^Tf | `[1]` | `[1]` | `[]` | `[2]` | `[2]` | `[16]` |

All sixteen entries are reproduced. They are extension results, not
layer counts: the standard C2 resolution has rank one in each degree, so
the E6 layers are single Z/2 or Z summands and the measured relations, in
the marked D/C/B/A bases, are

| Input `(s,omega)` | Degree | Relations | Group |
| --- | ---: | --- | --- |
| `(0,[1])` | 1 | `2D=0`, `2C=D` | Z/4 |
| `([1],0)` | 2 | `2D=0`, `2C=D`, `2B=C+2D` | Z/8 |
| `(0,0)` | 3 | `2D=0`, `2C=D`, `2B=C+D`; A free | Z ⊕ Z/8 |
| `([1],[1])` | 4 | `2D=0`, `2C=D`, `2B=C`, `2A=B` | Z/16 |

The free A generator at degree 3 is `H^0(BZ2;Z)=Z`, the `p+ip` layer; each
relation of degree 4 is measured in the target layer just below its
generator, so the joint presentation is the bidiagonal matrix with `2` on
the diagonal and `-1` below it, whose cokernel is Z/16.

## Wang and Gu, Table III: finite unitary groups

Table III (printed page 8) lists the classifications of finite unitary
groups in dimensions one to three without the intrinsic Kitaev and `p+ip`
layers. Four rows are in the fixture, with `G_b` the direct product of the
cyclic factors and `omega` the extension class of `G_f` (zero for the split
products, the cyclic carry on a `Z2` factor for `Z4^f × Z2`).

| Case | `G_f` | `s`, `omega` | Degree 2 | Degree 3 | Degree 4 |
| --- | --- | --- | --- | --- | --- |
| row 3 | Z2^f × Z4 | `0`, `0` | `[2]` | `[2,8]` | `[]` |
| row 5 | Z2^f × Z2 × Z4 | `0`, `0` | `[2,2,2]` | `[2,2,2,8,8]` | `[2,4]` |
| row 13 | Z4^f × Z2 | `0`, carry | `[2]` | `[4]` | `[2]` |
| row 9 | Z2^f × Z2 × Z2 × Z2 | `0`, `0` | `[2,2,2,2]` | `[2,4,4,4,8,8,8]` | `[2]^8` |

Rows 3, 5 and 13 are reproduced in every degree. For Z2^f × Z4 at degree
3 the layers with `p\ge1` are `B=H^1(Z4;Z2)=Z/2`, `C=H^2(Z4;Z2)=Z/2` and
`D=H^4(Z4;Z)=Z/4`, and the measured relations `4D=0`, `2C=D`, `2B=6D`
give Z/2 ⊕ Z/8, the printed `Z8 × Z2`; the full group adds the free `p+ip`
class, `[0,2,8]`. For Z4^f × Z2 at degree 3 the relations `2D=0`, `2C=D`
give the printed Z/4, and for Z2^f × Z2 × Z4 at degree 4 the three
generators `D1`, `D2` (order two) and `C` obey `2C=3D1+4D2=D1`, so the
group is Z/2 ⊕ Z/4 as printed.

Row 9 agrees at degrees 3 and 4 but not at degree 2, where the paper prints
`Z2^4`. The E6 layers with `p\ge1` of `G_b=Z2^3` at degree 2 are
`C=H^1(Z2^3;Z2)=Z2^3` and `D=H^3(Z2^3;Z)=Z2^3`; both survive to E6, the
relation matrix is diagonal, and their stacking group is `Z2^6`, of order
64. A group of order 16 cannot be the stacking group of these layers, so
the discrepancy is in the counting of one-dimensional layers, not in an
extension; the same happens for rows 10 and 11.

## Zhang, Ning, Qi and Gu, Table I: point groups

Table I (printed page 4) gives the groups of three-dimensional crystalline
topological superconductors for the 32 point groups at package degree 4.
The point group acts on space by the integral matrices of
`group_expression`, `s` is the determinant character, and `omega` is zero
for spin-1/2 fermions and the second Stiefel–Whitney class `w2` of the
spatial representation for spinless ones.

| Case | Matrices | `s`, `omega` | Printed | Computed |
| --- | --- | --- | --- | --- |
| Cs, spinless | one reflection | `[1]`, `[1]` | `[16]` | `[16]` |
| Cs, spin-1/2 | one reflection | `[1]`, `0` | `[]` | `[]` |
| D2, spin-1/2 | two rotations by π | `0`, `0` | `[2,2]` | `[2,2]` |
| S4, spin-1/2 | one four-fold rotoreflection | `[1]`, `0` | `[2,2]` | `[4]` |

The reflection group Cs is abstractly C2 with `s` the sign character, so
the spinless case has the inputs of the Table VII case Z4^Tf and reproduces
its Z/16 with the same four relations, while the spin-1/2 case has the
inputs of Z2^f × Z2^T and is zero at degree 4. For D2 the E6 line has only
the bosonic layer `H^5(D2;Z)=Z2^2`, so there is no extension to measure.
For S4, the cyclic group of order four generated by a rotoreflection, the
E6 line has `A=H^1(Z4;Z_s)=Z/2` and `C=H^3(Z4;Z2)=Z/2`, and the measured
relation `2A=C` makes the group Z/4, where Table I prints `Z2^2`; the
paper's Supplement S-3.8 derives the same nontrivial extension to Z/4, and
the fixture records this conflict in `source_conflicts_not_tested`.

## arXiv:2512.25069, Tables I–II: space groups

Tables I and II give the E6 layers `p+ip (1,0)`, `Kitaev (2,-1)`,
`complex fermion (3,-2)` and `bosonic (5,-4)` of the spin-1/2 electronic
case (`s=w1`, `omega=0`) of the 230 space groups at package degree 4,
without their extensions. Two groups are in the fixture, with the
resolutions of the SpaceGroupCohomology package.

| Group | Printed layers (p+ip; Kitaev; complex; bosonic) | Computed |
| --- | --- | --- |
| SG 2, P-1 | `Z^3`; `0`; `Z2^4`; `Z2` | the same |
| SG 7, Pc | `Z`; `Z2^3`; `Z2`; `0` | `Z ⊕ Z2`; `Z2^3`; `Z2`; `0` |

For Pc the package finds an extra Z/2 in the `p+ip` layer: the torsion
of `H^1(Pc;Z_s)`, generated by the orientation-reversing glide, survives
to E6 with the twist `s=w1`. Whether this class belongs to the paper's
classification depends on how the paper treats `omega+s^2`; the fixture
records the printed value and does not settle the question. The same
extra Z/2 appears for two dozen other space groups, and the other
layers agree for all but five groups.

## Interpreting results

Record the actual invariant list and the resolution and twist inputs of a
case, or its explicit unresolved status. An unresolved result is not zero
and cannot count as a match. Compare complete invariant factors, including
the free rank, rather than products of layer orders, and keep
`certified_ko=false`: these bounded comparisons do not establish
all-degree naturality or a complete identification of `ko` with the
classification of fermionic phases.
