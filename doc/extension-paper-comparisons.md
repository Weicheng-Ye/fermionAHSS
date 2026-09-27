# Full-group comparison fixtures

The fixture [extension-paper-samples.json](../data/extension-paper-samples.json)
records full invertible-phase groups from the literature, with their
inputs. There are four distinct C2 calculations; the crystalline
comparisons reuse their degree-four results. The expected groups below are
literature fixtures, not inputs to the solver.

The papers' spatial dimension is `d=p+q+2`, while the package degree is
`k=p+q+3=d+1`. Thus paper dimensions 0–3 correspond to package degrees
1–4: paper dimension `d` is `koFull(CyclicGroup(2),s,omega,d+1).invariants`,
or `full.invariants[d+3]` for `full := koFull_batch(CyclicGroup(2),s,omega,4)`.
The group argument is the bosonic quotient `CyclicGroup(2)`. Nonzero
twists `[1]` refer to the standard rank-one HAP cyclic resolution.

## Internal symmetry

[Wang–Gu, arXiv:1811.00536v3](https://arxiv.org/abs/1811.00536v3),
[Table VII, printed page 67](https://arxiv.org/pdf/1811.00536#page=67),
gives full invertible-phase groups. Tables III and VI instead omit
intrinsic Kitaev and `p+ip` phases in dimensions one and two. Use Table VII
when retaining all the package's filtration layers. Appendix E.2,
equations (E23)–(E24), specifies the nonzero twists.

| Fermionic symmetry | `s` | `omega` | Degree 1 (`d=0`) | Degree 2 (`d=1`) | Degree 3 (`d=2`) | Degree 4 (`d=3`) |
| --- | --- | --- | --- | --- | --- | --- |
| Z2^f × Z2 | `0` | `0` | `[2,2]` | `[2,2]` | `[0,8]` | `[]` |
| Z4^f | `0` | `[1]` | `[4]` | `[]` | `[0]` | `[]` |
| Z2^f × Z2^T | `[1]` | `0` | `[2]` | `[8]` | `[]` | `[]` |
| Z4^Tf | `[1]` | `[1]` | `[]` | `[2]` | `[2]` | `[16]` |

Here `[]` means the zero group, `[0]` means Z, and `[0,8]` means
Z ⊕ Z/8. These expected groups test extensions, not just phase counts.

## Crystalline symmetry

[Zhang–Ning–Qi–Gu, arXiv:2204.13558v2](https://arxiv.org/abs/2204.13558v2),
[Table I, printed page 4](https://arxiv.org/pdf/2204.13558#page=4),
gives crystalline topological-superconductor groups. The mapping in
[Section V, printed page 47](https://arxiv.org/pdf/2204.13558#page=47)
turns reflections into antiunitary symmetry and exchanges spinless and
spin-1/2 conventions. The resulting internal inputs are:

| Crystalline case | `s` | `omega` | Expected at package degree 4 (`d=3`) |
| --- | --- | --- | --- |
| Reflection Cs=C1h, spinless | `[1]` | `[1]` | `[16]` |
| Reflection Cs=C1h, spin-1/2 | `[1]` | `0` | `[]` |
| Rotation C2, spinless | `0` | `[1]` | `[]` |
| Rotation C2, spin-1/2 | `0` | `0` | `[]` |

[Supplement S-3.1–S-3.2, printed page 54](https://arxiv.org/pdf/2204.13558#page=54)
discusses these cases; the mirror result includes a nontrivial extension
between an order-two decoration and Z/8. The cases reuse the internal
C2 calculations above. Table II concerns charge-conserving insulators
and is not used here.

There is a source conflict outside this sample: spin-1/2 four-fold
rotoreflection S4 is Z/2 ⊕ Z/2 in Table I but Z/4 in the extension
calculation in [Supplement S-3.8, printed page 63](https://arxiv.org/pdf/2204.13558#page=63).
The fixture preserves both statements without selecting one and runs no
C4 case. This S4 is abstractly cyclic C4, not the symmetric group.

## Running the comparison

```sh
gap -q --quitonbreak examples/extension_papers.g
```

The runner computes each of the four distinct C2 cases with
`koFull_batch(CyclicGroup(2),s,omega,4)`, independently of the fixture,
and only then compares the actual invariant lists with the expected ones.
It prints a JSON report with the invariants, relation matrices, fixed flat
lifts and power witnesses of every labeled entry, and asserts that all 20
entries are computed and match. The crystalline entries reuse the
corresponding internal calculation, so these are 20 labeled comparisons,
not 20 independent runs.

## Measured relations

The groups come from measured stacking relations followed by Smith
reduction, not from classifications supplied to the solver. In the
marked D/C/B/A bases:

| Input `(s,omega)` | Degree | Relations | Group |
| --- | ---: | --- | --- |
| `(0,[1])` | 1 | `2D=0`, `2C=D` | Z/4 |
| `([1],0)` | 2 | `2D=0`, `2C=D`, `2B=C` | Z/8 |
| `(0,0)` | 3 | `2D=0`, `2C=D`, `2B=C+D`; A free | Z ⊕ Z/8 |
| `([1],[1])` | 4 | `2D=0`, `2C=D`, `2B=C+7D`, `2A=B+C-17D` | Z/16 |

For C2 with `s=omega=[1]` in package degree 4, all coordinates refer to
the same stored D/C/B/A columns. The joint presentation and its nonunit
Smith factor are

\[
R=\begin{pmatrix}
2&0&0&0\\
-1&2&0&0\\
-7&-1&2&0\\
17&-1&-1&2
\end{pmatrix},\qquad \operatorname{coker}(R)\cong\mathbb Z/16.
\]

The standard C2 resolution has rank one in each degree, so native states
are integer vectors `(A,B,C,D)`. The selected A lift is `(1,0,0,0)`, and
its stacked double is `(2,0,0,-42)`. The canonical lower product for the
measured coordinate vector `[-17,1,1]` is `(0,1,1,-18)`, and the
preceding-degree gauge `g=(-1,0,0,0)` satisfies

\[
(2,0,0,-42)=\operatorname{act}_R\bigl(g,(0,1,1,-18)\bigr)
\]

exactly, with a flat action result. The nonzero integral carries are
retained even when they have the same image modulo the lower relations.
Generator orders alone do not supply this comparison.

## Interpreting results

For each case, record the actual invariant list and resolution/twist
inputs, or an explicit unresolved/error status. An unresolved result is
not zero and cannot count as a match. Compare complete invariant factors,
including free rank, rather than products of E6 layer orders. Preserve
`certified_ko=false`: these bounded comparisons do not establish
all-degree naturality or a complete ko/SPT identification.
