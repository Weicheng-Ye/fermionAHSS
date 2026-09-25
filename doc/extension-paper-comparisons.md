# Full-group comparison fixtures

Reference extraction and sample calculations checked on 2026-09-25.
**All 20 labeled comparisons match, with no unresolved cases or
disagreements.** The initial low-degree run recorded 17 matches and 3
unresolved comparisons; its evidence is preserved below.
There are four distinct C2 calculations; the crystalline comparisons reuse
their degree-four results. The expected groups below are literature fixtures.
Machine-readable inputs and expectations are in
[extension-paper-samples.json](../data/extension-paper-samples.json).

The papers' spatial dimension is `d=p+q+2`, while the package degree is
`k=p+q+3=d+1`; see [the campaign conventions](../batch/README.md).
Thus paper dimensions 0–3 correspond to package degrees 1–4, and a
`koFull(...,4)` result uses entries `full.invariants[d+3]` for these cases.
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
The fixture preserves both statements without selecting one; no C4 run
is recorded. This S4 is abstractly cyclic C4, not the symmetric group.

## Historical low-degree results on 2026-09-25

The initial version of `examples/extension_papers.g` constructed each full
result independently of the fixture, then compared the actual invariant
lists. Its machine-readable evidence is
[extension-paper-results-20260925.json](../data/extension-paper-results-20260925.json).

| Input `(s,omega)` | Computed package degrees 1, 2, 3, 4 | Comparison |
| --- | --- | --- |
| `(0,0)` | `[2,2]`, `[2,2]`, unresolved, `[]` | 3 matches; degree 3 unresolved |
| `(0,[1])` | `[4]`, `[]`, `[0]`, `[]` | 4 matches |
| `([1],0)` | `[2]`, `[8]`, `[]`, `[]` | 4 matches |
| `([1],[1])` | `[]`, `[2]`, `[2]`, unresolved | 3 matches; degree 4 unresolved |

Thus 14 of the 16 Wang–Gu entries match. The three zero-group crystalline
entries match; spinless reflection remains unresolved, with expected `[16]`.
The two unresolved internal calculations need the general upper-lift and
lower-class reduction machinery beyond the then-supported stacking
adapter. The unresolved reflection comparison repeats the second of these
calculations; it is not an independent third obstruction.

The nonsplit `[4]` result comes from relations `2d=0, 2c=d`.
The `[8]` result comes from `2d=0, 2c=d, 2b=c`, with relation matrix
`[[2,0,0],[-1,2,0],[0,-1,2]]`. These are measured low-degree stacking
vectors followed by Smith reduction, not classifications supplied to the
solver. The expected `[0,8]` and `[16]` were not substituted for missing
computations.

## Complete-state follow-up on 2026-09-25

Run `gap -q --quitonbreak examples/extension_papers.g` from the package root.
The current run completed all four distinct C2 calculations. All 16
Wang–Gu entries and all four crystalline entries match, including the
previously unresolved Z ⊕ Z/8 and Z/16 groups. The crystalline entries
reuse the corresponding internal calculation, so these are 20 labeled
comparisons, not 20 independent runs.

The [complete-state results](../data/extension-paper-results-20260925-full.json)
record actual invariant lists, relation matrices, fixed flat lifts,
measured power vectors, exact gauge witnesses, finite multiplication
tables and source hashes. The solver never consults the expected groups;
the sample runner compares them only after `koFull` returns.

| Input `(s,omega)` | Computed package degrees 1, 2, 3, 4 | Comparison |
| --- | --- | --- |
| `(0,0)` | `[2,2]`, `[2,2]`, `[0,8]`, `[]` | 4 matches |
| `(0,[1])` | `[4]`, `[]`, `[0]`, `[]` | 4 matches |
| `([1],0)` | `[2]`, `[8]`, `[]`, `[]` | 4 matches |
| `([1],[1])` | `[]`, `[2]`, `[2]`, `[16]` | 4 matches |

The higher-degree implementation first solves the full defining equations
for every marked A, B, C and D generator, storing immutable flat tuples.
Every power and lower-coordinate calculation reuses those tuples. It
retains the exact boundary comparison with the canonical lower product
and performs a finite multiplication-table audit before completing a
higher-degree result. See [the implementation contract](extensions.md#complete-flat-representatives-and-boundary-comparisons).

For unitary C2 in package degree 3, the measured relations in the marked
D/C/B basis are

\[
2D=0,\qquad 2C=D,\qquad 2B=C+D.
\]

Their Smith form gives Z/8. The A generator is free, with its own stored
full flat representative, giving Z ⊕ Z/8 in the abelian abutment category.

For C2 with `s=omega=[1]` in package degree 4, the measured relations are

\[
2D=0,\qquad 2C=D,\qquad 2B=C+7D,\qquad 2A=B+C-17D.
\]

All coordinates refer to the same stored D/C/B/A columns. The joint
presentation and its nonunit Smith factor are

\[
R=\begin{pmatrix}
2&0&0&0\\
-1&2&0&0\\
-7&-1&2&0\\
17&-1&-1&2
\end{pmatrix},\qquad \operatorname{coker}(R)\cong\mathbb Z/16.
\]

In the complete normalized C2 bar basis, the selected A lift is
`(A,B,C,D)=(1,0,0,0)`. Its actual double is `(2,0,0,-42)`.
The canonical lower product for the measured coordinate vector `[-17,1,1]`
is `(0,1,1,-18)`. The stored preceding-degree gauge is `(-1,0,0,0)`,
with boundary `(2,1,1,-29)`. The implementation verifies the exact ordered
cochain equality

\[
(2,0,0,-42)=(2,1,1,-29)\mathbin{\times}(0,1,1,-18),
\]

and verifies that the boundary is flat. The nonzero integral carries
are retained even when they have the same image modulo the lower
relations. Generator orders alone do not supply this comparison.

The unitary degree-three torsion part passed all 64 ordered products of
its eight marked normal forms. The signed degree-four result passed all
256 ordered products of its 16 forms. Each comparison used literal
equality or a solved exact boundary witness; the resulting tables passed
identity, inverse, commutativity and associativity checks.

Fresh verification of this implementation passed:

- The full GAP package suite in GAP 4.15.1, including API/display
  compatibility, correlated extension vectors, forced B/C defining-choice
  changes, integral D primitives, fixed-lift reuse, actual C2 A-layer
  comparisons and finite quotient audits. A synthetic negative audit
  rejects an unsupported Z/2 relation instead of reporting a group.
- `python3 -m unittest discover -s python -p 'test_extension*.py' -v`:
  eight tests covering exact D carries, cache/input validation and the
  degree-six finite section. The nonzero off-shell A test checks
  `d(x × y)=d(x) × d(y)` and `d²(x)=0`, including the degree-seven
  endpoint; it is a direct operation test, not a nontrivial degree-six
  full-group classification.
- `koFull(TrivialGroup(),0,0,6)` returns
  `[[0],[],[2],[2],[0],[],[],[]]` across degrees -1 through 6.
- Code review, bundled-source checksum validation, documentation links
  and `git diff --check`.

## Interpreting results

For each case, record the actual invariant list and resolution/twist
inputs, or an explicit unresolved/error status. An unresolved result is
not zero and cannot count as a match. Compare complete invariant factors,
including free rank, rather than products of E6 layer orders. Preserve
`certified_ko=false`: these bounded comparisons do not establish
all-degree naturality or a complete ko/SPT identification.

## Historical implementation verification on 2026-09-25

The following checks passed in fresh processes:

- `gap -q --quitonbreak examples/c2.g` and `examples/twisted_c2.g`.
- `LoadPackage("fermionAHSS"); Assert(0,TestPackage("fermionAHSS"));`
  in GAP 4.15.1, including API/display, abstract extension, and production
  stacking tests. These cover correlated images, a Z/16 abstract chain,
  free summands, basis transformations, context reuse, and signed C2 × C2
  with two independent C-layer images and computed invariants `[4,8]`.
- `python3 -m unittest discover -s python -p test_stacking_low_phase.py -v`:
  four tests covering the payload checksum, every local signed integral
  carry, and explicit universal associativity/exchange boundary witnesses.
- `examples/extension_papers.g` through the standalone loader from a
  different working directory; the 20 recorded rows were reproduced.
- `koFull(TrivialGroup(),0,0,6)` covers degrees -1 through 6 and returns
  `[[0],[],[2],[2],[0],[],[],[]]`.
- `git diff --check` and code review of the implementation.

The results JSON records hashes of the runtime files, local phase data,
fixture and sample runner used for these comparisons. Existing page
calibration payloads and their hashes were not changed.
