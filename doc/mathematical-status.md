# Mathematical status and provenance

fermionAHSS is the fixed group-bar edition of the koAHSS research
implementation. Its page calculation covers rows `q=-4,-3,-2,-1,0`
through E6 and physical cutoff `p+q+3 <= k <= 6`. E6 is terminal within
this strip; it is not automatically the full ko E-infinity page. Neither
other ko rows nor abutment/stacking extensions are computed. Saved output
therefore retains `certified_ko: false`.

## Fixed secondary family

The runtime family is `chi7_tail` in input degrees zero through seven,
with epsilon `(1,0,0)` and eta `(1,0,1)`. Its word identities, right-product
suspension, and Thom calibration are the mathematical inputs to the
implementation. The exact common integral lift satisfies
`rho M'=psi'` and `d_s M'=2 tau'`. The square-zero assertion `J[Tau]=0`
means an integral coboundary, not a literally zero cochain.

The source-workspace provenance is
`note/extra/chi_suspension_degree7/`, including `tail_words.json`, and
`note/extra/secondary_operations_psi_note.tex`, appendix equations
`Th-squarezero` and theorem `global-squarezero`. These proof files are not
bundled. The packaged [ANF data](../data/chi-calibrated-degree7-anf.g) and
[helper reference](universal_helpers.md) fully specify runtime evaluation.
Historical `head`, `tail6`, and constrained degree-nine references have
different normalizations and cannot be interchanged with `chi7_tail`.
Local residual and injective-J diagnostics do not select production maps.

## Fixed tertiary family

The Danus universal helper and final T are implemented only in input degrees
zero through three. The general signed R0, R1, R2, and R3 family is used;
the compact even-input family is not silently spliced into it. R0 retains
full and quarter-input carries. Current R2 is
`R2sharp=R2old-A^cup3/4`; the R3 coefficient prescriptions still use
historical `V2fin` in `C0` and xi. Applying that cubic correction twice or
absorbing it into the binary L6 helper changes the operation.

The source supplements establish low-degree first-b descent, additivity,
and the specified suspension comparisons. These are separate claims from
residual solvability, exact boundary, and absolute normalization. The
rank-six `BPSO(6)=BPU(4)` comparison fixes the `q(omega)A` correction to
zero relative to Danus. The signed Euler calculation gives `mu_R=0`.
Thus final T adds `2 beta_3,s P^1_s rho_3,s`, which can contribute only
in input degree three here. The old coefficient-one correction belongs
only to legacy `TReference`.

The nonzero finite data are essential: low selectors `(0,1,0)`, V2 source
values `(0,3/4,0,3/4,1/4)`, R3 selectors `(1,0,1,1,1)`, xi `3/4`,
and suspension periods `(3/4,1/4,0,1/2,1/2)`. The R3 source calculation
has a 46-by-94 matrix and all 63 cycle-obstruction checks pass; source
denominators reach sixteen. The final rational phase may also have a
factor of three. See [tertiary formulas](tertiary_operations.md) and
[universal helpers](universal_helpers.md).

No all-degree universal R extension is asserted. The older degree-six
quarter/eighth-denominator obstructions concern a different range and do
not contradict the implemented low-degree family. The older polynomial
L6 search and the 14-word-helper/first-lift H3 problem are not solved by
renaming the current finite contractors.

GAP solves only the permitted defining cochains b and c. R is evaluated
from the fixed universal operators and calibrated source values, never
chosen by a local residual solve. A solution using the last nullhomotopy,
such as R=-U(c), would not define the required universal R(A,b).

## Evidence and packaging

Original research-workspace provenance (not bundled runtime dependencies):

- `notes/extra/tertiary_R_Danus/`: general R formulas and suspension
  arguments, including `r3.md` and `r2_cubic_normalization.md`.
- `notes/extra/tertiary_T_degree3/`: low-degree T descent, additivity,
  rank-six calibration, and normalization arguments.
- `notes/extra/tertiary_T_degree3/mu_verification/`: exact Euler calibration
  with `Xi(D)=V2(D)=13/4`, `kappa6 V2(Ss)(Tor)=3/2`, and
  `Uraw(Tor)=-1/2`; this is distinct from the production R3 computation.
- `notes/extra/tertiary_T_degree3/gap_verification_20260923/`: executed
  universal source/calibration checks, 56 Python tests, 20 invalid-input
  cases, and 20 bounded E2–E6 integrations. Those page integrations
  required only degree-zero T. Separate direct calls covered input degrees
  one through three, including signed C2 degree three with phase `-197/24`
  and T=0, and a C3 times C3 nonzero prime-three detector.
- `note/extra/tertiary_R_degree3/`: earlier polynomial data and older phase
  suspension fixtures, including 33/240 shuffle and 3/48 slant failures.
  Their failures remain historical evidence for that earlier normalization;
  new symbolic arguments do not mean those tests were rerun or repaired.

The flattened `python/` layout changes the package-data lookup in
`phase_eval.py`. Its packaged source hash and the dependent R3 source hash
are updated to match that path-only packaging change. Numeric source
coefficients, calibration values, and mathematical formulas are unchanged;
this is not a new calibration.

Finite executable checks support the adapter and arithmetic. They do not
replace the universal arguments or establish an all-degree theorem. Fresh
package checks are recorded separately from the historical evidence above.
