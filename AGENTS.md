# Working on fermionAHSS

## Scope and first reads

This is the independent GAP package for five-row twisted connective real
K-theory AHSS calculations through E6. Production operations use the fixed
normalized group-bar comparison. Read `README.md`, `doc/README.md`,
`doc/mathematical-status.md`, and `doc/verification.json` before changing it.
For formulas and calibration data read `doc/conventions.md`,
`doc/secondary_operations.md`, `doc/tertiary_operations.md`, and
`doc/universal_helpers.md`. For adapter work also read `doc/backends.md`.

This folder is a Git repository. Inspect `git status` and preserve unrelated
changes. Keep scratch calculations outside the package. Do not restore
Dold–Kan transport, a model-selection option, or a local R fallback without
an explicit change in the project's intended scope.

## Layout and loading

- `PackageInfo.g`, `init.g`, `read.g`, `load.g`: package metadata and loaders.
- `gap/koahss.gd`: public declarations; `gap/group_api.gi`: finite-group entry.
- `gap/pages.gi`: exact page quotients and unresolved dependency propagation.
- `gap/hap.gi`, `cochains.gi`, `natural_bar.gi`: cochains and bar transport.
- `gap/natural_words.gi`, `data/`: chi words and exact compiled ANF data.
- `gap/natural_secondary.gi`, `natural_tertiary.gi`: production operations.
- `python/`: flattened exact Danus kernel, worker, and calibration JSON.
- `doc/`: formulas, conventions, limits, and dated verification evidence.
- `examples/`: runnable small examples; `tst/smoke.tst`: package smoke check.

Use GAP >=4.12, Polycyclic >=2.16, HAP, GAP JSON, and Python >=3.10 as
`python3` on PATH. The Python kernel uses the standard library only.
Load with `LoadPackage("fermionAHSS")`,
`ReadPackage("fermionAHSS", "load.g")`, or an absolute `Read` of `load.g`.
A bare package name is not a valid `ReadPackage` argument. Repeated loading
is harmless but does not reload edited code: always test in a fresh process.
Do not load upstream koAHSS in the same process; both packages share globals.
The local `~/.gap/pkg/fermionAHSS` symlink points to this checkout.

## Public API and mathematical contracts

`koAHSS(group,s,omega,k[,n])` constructs a resolution for a finite GAP group.
For infinite groups or control of the resolution basis, explicitly construct
an integral HAP resolution and use `koAHSSHAPSpace(R,koAHSSNaturalOperations())`
with `koAHSSpages`. No transport option is needed.

- `k` is the largest displayed physical dimension `p+q+3`, in `[-1..6]`.
  Optional `n` counts pages starting at E2, in `[1..5]`. Without `n`, return
  one E6 table. With it, return E2 through E(n+1).
- Rows are `q=-4,-3,-2,-1,0`; GAP entry is `table[q+5][p+1]`. Their lengths
  are `max(0,k-q-2)`. Retain hidden targets of outgoing maps.
- Exact invariant lists: `[]` is zero, `[0]` is Z, `[2]` is Z/2. Unresolved
  records are not zero or classifications. Invalid identities remain errors.
- Both integral rows use Z_s. Twists are degree-one/two binary cocycles in
  the actual resolution basis; scalar `0` means zero. Named cohomology
  classes are not automatically coordinate vectors. A resolution models BG.
- Preserve exact integer/mod-two/rational arithmetic and representative
  lifts/projections; ranks alone lose torsion. Differentials act on row vectors.
- The fixed lower reference is chi7_tail with epsilon `(1,0,0)` and eta
  `(1,0,1)`. Preserve the whole secondary formula and matched integral lift.
- Production T uses Danus, zero rank correction, mu_R=0, and coefficient-two
  prime-three correction in input degree three. The legacy TReference
  coefficient one must not be applied to this T.
- Current R2 is R2old minus A^cup3/4. Do not apply that correction twice or
  change the historical V2fin used in the R3 coefficient prescriptions.
- Preserve raw rational phases, both integer carries, and corrections to b
  and c during bar transport. Only allowed defining cochains are solved locally;
  universal R values come from the fixed contractors and calibration data.
- Keep every nonzero source value and selector in the bundled JSON. T covers
  only input degrees 0–3, secondary helpers 0–7. Preserve `certified_ko:false`;
  finite checks do not prove all-degree naturality or solve abutment extensions.

## Editing and verification

Follow the existing GAP and Python style. Add declarations to `gap/koahss.gd`
and new modules to `load.g` in dependency order; `read.g` shares this loader.
Public APIs use `koAHSS...`; internal helpers generally use `KOAHSS_...`.
Keep package and standalone loading functional and idempotent.

`python/phase_eval.py` locates data one parent above `python/`; the GAP
worker path is `../python/worker.py`. `r3_source.json` records source hashes,
and `high_calibration.json` records its hash. Never update hashes merely to
silence failed validation. A formula change requires mathematical review and
recalibration; a proven path-only migration requires an explicit record that
the numerical payloads are unchanged. Preserve the MIT license attribution.

Run bounded checks from the package root:

```sh
gap -q --quitonbreak examples/c2.g
gap -q --quitonbreak examples/twisted_c2.g
gap -q --quitonbreak -c 'LoadPackage("fermionAHSS"); Assert(0,TestPackage("fermionAHSS")); QUIT;'
```

After mathematical changes, check the affected identities, torsion, lifts,
invalid inputs, and direct operation degrees. A small page run may never
invoke T1–T3. The initialization audit in `doc/verification.json` also used
upstream fixtures temporarily; they are not all included here. Do not claim
that running the bundled smoke check reruns that audit. Record fresh commands,
results, exclusions, and dependency skips. No full group catalogue is needed
for packaging changes. Run code review before committing code changes.

Report what changed, exact validation, saved evidence, and remaining scope.
Keep dated verification evidence historical rather than silently overwriting
it as proof of later calculations. Do not push or publish unless requested.
