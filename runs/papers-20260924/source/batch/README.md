The paper campaign runs the actual package in this directory, using the full
infinite affine space groups and exact rational Clifford lifts for
`w2+w1^2`. Group/twist constructors were adapted from the neighbouring
`fermionicSPT/batch` scripts; no legacy operation backend is loaded.

From the package root:

```sh
python3 batch/campaign.py prepare runs/papers-20260924
python3 batch/campaign.py submit runs/papers-20260924
python3 batch/report.py runs/papers-20260924
```

Defaults select all 230 space groups, separately for `(w1,0)` and
`(w1,w2+w1^2)`, and all fixed rows of Table III of arXiv:1811.00536.
Its infinite families are sampled at parameters 1–4 (`--family-max` changes
that explicit range). There are 489 isolated calculations: 29 finite cases
and 460 labelled space-group cases. Coincident space-group twist classes
remain separately labelled. Quaternion symmetry uses the bosonic quotient
`C2 x C2` with extension `x^2+xy+y^2`.

The paper's spatial dimension is `d=p+q+2`. The package display cutoff is
`k=p+q+3`, so the default **k=4**, E2–E6, covers paper dimensions 1–3.
In 3D the layer locations are p+ip `(1,0)`, Kitaev `(2,-1)`, complex
fermion `(3,-2)`, and bosonic `(5,-4)`. The finite-group phase counts omit
intrinsic p=0 Kitaev/p+ip layers in 1D/2D, as Table III does.

Each run freezes GAP/Python code, calibration data, scripts, and paper facts
in `source/`, with SHA256 hashes in `manifest.json`. The installed external
GAP packages remain dependencies; each space-group event records the exact
SpaceGroupCohomology module hash. A fresh run directory is required after
changes or for retries. Existing results are never silently overwritten.

The dispatcher maintains at most 150 total submitted jobs for this user,
including unrelated jobs and itself. It submits the next waiting task as
capacity becomes available. Slurm enforces the observed 100-running-job and
200-CPU account limits. Calculations request one CPU, `--mem=0` (all
allocatable memory of their normalx node), and up to 14 days. On the inspected
cluster this is 510000 or 768000 MiB per node; GAP's workspace limit uses 85%
of it, leaving the rest for Python and other subprocesses. Memory reservations
can limit concurrency until other jobs release their allocations. The
dispatcher runs on normal with 1 GiB and renews itself before its time limit.

Every completed page is flushed to `tasks/<case>/events.jsonl`. The task
directory also holds exact inputs, GAP log, GNU time output, and live/final
metrics. Measurements are:

- Wall seconds for setup, resolution, twist construction, and all page calls,
  excluding queue wait. Incremental page checkpointing may repeat prior-page
  work; this overhead is included.
- GNU time maximum RSS (KiB): largest process high-water mark, including
  inherited child resource usage, not the sum of simultaneous processes.
- Sampled process-tree peak RSS (KiB), every two seconds: concurrent GAP and
  descendants; short spikes can be missed and shared pages can be double-counted.
- Job-cgroup peak bytes when the cluster exposes a job-specific counter.
  Login-session counters are deliberately excluded.

The dispatcher refreshes `REPORT.md`, `summary.json`, `results.json`,
`comparison.json`, and `resources.csv` every minute. Each completed worker
also refreshes them. A missing terminal measurement after a Slurm kill stays
explicit; it is not a zero runtime or a successful calculation. The cluster's
Slurm accounting database was unavailable during setup, so measurements do
not depend on sacct.

The package does **not** solve stacking extensions or produce certified total
SPT groups. All records retain `certified_ko=false`. Table III comparisons
test phase counts only, never equality of a direct sum of layers to a stacked
classification. Space-group comparisons test individual E6 layers only for
the paper's effective electronic `(w1,0)` twist. The second twist is not
tabulated there. Blank SG210/219/228 cells remain unreported. Printed totals
and products of printed layers are kept separate, including inconsistencies.

Reference JSON was extracted from arXiv TeX, retaining raw table cells and
source hashes. Rebuild it from downloaded/extracted sources with:

```sh
python3 batch/extract_references.py main_classification.tex Gf.tex references
python3 -m unittest discover -s batch -p 'test_*.py' -v
gap -q --quitonbreak batch/test_paper_twists.g
gap -q --quitonbreak batch/test_spacegroup_twists.g
```

For a bounded local check use `prepare /tmp/fahss-check --groups table3-C2
--k 2 --pages 2`, then `work /tmp/fahss-check table3-C2-row1`.
