Resource policy update (2026-09-24): 53 still-pending jobs were resubmitted together with 6144 MiB and one CPU each; GAP workspace is capped at 5222 MiB. Already-running jobs were preserved. Details and old/new job IDs: `resource-change-2026-09-24T124555.111968_0000/policy.json`. All replacements have been admitted by Slurm; the original full-node policy below is historical.

This is the active consolidated campaign. `../latest` points here.
Dispatcher job: **390913**. All 489 cases are in `manifest.json`; the dispatcher
submits further cases as the per-user Slurm queue limit permits.

`REPORT.md` and `comparison.json` update automatically as calculations finish.
`resources.csv` contains runtime and memory measurements; raw measurements and
GAP page checkpoints are in `tasks/<case>/`. Unfinished calculations are not
classifications. The second space-group twist has no matching published table.

The initial run found that phoenix7 lacks `/usr/bin/time`. The second run used
portable kernel resource accounting, then exposed a Slurm environment issue:
full-node jobs inherited the dispatcher's 1024 MiB memory variable, giving GAP
an unintended 870 MiB workspace cap. Actual reservations were full-node memory.
The current runner queries the controller's node memory and strips inherited
Slurm resource variables when submitting jobs. Its limits were checked on a
compute node before launch and verified in actual calculation metrics:
433500 MiB GAP workspace on 510000 MiB nodes, and 652800 MiB on 768000 MiB nodes.

51 completed calculations from the earlier runs were retained, after comparing
all mathematical source hashes, reference facts, page settings, and exact task
inputs. Their original workspace limits and timings remain recorded. Their
incorrect reservation-size metadata was corrected using the submitted
`--mem=0` request and controller node memory; the original metric files are
preserved as `metrics-before-memory-correction.json`. These completed cases
used less memory than their original limit, and their mathematical results
were not changed. Unfinished calculations were restarted with the larger cap.

Reused task directories are symbolic links into the earlier run directories.
**Keep `papers-20260924` and `papers-20260924-v2`**: they retain the original
source snapshots, resource records, and calculation provenance.

The first full Table III comparison finished all 29 distinct group/twist
calculations. Of 90 dimension-by-source-row comparisons, 87 have matching phase
counts. Rows 9, 10, and 11 disagree in 1D: the computed associated-graded orders
are 64, 64, and 128 respectively, while each printed order is 16. Agreement
of phase counts does not resolve stacking extensions. Space-group comparisons
continue in the live report; disagreements are retained without modifying
either the computed output or the published table.

Validation evidence is in `verification/`: package smoke check, both bundled
examples, 69 finite-twist checks, 56 space-group twist checks, 15 Python tests,
and a full-node memory pilot. No mathematical package implementation was changed.

To refresh/read status:

```sh
python3 batch/report.py runs/papers-20260924-v3
squeue -u cliu
```
