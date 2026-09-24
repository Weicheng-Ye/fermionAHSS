#!/usr/bin/env python3
"""Resubmit only pending campaign tasks with a smaller per-job memory budget."""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

from campaign import now, queue, save, submission_environment


def resize(run, memory_mib):
    run = run.resolve()
    if memory_mib <= 0:
        raise ValueError("an explicit positive memory budget is required")
    manifest = json.loads((run/"manifest.json").read_text())
    active = queue()
    targets = []
    for task in manifest["tasks"]:
        directory = run/"tasks"/task["id"]
        submission = directory/"submission.json"
        if not submission.exists() or (directory/"metrics.json").exists():
            continue
        old = json.loads(submission.read_text())
        if active.get(old["job_id"]) == "PENDING":
            targets.append((task, directory, old))
    if not targets:
        print("No pending calculation jobs to resize.")
        return
    stamp = now().replace(":", "").replace("+", "_")
    control = run/("resource-change-"+stamp)
    control.mkdir()
    # Freeze the resource controller separately; the original mathematical
    # source snapshot and all already-running processes remain untouched.
    for name in ("campaign.py", "catalogue.py", "resize_pending.py"):
        shutil.copy2(Path(__file__).parent/name, control/name)
    policy = dict(created_utc=now(), memory_mib=memory_mib, cpus_per_task=1,
                  source_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in control.glob("*.py")},
                  mathematical_source=str(run/"source"),
                  old_pending_jobs=[old["job_id"] for _, _, old in targets],
                  replacements=[], skipped=[])
    save(control/"policy.json", policy)
    # State filtering protects jobs that started after the queue snapshot.
    subprocess.run(["scancel", "--quiet", "--state=PENDING",
                    ",".join(old["job_id"] for _, _, old in targets)], check=True)
    active = queue()
    for task, directory, old in targets:
        if old["job_id"] in active or (directory/"metrics-live.json").exists() or (directory/"events.jsonl").exists():
            policy["skipped"].append(dict(task=task["id"], job=old["job_id"], reason="started or still active during resize"))
            save(control/"policy.json", policy)
            continue
        save(directory/("submission-before-"+stamp+".json"), old)
        save(directory/"resource-request.json", dict(memory_mib=memory_mib, changed_utc=now(),
                                                     policy=str(control/"policy.json")))
        script = control/(task["id"]+".sh")
        script.write_text("#!/bin/bash\nset -eu\nexec "+shlex.join([
            sys.executable, str(control/"campaign.py"), "work", str(run), task["id"]])+"\n")
        command = ["sbatch", "--parsable", "--partition="+manifest["partition"],
                   "--nodes=1", "--ntasks=1", "--cpus-per-task=1", f"--mem={memory_mib}M",
                   "--time="+manifest["time_limit"], "--signal=B:TERM@120",
                   "--job-name=fahss-"+task["id"],
                   "--output="+str(run/"slurm"/(task["id"]+"-%j.log")), str(script)]
        result = subprocess.run(command, capture_output=True, text=True, check=True,
                                env=submission_environment())
        job_id = result.stdout.strip().split(";")[0]
        save(directory/"submission.json", dict(job_id=job_id, submitted_utc=now(), command=command,
                                               previous_job_id=old["job_id"], resource_policy=str(control/"policy.json")))
        policy["replacements"].append(dict(task=task["id"], old_job_id=old["job_id"], new_job_id=job_id))
        save(control/"policy.json", policy)
        print(task["id"], old["job_id"], "->", job_id, flush=True)
    save(run/"queue.json", dict(updated_utc=now(), jobs=queue()))
    print(f"Resubmitted {len(policy['replacements'])} tasks at {memory_mib} MiB; "
          f"left {len(policy['skipped'])} raced tasks untouched. Policy: {control}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--memory-mib", type=int, default=6144)
    args = parser.parse_args()
    resize(args.run, args.memory_mib)
