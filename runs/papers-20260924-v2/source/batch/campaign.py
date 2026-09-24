#!/usr/bin/env python3
"""Prepare, submit, monitor, and measure isolated paper calculations on Slurm."""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
import resource
from pathlib import Path
import shlex
import shutil
import signal
import socket
import subprocess
import sys
import time

from catalogue import build_catalogue

ROOT = Path(__file__).resolve().parents[1]


def now():
    return datetime.now(timezone.utc).isoformat()


def save(path, value):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def events(path):
    if not path.exists():
        return []
    lines = path.read_bytes().splitlines(keepends=True)
    result = []
    for i, line in enumerate(lines):
        try:
            result.append(json.loads(line))
        except (ValueError, UnicodeDecodeError):
            if i == len(lines)-1 and not line.endswith(b"\n"):
                break
            raise
    return result


def tasks_for(catalogue):
    tasks = []
    for spec in catalogue:
        if spec["kind"] == "spacegroup":
            for twist in spec["spacegroup_twists"]:
                family = twist["spacegroup_twist"]
                tasks.append(dict(id=spec["id"] + "-" + family,
                                  spec=spec, spacegroup_twist=family))
        else:
            # Deduplicate only identical named recipes, retaining source rows.
            recipes = defaultdict(list)
            for row in spec["source_rows"]:
                key = (row["paper_s"], row["paper_omega_kind"],
                       row.get("paper_omega_cyclic_order"))
                recipes[key].append(row)
            for rows in recipes.values():
                reduced = dict(spec, source_rows=rows)
                tasks.append(dict(id=spec["id"] + "-row" + str(rows[0]["row"]), spec=reduced))
    # Interleave finite and space-group examples so early results cover both
    # papers even when some of the finite product groups take a long time.
    finite = [t for t in tasks if t["spec"]["kind"] == "abelian"]
    spatial = [t for t in tasks if t["spec"]["kind"] == "spacegroup"]
    result = []
    for i in range(max(len(finite), len(spatial))):
        if i < len(finite):
            result.append(finite[i])
        if i < len(spatial):
            result.append(spatial[i])
    return result


def prepare(args):
    run = args.run.resolve()
    run.mkdir(parents=True, exist_ok=False)
    source = run / "source"
    source.mkdir()
    for name in ("batch", "gap", "python", "data", "references"):
        shutil.copytree(ROOT/name, source/name, ignore=shutil.ignore_patterns("__pycache__"))
    for name in ("load.g", "PackageInfo.g", "init.g", "read.g", "LICENSE"):
        shutil.copy2(ROOT/name, source/name)
    hashes = {str(p.relative_to(source)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(source.rglob("*")) if p.is_file()}
    reused = {}
    prior = None
    if args.reuse_run:
        prior = args.reuse_run.resolve()
        previous = json.loads((prior/"manifest.json").read_text())
        protected = [n for n in hashes if n.startswith(("gap/", "python/", "data/", "references/"))
                     or n.endswith(".g") or n == "batch/catalogue.py"]
        if any(previous["source_sha256"].get(n) != hashes[n] for n in protected):
            raise ValueError("cannot reuse results after changing mathematical code or reference inputs")
        if (previous["k"], previous["pages"]) != (args.k, args.pages):
            raise ValueError("cannot reuse results from a different page window")
        active = queue()
    catalogue = build_catalogue(range(1, args.family_max+1), range(1, args.family_max+1))
    tasks = tasks_for(catalogue)
    if args.groups:
        requested = set(args.groups.split(","))
        tasks = [t for t in tasks if t["spec"]["id"] in requested]
        if requested != {t["spec"]["id"] for t in tasks}:
            raise ValueError("unknown group selection")
    (run/"tasks").mkdir()
    (run/"slurm").mkdir()
    for task in tasks:
        directory = run/"tasks"/task["id"]
        if prior:
            old = prior/"tasks"/task["id"]
            metric = old/"metrics.json"
            submission = old/"submission.json"
            old_id = json.loads(submission.read_text())["job_id"] if submission.exists() else None
            status = json.loads(metric.read_text())["state"] if metric.exists() else None
            if status == "computed" or active.get(old_id) in ("RUNNING", "COMPLETING"):
                if json.loads((old/"job.json").read_text())["spec"] != task["spec"]:
                    raise ValueError("reused task inputs differ")
                directory.symlink_to(old, target_is_directory=True)
                reused[task["id"]] = str(old)
                continue
        directory.mkdir()
        job = dict(root=str(source), spec=task["spec"], results=str(directory/"events.jsonl"),
                   k=args.k, pages=args.pages, twists="paper", operations="chosen", seed=1,
                   max_twists=0, completed=[], previous_basis=None, previous_boundary=None)
        if "spacegroup_twist" in task:
            job["spacegroup_twist"] = task["spacegroup_twist"]
        save(directory/"job.json", job)
    manifest = dict(created_utc=now(), package="fermionAHSS", certified_ko=False,
                    k=args.k, pages=args.pages, paper_spatial_dimension_offset=1,
                    family_parameters=list(range(1, args.family_max+1)),
                    partition=args.partition, memory="all allocatable node memory (--mem=0)",
                    time_limit=args.time_limit, tasks=tasks, source_sha256=hashes,
                    gap=shutil.which("gap"), max_submitted=args.max_submitted,
                    reused_tasks=reused, prior_run=str(prior) if prior else None,
                    scope="Five-row E6 associated graded only; no stacking extensions")
    save(run/"manifest.json", manifest)
    for command, name in [(["scontrol", "show", "assoc_mgr", "users="+os.environ["USER"]], "limits"),
                          (["sinfo", "-N", "-l"], "nodes"),
                          (["scontrol", "show", "partition", args.partition], "partition")]:
        result = subprocess.run(command, capture_output=True, text=True)
        (run/(name+".txt")).write_text(result.stdout+result.stderr)
    print(f"Prepared {len(tasks)} isolated calculations in {run}")


def allocated_mib():
    # mem=0 may be exported literally as zero. Ask the controller for the
    # assigned node's configured memory rather than using host MemTotal.
    value = int(os.environ.get("SLURM_MEM_PER_NODE", "0"))
    if value:
        return value
    node = os.environ.get("SLURMD_NODENAME")
    if node:
        import re
        result = subprocess.run(["scontrol", "show", "node", node, "-o"],
                                capture_output=True, text=True, check=True)
        return int(re.search(r"RealMemory=(\d+)", result.stdout)[1])
    return 4096


def tree_rss(pid):
    """Sample simultaneous resident memory of GAP and all live descendants."""
    pending, seen, total = [pid], set(), 0
    while pending:
        current = pending.pop()
        if current in seen:
            continue
        seen.add(current)
        try:
            status = Path(f"/proc/{current}/status").read_text()
            total += next((int(s.split()[1]) for s in status.splitlines()
                           if s.startswith("VmRSS:")), 0)
            pending.extend(map(int, Path(f"/proc/{current}/task/{current}/children").read_text().split()))
        except (OSError, ProcessLookupError):
            pass
    return total


def cgroup_peak():
    if not os.environ.get("SLURM_JOB_ID"):
        return None  # A login-session cgroup includes unrelated processes.
    try:
        for line in Path("/proc/self/cgroup").read_text().splitlines():
            _, controllers, group = line.split(":", 2)
            if controllers == "":
                path = Path("/sys/fs/cgroup")/group.lstrip("/")/"memory.peak"
            elif "memory" in controllers.split(","):
                path = Path("/sys/fs/cgroup/memory")/group.lstrip("/")/"memory.max_usage_in_bytes"
            else:
                continue
            if path.exists() and ("job_" in str(path) or "job-" in str(path)):
                return dict(bytes=int(path.read_text()), path=str(path))
    except (OSError, ValueError):
        pass
    return None


def work(args):
    run = args.run.resolve()
    manifest = json.loads((run/"manifest.json").read_text())
    directory = run/"tasks"/args.task
    lock = (directory/".lock").open("w")
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    if (directory/"metrics.json").exists():
        raise RuntimeError("task already has a terminal result; prepare a fresh run")
    budget = allocated_mib()
    gap_limit = int(budget * .85)  # Leave room for Python, polymake, and OS.
    command = [manifest["gap"], "-q", "--quitonbreak", "-o", f"{gap_limit}m",
               "-K", f"{gap_limit}m", str(run/"source/batch/koahss_worker.g")]
    metrics = dict(started_utc=now(), host=socket.gethostname(), command=command,
                   slurm_job_id=os.environ.get("SLURM_JOB_ID"),
                   allocated_memory_mib=budget, gap_workspace_limit_mib=gap_limit,
                   sampled_tree_peak_rss_kib=0, sampling_interval_seconds=2,
                   state="running")
    save(directory/"metrics-live.json", metrics)
    start = time.monotonic()
    stopping = []
    def stop(signum, frame):
        stopping.append(signum)
    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)
    env = dict(os.environ, KOAHSS_JOB_FILE=str(directory/"job.json"),
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    with (directory/"gap.log").open("w") as log:
        usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
        process = subprocess.Popen(command, cwd=run/"source", env=env,
                                   stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        terminated = False
        while process.poll() is None:
            metrics["sampled_tree_peak_rss_kib"] = max(metrics["sampled_tree_peak_rss_kib"], tree_rss(process.pid))
            metrics["wall_seconds"] = round(time.monotonic()-start, 3)
            metrics["updated_utc"] = now()
            save(directory/"metrics-live.json", metrics)
            if stopping and not terminated:
                os.killpg(process.pid, signal.SIGTERM)
                terminated = True
            if terminated:
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
            try:
                process.wait(timeout=2)
            except subprocess.TimeoutExpired:
                pass
    metrics.update(ended_utc=now(), wall_seconds=round(time.monotonic()-start, 3),
                   returncode=process.returncode, cgroup_peak=cgroup_peak())
    # One measured GAP child per fresh worker. Kernel accounting propagates
    # descendants that GAP waits for and requires no node-local GNU time binary.
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    metrics["kernel_max_rss_kib"] = usage.ru_maxrss
    metrics["user_cpu_seconds"] = usage.ru_utime-usage_before.ru_utime
    metrics["system_cpu_seconds"] = usage.ru_stime-usage_before.ru_stime
    metrics["measurement_method"] = "getrusage(RUSAGE_CHILDREN), Linux ru_maxrss in KiB"
    history = events(directory/"events.jsonl")
    cases = [e for e in history if e.get("type") == "case"]
    if stopping:
        state = "interrupted"
    elif process.returncode or not history or history[-1].get("type") != "group_complete":
        state = "failed"
    elif not cases or any(c["status"] == "error" for c in cases):
        state = "calculation_error"
    elif any(c["status"] == "partial" for c in cases):
        state = "partial"
    else:
        state = "computed"
    metrics["state"] = state
    save(directory/"metrics.json", metrics)
    subprocess.run([sys.executable, str(run/"source/batch/report.py"), str(run)], check=True)
    return int(state not in ("computed", "partial"))


def queue():
    result = subprocess.run(["squeue", "-h", "-r", "-u", os.environ["USER"], "-o", "%i|%T"],
                            capture_output=True, text=True, check=True)
    return dict(line.split("|", 1) for line in result.stdout.splitlines() if "|" in line)


def submit_task(run, manifest, task):
    script = run/"slurm"/(task["id"]+".sh")
    script.write_text("#!/bin/bash\nset -eu\nexec " + shlex.join([
        sys.executable, str(run/"source/batch/campaign.py"), "work", str(run), task["id"]]) + "\n")
    command = ["sbatch", "--parsable", "--partition="+manifest["partition"],
               "--nodes=1", "--ntasks=1", "--cpus-per-task=1", "--mem=0",
               "--time="+manifest["time_limit"], "--signal=B:TERM@120",
               "--job-name=fahss-"+task["id"],
               "--output="+str(run/"slurm"/(task["id"]+"-%j.log")), str(script)]
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(result.stderr.strip())
    job_id = result.stdout.strip().split(";")[0]
    save(run/"tasks"/task["id"]/"submission.json",
         dict(job_id=job_id, submitted_utc=now(), command=command))
    return job_id


def dispatch(args):
    run = args.run.resolve()
    manifest = json.loads((run/"manifest.json").read_text())
    lock = (run/".dispatch.lock").open("w")
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    dispatch_start = time.monotonic()
    while True:
        try:
            active = queue()
            capacity = max(0, manifest["max_submitted"]-len(active))
            pending = [t for t in manifest["tasks"]
                       if not (run/"tasks"/t["id"]/"submission.json").exists()]
            for task in pending[:capacity]:
                job_id = submit_task(run, manifest, task)
                print(now(), task["id"], job_id, flush=True)
                active[job_id] = "PENDING"
            save(run/"queue.json", dict(updated_utc=now(), jobs=active))
            report = run/"source/batch/report.py"
            if report.exists():
                subprocess.run([sys.executable, str(report), str(run)], check=True)
            submitted = [json.loads(p.read_text())["job_id"] for p in (run/"tasks").glob("*/submission.json")]
            if len(submitted) == len(manifest["tasks"]) and not any(j in active for j in submitted):
                save(run/"dispatcher-complete.json", dict(ended_utc=now(), submitted=len(submitted)))
                return 0
            # A calculation can use 14 days; renew the dispatcher before its
            # own wall limit so later queued tasks still get submitted/reported.
            if time.monotonic()-dispatch_start > 13*86400 and len(active) < manifest["max_submitted"]:
                submit(args)
                return 0
        except (subprocess.SubprocessError, OSError, RuntimeError) as error:
            print(now(), "dispatcher retry:", str(error), flush=True)
        if args.once:
            return 0
        time.sleep(60)


def submit(args):
    run = args.run.resolve()
    script = run/"slurm/dispatch.sh"
    script.write_text("#!/bin/bash\nset -eu\nexec " + shlex.join([
        sys.executable, str(run/"source/batch/campaign.py"), "dispatch", str(run)]) + "\n")
    command = ["sbatch", "--parsable", "-p", "normal", "-c", "1", "--mem=1G",
               "--time=14-00:00:00", "--job-name=fahss-dispatch",
               "--output="+str(run/"slurm/dispatch-%j.log"), str(script)]
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    save(run/"dispatcher.json", dict(job_id=result.stdout.strip().split(";")[0],
                                    command=command, submitted_utc=now()))
    print(result.stdout.strip())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prep = commands.add_parser("prepare")
    prep.add_argument("run", type=Path)
    prep.add_argument("--k", type=int, default=4, choices=range(-1, 7))
    prep.add_argument("--pages", type=int, default=5, choices=range(1, 6))
    prep.add_argument("--family-max", type=int, default=4)
    prep.add_argument("--groups")
    prep.add_argument("--partition", default="normalx")
    prep.add_argument("--time-limit", default="14-00:00:00")
    prep.add_argument("--max-submitted", type=int, default=150)
    prep.add_argument("--reuse-run", type=Path, help="reuse completed/running tasks with identical mathematical sources")
    for name in ("submit", "dispatch", "work"):
        sub = commands.add_parser(name)
        sub.add_argument("run", type=Path)
        if name == "work":
            sub.add_argument("task")
        if name == "dispatch":
            sub.add_argument("--once", action="store_true")
    args = parser.parse_args()
    return globals()[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
