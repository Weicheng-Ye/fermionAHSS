#!/usr/bin/env python3
"""Export current results, resource measurements, and source-grounded comparisons."""
import argparse
from collections import Counter
import csv
import fcntl
import json
from math import prod
from pathlib import Path

from campaign import events, now, save


def primary(values):
    answer = []
    for value in values:
        if value == 0:
            answer.append(0)
            continue
        divisor = 2
        while divisor*divisor <= value:
            power = 1
            while value % divisor == 0:
                value //= divisor
                power *= divisor
            if power > 1:
                answer.append(power)
            divisor += 1
        if value > 1:
            answer.append(value)
    return sorted(answer)


def cell(table, p, q):
    if table is None or p < 0 or p >= len(table[q+4]):
        return None
    result = table[q+4][p]
    return result if isinstance(result, list) else None


def finite_expected(row, reference):
    number, params = row["row"], row["parameters"]
    if number == 2:
        return {"1": [], "2": [2*params["k"]+1], "3": []}
    if number == 3:
        k = params["k"]
        return {"1": [2], "2": [4*k, 2] if k % 2 == 0 else [8*k], "3": []}
    if number == 12:
        value = params["k"]//2
        return {"1": [], "2": [value] if value > 1 else [], "3": []}
    return reference["rows"][number-1]["groups_by_dimension"]


def compare_finite(task, table, reference):
    comparisons = []
    for source in task["spec"]["source_rows"]:
        for dimension, wanted in finite_expected(source, reference).items():
            d = int(dimension)
            pieces = []
            for q in (-4, -2, -1, 0):
                p = d-2-q  # paper dimension d = p+q+2, package k=d+1
                if p < 0 or (p == 0 and d in (1, 2)):
                    continue
                pieces.append(dict(p=p, q=q, invariants=cell(table, p, q)))
            entry = dict(task=task["id"], row=source["row"], parameters=source["parameters"],
                         dimension=d, paper_group=wanted, paper_order=prod(wanted),
                         computed_layers=pieces, classification_group_verified=False)
            if any(p["invariants"] is None for p in pieces):
                entry["status"] = "not_computed"
            else:
                values = sum((p["invariants"] for p in pieces), [])
                entry["computed_free_rank"] = values.count(0)
                entry["computed_graded_order"] = None if 0 in values else prod(values)
                entry["status"] = "order_matches_extensions_unresolved" if entry["computed_graded_order"] == prod(wanted) else "order_differs"
            comparisons.append(entry)
    return comparisons


def compare_space(task, table, reference):
    expected = reference["groups"][task["spec"]["number"]-1]
    result = dict(task=task["id"], group=task["spec"]["id"], twist=task["spacegroup_twist"],
                  paper_printed_total=expected["printed_total"],
                  paper_total_consistent=expected["printed_total_matches_layers"], layers={})
    for layer, (p, q) in reference["ahss_mapping"].items():
        actual = cell(table, p, q)
        wanted = expected["layers"][layer]
        if task["spacegroup_twist"] != "zero":
            wanted, status = None, "twist_not_tabulated"
        elif wanted is None:
            status = "unreported_by_paper"
        elif actual is None:
            status = "not_computed"
        else:
            status = "matches" if primary(actual) == primary(wanted) else "differs"
        result["layers"][layer] = dict(computed=actual, paper=wanted, status=status)
    statuses = [entry["status"] for entry in result["layers"].values()]
    if task["spacegroup_twist"] != "zero":
        result["status"] = "twist_not_tabulated"
    elif "differs" in statuses:
        result["status"] = "layer_discrepancy"
    elif "not_computed" in statuses:
        result["status"] = "not_computed"
    elif "unreported_by_paper" in statuses:
        result["status"] = "paper_incomplete"
    else:
        result["status"] = "layers_match"
    return result


def report(run):
    lock = (run/".report.lock").open("w")
    fcntl.flock(lock, fcntl.LOCK_EX)
    manifest = json.loads((run/"manifest.json").read_text())
    references = run/"source/references"
    finite_ref = json.loads((references/"table-iii.json").read_text())
    space_ref = json.loads((references/"spacegroups.json").read_text())
    queue = json.loads((run/"queue.json").read_text()).get("jobs", {}) if (run/"queue.json").exists() else {}
    records, finite, space = [], [], []
    for task in manifest["tasks"]:
        directory = run/"tasks"/task["id"]
        history = events(directory/"events.jsonl")
        terminal = directory/"metrics.json"
        live = directory/"metrics-live.json"
        metrics = json.loads((terminal if terminal.exists() else live).read_text()) if terminal.exists() or live.exists() else {}
        submission = json.loads((directory/"submission.json").read_text()) if (directory/"submission.json").exists() else {}
        job = submission.get("job_id")
        state = metrics.get("state", "not_submitted")
        if not terminal.exists() and job:
            state = queue.get(job, "no_longer_in_queue_without_terminal_metrics")
        pages = {e["page"]: e["table"] for e in history if e.get("type") == "page_complete"}
        last_start = next((e.get("page") for e in reversed(history) if e.get("type") == "page_start"), None)
        record = dict(task=task["id"], group=task["spec"]["id"], job_id=job, state=state,
                      completed_pages=sorted(pages), current_page=last_start,
                      metrics=metrics, tables=pages, certified_ko=False)
        records.append(record)
        if task["spec"]["kind"] == "spacegroup":
            space.append(compare_space(task, pages.get(6), space_ref))
        else:
            finite.extend(compare_finite(task, pages.get(6), finite_ref))
    counts = dict(Counter(r["state"] for r in records))
    summary = dict(updated_utc=now(), counts=counts, total_tasks=len(records),
                   completed_e6=sum(6 in r["completed_pages"] for r in records), certified_ko=False,
                   finite_comparison_counts=dict(Counter(r["status"] for r in finite)),
                   space_comparison_counts=dict(Counter(r["status"] for r in space)))
    save(run/"summary.json", summary)
    save(run/"results.json", records)
    save(run/"comparison.json", dict(table_iii=finite, spacegroups=space))
    with (run/"resources.csv").open("w") as stream:
        fields = ["task", "group", "job_id", "state", "host", "wall_seconds",
                  "allocated_memory_mib", "gap_workspace_limit_mib", "gnu_time_max_rss_kib",
                  "sampled_tree_peak_rss_kib", "cgroup_peak_bytes"]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for record in records:
            value = {**record, **record["metrics"]}
            value["cgroup_peak_bytes"] = (record["metrics"].get("cgroup_peak") or {}).get("bytes")
            writer.writerow({key: value.get(key) for key in fields})
    lines = [f"Updated {summary['updated_utc']}.", "",
             f"{summary['completed_e6']} / {len(records)} calculations have saved E6 tables.",
             "Task states: " + json.dumps(counts, sort_keys=True), "",
             "These are five-row associated-graded calculations, certified_ko=false. Stacking extensions are unresolved.",
             "Paper spatial dimension d uses p+q=d-2 (package cutoff k=d+1).",
             "Table III omits intrinsic p=0 Kitaev/p+ip factors in 1D/2D; comparisons test phase counts only.", "",
             "Table III comparison: " + json.dumps(summary["finite_comparison_counts"], sort_keys=True),
             "Space-group comparison: " + json.dumps(summary["space_comparison_counts"], sort_keys=True), "",
             "Sources: [Table III](https://arxiv.org/pdf/1811.00536#page=8), "
             "[Tables I–II](https://arxiv.org/pdf/2512.25069#page=6).",
             "The second space-group twist has no reference in these tables. SG210, SG219, SG228 contain unreported cells.",
             "Paper rows whose printed total differs from their printed layers: " +
             ", ".join(str(g["number"]) for g in space_ref["groups"] if g["printed_total_matches_layers"] is False)+".", "",
             "Measured wall time includes setup, resolution, twists, and pages; it excludes queue wait. "
             "GAP event runtime_ms is cumulative GAP CPU time per case, excluding external-worker CPU. "
             "GNU time RSS is the largest process high-water mark; sampled tree RSS sums live descendants every 2 s. "
             "Cgroup peak, when available, measures the job cgroup and includes process overhead.", "",
             "Completed results and discrepancies:", ""]
    for entry in finite:
        if entry["status"] != "not_computed":
            lines.append(f"- {entry['task']}, row {entry['row']} {entry['parameters']}, {entry['dimension']}D: "
                         f"{entry['status']}; graded order {entry['computed_graded_order']}, paper order {entry['paper_order']}.")
    for entry in space:
        if any(layer["computed"] is not None for layer in entry["layers"].values()):
            lines.append(f"- {entry['task']}: {entry['status']}; " + "; ".join(
                f"{name}={layer['computed']} (paper {layer['paper']}, {layer['status']})"
                for name, layer in entry["layers"].items()))
    if not any(r["completed_pages"] for r in records):
        lines.append("No page results yet.")
    (run/"REPORT.md").write_text("\n".join(lines)+"\n")
    print(json.dumps(summary, sort_keys=True), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    report(parser.parse_args().run.resolve())
