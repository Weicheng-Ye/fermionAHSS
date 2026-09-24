import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import subprocess

from campaign import allocated_mib, build_catalogue, events, tasks_for, submission_environment
from report import cell, compare_finite, compare_space, finite_expected, primary

ROOT = Path(__file__).resolve().parents[1]


class CampaignTests(unittest.TestCase):
    def setUp(self):
        self.tasks = tasks_for(build_catalogue())
        self.finite_ref = json.loads((ROOT/"references/table-iii.json").read_text())
        self.space_ref = json.loads((ROOT/"references/spacegroups.json").read_text())

    def test_every_case_and_source_row_preserved(self):
        self.assertEqual(len(self.tasks), 489)
        self.assertEqual(len({t["id"] for t in self.tasks}), 489)
        finite = [t for t in self.tasks if t["spec"]["kind"] == "abelian"]
        self.assertEqual(len(finite), 29)
        self.assertEqual({r["row"] for t in finite for r in t["spec"]["source_rows"]}, set(range(1, 22)))
        for n in range(1, 231):
            self.assertEqual({t["spacegroup_twist"] for t in self.tasks if t["spec"].get("number") == n},
                             {"zero", "w2_plus_w1_squared"})

    def test_comparison_does_not_claim_extensions(self):
        task = next(t for t in self.tasks if t["id"] == "table3-C2-row1")
        table = [[[] for _ in range(n)] for n in (6, 5, 4, 3, 2)]
        # 2D C2: three Z2 layers have order 8 but are not proof of Z8.
        table[0][4] = [2]
        table[2][2] = [2]
        table[3][1] = [2]
        table[4][0] = [0]  # intrinsic p+ip layer must be omitted
        result = next(r for r in compare_finite(task, table, self.finite_ref)
                      if r["row"] == 1 and r["dimension"] == 2)
        self.assertEqual(result["status"], "order_matches_extensions_unresolved")
        self.assertFalse(result["classification_group_verified"])
        table[2][2] = {"status": "unresolved"}
        result = next(r for r in compare_finite(task, table, self.finite_ref)
                      if r["row"] == 1 and r["dimension"] == 2)
        self.assertEqual(result["status"], "not_computed")

    def test_parameterized_reference(self):
        def row(n, k):
            return dict(row=n, parameters={"k": k})
        self.assertEqual(finite_expected(row(3, 2), self.finite_ref)["2"], [8, 2])
        self.assertEqual(finite_expected(row(3, 3), self.finite_ref)["2"], [24])
        self.assertEqual(finite_expected(row(12, 2), self.finite_ref)["2"], [])

    def test_paper_coverage_and_twists(self):
        self.assertEqual(len(self.space_ref["groups"]), 230)
        missing = {g["number"] for g in self.space_ref["groups"] if any(v is None for v in g["layers"].values())}
        self.assertEqual(missing, {210, 219, 228})
        task = next(t for t in self.tasks if t["id"] == "SG001-w2_plus_w1_squared")
        self.assertEqual(compare_space(task, None, self.space_ref)["status"], "twist_not_tabulated")
        self.assertEqual(primary([6, 4]), primary([2, 3, 4]))

    def test_interrupted_events_not_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/"events.jsonl"
            p.write_text('{"type":"page_start"}\n{"type":')
            self.assertEqual(events(p), [{"type": "page_start"}])
        self.assertIsNone(cell(None, 1, 0))

    def test_full_node_memory_ignores_dispatcher_environment(self):
        with patch.dict("os.environ", {"SLURM_MEM_PER_NODE": "1024", "SLURMD_NODENAME": "phoenix13"}):
            with patch("campaign.subprocess.run", return_value=subprocess.CompletedProcess([], 0, "NodeName=phoenix13 RealMemory=768000 AllocMem=768000")):
                self.assertEqual(allocated_mib(), 768000)
            self.assertNotIn("SLURM_MEM_PER_NODE", submission_environment())


if __name__ == "__main__":
    unittest.main()
