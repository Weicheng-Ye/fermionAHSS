"""Tests for quotient groups and provenance in the batch catalogue."""

import json
import unittest

try:
    from .catalogue import build_catalogue
except ImportError:
    from catalogue import build_catalogue


class CatalogueTests(unittest.TestCase):
    def test_default_rows_and_space_group_numbers(self):
        catalogue = build_catalogue()
        finite = [entry for entry in catalogue if entry["kind"] == "abelian"]
        space_groups = [entry for entry in catalogue
                        if entry["kind"] == "spacegroup"]
        self.assertEqual(len(finite), 17)
        self.assertEqual(
            {row["row"] for entry in finite for row in entry["source_rows"]},
            set(range(1, 22)),
        )
        self.assertEqual([entry["number"] for entry in space_groups],
                         list(range(1, 231)))
        self.assertEqual([entry["id"] for entry in space_groups],
                         [f"SG{number:03d}" for number in range(1, 231)])
        self.assertEqual(len({entry["id"] for entry in catalogue}), len(catalogue))
        self.assertEqual(json.loads(json.dumps(catalogue)), catalogue)

    def test_selected_families_retain_every_parameter(self):
        catalogue = build_catalogue(k_values=[1, 5], n_values=[1, 5])
        rows = [(entry, row) for entry in catalogue if entry["kind"] == "abelian"
                for row in entry["source_rows"]]
        self.assertEqual({row["row"] for _, row in rows}, set(range(1, 22)))
        self.assertEqual(
            [(row["parameters"]["k"], entry["factors"])
             for entry, row in rows if row["row"] == 2],
            [(1, [3]), (5, [11])],
        )
        self.assertEqual(
            [(row["parameters"]["k"], entry["factors"])
             for entry, row in rows if row["row"] == 3],
            [(1, [2]), (5, [10])],
        )
        cyclic_fermionic = [(entry, row) for entry, row in rows if row["row"] == 12]
        self.assertEqual(
            [(row["parameters"], entry["factors"], row["fermionic_group"])
             for entry, row in cyclic_fermionic],
            [({"n": 1, "k": 2}, [2], "C4^f"),
             ({"n": 5, "k": 32}, [32], "C64^f")],
        )

    def test_quotient_deduplication_preserves_distinct_extensions(self):
        catalogue = build_catalogue([1, 1], [1, 1])
        entries = {entry["id"]: entry for entry in catalogue}
        c2_rows = entries["table3-C2"]["source_rows"]
        self.assertEqual([row["row"] for row in c2_rows], [1, 3, 12, 19, 20])
        self.assertEqual([(row["paper_s"], row["paper_omega"])
                          for row in c2_rows if row["row"] in (19, 20)],
                         [("x", "0"), ("x", "x^2")])
        klein_rows = entries["table3-C2xC2"]["source_rows"]
        self.assertEqual([row["row"] for row in klein_rows], [4, 13, 21])
        quaternion = next(row for row in klein_rows if row["row"] == 21)
        self.assertEqual(quaternion["fermionic_group"], "Q8^f")
        self.assertEqual(quaternion["paper_omega"], "x^2 + x*y + y^2")

    def test_carry_factor_survives_sorting(self):
        entry = next(entry for entry in build_catalogue()
                     if entry["id"] == "table3-C2xC4")
        rows = {row["row"]: row for row in entry["source_rows"]}
        self.assertEqual(entry["factors"], [2, 4])
        self.assertEqual(rows[14]["quotient_factors_in_paper_order"], [2, 4])
        self.assertEqual(rows[14]["paper_omega_cyclic_order"], 2)
        self.assertEqual(rows[15]["quotient_factors_in_paper_order"], [4, 2])
        self.assertEqual(rows[15]["paper_omega_cyclic_order"], 4)
        self.assertIn("first paper-order C4 factor", rows[15]["paper_omega"])
        self.assertIn("floor((a0+b0)/4)", rows[15]["paper_omega"])

    def test_paper_metadata_contains_only_listed_twist_recipes(self):
        entries = {entry["id"]: entry for entry in build_catalogue()}
        self.assertEqual([row["row"] for row in entries["table3-C6"]["source_rows"]], [3])
        self.assertEqual([row["row"] for row in entries["table3-C2xC4"]["source_rows"]],
                         [5, 14, 15])
        rows = [row for entry in entries.values() if entry["kind"] == "abelian"
                for row in entry["source_rows"]]
        self.assertEqual({row["row"] for row in rows if row["paper_s"] != "0"}, {19, 20})
        self.assertEqual({row["row"] for row in rows
                          if row["paper_omega_kind"] == "cyclic_carry"},
                         set(range(12, 19)) | {20})
        self.assertEqual({row["row"] for row in rows
                          if row["paper_omega_kind"] == "quaternion"}, {21})
        self.assertTrue(all(row["paper_omega_kind"] in ("zero", "cyclic_carry", "quaternion")
                            for row in rows))

    def test_space_groups_have_exactly_the_two_requested_labeled_twists(self):
        groups = [entry for entry in build_catalogue() if entry["kind"] == "spacegroup"]
        for entry in groups:
            with self.subTest(group=entry["id"]):
                twists = entry["spacegroup_twists"]
                self.assertEqual([(twist["spacegroup_twist"], twist["s"], twist["omega"])
                                  for twist in twists],
                                 [("zero", "w1", "0"),
                                  ("w2_plus_w1_squared", "w1", "w2+w1^2")])
        groups[0]["spacegroup_twists"][0]["s"] = "changed"
        self.assertEqual(groups[1]["spacegroup_twists"][0]["s"], "w1")
        fresh = next(entry for entry in build_catalogue() if entry["id"] == "SG001")
        self.assertEqual(fresh["spacegroup_twists"][0]["s"], "w1")

    def test_parameter_validation(self):
        for value in (0, -1, True, 1.0, "1", None):
            for parameter in ("k_values", "n_values"):
                with self.subTest(parameter=parameter, value=value):
                    with self.assertRaises(ValueError):
                        build_catalogue(**{parameter: [value]})
        for parameter in ("k_values", "n_values"):
            with self.subTest(parameter=parameter):
                with self.assertRaises(TypeError):
                    build_catalogue(**{parameter: 1})

    def test_empty_families_are_explicitly_omitted(self):
        rows = {row["row"] for entry in build_catalogue([], [])
                if entry["kind"] == "abelian" for row in entry["source_rows"]}
        self.assertEqual(rows, set(range(1, 22)) - {2, 3, 12})

    def test_each_call_returns_independent_metadata(self):
        catalogue = build_catalogue()
        catalogue[0]["factors"].append(99)
        catalogue[0]["source_rows"][0]["parameters"]["changed"] = True
        fresh = build_catalogue()
        self.assertNotIn(99, fresh[0]["factors"])
        self.assertNotIn("changed", fresh[0]["source_rows"][0]["parameters"])


if __name__ == "__main__":
    unittest.main()
