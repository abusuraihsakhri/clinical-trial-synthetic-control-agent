import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from synthetic_control_arm.cli import (
    generate_benchmark_synthetic_cohort,
    load_cohort_from_csv,
    main,
    run_batch_matching,
)


class TestCLI(unittest.TestCase):
    def test_demo_is_deterministic(self):
        first = generate_benchmark_synthetic_cohort(3, 4)
        second = generate_benchmark_synthetic_cohort(3, 4)
        self.assertEqual(
            [(s.subject_id, s.covariates) for s in first],
            [(s.subject_id, s.covariates) for s in second],
        )

    def test_sample_csv_loads(self):
        subjects = load_cohort_from_csv(Path(__file__).parents[1] / "sample.csv")
        self.assertEqual(len(subjects), 15)
        self.assertEqual(sum(subject.is_treated for subject in subjects), 5)

    def test_batch_output_contains_real_match_metadata(self):
        sample = Path(__file__).parents[1] / "sample.csv"
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "results.csv"
            self.assertEqual(run_batch_matching(sample, output, quiet=True), 0)
            with output.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            matched = [row for row in rows if row["match_status"] == "MATCHED"]
            self.assertTrue(matched)
            self.assertTrue(all(row["matched_pair_id"] for row in matched))
            self.assertTrue(all(row["matched_subject_id"] for row in matched))
            self.assertTrue(all(row["logit_ps_distance"] for row in matched))

    def test_invalid_boolean_in_csv_returns_clear_error(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.csv"
            path.write_text(
                "subject_id,is_treated,age,time_to_event_months,event_observed,response_achieved\n"
                "T1,maybe,60,10,true,false\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "is_treated"):
                load_cohort_from_csv(path)

    def test_package_module_help_works(self):
        result = subprocess.run(
            [sys.executable, "-m", "synthetic_control_arm.cli", "--help"],
            cwd=Path(__file__).parents[1],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Propensity-score matching", result.stdout)

    def test_json_demo_is_valid_json(self):
        result = subprocess.run(
            [sys.executable, "-m", "synthetic_control_arm.cli", "--demo", "--json"],
            cwd=Path(__file__).parents[1],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        payload = json.loads(result.stdout)
        self.assertIn("balance_assessment", payload)
        self.assertIn("matched_pairs", payload)

    def test_main_rejects_unknown_file_type(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "cohort.txt"
            path.write_text("x", encoding="utf-8")
            self.assertEqual(main(["--file", str(path)]), 2)


    def test_duplicate_normalized_csv_headers_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "duplicate.csv"
            path.write_text(
                "subject_id,is_treated,age, age,time_to_event_months,event_observed,response_achieved\n"
                "T1,true,60,61,10,true,false\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "unique after trimming"):
                load_cohort_from_csv(path)

    def test_blank_csv_header_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "blank-header.csv"
            path.write_text(
                "subject_id,is_treated,,time_to_event_months,event_observed,response_achieved\n"
                "T1,true,60,10,true,false\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "non-empty"):
                load_cohort_from_csv(path)


if __name__ == "__main__":
    unittest.main()
