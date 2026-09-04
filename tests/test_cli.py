"""
Integration and CLI Tests for Clinical Trial Synthetic Control Agent
=====================================================================
Tests CLI batch processing, CSV loading, demo generation, and report formatting.
"""

import unittest
import tempfile
import os
import csv
import sys

# Ensure root directory is in sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from cli import (
    load_cohort_from_csv,
    run_batch_matching,
    main,
    generate_benchmark_synthetic_cohort,
    format_sca_report,
)
from synthetic_control_arm import (
    SubjectRecord,
    SyntheticControlAgentEngine,
    RegulatoryAdequacyTier,
)


class TestCLIBatchAndWorkflows(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.csv_input = os.path.join(self.temp_dir.name, "test_cohort.csv")
        self.csv_output = os.path.join(self.temp_dir.name, "test_results.csv")

        # Create realistic test CSV
        with open(self.csv_input, mode="w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                "subject_id", "is_treated", "age", "ecog", "prior_lines", "ldh",
                "time_to_event_months", "event_observed", "response_achieved"
            ])
            # Treated
            for i in range(5):
                writer.writerow([f"TRIAL-{i+1:03d}", "True", 58.0 + i, 0.0, 1.0, 210.0 + i*5, 22.0 + i, "True", "True"])
            # Control
            for j in range(12):
                writer.writerow([f"RWD-{j+1:03d}", "False", 60.0 + (j%4), float(j%2), 1.0 + float(j%3), 220.0 + j*8, 15.0 + j, "True", "False"])

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_load_cohort_from_csv(self):
        subjects = load_cohort_from_csv(self.csv_input)
        self.assertEqual(len(subjects), 17)
        treated = [s for s in subjects if s.is_treated]
        controls = [s for s in subjects if not s.is_treated]
        self.assertEqual(len(treated), 5)
        self.assertEqual(len(controls), 12)
        self.assertIn("age", subjects[0].covariates)
        self.assertIn("ldh", subjects[0].covariates)

    def test_run_batch_matching_execution(self):
        ret = run_batch_matching(self.csv_input, self.csv_output, caliper=0.2)
        self.assertEqual(ret, 0)
        self.assertTrue(os.path.exists(self.csv_output))

        with open(self.csv_output, mode="r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        self.assertEqual(len(rows), 17)
        self.assertIn("cohort_arm", rows[0])
        self.assertIn("propensity_score", rows[0])
        self.assertIn("hazard_ratio", rows[0])
        self.assertIn("regulatory_tier", rows[0])

    def test_cli_main_batch_subcommand(self):
        ret = main(["batch", "-i", self.csv_input, "-o", self.csv_output])
        self.assertEqual(ret, 0)
        self.assertTrue(os.path.exists(self.csv_output))

    def test_cli_main_demo(self):
        ret = main(["--demo"])
        self.assertEqual(ret, 0)

    def test_format_sca_report(self):
        subjects = generate_benchmark_synthetic_cohort(10, 30)
        res = SyntheticControlAgentEngine.generate_synthetic_control_arm(subjects)
        report = format_sca_report(res)
        self.assertIn("SYNTHETIC CONTROL ARM (SCA)", report)
        self.assertIn("COVARIATE BALANCE DIAGNOSTICS", report)
        self.assertIn("COMPARATIVE CLINICAL EFFICACY", report)


if __name__ == "__main__":
    unittest.main()
