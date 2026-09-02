"""
End-to-End Regulatory Adequacy and Propensity Matching Tests
"""

import unittest
import os
import sys

# Ensure root directory is in sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from synthetic_control_arm import (
    MatchingMethod,
    RegulatoryAdequacyTier,
    SubjectRecord,
    MatchedPair,
    CovariateBalance,
    SyntheticControlAgentEngine,
    parse_synthetic_cohort_dict,
)


class TestSyntheticControlAgentFull(unittest.TestCase):

    def test_multi_covariate_matching(self):
        treated = [
            SubjectRecord(f"T-{i}", True, {"age": 60.0 + i, "ecog": 1.0, "prior_lines": 2.0}, 18.0, True, True)
            for i in range(15)
        ]
        control = [
            SubjectRecord(f"C-{j}", False, {"age": 55.0 + (j % 15), "ecog": float(j % 2), "prior_lines": float((j % 3) + 1)}, 10.0, True, False)
            for j in range(45)
        ]
        res = SyntheticControlAgentEngine.generate_synthetic_control_arm(treated + control)
        self.assertEqual(res.trial_arm_size, 15)
        self.assertEqual(res.rwd_pool_size, 45)
        self.assertGreater(res.matched_pairs_count, 0)
        self.assertEqual(len(res.covariate_balance), 3)

    def test_strong_regulatory_grade_cohort(self):
        # Closely matched distribution
        treated = [
            SubjectRecord(f"T-{i}", True, {"age": 60.0, "ecog": 1.0}, 24.0, True, True)
            for i in range(20)
        ]
        control = [
            SubjectRecord(f"C-{j}", False, {"age": 60.0, "ecog": 1.0}, 12.0, True, False)
            for j in range(40)
        ]
        res = SyntheticControlAgentEngine.generate_synthetic_control_arm(treated + control)
        self.assertEqual(res.regulatory_adequacy_tier, RegulatoryAdequacyTier.STRONG_REGULATORY_GRADE)
        self.assertTrue(res.all_covariates_balanced)
        self.assertLess(res.mean_absolute_smd_post, 0.05)

    def test_high_confounding_risk_tier(self):
        # Drastically disparate distributions with little overlap
        treated = [
            SubjectRecord(f"T-{i}", True, {"age": 30.0, "ecog": 0.0}, 30.0, True, True)
            for i in range(10)
        ]
        control = [
            SubjectRecord(f"C-{j}", False, {"age": 85.0, "ecog": 3.0}, 5.0, True, False)
            for j in range(30)
        ]
        res = SyntheticControlAgentEngine.generate_synthetic_control_arm(treated + control, caliper_sd_multiplier=0.01)
        self.assertIn(res.regulatory_adequacy_tier, [
            RegulatoryAdequacyTier.HIGH_CONFOUNDING_RISK,
            RegulatoryAdequacyTier.ACCEPTABLE_SUPPORTIVE
        ])

    def test_zero_variance_covariates(self):
        treated = [
            SubjectRecord(f"T-{i}", True, {"constant_val": 1.0}, 15.0, True)
            for i in range(5)
        ]
        control = [
            SubjectRecord(f"C-{j}", False, {"constant_val": 1.0}, 10.0, True)
            for j in range(10)
        ]
        res = SyntheticControlAgentEngine.generate_synthetic_control_arm(treated + control)
        self.assertEqual(res.covariate_balance[0].post_smd, 0.0)


if __name__ == "__main__":
    unittest.main()
