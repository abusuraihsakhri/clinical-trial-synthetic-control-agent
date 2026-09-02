"""
Unit and Integration Test Suite for Synthetic Control Arm (SCA) Engine
=====================================================================
Tests Propensity Score Modeling, Caliper Matching, Covariate Balance (SMD),
Kaplan-Meier Survival Estimation, Log-Rank Tests, and Regulatory Adequacy.
"""

import unittest
import json
import os
import sys
import math

# Ensure root directory is in sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from synthetic_control_arm import (
    MatchingMethod,
    RegulatoryAdequacyTier,
    SubjectRecord,
    MatchedPair,
    CovariateBalance,
    SurvivalCurvePoint,
    SyntheticControlAnalysisResult,
    BiostatisticalMath,
    PropensityScoreEngine,
    SurvivalAnalysisEngine,
    SyntheticControlAgentEngine,
    parse_synthetic_cohort_dict,
)


class TestBiostatisticalMath(unittest.TestCase):
    """Test statistical math calculations."""

    def test_mean_and_variance(self):
        data = [10.0, 20.0, 30.0, 40.0, 50.0]
        self.assertAlmostEqual(BiostatisticalMath.mean(data), 30.0)
        self.assertAlmostEqual(BiostatisticalMath.variance(data), 250.0)
        self.assertAlmostEqual(BiostatisticalMath.std_dev(data), math.sqrt(250.0))

    def test_empty_mean(self):
        self.assertEqual(BiostatisticalMath.mean([]), 0.0)
        self.assertEqual(BiostatisticalMath.variance([5.0]), 0.0)

    def test_smd_identical_distributions(self):
        t = [50.0, 55.0, 60.0]
        c = [50.0, 55.0, 60.0]
        smd = BiostatisticalMath.standardized_mean_difference(t, c)
        self.assertAlmostEqual(smd, 0.0)

    def test_smd_divergent_distributions(self):
        t = [70.0, 75.0, 80.0]
        c = [50.0, 55.0, 60.0]
        smd = BiostatisticalMath.standardized_mean_difference(t, c)
        self.assertGreater(smd, 2.0)

    def test_smd_empty_cases(self):
        self.assertEqual(BiostatisticalMath.standardized_mean_difference([], [1.0]), 0.0)
        self.assertEqual(BiostatisticalMath.standardized_mean_difference([1.0], []), 0.0)

    def test_sigmoid_limits(self):
        self.assertAlmostEqual(BiostatisticalMath.sigmoid(0.0), 0.5)
        self.assertEqual(BiostatisticalMath.sigmoid(50.0), 1.0)
        self.assertEqual(BiostatisticalMath.sigmoid(-50.0), 0.0)

    def test_chi2_survival_function(self):
        p_val = BiostatisticalMath.chi2_sf_1df(3.84146)
        self.assertAlmostEqual(p_val, 0.05, places=3)
        self.assertEqual(BiostatisticalMath.chi2_sf_1df(0.0), 1.0)
        self.assertEqual(BiostatisticalMath.chi2_sf_1df(-5.0), 1.0)


class TestPropensityScoreAndMatching(unittest.TestCase):
    """Test propensity score modeling and nearest neighbor matching."""

    def setUp(self):
        self.subjects: List[SubjectRecord] = []
        for i in range(10):
            self.subjects.append(
                SubjectRecord(
                    subject_id=f"T-{i}",
                    is_treated=True,
                    covariates={"age": 55.0 + i, "ecog": 1.0},
                    time_to_event_months=18.0,
                    event_observed=True,
                    response_achieved=True,
                )
            )
        for j in range(30):
            self.subjects.append(
                SubjectRecord(
                    subject_id=f"C-{j}",
                    is_treated=False,
                    covariates={"age": 52.0 + j, "ecog": 1.0 if j % 2 == 0 else 2.0},
                    time_to_event_months=12.0,
                    event_observed=True,
                    response_achieved=False,
                )
            )

    def test_propensity_score_bounds(self):
        PropensityScoreEngine.fit_and_predict_propensity_scores(self.subjects, ["age", "ecog"])
        for s in self.subjects:
            self.assertIsNotNone(s.propensity_score)
            self.assertGreater(s.propensity_score, 0.0)
            self.assertLess(s.propensity_score, 1.0)

    def test_matching_generates_pairs(self):
        res = SyntheticControlAgentEngine.generate_synthetic_control_arm(
            self.subjects, covariate_names=["age", "ecog"], caliper_sd_multiplier=0.5
        )
        self.assertGreater(res.matched_pairs_count, 0)
        self.assertEqual(res.trial_arm_size, 10)
        self.assertEqual(res.rwd_pool_size, 30)

    def test_post_matching_smd_improvement(self):
        res = SyntheticControlAgentEngine.generate_synthetic_control_arm(
            self.subjects, covariate_names=["age", "ecog"], caliper_sd_multiplier=0.5
        )
        self.assertLessEqual(res.mean_absolute_smd_post, res.mean_absolute_smd_pre + 0.1)

    def test_fit_and_predict_empty_covariates(self):
        subs = [SubjectRecord("T1", True, {}, 10.0, True)]
        out = PropensityScoreEngine.fit_and_predict_propensity_scores(subs, [])
        self.assertEqual(out[0].propensity_score, 0.5)


class TestSurvivalAnalysisEngine(unittest.TestCase):
    """Test Kaplan-Meier survival curves and Log-Rank tests."""

    def test_kaplan_meier_monotone_decreasing(self):
        subs = [
            SubjectRecord("S1", True, {}, 4.0, True),
            SubjectRecord("S2", True, {}, 8.0, True),
            SubjectRecord("S3", True, {}, 12.0, False),  # Censored
            SubjectRecord("S4", True, {}, 16.0, True),
        ]
        curve, med = SurvivalAnalysisEngine.calculate_kaplan_meier(subs)
        self.assertEqual(curve[0].survival_prob, 1.0)
        for i in range(1, len(curve)):
            self.assertLessEqual(curve[i].survival_prob, curve[i - 1].survival_prob)

    def test_kaplan_meier_empty(self):
        curve, med = SurvivalAnalysisEngine.calculate_kaplan_meier([])
        self.assertEqual(curve, [])
        self.assertEqual(med, 0.0)

    def test_log_rank_superiority(self):
        treated = [SubjectRecord(f"T{i}", True, {}, 24.0, True) for i in range(15)]
        control = [SubjectRecord(f"C{i}", False, {}, 6.0, True) for i in range(15)]

        hr, ci_l, ci_h, chi2, p_val = SurvivalAnalysisEngine.log_rank_test_and_hazard_ratio(treated, control)
        self.assertLess(hr, 1.0)
        self.assertLess(p_val, 0.05)
        self.assertGreater(chi2, 3.84)

    def test_log_rank_no_events(self):
        treated = [SubjectRecord(f"T{i}", True, {}, 24.0, False) for i in range(5)]
        control = [SubjectRecord(f"C{i}", False, {}, 12.0, False) for i in range(5)]
        hr, ci_l, ci_h, chi2, p_val = SurvivalAnalysisEngine.log_rank_test_and_hazard_ratio(treated, control)
        self.assertEqual(hr, 1.0)
        self.assertEqual(p_val, 1.0)


class TestEndToEndSCAAnalysis(unittest.TestCase):
    """Test end-to-end synthetic control arm pipeline and edge cases."""

    def test_complete_sca_analysis(self):
        subjects = []
        for i in range(20):
            subjects.append(
                SubjectRecord(
                    subject_id=f"T-{i}",
                    is_treated=True,
                    covariates={"age": 60.0 + (i % 5), "ecog": float(i % 2), "ldh": 200.0 + (i * 2)},
                    time_to_event_months=20.0 + i,
                    event_observed=True,
                    response_achieved=True,
                )
            )
        for j in range(60):
            subjects.append(
                SubjectRecord(
                    subject_id=f"C-{j}",
                    is_treated=False,
                    covariates={"age": 58.0 + (j % 10), "ecog": float(j % 3), "ldh": 220.0 + (j * 3)},
                    time_to_event_months=12.0 + j,
                    event_observed=True,
                    response_achieved=(j % 3 == 0),
                )
            )

        res = SyntheticControlAgentEngine.generate_synthetic_control_arm(subjects)
        self.assertIn(res.regulatory_adequacy_tier, [
            RegulatoryAdequacyTier.STRONG_REGULATORY_GRADE,
            RegulatoryAdequacyTier.ACCEPTABLE_SUPPORTIVE,
            RegulatoryAdequacyTier.HIGH_CONFOUNDING_RISK,
        ])
        self.assertGreater(len(res.covariate_balance), 0)
        self.assertGreater(len(res.recommendations), 0)

    def test_empty_treated_raises_value_error(self):
        controls_only = [SubjectRecord("C1", False, {"age": 60.0}, 10.0, True)]
        with self.assertRaises(ValueError):
            SyntheticControlAgentEngine.generate_synthetic_control_arm(controls_only)

    def test_parse_synthetic_cohort_dict(self):
        data = {
            "subjects": [
                {
                    "subject_id": "PT-01",
                    "is_treated": True,
                    "covariates": {"age": 58.0, "ecog": 0.0},
                    "time_to_event_months": 15.4,
                    "event_observed": True,
                    "response_achieved": True,
                }
            ]
        }
        subs = parse_synthetic_cohort_dict(data)
        self.assertEqual(len(subs), 1)
        self.assertEqual(subs[0].subject_id, "PT-01")
        self.assertTrue(subs[0].is_treated)
        self.assertEqual(subs[0].covariates["age"], 58.0)

    def test_json_export_structure(self):
        subjects = [
            SubjectRecord("T1", True, {"age": 55.0}, 15.0, True),
            SubjectRecord("C1", False, {"age": 56.0}, 10.0, True),
        ]
        res = SyntheticControlAgentEngine.generate_synthetic_control_arm(subjects)
        d = res.to_dict()
        self.assertIn("trial_arm_size", d)
        self.assertIn("efficacy_outcomes", d)
        self.assertIn("covariate_balance", d)
        self.assertIn("regulatory_adequacy_tier", d)


if __name__ == "__main__":
    unittest.main()
