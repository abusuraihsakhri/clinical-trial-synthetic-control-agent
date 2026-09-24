import json
import math
import unittest

from synthetic_control_arm import (
    BiostatisticalMath,
    MatchingMethod,
    PropensityScoreEngine,
    RegulatoryAdequacyTier,
    SubjectRecord,
    SurvivalAnalysisEngine,
    SyntheticControlAgentEngine,
    parse_synthetic_cohort_dict,
)


class TestBiostatisticalMath(unittest.TestCase):
    def test_mean_variance_and_chi_square_tail(self):
        values = [10.0, 20.0, 30.0, 40.0, 50.0]
        self.assertAlmostEqual(BiostatisticalMath.mean(values), 30.0)
        self.assertAlmostEqual(BiostatisticalMath.variance(values), 250.0)
        self.assertAlmostEqual(
            BiostatisticalMath.chi2_sf_1df(3.84146), 0.05, places=3
        )

    def test_smd_identical_constant_groups_is_zero(self):
        self.assertEqual(
            BiostatisticalMath.standardized_mean_difference(
                [60.0, 60.0], [60.0, 60.0]
            ),
            0.0,
        )

    def test_smd_separated_constant_groups_is_infinite(self):
        value = BiostatisticalMath.standardized_mean_difference(
            [30.0, 30.0], [85.0, 85.0]
        )
        self.assertTrue(math.isinf(value))

    def test_smd_empty_group_is_not_interpretable(self):
        self.assertTrue(
            math.isinf(
                BiostatisticalMath.standardized_mean_difference([], [1.0])
            )
        )


class TestPropensityAndMatching(unittest.TestCase):
    @staticmethod
    def overlapping_subjects():
        subjects = []
        for i in range(10):
            subjects.append(
                SubjectRecord(
                    f"T-{i}",
                    True,
                    {"age": 55.0 + i, "ecog": float(i % 2)},
                    18.0 + i,
                    i % 3 != 0,
                    i % 2 == 0,
                )
            )
        for i in range(30):
            subjects.append(
                SubjectRecord(
                    f"C-{i}",
                    False,
                    {"age": 52.0 + (i % 16), "ecog": float(i % 2)},
                    11.0 + i / 2,
                    i % 4 != 0,
                    i % 4 == 0,
                )
            )
        return subjects

    def test_propensity_scores_are_bounded(self):
        subjects = self.overlapping_subjects()
        PropensityScoreEngine.fit_and_predict_propensity_scores(
            subjects, ["age", "ecog"]
        )
        for subject in subjects:
            self.assertIsNotNone(subject.propensity_score)
            self.assertGreater(subject.propensity_score, 0.0)
            self.assertLess(subject.propensity_score, 1.0)

    def test_matching_distance_is_on_logit_scale(self):
        result = SyntheticControlAgentEngine.generate_synthetic_control_arm(
            self.overlapping_subjects(), caliper_sd_multiplier=0.5
        )
        self.assertGreater(result.matched_pairs_count, 0)
        for pair in result.matched_pairs:
            expected = abs(
                BiostatisticalMath.logit(pair.treated_ps)
                - BiostatisticalMath.logit(pair.control_ps)
            )
            self.assertAlmostEqual(pair.ps_distance, expected, places=10)
            self.assertLessEqual(pair.ps_distance, result.caliper_width + 1e-12)

    def test_duplicate_ids_are_rejected(self):
        subjects = [
            SubjectRecord("A", True, {"age": 60.0}, 10.0, True),
            SubjectRecord("A", False, {"age": 60.0}, 10.0, True),
        ]
        with self.assertRaisesRegex(ValueError, "unique"):
            SyntheticControlAgentEngine.generate_synthetic_control_arm(subjects)

    def test_missing_covariate_is_rejected_when_explicitly_requested(self):
        subjects = [
            SubjectRecord("T", True, {"age": 60.0}, 10.0, True),
            SubjectRecord("C", False, {"ecog": 1.0}, 10.0, True),
        ]
        with self.assertRaisesRegex(ValueError, "missing covariates"):
            SyntheticControlAgentEngine.generate_synthetic_control_arm(
                subjects, covariate_names=["age"]
            )

    def test_invalid_caliper_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "greater than zero"):
            SyntheticControlAgentEngine.generate_synthetic_control_arm(
                self.overlapping_subjects(), caliper_sd_multiplier=0
            )

    def test_unimplemented_matching_methods_do_not_silently_fall_back(self):
        with self.assertRaises(NotImplementedError):
            SyntheticControlAgentEngine.generate_synthetic_control_arm(
                self.overlapping_subjects(),
                matching_method=MatchingMethod.IPTW_WEIGHTING,
            )

    def test_no_matches_raises_instead_of_using_unmatched_controls(self):
        subjects = [
            SubjectRecord(f"T-{i}", True, {"age": 20.0 + i}, 20.0, True)
            for i in range(4)
        ] + [
            SubjectRecord(f"C-{i}", False, {"age": 80.0 + i}, 8.0, True)
            for i in range(4)
        ]
        with self.assertRaisesRegex(ValueError, "No treated-control pairs"):
            SyntheticControlAgentEngine.generate_synthetic_control_arm(
                subjects, caliper_sd_multiplier=1e-6
            )

    def test_infinite_smd_serializes_as_json_null(self):
        subjects = [
            SubjectRecord(f"T-{i}", True, {"age": 30.0}, 20.0, True)
            for i in range(3)
        ] + [
            SubjectRecord(f"C-{i}", False, {"age": 85.0}, 8.0, True)
            for i in range(3)
        ]
        result = SyntheticControlAgentEngine.generate_synthetic_control_arm(
            subjects, caliper_sd_multiplier=100.0
        )
        payload = result.to_dict()
        self.assertIsNone(payload["mean_absolute_smd_post"])
        json.dumps(payload, allow_nan=False)

    def test_result_uses_neutral_balance_label(self):
        result = SyntheticControlAgentEngine.generate_synthetic_control_arm(
            self.overlapping_subjects(), caliper_sd_multiplier=0.5
        )
        self.assertIn(
            result.balance_assessment,
            {"WELL_BALANCED", "PARTIAL_BALANCE", "POOR_BALANCE_OR_OVERLAP"},
        )
        self.assertIn(
            result.regulatory_adequacy_tier,
            set(RegulatoryAdequacyTier),
        )


class TestSurvivalAnalysis(unittest.TestCase):
    def test_kaplan_meier_median_not_reached_is_none(self):
        subjects = [
            SubjectRecord("S1", True, {"x": 1.0}, 12.0, True),
            SubjectRecord("S2", True, {"x": 1.0}, 20.0, False),
            SubjectRecord("S3", True, {"x": 1.0}, 30.0, False),
        ]
        _, median = SurvivalAnalysisEngine.calculate_kaplan_meier(subjects)
        self.assertIsNone(median)

    def test_no_events_returns_non_estimable_hr(self):
        treated = [SubjectRecord("T", True, {"x": 1.0}, 10.0, False)]
        control = [SubjectRecord("C", False, {"x": 1.0}, 10.0, False)]
        hr, low, high, chi2, p_value = (
            SurvivalAnalysisEngine.log_rank_test_and_hazard_ratio(treated, control)
        )
        self.assertIsNone(hr)
        self.assertIsNone(low)
        self.assertIsNone(high)
        self.assertEqual(chi2, 0.0)
        self.assertEqual(p_value, 1.0)


class TestJSONParsing(unittest.TestCase):
    def test_string_false_is_parsed_as_false(self):
        subjects = parse_synthetic_cohort_dict(
            {
                "subjects": [
                    {
                        "subject_id": "C-1",
                        "is_treated": "false",
                        "covariates": {"age": 62},
                        "time_to_event_months": 12,
                        "event_observed": "false",
                        "response_achieved": "0",
                    }
                ]
            }
        )
        self.assertFalse(subjects[0].is_treated)
        self.assertFalse(subjects[0].event_observed)
        self.assertFalse(subjects[0].response_achieved)

    def test_invalid_boolean_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "boolean"):
            parse_synthetic_cohort_dict(
                {
                    "subjects": [
                        {
                            "subject_id": "X",
                            "is_treated": "sometimes",
                            "covariates": {"age": 60},
                            "time_to_event_months": 12,
                            "event_observed": True,
                        }
                    ]
                }
            )


if __name__ == "__main__":
    unittest.main()
