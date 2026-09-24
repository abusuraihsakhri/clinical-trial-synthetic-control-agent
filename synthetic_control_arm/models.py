"""Data models for synthetic-control matching and outcome diagnostics."""

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional


class MatchingMethod(str, Enum):
    """Matching methods exposed by the public API."""

    NEAREST_NEIGHBOR_CALIPER = "NEAREST_NEIGHBOR_CALIPER"
    EXACT_AND_CALIPER = "EXACT_AND_CALIPER"
    IPTW_WEIGHTING = "IPTW_WEIGHTING"


class RegulatoryAdequacyTier(str, Enum):
    """Legacy compatibility values for the tool's heuristic balance assessment.

    The names are retained to avoid breaking existing consumers. They are not FDA,
    EMA, or other regulator-defined grades and must not be interpreted as approval,
    validation, or regulatory acceptability.
    """

    STRONG_REGULATORY_GRADE = "STRONG_REGULATORY_GRADE"
    ACCEPTABLE_SUPPORTIVE = "ACCEPTABLE_SUPPORTIVE"
    HIGH_CONFOUNDING_RISK = "HIGH_CONFOUNDING_RISK"


BALANCE_ASSESSMENT_LABELS = {
    RegulatoryAdequacyTier.STRONG_REGULATORY_GRADE: "WELL_BALANCED",
    RegulatoryAdequacyTier.ACCEPTABLE_SUPPORTIVE: "PARTIAL_BALANCE",
    RegulatoryAdequacyTier.HIGH_CONFOUNDING_RISK: "POOR_BALANCE_OR_OVERLAP",
}


@dataclass
class SubjectRecord:
    subject_id: str
    is_treated: bool
    covariates: Dict[str, float]
    time_to_event_months: float
    event_observed: bool
    response_achieved: bool = False
    propensity_score: Optional[float] = None
    iptw_weight: Optional[float] = None


@dataclass
class MatchedPair:
    treated_id: str
    control_id: str
    ps_distance: float
    treated_ps: float
    control_ps: float


@dataclass
class CovariateBalance:
    covariate_name: str
    pre_treated_mean: float
    pre_control_mean: float
    pre_smd: float
    post_treated_mean: float
    post_control_mean: float
    post_smd: float
    is_balanced: bool


@dataclass
class SurvivalCurvePoint:
    time_months: float
    n_at_risk: int
    events: int
    censored: int
    survival_prob: float


@dataclass
class SyntheticControlAnalysisResult:
    trial_arm_size: int
    rwd_pool_size: int
    matched_pairs_count: int
    matching_method: MatchingMethod
    caliper_width: float
    covariate_balance: List[CovariateBalance]
    all_covariates_balanced: bool
    mean_absolute_smd_pre: float
    mean_absolute_smd_post: float
    hazard_ratio: Optional[float]
    hazard_ratio_ci_low: Optional[float]
    hazard_ratio_ci_high: Optional[float]
    log_rank_test_statistic: float
    log_rank_p_value: float
    median_survival_treated_months: Optional[float]
    median_survival_synthetic_control_months: Optional[float]
    orr_treated_pct: float
    orr_synthetic_control_pct: float
    att_orr_diff_pct: float
    regulatory_adequacy_tier: RegulatoryAdequacyTier
    recommendations: List[str]
    matched_pairs: List[MatchedPair] = field(default_factory=list)
    matching_retention_pct: float = 0.0

    @staticmethod
    def _round_optional(value: Optional[float], digits: int) -> Optional[float]:
        if value is None or not math.isfinite(value):
            return None
        return round(value, digits)

    @property
    def balance_assessment(self) -> str:
        return BALANCE_ASSESSMENT_LABELS[self.regulatory_adequacy_tier]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "trial_arm_size": self.trial_arm_size,
            "rwd_pool_size": self.rwd_pool_size,
            "matched_pairs_count": self.matched_pairs_count,
            "matching_retention_pct": round(self.matching_retention_pct, 1),
            "matching_method": self.matching_method.value,
            "caliper_width_logit_sd": round(self.caliper_width, 4),
            "mean_absolute_smd_pre": self._round_optional(self.mean_absolute_smd_pre, 3),
            "mean_absolute_smd_post": self._round_optional(self.mean_absolute_smd_post, 3),
            "all_covariates_balanced": self.all_covariates_balanced,
            "balance_assessment": self.balance_assessment,
            "regulatory_adequacy_tier": self.regulatory_adequacy_tier.value,
            "efficacy_outcomes": {
                "hazard_ratio_logrank_approx": self._round_optional(self.hazard_ratio, 3),
                "hazard_ratio_95_ci": [
                    self._round_optional(self.hazard_ratio_ci_low, 3),
                    self._round_optional(self.hazard_ratio_ci_high, 3),
                ],
                "log_rank_test_statistic": round(self.log_rank_test_statistic, 3),
                "log_rank_p_value": round(self.log_rank_p_value, 4),
                "median_os_treated_months": self._round_optional(self.median_survival_treated_months, 1),
                "median_os_synthetic_control_months": self._round_optional(self.median_survival_synthetic_control_months, 1),
                "orr_treated_pct": round(self.orr_treated_pct, 1),
                "orr_synthetic_control_pct": round(self.orr_synthetic_control_pct, 1),
                "att_orr_difference_pct": round(self.att_orr_diff_pct, 1),
            },
            "covariate_balance": [
                {
                    "covariate": cb.covariate_name,
                    "pre_smd": self._round_optional(cb.pre_smd, 3),
                    "post_smd": self._round_optional(cb.post_smd, 3),
                    "pre_treated_mean": round(cb.pre_treated_mean, 2),
                    "pre_control_mean": round(cb.pre_control_mean, 2),
                    "post_treated_mean": round(cb.post_treated_mean, 2),
                    "post_control_mean": round(cb.post_control_mean, 2),
                    "is_balanced": cb.is_balanced,
                }
                for cb in self.covariate_balance
            ],
            "matched_pairs": [
                {
                    "treated_id": pair.treated_id,
                    "control_id": pair.control_id,
                    "logit_ps_distance": round(pair.ps_distance, 5),
                    "treated_ps": round(pair.treated_ps, 5),
                    "control_ps": round(pair.control_ps, 5),
                }
                for pair in self.matched_pairs
            ],
            "recommendations": self.recommendations,
        }
