"""
Data Models & Definitions for Synthetic Control Arm (SCA) Generation.
Domain: Real-World Evidence (RWE), Propensity Score Matching & Survival Analysis
Standards: FDA / EMA RWD Guidance, ISPOR-ISPE Good Research Practices
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Any, Tuple


class MatchingMethod(str, Enum):
    NEAREST_NEIGHBOR_CALIPER = "NEAREST_NEIGHBOR_CALIPER"
    EXACT_AND_CALIPER = "EXACT_AND_CALIPER"
    IPTW_WEIGHTING = "IPTW_WEIGHTING"


class RegulatoryAdequacyTier(str, Enum):
    STRONG_REGULATORY_GRADE = "STRONG_REGULATORY_GRADE"  # All SMD < 0.10, Rubin B < 0.5, large overlap
    ACCEPTABLE_SUPPORTIVE = "ACCEPTABLE_SUPPORTIVE"       # SMD < 0.20, moderate overlap
    HIGH_CONFOUNDING_RISK = "HIGH_CONFOUNDING_RISK"       # Residual imbalances SMD >= 0.20


@dataclass
class SubjectRecord:
    subject_id: str
    is_treated: bool  # True = Single-Arm Trial Patient, False = Historical / RWD Registry Patient
    covariates: Dict[str, float]  # e.g., {'age': 62.0, 'ecog': 1.0, 'prior_lines': 2.0, 'ldh': 240.0}
    time_to_event_months: float  # Follow-up time to OS event or censoring
    event_observed: bool  # True = Event (Death/Progression), False = Censored
    response_achieved: bool = False  # Objective Response (ORR: Complete / Partial Response)
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
    is_balanced: bool  # True if post_smd < 0.10


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
    hazard_ratio: float
    hazard_ratio_ci_low: float
    hazard_ratio_ci_high: float
    log_rank_test_statistic: float
    log_rank_p_value: float
    median_survival_treated_months: float
    median_survival_synthetic_control_months: float
    orr_treated_pct: float
    orr_synthetic_control_pct: float
    att_orr_diff_pct: float
    regulatory_adequacy_tier: RegulatoryAdequacyTier
    recommendations: List[str]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "trial_arm_size": self.trial_arm_size,
            "rwd_pool_size": self.rwd_pool_size,
            "matched_pairs_count": self.matched_pairs_count,
            "matching_method": self.matching_method.value,
            "caliper_width": round(self.caliper_width, 4),
            "mean_absolute_smd_pre": round(self.mean_absolute_smd_pre, 3),
            "mean_absolute_smd_post": round(self.mean_absolute_smd_post, 3),
            "all_covariates_balanced": self.all_covariates_balanced,
            "efficacy_outcomes": {
                "hazard_ratio": round(self.hazard_ratio, 3),
                "hazard_ratio_95_ci": [round(self.hazard_ratio_ci_low, 3), round(self.hazard_ratio_ci_high, 3)],
                "log_rank_p_value": round(self.log_rank_p_value, 4),
                "median_os_treated_months": round(self.median_survival_treated_months, 1),
                "median_os_synthetic_control_months": round(self.median_survival_synthetic_control_months, 1),
                "orr_treated_pct": round(self.orr_treated_pct, 1),
                "orr_synthetic_control_pct": round(self.orr_synthetic_control_pct, 1),
                "att_orr_difference_pct": round(self.att_orr_diff_pct, 1),
            },
            "regulatory_adequacy_tier": self.regulatory_adequacy_tier.value,
            "covariate_balance": [
                {
                    "covariate": cb.covariate_name,
                    "pre_smd": round(cb.pre_smd, 3),
                    "post_smd": round(cb.post_smd, 3),
                    "pre_treated_mean": round(cb.pre_treated_mean, 2),
                    "pre_control_mean": round(cb.pre_control_mean, 2),
                    "post_treated_mean": round(cb.post_treated_mean, 2),
                    "post_control_mean": round(cb.post_control_mean, 2),
                    "is_balanced": cb.is_balanced,
                }
                for cb in self.covariate_balance
            ],
            "recommendations": self.recommendations,
        }
