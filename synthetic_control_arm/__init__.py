"""
Synthetic Control Arm Package - Real-World Evidence & Propensity Matching Engine
"""

from .models import (
    MatchingMethod,
    RegulatoryAdequacyTier,
    SubjectRecord,
    MatchedPair,
    CovariateBalance,
    SurvivalCurvePoint,
    SyntheticControlAnalysisResult,
)
from .engine import (
    BiostatisticalMath,
    PropensityScoreEngine,
    SurvivalAnalysisEngine,
    SyntheticControlAgentEngine,
    parse_synthetic_cohort_dict,
)

__all__ = [
    "MatchingMethod",
    "RegulatoryAdequacyTier",
    "SubjectRecord",
    "MatchedPair",
    "CovariateBalance",
    "SurvivalCurvePoint",
    "SyntheticControlAnalysisResult",
    "BiostatisticalMath",
    "PropensityScoreEngine",
    "SurvivalAnalysisEngine",
    "SyntheticControlAgentEngine",
    "parse_synthetic_cohort_dict",
]
