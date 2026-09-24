"""Core propensity-score matching and matched-cohort outcome diagnostics."""

import math
from typing import Any, Dict, List, Optional, Tuple

from .models import (
    CovariateBalance,
    MatchedPair,
    MatchingMethod,
    RegulatoryAdequacyTier,
    SubjectRecord,
    SurvivalCurvePoint,
    SyntheticControlAnalysisResult,
)


class BiostatisticalMath:
    """Small dependency-free statistical helpers used by the engine."""

    @staticmethod
    def mean(values: List[float]) -> float:
        return 0.0 if not values else sum(values) / len(values)

    @staticmethod
    def variance(values: List[float], sample: bool = True) -> float:
        n = len(values)
        if n < 2:
            return 0.0
        mean_value = BiostatisticalMath.mean(values)
        denominator = n - 1 if sample else n
        return sum((value - mean_value) ** 2 for value in values) / denominator

    @staticmethod
    def std_dev(values: List[float], sample: bool = True) -> float:
        return math.sqrt(BiostatisticalMath.variance(values, sample))

    @staticmethod
    def standardized_mean_difference(
        treated: List[float], control: List[float]
    ) -> float:
        """Return the absolute SMD using the pooled within-group SD.

        If both groups are constant but have different means, the standardized
        difference is unbounded rather than zero. Empty groups are likewise not
        interpretable and return infinity.
        """
        if not treated or not control:
            return math.inf
        mean_t = BiostatisticalMath.mean(treated)
        mean_c = BiostatisticalMath.mean(control)
        var_t = BiostatisticalMath.variance(treated)
        var_c = BiostatisticalMath.variance(control)
        pooled_sd = math.sqrt((var_t + var_c) / 2.0)
        if pooled_sd < 1e-12:
            return 0.0 if math.isclose(mean_t, mean_c, abs_tol=1e-12) else math.inf
        return abs(mean_t - mean_c) / pooled_sd

    @staticmethod
    def sigmoid(z: float) -> float:
        if z > 35.0:
            return 1.0
        if z < -35.0:
            return 0.0
        return 1.0 / (1.0 + math.exp(-z))

    @staticmethod
    def logit(probability: float) -> float:
        p = max(1e-12, min(1.0 - 1e-12, probability))
        return math.log(p / (1.0 - p))

    @staticmethod
    def chi2_sf_1df(value: float) -> float:
        if value <= 0:
            return 1.0
        return math.erfc(math.sqrt(value / 2.0))


class PropensityScoreEngine:
    """Regularized logistic propensity-score model implemented without dependencies."""

    @classmethod
    def fit_and_predict_propensity_scores(
        cls,
        subjects: List[SubjectRecord],
        covariate_names: List[str],
        iterations: int = 300,
        lr: float = 0.05,
    ) -> List[SubjectRecord]:
        if not subjects or not covariate_names:
            for subject in subjects:
                subject.propensity_score = 0.5
            return subjects

        means: Dict[str, float] = {}
        scales: Dict[str, float] = {}
        for covariate in covariate_names:
            values = [subject.covariates[covariate] for subject in subjects]
            means[covariate] = BiostatisticalMath.mean(values)
            scale = BiostatisticalMath.std_dev(values)
            scales[covariate] = scale if scale > 1e-12 else 1.0

        design_matrix: List[List[float]] = []
        outcomes: List[float] = []
        for subject in subjects:
            row = [1.0]
            row.extend(
                (subject.covariates[covariate] - means[covariate])
                / scales[covariate]
                for covariate in covariate_names
            )
            design_matrix.append(row)
            outcomes.append(1.0 if subject.is_treated else 0.0)

        n_samples = len(subjects)
        n_features = len(covariate_names) + 1
        weights = [0.0] * n_features
        treated_fraction = max(0.01, min(0.99, sum(outcomes) / n_samples))
        weights[0] = math.log(treated_fraction / (1.0 - treated_fraction))

        for _ in range(iterations):
            gradients = [0.0] * n_features
            for i, row in enumerate(design_matrix):
                probability = BiostatisticalMath.sigmoid(
                    sum(weights[j] * row[j] for j in range(n_features))
                )
                error = probability - outcomes[i]
                for j in range(n_features):
                    gradients[j] += error * row[j]

            for j in range(n_features):
                l2_penalty = 0.01 * weights[j] if j > 0 else 0.0
                weights[j] -= lr * (gradients[j] / n_samples + l2_penalty)

        for subject, row in zip(subjects, design_matrix):
            score = BiostatisticalMath.sigmoid(
                sum(weights[j] * row[j] for j in range(n_features))
            )
            subject.propensity_score = max(0.001, min(0.999, score))
        return subjects


class SurvivalAnalysisEngine:
    """Kaplan-Meier estimation and log-rank score-test diagnostics."""

    @classmethod
    def calculate_kaplan_meier(
        cls, subjects: List[SubjectRecord]
    ) -> Tuple[List[SurvivalCurvePoint], Optional[float]]:
        if not subjects:
            return [], None

        ordered = sorted(subjects, key=lambda subject: subject.time_to_event_months)
        times = sorted({subject.time_to_event_months for subject in ordered})
        curve = [SurvivalCurvePoint(0.0, len(ordered), 0, 0, 1.0)]
        survival = 1.0
        n_at_risk = len(ordered)
        median_time: Optional[float] = None

        for time in times:
            events = sum(
                1
                for subject in ordered
                if subject.time_to_event_months == time and subject.event_observed
            )
            censored = sum(
                1
                for subject in ordered
                if subject.time_to_event_months == time and not subject.event_observed
            )
            if n_at_risk > 0 and events:
                survival *= 1.0 - events / n_at_risk
            curve.append(
                SurvivalCurvePoint(time, n_at_risk, events, censored, survival)
            )
            if median_time is None and survival <= 0.5:
                median_time = time
            n_at_risk -= events + censored

        return curve, median_time

    @classmethod
    def log_rank_test_and_hazard_ratio(
        cls, treated: List[SubjectRecord], control: List[SubjectRecord]
    ) -> Tuple[Optional[float], Optional[float], Optional[float], float, float]:
        """Return a log-rank test and score-test approximation to the hazard ratio.

        The hazard ratio is an approximation derived from the log-rank score and
        information, not a fitted Cox proportional-hazards model.
        """
        all_subjects = treated + control
        if not treated or not control or not any(
            subject.event_observed for subject in all_subjects
        ):
            return None, None, None, 0.0, 1.0

        event_times = sorted(
            {
                subject.time_to_event_months
                for subject in all_subjects
                if subject.event_observed
            }
        )
        observed_treated = 0.0
        expected_treated = 0.0
        variance = 0.0

        for time in event_times:
            n_treated = sum(
                1 for subject in treated if subject.time_to_event_months >= time
            )
            n_control = sum(
                1 for subject in control if subject.time_to_event_months >= time
            )
            n_total = n_treated + n_control
            if n_total <= 1:
                continue

            events_treated = sum(
                1
                for subject in treated
                if subject.time_to_event_months == time and subject.event_observed
            )
            events_control = sum(
                1
                for subject in control
                if subject.time_to_event_months == time and subject.event_observed
            )
            events_total = events_treated + events_control
            if not events_total:
                continue

            observed_treated += events_treated
            expected_treated += (n_treated / n_total) * events_total
            variance += (
                n_treated
                * n_control
                * events_total
                * (n_total - events_total)
                / ((n_total**2) * (n_total - 1))
            )

        if variance <= 1e-12:
            return None, None, None, 0.0, 1.0

        score = observed_treated - expected_treated
        chi2 = (score**2) / variance
        p_value = BiostatisticalMath.chi2_sf_1df(chi2)

        log_hr = score / variance
        standard_error = math.sqrt(1.0 / variance)
        hazard_ratio = math.exp(log_hr)
        ci_low = math.exp(log_hr - 1.96 * standard_error)
        ci_high = math.exp(log_hr + 1.96 * standard_error)
        return hazard_ratio, ci_low, ci_high, chi2, p_value


class SyntheticControlAgentEngine:
    """End-to-end propensity-score matching and matched-cohort analysis."""

    @staticmethod
    def _validate_subjects(
        subjects: List[SubjectRecord], covariate_names: List[str]
    ) -> None:
        if not subjects:
            raise ValueError("At least one subject is required.")
        identifiers = [subject.subject_id.strip() for subject in subjects]
        if any(not identifier for identifier in identifiers):
            raise ValueError("subject_id must be non-empty for every subject.")
        if len(set(identifiers)) != len(identifiers):
            raise ValueError("subject_id values must be unique.")
        if not covariate_names:
            raise ValueError("At least one numeric baseline covariate is required.")

        for subject in subjects:
            if not math.isfinite(subject.time_to_event_months) or subject.time_to_event_months < 0:
                raise ValueError(
                    f"Invalid time_to_event_months for subject {subject.subject_id}."
                )
            missing = [name for name in covariate_names if name not in subject.covariates]
            if missing:
                raise ValueError(
                    f"Subject {subject.subject_id} is missing covariates: {', '.join(missing)}."
                )
            for name in covariate_names:
                value = subject.covariates[name]
                if not math.isfinite(value):
                    raise ValueError(
                        f"Covariate {name!r} is not finite for subject {subject.subject_id}."
                    )

    @classmethod
    def generate_synthetic_control_arm(
        cls,
        subjects: List[SubjectRecord],
        covariate_names: Optional[List[str]] = None,
        caliper_sd_multiplier: float = 0.2,
        matching_method: MatchingMethod = MatchingMethod.NEAREST_NEIGHBOR_CALIPER,
    ) -> SyntheticControlAnalysisResult:
        treated_pool = [subject for subject in subjects if subject.is_treated]
        control_pool = [subject for subject in subjects if not subject.is_treated]
        if not treated_pool or not control_pool:
            raise ValueError(
                "Input population must contain both treated and control subjects."
            )
        if matching_method != MatchingMethod.NEAREST_NEIGHBOR_CALIPER:
            raise NotImplementedError(
                f"Matching method {matching_method.value} is not implemented."
            )
        if not math.isfinite(caliper_sd_multiplier) or caliper_sd_multiplier <= 0:
            raise ValueError("caliper_sd_multiplier must be a finite value greater than zero.")

        if covariate_names is None:
            covariate_names = sorted(
                set.intersection(*(set(subject.covariates) for subject in subjects))
            ) if subjects else []
        else:
            covariate_names = list(dict.fromkeys(covariate_names))
        cls._validate_subjects(subjects, covariate_names)

        PropensityScoreEngine.fit_and_predict_propensity_scores(
            subjects, covariate_names
        )

        logits = [
            BiostatisticalMath.logit(subject.propensity_score or 0.5)
            for subject in subjects
        ]
        sd_logit_ps = BiostatisticalMath.std_dev(logits)
        caliper_width = caliper_sd_multiplier * sd_logit_ps
        if caliper_width <= 1e-12:
            caliper_width = 1e-12

        matched_pairs: List[MatchedPair] = []
        available_controls = list(control_pool)
        for treated in sorted(
            treated_pool,
            key=lambda subject: BiostatisticalMath.logit(
                subject.propensity_score or 0.5
            ),
        ):
            treated_ps = treated.propensity_score or 0.5
            treated_logit = BiostatisticalMath.logit(treated_ps)
            best_control: Optional[SubjectRecord] = None
            best_distance = math.inf
            for control in available_controls:
                control_ps = control.propensity_score or 0.5
                distance = abs(
                    treated_logit - BiostatisticalMath.logit(control_ps)
                )
                if distance <= caliper_width and distance < best_distance:
                    best_distance = distance
                    best_control = control
            if best_control is not None:
                matched_pairs.append(
                    MatchedPair(
                        treated_id=treated.subject_id,
                        control_id=best_control.subject_id,
                        ps_distance=best_distance,
                        treated_ps=treated_ps,
                        control_ps=best_control.propensity_score or 0.5,
                    )
                )
                available_controls.remove(best_control)

        if not matched_pairs:
            raise ValueError(
                "No treated-control pairs met the propensity-score caliper. "
                "Review overlap, covariates, or the caliper setting."
            )

        subject_by_id = {subject.subject_id: subject for subject in subjects}
        matched_treated = [subject_by_id[pair.treated_id] for pair in matched_pairs]
        matched_control = [subject_by_id[pair.control_id] for pair in matched_pairs]

        balance_results: List[CovariateBalance] = []
        pre_smd_values: List[float] = []
        post_smd_values: List[float] = []
        for covariate in covariate_names:
            pre_treated = [subject.covariates[covariate] for subject in treated_pool]
            pre_control = [subject.covariates[covariate] for subject in control_pool]
            post_treated = [subject.covariates[covariate] for subject in matched_treated]
            post_control = [subject.covariates[covariate] for subject in matched_control]
            pre_smd = BiostatisticalMath.standardized_mean_difference(
                pre_treated, pre_control
            )
            post_smd = BiostatisticalMath.standardized_mean_difference(
                post_treated, post_control
            )
            pre_smd_values.append(pre_smd)
            post_smd_values.append(post_smd)
            balance_results.append(
                CovariateBalance(
                    covariate_name=covariate,
                    pre_treated_mean=BiostatisticalMath.mean(pre_treated),
                    pre_control_mean=BiostatisticalMath.mean(pre_control),
                    pre_smd=pre_smd,
                    post_treated_mean=BiostatisticalMath.mean(post_treated),
                    post_control_mean=BiostatisticalMath.mean(post_control),
                    post_smd=post_smd,
                    is_balanced=post_smd < 0.10,
                )
            )

        all_balanced = all(item.is_balanced for item in balance_results)
        mean_pre_smd = BiostatisticalMath.mean(pre_smd_values)
        mean_post_smd = BiostatisticalMath.mean(post_smd_values)
        retention_pct = len(matched_pairs) / len(treated_pool) * 100.0

        _, median_treated = SurvivalAnalysisEngine.calculate_kaplan_meier(matched_treated)
        _, median_control = SurvivalAnalysisEngine.calculate_kaplan_meier(matched_control)
        hr, ci_low, ci_high, chi2, p_value = (
            SurvivalAnalysisEngine.log_rank_test_and_hazard_ratio(
                matched_treated, matched_control
            )
        )

        orr_treated = (
            sum(subject.response_achieved for subject in matched_treated)
            / len(matched_treated)
            * 100.0
        )
        orr_control = (
            sum(subject.response_achieved for subject in matched_control)
            / len(matched_control)
            * 100.0
        )

        if all_balanced and mean_post_smd <= 0.08 and retention_pct >= 75.0:
            tier = RegulatoryAdequacyTier.STRONG_REGULATORY_GRADE
        elif mean_post_smd <= 0.15 and retention_pct >= 50.0:
            tier = RegulatoryAdequacyTier.ACCEPTABLE_SUPPORTIVE
        else:
            tier = RegulatoryAdequacyTier.HIGH_CONFOUNDING_RISK

        recommendations = cls._generate_recommendations(
            tier=tier,
            all_balanced=all_balanced,
            retention_pct=retention_pct,
            hazard_ratio=hr,
            p_value=p_value,
        )

        return SyntheticControlAnalysisResult(
            trial_arm_size=len(treated_pool),
            rwd_pool_size=len(control_pool),
            matched_pairs_count=len(matched_pairs),
            matching_method=matching_method,
            caliper_width=caliper_width,
            covariate_balance=balance_results,
            all_covariates_balanced=all_balanced,
            mean_absolute_smd_pre=mean_pre_smd,
            mean_absolute_smd_post=mean_post_smd,
            hazard_ratio=hr,
            hazard_ratio_ci_low=ci_low,
            hazard_ratio_ci_high=ci_high,
            log_rank_test_statistic=chi2,
            log_rank_p_value=p_value,
            median_survival_treated_months=median_treated,
            median_survival_synthetic_control_months=median_control,
            orr_treated_pct=orr_treated,
            orr_synthetic_control_pct=orr_control,
            att_orr_diff_pct=orr_treated - orr_control,
            regulatory_adequacy_tier=tier,
            recommendations=recommendations,
            matched_pairs=matched_pairs,
            matching_retention_pct=retention_pct,
        )

    @staticmethod
    def _generate_recommendations(
        tier: RegulatoryAdequacyTier,
        all_balanced: bool,
        retention_pct: float,
        hazard_ratio: Optional[float],
        p_value: float,
    ) -> List[str]:
        recommendations = [
            "Balance labels are tool-defined diagnostics, not regulatory determinations."
        ]
        if tier == RegulatoryAdequacyTier.STRONG_REGULATORY_GRADE:
            recommendations.append(
                "Matched covariates meet the conventional SMD < 0.10 balance benchmark with good treated-cohort retention."
            )
        elif tier == RegulatoryAdequacyTier.ACCEPTABLE_SUPPORTIVE:
            recommendations.append(
                "Residual imbalance or overlap limitations remain; inspect individual SMDs and consider sensitivity analyses."
            )
        else:
            recommendations.append(
                "Poor balance or limited overlap materially constrains interpretation; revise the cohort or adjustment strategy."
            )
        if retention_pct < 80.0:
            recommendations.append(
                f"Only {retention_pct:.1f}% of treated subjects were retained after caliper matching."
            )
        if not all_balanced:
            recommendations.append(
                "At least one measured covariate remains above the conventional absolute SMD threshold of 0.10."
            )
        if hazard_ratio is not None:
            recommendations.append(
                f"Matched-sample log-rank HR approximation = {hazard_ratio:.2f} (p={p_value:.4f}); this is not a fitted Cox model or proof of causality."
            )
        return recommendations


def _parse_bool(value: Any, field_name: str) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and value in (0, 1):
        return bool(value)
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "1", "t", "yes", "y"}:
            return True
        if normalized in {"false", "0", "f", "no", "n"}:
            return False
    raise ValueError(f"{field_name} must be a boolean value.")


def parse_synthetic_cohort_dict(data: Dict[str, Any]) -> List[SubjectRecord]:
    raw_subjects = data.get("subjects")
    if not isinstance(raw_subjects, list):
        raise ValueError("JSON input must contain a 'subjects' array.")

    subjects: List[SubjectRecord] = []
    for index, raw in enumerate(raw_subjects, start=1):
        if not isinstance(raw, dict):
            raise ValueError(f"Subject {index} must be an object.")
        covariates = raw.get("covariates", {})
        if not isinstance(covariates, dict):
            raise ValueError(f"Subject {index} covariates must be an object.")
        subjects.append(
            SubjectRecord(
                subject_id=str(raw.get("subject_id", "")).strip(),
                is_treated=_parse_bool(raw.get("is_treated"), "is_treated"),
                covariates={str(key): float(value) for key, value in covariates.items()},
                time_to_event_months=float(raw.get("time_to_event_months")),
                event_observed=_parse_bool(
                    raw.get("event_observed"), "event_observed"
                ),
                response_achieved=_parse_bool(
                    raw.get("response_achieved", False), "response_achieved"
                ),
            )
        )
    return subjects
