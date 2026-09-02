"""
Core Biostatistical Engine for Synthetic Control Arm (SCA) Generation
Domain: Real-World Evidence (RWE), Propensity Score Matching & Survival Analysis
Standards: FDA / EMA RWD Guidance, ISPOR-ISPE Good Research Practices
"""

import math
from typing import Dict, List, Optional, Any, Tuple
from .models import (
    MatchingMethod,
    RegulatoryAdequacyTier,
    SubjectRecord,
    MatchedPair,
    CovariateBalance,
    SurvivalCurvePoint,
    SyntheticControlAnalysisResult,
)


class BiostatisticalMath:
    """Standard statistical math routines for RWE analytics."""

    @staticmethod
    def mean(values: List[float]) -> float:
        if not values:
            return 0.0
        return sum(values) / len(values)

    @staticmethod
    def variance(values: List[float], sample: bool = True) -> float:
        n = len(values)
        if n < 2:
            return 0.0
        m = BiostatisticalMath.mean(values)
        ss = sum((x - m) ** 2 for x in values)
        return ss / (n - 1 if sample else n)

    @staticmethod
    def std_dev(values: List[float], sample: bool = True) -> float:
        return math.sqrt(BiostatisticalMath.variance(values, sample))

    @staticmethod
    def standardized_mean_difference(treated: List[float], control: List[float]) -> float:
        """
        Calculate Standardized Mean Difference (SMD / Cohen's d):
        SMD = (Mean_T - Mean_C) / sqrt((Var_T + Var_C) / 2)
        """
        if not treated or not control:
            return 0.0
        m_t = BiostatisticalMath.mean(treated)
        m_c = BiostatisticalMath.mean(control)
        v_t = BiostatisticalMath.variance(treated)
        v_c = BiostatisticalMath.variance(control)
        pooled_sd = math.sqrt((v_t + v_c) / 2.0)
        if pooled_sd < 1e-9:
            return 0.0
        return abs(m_t - m_c) / pooled_sd

    @staticmethod
    def sigmoid(z: float) -> float:
        if z > 35.0:
            return 1.0
        if z < -35.0:
            return 0.0
        return 1.0 / (1.0 + math.exp(-z))

    @staticmethod
    def chi2_sf_1df(x: float) -> float:
        """Survival function (p-value) for Chi-Square distribution with 1 degree of freedom."""
        if x <= 0:
            return 1.0
        # P(Chi2(1) > x) = 2 * (1 - Phi(sqrt(x))) = erfc(sqrt(x/2))
        return math.erfc(math.sqrt(x / 2.0))


class PropensityScoreEngine:
    """
    Multivariate Propensity Score estimation via regularized logistic regression.
    """

    @classmethod
    def fit_and_predict_propensity_scores(
        cls, subjects: List[SubjectRecord], covariate_names: List[str], iterations: int = 150, lr: float = 0.05
    ) -> List[SubjectRecord]:
        """
        Fits logistic regression model P(T=1 | X) using batch gradient descent
        and populates `propensity_score` on each subject.
        """
        if not subjects or not covariate_names:
            for s in subjects:
                s.propensity_score = 0.5
            return subjects

        # Standardize covariates for stable gradient descent
        mins: Dict[str, float] = {}
        maxs: Dict[str, float] = {}
        for c in covariate_names:
            vals = [s.covariates.get(c, 0.0) for s in subjects]
            min_v = min(vals)
            max_v = max(vals)
            spread = max_v - min_v if max_v != min_v else 1.0
            mins[c] = min_v
            maxs[c] = spread

        # Normalize features
        X_norm = []
        Y = []
        for s in subjects:
            row = [1.0]  # Intercept
            for c in covariate_names:
                v = s.covariates.get(c, 0.0)
                norm_v = (v - mins[c]) / maxs[c]
                row.append(norm_v)
            X_norm.append(row)
            Y.append(1.0 if s.is_treated else 0.0)

        n_features = len(covariate_names) + 1
        n_samples = len(subjects)
        weights = [0.0] * n_features

        # Initialize weights with prior log-odds
        n_treated = sum(Y)
        p_t = max(0.01, min(0.99, n_treated / n_samples))
        weights[0] = math.log(p_t / (1.0 - p_t))

        # Gradient Descent
        for _ in range(iterations):
            gradients = [0.0] * n_features
            for i in range(n_samples):
                z = sum(weights[j] * X_norm[i][j] for j in range(n_features))
                p = BiostatisticalMath.sigmoid(z)
                err = p - Y[i]
                for j in range(n_features):
                    gradients[j] += err * X_norm[i][j]

            # Update weights with L2 regularization
            for j in range(n_features):
                reg = 0.01 * weights[j] if j > 0 else 0.0
                weights[j] -= lr * (gradients[j] / n_samples + reg)

        # Predict scores
        for i, s in enumerate(subjects):
            z = sum(weights[j] * X_norm[i][j] for j in range(n_features))
            ps = BiostatisticalMath.sigmoid(z)
            # Bound propensity score away from 0 and 1
            s.propensity_score = max(0.001, min(0.999, ps))

        return subjects


class SurvivalAnalysisEngine:
    """
    Kaplan-Meier survival estimation, Log-Rank testing, and Hazard Ratio modeling.
    """

    @classmethod
    def calculate_kaplan_meier(cls, subjects: List[SubjectRecord]) -> Tuple[List[SurvivalCurvePoint], float]:
        """
        Computes Kaplan-Meier survival curve and median survival time (months).
        """
        if not subjects:
            return [], 0.0

        # Sort by time ascending
        sorted_sub = sorted(subjects, key=lambda s: s.time_to_event_months)
        unique_times = sorted(list(set(s.time_to_event_months for s in sorted_sub)))

        curve = [SurvivalCurvePoint(0.0, len(sorted_sub), 0, 0, 1.0)]
        current_surv = 1.0
        n_at_risk = len(sorted_sub)
        median_time = sorted_sub[-1].time_to_event_months

        for t in unique_times:
            events = sum(1 for s in sorted_sub if s.time_to_event_months == t and s.event_observed)
            censored = sum(1 for s in sorted_sub if s.time_to_event_months == t and not s.event_observed)

            if n_at_risk > 0 and events > 0:
                step = 1.0 - (events / n_at_risk)
                current_surv *= step

            curve.append(SurvivalCurvePoint(t, n_at_risk, events, censored, current_surv))

            if current_surv <= 0.5 and median_time == sorted_sub[-1].time_to_event_months:
                median_time = t

            n_at_risk -= (events + censored)

        return curve, median_time

    @classmethod
    def log_rank_test_and_hazard_ratio(
        cls, treated: List[SubjectRecord], control: List[SubjectRecord]
    ) -> Tuple[float, float, float, float, float]:
        """
        Performs Mantel-Haenszel Log-Rank test and calculates Hazard Ratio (HR) with 95% CI.
        Returns (HR, CI_low, CI_high, chi2_statistic, p_value).
        """
        all_subjects = treated + control
        if not treated or not control or not any(s.event_observed for s in all_subjects):
            return 1.0, 0.5, 2.0, 0.0, 1.0

        all_times = sorted(list(set(s.time_to_event_months for s in all_subjects if s.event_observed)))

        O_T = 0.0  # Observed events treated
        E_T = 0.0  # Expected events treated
        Var_T = 0.0

        O_C = 0.0  # Observed events control
        E_C = 0.0  # Expected events control

        for t in all_times:
            # Number at risk right before time t
            n_t = sum(1 for s in treated if s.time_to_event_months >= t)
            n_c = sum(1 for s in control if s.time_to_event_months >= t)
            n_total = n_t + n_c

            if n_total <= 1:
                continue

            # Number of events at time t
            d_t = sum(1 for s in treated if s.time_to_event_months == t and s.event_observed)
            d_c = sum(1 for s in control if s.time_to_event_months == t and s.event_observed)
            d_total = d_t + d_c

            if d_total == 0:
                continue

            O_T += d_t
            O_C += d_c

            e_t = (n_t / n_total) * d_total
            e_c = (n_c / n_total) * d_total
            E_T += e_t
            E_C += e_c

            var_t = (n_t * n_c * d_total * (n_total - d_total)) / ((n_total ** 2) * (n_total - 1))
            Var_T += var_t

        if Var_T > 1e-9:
            chi2 = ((O_T - E_T) ** 2) / Var_T
            p_val = BiostatisticalMath.chi2_sf_1df(chi2)
        else:
            chi2 = 0.0
            p_val = 1.0

        # Hazard ratio approximation: HR = (O_T / E_T) / (O_C / E_C)
        if E_T > 0 and E_C > 0 and O_C > 0 and O_T > 0:
            hr = (O_T / E_T) / (O_C / E_C)
            se_ln_hr = math.sqrt((1.0 / E_T) + (1.0 / E_C))
            ci_low = math.exp(math.log(hr) - 1.96 * se_ln_hr)
            ci_high = math.exp(math.log(hr) + 1.96 * se_ln_hr)
        else:
            hr = 1.0
            ci_low = 0.5
            ci_high = 2.0

        return hr, ci_low, ci_high, chi2, p_val


class SyntheticControlAgentEngine:
    """
    Main Orchestrator for Propensity Score Matching and Synthetic Control Arm Analysis.
    """

    @classmethod
    def generate_synthetic_control_arm(
        cls,
        subjects: List[SubjectRecord],
        covariate_names: Optional[List[str]] = None,
        caliper_sd_multiplier: float = 0.2,
        matching_method: MatchingMethod = MatchingMethod.NEAREST_NEIGHBOR_CALIPER,
    ) -> SyntheticControlAnalysisResult:
        """
        Executes end-to-end SCA creation:
        1. Propensity score modeling
        2. Caliper-based 1:1 Nearest-Neighbor matching
        3. Pre/Post SMD balance verification
        4. Kaplan-Meier and Cox/Log-Rank comparative survival analytics
        5. Objective Response Rate (ORR) comparative analytics
        6. Regulatory adequacy determination
        """
        treated_pool = [s for s in subjects if s.is_treated]
        control_pool = [s for s in subjects if not s.is_treated]

        if not treated_pool or not control_pool:
            raise ValueError("Input population must contain both treated (trial) and control (RWD) subjects.")

        if covariate_names is None:
            # Discover common covariate names
            covariate_names = sorted(list(set(k for s in subjects for k in s.covariates.keys())))

        # 1. Fit Propensity Scores
        PropensityScoreEngine.fit_and_predict_propensity_scores(subjects, covariate_names)

        # 2. Calculate caliper width based on logit of PS
        ps_logits = [math.log(s.propensity_score / (1.0 - s.propensity_score)) for s in subjects if s.propensity_score is not None]
        sd_logit_ps = BiostatisticalMath.std_dev(ps_logits)
        caliper_width = caliper_sd_multiplier * sd_logit_ps if sd_logit_ps > 0 else 0.1

        # 3. Perform 1:1 Nearest Neighbor Matching on PS
        matched_pairs: List[MatchedPair] = []
        available_control = list(control_pool)

        # Pre-match covariate statistics
        balance_results: List[CovariateBalance] = []
        pre_smd_list = []
        post_smd_list = []

        # Sort treated subjects
        for t in sorted(treated_pool, key=lambda s: s.propensity_score or 0.5):
            t_ps = t.propensity_score or 0.5
            best_c = None
            min_dist = float("inf")

            for c in available_control:
                c_ps = c.propensity_score or 0.5
                dist = abs(t_ps - c_ps)
                if dist <= caliper_width and dist < min_dist:
                    min_dist = dist
                    best_c = c

            if best_c is not None:
                matched_pairs.append(
                    MatchedPair(
                        treated_id=t.subject_id,
                        control_id=best_c.subject_id,
                        ps_distance=min_dist,
                        treated_ps=t_ps,
                        control_ps=best_c.propensity_score or 0.5,
                    )
                )
                available_control.remove(best_c)

        matched_treated = [s for s in treated_pool if any(p.treated_id == s.subject_id for p in matched_pairs)]
        matched_control = [s for s in control_pool if any(p.control_id == s.subject_id for p in matched_pairs)]

        # If matching yielded 0 pairs, fallback to top available
        if not matched_control:
            matched_control = control_pool[:len(treated_pool)]
            matched_treated = treated_pool[:len(matched_control)]

        # 4. Covariate Balance Diagnostics
        for cov in covariate_names:
            pre_t_vals = [s.covariates.get(cov, 0.0) for s in treated_pool]
            pre_c_vals = [s.covariates.get(cov, 0.0) for s in control_pool]
            post_t_vals = [s.covariates.get(cov, 0.0) for s in matched_treated]
            post_c_vals = [s.covariates.get(cov, 0.0) for s in matched_control]

            pre_smd = BiostatisticalMath.standardized_mean_difference(pre_t_vals, pre_c_vals)
            post_smd = BiostatisticalMath.standardized_mean_difference(post_t_vals, post_c_vals)

            pre_smd_list.append(pre_smd)
            post_smd_list.append(post_smd)

            balance_results.append(
                CovariateBalance(
                    covariate_name=cov,
                    pre_treated_mean=BiostatisticalMath.mean(pre_t_vals),
                    pre_control_mean=BiostatisticalMath.mean(pre_c_vals),
                    pre_smd=pre_smd,
                    post_treated_mean=BiostatisticalMath.mean(post_t_vals),
                    post_control_mean=BiostatisticalMath.mean(post_c_vals),
                    post_smd=post_smd,
                    is_balanced=(post_smd < 0.10),
                )
            )

        all_balanced = all(cb.is_balanced for cb in balance_results)
        mean_pre_smd = BiostatisticalMath.mean(pre_smd_list)
        mean_post_smd = BiostatisticalMath.mean(post_smd_list)

        # 5. Comparative Survival Outcomes
        _, med_os_t = SurvivalAnalysisEngine.calculate_kaplan_meier(matched_treated)
        _, med_os_c = SurvivalAnalysisEngine.calculate_kaplan_meier(matched_control)
        hr, ci_l, ci_h, chi2, p_val = SurvivalAnalysisEngine.log_rank_test_and_hazard_ratio(
            matched_treated, matched_control
        )

        # 6. Binary Efficacy (ORR)
        orr_t = (sum(1 for s in matched_treated if s.response_achieved) / len(matched_treated) * 100.0) if matched_treated else 0.0
        orr_c = (sum(1 for s in matched_control if s.response_achieved) / len(matched_control) * 100.0) if matched_control else 0.0
        att_orr = orr_t - orr_c

        # 7. Regulatory Adequacy Tier
        if all_balanced and mean_post_smd <= 0.08 and len(matched_pairs) >= int(0.75 * len(treated_pool)):
            reg_tier = RegulatoryAdequacyTier.STRONG_REGULATORY_GRADE
        elif mean_post_smd <= 0.15:
            reg_tier = RegulatoryAdequacyTier.ACCEPTABLE_SUPPORTIVE
        else:
            reg_tier = RegulatoryAdequacyTier.HIGH_CONFOUNDING_RISK

        recs = cls._generate_recommendations(reg_tier, all_balanced, len(matched_pairs), len(treated_pool), hr, p_val)

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
            hazard_ratio_ci_low=ci_l,
            hazard_ratio_ci_high=ci_h,
            log_rank_test_statistic=chi2,
            log_rank_p_value=p_val,
            median_survival_treated_months=med_os_t,
            median_survival_synthetic_control_months=med_os_c,
            orr_treated_pct=orr_t,
            orr_synthetic_control_pct=orr_c,
            att_orr_diff_pct=att_orr,
            regulatory_adequacy_tier=reg_tier,
            recommendations=recs,
        )

    @classmethod
    def _generate_recommendations(
        cls, tier: RegulatoryAdequacyTier, all_balanced: bool, matched_n: int, treated_n: int, hr: float, p_val: float
    ) -> List[str]:
        recs = []
        retention_pct = (matched_n / treated_n * 100.0) if treated_n > 0 else 0.0

        if tier == RegulatoryAdequacyTier.STRONG_REGULATORY_GRADE:
            recs.append("Covariate balance meets stringent FDA/EMA RWE benchmark (SMD < 0.10 across all key confounders).")
        elif tier == RegulatoryAdequacyTier.ACCEPTABLE_SUPPORTIVE:
            recs.append("Acceptable supportive evidence; consider doubly robust estimation (IPTW + outcome regression) to address mild residual variance.")
        else:
            recs.append("High risk of residual unmeasured confounding. Expand RWD cohort or calibrate wider inclusion criteria.")

        if retention_pct < 80.0:
            recs.append(f"Matching caliper excluded {100.0 - retention_pct:.1f}% of trial cohort. Evaluate caliper expansion or full optimal matching.")

        if hr < 0.70 and p_val < 0.05:
            recs.append(f"Statistically significant OS benefit observed (HR={hr:.2f}, p={p_val:.4f}) versus Synthetic Control Arm.")
        else:
            recs.append(f"Treatment effect estimate (HR={hr:.2f}) does not demonstrate definitive statistical superiority against historical controls.")

        return recs


def parse_synthetic_cohort_dict(data: Dict[str, Any]) -> List[SubjectRecord]:
    """Parse JSON/dictionary into List[SubjectRecord]."""
    subjects = []
    for raw in data.get("subjects", []):
        subjects.append(
            SubjectRecord(
                subject_id=str(raw.get("subject_id", "SUBJ-UNKNOWN")),
                is_treated=bool(raw.get("is_treated", False)),
                covariates={str(k): float(v) for k, v in raw.get("covariates", {}).items()},
                time_to_event_months=float(raw.get("time_to_event_months", 12.0)),
                event_observed=bool(raw.get("event_observed", True)),
                response_achieved=bool(raw.get("response_achieved", False)),
            )
        )
    return subjects
