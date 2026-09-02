#!/usr/bin/env python3
"""
Command-Line Interface for Clinical Trial Synthetic Control Agent
=================================================================
Performs Propensity Score Matching (PSM), Covariate Balance Diagnostics (SMD),
Kaplan-Meier Survival Analysis, and Regulatory Adequacy Assessment for Single-Arm Trials.
"""

import sys
import os
import json
import argparse
import random
from typing import Dict, List, Any

# Ensure project path is accessible
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from synthetic_control_arm import (
    SyntheticControlAgentEngine,
    SubjectRecord,
    SyntheticControlAnalysisResult,
    parse_synthetic_cohort_dict,
)


def generate_benchmark_synthetic_cohort(n_treated: int = 40, n_control: int = 150) -> List[SubjectRecord]:
    """Generates a synthetic oncology cohort (Trial Arm vs RWD Registry)."""
    random.seed(42)
    subjects: List[SubjectRecord] = []

    # Treated group (Single-Arm Trial: younger, slightly better ECOG, high response rate)
    for i in range(n_treated):
        age = random.gauss(58.0, 8.0)
        ecog = random.choices([0, 1, 2], weights=[0.6, 0.35, 0.05])[0]
        prior_lines = random.choices([1, 2, 3], weights=[0.5, 0.4, 0.1])[0]
        ldh = random.gauss(210.0, 45.0)

        # Survival time in months (better for treated)
        time_os = random.expovariate(1.0 / 22.0) + 4.0
        event = random.random() < 0.65
        orr = random.random() < 0.55

        subjects.append(
            SubjectRecord(
                subject_id=f"TRIAL-PT-{i+1:03d}",
                is_treated=True,
                covariates={"age": round(age, 1), "ecog": float(ecog), "prior_lines": float(prior_lines), "ldh": round(ldh, 1)},
                time_to_event_months=round(time_os, 1),
                event_observed=event,
                response_achieved=orr,
            )
        )

    # Control pool (RWD Registry: older, more prior lines, higher LDH, worse survival)
    for j in range(n_control):
        age = random.gauss(64.0, 10.0)
        ecog = random.choices([0, 1, 2], weights=[0.3, 0.5, 0.2])[0]
        prior_lines = random.choices([1, 2, 3, 4], weights=[0.2, 0.4, 0.3, 0.1])[0]
        ldh = random.gauss(260.0, 60.0)

        time_os = random.expovariate(1.0 / 13.0) + 2.0
        event = random.random() < 0.85
        orr = random.random() < 0.25

        subjects.append(
            SubjectRecord(
                subject_id=f"RWD-REG-{j+1:03d}",
                is_treated=False,
                covariates={"age": round(age, 1), "ecog": float(ecog), "prior_lines": float(prior_lines), "ldh": round(ldh, 1)},
                time_to_event_months=round(time_os, 1),
                event_observed=event,
                response_achieved=orr,
            )
        )

    return subjects


def format_sca_report(res: SyntheticControlAnalysisResult) -> str:
    lines = []
    lines.append("=" * 82)
    lines.append(" SYNTHETIC CONTROL ARM (SCA) & PROPENSITY SCORE MATCHING REPORT")
    lines.append(" Regulatory Benchmark: FDA / EMA Real-World Evidence Guidance")
    lines.append("=" * 82)
    lines.append(f"Trial Arm Sample (Treated)   : {res.trial_arm_size} patients")
    lines.append(f"RWD Control Pool Size        : {res.rwd_pool_size} patients")
    lines.append(f"Matched Pairs Generated (1:1): {res.matched_pairs_count} pairs (Caliper = {res.caliper_width:.4f} logit SD)")
    lines.append(f"Regulatory Adequacy Tier     : [{res.regulatory_adequacy_tier.value}]")
    lines.append("-" * 82)
    lines.append("COVARIATE BALANCE DIAGNOSTICS (Standardized Mean Difference - SMD):")
    lines.append(f"{'Covariate':<15} | {'Pre-T Mean':<10} | {'Pre-C Mean':<10} | {'Pre-SMD':<8} | {'Post-T Mean':<11} | {'Post-C Mean':<11} | {'Post-SMD':<8} | {'Balance'}")
    lines.append("-" * 82)
    for cb in res.covariate_balance:
        bal_mark = "[PASS <0.10]" if cb.is_balanced else "[FAIL >=0.10]"
        lines.append(f"{cb.covariate_name:<15} | {cb.pre_treated_mean:<10.1f} | {cb.pre_control_mean:<10.1f} | {cb.pre_smd:<8.3f} | {cb.post_treated_mean:<11.1f} | {cb.post_control_mean:<11.1f} | {cb.post_smd:<8.3f} | {bal_mark}")
    lines.append("-" * 82)
    lines.append(f"Mean Absolute SMD: Pre-Matching = {res.mean_absolute_smd_pre:.3f}  ==>  Post-Matching = {res.mean_absolute_smd_post:.3f}")
    lines.append(f"All Covariates Balanced (SMD < 0.10): {'YES (Optimal)' if res.all_covariates_balanced else 'NO (Residual Variance)'}")
    lines.append("-" * 82)
    lines.append("COMPARATIVE CLINICAL EFFICACY & SURVIVAL OUTCOMES:")
    lines.append(f"  * Hazard Ratio (HR, OS)    : {res.hazard_ratio:.3f} (95% CI: {res.hazard_ratio_ci_low:.3f} - {res.hazard_ratio_ci_high:.3f})")
    lines.append(f"  * Log-Rank Test p-value    : {res.log_rank_p_value:.4f} (Chi2 = {res.log_rank_test_statistic:.2f})")
    lines.append(f"  * Median Overall Survival  : Treated = {res.median_survival_treated_months:.1f} mos vs Synthetic Control = {res.median_survival_synthetic_control_months:.1f} mos")
    lines.append(f"  * Objective Response (ORR) : Treated = {res.orr_treated_pct:.1f}% vs Synthetic Control = {res.orr_synthetic_control_pct:.1f}% (ATT Diff: +{res.att_orr_diff_pct:.1f}%)")
    lines.append("-" * 82)
    lines.append("REGULATORY & METHODOLOGICAL RECOMMENDATIONS:")
    for idx, rec in enumerate(res.recommendations, start=1):
        lines.append(f"  {idx}. {rec}")
    lines.append("=" * 82)
    return "\n".join(lines)


def interactive_mode():
    print("\n--- Interactive Synthetic Control Arm Builder ---")
    n_treated = int(input("Number of Single-Arm Trial Patients [e.g. 30]: ").strip() or "30")
    n_control = int(input("Number of RWD Registry Patients [e.g. 100]: ").strip() or "100")
    caliper = float(input("Caliper Multiplier (SD of logit PS) [default 0.2]: ").strip() or "0.2")

    print("\nGenerating realistic baseline cohort...")
    subjects = generate_benchmark_synthetic_cohort(n_treated, n_control)
    res = SyntheticControlAgentEngine.generate_synthetic_control_arm(subjects, caliper_sd_multiplier=caliper)
    print("\n" + format_sca_report(res))


def main():
    parser = argparse.ArgumentParser(
        description="Clinical Trial Synthetic Control Agent - Propensity Score Matching & Survival Analysis"
    )
    parser.add_argument("--demo", action="store_true", help="Run benchmark matching analysis with synthetic oncology cohort")
    parser.add_argument("--file", type=str, help="Path to JSON file containing patient cohort records")
    parser.add_argument("--caliper", type=float, default=0.2, help="Caliper width multiplier (default: 0.2 SD of logit PS)")
    parser.add_argument("--json", action="store_true", help="Output analysis results in JSON format")
    parser.add_argument("--interactive", action="store_true", help="Run interactive cohort matching wizard")

    args = parser.parse_args()

    if args.interactive:
        interactive_mode()
        return

    if args.file:
        with open(args.file, "r", encoding="utf-8") as f:
            data = json.load(f)
        subjects = parse_synthetic_cohort_dict(data)
    else:
        subjects = generate_benchmark_synthetic_cohort()

    res = SyntheticControlAgentEngine.generate_synthetic_control_arm(
        subjects, caliper_sd_multiplier=args.caliper
    )

    if args.json:
        print(json.dumps(res.to_dict(), indent=2))
    else:
        print(format_sca_report(res))


if __name__ == "__main__":
    main()
