"""Command-line interface for the synthetic-control matching engine."""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

from .engine import SyntheticControlAgentEngine, parse_synthetic_cohort_dict
from .models import MatchedPair, SubjectRecord, SyntheticControlAnalysisResult

VERSION = "2.1.0"
_TRUE_VALUES = {"true", "1", "t", "yes", "y"}
_FALSE_VALUES = {"false", "0", "f", "no", "n"}
_REQUIRED_CSV_FIELDS = {
    "subject_id",
    "is_treated",
    "time_to_event_months",
    "event_observed",
    "response_achieved",
}
_NON_COVARIATE_FIELDS = _REQUIRED_CSV_FIELDS | {
    "propensity_score",
    "iptw_weight",
    "cohort_arm",
    "match_status",
    "matched_pair_id",
    "logit_ps_distance",
    "balance_assessment",
    "regulatory_tier",
    "mean_post_smd",
    "hazard_ratio",
    "log_rank_p_value",
}


def _parse_bool(value: Any, field_name: str, row_number: Optional[int] = None) -> bool:
    normalized = str(value).strip().lower()
    if normalized in _TRUE_VALUES:
        return True
    if normalized in _FALSE_VALUES:
        return False
    where = f" on row {row_number}" if row_number is not None else ""
    raise ValueError(f"{field_name}{where} must be true/false or 1/0.")


def _parse_finite_float(value: Any, field_name: str, row_number: int) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field_name} on row {row_number} must be numeric.") from exc
    if not math.isfinite(number):
        raise ValueError(f"{field_name} on row {row_number} must be finite.")
    return number


def generate_benchmark_synthetic_cohort(
    n_treated: int = 40, n_control: int = 150, seed: int = 42
) -> List[SubjectRecord]:
    """Generate a deterministic synthetic oncology-like demonstration cohort."""
    if n_treated < 1 or n_control < 1:
        raise ValueError("Demo cohort sizes must both be at least 1.")
    random.seed(seed)
    subjects: List[SubjectRecord] = []

    for i in range(n_treated):
        subjects.append(
            SubjectRecord(
                subject_id=f"TRIAL-PT-{i + 1:03d}",
                is_treated=True,
                covariates={
                    "age": round(random.gauss(58.0, 8.0), 1),
                    "ecog": float(
                        random.choices([0, 1, 2], weights=[0.6, 0.35, 0.05])[0]
                    ),
                    "prior_lines": float(
                        random.choices([1, 2, 3], weights=[0.5, 0.4, 0.1])[0]
                    ),
                    "ldh": round(random.gauss(210.0, 45.0), 1),
                },
                time_to_event_months=round(random.expovariate(1.0 / 22.0) + 4.0, 1),
                event_observed=random.random() < 0.65,
                response_achieved=random.random() < 0.55,
            )
        )

    for i in range(n_control):
        subjects.append(
            SubjectRecord(
                subject_id=f"RWD-REG-{i + 1:03d}",
                is_treated=False,
                covariates={
                    "age": round(random.gauss(64.0, 10.0), 1),
                    "ecog": float(
                        random.choices([0, 1, 2], weights=[0.3, 0.5, 0.2])[0]
                    ),
                    "prior_lines": float(
                        random.choices([1, 2, 3, 4], weights=[0.2, 0.4, 0.3, 0.1])[0]
                    ),
                    "ldh": round(random.gauss(260.0, 60.0), 1),
                },
                time_to_event_months=round(random.expovariate(1.0 / 13.0) + 2.0, 1),
                event_observed=random.random() < 0.85,
                response_achieved=random.random() < 0.25,
            )
        )
    return subjects


def load_cohort_from_csv(csv_path: str | Path) -> List[SubjectRecord]:
    """Load and strictly validate a patient-level cohort CSV."""
    path = Path(csv_path)
    if not path.is_file():
        raise ValueError(f"Input CSV does not exist: {path}")

    subjects: List[SubjectRecord] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        raw_fieldnames = reader.fieldnames or []
        if any(name is None or not name.strip() for name in raw_fieldnames):
            raise ValueError("CSV column names must be non-empty.")
        fieldnames = [name.strip() for name in raw_fieldnames]
        if len(set(fieldnames)) != len(fieldnames):
            raise ValueError(
                "CSV column names must be unique after trimming whitespace."
            )
        missing = sorted(_REQUIRED_CSV_FIELDS - set(fieldnames))
        if missing:
            raise ValueError(
                "CSV is missing required columns: " + ", ".join(missing)
            )
        covariate_fields = [
            name for name in fieldnames if name not in _NON_COVARIATE_FIELDS
        ]
        if not covariate_fields:
            raise ValueError("CSV must contain at least one numeric baseline covariate.")

        for row_number, raw_row in enumerate(reader, start=2):
            row = {str(key).strip(): value for key, value in raw_row.items() if key}
            subject_id = str(row.get("subject_id", "")).strip()
            if not subject_id:
                raise ValueError(f"subject_id is required on row {row_number}.")
            time_value = _parse_finite_float(
                row.get("time_to_event_months"), "time_to_event_months", row_number
            )
            if time_value < 0:
                raise ValueError(
                    f"time_to_event_months on row {row_number} must be >= 0."
                )
            covariates: Dict[str, float] = {}
            for field in covariate_fields:
                raw_value = row.get(field, "")
                if raw_value is None or str(raw_value).strip() == "":
                    raise ValueError(
                        f"Covariate {field!r} is missing on row {row_number}."
                    )
                covariates[field] = _parse_finite_float(
                    raw_value, field, row_number
                )

            subjects.append(
                SubjectRecord(
                    subject_id=subject_id,
                    is_treated=_parse_bool(
                        row.get("is_treated"), "is_treated", row_number
                    ),
                    covariates=covariates,
                    time_to_event_months=time_value,
                    event_observed=_parse_bool(
                        row.get("event_observed"), "event_observed", row_number
                    ),
                    response_achieved=_parse_bool(
                        row.get("response_achieved"),
                        "response_achieved",
                        row_number,
                    ),
                )
            )
    if not subjects:
        raise ValueError("CSV contains no data rows.")
    return subjects


def _fmt_optional(value: Optional[float], digits: int = 2) -> str:
    return "Not reached/estimable" if value is None else f"{value:.{digits}f}"


def format_sca_report(result: SyntheticControlAnalysisResult) -> str:
    lines = [
        "=" * 84,
        " SYNTHETIC CONTROL MATCHING & OUTCOME DIAGNOSTICS",
        " Heuristic balance labels only; not a regulatory determination",
        "=" * 84,
        f"Treated cohort              : {result.trial_arm_size}",
        f"External/RWD control pool   : {result.rwd_pool_size}",
        f"Matched pairs               : {result.matched_pairs_count}",
        f"Treated retention           : {result.matching_retention_pct:.1f}%",
        f"Logit-PS caliper width      : {result.caliper_width:.4f}",
        f"Balance assessment          : {result.balance_assessment}",
        "-" * 84,
        "COVARIATE BALANCE (absolute standardized mean difference)",
    ]
    for item in result.covariate_balance:
        status = "PASS" if item.is_balanced else "REVIEW"
        lines.append(
            f"  {item.covariate_name:<18} pre={item.pre_smd:>7.3f}  "
            f"post={item.post_smd:>7.3f}  {status}"
        )
    lines.extend(
        [
            "-" * 84,
            "MATCHED-SAMPLE OUTCOMES",
            f"  Log-rank HR approximation : {_fmt_optional(result.hazard_ratio, 3)}",
            f"  95% CI                    : "
            f"{_fmt_optional(result.hazard_ratio_ci_low, 3)} to "
            f"{_fmt_optional(result.hazard_ratio_ci_high, 3)}",
            f"  Log-rank p-value          : {result.log_rank_p_value:.4f}",
            f"  Median survival treated   : {_fmt_optional(result.median_survival_treated_months, 1)} months",
            f"  Median survival control   : {_fmt_optional(result.median_survival_synthetic_control_months, 1)} months",
            f"  ORR treated/control       : {result.orr_treated_pct:.1f}% / {result.orr_synthetic_control_pct:.1f}%",
            f"  ORR difference            : {result.att_orr_diff_pct:+.1f} percentage points",
            "-" * 84,
        ]
    )
    lines.extend(f"  - {recommendation}" for recommendation in result.recommendations)
    lines.append("=" * 84)
    return "\n".join(lines)


def _pair_maps(
    pairs: List[MatchedPair],
) -> tuple[Dict[str, tuple[str, str, float]], Dict[str, tuple[str, str, float]]]:
    treated_map: Dict[str, tuple[str, str, float]] = {}
    control_map: Dict[str, tuple[str, str, float]] = {}
    for index, pair in enumerate(pairs, start=1):
        pair_id = f"PAIR-{index:04d}"
        treated_map[pair.treated_id] = (pair_id, pair.control_id, pair.ps_distance)
        control_map[pair.control_id] = (pair_id, pair.treated_id, pair.ps_distance)
    return treated_map, control_map


def write_batch_results(
    subjects: List[SubjectRecord],
    result: SyntheticControlAnalysisResult,
    output_path: str | Path,
) -> None:
    covariates = sorted({name for subject in subjects for name in subject.covariates})
    fieldnames = [
        "subject_id",
        "cohort_arm",
        "match_status",
        "matched_pair_id",
        "matched_subject_id",
        "logit_ps_distance",
        "propensity_score",
        "time_to_event_months",
        "event_observed",
        "response_achieved",
        "balance_assessment",
        "mean_post_smd",
        "hazard_ratio_logrank_approx",
        "log_rank_p_value",
        *covariates,
    ]
    treated_map, control_map = _pair_maps(result.matched_pairs)
    rows: List[Dict[str, Any]] = []
    for subject in subjects:
        match_info = (treated_map if subject.is_treated else control_map).get(
            subject.subject_id
        )
        row: Dict[str, Any] = {
            "subject_id": subject.subject_id,
            "cohort_arm": "TRIAL_TREATED" if subject.is_treated else "RWD_CONTROL",
            "match_status": "MATCHED" if match_info else "UNMATCHED",
            "matched_pair_id": match_info[0] if match_info else "",
            "matched_subject_id": match_info[1] if match_info else "",
            "logit_ps_distance": round(match_info[2], 6) if match_info else "",
            "propensity_score": round(subject.propensity_score or 0.0, 6),
            "time_to_event_months": subject.time_to_event_months,
            "event_observed": subject.event_observed,
            "response_achieved": subject.response_achieved,
            "balance_assessment": result.balance_assessment,
            "mean_post_smd": round(result.mean_absolute_smd_post, 4),
            "hazard_ratio_logrank_approx": ""
            if result.hazard_ratio is None
            else round(result.hazard_ratio, 4),
            "log_rank_p_value": round(result.log_rank_p_value, 5),
        }
        row.update(subject.covariates)
        rows.append(row)

    with Path(output_path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_batch_matching(
    input_path: str | Path,
    output_path: str | Path,
    caliper: float = 0.2,
    quiet: bool = False,
) -> int:
    subjects = load_cohort_from_csv(input_path)
    result = SyntheticControlAgentEngine.generate_synthetic_control_arm(
        subjects, caliper_sd_multiplier=caliper
    )
    write_batch_results(subjects, result, output_path)
    if not quiet:
        print(
            f"Processed {len(subjects)} subjects -> {output_path} "
            f"({result.matched_pairs_count} pairs; {result.balance_assessment})"
        )
    return 0


def interactive_mode() -> int:
    print("\nInteractive synthetic-control demonstration")
    n_treated = int(input("Treated subjects [30]: ").strip() or "30")
    n_control = int(input("Control subjects [100]: ").strip() or "100")
    caliper = float(input("Caliper multiplier [0.2]: ").strip() or "0.2")
    subjects = generate_benchmark_synthetic_cohort(n_treated, n_control)
    result = SyntheticControlAgentEngine.generate_synthetic_control_arm(
        subjects, caliper_sd_multiplier=caliper
    )
    print("\n" + format_sca_report(result))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="clinical-trial-synthetic-control-agent",
        description=(
            "Propensity-score matching, balance diagnostics, and matched-cohort "
            "survival/response summaries for external-control analyses."
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    subparsers = parser.add_subparsers(dest="subcommand")

    batch = subparsers.add_parser("batch", help="Analyze a patient-level CSV cohort")
    batch.add_argument("-i", "--input", required=True, help="Input cohort CSV")
    batch.add_argument("-o", "--output", default="results.csv", help="Output CSV")
    batch.add_argument(
        "--caliper",
        type=float,
        default=0.2,
        help="Caliper multiplier on the SD of the logit propensity score",
    )

    parser.add_argument("--demo", action="store_true", help="Run deterministic demo data")
    parser.add_argument("--file", help="Analyze a JSON or CSV cohort file")
    parser.add_argument(
        "--caliper",
        type=float,
        default=0.2,
        help="Caliper multiplier on the SD of the logit propensity score",
    )
    parser.add_argument("--json", action="store_true", help="Print structured JSON")
    parser.add_argument("--interactive", action="store_true", help="Run interactive demo")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.subcommand == "batch":
            return run_batch_matching(args.input, args.output, caliper=args.caliper)
        if args.interactive:
            return interactive_mode()
        if args.file:
            path = Path(args.file)
            if path.suffix.lower() == ".csv":
                subjects = load_cohort_from_csv(path)
            elif path.suffix.lower() == ".json":
                with path.open("r", encoding="utf-8") as handle:
                    subjects = parse_synthetic_cohort_dict(json.load(handle))
            else:
                raise ValueError("--file must point to a .csv or .json file.")
        else:
            subjects = generate_benchmark_synthetic_cohort()

        result = SyntheticControlAgentEngine.generate_synthetic_control_arm(
            subjects, caliper_sd_multiplier=args.caliper
        )
        if args.json:
            print(json.dumps(result.to_dict(), indent=2, allow_nan=False))
        else:
            print(format_sca_report(result))
        return 0
    except (OSError, ValueError, NotImplementedError, json.JSONDecodeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
