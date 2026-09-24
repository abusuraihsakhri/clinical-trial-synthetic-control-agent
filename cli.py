#!/usr/bin/env python3
"""Backward-compatible repository-root CLI wrapper."""

from synthetic_control_arm.cli import (
    VERSION,
    build_parser,
    format_sca_report,
    generate_benchmark_synthetic_cohort,
    interactive_mode,
    load_cohort_from_csv,
    main,
    run_batch_matching,
    write_batch_results,
)

__all__ = [
    "VERSION",
    "build_parser",
    "format_sca_report",
    "generate_benchmark_synthetic_cohort",
    "interactive_mode",
    "load_cohort_from_csv",
    "main",
    "run_batch_matching",
    "write_batch_results",
]

if __name__ == "__main__":
    raise SystemExit(main())
