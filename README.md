# Clinical Trial Synthetic Control Agent

### [Open the Live Application →](https://abusuraihsakhri.github.io/clinical-trial-synthetic-control-agent/)

A dependency-free Python tool for propensity-score matching and matched-cohort diagnostics in external-control analyses.

## What it does

The project compares a treated cohort with an observational or historical control pool using a compact, reproducible workflow:

- regularized logistic propensity-score estimation from numeric baseline covariates;
- 1:1 nearest-neighbor matching without replacement;
- a caliper defined on the **logit propensity-score** scale;
- pre/post absolute standardized mean differences (SMDs);
- Kaplan–Meier median survival estimation;
- log-rank testing with a score-test hazard-ratio approximation;
- matched-sample objective response rate (ORR) summaries; and
- CSV and JSON output for downstream review.

The browser interface runs the same Python engine client-side through Pyodide. Cohort files are processed in the browser and are not uploaded to this repository or to an application server.

## Scope and interpretation

This is a research and educational utility, not a validated statistical package, medical device, or regulatory decision system. The current engine treats covariates as numeric variables; categorical variables require suitable encoding before use. Matching addresses measured covariates only and does not remove bias from unmeasured confounding, data-quality problems, selection mechanisms, or endpoint incompatibility.

The public API retains the historical field `regulatory_adequacy_tier` for backward compatibility, but its enum names are **tool-defined legacy labels**. New interfaces expose the neutral `balance_assessment` labels `WELL_BALANCED`, `PARTIAL_BALANCE`, and `POOR_BALANCE_OR_OVERLAP`. These labels are not FDA, EMA, or other regulator-defined grades.

The default caliper multiplier is `0.20 × SD(logit propensity score)`, a commonly cited methodological recommendation for propensity-score matching rather than a regulatory requirement. The tool uses absolute SMD `< 0.10` as a conventional measured-covariate balance flag. Interpretation should remain study-specific.

## Browser application

The static application is under `web/` and is deployed through GitHub Pages. It supports:

- local CSV selection or drag-and-drop;
- the included demonstration cohort;
- adjustable caliper multiplier;
- matched-pair retention and balance summaries;
- matched-sample survival and ORR diagnostics;
- downloadable annotated CSV output;
- light and dark themes; and
- responsive keyboard-accessible controls.

Python runs in WebAssembly through Pyodide, so no Flask, FastAPI, or server-side Python process is required.

## CLI

Python 3.10 or newer is required. The runtime package has no third-party Python dependencies.

```bash
python -m pip install -e .
clinical-trial-synthetic-control-agent --demo
```

Analyze the included CSV and write an annotated output file:

```bash
clinical-trial-synthetic-control-agent batch \
  --input sample.csv \
  --output matched-results.csv \
  --caliper 0.20
```

Structured demo output:

```bash
clinical-trial-synthetic-control-agent --demo --json
```

The legacy root command remains available:

```bash
python cli.py --demo
```

## CSV format

Required columns:

| Column | Meaning |
| --- | --- |
| `subject_id` | Unique subject identifier within the input file |
| `is_treated` | `true/false`, `1/0`, `yes/no`, or equivalent supported boolean form |
| `time_to_event_months` | Non-negative follow-up/event time |
| `event_observed` | Whether the time-to-event endpoint occurred |
| `response_achieved` | Binary response indicator |

Add at least one numeric baseline covariate column, for example `age`, `ecog`, `prior_lines`, or `ldh`. Header names must be non-empty and unique after surrounding whitespace is removed. Every analyzed subject must have a finite value for every selected covariate. The CSV and JSON parsers fail explicitly on missing, non-finite, or invalid values instead of silently coercing them to zero.

## Statistical notes

The propensity model is an L2-regularized logistic regression implemented with batch gradient descent and z-standardized predictors. Matching is greedy 1:1 nearest-neighbor matching without replacement on `logit(PS)`. The reported caliper is in the same logit-PS units used for the match distance.

SMD is computed as the absolute mean difference divided by the square root of the mean within-group variances. If both groups are constant but have different means, the SMD is treated as unbounded rather than zero. Kaplan–Meier median survival is reported as not reached when the estimated survival curve never falls to 0.5 or below.

The reported hazard ratio is a **log-rank score-test approximation**. It is not a fitted Cox proportional-hazards model and should not be described or interpreted as one.

## Development and testing

```bash
python -m pip install -e .
python -m pip install pytest ruff
ruff check .
python -m pytest -q
python -m synthetic_control_arm.cli --help
python cli.py batch -i sample.csv -o out_smoke.csv
```

CI runs the test suite, lint checks, package/console-script smoke tests, and a batch export smoke test on Python 3.10–3.12.

## Privacy and security

The CLI reads local files only. The browser application processes cohort data in browser memory after static assets and the Pyodide runtime have loaded. It does not send cohort contents to a project backend. Avoid using direct identifiers or protected health information unless your local workflow, governance, and authorization permit it.

The repository contains no application secrets or client-side API credentials.

## Browser compatibility

The browser interface targets current desktop and mobile versions of Chromium-based browsers, Firefox, and Safari with WebAssembly support. Initial load requires network access to retrieve the version-pinned Pyodide runtime from jsDelivr; subsequent runtime caching is browser-dependent.

## Method references

- Austin PC. Optimal caliper widths for propensity-score matching when estimating differences in means and differences in proportions in observational studies. *Pharmaceutical Statistics*. 2011;10(2):150–161. doi:10.1002/pst.433.
- U.S. Food and Drug Administration. Considerations for the Design and Conduct of Externally Controlled Trials for Drug and Biological Products. Draft Guidance for Industry, February 2023.

## License

MIT License. See [LICENSE](LICENSE).
