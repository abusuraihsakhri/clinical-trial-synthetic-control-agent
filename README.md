# Clinical Trial Synthetic Control Agent

> **Domain:** Real-World Evidence (RWE), Single-Arm Oncology Trials, & Observational Causal Inference  
> **Regulatory Guidelines & Standards:** FDA Real-World Evidence Guidance (2023/2024), EMA Guideline on Registry-Based Studies, ISPOR-ISPE Good Research Practices

---

## 📖 Overview

The **Clinical Trial Synthetic Control Agent** is a biostatistical computing engine designed to construct regulatory-grade **Synthetic Control Arms (SCA)** from Real-World Data (RWD) such as electronic health records (EHR), patient registries, and historical control cohorts.

In single-arm Phase II or pivotal oncology trials where randomized concurrent control groups are ethically unfeasible or logistically prohibitive (e.g., rare biomarker-driven oncology or orphan indications), this system implements causal inference architectures:
- **Propensity Score Modeling**: Regularized multivariable logistic regression modeling probability of treatment assignment given high-dimensional baseline confounders.
- **Nearest-Neighbor Caliper Matching**: Greedy 1:1 matching on the logit of the propensity score with customizable caliper restrictions.
- **Covariate Balance Diagnostics**: Standardized Mean Difference (SMD / Cohen's $d$) profiling before and after matching.
- **Survival Analytics & Causal Estimators**: Kaplan-Meier survival curves, Mantel-Haenszel log-rank hypothesis testing, and Cox-approximated Hazard Ratios (HR) with 95% confidence intervals.
- **Binary Response Comparative Diagnostics**: Objective Response Rate (ORR) and Average Treatment Effect on the Treated (ATT) estimation.
- **Automated Regulatory Adequacy Tiering**: Systematic grading against FDA/EMA clinical adequacy benchmarks.

---

## 📐 Methodological & Biostatistical Formulations

### 1. Propensity Score Specification
For patient $i$ with baseline covariate vector $\mathbf{X}_i = [X_{i1}, X_{i2}, \dots, X_{ip}]^T$, the propensity score $e(\mathbf{X}_i)$ denotes the conditional probability of assignment to the experimental trial arm ($T_i = 1$) versus the observational RWD registry pool ($T_i = 0$):

$$e(\mathbf{X}_i) = P(T_i = 1 \mid \mathbf{X}_i) = \sigma(\boldsymbol{\beta}^T \mathbf{X}_i) = \frac{1}{1 + \exp\left(-\left(\beta_0 + \sum_{j=1}^p \beta_j X_{ij}\right)\right)}$$

Optimization is conducted using regularized maximum likelihood via batch gradient descent with L2 penalty:

$$\mathcal{L}(\boldsymbol{\beta}) = -\frac{1}{N}\sum_{i=1}^N \left[ T_i \ln(e(\mathbf{X}_i)) + (1 - T_i)\ln(1 - e(\mathbf{X}_i)) \right] + \frac{\lambda}{2}\|\boldsymbol{\beta}_{1:p}\|_2^2$$

### 2. Caliper Restriction on the Logit Propensity Score
Matching is executed on the logit of the propensity score ($q_i = \text{logit}(e(\mathbf{X}_i)) = \ln\frac{e(\mathbf{X}_i)}{1 - e(\mathbf{X}_i)}$) within a maximum allowable distance $\delta$:

$$\delta = \kappa \cdot \text{SD}\left(\text{logit}(e(\mathbf{X}))\right)$$

where standard regulatory guidelines recommend $\kappa \in [0.1, 0.2]$ (Austin 2011). A treated subject $i$ is matched to historical subject $j$ satisfying:

$$|q_i - q_j| \le \delta \quad \text{and} \quad j = \arg\min_{k \in \mathcal{C}_{\text{available}}} |q_i - q_k|$$

### 3. Covariate Balance Metric (Standardized Mean Difference - SMD)
For each covariate $j$, the pre- and post-match balance is assessed using Cohen's $d$:

$$\text{SMD}_j = \frac{|\bar{X}_{T, j} - \bar{X}_{C, j}|}{\sqrt{\frac{s_{T, j}^2 + s_{C, j}^2}{2}}}$$

Where:
- $\bar{X}_{T, j}, \bar{X}_{C, j}$ are the sample means in the treated and control cohorts.
- $s_{T, j}^2, s_{C, j}^2$ are the sample variances.
- An $\text{SMD}_j < 0.10$ indicates adequate balance under FDA/EMA RWE guidelines.

### 4. Mantel-Haenszel Log-Rank Statistic & Hazard Ratio
Across ordered unique event times $t_1 < t_2 < \dots < t_K$:

$$E_{T, k} = n_{T, k} \cdot \left(\frac{d_k}{n_k}\right), \quad \text{Var}(O_{T, k}) = \frac{n_{T, k} n_{C, k} d_k (n_k - d_k)}{n_k^2 (n_k - 1)}$$

$$\chi^2 = \frac{\left( \sum_{k=1}^K (O_{T, k} - E_{T, k}) \right)^2}{\sum_{k=1}^K \text{Var}(O_{T, k})}, \quad p = \text{erfc}\left(\sqrt{\frac{\chi^2}{2}}\right)$$

The approximated hazard ratio (HR) and its 95% confidence bounds are derived as:

$$\widehat{\text{HR}} \approx \frac{\sum O_{T, k} / \sum E_{T, k}}{\sum O_{C, k} / \sum E_{C, k}}, \quad 95\% \text{ CI} = \exp\left(\ln(\widehat{\text{HR}}) \pm 1.96 \cdot \sqrt{\frac{1}{\sum E_{T, k}} + \frac{1}{\sum E_{C, k}}}\right)$$

### 5. Average Treatment Effect on the Treated (ATT) for Response
For binary endpoint (Objective Response Rate - ORR):

$$\widehat{\text{ATT}}_{\text{ORR}} = \frac{1}{N_{\text{matched}}} \sum_{i \in \mathcal{M}_T} Y_i - \frac{1}{N_{\text{matched}}} \sum_{j \in \mathcal{M}_C} Y_j$$

---

## 📊 Regulatory Adequacy Tiering Criteria

| Tier | Covariate Balance (SMD) | Cohort Overlap / Retention | Regulatory Use Case |
|:-----|:------------------------|:---------------------------|:--------------------|
| **`STRONG_REGULATORY_GRADE`** | All covariates $\text{SMD} < 0.10$, mean $\text{SMD} \le 0.08$ | $\ge 75\%$ treated arm matched | Pivotal filing supportive evidence (FDA/EMA submission) |
| **`ACCEPTABLE_SUPPORTIVE`** | Mean $\text{SMD} \le 0.15$, minor residual imbalance | $\ge 50\%$ treated arm matched | Phase II signal detection, secondary efficacy endpoint |
| **`HIGH_CONFOUNDING_RISK`** | Any key confounder $\text{SMD} > 0.15$ or mean $> 0.15$ | $< 50\%$ treated arm matched | Exploratory internal discovery; requires cohort expansion |

---

## 💻 CLI Quickstart & Usage

### 1. Batch Matching with CSV Files
Run batch matching on patient cohorts formatted in CSV:

```bash
python cli.py batch -i sample.csv -o out_matched.csv
```
Or using full parameter flags:
```bash
python cli.py batch --input sample.csv --output out_matched.csv --caliper 0.2
```

### 2. Demonstration Analysis (Synthetic Oncology Cohort)
Run an automated benchmark analysis comparing a single-arm trial arm ($N=40$) against an RWD registry pool ($N=150$):

```bash
python cli.py --demo
```

### 3. JSON Output Mode
Export structured causal matching diagnostics and survival metrics as JSON:

```bash
python cli.py --demo --json
```

### 4. Interactive Cohort Matching Wizard
Interactively define sample sizes, caliper widths, and inspect matched survival results:

```bash
python cli.py --interactive
```

---

## 📋 Input Data Format (`sample.csv`)

Input CSV files support patient-level trial and observational records. Common clinical baseline covariates (`age`, `ecog`, `prior_lines`, `ldh`, etc.) are automatically extracted and balanced:

```csv
subject_id,is_treated,age,ecog,prior_lines,ldh,time_to_event_months,event_observed,response_achieved
TRIAL-PT-001,True,58.2,0.0,1.0,212.4,24.5,True,True
TRIAL-PT-002,True,62.1,1.0,2.0,230.1,18.2,False,True
RWD-REG-001,False,60.1,0.0,1.0,215.0,20.4,True,True
RWD-REG-002,False,63.5,1.0,2.0,238.2,16.1,True,False
```

### Field Definitions:
- `subject_id`: Unique patient identifier (e.g. `TRIAL-PT-001` or `RWD-REG-001`).
- `is_treated`: Boolean (`True` for single-arm trial experimental cohort; `False` for RWD observational control pool).
- `time_to_event_months`: Total follow-up duration (in months) to event or right-censoring.
- `event_observed`: Boolean (`True` if target clinical event occurred; `False` if right-censored).
- `response_achieved`: Boolean (`True` for confirmed objective response [RECIST CR/PR]; `False` otherwise).
- *(Additional columns)*: Numerical baseline covariates automatically standardized and included in the propensity model.

---

## 🧪 Testing & Quality Assurance

Run the automated test suite with pytest:

```bash
python -m pytest -p no:zarr
```

Run the end-to-end CLI batch smoke test:

```bash
python cli.py batch -i sample.csv -o out_smoke.csv
python -c "import os; assert os.path.exists('out_smoke.csv'); os.remove('out_smoke.csv')"
```

---

## 📄 License

MIT License. See [LICENSE](LICENSE) for details.
