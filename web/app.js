"use strict";

const PYODIDE_INDEX = "https://cdn.jsdelivr.net/pyodide/v314.0.7/full/";
const MAX_FILE_BYTES = 10 * 1024 * 1024;
const PYTHON_FILES = [
  "synthetic_control_arm/__init__.py",
  "synthetic_control_arm/models.py",
  "synthetic_control_arm/engine.py",
  "synthetic_control_arm/cli.py",
];

const elements = {
  runtime: document.getElementById("runtimeStatus"),
  file: document.getElementById("csvFile"),
  dropzone: document.getElementById("dropzone"),
  fileMeta: document.getElementById("fileMeta"),
  sample: document.getElementById("sampleButton"),
  clear: document.getElementById("clearButton"),
  analyze: document.getElementById("analyzeButton"),
  caliper: document.getElementById("caliper"),
  empty: document.getElementById("emptyState"),
  error: document.getElementById("errorState"),
  errorText: document.getElementById("errorText"),
  result: document.getElementById("resultContent"),
  subtitle: document.getElementById("resultsSubtitle"),
  metrics: document.getElementById("metricsGrid"),
  balanceBody: document.getElementById("balanceBody"),
  badge: document.getElementById("assessmentBadge"),
  recommendations: document.getElementById("recommendations"),
  download: document.getElementById("downloadButton"),
  theme: document.getElementById("themeToggle"),
};

let pyodide = null;
let csvText = "";
let csvName = "";
let outputCsv = "";

function setRuntime(message, state = "") {
  elements.runtime.className = `runtime-pill ${state}`.trim();
  elements.runtime.querySelector("span:last-child").textContent = message;
}

function showError(message) {
  elements.empty.hidden = true;
  elements.result.hidden = true;
  elements.error.hidden = false;
  elements.errorText.textContent = message;
  elements.subtitle.textContent = "Review the input and try again.";
  elements.download.disabled = true;
}

function clearError() {
  elements.error.hidden = true;
  elements.errorText.textContent = "";
}

function setCsv(text, name) {
  csvText = text;
  csvName = name;
  outputCsv = "";
  elements.fileMeta.textContent = `${name} · ${(new Blob([text]).size / 1024).toFixed(1)} KB`;
  elements.clear.disabled = false;
  elements.analyze.disabled = !pyodide;
  elements.download.disabled = true;
  clearError();
  elements.result.hidden = true;
  elements.empty.hidden = false;
  elements.subtitle.textContent = "Ready to analyze.";
}

function resetData() {
  csvText = "";
  csvName = "";
  outputCsv = "";
  elements.file.value = "";
  elements.fileMeta.textContent = "Maximum 10 MB";
  elements.clear.disabled = true;
  elements.analyze.disabled = true;
  elements.download.disabled = true;
  clearError();
  elements.result.hidden = true;
  elements.empty.hidden = false;
  elements.subtitle.textContent = "Load a valid CSV to begin.";
}

async function readFile(file) {
  if (!file) return;
  if (file.size > MAX_FILE_BYTES) {
    showError("The selected CSV exceeds the 10 MB browser limit.");
    return;
  }
  if (!file.name.toLowerCase().endsWith(".csv")) {
    showError("Select a .csv file.");
    return;
  }
  setCsv(await file.text(), file.name);
}

function metric(label, value, subtext = "") {
  const card = document.createElement("div");
  card.className = "metric-card";
  const labelNode = document.createElement("div");
  labelNode.className = "metric-label";
  labelNode.textContent = label;
  const valueNode = document.createElement("div");
  valueNode.className = "metric-value";
  valueNode.textContent = value;
  card.append(labelNode, valueNode);
  if (subtext) {
    const sub = document.createElement("div");
    sub.className = "metric-sub";
    sub.textContent = subtext;
    card.append(sub);
  }
  return card;
}

function numberOrNA(value, digits = 2) {
  return value === null || value === undefined ? "N/A" : Number(value).toFixed(digits);
}

function renderResults(data) {
  elements.metrics.replaceChildren();
  const outcomes = data.efficacy_outcomes || {};
  elements.metrics.append(
    metric("Matched pairs", String(data.matched_pairs_count), `${numberOrNA(data.matching_retention_pct, 1)}% treated retained`),
    metric("Post-match SMD", numberOrNA(data.mean_absolute_smd_post, 3), "Mean absolute SMD"),
    metric("Log-rank HR", numberOrNA(outcomes.hazard_ratio_logrank_approx, 3), "Approximation; not Cox"),
    metric("Log-rank p", numberOrNA(outcomes.log_rank_p_value, 4)),
    metric("Median OS", `${numberOrNA(outcomes.median_os_treated_months, 1)} / ${numberOrNA(outcomes.median_os_synthetic_control_months, 1)}`, "Treated / control, months"),
    metric("ORR difference", `${numberOrNA(outcomes.att_orr_difference_pct, 1)} pp`, "Treated − control")
  );

  elements.balanceBody.replaceChildren();
  for (const item of data.covariate_balance || []) {
    const row = document.createElement("tr");
    const name = document.createElement("td");
    const pre = document.createElement("td");
    const post = document.createElement("td");
    const status = document.createElement("td");
    name.textContent = item.covariate;
    pre.textContent = numberOrNA(item.pre_smd, 3);
    post.textContent = numberOrNA(item.post_smd, 3);
    status.textContent = item.is_balanced ? "Balanced" : "Review";
    status.className = item.is_balanced ? "balance-pass" : "balance-review";
    row.append(name, pre, post, status);
    elements.balanceBody.append(row);
  }

  elements.badge.textContent = (data.balance_assessment || "UNCLASSIFIED").replaceAll("_", " ");
  elements.recommendations.replaceChildren();
  for (const note of data.recommendations || []) {
    const item = document.createElement("li");
    item.textContent = note;
    elements.recommendations.append(item);
  }

  elements.empty.hidden = true;
  elements.error.hidden = true;
  elements.result.hidden = false;
  elements.subtitle.textContent = `${data.trial_arm_size} treated · ${data.rwd_pool_size} controls · ${data.matched_pairs_count} matched pairs`;
  elements.download.disabled = !outputCsv;
}

async function initializePython() {
  try {
    setRuntime("Loading Python…");
    pyodide = await loadPyodide({ indexURL: PYODIDE_INDEX });
    pyodide.FS.mkdirTree("/home/pyodide/synthetic_control_arm");
    for (const path of PYTHON_FILES) {
      const response = await fetch(path, { cache: "no-cache" });
      if (!response.ok) throw new Error(`Could not load ${path}.`);
      pyodide.FS.writeFile(`/home/pyodide/${path}`, await response.text(), { encoding: "utf8" });
    }
    pyodide.runPython("import sys\nsys.path.insert(0, '/home/pyodide')\nimport synthetic_control_arm");
    setRuntime("Python ready", "ready");
    elements.analyze.disabled = !csvText;
  } catch (error) {
    console.error(error);
    setRuntime("Runtime failed", "error");
    showError("Python runtime failed to load. Check your network connection and reload the page.");
  }
}

async function analyze() {
  clearError();
  if (!pyodide) {
    showError("Python is still loading.");
    return;
  }
  if (!csvText) {
    showError("Choose a CSV or load the sample first.");
    return;
  }
  const caliper = Number(elements.caliper.value);
  if (!Number.isFinite(caliper) || caliper <= 0 || caliper > 2) {
    showError("Caliper multiplier must be greater than 0 and no more than 2.");
    return;
  }

  elements.analyze.disabled = true;
  elements.analyze.textContent = "Analyzing…";
  elements.subtitle.textContent = "Running Python matching engine…";
  try {
    pyodide.globals.set("browser_csv_text", csvText);
    pyodide.globals.set("browser_caliper", caliper);
    const serialized = pyodide.runPython(`
from pathlib import Path
import json
from synthetic_control_arm import SyntheticControlAgentEngine
from synthetic_control_arm.cli import load_cohort_from_csv, write_batch_results

_input = Path('/tmp/browser_input.csv')
_output = Path('/tmp/browser_output.csv')
_input.write_text(browser_csv_text, encoding='utf-8')
_subjects = load_cohort_from_csv(_input)
_result = SyntheticControlAgentEngine.generate_synthetic_control_arm(
    _subjects, caliper_sd_multiplier=float(browser_caliper)
)
write_batch_results(_subjects, _result, _output)
json.dumps({
    'result': _result.to_dict(),
    'csv': _output.read_text(encoding='utf-8')
}, allow_nan=False)
    `);
    const payload = JSON.parse(serialized);
    outputCsv = payload.csv;
    renderResults(payload.result);
  } catch (error) {
    console.error(error);
    const message = String(error).replace(/^PythonError:\s*/, "").split("\n").slice(-1)[0];
    showError(message || "Unknown analysis error.");
  } finally {
    elements.analyze.textContent = "Analyze cohort";
    elements.analyze.disabled = !pyodide || !csvText;
  }
}

function downloadResults() {
  if (!outputCsv) return;
  const blob = new Blob([outputCsv], { type: "text/csv;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  const stem = (csvName || "cohort").replace(/\.csv$/i, "").replace(/[^a-zA-Z0-9._-]+/g, "-");
  link.download = `${stem}-matched-results.csv`;
  document.body.append(link);
  link.click();
  link.remove();
  URL.revokeObjectURL(url);
}

function initializeTheme() {
  const saved = localStorage.getItem("sc-theme");
  const preferredDark = window.matchMedia("(prefers-color-scheme: dark)").matches;
  const theme = saved || (preferredDark ? "dark" : "light");
  document.documentElement.dataset.theme = theme;
  elements.theme.setAttribute("aria-label", theme === "dark" ? "Switch to light theme" : "Switch to dark theme");
}

function toggleTheme() {
  const next = document.documentElement.dataset.theme === "dark" ? "light" : "dark";
  document.documentElement.dataset.theme = next;
  localStorage.setItem("sc-theme", next);
  elements.theme.setAttribute("aria-label", next === "dark" ? "Switch to light theme" : "Switch to dark theme");
}

elements.file.addEventListener("change", () => readFile(elements.file.files[0]));
elements.sample.addEventListener("click", async () => {
  try {
    const response = await fetch("sample.csv", { cache: "no-cache" });
    if (!response.ok) throw new Error("Sample CSV is unavailable.");
    setCsv(await response.text(), "sample.csv");
  } catch (error) {
    showError(String(error.message || error));
  }
});
elements.clear.addEventListener("click", resetData);
elements.analyze.addEventListener("click", analyze);
elements.download.addEventListener("click", downloadResults);
elements.theme.addEventListener("click", toggleTheme);

for (const eventName of ["dragenter", "dragover"]) {
  elements.dropzone.addEventListener(eventName, (event) => {
    event.preventDefault();
    elements.dropzone.classList.add("dragover");
  });
}
for (const eventName of ["dragleave", "drop"]) {
  elements.dropzone.addEventListener(eventName, (event) => {
    event.preventDefault();
    elements.dropzone.classList.remove("dragover");
  });
}
elements.dropzone.addEventListener("drop", (event) => {
  const [file] = event.dataTransfer.files;
  readFile(file);
});

initializeTheme();
initializePython();
