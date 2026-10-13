# Heat and Humidity Reshape Urban Economic Activity at Night

R code accompanying the study *Heat and Humidity Reshape Urban Economic Activity at Night*. The repository contains the main city–date statistical model, a separate analysis of within-night transaction timing, and synthetic examples for inspecting the code without access to restricted transaction data.

**Contact:** chenchen1107@connect.hku.hk  
**License:** [MIT](LICENSE)

## Quick start

Install R and the required packages:

```bash
Rscript -e "install.packages(c('mgcv', 'dplyr', 'data.table', 'Matrix'), repos='https://cloud.r-project.org')"
```

From the repository root, run the main-model demonstration:

```bash
Rscript examples/run_demo.R
```

The included [synthetic dataset](examples/main_panel_synthetic.csv) has 360 city–date observations (3 cities × 120 days). The demo fits the transaction-count (TC), total-spending (TS), and average-spending (AS) models and runs **4 bootstrap replicates per outcome** on one CPU worker. It writes results to `examples/output_demo/`, including `ewme_bootstrap_CI.csv`, model diagnostics, and a summary report. Existing files with the same names in that directory may be overwritten. A recorded Windows test completed the analysis in approximately **4 seconds**; runtime varies by machine. Demo estimates are illustrative, not study findings.

To regenerate the synthetic input or run the separate Timing checks:

```bash
Rscript examples/generate_demo_data.R
Rscript Timing/validate_timing.R
```

The Timing validation checks data matching, design construction, parameter recovery, and counterfactual handling; it is **not** an end-to-end run of the full Timing analysis. See [demo verification details](examples/TESTING.md).

## Code and environment

| Location | Purpose |
| --- | --- |
| [`Main_model.R`](Main_model.R) | Main generalized additive mixed models; exposure-weighted mean effects (EWMEs); date-block bootstrap inference. |
| [`Timing/`](Timing/) | Conditional model of the hourly distribution of nighttime activity; penalty selection, estimation, bootstrap, and grid confidence intervals. |
| [`examples/`](examples/) | Fully synthetic main-model input, generator, demonstration script, and sample outputs. |

The study analyses used **Python 3.12.9**, **R 4.4.1**, and **mgcv 1.9-3**. Python is not required to run the R scripts in this repository. The public main-model demo was tested on **Windows 11** with **R 4.5.1**, **mgcv 1.9-3**, and **dplyr 1.1.4**; this test does not establish compatibility with every research-server configuration. The Timing code additionally uses `data.table` and `Matrix`. Base R provides `parallel` and `splines`. No specialized hardware is required for the small demo; full research runs are substantially more demanding.

## Data access and required input

The original Alipay transaction-derived data and prepared research panels are **restricted and are not distributed**. The synthetic examples contain no real transaction observations and cannot reproduce the numerical estimates reported in the paper.

For an authorized main-model run, place `main_panel.csv` in `data/` next to `Main_model.R` (or edit `PANEL_FILE` in the script). Each row represents one city–date observation, with unique `city,DATE` keys:

| Field | Description |
| --- | --- |
| `city`, `DATE` | City name in `MAIN_CITIES`; calendar date (`YYYY-MM-DD`). |
| `TEMP`, `RH` | Temperature (°C) and relative humidity (%; 0–100). |
| `WS`, `PRECIP`, `PM25` | Wind speed (m/s), nonnegative precipitation (mm), and PM2.5 (µg/m³). |
| `calendar_type` | Monday–Sunday, `Statutory holiday`, or `Make-up workday`. |
| `model_cnt` | Transaction count (TC). |
| `model_amt` | Total spending (TS), in original monetary units. |
| `model_amt_per_cnt` | Average spending (AS): `model_amt / model_cnt`. |

The script excludes rows with invalid or missing required fields. After filtering, the date-block bootstrap requires a balanced city–date panel and adequate covariate variation for spline fitting. Supply **unscaled outcomes**: the code divides TS by 100,000 and AS by 10 before fitting.

## Main-model research run

With an authorized input panel:

```bash
Rscript Main_model.R
```

The main model uses negative-binomial (TC) and Gamma-log (TS and AS) families, smooth temperature and humidity terms, their interaction, city-specific exposure responses, temporal and meteorological controls, and AR(1) correlation. Default inference uses **1,000 circular 7-day calendar-date block-bootstrap replicates per outcome**, with the same date draws across cities and fixed smoothing parameters during bootstrap refits.

Edit `PANEL_FILE`, `OUTPUT_DIR`, and `main_model_defaults()` in `Main_model.R` to change paths or settings. The default run is serial (`n_cores = 1`); parallel execution is supported on Unix. The full study-scale run may require substantial CPU time and memory; no general production runtime is claimed. Use a new output directory or set `force_refit = TRUE` when changing inputs or settings to avoid reusing incompatible bootstrap caches.

Key outputs are written under `output_main/`:

| Output | Contents |
| --- | --- |
| `ewme_bootstrap_CI.csv` | EWME point estimates and bootstrap confidence intervals. |
| `model_diagnostics.csv`, `model_run_status.csv` | Model fit and run diagnostics. |
| `{TC,TS,AS}_bootstrap_status.csv` | Bootstrap convergence and error records. |
| `{TC,TS,AS}_bootstrap_ewme_replicates.csv` | Replicate-level EWME estimates. |
| `models/*.rds`, `bootstrap_bundles/` | Fitted models and cached bootstrap results. |
| `analysis_data_model_input.csv` | Filtered and scaled model input. |

**Do not publish outputs derived from restricted data without authorization**, including filtered analysis inputs, fitted models, and bootstrap artifacts.

## Timing analysis

The [`Timing/README.md`](Timing/README.md) documents the separate **00–04** workflow, its hourly input files, and configuration. The conditional model includes city–hour baselines, temperature/humidity splines and city-specific period slopes, sunset/calendar controls, and hourly PM2.5, precipitation, and wind-speed controls. Running the complete Timing workflow requires authorized hourly data and locally configured input paths; `Timing/validate_timing.R` provides the public synthetic checks.
