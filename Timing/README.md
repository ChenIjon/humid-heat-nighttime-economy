# Heat and Humidity Reshape Urban Economic Activity at Night

Contact: chenchen1107@connect.hku.hk

The timing analysis estimates how temperature and humidity reshape the hourly distribution of nighttime economic activity using a penalized conditional model with city-specific effects.

Real transaction-derived hourly panels and generated artifacts are restricted and are not distributed. Aggregation does not establish release permission. For the public synthetic main-model demonstration and data release restrictions, see [../README.md](../README.md). The Timing synthetic check below is separate from the server research runs and does not reproduce study estimates.

```bash
Rscript 00_Timing_prepare_panel.R
Rscript 01_Timing_penalty.R
Rscript 02_Timing_model.R
Rscript 03_Timing_inference.R
Rscript 04_Timing_grid_ci.R
```

Run 00-04 in order: prepare the adjusted hourly panel, select penalties, fit the model and hourly profile effects, run the date-block bootstrap, then calculate confidence intervals on the temperature-humidity grid. `timing_helpers.R` contains the shared paths and model functions; `hourly_weather_controls.R` contains the hourly covariate configuration and functions.

The original model structure is retained: shared temperature and humidity natural splines (df = 3) by hour, city-hour baselines, city-specific temperature and humidity slopes for four nighttime periods, and sunset/calendar controls. The conditional likelihood, Newton solver, Laplace penalty selection and circular 7-day date-block bootstrap are retained.

PM2.5, precipitation and wind speed are included as shared linear controls. Each variable is decomposed into its 12-hour nighttime mean and hourly deviation from that mean. Nighttime means interact with hours 1-11 (hour 0 is the reference); hourly deviations interact with hours 0-11, giving 69 candidate control columns. Redundant controls are removed using within-night design differences. Temperature/humidity counterfactuals keep all three controls at their observed values.

Inputs are configured at the top of `hourly_weather_controls.R`:

- Original panel: `data/hourly_profile_panel.csv` under this Timing directory.
- Population-weighted ERA5 hourly data: `/data/ccpeng/NEI/ERA5_population_weighted_hourly_long.csv`, with `city,date,hour,precipitation_mm,wind_speed_mps` (mm per hour; m/s).
- PM2.5 hourly data: `/data/ccpeng/NEI/PM25_city_hourly_2024.csv`, with `city,date,hour,PM25` (µg/m³).

Hourly inputs use Beijing-time calendar dates and clock hours 0-23, with unique city/date/hour keys. Dates may be `YYYYMMDD` or `YYYY-MM-DD`. Matching accounts for midnight and city-specific night start times. A missing, nonfinite or negative hourly covariate excludes the entire night; no imputation is performed. Outcome/exposure filtering also retains complete 12-hour nights and checks that hourly outcomes sum to the nighttime total.

The adjusted panel is written to `data/hourly_profile_panel_adjusted.csv`; model outputs are written to `output_timing_model_adjusted/`. Repeated runs overwrite the same adjusted outputs. The original panel and previous model outputs are preserved. Input hashes and model-version checks reject stale panels and incompatible penalties/models; rerun penalty selection for the adjusted specification.

- `EXPOSURE_OUTCOME`: `cnt` (default) models transaction counts (TC); `amt` models total spending (TS). Average spending per transaction (AS) equals total spending divided by transaction counts (TS/TC). Use a separate result directory for each outcome.
- Set `N_CORES` for your machine. Windows bootstrap/grid evaluation is serial.
- Set `N_BOOT` and `GRID_CI_BOOT` to the same value (default 1000).

Required packages: data.table, mgcv and Matrix; splines and parallel are bundled with R.

Lightweight validation: `Rscript validate_timing.R`. This uses synthetic inputs in a temporary directory to check hourly matching, complete-night filtering, control design and rank, parameter recovery, consistency across all four stages, and counterfactual controls.

The main analyses were conducted using Python (version 3.12.9) and R (version 4.4.1), with the primary statistical models implemented using the mgcv package (version 1.9-3).

This statement describes the server research environment. Demonstration/test runtimes are recorded separately and do not redefine those research versions.
