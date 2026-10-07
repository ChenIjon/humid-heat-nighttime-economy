# Heat and Humidity Reshape Urban Economic Activity at Night

Contact: chenchen1107@connect.hku.hk

The timing analysis estimates how temperature and humidity reshape the hourly distribution of nighttime economic activity using a penalized conditional model with city-specific effects.

```bash
Rscript 01_Timing_penalty.R
Rscript 02_Timing_model.R
Rscript 03_Timing_inference.R
Rscript 04_Timing_grid_ci.R
```

Run 01-04 in order: penalty selection, model fitting and hourly profile effects, date-block bootstrap, then confidence intervals on the temperature-humidity grid. `timing_helpers.R` contains the shared paths and functions.

- `EXPOSURE_OUTCOME`: `cnt` (default) models transaction counts (TC); `amt` models total spending (TS). Average spending per transaction (AS) equals total spending divided by transaction counts (TS/TC). Use a separate result directory for each outcome.
- Set `N_CORES` for your machine. Windows bootstrap/grid evaluation is serial.
- Set `N_BOOT` and `GRID_CI_BOOT` to the same value (default 1000).

Dependencies: R 4.5.1, data.table 1.17.8, mgcv 1.9.3, Matrix 1.7.3; splines and parallel are bundled with R.
