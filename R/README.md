# Heat and Humidity Reshape Urban Economic Activity at Night

Contact: chenchen1107@connect.hku.hk

`Main_model.R` fits generalized additive mixed models for daily economic activity, with temperature-humidity interactions, city-specific effects, calendar controls and AR(1) correlation. It estimates exposure-weighted marginal effects (EWME) and confidence intervals using a date-block bootstrap.

```bash
Rscript Main_model.R
```

Dependencies: R 4.5.1, mgcv 1.9.3, dplyr 1.1.4; parallel is bundled with R.

