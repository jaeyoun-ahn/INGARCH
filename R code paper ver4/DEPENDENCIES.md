# Dependency manifest

## Direct R package dependencies

| Package | Version constraint | Used by |
|---|---:|---|
| RTMB | >= 1.9 | Model fitting, automatic differentiation, simulations, bootstraps, diagnostics |
| data.table | any current CRAN version | Descriptive summaries |
| ggplot2 | any current CRAN version | Coverage plots |
| lme4 | any current CRAN version | Descriptive mixed-model calculation |

The development copy of RTMB was version 1.9, built under R 4.6.1 on
64-bit Windows. RTMB declares compiled dependencies including TMB and Rcpp;
these should be installed from the package repository rather than committed as
a platform-specific `Rlib` directory.

## Base and recommended R packages

The code also uses `parallel`, `grid`, `stats`, `utils`, `graphics`, `Matrix`,
and `MASS`. These are included with, or installed as recommended dependencies
of, a standard R distribution and RTMB installation.

## Local file dependency graph

```text
data.RData + dataout.RData
  -> aux_data_manage2.R
     -> aux_data_manage3.R
        -> aux_data_manage4.R
           -> Data_analysis_INGARCH2.R
           -> run_NB_bootstrap.R
           -> run_NB_parametric_bootstrap.R
           -> heterogeneous_LGPIF_diagnostics.R

funs_INGARCH.R
  -> Data_analysis_INGARCH2.R
  -> run_NB_bootstrap.R
  -> run_NB_parametric_bootstrap.R
  -> sim_study_INGARCH.R

Data_analysis_INGARCH2.R
  -> run_model_loss_bootstrap.R

sim_study_INGARCH.R
  -> sim_study_INGARCH_missing.R

data.RData + dataout.RData + aux_data_manage2.R
  -> Data_summary_INGARCH.R
```
