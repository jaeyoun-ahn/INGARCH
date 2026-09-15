# Reproducibility code for INGARCH paper, version 4

This folder accompanies `INGARCH_revision_ver4.pdf` (September 14, 2026).
Run all commands from this directory because the scripts use relative paths.

## Requirements

- R 4.6.0 or later (the development environment used R 4.6.1 on Windows)
- CRAN packages: `RTMB` 1.9 or later, `data.table`, `ggplot2`, and `lme4`

Install direct package dependencies with:

```r
install.packages(c("RTMB", "data.table", "ggplot2", "lme4"))
```

`RTMB` installs its compiled dependencies (including `TMB` and `Rcpp`). The
scripts prefer a local `Rlib` directory when one exists but also search the
normal R library paths, so the binary library used by the authors is not
required or portable across operating systems.

## Included inputs

- `data.RData` and `dataout.RData`: data inputs used by the empirical analysis
- `aux_data_manage2.R`, `aux_data_manage3.R`, and `aux_data_manage4.R`: data
  preparation for the time-varying coverage analysis
- `funs_INGARCH.R`: model likelihoods, prediction, and bootstrap helpers

## Main scripts and manuscript outputs

1. `sim_study_INGARCH.R` reproduces the complete-data simulation underlying
   Table 1. Set `INGARCH_SIM_REPLICATES` to change the default of 1,000
   replicates; the default run can take a long time.
2. `sim_study_INGARCH_missing.R` reproduces the missing-data simulation in
   Table 2. Set `INGARCH_MISSING_SIM_REPLICATES` to change the default of
   1,000 replicates.
3. `Data_summary_INGARCH.R` produces the summaries in Tables 3-4 and Figure 1.
4. `Data_analysis_INGARCH2.R` fits the proposed and benchmark models and
   produces the values in Table 6.
5. `run_NB_bootstrap.R` produces the nonparametric bootstrap estimates and
   intervals in Table 5. Set `NB_BOOTSTRAP_REPLICATES` to change the default
   of 1,000 replicates.
6. `run_NB_parametric_bootstrap.R` performs the likelihood-ratio parametric
   bootstrap discussed in Section 5.4. Set `NB_PARAMETRIC_BOOTSTRAP_REPLICATES`
   to change the default of 1,000 replicates.
7. `heterogeneous_LGPIF_diagnostics.R` produces the diagnostics in Figure 2.
8. `run_model_loss_bootstrap.R` produces the paired model-comparison bootstrap
   in Table 7. Set `MODEL_LOSS_BOOTSTRAP_REPLICATES` to change the default of
   1,000 replicates.

The bootstrap and simulation scripts use 11 parallel workers by default.
Reduce `workers` in the relevant script on smaller machines.

## Quick smoke test

Before launching the long simulations, verify that the analysis loads:

```r
source("Data_analysis_INGARCH2.R")
```

The large `.rds`, `.csv`, `.png`, and `.tex` result files in the authors'
working directory are generated outputs, not inputs required to run the code.
