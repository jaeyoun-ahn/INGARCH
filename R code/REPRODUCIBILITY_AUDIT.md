# Reproducibility audit

Audit date: 2026-09-15

Manuscript checked: `INGARCH_revision_ver4.pdf`, dated 2026-09-14.

## Coverage

All empirical and simulation results in the current paper have an identified
generating script:

| Manuscript item | Generating code |
|---|---|
| Table 1 | `sim_study_INGARCH.R` |
| Table 2 | `sim_study_INGARCH_missing.R` |
| Tables 3-4 and Figure 1 | `Data_summary_INGARCH.R` |
| Table 5 | `run_NB_bootstrap.R` |
| Figure 2 | `heterogeneous_LGPIF_diagnostics.R` |
| Table 6 | `Data_analysis_INGARCH2.R` |
| Table 7 | `run_model_loss_bootstrap.R` |
| Parametric-bootstrap result discussed in Section 5.4 | `run_NB_parametric_bootstrap.R` |

## Corrections made during the audit

- Added the missing `aux_data_manage3.R` local dependency.
- Made `Data_summary_INGARCH.R` source the minimal preprocessing dependency
  (`aux_data_manage2.R`) that it actually uses.
- Changed the default complete-data simulation count from 500 to the 1,000
  replicates reported in the manuscript.
- Removed a hard-coded read of the 1,000-replicate parametric-bootstrap CSV;
  the diagnostic code now uses the bootstrap object generated in the same run,
  so alternate replicate counts work.
- Added run instructions and a package dependency manifest.

## Paper/result discrepancies to resolve before public release

1. The current generated output for the INAR benchmark in Table 6 does not
   exactly match the manuscript. `model_comparison_varying2.csv` gives MSE
   0.222237, MAE 0.131661, and average NLL 0.218929; the manuscript gives
   0.2234, 0.1304, and 0.2193. The other Table 6 rows match the stored output
   after rounding.
2. In Table 7, the manuscript appears to interchange the MSE and MAE row
   labels. In `model_loss_bootstrap_1000_summary.csv`, the intervals
   [-0.0042, 0.0004] and [-0.0039, 0.0028] belong to MAE, while
   [-0.0265, 0.0031] and [-0.0043, 0.0362] belong to MSE.
3. Confirm that the licenses/terms for redistributing `data.RData` and
   `dataout.RData` allow them to be placed in a public repository, and add the
   repository's chosen open-source license for the code.

## Verification limitation

The audit checked the complete manuscript application/simulation sections,
the local source graph, declared package calls, input files, and stored result
tables. R is not installed on the audit host's command path, so an executable
end-to-end smoke test was not run here.
