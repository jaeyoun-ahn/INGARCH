source("sim_study_INGARCH.R")

run_missingness_study <- function(Delta=c(.5, .8), n=1000,
                                  missing=c(0, .05, .10, .20),
                                  replicates=1000, workers=11, seed=200){
  design <- expand.grid(replicate=seq_len(replicates), n=n, Delta=Delta)
  cluster <- parallel::makeCluster(min(workers, nrow(design)))
  parallel::clusterExport(
    cluster,
    c("design", "missing", "beta_true", "gamma_true",
      "sample_observation_patterns", "simulate_nb_ingarch", "fit_simulation",
      "rtmb_state_nb_nll", "make_state_nb_objective"),
    envir=environment()
  )
  invisible(parallel::clusterEvalQ(cluster, {
    .libPaths(c("Rlib", .libPaths()))
    library(RTMB)
  }))
  parallel::clusterSetRNGStream(cluster, seed)
  result <- parallel::parLapply(cluster, seq_len(nrow(design)), function(i){
    setting <- design[i, ]
    complete <- simulate_nb_ingarch(setting$Delta, setting$n, 0)
    complete_fit <- fit_simulation(complete, keep_fit=TRUE)

    lapply(missing, function(rate){
      data <- complete
      data$W <- sample_observation_patterns(setting$n, rate)
      data$Y_train[!data$W] <- NA
      fit <- if(rate==0) complete_fit else
        fit_simulation(data, complete_fit$par, keep_fit=TRUE)
      transform(fit$result, replicate=setting$replicate, n=setting$n,
                Delta=setting$Delta, missing=rate,
                realized_missing=1-mean(data$W))
    })
  })
  parallel::stopCluster(cluster)
  estimates <- do.call(rbind, unlist(result, recursive=FALSE))
  truth <- c(beta_true, delta=NA_real_, gamma_1_0=gamma_true)
  estimates$truth <- truth[estimates$parameter]
  estimates$truth[estimates$parameter=="delta"] <-
    estimates$Delta[estimates$parameter=="delta"]
  estimates$covered <- abs(estimates$estimate-estimates$truth) <=
    1.96*estimates$SE

  complete_rows <- estimates[estimates$missing==0, ]
  complete_estimate <- complete_rows$estimate[
    match(paste(estimates$replicate, estimates$n, estimates$Delta,
                estimates$parameter),
          paste(complete_rows$replicate, complete_rows$n,
                complete_rows$Delta, complete_rows$parameter))
  ]
  estimates$change_from_complete <- estimates$estimate-complete_estimate

  summary <- do.call(rbind, lapply(
    split(estimates, list(estimates$n, estimates$missing, estimates$Delta,
                          estimates$parameter)),
    function(x){
      valid <- is.finite(x$SE)
      centered <- x$estimate-mean(x$estimate)
      spread <- sqrt(mean(centered^2))
      data.frame(
        n=x$n[1], Delta=x$Delta[1], missing=x$missing[1],
        realized_missing=mean(x$realized_missing), parameter=x$parameter[1],
        truth=x$truth[1], mean=mean(x$estimate),
        bias=mean(x$estimate-x$truth), variance=var(x$estimate),
        empirical_SD=sd(x$estimate), MSE=mean((x$estimate-x$truth)^2),
        mean_SE=mean(x$SE[valid]),
        SE_to_SD=mean(x$SE[valid])/sd(x$estimate),
        mean_Wald_width=mean(2*1.96*x$SE[valid]),
        coverage=mean(x$covered[valid]), valid_Wald_rate=mean(valid),
        mean_change_from_complete=mean(x$change_from_complete),
        SD_change_from_complete=sd(x$change_from_complete),
        RMSE_change_from_complete=sqrt(mean(x$change_from_complete^2)),
        skewness=mean(centered^3)/spread^3,
        excess_kurtosis=mean(centered^4)/spread^4-3,
        Shapiro_p=shapiro.test(x$estimate)$p.value,
        convergence_rate=mean(x$convergence==0)
      )
    }
  ))
  rownames(summary) <- NULL

  baseline <- summary[summary$missing==0,
                      c("n", "Delta", "parameter", "variance", "MSE",
                        "mean_Wald_width", "coverage")]
  names(baseline)[4:7] <- paste0(names(baseline)[4:7], "_complete")
  summary <- merge(summary, baseline,
                   by=c("n", "Delta", "parameter"), sort=FALSE)
  summary$variance_inflation <- summary$variance/summary$variance_complete
  summary$relative_efficiency <- summary$variance_complete/summary$variance
  summary$MSE_inflation <- summary$MSE/summary$MSE_complete
  summary$Wald_width_ratio <- summary$mean_Wald_width/
    summary$mean_Wald_width_complete
  summary$coverage_change <- summary$coverage-summary$coverage_complete

  list(summary=summary, estimates=estimates, design=design,
       missing=missing, seed=seed)
}

replicates <- as.integer(Sys.getenv("INGARCH_MISSING_SIM_REPLICATES", "1000"))
simulation_time <- system.time(simulation <- run_missingness_study(
  Delta=c(.5, .8), n=1000, missing=c(0, .05, .10, .20),
  replicates=replicates
))
output <- paste0("INGARCH_simulation_missing_", replicates)
saveRDS(simulation, paste0(output, ".rds"))
write.csv(simulation$summary, paste0(output, ".csv"), row.names=FALSE)
write.csv(simulation$estimates, paste0(output, "_estimates.csv"), row.names=FALSE)
print(simulation$summary)
print(simulation_time)
