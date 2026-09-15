.libPaths(c("Rlib", .libPaths()))
library(RTMB)
source("funs_INGARCH.R")

beta_true <- c("(Intercept)"=-3.5, Type2=.3, Type3=-.3, Continuous=1.3)
gamma_true <- .7

sample_observation_patterns <- function(n, missing=0){
  weights <- c(.6, .3, .1)
  scale <- 5*missing/sum((1:3)*weights)
  gaps <- sample(0:3, n, replace=TRUE,
                 prob=c(1-scale*sum(weights), scale*weights))
  W <- matrix(TRUE, n, 5)
  for(i in seq_len(n)){
    times <- switch(gaps[i]+1, integer(0), sample(2:4, 1),
                    {first <- sample(2:3, 1); first:(first+1)}, 2:4)
    W[i, times] <- FALSE
  }
  W
}

simulate_nb_ingarch <- function(Delta, n, missing=0,
                                beta=beta_true, gamma=gamma_true){
  category <- sample(rep(1:3, length.out=n))
  individual_mean <- rnorm(n)
  continuous <- matrix(rnorm(6*n, rep(individual_mean, 6), .125), n, 6)
  XX <- array(1, c(n, 6, length(beta)))
  XX[, , 2] <- category==2
  XX[, , 3] <- category==3
  XX[, , 4] <- continuous
  dimnames(XX)[[3]] <- names(beta)
  lambda <- exp(matrix(matrix(XX, 6*n, length(beta)) %*% beta, n, 6))

  W <- sample_observation_patterns(n, missing)
  Y <- matrix(NA_real_, n, 5)
  A <- B <- rep(gamma, n)
  for(t in 1:5){
    observed <- W[, t]
    Y[observed, t] <- rnbinom(sum(observed), size=A[observed],
                              mu=A[observed]*lambda[observed, t]/B[observed])
    A_post <- ifelse(observed, A+Y[, t], A)
    B_post <- ifelse(observed, B+lambda[, t], B)
    q <- 1/(Delta^2+(1-Delta^2)*B_post/gamma)
    B <- q*B_post
    A <- Delta*q*A_post+(1-Delta)*B
  }
  Y_test <- rnbinom(n, size=A, mu=A*lambda[, 6]/B)
  list(Y_train=Y, Y_test=Y_test, W=W,
       XX_train=XX[, 1:5, , drop=FALSE], XX_test=XX[, 6, ],
       category=category)
}

fit_simulation <- function(data, start=NULL, keep_fit=FALSE){
  np <- dim(data$XX_train)[3]
  if(is.null(start))
    start <- c(log(mean(data$Y_train, na.rm=TRUE)), rep(0, np-1),
               .1, log(gamma_true))
  obj <- make_state_nb_objective(TRUE, data$Y_train, data$W,
                                 data$XX_train, start)
  fit <- optim(obj$par, obj$fn, obj$gr, method="L-BFGS-B",
               lower=c(rep(-15, np), 0, log(.02)),
               upper=c(rep( 15, np), .999, log(50)),
               control=list(maxit=5000, factr=1e7))
  estimate <- c(fit$par[seq_len(np)], 1-fit$par[np+1],
                exp(fit$par[np+2]))
  names(estimate) <- c(dimnames(data$XX_train)[[3]], "delta", "gamma_1_0")
  information <- obj$he(fit$par)*sum(data$W)
  derivative <- c(rep(1, np), -1, estimate[np+2])
  vcov <- diag(derivative) %*% solve(information) %*% diag(derivative)
  SE <- sqrt(ifelse(diag(vcov)>0, diag(vcov), NA_real_))
  result <- data.frame(
    parameter=c(dimnames(data$XX_train)[[3]], "delta", "gamma_1_0"),
    estimate=unname(estimate), SE=SE, convergence=fit$convergence
  )
  if(keep_fit) list(result=result, par=fit$par) else result
}

run_simulation_study <- function(Delta=c(.5, .8), n=c(200, 500, 1000),
                                 missing=0, replicates=500,
                                 workers=11, seed=100){
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
    fit <- fit_simulation(simulate_nb_ingarch(
      setting$Delta, setting$n, missing
    ))
    transform(fit, replicate=setting$replicate,
              n=setting$n, Delta=setting$Delta, missing=missing)
  })
  parallel::stopCluster(cluster)
  estimates <- do.call(rbind, result)
  truth <- c(beta_true, delta=NA_real_, gamma_1_0=gamma_true)
  estimates$truth <- truth[estimates$parameter]
  estimates$truth[estimates$parameter=="delta"] <-
    estimates$Delta[estimates$parameter=="delta"]
  estimates$covered <- abs(estimates$estimate-estimates$truth) <= 1.96*estimates$SE

  summary <- do.call(rbind, lapply(
    split(estimates, list(estimates$n, estimates$Delta, estimates$parameter)),
    function(x){
      centered <- x$estimate-mean(x$estimate)
      spread <- sqrt(mean(centered^2))
      data.frame(n=x$n[1], Delta=x$Delta[1], missing=x$missing[1],
                 parameter=x$parameter[1], truth=x$truth[1],
                 mean=mean(x$estimate), bias=mean(x$estimate-x$truth),
                 variance=var(x$estimate),
                 MSE=mean((x$estimate-x$truth)^2),
                 mean_SE=mean(x$SE, na.rm=TRUE),
                 coverage=mean(x$covered, na.rm=TRUE),
                 valid_Wald_rate=mean(is.finite(x$SE)),
                 skewness=mean(centered^3)/spread^3,
                 excess_kurtosis=mean(centered^4)/spread^4-3,
                 Shapiro_p=shapiro.test(x$estimate)$p.value,
                 convergence_rate=mean(x$convergence==0))
    }
  ))
  rownames(summary) <- NULL
  list(summary=summary, estimates=estimates, design=design, seed=seed)
}

if(sys.nframe()==0){
  replicates <- as.integer(Sys.getenv("INGARCH_SIM_REPLICATES", "1000"))
  simulation_time <- system.time(simulation <- run_simulation_study(
    Delta=c(.5, .8), n=c(200, 500, 1000), missing=0, replicates=replicates
  ))
  output <- paste0("INGARCH_simulation_complete_", replicates)
  saveRDS(simulation, paste0(output, ".rds"))
  write.csv(simulation$summary, paste0(output, ".csv"), row.names=FALSE)
  write.csv(simulation$estimates, paste0(output, "_estimates.csv"), row.names=FALSE)
  print(simulation$summary)
  print(simulation_time)
}
