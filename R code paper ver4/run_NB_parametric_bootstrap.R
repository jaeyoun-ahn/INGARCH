source("aux_data_manage4.R")

.libPaths(c("Rlib", .libPaths()))
library(RTMB)

W_train <- !is.na(Y_train)
p <- dim(XX_train)[3]
mean_y <- mean(Y_train, na.rm=TRUE)
source("funs_INGARCH.R")

fit_null_model <- function(Y, W, XX, start){
  obj <- make_state_nb_objective(FALSE, Y, W, XX, start)
  optim(obj$par, obj$fn, obj$gr, method="L-BFGS-B",
        control=list(maxit=5000, factr=1e7))
}

fit_alternative_model <- function(Y, W, XX, null_fit){
  start <- c(null_fit$par[seq_len(p)], 0, null_fit$par[p+1])
  obj <- make_state_nb_objective(TRUE, Y, W, XX, start)
  fits <- lapply(c(0, .05), function(rho){
    start[p+1] <- rho
    optim(start, obj$fn, obj$gr, method="L-BFGS-B",
          lower=c(rep(-Inf, p), 0, -Inf),
          upper=c(rep( Inf, p), .999, Inf),
          control=list(maxit=5000, factr=1e7))
  })
  fits[[which.min(vapply(fits, function(x) x$value, numeric(1)))]]
}

null_fit <- fit_null_model(
  Y_train, W_train, XX_train,
  c(log(mean_y), rep(0, p-1), 0)
)
alternative_fit <- fit_alternative_model(
  Y_train, W_train, XX_train, null_fit
)

n_observed <- sum(W_train)
observed <- data.frame(
  LR=max(0, 2*n_observed*(null_fit$value-alternative_fit$value)),
  delta=1-alternative_fit$par[p+1]
)

B <- as.integer(Sys.getenv("NB_PARAMETRIC_BOOTSTRAP_REPLICATES", "1000"))
workers <- min(11, B)
seed <- 20260910

cluster <- parallel::makeCluster(workers)
parallel::clusterExport(
  cluster,
  c("Y_train", "W_train", "XX_train", "p", "null_fit",
    "rtmb_state_nb_nll", "make_state_nb_objective",
    "fit_null_model", "fit_alternative_model"),
  envir=.GlobalEnv
)
invisible(parallel::clusterEvalQ(cluster, {
  .libPaths(c("Rlib", .libPaths()))
  library(RTMB)
}))
parallel::clusterSetRNGStream(cluster, seed)

bootstrap_time <- system.time(
  bootstrap <- parallel::parLapply(cluster, seq_len(B), function(b){
    beta <- null_fit$par[seq_len(p)]
    gamma <- exp(null_fit$par[p+1])
    lambda <- exp(matrix(
      matrix(XX_train, nrow(Y_train)*ncol(Y_train), p) %*% beta,
      nrow(Y_train), ncol(Y_train)
    ))
    random_effect <- rgamma(nrow(Y_train), gamma, gamma)
    Y <- matrix(NA_real_, nrow(Y_train), ncol(Y_train))
    Y[W_train] <- rpois(sum(W_train), (random_effect*lambda)[W_train])

    fit0 <- fit_null_model(Y, W_train, XX_train, null_fit$par)
    fit1 <- fit_alternative_model(Y, W_train, XX_train, fit0)
    data.frame(
      replicate=b,
      LR=max(0, 2*sum(W_train)*(fit0$value-fit1$value)),
      delta=1-fit1$par[p+1],
      null_NLL=fit0$value,
      alternative_NLL=fit1$value,
      null_convergence=fit0$convergence,
      alternative_convergence=fit1$convergence
    )
  })
)
parallel::stopCluster(cluster)

bootstrap <- do.call(rbind, bootstrap)
p_value <- data.frame(
  LR=(1+sum(bootstrap$LR >= observed$LR))/(B+1),
  delta=(1+sum(bootstrap$delta <= observed$delta))/(B+1)
)

results <- list(observed=observed, bootstrap=bootstrap,
                p_value=p_value, null_fit=null_fit,
                alternative_fit=alternative_fit, seed=seed,
                bootstrap_time=bootstrap_time)
output <- paste0("NB_INGARCH_parametric_bootstrap_", B)
saveRDS(results, paste0(output, ".rds"))
write.csv(bootstrap, paste0(output, ".csv"), row.names=FALSE)

print(observed)
print(bootstrap)
print(p_value)
print(bootstrap_time)


dat <- bootstrap
#QQ-plot
myqvalues <- dat$LR[dat$delta<1]
chisqq <- qchisq(seq(1, length(myqvalues), 1)/(length(myqvalues)+1), df=1)
plot(sort(myqvalues), chisqq)
abline(a=0, b=1)
names(dat)
pchisq(2.3050, df=1, lower.tail = FALSE)/2
(1+sum(dat$LR>2.3050))/(nrow(dat)+1)
