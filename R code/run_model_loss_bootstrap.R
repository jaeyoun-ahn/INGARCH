source("Data_analysis_INGARCH2.R")

B <- as.integer(Sys.getenv("MODEL_LOSS_BOOTSTRAP_REPLICATES", "1000"))
workers <- 11
seed <- 20260910

type_columns <- grep("^Type", dimnames(XX_train)[[3]])
type_design <- XX_train[, 1, type_columns, drop=FALSE]
dim(type_design) <- c(nrow(Y_train), length(type_columns))
strata <- max.col(cbind(1-rowSums(type_design), type_design))
test_rows <- match(ID_test, ID_train)
training_only <- setdiff(seq_len(nrow(Y_train)), test_rows)

sample_strata <- function(rows)
  unlist(lapply(split(rows, strata[rows]), function(x)
    sample(x, length(x), replace=TRUE)), use.names=FALSE)

set.seed(seed)
samples <- replicate(B, list(test=sample_strata(test_rows),
                             training_only=sample_strata(training_only)),
                     simplify=FALSE)

starts <- list(state=fit_state_nb$par, random=fit_random$par,
               xu=fit_xu$par, loglinear=fit_loglinear$par,
               inar=fit_inar$par, ssm=fit_ssm_poisson$par)

fit_bootstrap_models <- function(sample){
  rows <- c(sample$test, sample$training_only)
  test_position <- seq_along(sample$test)
  test_index <- match(ID_train[sample$test], ID_test)
  Y <- Y_train[rows, , drop=FALSE]
  W <- W_train[rows, , drop=FALSE]
  XX <- XX_train[rows, , , drop=FALSE]
  Yt <- Y_test[test_index]
  XXt <- XX_test[test_index, , drop=FALSE]
  predictor <- function(beta)
    matrix(matrix(XX, nrow(Y)*ncol(Y), p) %*% beta, nrow(Y), ncol(Y))

  obj1 <- make_state_nb_objective(TRUE, Y, W, XX, starts$state)
  fit1 <- optim(obj1$par, obj1$fn, obj1$gr, method="L-BFGS-B",
                lower=c(rep(-Inf, p), 0, -Inf),
                upper=c(rep( Inf, p), .999, Inf),
                control=list(maxit=5000, factr=1e7))
  u1 <- unpack_state_nb(fit1$par)
  f1 <- gamma_filter(Y, W, exp(predictor(u1$beta)), u1$delta, u1$shape)
  A1 <- f1$A[test_position, ncol(Y)+1]
  B1 <- f1$B[test_position, ncol(Y)+1]
  pred1 <- A1/B1*exp(drop(XXt %*% u1$beta))
  logp1 <- dnbinom(Yt, size=A1, mu=pred1, log=TRUE)

  obj2 <- make_state_nb_objective(FALSE, Y, W, XX, starts$random)
  fit2 <- optim(obj2$par, obj2$fn, obj2$gr, method="L-BFGS-B",
                control=list(maxit=5000, factr=1e7))
  u2 <- unpack_state_nb(fit2$par, FALSE)
  f2 <- gamma_filter(Y, W, exp(predictor(u2$beta)), 1, u2$shape)
  A2 <- f2$A[test_position, ncol(Y)+1]
  B2 <- f2$B[test_position, ncol(Y)+1]
  pred2 <- A2/B2*exp(drop(XXt %*% u2$beta))
  logp2 <- dnbinom(Yt, size=A2, mu=pred2, log=TRUE)

  obj3 <- make_xu_objective(starts$xu, Y, W)
  fit3 <- optim(obj3$par, obj3$fn, obj3$gr, method="L-BFGS-B",
                control=list(maxit=5000, factr=1e7))
  u3 <- unpack_xu(fit3$par)
  pred3 <- xu_means(Y, W, u3)[test_position, ncol(Y)+1]
  logp3 <- dnbinom(Yt, size=pred3*u3$precision, mu=pred3, log=TRUE)

  fit4 <- fit_loglinear_rtmb(list(starts$loglinear), Y, W, XX)
  u4 <- unpack_loglinear(fit4$par)
  forecast_xb <- rep(NA_real_, nrow(Y))
  forecast_xb[test_position] <- drop(XXt %*% u4$beta)
  pred4 <- exp(loglinear_means(Y, W, u4, forecast_xb, XX)$prediction[
    test_position])
  logp4 <- dpois(Yt, pred4, log=TRUE)

  groups <- make_inar_groups(Y, W)
  obj5 <- make_inar_objective(starts$inar, Y, XX, groups$transitions,
                              groups$pairs, groups$first, groups$nobs)
  fit5 <- optim(obj5$par, obj5$fn, obj5$gr, method="L-BFGS-B",
                control=list(maxit=5000, factr=1e7))
  u5 <- unpack_inar(fit5$par)
  last_y <- vapply(test_position, function(i)
    tail(Y[i, W[i, ]], 1), numeric(1))
  innovation5 <- exp(drop(XXt %*% u5$beta))
  pred5 <- u5$alpha*last_y+innovation5
  logp5 <- vapply(seq_along(Yt), function(i)
    inar_logpmf(Yt[i], last_y[i], innovation5[i], u5$alpha), numeric(1))

  obj6 <- make_ssm_poisson_objective(starts$ssm, Y, W, XX,
                                     gh_nodes, gh_weights)
  fit6 <- optim(obj6$par, obj6$fn, obj6$gr, method="L-BFGS-B",
                control=list(maxit=5000, factr=1e7))
  f6 <- ssm_poisson_filter(fit6$par, Y, W, XX, gh_nodes, gh_weights)
  u6 <- f6$parameters
  last_weights <- f6$weights[test_position, ]
  xb6 <- drop(XXt %*% u6$beta)
  conditional_means <- exp(outer(xb6, u6$phi*f6$states, "+")+u6$sigma^2/2)
  pred6 <- rowSums(last_weights*conditional_means)
  logp6 <- vapply(seq_along(Yt), function(i){
    rates <- exp(xb6[i]+outer(u6$phi*f6$states,
                              u6$sigma*gh_nodes, "+"))
    log(sum(last_weights[i, ]*rowSums(sweep(
      dpois(Yt[i], rates), 2, gh_weights, "*"))))
  }, numeric(1))

  predictions <- list(pred1, pred2, pred3, pred4, pred5, pred6)
  log_probabilities <- list(logp1, logp2, logp3, logp4, logp5, logp6)
  fits <- list(fit1, fit2, fit3, fit4, fit5, fit6)
  data.frame(
    model=names(models),
    test_MSE=vapply(predictions, function(x) mean((Yt-x)^2), numeric(1)),
    test_MAE=vapply(predictions, function(x) mean(abs(Yt-x)), numeric(1)),
    test_NLL=vapply(log_probabilities, function(x) -sum(x), numeric(1)),
    test_mean_NLL=vapply(log_probabilities, function(x) -mean(x), numeric(1)),
    convergence=vapply(fits, `[[`, numeric(1), "convergence")
  )
}

cluster <- parallel::makeCluster(workers)
parallel::clusterExport(
  cluster,
  c("samples", "starts", "Y_train", "W_train", "XX_train", "Y_test",
    "XX_test", "ID_train", "ID_test", "p", "models", "gh_nodes",
    "gh_weights", "fit_bootstrap_models", "rtmb_state_nb_nll",
    "make_state_nb_objective", "unpack_state_nb", "gamma_filter",
    "inv_logit",
    "rtmb_xu_nll", "make_xu_objective", "unpack_xu", "xu_means",
    "rtmb_loglinear_nll", "make_loglinear_objective", "fit_loglinear_rtmb",
    "unpack_loglinear", "loglinear_means", "make_inar_groups",
    "rtmb_inar_nll", "make_inar_objective", "unpack_inar", "inar_logpmf",
    "log_sum_exp", "rtmb_ssm_poisson_nll", "make_ssm_poisson_objective",
    "unpack_ssm_poisson", "ssm_poisson_filter"),
  envir=.GlobalEnv
)
invisible(parallel::clusterEvalQ(cluster, {
  .libPaths(c("Rlib", .libPaths()))
  library(RTMB)
}))
parallel::clusterSetRNGStream(cluster, seed)
bootstrap_time <- system.time(result <- parallel::parLapply(
  cluster, samples, fit_bootstrap_models
))
parallel::stopCluster(cluster)

losses <- do.call(rbind, Map(function(x, b) transform(x, replicate=b),
                             result, seq_along(result)))
losses <- losses[order(losses$replicate, match(losses$model, names(models))), ]
rownames(losses) <- NULL

pairs <- combn(names(models), 2, simplify=FALSE)
metrics <- c("test_MSE", "test_MAE", "test_NLL")
pairwise <- do.call(rbind, lapply(metrics, function(metric){
  wide <- matrix(losses[[metric]], B, length(models), byrow=TRUE,
                 dimnames=list(NULL, names(models)))
  do.call(rbind, lapply(pairs, function(pair)
    data.frame(replicate=seq_len(B), metric=metric,
               model_A=pair[1], model_B=pair[2],
               improvement_A=wide[, pair[2]]-wide[, pair[1]])))
}))

observed <- setNames(lapply(metrics, function(metric)
  setNames(model_comparison[[metric]], model_comparison$model)), metrics)
pairwise_summary <- do.call(rbind, lapply(
  split(pairwise, list(pairwise$metric, pairwise$model_A, pairwise$model_B),
        drop=TRUE),
  function(x){
    difference <- x$improvement_A
    observed_difference <- observed[[x$metric[1]]][x$model_B[1]]-
      observed[[x$metric[1]]][x$model_A[1]]
    data.frame(metric=x$metric[1], model_A=x$model_A[1], model_B=x$model_B[1],
               observed_improvement_A=observed_difference,
               bootstrap_mean=mean(difference), bootstrap_SE=sd(difference),
               CI_lower=quantile(difference, .025),
               CI_upper=quantile(difference, .975),
               probability_A_better=mean(difference>0))
  }
))
rownames(pairwise_summary) <- NULL

convergence <- aggregate(convergence~model, losses,
                         function(x) mean(x==0))
names(convergence)[2] <- "convergence_rate"

output <- paste0("model_loss_bootstrap_", B)
bootstrap <- list(losses=losses, pairwise=pairwise,
                  pairwise_summary=pairwise_summary,
                  convergence=convergence, samples=samples,
                  seed=seed, bootstrap_time=bootstrap_time)
saveRDS(bootstrap, paste0(output, ".rds"))
write.csv(losses, paste0(output, "_losses.csv"), row.names=FALSE)
write.csv(pairwise, paste0(output, "_pairwise.csv"), row.names=FALSE)
write.csv(pairwise_summary, paste0(output, "_summary.csv"), row.names=FALSE)
print(pairwise_summary)
print(convergence)
print(bootstrap_time)
