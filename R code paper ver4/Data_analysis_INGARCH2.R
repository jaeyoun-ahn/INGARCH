source("aux_data_manage4.R")

.libPaths(c("Rlib", .libPaths()))
library(RTMB)

W_train <- !is.na(Y_train)
test_rows <- match(ID_test, ID_train)
p <- dim(XX_train)[3]
mean_y <- mean(Y_train, na.rm=TRUE)

source("funs_INGARCH.R")

inar_transitions <- do.call(rbind, lapply(seq_len(nrow(Y_train)), function(i){
  observed <- which(W_train[i, ])
  if(length(observed)<2) return(NULL)
  data.frame(policy=i, time=tail(observed, -1),
             previous=Y_train[i, head(observed, -1)],
             y=Y_train[i, tail(observed, -1)])
}))
inar_pairs <- split(seq_len(nrow(inar_transitions)),
                    paste(inar_transitions$previous, inar_transitions$y))
inar_first <- cbind(seq_len(nrow(Y_train)), apply(W_train, 1, which.max))
inar_nobs <- nrow(inar_transitions)+nrow(inar_first)

gh_matrix <- matrix(0, 11, 11)
gh_matrix[row(gh_matrix)==col(gh_matrix)+1] <- sqrt((1:10)/2)
gh <- eigen(gh_matrix+t(gh_matrix), symmetric=TRUE)
gh_weights <- gh$vectors[1, ]^2
gh_nodes <- sqrt(2)*gh$values

# 1. Heterogeneous state-space NB-INGARCH -------------------------------


time1 <- system.time({
  obj_state_nb <- make_state_nb_objective()
  fit_state_nb <- optim(obj_state_nb$par, obj_state_nb$fn, obj_state_nb$gr,
                      method="L-BFGS-B",
                      lower=c(rep(-Inf, p), 0, -Inf),
                      upper=c(rep( Inf, p), .999, Inf),
                      control=list(maxit=5000, factr=1e7))
})
u_state <- unpack_state_nb(fit_state_nb$par)
lambda_state <- exp(train_predictor(u_state$beta))
f_state <- gamma_filter(
  Y_train, W_train, lambda_state,
  u_state$delta, u_state$shape
)
A_state <- f_state$A[test_rows, ncol(Y_train)+1]
B_state <- f_state$B[test_rows, ncol(Y_train)+1]
pred_state <- A_state/B_state*exp(test_predictor(u_state$beta))
logp_state <- dnbinom(Y_test, size=A_state, mu=pred_state, log=TRUE)


# 2. Poisson-Gamma random-effects model ---------------------------------

time2 <- system.time({
  obj_random <- make_state_nb_objective(FALSE)
  fit_random <- optim(obj_random$par, obj_random$fn, obj_random$gr,
                    method="L-BFGS-B",
                    control=list(maxit=5000, factr=1e7))
})
u_random <- unpack_state_nb(fit_random$par, estimate_delta=FALSE)
lambda_random <- exp(train_predictor(u_random$beta))
f_random <- gamma_filter(
  Y_train, W_train, lambda_random,
  1, u_random$shape
)
A_random <- f_random$A[test_rows, ncol(Y_train)+1]
B_random <- f_random$B[test_rows, ncol(Y_train)+1]
pred_random <- A_random/B_random*exp(test_predictor(u_random$beta))
logp_random <- dnbinom(Y_test, size=A_random, mu=pred_random, log=TRUE)

inference_state_nb <- nb_inference(fit_state_nb)
inference_random <- nb_inference(fit_random, estimate_delta=FALSE)


# 3. Xu et al. homogeneous NB-INGARCH(1,1) -----------------------------


alpha0 <- 0.1
beta0 <- 0.7
rest0 <- 1-alpha0-beta0
xu_start <- c(log(mean_y*rest0), log(alpha0/rest0),
              log(beta0/rest0), 0)
time3 <- system.time({
  obj_xu <- make_xu_objective(xu_start)
  fit_xu <- optim(obj_xu$par, obj_xu$fn, obj_xu$gr, method="L-BFGS-B",
                  control=list(maxit=5000, factr=1e7))
})
u_xu <- unpack_xu(fit_xu$par)
pred_xu <- xu_means(Y_train, W_train, u_xu)[test_rows, ncol(Y_train)+1]
logp_xu <- dnbinom(Y_test, size=pred_xu*u_xu$precision,
                   mu=pred_xu, log=TRUE)


# 4. Log-linear Poisson autoregression ----------------------------------


ab_starts <- rbind(c(.5, .2), c(.2, .5), c(-.2, .4),
                   c(.4, -.2), c(.05, .05))
time4 <- system.time(
  fit_loglinear <- fit_loglinear_rtmb(
  lapply(1:nrow(ab_starts), function(i)
    loglinear_start(ab_starts[i, 1], ab_starts[i, 2]))
)
)
u_loglinear <- unpack_loglinear(fit_loglinear$par)
forecast_xb <- rep(NA_real_, nrow(Y_train))
forecast_xb[test_rows] <- test_predictor(u_loglinear$beta)
L_loglinear <- loglinear_means(
  Y_train, W_train, u_loglinear, forecast_xb
)
pred_loglinear <- exp(L_loglinear$prediction[test_rows])
logp_loglinear <- dpois(Y_test, pred_loglinear, log=TRUE)


# 5. Log-linked Poisson INAR(1) -----------------------------------------


inar_starts <- lapply(c(.05, .25, .50), function(alpha)
  c(log(mean_y*(1-alpha)), rep(0, p-1), qlogis(alpha)))
time5 <- system.time({
  obj_inar <- make_inar_objective(inar_starts[[1]])
  fit_inar <- best_fit(inar_starts, obj_inar$fn, gr=obj_inar$gr)
})
u_inar <- unpack_inar(fit_inar$par)
last_y <- vapply(test_rows, function(i)
  tail(Y_train[i, W_train[i, ]], 1), numeric(1))
mu_inar_test <- exp(test_predictor(u_inar$beta))
pred_inar <- u_inar$alpha*last_y+mu_inar_test
logp_inar <- vapply(seq_along(Y_test), function(i)
  inar_logpmf(Y_test[i], last_y[i], mu_inar_test[i], u_inar$alpha),
  numeric(1))


# 6. Poisson Gaussian AR(1) state-space model --------------------------

ssm_starts <- list(c(u_loglinear$beta, qlogis((.8-.001)/.998), log(.3)))
time6 <- system.time({
  obj_ssm_poisson <- make_ssm_poisson_objective(ssm_starts[[1]])
  fit_ssm_poisson <- best_fit(ssm_starts, obj_ssm_poisson$fn,
                              gr=obj_ssm_poisson$gr)
})
f_ssm_poisson <- ssm_poisson_filter(fit_ssm_poisson$par)
last_weights <- f_ssm_poisson$weights[test_rows, ]
u_ssm_poisson <- f_ssm_poisson$parameters
xb_test <- test_predictor(u_ssm_poisson$beta)
conditional_means <- exp(outer(xb_test, u_ssm_poisson$phi*f_ssm_poisson$states,
                               "+") + u_ssm_poisson$sigma^2/2)
pred_ssm_poisson <- rowSums(last_weights*conditional_means)
logp_ssm_poisson <- vapply(seq_along(Y_test), function(i){
  rates <- exp(xb_test[i]+outer(u_ssm_poisson$phi*f_ssm_poisson$states,
                                u_ssm_poisson$sigma*gh_nodes, "+"))
  log(sum(last_weights[i, ]*rowSums(sweep(
    dpois(Y_test[i], rates), 2, gh_weights, "*"))))
}, numeric(1))


# Out-of-sample comparison -----------------------------------------------



models <- list(
  "Heterogeneous NB-INGARCH" =
    list(fit=fit_state_nb, prediction=pred_state, logp=logp_state),
  "Poisson-Gamma random effects" =
    list(fit=fit_random, prediction=pred_random, logp=logp_random),
  "Xu et al. homogeneous NB-INGARCH" =
    list(fit=fit_xu, prediction=pred_xu, logp=logp_xu),
  "Log-linear Poisson autoregression with covariates" =
    list(fit=fit_loglinear, prediction=pred_loglinear, logp=logp_loglinear),
  "Log-linked Poisson INAR(1) with covariates" =
    list(fit=fit_inar, prediction=pred_inar, logp=logp_inar),
  "Poisson Gaussian AR(1) state-space" =
    list(fit=fit_ssm_poisson, prediction=pred_ssm_poisson,
         logp=logp_ssm_poisson)
)

time_elapsed <- c(time1['elapsed'], time2['elapsed'], 
                  time3['elapsed'], time4['elapsed'], 
                  time5['elapsed'], time6['elapsed'])

model_comparison <- data.frame(
  model=names(models),
  n_parameters=vapply(models, function(x) length(x$fit$par), numeric(1)),
  test_MSE=vapply(models, function(x) mean((Y_test-x$prediction)^2), numeric(1)),
  test_MAE=vapply(models, function(x) mean(abs(Y_test-x$prediction)), numeric(1)),
  test_NLL=vapply(models, function(x) -sum(x$logp), numeric(1)),
  time_elapsed = time_elapsed, 
  convergence=vapply(models, function(x) x$fit$convergence, numeric(1)),
  row.names=NULL
)
model_comparison$MSE_rank <- rank(model_comparison$test_MSE)
model_comparison$MAE_rank <- rank(model_comparison$test_MAE)
model_comparison$NLL_rank <- rank(model_comparison$test_NLL)
model_comparison$average_improvement <-
  (model_comparison$test_NLL-model_comparison$test_NLL[1])/length(Y_test)
model_comparison$average_improvement[1] <- NA


print(model_comparison)
print(inference_state_nb$estimates)
print(inference_random$estimates)
write.csv(model_comparison, "model_comparison_varying2.csv")
