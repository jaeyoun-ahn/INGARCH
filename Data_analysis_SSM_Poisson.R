source("aux_data_manage.R")
library(nimble)

###############
## Data prep
###############

Y_mat    <- Y_train
Y_mat_na <- Y_mat
Y_mat_na[is.na(Y_mat_na)] <- 0L

nn <- nrow(Y_mat)
tt <- ncol(Y_mat)
kk <- ncol(XX_train)

zero_mean <- rep(0, kk)
diag_cov  <- diag(rep(10.0, kk)) + 5

constants <- list(
  nn        = nn,
  tt        = tt,
  kk        = kk,
  X         = XX_train,
  zero_mean = zero_mean,
  diag_cov  = diag_cov,
  M         = IDXint,
  num_obs   = num_obs1
)

data_list <- list(
  Y = Y_mat_na
)

#######################
## Model
#######################

code_ss_pois <- nimbleCode({
  
  ## AR(1) coefficient
  phi ~ dunif(0.001, 0.999)
  
  ## innovation sd for latent state
  sigma_state ~ dgamma(shape=0.05, rate=0.05)
  tau_state <- 1 / (sigma_state * sigma_state)
  
  ## stationary precision for S[,1]
  tau_init <- (1 - phi * phi) / (sigma_state * sigma_state)
  
  for (i in 1:nn) {
    
    ## latent Gaussian AR(1) state
    S[i,1] ~ dnorm(0, tau_init)
    
    for (t in 1:(tt-1)) {
      S[i,t+1] ~ dnorm(phi * S[i,t], tau_state)
    }
    
    ## linear predictor
    for (t in 1:tt) {
      eta[i,t] <- inprod(beta[1:kk], X[i,1:kk]) + S[i,t]
      lambda[i,t] <- exp(eta[i,t])
    }
    
    ## observed likelihood only at observed times
    for (j in 1:num_obs[i]) {
      Y[i, M[i,j]] ~ dpois(lambda[i, M[i,j]])
    }
    
    ## one-step-ahead prediction
    S_fore[i] ~ dnorm(phi * S[i,tt], tau_state)
    eta_fore[i] <- inprod(beta[1:kk], X[i,1:kk]) + S_fore[i]
    lambda_fore[i] <- exp(eta_fore[i])
    y6[i] ~ dpois(lambda_fore[i])
  }
  
  ## prior for beta
  beta[1:kk] ~ dmnorm(zero_mean[1:kk], diag_cov[1:kk,1:kk])
})

#######################
## Initial values
#######################

make_inits_ss_pois <- function() {
  list(
    beta        = rnorm(kk, 0, 0.1),
    phi         = runif(1, 0.2, 0.9),
    sigma_state = runif(1, 0.2, 1.0),
    S           = matrix(rnorm(nn * tt, 0, 0.3), nrow = nn, ncol = tt),
    S_fore      = rnorm(nn, 0, 0.3),
    y6          = rep(0L, nn)
  )
}

#######################
## Run MCMC
#######################

set.seed(101)

samples_model_ss_pois <- nimbleMCMC(
  code      = code_ss_pois,
  constants = constants,
  data      = data_list,
  inits     = list(
    make_inits_ss_pois(),
    make_inits_ss_pois(),
    make_inits_ss_pois()
  ),
  monitors  = c("beta", "phi", "sigma_state", "S", "S_fore", "y6"),
  niter     = 50000,
  nburnin   = 10000,
  thin      = 5,
  nchains   = 3,
  WAIC      = TRUE
)

library(MCMCvis)
summary_mcmc_ss_pois <- MCMCsummary(
  object = samples_model_ss_pois$samples,
  round  = 5,
  params = c("beta", "phi", "sigma_state")
)
print(summary_mcmc_ss_pois)


library(coda)

samps_all_ss <- rbind(
  samples_model_ss_pois$samples$chain1,
  samples_model_ss_pois$samples$chain2,
  samples_model_ss_pois$samples$chain3
)

n_samp <- nrow(samps_all_ss)

idx_beta <- grep("^beta\\[", colnames(samps_all_ss))
beta_samps <- samps_all_ss[, idx_beta, drop = FALSE]

phi_samps   <- samps_all_ss[, "phi"]
sig_samps   <- samps_all_ss[, "sigma_state"]

nn_test <- length(Y_test)
log_pred_i_ss_pois <- numeric(nn_test)

test_list <- which(ID_train %in% ID_test)

for (i in 1:nn_test) {
  
  y_obs_i <- Y_test[i]
  logw <- numeric(n_samp)
  
  for (m in 1:n_samp) {
    cc_i <- sum(beta_samps[m, ] * XX_test[i, ])
    
    s_name <- sprintf("S[%d, %d]", test_list[i], tt)
    s_last <- samps_all_ss[m, s_name]
    
    mu_fore <- phi_samps[m] * s_last
    s_fore  <- rnorm(1, mean = mu_fore, sd = sig_samps[m])
    
    lambda_fore <- exp(cc_i + s_fore)
    
    logw[m] <- dpois(y_obs_i, lambda = lambda_fore, log = TRUE)
  }
  
  maxlog <- max(logw)
  log_pred_i_ss_pois[i] <- maxlog + log(mean(exp(logw - maxlog)))
}

test_loglik_postpred_ss_pois <- sum(log_pred_i_ss_pois)
test_loglik_postpred_ss_pois
# [1] -217.9274

## MSE
y6_idx <- grep("^y6\\[", colnames(samps_all_ss))
y6_samples <- samps_all_ss[, y6_idx, drop = FALSE]
y6_mean <- colMeans(y6_samples)

mean((y6_mean[ID_train %in% ID_test] - Y_test)^2)
# [1] 1.442389

## MAE
mean(abs(y6_mean[ID_train %in% ID_test] - Y_test))
# [1] 0.6403126

## WAIC
samples_model_ss_pois$WAIC

# nimbleList object of type waicNimbleList
# Field "WAIC":
#   [1] 2657.118
# Field "lppd":
#   [1] -1135.999
# Field "pWAIC":
#   [1] 192.5601

