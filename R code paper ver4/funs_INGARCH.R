train_predictor <- function(beta)
  matrix(matrix(XX_train, nrow(Y_train)*ncol(Y_train), p) %*% beta,
         nrow(Y_train), ncol(Y_train))

test_predictor <- function(beta) drop(XX_test %*% beta)

inv_logit <- function(x) 1/(1+exp(-x))

log_sum_exp <- function(x){
  m <- max(x)
  m + log(sum(exp(x-m)))
}

best_fit <- function(starts, objective, ...){
  fits <- lapply(starts, function(start)
    optim(start, objective, ..., method="L-BFGS-B",
          control=list(maxit=5000, factr=1e7)))
  fits[[which.min(vapply(fits, function(x) x$value, numeric(1)))]]
}

# 1. Heterogeneous state-space NB-INGARCH -------------------------------

unpack_state_nb <- function(par, estimate_delta=TRUE){
  beta <- par[1:p]
  if(estimate_delta){
    delta <- 1-par[p+1]
    shape <- exp(par[p+2])
  } else {
    delta <- 1
    shape <- exp(par[p+1])
  }
  list(beta=beta, delta=delta, shape=shape)
}

gamma_filter <- function(Y, W, lambda, delta, shape){
  A <- B <- matrix(NA_real_, nrow(Y), ncol(Y)+1)
  A[, 1] <- B[, 1] <- shape
  for(t in 1:ncol(Y)){
    A_post <- ifelse(W[, t], A[, t]+Y[, t], A[, t])
    B_post <- ifelse(W[, t], B[, t]+lambda[, t], B[, t])
    q <- 1/(delta^2+(1-delta^2)*B_post/shape)
    B[, t+1] <- q*B_post
    A[, t+1] <- delta*q*A_post+(1-delta)*B[, t+1]
  }
  list(A=A, B=B)
}

state_nb_nll <- function(par, estimate_delta=TRUE){
  u <- unpack_state_nb(par, estimate_delta)
  lambda <- exp(train_predictor(u$beta))
  f <- gamma_filter(Y_train, W_train, lambda, u$delta, u$shape)
  A <- f$A[, 1:ncol(Y_train)]
  B <- f$B[, 1:ncol(Y_train)]
  ll <- ifelse(W_train,
               dnbinom(Y_train, size=A, mu=A*lambda/B, log=TRUE), 0)
  -sum(ll)/sum(W_train)
}

rtmb_state_nb_nll <- function(par, Y, W, XX, estimate_delta=TRUE){
  beta <- par$beta
  delta <- if(estimate_delta) 1-par$rho else 1
  shape <- exp(par$log_shape)
  lambda <- exp(matrix(matrix(XX, nrow(Y)*ncol(Y), length(beta)) %*% beta,
                       nrow(Y), ncol(Y)))
  A <- B <- rep(shape, nrow(Y))
  ll <- 0
  
  for(t in seq_len(ncol(Y))){
    observed <- W[, t]
    mu <- A[observed]*lambda[observed, t]/B[observed]
    ll <- ll + sum(RTMB::dnbinom2(Y[observed, t], mu,
                                  mu+mu^2/A[observed], TRUE))
    A_post <- ifelse(observed, A+Y[, t], A)
    B_post <- ifelse(observed, B+lambda[, t], B)
    q <- 1/(delta^2+(1-delta^2)*B_post/shape)
    B <- q*B_post
    A <- delta*q*A_post+(1-delta)*B
  }
  -ll/sum(W)
}

make_state_nb_objective <- function(estimate_delta=TRUE, Y=Y_train,
                                    W=W_train, XX=XX_train, start=NULL){
  np <- dim(XX)[3]
  if(is.null(start))
    start <- c(log(mean(Y, na.rm=TRUE)), rep(0, np-1),
               if(estimate_delta) .1, 0)
  parameters <- list(beta=start[seq_len(np)], log_shape=start[np+1])
  if(estimate_delta)
    parameters <- list(beta=start[seq_len(np)], rho=start[np+1],
                       log_shape=start[np+2])
  RTMB::MakeADFun(function(par)
    rtmb_state_nb_nll(par, Y, W, XX, estimate_delta), parameters, silent=TRUE)
}

nb_inference <- function(fit, estimate_delta=TRUE){
  u <- unpack_state_nb(fit$par, estimate_delta)
  estimate <- c(u$beta, if(estimate_delta) delta=u$delta,
                gamma_1_0=u$shape)
  derivative <- c(rep(1, p), if(estimate_delta) -1,
                  u$shape)
  obj <- make_state_nb_objective(estimate_delta, start=fit$par)
  information_raw <- obj$he(fit$par)*sum(W_train)
  vcov_raw <- solve(information_raw)
  jacobian <- diag(derivative)
  vcov <- jacobian %*% vcov_raw %*% jacobian
  names(estimate) <- c(dimnames(XX_train)[[3]],
                       if(estimate_delta) "delta", "gamma_1_0")
  dimnames(vcov) <- list(names(estimate), names(estimate))
  list(estimates=data.frame(estimate, standard_error=sqrt(diag(vcov))),
       information=solve(vcov), vcov=vcov,
       information_raw=information_raw, vcov_raw=vcov_raw)
}

# 2. Poisson-Gamma random-effects model ---------------------------------

bootstrap_state_nb <- function(B=125, workers=11, seed=20260909){
  set.seed(seed)
  type_columns <- grep("^Type", dimnames(XX_train)[[3]])
  type_design <- XX_train[, 1, type_columns, drop=FALSE]
  dim(type_design) <- c(nrow(Y_train), length(type_columns))
  strata <- max.col(cbind(1-rowSums(type_design), type_design))
  groups <- split(seq_len(nrow(Y_train)), strata)
  samples <- replicate(B, unlist(lapply(groups, function(x)
    sample(x, length(x), replace=TRUE))), simplify=FALSE)
  cluster <- parallel::makeCluster(workers)
  on.exit(parallel::stopCluster(cluster))
  
  parallel::clusterExport(
    cluster,
    c("Y_train", "W_train", "XX_train", "p", "fit_state_nb", "fit_random",
      "rtmb_state_nb_nll", "make_state_nb_objective"),
    envir=.GlobalEnv
  )
  parallel::clusterEvalQ(cluster, {
    .libPaths(c("Rlib", .libPaths()))
    library(RTMB)
  })
  
  fits <- parallel::parLapply(cluster, samples, function(rows){
    Y <- Y_train[rows, , drop=FALSE]
    W <- W_train[rows, , drop=FALSE]
    XX <- XX_train[rows, , , drop=FALSE]
    
    obj_state <- make_state_nb_objective(TRUE, Y, W, XX, fit_state_nb$par)
    state <- optim(obj_state$par, obj_state$fn, obj_state$gr,
                   method="L-BFGS-B",
                   lower=c(rep(-Inf, p), 0, -Inf),
                   upper=c(rep( Inf, p), .999, Inf),
                   control=list(maxit=5000, factr=1e7))
    
    obj_random <- make_state_nb_objective(FALSE, Y, W, XX, fit_random$par)
    random <- optim(obj_random$par, obj_random$fn, obj_random$gr,
                    method="L-BFGS-B",
                    control=list(maxit=5000, factr=1e7))
    
    list(
      state=c(state$par[seq_len(p)], delta=1-state$par[p+1],
              gamma_1_0=exp(state$par[p+2])),
      random=c(random$par[seq_len(p)],
               gamma_1_0=exp(random$par[p+1])),
      convergence=c(state=state$convergence, random=random$convergence)
    )
  })
  
  parameter_names <- dimnames(XX_train)[[3]]
  state <- do.call(rbind, lapply(fits, `[[`, "state"))
  random <- do.call(rbind, lapply(fits, `[[`, "random"))
  colnames(state) <- c(parameter_names, "delta", "gamma_1_0")
  colnames(random) <- c(parameter_names, "gamma_1_0")
  convergence <- do.call(rbind, lapply(fits, `[[`, "convergence"))
  
  summarize <- function(draws, estimate, model, convergence_name){
    keep <- convergence[, convergence_name] == 0
    data.frame(
      model=model,
      parameter=colnames(draws),
      estimate=estimate,
      bootstrap_SE=apply(draws[keep, , drop=FALSE], 2, sd),
      CI_lower=apply(draws[keep, , drop=FALSE], 2, quantile, .025),
      CI_upper=apply(draws[keep, , drop=FALSE], 2, quantile, .975),
      successful_replicates=sum(keep)
    )
  }
  
  state_estimate <- c(u_state$beta, delta=u_state$delta,
                      gamma_1_0=u_state$shape)
  random_estimate <- c(u_random$beta, gamma_1_0=u_random$shape)
  summary <- rbind(
    summarize(state, state_estimate, "heterogeneous", "state"),
    summarize(random, random_estimate, "random_effects", "random")
  )
  list(summary=summary, state=state, random=random,
       convergence=convergence, samples=samples, strata=table(strata),
       seed=seed)
}

# 3. Xu et al. homogeneous NB-INGARCH(1,1) -----------------------------

unpack_xu <- function(par){
  weights <- exp(c(0, par[2], par[3]))
  weights <- weights/sum(weights)
  list(omega=exp(par[1]), alpha=weights[2], beta=weights[3],
       precision=exp(par[4]))
}

xu_means <- function(Y, W, u){
  M <- matrix(NA_real_, nrow(Y), ncol(Y)+1)
  M[, 1] <- u$omega/(1-u$alpha-u$beta)
  for(t in 1:ncol(Y)){
    M[, t+1] <- ifelse(W[, t],
                       u$omega+u$alpha*Y[, t]+u$beta*M[, t], M[, t])
  }
  M
}

xu_nll <- function(par){
  u <- unpack_xu(par)
  M <- xu_means(Y_train, W_train, u)[, 1:ncol(Y_train)]
  ll <- ifelse(W_train,
               dnbinom(Y_train, size=M*u$precision, mu=M, log=TRUE), 0)
  -sum(ll)/sum(W_train)
}

rtmb_xu_nll <- function(par, Y, W){
  weights <- exp(c(0, par$logit_weights))
  weights <- weights/sum(weights)
  omega <- exp(par$log_omega)
  alpha <- weights[2]
  beta <- weights[3]
  precision <- exp(par$log_precision)
  M <- rep(omega/(1-alpha-beta), nrow(Y))
  ll <- 0
  for(t in seq_len(ncol(Y))){
    observed <- W[, t]
    mu <- M[observed]
    ll <- ll+sum(RTMB::dnbinom2(Y[observed, t], mu,
                                 mu+mu/precision, TRUE))
    M <- ifelse(observed, omega+alpha*Y[, t]+beta*M, M)
  }
  -ll/sum(W)
}

make_xu_objective <- function(start, Y=Y_train, W=W_train){
  parameters <- list(log_omega=start[1], logit_weights=start[2:3],
                     log_precision=start[4])
  RTMB::MakeADFun(function(par) rtmb_xu_nll(par, Y, W),
                  parameters, silent=TRUE)
}

# 4. Log-linear Poisson autoregression ----------------------------------

unpack_loglinear <- function(par){
  raw <- par[(p+1):(p+2)]
  length_raw <- sqrt(sum(raw^2))
  if(length_raw==0) return(list(beta=par[1:p], a=0, b=0))
  direction <- raw/length_raw
  boundary <- if(direction[1]*direction[2]>=0)
    1/sum(abs(direction)) else 1
  ab <- boundary*tanh(length_raw)*direction
  list(beta=par[1:p], a=ab[1], b=ab[2])
}

loglinear_means <- function(Y, W, u, forecast_xb=NULL, XX=XX_train){
  L <- matrix(NA_real_, nrow(Y), ncol(Y))
  prediction <- rep(NA_real_, nrow(Y))
  xb <- matrix(matrix(XX, nrow(Y)*ncol(Y), length(u$beta)) %*% u$beta,
               nrow(Y), ncol(Y))
  for(i in 1:nrow(Y)){
    observed <- which(W[i, ])
    L[i, observed[1]] <- xb[i, observed[1]]
    if(length(observed)>1){
      for(j in 2:length(observed)){
        current <- observed[j]
        previous <- observed[j-1]
        L[i, current] <- xb[i, current]+u$a*L[i, previous]+u$b*log1p(Y[i, previous])
      }
    }
    if(!is.null(forecast_xb)){
      last <- tail(observed, 1)
      prediction[i] <- forecast_xb[i]+u$a*L[i, last]+u$b*log1p(Y[i, last])
    }
  }
  list(training=L, prediction=prediction)
}

loglinear_nll <- function(par){
  L <- loglinear_means(Y_train, W_train, unpack_loglinear(par))$training
  ll <- ifelse(W_train,
               dpois(Y_train, exp(L), log=TRUE), 0)
  -sum(ll)/sum(W_train)
}

rtmb_loglinear_nll <- function(par, Y, W, XX, previous, same_sign){
  raw <- par$raw
  raw_length <- sqrt(sum(raw^2)+1e-20)
  direction <- raw/raw_length
  boundary <- if(same_sign) 1/sum(abs(direction)) else 1
  ab <- boundary*tanh(raw_length)*direction
  xb <- matrix(matrix(XX, nrow(Y)*ncol(Y), length(par$beta)) %*% par$beta,
               nrow(Y), ncol(Y))
  eta <- matrix(0, nrow(Y), ncol(Y))
  ll <- 0
  for(t in seq_len(ncol(Y))){
    observed <- which(W[, t])
    first <- observed[previous[observed, t]==0]
    later <- observed[previous[observed, t]>0]
    if(length(first)) eta[first, t] <- xb[first, t]
    if(length(later)){
      index <- cbind(later, previous[later, t])
      eta[later, t] <- xb[later, t]+ab[1]*eta[index]+ab[2]*log1p(Y[index])
    }
    ll <- ll+sum(RTMB::dpois(Y[observed, t], exp(eta[observed, t]), TRUE))
  }
  -ll/sum(W)
}

make_loglinear_objective <- function(start, Y=Y_train, W=W_train,
                                     XX=XX_train){
  np <- dim(XX)[3]
  same_sign <- start[np+1]*start[np+2]>=0
  previous <- matrix(0L, nrow(Y), ncol(Y))
  for(i in seq_len(nrow(Y))){
    observed <- which(W[i, ])
    if(length(observed)>1)
      previous[cbind(i, tail(observed, -1))] <- head(observed, -1)
  }
  parameters <- list(beta=start[seq_len(np)], raw=start[np+(1:2)])
  RTMB::MakeADFun(function(par)
                    rtmb_loglinear_nll(par, Y, W, XX, previous, same_sign),
                  parameters, silent=TRUE)
}

fit_loglinear_rtmb <- function(starts, Y=Y_train, W=W_train, XX=XX_train){
  np <- dim(XX)[3]
  fits <- lapply(starts, function(start){
    obj <- make_loglinear_objective(start, Y, W, XX)
    positive <- start[np+(1:2)]>0
    lower <- c(rep(-Inf, np), ifelse(positive, 0, -Inf))
    upper <- c(rep( Inf, np), ifelse(positive, Inf, 0))
    optim(obj$par, obj$fn, obj$gr, method="L-BFGS-B",
          lower=lower, upper=upper,
          control=list(maxit=5000, factr=1e7))
  })
  fits[[which.min(vapply(fits, function(x) x$value, numeric(1)))]]
}

loglinear_start <- function(a, b){
  direction <- c(a, b)/sqrt(a^2+b^2)
  boundary <- if(a*b>=0) 1/sum(abs(direction)) else 1
  raw <- atanh(sqrt(a^2+b^2)/boundary)*direction
  intercept <- (1-a)*log(mean_y)-b*log1p(mean_y)
  c(intercept, rep(0, p-1), raw)
}

# 5. Log-linked Poisson INAR(1) -----------------------------------------

make_inar_groups <- function(Y, W){
  transitions <- do.call(rbind, lapply(seq_len(nrow(Y)), function(i){
    observed <- which(W[i, ])
    if(length(observed)<2) return(NULL)
    data.frame(policy=i, time=tail(observed, -1),
               previous=Y[i, head(observed, -1)],
               y=Y[i, tail(observed, -1)])
  }))
  list(
    transitions=transitions,
    pairs=split(seq_len(nrow(transitions)),
                paste(transitions$previous, transitions$y)),
    first=cbind(seq_len(nrow(Y)), apply(W, 1, which.max)),
    nobs=nrow(transitions)+nrow(Y)
  )
}

unpack_inar <- function(par){
  list(beta=par[1:p], alpha=inv_logit(par[p+1]))
}

inar_logpmf <- function(y, previous_y, innovation_mean, alpha){
  survivors <- 0:min(y, previous_y)
  log_sum_exp(
    dbinom(survivors, previous_y, alpha, log=TRUE) +
      dpois(y-survivors, innovation_mean, log=TRUE)
  )
}

inar_nll <- function(par){
  u <- unpack_inar(par)
  mu <- exp(train_predictor(u$beta))
  ll <- sum(dpois(Y_train[inar_first], mu[inar_first], log=TRUE))
  for(rows in inar_pairs){
    previous <- inar_transitions$previous[rows[1]]
    y <- inar_transitions$y[rows[1]]
    survivors <- 0:min(previous, y)
    log_terms <- matrix(
      dpois(matrix(y-survivors, nrow=length(rows), ncol=length(survivors),
                   byrow=TRUE), mu[cbind(inar_transitions$policy[rows],
                                         inar_transitions$time[rows])], log=TRUE),
      nrow=length(rows)
    )
    log_terms <- sweep(log_terms, 2,
                       dbinom(survivors, previous, u$alpha, log=TRUE), "+")
    largest <- apply(log_terms, 1, max)
    ll <- ll+sum(largest+log(rowSums(exp(log_terms-largest))))
  }
  -ll/inar_nobs
}

rtmb_inar_nll <- function(par, Y, XX, transitions, pairs, first, nobs){
  alpha <- 1/(1+exp(-par$logit_alpha))
  mu <- exp(matrix(matrix(XX, nrow(Y)*ncol(Y), length(par$beta)) %*%
                     par$beta, nrow(Y), ncol(Y)))
  ll <- sum(RTMB::dpois(Y[first], mu[first], TRUE))
  for(rows in pairs){
    previous <- transitions$previous[rows[1]]
    y <- transitions$y[rows[1]]
    innovation <- mu[cbind(transitions$policy[rows],
                           transitions$time[rows])]
    logp <- RTMB::dbinom(0, previous, alpha, TRUE)+
      RTMB::dpois(y, innovation, TRUE)
    if(min(previous, y)>0) for(survivors in seq_len(min(previous, y))){
      term <- RTMB::dbinom(survivors, previous, alpha, TRUE)+
        RTMB::dpois(y-survivors, innovation, TRUE)
      logp <- RTMB::logspace_add(logp, term)
    }
    ll <- ll+sum(logp)
  }
  -ll/nobs
}

make_inar_objective <- function(start, Y=Y_train, XX=XX_train,
                                transitions=inar_transitions,
                                pairs=inar_pairs, first=inar_first,
                                nobs=inar_nobs){
  np <- dim(XX)[3]
  parameters <- list(beta=start[seq_len(np)], logit_alpha=start[np+1])
  RTMB::MakeADFun(function(par)
    rtmb_inar_nll(par, Y, XX, transitions, pairs, first, nobs),
    parameters, silent=TRUE)
}


# 6. Poisson Gaussian AR(1) state-space model --------------------------

unpack_ssm_poisson <- function(par){
  list(beta=par[1:p], phi=.001+.998*inv_logit(par[p+1]),
       sigma=exp(par[p+2]))
}

ssm_poisson_filter <- function(par, Y=Y_train, W=W_train, XX=XX_train,
                               nodes=gh_nodes, weights=gh_weights){
  u <- unpack_ssm_poisson(par)
  state_sd <- u$sigma/sqrt(1-u$phi^2)
  states <- state_sd*nodes
  transition <- outer(states, states,
                      function(previous, current)
                        dnorm(current, u$phi*previous, u$sigma)/
                        dnorm(current, 0, state_sd))
  transition <- sweep(transition, 2, weights, "*")
  transition <- transition/rowSums(transition)
  xb <- matrix(matrix(XX, nrow(Y)*ncol(Y), length(u$beta)) %*% u$beta,
               nrow(Y), ncol(Y))
  alpha <- matrix(weights, nrow(Y), length(states), byrow=TRUE)
  loglik <- 0
  for(t in 1:ncol(Y)){
    if(t>1) alpha <- alpha %*% transition
    observed <- which(W[, t])
    alpha[observed, ] <- alpha[observed, ]*
      dpois(Y[observed, t], exp(outer(xb[observed, t], states, "+")))
    scale <- rowSums(alpha)
    if(any(!is.finite(scale)) || any(scale<=0)) return(NULL)
    loglik <- loglik+sum(log(scale))
    alpha <- alpha/scale
  }
  list(loglik=loglik, weights=alpha, states=states,
       transition=transition, parameters=u)
}

ssm_poisson_nll <- function(par){
  f <- ssm_poisson_filter(par)
  if(is.null(f)) return(1e100)
  -f$loglik/sum(W_train)
}

rtmb_ssm_poisson_nll <- function(par, Y, W, XX, nodes, weights){
  phi <- .001+.998/(1+exp(-par$logit_phi))
  sigma <- exp(par$log_sigma)
  state_sd <- sigma/sqrt(1-phi^2)
  states <- state_sd*nodes
  K <- length(states)
  previous <- matrix(nodes, K, K)
  current <- matrix(nodes, K, K, byrow=TRUE)
  transition <- exp(-.5*log(1-phi^2)+.5*current^2-
                      .5*(current-phi*previous)^2/(1-phi^2))
  transition <- sweep(transition, 2, weights, "*")
  transition <- transition/rowSums(transition)
  xb <- matrix(matrix(XX, nrow(Y)*ncol(Y), length(par$beta)) %*% par$beta,
               nrow(Y), ncol(Y))
  alpha <- matrix(weights, nrow(Y), K, byrow=TRUE)
  ll <- 0
  for(t in seq_len(ncol(Y))){
    if(t>1) alpha <- alpha %*% transition
    observed <- which(W[, t])
    rate <- exp(outer(xb[observed, t], states, "+"))
    alpha[observed, ] <- alpha[observed, ]*
      RTMB::dpois(Y[observed, t], rate)
    scale <- rowSums(alpha)
    ll <- ll+sum(log(scale))
    alpha <- alpha/scale
  }
  -ll/sum(W)
}

make_ssm_poisson_objective <- function(start, Y=Y_train, W=W_train,
                                       XX=XX_train, nodes=gh_nodes,
                                       weights=gh_weights){
  np <- dim(XX)[3]
  parameters <- list(beta=start[seq_len(np)], logit_phi=start[np+1],
                     log_sigma=start[np+2])
  RTMB::MakeADFun(function(par)
    rtmb_ssm_poisson_nll(par, Y, W, XX, nodes, weights),
    parameters, silent=TRUE)
}
