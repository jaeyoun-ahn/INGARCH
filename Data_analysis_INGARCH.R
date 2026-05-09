source("aux_data_manage.R")

W_train <- 1*!is.na(Y_train)

test_ind <- which(ID_train %in% ID_test)


#utilities
inv_logit <- function(x) 1/(1+exp(-x))


unpack1 <- function(params){
  beta <- params[1:p]
  delta <- inv_logit(params[p+1])
  mu <- exp(params[p+2])
  list(beta=beta, delta=delta, mu=mu)
}


#calculates mu's
#once mu's are calculated, dpois can be evaluated on 
#desired mu's


mu_recursion <- function(z, w, probs, Delta, mu, 
                         stationary =TRUE){
  n_pol <- nrow(z)
  n_col <- ncol(z)
  mu_mat <- matrix(NA, nrow = nrow(z), ncol = ncol(z)+1)
  mu_mat[, 1] <- mu
  for(j in 1:(n_col)){
    
    
    currentmu <- mu_mat[, j]
    currentz <- z[, j]
    currentw <- w[, j]
    if(stationary){
      currentDelta <- ifelse(currentw, 
                             sqrt(Delta^2*probs[, 1]/(Delta^2*probs[, 1]+probs[, j])),
                             1)
    } else {
      currentDelta <- rep(Delta, n_pol)
    }
    currentprobs <- probs[, j]
    currentlam <- currentmu * currentw * currentprobs
    currentmu_filt <- ifelse(currentw==0, 
                             currentmu, 
                             currentz + (1-currentprobs)*currentmu)
    nextmu <- currentDelta*currentmu_filt + (1-currentDelta)*mu
    mu_mat[, j+1] <- nextmu
  }
  mu_mat
}





#evaluates negative log likelihood; the function to optimize
poisLogLikvec <- function(params, z, w, XX, stationary=TRUE){
  params_unpacked <- unpack1(params)
  n_time <- ncol(z)
  lin_comps <- XX %*% params_unpacked$beta
  probs <- matrix(rep(inv_logit(lin_comps), n_time), ncol = n_time, 
                  nrow = nrow(XX))
  mu_recurred <- mu_recursion(z=z, w=w, probs=probs, 
                              mu = params_unpacked$mu, 
                              Delta= params_unpacked$delta, 
                              stationary=stationary)[, 1:ncol(z)]
  myloglik <- ifelse(w, 
                     dpois(z, lambda = mu_recurred*probs, log=TRUE), 
                     0)
  -myloglik
}

tryloglik <- poisLogLikvec(rep(0, 10), z=Y_train, w=W_train, XX=XX_train)
mean(tryloglik)
head(tryloglik)

poisLogLik <- function(params){
  mynll <- poisLogLikvec(params=params, z=Y_train, w=W_train, XX=XX_train, 
                         stationary=TRUE)
  sum(mynll)/sum(W_train)
}



logsumexp <- function(x){
  m <- max(x)
  m + log(sum(exp(x-m)))
}

logsumexp_mat <- function(x){
  m <-do.call(pmax, x)
  x_mat <- do.call(cbind, x)
  m + log(rowSums(exp(x_mat-m)))
}

#for debugging
# params <-   c(rep(0, 9), log(2*mean(Y_train, na.rm=TRUE)))
# z <- Y_train
# w <- W_train
# XX <- XX_train

#evaluates negative log likelihood; the function to optimize
poisLogLikIntvec <- function(params, z, w, XX, stationary=TRUE){
  params_unpacked <- unpack1(params)
  mu <- params_unpacked$mu
  n_time <- ncol(z)
  n <- nrow(z)
  probs <- inv_logit(XX %*% params_unpacked$beta)
  max_z <- max(qpois(1e-12, lambda = probs*mu, lower.tail = FALSE))
  logliklist <- rep(list(NULL), max_z+1)
  new_w <- cbind(1, w)
  newprobs <- matrix(rep(probs, n_time+1), nrow = nrow(z), 
                     ncol = n_time+1)
  for(i in 1:(max_z+1)){
    currentz <- i-1
    new_z <- cbind(currentz, z)
    logliklist[[i]] <- -rowSums(poisLogLikvec(params = params, z=new_z, 
                                              w=new_w, XX = XX, 
                                              stationary=stationary))
    
  }
  #logliklist
  loglikvec_added <- logsumexp_mat(logliklist)
  -loglikvec_added
}

#tryloglikint <- poisLogLikIntvec(rep(0, 10), z=Y_train, w=W_train, XX=XX_train)

poisLogLikInt <- function(params){
  mynll <- poisLogLikIntvec(params=params, z=Y_train, w=W_train, XX=XX_train, 
                            stationary=TRUE)
  sum(mynll)/sum(W_train)
}


mypois_int <- optim(c(rep(0, 9), log(2*mean(Y_train, na.rm=TRUE))), poisLogLikInt, 
                    control = list(maxit = 5000), 
                    method = "L-BFGS-B")
unpacked1_int <- unpack1(mypois_int$par)
mypois_int$convergence
mypois_int$counts
unpacked1_int
beta1_int <- unpacked1_int$beta
delta1_int <- unpacked1_int$delta
mu1_int <- unpacked1_int$mu

probs1_int <- inv_logit(XX_train %*% beta1_int)
probs1_int_mat <-  matrix(rep(probs1_int, ncol(Y_train)+1), 
                          nrow = nrow(Y_train), ncol = ncol(Y_train)+1)
max_z <- max(qpois(1e-12, lambda = probs1_int*mu1_int, lower.tail = FALSE))
new_w <- cbind(1, W_train)
loglikvec_pois <- rep(list(NULL), max_z+1)
msevec_pois <- rep(NA, max_z+1)
n_time <- ncol(Y_train)
for(i in 1:(max_z+1)){
  currentz <- i-1
  dens_z <- dpois(currentz, lambda = probs1_int[test_ind]*mu1_int, log=TRUE)
  new_z <- cbind(currentz, Y_train)
  mu_recurred_future <- mu_recursion(z=new_z, w=new_w, probs=probs1_int_mat, 
                                     mu = mu1_int, 
                                     Delta= delta1_int)[test_ind, n_time+2]
  preds1_int <- mu_recurred_future*probs1_int[test_ind]
  loglikvec_pois[[i]] <- dens_z + dpois(Y_test, lambda = preds1_int, log = TRUE)
  msevec_pois[i] <- mean(exp(dens_z) * (Y_test - preds1_int)^2)
}

#utilities
inv_logit <- function(x) 1/(1+exp(-x))


unpack_nb <- function(params){
  beta <- params[1:p]
  delta <- inv_logit(params[p+1])
  a <- exp(params[p+2])
  list(beta=beta, delta=delta, a=a)
}

unpack_nb_rand <- function(params){
  beta <- params[1:p]
  a <- exp(params[p+1])
  list(beta=beta, a=a)
}

ldnbinom2 <- function(x, a, b, lam){
  dnbinom(x=x, size = a, mu = a*lam/b, log=TRUE)
}




#calculates ab's
#once ab's are calculated, dnbinom can be evaluated on 
#desired a's, b's and lambda's
ab_recursion <- function(z, w, lambda, Delta, a){
  n_pol <- nrow(z)
  n_col <- ncol(z)
  a_mat <- matrix(NA, nrow = n_pol, ncol = n_col+1)
  b_mat <- matrix(NA, nrow = n_pol, ncol = n_col+1)
  a_mat[, 1] <- a
  b_mat[, 1] <- a
  for(j in 1:(n_col)){
    currenta <- a_mat[, j]; currentb <- b_mat[, j]
    currentz <- z[, j];     currentw <- w[, j]
    currentlam <- lambda[, j]
    currenta_filt <- ifelse(currentw==0, 
                            currenta, 
                            currenta + currentz)
    currentb_filt <- ifelse(currentw==0, 
                            currentb, 
                            currentb + currentlam)
    q <- 1/(Delta^2 +(1-Delta^2)*(currentb_filt)/a)
    nextb <- q*currentb_filt
    nexta <- Delta*q*currenta_filt + (1-Delta)*nextb
    a_mat[, j+1] <- nexta
    b_mat[, j+1] <- nextb
  }
  list(a_mat = a_mat, b_mat = b_mat)
}


#evaluates negative log likelihood; the function to optimize
negbinLogLik <- function(params, z, w, XX){
  params_unpacked <- unpack_nb(params)
  n_time <- ncol(z)
  lin_comps <- XX %*% params_unpacked$beta
  lam <- matrix(rep(exp(lin_comps), n_time), ncol = n_time, 
                nrow = nrow(XX))
  ab_recurred <- ab_recursion(z=z, w=w, lambda = lam, 
                              a = params_unpacked$a, 
                              Delta= params_unpacked$delta)
  a_recurred <- ab_recurred$a_mat[, 1:ncol(z)]
  b_recurred <- ab_recurred$b_mat[, 1:ncol(z)]
  myloglik <- ifelse(w, 
                     ldnbinom2(z, a = a_recurred, b=b_recurred, lam=lam), 
                     0)
  sum(-myloglik)/sum(w)
}

negbinLogLik(rep(0, 10), z=Y_train, w=W_train, XX=XX_train)




# negbinLogLik_rand <- function(params, z, w, XX){
#   params_unpacked <- unpack_nb_rand(params)
#   newparams <- c(params_unpacked$beta, Inf, log(params_unpacked$a))
#   negbinLogLik(params = newparams, 
#                z=z, w=w, XX=XX)
# }


negbinLogLik_rand <- function(params, z, w, XX){
  params_unpacked <- unpack_nb_rand(params)
  n_time <- ncol(z)
  lin_comps <- XX %*% params_unpacked$beta
  lam <- matrix(rep(exp(lin_comps), n_time), ncol = n_time, 
                nrow = nrow(XX))
  ab_recurred <- ab_recursion(z=z, w=w, lambda = lam, 
                              a = params_unpacked$a, 
                              Delta= 1)
  a_recurred <- ab_recurred$a_mat[, 1:ncol(z)]
  b_recurred <- ab_recurred$b_mat[, 1:ncol(z)]
  myloglik <- ifelse(w, 
                     ldnbinom2(z, a = a_recurred, b=b_recurred, lam=lam), 
                     0)
  sum(-myloglik)/sum(w)
}


mynegbin <- optim(c(log(mean(Y_train, na.rm=TRUE)), rep(0, 9)), negbinLogLik, z=Y_train, 
                  w=W_train, XX = XX_train, 
                  control = list(maxit = 5000, factr = 1e7), 
                  method = "L-BFGS-B")
unpacked_nb <- unpack_nb(mynegbin$par)
beta_nb <- unpacked_nb$beta
delta_nb <- unpacked_nb$delta
a_nb <- unpacked_nb$a

lambda_nb <- exp(XX_train %*% beta_nb)
lambda_nb_mat <-  matrix(rep(lambda_nb, ncol(Y_train)), 
                         nrow = nrow(Y_train), ncol = ncol(Y_train))

ab_future <- ab_recursion(z=Y_train, w=W_train, 
                          lambda = lambda_nb_mat, Delta=delta_nb, a=a_nb)
a_future <- ab_future$a_mat[test_ind, 6]
b_future <- ab_future$b_mat[test_ind, 6]
preds_nb <- a_future/b_future * lambda_nb[test_ind]











mynegbin_rand <- optim(c(log(mean(Y_train, na.rm=TRUE)), rep(0, 8)), negbinLogLik_rand, z=Y_train, 
                       w=W_train, XX = XX_train, 
                       control = list(maxit=5000, factr = 1e7), 
                       method = "L-BFGS-B")
unpacked_nb_rand <- unpack_nb_rand(mynegbin_rand$par)

beta_nb_rand <- unpacked_nb_rand$beta
Delta_nb_rand <- 1
a_rand <- unpacked_nb_rand$a
lambda_rand <- exp(XX_train %*% beta_nb_rand)
lambda_rand_mat <- matrix(rep(lambda_rand, ncol(Y_train)), ncol = ncol(Y_train), 
                          nrow = nrow(Y_train))
ab_rand_future <- ab_recursion(z=Y_train, 
                               w = W_train, 
                               lambda = lambda_rand_mat, 
                               Delta = Delta_nb_rand, 
                               a = a_rand)
a_rand_future <- ab_rand_future$a_mat[test_ind, 6]
b_rand_future <- ab_rand_future$b_mat[test_ind, 6]
preds_nb_rand <- a_rand_future/b_rand_future * lambda_rand[test_ind]


mypois_int$convergence
mypois_int$counts
unpacked1_int


sum(logsumexp_mat(loglikvec_pois))
sum(msevec_pois)


mynegbin_rand$convergence
mynegbin_rand$counts

unpacked_nb_rand
mean((preds_nb_rand-Y_test)^2)
sum(ldnbinom2(Y_test, 
              a = a_rand_future, 
              b = b_rand_future, 
              lam = lambda_rand[test_ind]))



mynegbin$convergence
mynegbin$counts
unpacked_nb

mean((preds_nb-Y_test)^2)
sum(ldnbinom2(Y_test, a = a_future, b=b_future, 
              lam = lambda_nb[test_ind]))


