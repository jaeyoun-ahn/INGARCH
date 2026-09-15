source("aux_data_manage4.R")

.libPaths(c("Rlib", .libPaths()))
library(RTMB)

W_train <- !is.na(Y_train)
p <- dim(XX_train)[3]
mean_y <- mean(Y_train, na.rm=TRUE)
source("funs_INGARCH.R")

obj_state_nb <- make_state_nb_objective()
fit_state_nb <- optim(
  obj_state_nb$par, obj_state_nb$fn, obj_state_nb$gr,
  method="L-BFGS-B", lower=c(rep(-Inf, p), 0, -Inf),
  upper=c(rep(Inf, p), .999, Inf)
)
u_state <- unpack_state_nb(fit_state_nb$par)

obj_random <- make_state_nb_objective(FALSE)
fit_random <- optim(obj_random$par, obj_random$fn, obj_random$gr,
                    method="L-BFGS-B")
u_random <- unpack_state_nb(fit_random$par, FALSE)

B <- as.integer(Sys.getenv("NB_BOOTSTRAP_REPLICATES", "1000"))
bootstrap_time <- system.time(
  NB_bootstrap <- bootstrap_state_nb(B=B, workers=11, seed=20260909)
)

output <- paste0("NB_INGARCH_bootstrap_City_reference_", B)
saveRDS(NB_bootstrap, paste0(output, ".rds"))
write.csv(NB_bootstrap$summary, paste0(output, ".csv"),
          row.names=FALSE)

print(NB_bootstrap$summary)
print(table(NB_bootstrap$convergence[, "state"]))
print(table(NB_bootstrap$convergence[, "random"]))
print(bootstrap_time)
