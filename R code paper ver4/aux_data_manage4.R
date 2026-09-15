source("aux_data_manage3.R")

data.train$lnCoverageIM <- log(data.train$CoverageIM)
data.valid$lnCoverageIM <- log(data.valid$CoverageIM)

log_coverage <- matrix(NA_real_, n_pol, n_time)
for(r in seq_len(nrow(data.train)))
  log_coverage[policy_row[r], data.train$Year[r]-2005] <-
    data.train$lnCoverageIM[r]

# Keep entity type fixed and insert the coverage observed in each period.
for(t in seq_len(n_time))
  XX_train[, t, ] <- cbind(entity_train, log_coverage[, t])
dimnames(XX_train)[[3]] <- c(colnames(entity_train), "lnCoverageIM")
XX_test <- cbind(entity_test, lnCoverageIM=data.valid$lnCoverageIM)

X_train <- XX_train
X_test <- XX_test
