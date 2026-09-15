source("aux_data_manage2.R")

policies <- unique(data.train$PolicyNum)
first_row <- match(policies, data.train$PolicyNum)
policy_row <- match(data.train$PolicyNum, policies)
test_policy_row <- match(data.valid$PolicyNum, policies)

coverage_first <- data.train$CoverageIM[first_row]
data.train$CoverageIM1 <- coverage_first[policy_row]
data.valid$CoverageIM1 <- coverage_first[test_policy_row]
data.train$lnCoverageIM1 <- log(data.train$CoverageIM1)
data.valid$lnCoverageIM1 <- log(data.valid$CoverageIM1)

n_pol <- length(policies)
n_time <- 5
Y_train <- matrix(NA_real_, n_pol, n_time)
for(r in seq_len(nrow(data.train)))
  Y_train[policy_row[r], data.train$Year[r]-2005] <- data.train$n[r]

ID_train <- policies
ID_test <- data.valid$PolicyNum
Y_test <- data.valid$n
Length_observed_train <- rowSums(!is.na(Y_train))

entity_train <- model.matrix(
  ~ TypeCounty+TypeMisc+TypeSchool+TypeTown+TypeVillage,
  data=data.train[first_row, ]
)
entity_test <- model.matrix(
  ~ TypeCounty+TypeMisc+TypeSchool+TypeTown+TypeVillage,
  data=data.valid
)

# Repeat the first-period log coverage over all five training periods.
XX_train <- array(NA_real_, c(n_pol, n_time, ncol(entity_train)+1))
for(t in seq_len(n_time))
  XX_train[, t, ] <- cbind(entity_train, log(coverage_first))
dimnames(XX_train)[[3]] <- c(colnames(entity_train), "lnCoverageIM1")
XX_test <- cbind(entity_test, lnCoverageIM1=data.valid$lnCoverageIM1)

X_train <- XX_train
X_test <- XX_test
IDX <- !is.na(Y_train)
num_obs <- rowSums(IDX)
n <- nrow(Y_train)
p <- dim(XX_train)[3]
k <- p
betaMean <- rep(0, p)
betaCov <- diag(p)

IDXint <- matrix(NA_integer_, n, n_time)
for(i in seq_len(n)) IDXint[i, seq_len(num_obs[i])] <- which(IDX[i, ])
num_obs1 <- num_obs
IDXint1 <- IDXint

I_train <- Y_train>0
N_train <- Y_train-1
N_train[N_train==-1] <- NA
IDX2 <- !is.na(N_train)
IDX3 <- which(rowSums(IDX2)>0)
num_obs2 <- rowSums(IDX2)
IDXint2 <- matrix(NA_integer_, n, n_time)
for(i in IDX3) IDXint2[i, seq_len(num_obs2[i])] <- which(IDX2[i, ])

I_mat <- I_train
N_mat <- N_train
wI_mat <- 1L*!is.na(I_mat)
wN_mat <- 1L*!is.na(N_mat)
I_mat[is.na(I_mat)] <- 0L
N_mat[is.na(N_mat)] <- 0L
T_I <- ncol(I_mat)
T_N <- ncol(N_mat)

Y_mat <- Y_train
wY_mat <- 1L*!is.na(Y_mat)
Y_mat[!wY_mat] <- 0L
