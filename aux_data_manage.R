source("aux_data_manage2.R")

#Inspect the data
head(data.train, 1)           #training data 
head(data.valid,1)        #test data
table(data.train$col.Cov.Idx) 
# n:frequency, s:aggregate loss, Type: one-hot encoding, col.Cov.Idx: 0,1,2
# id: PolicyNum
mean(data.train$PolicyNum %in%  data.valid$PolicyNum) # not all training data is in test data
mean(data.valid$PolicyNum %in% data.train$PolicyNum) # all test data is in training data


unique_policy = unique(data.train$PolicyNum)
n_pol = length(unique_policy)
Y_train = matrix(NA, nrow=n_pol, ncol=5)
ID_train = rep(NA, n_pol)
X_train = matrix(NA, nrow=n_pol, ncol=7)
Length_observed_train = rep(NA, n_pol)

for(i in 1:n_pol){
  temp_idx = which(unique_policy[i]== data.train$PolicyNum  )
  year_idx = data.train$Year[temp_idx]-2005
  len_temp = length(temp_idx)
  for(t in 1:5){
    Y_train[i, year_idx] = c(data.train$n[temp_idx])
    #Y_train[i, 1:5] = c(data.train$n[temp_idx], rep(NA, 5-len_temp))
  }
  X_train[i, 1:6] = with(c(TypeCity[temp_idx[1]], TypeCounty[temp_idx[1]], TypeSchool[temp_idx[1]], TypeTown[temp_idx[1]], 
                           TypeVillage[temp_idx[1]], col.Cov.Idx[temp_idx[1]]), data=data.train)
  ID_train[i] = unique_policy[i]
  Length_observed_train[i] = len_temp
}

X_train <- as.data.frame(X_train, stringsAsFactors = FALSE)
colnames(X_train) <- c(
  "CityType",
  "CountyType",
  "SchoolType",
  "TownType",
  "VillageType",
  "OtherCov"
)
XX_train = model.matrix(~ CityType+CountyType+SchoolType+TownType+VillageType+factor(OtherCov), data=X_train)
head(XX_train)





n_test = length(data.valid$PolicyNum)
X_test = matrix(NA, nrow=n_pol, ncol=6)
Y_test = rep(NA, nrow=n_test)
ID_test = rep(NA, nrow=n_test)
i=100
which(data.valid$PolicyNum[i]== unique_policy)

for(i in 1:n_test){
  temp_idx = which(data.valid$PolicyNum[i]== unique_policy  )
  Y_test[i] = data.valid$n[i]
  ID_test[i] = unique_policy[temp_idx]
  X_test[i, 1:6] = with(c(TypeCity[i], TypeCounty[i], TypeSchool[i], TypeTown[i], 
                          TypeVillage[i], col.Cov.Idx[i]), data=data.valid)
}

X_test <- as.data.frame(X_test, stringsAsFactors = FALSE)
colnames(X_test) <- c(
  "CityType",
  "CountyType",
  "SchoolType",
  "TownType",
  "VillageType",
  "OtherCov"
)
XX_test = model.matrix(~ CityType+CountyType+SchoolType+TownType+VillageType+factor(OtherCov), data=X_test)
head(XX_test, 1)

# Input DATA


head(Y_train)
head(XX_train)
IDX = !is.na(Y_train)
num_obs= rowSums(!is.na(Y_train))
n = dim(Y_train)[1]
p = dim(XX_train)[2]
k = dim(XX_train)[2]

betaMean  <- rep(0, p)        
betaCov   <- diag(1, p)
dim(XX_train)
IDXint = matrix(NA, nrow=n, ncol=5)
for(i in 1:n){
  obsCols    <- which(IDX[i,])
  IDXint[i, 1:num_obs[i]] = obsCols
}
length(betaMean)
dim(betaCov)
















num_obs1 = num_obs
IDXint1 = IDXint

I_train = (Y_train>0)
N_train = Y_train-1
N_train[N_train==-1] <- NA

IDX2 = !is.na(N_train) 
IDXint2 = matrix(NA, nrow=n, ncol=5)

IDX3 = which(rowSums(!is.na(N_train))>0)
length(IDX3)

num_obs2= rowSums(!is.na(N_train))
for(i in IDX3){
  # print(i)
  obsCols2    <- which(IDX2[i,])
  IDXint2[i, 1:num_obs2[i]] = obsCols2
}



###############################
#### Final data to be used ####
###############################

## Data for classical zero-inflated model
I_mat  <- I_train
N_mat  <- N_train
wI_mat <- !is.na(I_mat) * 1L     # 1 = 관측, 0 = NA
wN_mat <- !is.na(N_mat) * 1L
I_mat[is.na(I_mat)] <- 0L        # NA → 0 (더미)
N_mat[is.na(N_mat)] <- 0L
T_I <- ncol(I_mat);  T_N <- ncol(N_mat)



## Data for state-space copula model
n   <- nrow(Y_train)
k   <- ncol(XX_train)

Y_mat <- Y_train        
wY_mat <- !is.na(Y_mat) * 1L
Y_mat[!wY_mat] <- 0L
ID_train
#############
# Train data: 
#############
# Y_mat
# n \times 5 matrix of counting observations

# XX_train
# n \times k matrix of explanatory variable including the intercept.
# Assume there are k-1 explanatory variables and one intercept so that it has k columns
# Assume that explanatory variable does not change with time so that it is not a 3 dim array.

# wY_mat
# n \times 5 matrix of 0,1 exposures

# ID_train
# length n vector of ID to be linked with test data


#############
# Test data: 
#############

Y_test
XX_test
ID_test

# Y_test
# n_test vector of counting observations

# XX_test
# n_test \times k matrix of explanatory variable including the intercept.
# Assume there are k-1 explanatory variables and one intercept so that it has k columns
# Assume that explanatory variable does not change with time so that it is not a 3 dim array.

# ID_test
# length n_test vector of ID to be linked with train data

















