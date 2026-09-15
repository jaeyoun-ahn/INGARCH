load("data.RData")
head(data)

dim(data)
length( unique(data$PolicyNum) )
data$col.freq<-data$FreqIM      #Collision old and new
data$col.Cov<-data$CoverageIM #Coverage old and new

index<-(data$CoverageIM>2e-6)
data.train<-data[index,]
dim(data.train)
length(unique(data.train$PolicyNum))

mycov<-sort(data.train$col.Cov)
length(mycov)
#data.train$col.Cov.Idx<- as.factor((data.train$col.Cov>0.635)  +1)
4591/3
4591*2/3
data.train$col.Cov.Idx<- as.factor((data.train$col.Cov>quantile(mycov, 1/3))+(data.train$col.Cov>quantile(mycov, 2/3))  )
unique(data.train$col.Cov.Idx)

#Define frequency and severity of collision of old and new
data.train$n <-data.train$col.freq
data.train$s <- data.train$yAvgIM * data.train$FreqIM


#Define frequency and severity of collision of old and new
index2<-(data.train$col.freq>0)
data.train.sev<-data.train[index2,]
data.train.sev$m<- data.train.sev$yAvgIM * data.train.sev$FreqIM


load("dataout.RData")
head(dataout)

# Check whether the data is is the set
x<-c(1,2,3,4,5)
y<-c(2,4,6, 2)
# match(y, x)
y %in% x

# Out of test data, use index which is from the training data & coverage>0
out.idx1 <- dataout$PolicyNum %in% data.train$PolicyNum
out.idx2 <- dataout$CoverageIM>2e-6
out.idx <- out.idx1 & out.idx2
length( unique(data.train$PolicyNum) )
sum(out.idx)
data.valid<-dataout[out.idx,] #test data


# Now, define explanatory variable in test data. 
# (I am not sure the explanatory variable will be ever used or not.)
out.idx.sev<- data.valid$PolicyNum %in% data.train$PolicyNum
data.valid$col.freq<-data.valid$FreqIM
data.valid$col.Cov<-data.valid$CoverageIM
data.valid$col.Cov.Idx<-as.factor((data.valid$col.Cov>mycov[1540])+(data.valid$col.Cov>mycov[3081])  )
#data.valid$col.Cov.Idx<- as.factor((data.valid$col.Cov>0.635)  +1)
data.valid$s <- data.valid$yAvgIM * data.valid$FreqIM 
data.valid$n <- data.valid$FreqIM

length(data.valid$n )

