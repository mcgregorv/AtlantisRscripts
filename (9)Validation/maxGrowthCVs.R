# read in the cvs estimated in compareGrowthByAgeClass.R and check which groups have what maxs

thisData<-read.csv(paste(DIR$'Tables', "sizeAtAgeCVs.csv", sep=""))

test<-apply(thisData[,-1], 1, sum, na.rm=TRUE)
index<-test>0
thisData<-thisData[index,]

test<-apply(thisData[,-1],1, max, na.rm=TRUE)

# test[test>0.2]
# 6        27        33        40        41        42 
# 0.2971678 0.2069821 0.4050623 0.2010993 0.2522480 0.3802589 
# 
# thisData[test>0.24,1]
# [1] BEE LIN PFM PFS

