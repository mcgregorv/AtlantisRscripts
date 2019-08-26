# analyse stability using tracers rel to base runs - tracers were read in and stored in readInChaosRunsAndSaveForAnalyses.R
# use allRankings as possible explanatory variables - this table was created in network analyses
# have 2 main sets of runs so far - SampleA and SampleKEY - both have fished and non-fished version
# Use DChaosBASE and DChaosFISH to compare
basePath <-paste(DIR$'Base',"AtlantisModels\\chaos\\",sep="")
thisChaosVersion <- "SampleA";
thisChaosVersion <- "FISHSampleA";
thisChaosVersion <- "SampleKEY";
thisChaosVersion <- "FISHSampleKEYPlusBP";
thisChaosVersion <- "SampleKEYPlusBP"
nChaosRuns <- 35

#read in all the rankings - TL, lifespan, keystoneness,...
allTheRankings <- read.csv(paste(DIR$'Tables', "allTheRankings.csv", sep="")); nrankings <- dim(allTheRankings)[2]
groupsDF<- read.csv(paste(basePath,"..\\CRAM_groups.csv", sep="")); ng<-dim(groupsDF)[1]

## add in species growth rates (for age-structured)

## add informance - how well defined the group is combined with how well they performed in the historical model
rankDF<-read.csv(paste(DIR$'Tables',"interactionEffectsRANKINGs.csv",sep=""))

thisCex=1.65

plotPath <- paste(DIR$'Reports',"Chaos\\Figures\\", sep="") ## paper version
plotPath <- paste(DIR$'Figures',"Chaos\\", sep="")

chaosVersions <- c("SampleKEYPlusBP"); nchaosVersions <- length(chaosVersions)
chaosRunDescr <- c("A: Weighted keystone - with fishing")

nruns <- 35

storeRelResponse <- array(NA, dim=c(nchaosVersions, nruns, ng, nts))
for(c in 1:nchaosVersions){
  thisChaosVersion <- chaosVersions[c]
  
  # nChaosRuns <- 35
  groupsDF <- read.csv(paste(basePath, "..\\CRAM_groups.csv",sep="")); ng <- dim(groupsDF)[1]
  nts <- 151; nlayers <-6
  
  # load baseArray - has tracers from Chaos BASE and FISH runs (c("outputDChaosFISH", "outputDChaosBASE"))
  load(paste(basePath,"ChaosNtracersBASEandBaseFISH_baseAarray",sep=""))
  
  # bring in storeNarray
  load(paste(basePath,"ChaosNtracers",thisChaosVersion,sep=""))
  
  #summarise the largest proportional differences - perhaps from 50 years, 100 years, and at the end (151 years)
  allRelTracers <- 0*storeNarray
  for(g in 1:ng){
    thisCode<-as.character(groupsDF$Code[g]); 
    thisBaseTracer <- baseArray[1,g,]
    theseChaosTracers <- storeNarray[,g,]
    
    relTracers <- 0*theseChaosTracers; nruns <- dim(theseChaosTracers)[1]
    for(r in 1:nruns){
      allRelTracers[r,g,]<-(theseChaosTracers[r,] - thisBaseTracer)/thisBaseTracer
    }
  }
  storeRelResponse[c,,,]<-allRelTracers
  
}

# turn into df with response and explanatory vars
this_timeStep <- 150
# test for actual nts
test <- apply(storeRelResponse,4, sum, na.rm=TRUE)
this_nts <- max(seq(1,nts)[test>0])
this_timeStep <- this_nts
allData2fit <- data.frame(matrix(NA, ncol=(nrankings+6), nrow=(nruns * nchaosVersions * ng)))
colnames(allData2fit)<- c("ChaosRun", "relDiff", "relDiffLast20years", "maxDiff", "medDiff", "Informance", colnames(allTheRankings))

for(c in 1:nchaosVersions){
  for(r in 1:nruns){
    for(g in 1:ng){
      thisRow <- ng*nruns*(c-1)  + ng * (r-1) + g
      allData2fit$ChaosRun[thisRow] <- c; 
      allData2fit$relDiff[thisRow] <- storeRelResponse[c,r,g,this_timeStep]
      allData2fit$relDiffLast20years[thisRow] <- mean(abs(storeRelResponse[c,r,g,(this_nts-20):this_nts]), na.rm=TRUE)
      allData2fit$maxDiff[thisRow] <- max(abs(storeRelResponse[c,r,g,35:this_nts]), na.rm=TRUE)
      allData2fit$medDiff[thisRow] <- median(abs(storeRelResponse[c,r,g,35:this_nts]), na.rm=TRUE)
      thisCode <- as.character(groupsDF$Code[g])
      allData2fit[thisRow,c(colnames(allTheRankings))] <- allTheRankings[allTheRankings$Code == thisCode,]
      allData2fit$Code[thisRow] <- thisCode
      thisNumCohorts <- groupsDF$NumCohorts[g]
      if(thisNumCohorts==1){
        allData2fit$Keystone[thisRow] <- 0
      }
      thisInf <- rankDF$informPerformRank[rankDF$Code==thisCode]
      if(length(thisInf)==0){thisInf <- 16}
      allData2fit$Informance[thisRow] <- thisInf
    }
  }
}

test<- is.na(allData2fit$relDiff)
table(test) # some runs did not complete - take them out for this analyses as it  won't be clear which species group it was because of

data2fit <- allData2fit[!test,]

# some of the response vars need to be factors - the categorical or ordered ones. NO, they need their order
# data2fit$Keystone <- factor(data2fit$Keystone, levels=0:max(data2fit$Keystone))
# data2fit$Response <- factor(data2fit$Response, levels=unique(sort(data2fit$Response)))
head(data2fit)
# ChaosRun      relDiff    maxDiff Code       TL Keystone Response NumL1cons Lifespan   propAdM propJuvM       B0 PropByTopPrey
# 56        1 -0.297068344 0.29808869  ASQ 4.751959       19       17        14        2        NA       NA 16255.94     0.7653883
# 57        1  0.008208659 0.05560994  BAL 4.405262       25       27         1       80        NA       NA  6000.22     0.9610216
# 58        1  0.145957021 0.17490656   BB 0.100000        0       NA         1        1        NA       NA 78915.15            NA
# 59        1 -0.067463625 0.10121014   BC 1.100000        0       NA        20        1        NA       NA 25075.40     0.9573303
# 60        1  0.092585810 0.10416888   BD 1.104599        0       NA        17        1        NA       NA 41715.70     0.9949211
# 61        1 -0.178250278 0.17825028  BEE 3.597647       18       22        15       30 0.7961681        0 20907.44     0.6185041


getBestLModel <- function(thisVar, thisFactors, cutoff, data2fit){
  this_df <- data.frame(matrix(NA, ncol=2, nrow=length(thisFactors))); colnames(this_df)<-c("Factor", "Rsq")
  existingRsq <- 0
  modelFactors <- c()
  leftFactors<-thisFactors
  keepGoing <- TRUE
  i <- 1
  while(keepGoing==TRUE){
    testRsq <- c(); keep4nowFactors<-c()
    for(f in 1:length(leftFactors)){
      testFactors <- c(modelFactors, leftFactors[f])
      testFormula <- as.formula(paste(thisVar,"~", paste(testFactors, collapse="+")))  
      testModel <- glm(testFormula, data=data2fit)
      thisRsq <- 1-(summary(testModel)$"deviance"/summary(testModel)$"null.deviance")
      if((thisRsq-existingRsq) >cutoff){
        keep4nowFactors <- c(keep4nowFactors, leftFactors[f])
        testRsq <- c(testRsq, thisRsq)
      }
    }
    if(length(keep4nowFactors)>0){
      if((max(testRsq)-existingRsq)>cutoff){
        # then keep the top one
        index <- testRsq == max(testRsq)
        thisFactor <- keep4nowFactors[index]; thisRsq <- testRsq[index]
        thisRsq_improvement <- thisRsq - existingRsq
        this_df[i,] <- c(thisFactor, thisRsq)
        leftFactors <- leftFactors[leftFactors != thisFactor]
        existingRsq <- thisRsq
        modelFactors <- c(modelFactors, thisFactor)
        i <- i +1
      }else{
        keepGoing <- FALSE
      }
    } else{
      keepGoing=FALSE
    }
  }
  index<-!is.na(this_df$Factor)
  return(this_df[index,])
}

thisFactors <- c("TL", "Keystone", "Response", "NumL1cons", "Lifespan",   "propAdM", "propJuvM",       "B0", "PropByTopPrey")
nfactors <- length(thisFactors)

getBestLModel(thisVar="maxDiff", thisFactors = thisFactors, cutoff=0.01, data2fit)

getBestLModel(thisVar="medDiff", thisFactors = thisFactors, cutoff=0.01, data2fit)

getBestLModel(thisVar="relDiff", thisFactors = thisFactors, cutoff=0.01, data2fit)

getBestLModel(thisVar="relDiffLast20years", thisFactors = thisFactors, cutoff=0.01, data2fit)

# take out any zeros
this_data2fit <- data2fit[abs(data2fit$relDiffLast20years) >0,]

this_data2fit$logRelDiff <- log(abs(this_data2fit$relDiffLast20years))
this_data2fit$cubeRoot <- (abs(this_data2fit$relDiffLast20years))^{1/3}

thisVar <- "cubeRoot"
# thisVar <- "logRelDiff"
# thisVar <- "relDiffLast20years"

this_df <- getBestLModel(thisVar=thisVar, thisFactors = thisFactors, cutoff=0.01, data2fit=this_data2fit)
thisFormula <- as.formula(paste(thisVar, "~", paste(this_df$Factor, collapse=" + "), sep=""))
this_df 

test<-apply(this_data2fit[,c("relDiffLast20years", this_df$Factor)], 1, sum)
index <- is.na(test) # need all vars to be not na
this_data2fit <- this_data2fit[!index,]

this_model <- glm(thisFormula, data=this_data2fit)

testFit <- predict.glm(this_model)

thisMin <- min(c(testFit, this_data2fit[,c(thisVar)]), na.rm=TRUE)
thisMax <- max(c(testFit, this_data2fit[,c(thisVar)]), na.rm=TRUE)
thisLine <- seq(thisMin, thisMax, length.out = 100)
plot(x=this_data2fit[,c(thisVar)], y=testFit, pch=20, ylim=c(thisMin, thisMax), xlim=c(thisMin, thisMax))
points(x=thisLine, y=thisLine, type="l", col="red", lwd=2, lty=2)

########################################
## keep only factors that work for biomass pool groups as well
bp_factors <- c("TL",     "Keystone",   "NumL1cons",     "Lifespan" ,  "B0", "Informance")

bp_data2fit <- this_data2fit
test<-apply(bp_data2fit[,c("relDiffLast20years", bp_factors)], 1, sum)
index <- is.na(test) # need all vars to be not na
bp_data2fit <- bp_data2fit[!index,]
# check for any zeros in resp var
test <- bp_data2fit$relDiffLast20years==0
bp_data2fit<-bp_data2fit[!test,]

bp_data2fit$logRelDiff <- log(bp_data2fit$relDiffLast20years)
bp_data2fit$cubeRoot <- unlist(lapply(bp_data2fit$relDiffLast20years, FUN=function(x){x^(1/3)}))


thisVar <- "cubeRoot"
thisVar <- "logRelDiff"
this_df <- getBestLModel(thisVar=thisVar, thisFactors = bp_factors, cutoff=0.01, data2fit=bp_data2fit)
thisFormula <- as.formula(paste(thisVar, "~", paste(this_df$Factor, collapse=" + "), sep=""))
this_df

this_model <- glm(thisFormula, data=bp_data2fit)

testFit <- predict.glm(this_model)

thisMin <- min(c(testFit, bp_data2fit[,c(thisVar)]), na.rm=TRUE)
thisMax <- max(c(testFit, bp_data2fit[,c(thisVar)]), na.rm=TRUE)
thisLine <- seq(thisMin, thisMax, length.out = 100)

plot(x=bp_data2fit[,c(thisVar)], y=testFit, pch=20, ylim=c(thisMin, thisMax), xlim=c(thisMin, thisMax))
points(x=seq(-1,1), y=seq(-1,1), type="l", col="red", lwd=2, lty=2)
points(x=thisLine, y=thisLine, type="l", col="red", lwd=2, lty=2)
# 
# plot(this_model)
# 
# timeIndex <- 1:151
# thisYmax <- 100*max(allRelTracers[,,timeIndex], na.rm=TRUE); thisYmin <- 100*min(allRelTracers[,,timeIndex], na.rm=TRUE)
# plot(timeIndex, type="n", ylim=c(thisYmin, thisYmax), ylab="Relative difference (%)", xlab="")
# for(t in timeIndex){
#   thisData<-100*allRelTracers[,,t]
#   boxplot(as.double(thisData), at=t, outline=TRUE, col=myGrey_trans, add=TRUE, pch=20, cex=0.8, border=myGrey, yaxt="n")  
# }
# # par(new=TRUE)
# plot(timeIndex, type="n", ylim=c(-10, 10), ylab="", xlab="")
# for(t in timeIndex){
#   thisData<-100*allRelTracers[chaosDescIndex,,t]
#   boxplot(as.double(thisData), at=t, outline=FALSE, col=myBlue_trans, add=TRUE, pch=20, cex=0.8, border=myBlue, yaxt="n")  
# }


## test for correlation between explan, vars
nfactors <- length(thisFactors)
allCors <- array(NA, dim=c(nfactors, nfactors))
for(i in 1:nfactors){
  for(j in 1:nfactors){
    # if (i != j){
    y1 <- allData2fit[,c(thisFactors[i])]; y2 <- allData2fit[,c(thisFactors[j])]
    thisIndex <- !is.na(y1) & !is.na(y2)
    thisCor <- cor(x=y1[thisIndex], y=y2[thisIndex], method="spearman")
    allCors[i,j]<-thisCor
  }
  # }
}



# plot(this_model)


##############################################

# if don't use the fishing ones:
getIndexfromC <- function(c){
  index <- c()
  for(r in 1:nruns){
    for(g in 1:ng){
      this_i <- ng * nruns * (c-1) + ng * (r-1) + g
      index <- c(index, this_i)
    }
  }
  return(index)
}
# c1index <- getIndexfromC(c=1); c3index <- getIndexfromC(c=3)
chaosDescIndex <- c(1,3); chaosRunIndex <- unlist(lapply(chaosDescIndex, getIndexfromC))

this_data2fit <- allData2fit[chaosRunIndex,]
index <- is.na(this_data2fit$relDiff)

this_data2fit<-this_data2fit[!index,]

getBestLModel(thisVar="relDiff", thisFactors = thisFactors, cutoff=0.01, data2fit=this_data2fit)

getBestLModel(thisVar="maxDiff", thisFactors = thisFactors, cutoff=0.01, data2fit=this_data2fit)
getBestLModel(thisVar="medDiff", thisFactors = thisFactors, cutoff=0.01, data2fit=this_data2fit)
getBestLModel(thisVar="relDiffLast20years", thisFactors = thisFactors, cutoff=0.01, data2fit=this_data2fit)

# set the best one up as a model
this_data2fit$logRelDiff <- log(this_data2fit$relDiffLast20years)
this_data2fit$cubeRoot <- this_data2fit$relDiffLast20years^(1/3)
thisVar <- "cubeRoot"
thisVar <- "logRelDiff"
this_df <- getBestLModel(thisVar=thisVar, thisFactors = thisFactors, cutoff=0.01, data2fit=this_data2fit)
thisFormula <- as.formula(paste(thisVar, "~", paste(this_df$Factor, collapse=" + "), sep=""))
this_df 

test<-apply(this_data2fit[,c("relDiffLast20years", this_df$Factor)], 1, sum)
index <- is.na(test) # need all vars to be not na
this_data2fit <- this_data2fit[!index,]

this_model <- glm(thisFormula, data=this_data2fit)

testFit <- predict.glm(this_model)

plot(x=this_data2fit$relDiffLast20years, y=testFit, pch=20)
points(x=seq(-1,1), y=seq(-1,1), type="l", col="red", lwd=2, lty=2)

plot(this_model)

