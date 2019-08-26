# analyse stability using tracers rel to base runs - tracers were read in and stored in readInChaosRunsAndSaveForAnalyses.R
# use allRankings as possible explanatory variables - this table was created in network analyses
# have 2 main sets of runs so far - SampleA and SampleKEY - both have fished and non-fished version
# Use DChaosBASE and DChaosFISH to compare
basePath <-paste(DIR$'Base',"AtlantisModels\\chaos\\",sep="")
thisChaosVersion <- "SampleA";
thisChaosVersion <- "FISHSampleA";
# thisChaosVersion <- "SampleKEY";
thisChaosVersion <- "FISHSampleKEY";
thisChaosVersion <- "SampleKEYPlusBP"
# nChaosRuns <- 35

getPredVars <- function(this_df){
  temp<- this_df$Factor
  temp2 <- sort(unique(unlist(lapply(temp, FUN=function(x){unlist(str_split(x,":"))}))))
  return(temp2)
}
getVarValues <- function(thisVar){
  thisData <- bp_data2fit[,c(thisVar)]
  thisMin <- min(thisData, na.rm=TRUE); thisMax <- max(thisData, na.rm=TRUE)
  thisValues<- pretty (seq(thisMin  , thisMax))
  return(thisValues)
}
getVarValues <- function(thisVar){
  thisTab <- table(as_data2fit[,c(thisVar)])
  thisValues<- names(thisTab)
  if(length(thisValues)>10){
    thisValues <- pretty(as.double(thisValues), n=10)
  }
  return(thisValues)
}

getPercentile<-function(x, k){
  y <- sort(x[!is.na(x)]); ny <- length(y)
  i1<-round(ny*(k/2)); i2 <- round(ny*(1-(k/2)))
  return(y[c(i1, i2)])
}
getColor <- function(x, log=FALSE, base=10){
  thisCol<-"white"
  if(!is.na(x)){
    thisIndex <- round((x-thisMin)/(thisMax-thisMin),2)*100 +1
    # if(x<1){
    #   thisIndex <- round((x-thisMin)/(1-thisMin),2)*100 +1
    #   thisCol<-colByLess[thisIndex]
    # } else if(x>=1){
    #   thisIndex <- round((x-1)/(thisMax-1),2)*100 +1
    #   thisCol<-colByMore[thisIndex]
    # }
    if(log==TRUE){
      if(x<1){x<-1}
      thisIndex <- round((log(x, base=base)/log(thisMax, base=base)),2)*100 +1
    }
    if(length(thisIndex)>0){
      if(thisIndex > 101){
        thisIndex <- 101
      }
      if(sum(x, na.rm=TRUE)>0 & thisIndex>0){
        thisCol<-absColorRamp[thisIndex]
      }
    }
  }
  return(thisCol)
}
plotGrid<-function(x,y){
  thisX<-c(x-0.5,x-0.5,x+0.5,x+0.5); thisY<-c(y-0.5,y+0.5,y+0.5,y-0.5)
  index <- tempDF$x==x & tempDF$y==y
  thisCol<-plotColour[index]
  if(length(thisCol)>0){
    polygon(x=thisX,y=thisY,col=thisCol,border=NA)
  }
  return(NULL)
}
#read in all the rankings - TL, lifespan, keystoneness,...
allTheRankings <- read.csv(paste(DIR$'Tables', "allTheRankings.csv", sep="")); nrankings <- dim(allTheRankings)[2]
groupsDF<- read.csv(paste(basePath,"..\\CRAM_groups.csv", sep="")); ng<-dim(groupsDF)[1]

# load baseArray - has tracers from Chaos BASE and FISH runs (c("outputDChaosFISH", "outputDChaosBASE"))
load(paste(basePath,"ChaosNtracersBASEandBaseFISH_baseAarray",sep=""))

## add informance - how well defined the group is combined with how well they performed in the historical model
rankDF<-read.csv(paste(DIR$'Tables',"interactionEffectsRANKINGs.csv",sep=""))
thisGroupsDF <- groupsDF[groupsDF$Code != "DC",]; ng <- dim(thisGroupsDF)[1]
keepGindex <- grep("DC", groupsDF$Code, invert = TRUE)
thisCex=1.65

plotPath <- paste(DIR$'Reports',"Chaos\\Figures\\", sep="") ## paper version
plotPath <- paste(DIR$'Figures',"Chaos\\", sep="")

chaosVersions <- c("DChaosNtracers", "DChaosNtracersFISH", paste("ChaosNtracers",c("SampleInformanceAll", "FISHSampleInformanceAll", "SampleKEYPlusBP", "FISHSampleKEYPlusBP"), sep="")); nchaosVersions <- length(chaosVersions)
chaosRunDescr <- c("A: All up or down - no fishing", "B: All up or down - with fishing", "C: Weighted randon - no fishing", "D: Weighted random - with fishing", 
                   "E: Weighted Keystone - no fishing", "F: Weighted keystone - with fishing")

nruns <- 35
nts <- 151; nlayers <-6
storeRelResponse <- array(NA, dim=c(nchaosVersions, nruns, ng, nts))
storeRatioResponse <- storeRelResponse
for(c in 1:nchaosVersions){
  rm(storeNarray, tracersArray, thisBaseArray)
  thisChaosVersion <- chaosVersions[c]

  if(length(grep("with fishing",chaosRunDescr[c]))>0){
    thisBaseArray <- baseArray[1,,]
  } else{
    thisBaseArray <- baseArray[2,,]
  }
  # bring in storeNarray
  load(paste(basePath,thisChaosVersion,sep=""))
  # if(!exists("storeNarray")){
  #   storeNarray <- apply(tracersArray, c(1,2,3), FUN=function(x){x*mg_2_tonne*X_CN})
  #   thisBaseArray <- tracersArray[5,,]
  #   storeNarray <- tracersArray[-5,,]
  # } 
  if(dim(storeNarray)[2]==(ng+1)){storeNarray<-storeNarray[,keepGindex,]}
  if(dim(thisBaseArray)[1]==(ng+1)){thisBaseArray<-thisBaseArray[keepGindex,]}

  # # 
  #summarise the largest proportional differences - perhaps from 50 years, 100 years, and at the end (151 years)
  allRelTracers <- 0*storeNarray; allRatioTracers <- allRelTracers
  for(g in 1:ng){
    thisCode<-as.character(thisGroupsDF$Code[g]); 
    thisBaseTracer <- thisBaseArray[g,]
    theseChaosTracers <- storeNarray[,g,]
    # plot(thisBaseTracer, type="l", col="red")
    relTracers <- 0*theseChaosTracers; this_nruns <- dim(theseChaosTracers)[1]
    for(r in 1:this_nruns){
      allRelTracers[r,g,]<-(theseChaosTracers[r,] - thisBaseTracer)/thisBaseTracer
      allRatioTracers[r,g,] <- (theseChaosTracers[r,])/thisBaseTracer
      # points(theseChaosTracers[r,], type="l", col=myGreen)
    }
  }
  storeRelResponse[c,1:this_nruns,,]<-allRelTracers
  storeRatioResponse[c,1:this_nruns,,]<-allRatioTracers

  rm(storeNarray, tracersArray)
}

par(mfrow=c(3,2))
for(c in 1:nchaosVersions){
  hist(storeRatioResponse[c,,,50])
  hist(storeRatioResponse[c,,,1], col=myOrange, add=TRUE)
  hist(storeRatioResponse[c,,,151], col=myGrey_trans, add=TRUE)
}
# turn into df with response and explanatory vars
this_timeStep <- 150; half_timeStep <- 50 # not actually half, but matches the 50 years at which keystoneness and such were assessed in prev paper
allData2fit <- data.frame(matrix(NA, ncol=(nrankings+10), nrow=(nruns * nchaosVersions * ng)))
thisColnames <- c("ChaosRun", "relDiff", "half_relDiff", "relDiffLast20years", "ratioDiff", "half_ratioDiff", "ratioDiffLast20years", "maxDiff", "medDiff", "Informance", colnames(allTheRankings))
colnames(allData2fit)<- thisColnames
chaosLetters <- c("A", "B", "C", "D", "E", "F")
for(c in 1:nchaosVersions){
  for(r in 1:nruns){
    for(g in 1:ng){
      thisRow <- ng*nruns*(c-1)  + ng * (r-1) + g
      allData2fit$ChaosRun[thisRow] <- chaosLetters[c]; 
      allData2fit$relDiff[thisRow] <- storeRelResponse[c,r,g,this_timeStep]
      allData2fit$half_relDiff[thisRow] <- storeRelResponse[c,r,g,half_timeStep]
      allData2fit$relDiffLast20years[thisRow] <- mean((storeRelResponse[c,r,g,130:150]))
      allData2fit$ratioDiff[thisRow] <- storeRatioResponse[c,r,g,this_timeStep]
      allData2fit$half_ratioDiff[thisRow] <- storeRatioResponse[c,r,g,half_timeStep]
      allData2fit$ratioDiffLast20years[thisRow] <- mean((storeRatioResponse[c,r,g,130:150]))
      
      allData2fit$maxDiff[thisRow] <- max(abs(storeRelResponse[c,r,g,35:nts]))
      allData2fit$medDiff[thisRow] <- median(abs(storeRelResponse[c,r,g,35:nts]))
      thisCode <- as.character(thisGroupsDF$Code[g])
      allData2fit[thisRow,c(colnames(allTheRankings))] <- allTheRankings[allTheRankings$Code == thisCode,]
      allData2fit$Code[thisRow] <- thisCode
      thisNumCohorts <- thisGroupsDF$NumCohorts[g]
      if(thisNumCohorts==1){
        allData2fit$Keystone[thisRow] <- 0
      }
      thisInf <- rankDF$informPerformRank[rankDF$Code==thisCode]
      if(length(thisInf)==0){thisInf <- 16}
      allData2fit$Informance[thisRow] <- thisInf
    }
  }
}

# write out data2fit, so can call it in elsewhere
# take the X, X.1 columns out first
index <- grep("^X", colnames(allData2fit), invert = TRUE)
allData2fit <- allData2fit[,index]
write.csv(allData2fit, paste(DIR$'Tables',"Chaos_allData2fit.csv", sep=""), row.names = FALSE)


plotVar <- "half_ratioDiff"
thisMax <- max(allData2fit[c(plotVar)], na.rm=TRUE)
par(mfrow=c(2,2), mar=c(3,3,1,1))
for(g in 1:ng){
  thisCode <- thisGroupsDF$Code[g]
  hist(allData2fit[allData2fit$Code==thisCode,c(plotVar)], xlim=c(0, thisMax), ylab="", xlab="", main=""); mtext(thisCode, side=3, adj=0, line=-1)
  hist(allData2fit[allData2fit$Code==thisCode,c("ratioDiff")], add=TRUE, col=myOrange); 
  
}

test<- is.na(allData2fit$relDiff)
table(test) # some runs did not complete - take them out for this analyses as it  won't be clear which species group it was because of
test<- is.na(allData2fit$ratioDiff)
table(test)

data2fit <- allData2fit[!test,]
# turn B0 to tonnes . nope, already are
# data2fit$B0 <- data2fit$B0 * X_CN * mg_2_tonne
# change to 1000's though - easier later to plot, might as well do it here
data2fit$B0 <- data2fit$B0/1000

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


thisFactors <- c("Informance", "TL", "Keystone", "Response", "NumL1cons", "Lifespan", "propAdM", "propJuvM", "B0", "PropByTopPrey", "Linf")
# take out Linf as cuts out whales and such
thisFactors <- c("Informance", "TL", "Keystone", "Response", "NumL1cons", "Lifespan", "propAdM", "propJuvM", "B0", "PropByTopPrey","ChaosRun")

## keep only factors that work for biomass pool groups as well
bp_ind_factors <- c("TL",     "Keystone",   "NumL1cons",     "Lifespan" ,  "B0", "Informance", "Linf")

allInts <- c()
for(i in 1:length(bp_ind_factors)){
  for(j in 1:length(bp_ind_factors)){
    if(i != j){ 
      thisInt <- paste(bp_ind_factors[i],":", bp_ind_factors[j], sep="")
      checkIntRev <- paste(bp_ind_factors[j],":", bp_ind_factors[i], sep="")
      if(!(checkIntRev %in% allInts)){
        allInts <- c(allInts, thisInt)
      }
    }
  }
}
bp_factors<-c(bp_ind_factors, allInts)

test <- data2fit$relDiffLast20years==0
data2fit<-data2fit[!test,]

data2fit$logRatio <- unlist(lapply(data2fit$ratioDiffLast20years, FUN=function(x){log(x)}))
data2fit$cubeRat <- data2fit$ratioDiffLast20years^(1/3)
data2fit$logRatio <- log(data2fit$ratioDiffLast20years)
data2fit$cubeRel <- unlist(lapply(data2fit$relDiffLast20years, FUN=function(x){x^(1/3)}))
data2fit$logRel <- unlist(lapply(data2fit$relDiffLast20years, FUN=function(x){log((x)^2)}))
data2fit$logMax <- unlist(lapply(data2fit$maxDiff, FUN= function(x){log((x)+1)}))
data2fit$cubeMax <- unlist(lapply(data2fit$maxDiff, FUN= function(x){x^(1/3)}))

par(mfrow=c(2,2), mar=c(4,4,1,1))
plot(y=data2fit$cubeRel, x=data2fit$relDiff)
plot(y=data2fit$logRel, x=data2fit$relDiffLast20years)
# some of the predictor variables should be factors, not numbers - keystoneness and responsiveness are rankings
# actualFactors <- c("Informance", "Keystone", "Response", "ChaosRun"); naf <- length(actualFactors)
# actualFactors <- c("Informance", "ChaosRun"); naf <- length(actualFactors)
# actualFactors <- c("ChaosRun")
# for(f in 1:naf){
#   thisFactor <- actualFactors[f]
#   thisLevels <- sort(unique(data2fit[,c(thisFactor)]))
#   data2fit[,c(thisFactor)]<- factor(data2fit[,c(thisFactor)], levels=thisLevels)
# }
# check for any zeros in resp var
# bp_data2fit<-data2fit[!test,]
# table(test)

getBestLModel(thisVar="logRatio", thisFactors = bp_factors, cutoff=0.01, data2fit)
getBestLModel(thisVar="cubeRel", thisFactors = bp_factors, cutoff=0.01, data2fit)
getBestLModel(thisVar="logRel", thisFactors = bp_factors, cutoff=0.01, data2fit)
getBestLModel(thisVar="logMax", thisFactors = bp_factors, cutoff=0.01, data2fit)
getBestLModel(thisVar="cubeRat", thisFactors = bp_factors, cutoff=0.01, data2fit)
###########################################################################
# getBestLModel(thisVar="maxDiff", thisFactors = bp_factors, cutoff=0.01, data2fit)
# getBestLModel(thisVar="medDiff", thisFactors = thisFactors, cutoff=0.01, data2fit)
# getBestLModel(thisVar="relDiff", thisFactors = thisFactors, cutoff=0.01, data2fit)
# getBestLModel(thisVar="relDiffLast20years", thisFactors = thisFactors, cutoff=0.01, data2fit)
# getBestLModel(thisVar="ratioDiff", thisFactors = thisFactors, cutoff=0.01, data2fit)
# getBestLModel(thisVar="ratioDiffLast20years", thisFactors = thisFactors, cutoff=0.01, data2fit)
# 
# ##############################################################################################################################################################
###############################################################################
###############################################################################
###############################################################################
##############################################################################################################################################################
single_factors <- c("TL", "Keystone",  "Response", "Informance", "NumL1cons", "Lifespan", "propAdM",  "propJuvM",  "B0" ,   "PropByTopPrey", "Linf" ,"ChaosRun")
# no Linf as cuts out whales and such
single_factors <- c("TL", "Response", "Keystone", "Informance", "NumL1cons", "Lifespan", "propAdM",  "propJuvM",  "B0" ,   "PropByTopPrey","ChaosRun")
as_allInts <- c()
for(i in 1:length(single_factors)){
  for(j in 1:length(single_factors)){
    if(i != j){ 
      thisInt <- paste(single_factors[i],":", single_factors[j], sep="")
      checkIntRev <- paste(single_factors[j],":", single_factors[i], sep="")
      if(!(checkIntRev %in% as_allInts)){
        as_allInts <- c(as_allInts, thisInt)
      }
    }
  }
}
as_factors<-c(single_factors, as_allInts)

# have a look at only the age-structured groups
asCodes <- as.character(thisGroupsDF$Code[thisGroupsDF$NumCohorts>1])
as_index <- data2fit$Code %in% asCodes
as_data2fit <- data2fit[as_index,]
# do the runs look different..? NO
thisMax <- max(as_data2fit$ratioDiffLast20years)*100
par(mfrow=c(3,2))
for(c in 1:nchaosVersions){
  thisChaos <- levels(as_data2fit$ChaosRun)[c]
  testData <- 100*as_data2fit$relDiffLast20years[as_data2fit$ChaosRun==thisChaos]
  # hist(testData, xlim=c(0, thisMax), col="red", border="red")
  boxplot(testData, outline=TRUE)
}

# thisRespVar <-  "logRatio"
# this_k <- (-0.25)*min(as_data2fit[,c(thisRespVar)], na.rm=TRUE)
# check for all outliers
check_factors <- c(single_factors, c("relDiffLast20years", "ratioDiffLast20years"))
check_factors <- c( c("relDiffLast20years", "ratioDiffLast20years"))
check_factors <- c( "relDiffLast20years")

# check_factors <- c("logRatio", "logRelDiff")
nvars <- length(check_factors)
for(v in 1:nvars){
  testVar <- check_factors[v]
  # testVar <- "ratioDiffLast20years"
  thisLimits <- getPercentile(as_data2fit[,c(testVar)], 0.01)
  outIndex <- as_data2fit[,c(testVar)]< min(thisLimits) | as_data2fit[,c(testVar)]> max(thisLimits)
  x <- grep("TRUE", names(table(outIndex)))
  if(length(x)>0){
    numOutliers <- table(outIndex)[x]
    cat(paste("losing ", round(numOutliers/dim(as_data2fit)[1],2), "from ", check_factors[v]))
    as_data2fit <- as_data2fit[!outIndex,]
    cat("Num groups:", length(unique(as_data2fit$Code)),"\n")
  }
}
table(as_data2fit$Code)


# this_k <-0
# as_data2fit$logRatio <- unlist(lapply(as_data2fit$ratioDiffLast20years, FUN=function(x){log(x+this_k, base=exp(1))}))
# as_data2fit$logRelDiff <- unlist(lapply(100*as_data2fit$relDiffLast20years, FUN=function(x){log(x+this_k, base=exp(1))}))
# as_data2fit$cubeRatio <- unlist(lapply(as_data2fit$ratioDiffLast20years, FUN=function(x){x^{1/4}}))
# as_data2fit$cubeMax <- unlist(lapply(as_data2fit$maxDiff, FUN=function(x){x^{1/10}}))
as_data2fit$testCubeRel <- abs(100*as_data2fit$relDiffLast20years)^(1/6)
as_data2fit$testLogRel <- log(as_data2fit$relDiffLast20years^3)

# thisVar <- "logRatio"
thisVar <- "cubeRel"
thisVar <- "testCubeRel"
# thisVar <- "testLogRel"
# thisVar <- "cubeMax"
# thisVar <- "relDiff"
thisRespVar<-thisVar
plot(x=as_data2fit$relDiffLast20years, y=as_data2fit[,c(thisRespVar)])
#check for outliers
thisLimits <- getPercentile(as_data2fit[,c(thisRespVar)], k=0.01)
hist(as_data2fit[,c(thisRespVar)]); abline(v=thisLimits, col="red")
thisIndex <- as_data2fit[,c(thisRespVar)]<min(thisLimits) | as_data2fit[,c(thisRespVar)]>max(thisLimits)
hist(as_data2fit[!thisIndex, c(thisRespVar)])
this_data2fit <- as_data2fit[!thisIndex,]

thisChaosRuns <- c("A","B")
thisChaosRuns<-c("A", "B", "C", "D", "E", "F")
this_data2fit<-this_data2fit[(this_data2fit$ChaosRun %in% thisChaosRuns),]
this_data2fit$ChaosRun <- factor(this_data2fit$ChaosRun, levels=thisChaosRuns)

this_df <- getBestLModel(thisVar=thisRespVar, thisFactors = as_factors, cutoff=0.02, data2fit=this_data2fit)
thisFormula <- as.formula(paste(thisRespVar, "~", paste(this_df$Factor, collapse=" + "), sep=""))
this_df

this_model <- glm(thisFormula, data=this_data2fit)

testFit <- predict.glm(this_model, newdata = this_data2fit)
hist(testFit)
thisMin <- min(c(testFit, this_data2fit[,c(thisRespVar)]), na.rm=TRUE)
thisMax <- max(c(testFit, this_data2fit[,c(thisRespVar)]), na.rm=TRUE)
thisLine <- seq(thisMin, thisMax, length.out = 100)
# par(mfrow=c(3,2), mar=c(2,4,1,1))
plot(x=this_data2fit[,c(thisRespVar)], y=testFit, pch=20, ylim=c(thisMin, thisMax), xlim=c(thisMin, thisMax))
points(x=thisLine, y=thisLine, type="l", col="red", lwd=2, lty=2)
mtext(this_k, side=3, adj=0, font=2)
par(mfrow=c(2,2))
# plot(this_model)
#######################################################

thisPredVars <- getPredVars(this_df); npvars <- length(thisPredVars)
thisPredVarLengths <- rep(NA, npvars)
thisPredVarValues <- NULL
for(v in 1:npvars){
  thisVar <- thisPredVars[v]
  tempVarValues<-getVarValues(thisVar)
  ## fix the odd things - don't bother with zero biomass, lifespan, Linf, 
  # and keystonennes and responsiveness should just have unique values
  if(thisVar %in% c("B0", "Lifespan", "Linf", "NumL1cons")){
    tempVarValues <- tempVarValues[tempVarValues!=0]
  # } else if(thisVar %in% c("Keystone", "Informance")){
  } else if(thisVar %in% c("Informance", "Keystone","Response")){
    thisUnique <- sort(unique(as_data2fit[,c(thisVar)]))
    tempVarValues <- thisUnique
  }
  thisPredVarValues[[v]] <- tempVarValues
  thisPredVarLengths[v]<-length(thisPredVarValues[[v]])
}


thisSinglePredVars <- this_df$Factor[grep(":", this_df$Factor,invert = TRUE)]
if(length(thisSinglePredVars)>0){
par(mfrow=c(3,2))
for(v in 1:length(thisSinglePredVars)){
  thisVar <- thisSinglePredVars[v]
  thisVarValues <- sort(unique(this_data2fit[,c(thisVar)])); nvv <- length(thisVarValues)
    testMed<-data.frame(matrix(0, ncol=dim(this_model$data)[2], nrow=nvv)); colnames(testMed)<-colnames(this_model$data); 
     for(i in 1:length(single_factors)){
      fillVar <- single_factors[i]
      if( is.factor(this_data2fit[,fillVar])==FALSE){
        temp<- boxplot(this_model$data[,c(fillVar)], plot=FALSE)$stats
        fillValueMed <- temp[3]; 
        testMed[,c(fillVar)]<- fillValueMed; 
      } else{
        testMed[,c(fillVar)]<-names(sort(table(this_data2fit[,fillVar]),decreasing = TRUE)[1])
        testMed[,c(fillVar)] <- factor(testMed[,c(fillVar)], levels=sort(unique(as_data2fit[,c(fillVar)])))
      }
    }
    testMed[,c(thisVar)]<-thisVarValues;
    thisPred <- predict.glm(this_model,newdata = testMed)
    transPred <- unlist(lapply(thisPred, FUN= function(x){x^3}))
     thisYmin <- min(c(exp(c( thisPred)),0), na.rm=TRUE); thisYmax <- max(exp(c(thisPred)), na.rm=TRUE)
    plot(x=thisVarValues, y=transPred, type="l", xlab=thisVar, ylab="Relative change in biomass") 
    # points(exp(thisPredUQ), type="l", lty=4, col=myBlue); points(exp(thisPredLQ), type="l", lty=4, col="red")
}
}

# interaction effects as heat map
absColorRamp <- colorRampPalette(colors=c(myLightAqua, myAqua, myBlue,"midnightblue"))(101)

varDesc <- data.frame(matrix(NA, nrow=length(single_factors), ncol=2))
colnames(varDesc)<- c("Variable", "Decr")
varDesc$Variable <- single_factors
# varDesc$Decr <- c("Trophic level", "Keystone ranking", "Responsive ranking", "Level 1 connections", "Lifespan (years)", 
#                   "Additional adult mortality", "Additional juvenile mortality", expression(B[0]~"('000 tonnes)"),
#                   "Diet proportion by top prey", expression(L[infinity]~"(cm)"))

varDesc$Decr <- c("Trophic level", "Responsive ranking", "Keystone ranking", "Informance", "Level 1 connections", "Lifespan (years)", 
                  "Additional adult mortality", "Additional juvenile mortality", expression(B[0]~"('000 tonnes)"),
                  "Diet proportion by top prey","Simulation set")

x <- grep(":", this_df$Factor); intVars <- this_df$Factor[x]; nivars <- length(intVars)
lab_cex <- 0.8
thisMax <-10; thisMin <- 0

# pdf(paste(plotPath, "IntEffects_",thisRespVar ,".pdf", sep=""), height=4, width=7)
par(mfrow=c(2,2), mar=c(4,4,1,1), oma=c(0,0,0,6))
for(i in 1:nivars){
  thisIntVar <- intVars[i]; thisIntVars <- (unlist(str_split(thisIntVar,":")))
  thisValues1 <- as.double(thisPredVarValues[[grep(thisIntVars[1], thisPredVars)]])
  thisValues2 <- as.double(thisPredVarValues[[grep(thisIntVars[2], thisPredVars)]])
  if(is.factor(this_model$data[,c(thisIntVars[1])])){
    thisValues1<-thisPredVarValues[[grep(thisIntVars[1], thisPredVars)]]
  }
  if(is.factor(this_model$data[,c(thisIntVars[2])])){
    thisValues2<-thisPredVarValues[[grep(thisIntVars[2], thisPredVars)]]
  }
  test_df <- data.frame(matrix(0, ncol=dim(this_model$data)[2], nrow=(length(thisValues1) * length(thisValues2))))
  colnames(test_df)<- colnames(this_model$data)
   for(s in 1:length(single_factors)){
    fillVar <- single_factors[s]
    if( is.factor(this_data2fit[,fillVar])==FALSE){
      temp<- boxplot(this_model$data[,c(fillVar)], plot=FALSE)$stats
      fillValueMed <- temp[3]; 
      test_df[,c(fillVar)]<- fillValueMed; 
    } else{
      test_df[,c(fillVar)]<-names(sort(table(as_data2fit[,fillVar]),decreasing = TRUE)[1])
      test_df[,c(fillVar)] <- factor(test_df[,c(fillVar)], levels=sort(unique(as_data2fit[,c(fillVar)])))
    }
  }
  test_df[,c(thisIntVars[1])] <- sort(rep(thisValues1, length(thisValues2)));
  test_df[,c(thisIntVars[2])] <-  rep(thisValues2, length(thisValues1))
  if(is.factor(this_model$data[,c(thisIntVars[1])])){
    test_df[,c(thisIntVars[1])] <- factor(test_df[,c(thisIntVars[1])], levels=levels(this_model$data[,c(thisIntVars[1])]))
  }
  if(is.factor(this_model$data[,c(thisIntVars[2])])){
    test_df[,c(thisIntVars[2])] <- factor(test_df[,c(thisIntVars[2])], levels=levels(this_model$data[,c(thisIntVars[2])]))
  }
  # test_df$ChaosRun <- factor(test_df$ChaosRun, levels=thisChaosRuns)
  test_df$Predicted <- ((predict.glm(this_model, newdata = test_df)))^6
  # thisMax <- max(test_df$Predicted, na.rm=TRUE); thisMin <- min(test_df$Predicted, na.rm=TRUE)
  cat(paste("-- ",thisMax, thisMin, sep=", "))
  plotColour <- unlist(lapply(test_df$Predicted, getColor, log=TRUE))
  tempDF<-data.frame(cbind("x"=test_df[,c(thisIntVars[1])]),"y"=test_df[,c(thisIntVars[2])])
  x1<-sort(rep(seq(1,length(thisValues1)), length(thisValues2))); y1 <- rep(seq(1,length(thisValues2)), length(thisValues1))
  tempDF <- data.frame(cbind("x"=x1,"y"=y1))
  # pdf(paste(plotPath,"KemptonsQ_SSRsens.pdf",sep=""),height=4,width=10)
  plot(x=tempDF$x,y=tempDF$y,type="n",xlab="",ylab="",yaxt="n", xaxt="n")
  temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
  par(las=1)
  xx <- thisValues1; xxx <- seq(1, length(thisValues1))
  if(length(xx)>10){
    aIndex<-seq(1,length(xx), by=ceiling(length(xx)/10))
    axis1Labs  <- xx[aIndex]; axis1At <- xxx[aIndex]
  }else{
    axis1Labs  <- xx; axis1At <- xxx
  }
  xx <- thisValues2; xxx <- seq(1, length(thisValues2))
  if(length(xx)>10){
    aIndex<-seq(1,length(xx), by=ceiling(length(xx)/10))
    axis2Labs  <- xx[aIndex]; axis2At <- xxx[aIndex]
  } else{
    axis2Labs  <- xx; axis2At <- xxx
  }
  axis(at=axis1At, labels=axis1Labs, side=1)
  axis(at=axis2At, labels=axis2Labs, side=2)
  par(las=0)
  mtext(varDesc$Decr[varDesc$Variable==thisIntVars[1]], side=1, adj=0.5, line=2.5, cex=lab_cex); mtext(varDesc$Decr[varDesc$Variable==thisIntVars[2]], side=2, adj=0.5, line=2.5, cex=lab_cex)
  box()
}

predLegendValues <- pretty(seq(0,thisMax, length.out=5))
predLegendValues <- c(1,2,3,5,8)
legendValues <- predLegendValues # absoulte, and percentage


legendCols <- unlist(lapply(predLegendValues, getColor, log=TRUE))
# makeBlankPlot()
par(xpd=NA)
legend(legend=legendValues, col=legendCols, x="right", bty="n", pch=15, pt.cex=1.5, title="Absolute relative\ndifference (%)", inset=-0.5)
# 
# par(xpd=TRUE)
# par(mar=c(0,0,0,0))
# makeBlankPlot()
# legend(legend=chaosRunDescr, col="white", x="center", bty="n")
# dev.off()

thisPlotMax <- 0.1; thisMin <- 0

if(nivars>4){
  pdf(paste(plotPath, "IntEffects_",thisRespVar ,".pdf", sep=""), height=4, width=7)
  par(mfrow=c(2,3), mar=c(4,4,1.2,1))
  for(i in 1:nivars){
    thisIntVar <- intVars[i]; thisIntVars <- (unlist(str_split(thisIntVar,":")))
    thisValues1 <- as.double(thisPredVarValues[[grep(thisIntVars[1], thisPredVars)]])
    thisValues2 <- as.double(thisPredVarValues[[grep(thisIntVars[2], thisPredVars)]])
    if(is.factor(this_model$data[,c(thisIntVars[1])])){
      thisValues1<-thisPredVarValues[[grep(thisIntVars[1], thisPredVars)]]
    }
    if(is.factor(this_model$data[,c(thisIntVars[2])])){
      thisValues2<-thisPredVarValues[[grep(thisIntVars[2], thisPredVars)]]
    }
    test_df <- data.frame(matrix(0, ncol=dim(this_model$data)[2], nrow=(length(thisValues1) * length(thisValues2))))
    colnames(test_df)<- colnames(this_model$data)
    for(s in 1:length(single_factors)){
      fillVar <- single_factors[s]
      if( is.factor(this_data2fit[,fillVar])==FALSE){
        temp<- boxplot(this_model$data[,c(fillVar)], plot=FALSE)$stats
        fillValueMed <- temp[3]; 
        test_df[,c(fillVar)]<- fillValueMed; 
      } else{
        test_df[,c(fillVar)]<-names(sort(table(as_data2fit[,fillVar]),decreasing = TRUE)[1])
        test_df[,c(fillVar)] <- factor(test_df[,c(fillVar)], levels=sort(unique(as_data2fit[,c(fillVar)])))
      }
    }
    test_df[,c(thisIntVars[1])] <- sort(rep(thisValues1, length(thisValues2)));
    test_df[,c(thisIntVars[2])] <-  rep(thisValues2, length(thisValues1))
    if(is.factor(this_model$data[,c(thisIntVars[1])])){
      test_df[,c(thisIntVars[1])] <- factor(test_df[,c(thisIntVars[1])], levels=levels(this_model$data[,c(thisIntVars[1])]))
    }
    if(is.factor(this_model$data[,c(thisIntVars[2])])){
      test_df[,c(thisIntVars[2])] <- factor(test_df[,c(thisIntVars[2])], levels=levels(this_model$data[,c(thisIntVars[2])]))
    }
    # test_df$ChaosRun <- factor(test_df$ChaosRun, levels=thisChaosRuns)
    test_df$Predicted <- ((predict.glm(this_model, newdata = test_df)))^6
    # thisMax <- max(test_df$Predicted, na.rm=TRUE); thisMin <- min(test_df$Predicted, na.rm=TRUE)
    cat(paste("-- ",thisMax, thisMin, sep=", "))
    plotColour <- unlist(lapply(test_df$Predicted, getColor, log=TRUE))
    tempDF<-data.frame(cbind("x"=test_df[,c(thisIntVars[1])]),"y"=test_df[,c(thisIntVars[2])])
    x1<-sort(rep(seq(1,length(thisValues1)), length(thisValues2))); y1 <- rep(seq(1,length(thisValues2)), length(thisValues1))
    tempDF <- data.frame(cbind("x"=x1,"y"=y1))
    # pdf(paste(plotPath,"KemptonsQ_SSRsens.pdf",sep=""),height=4,width=10)
    plot(x=tempDF$x,y=tempDF$y,type="n",xlab="",ylab="",yaxt="n", xaxt="n")
    temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
    par(las=1)
    xx <- thisValues1; xxx <- seq(1, length(thisValues1))
    if(length(xx)>10){
      aIndex<-seq(1,length(xx), by=ceiling(length(xx)/10))
      axis1Labs  <- xx[aIndex]; axis1At <- xxx[aIndex]
    }else{
      axis1Labs  <- xx; axis1At <- xxx
    }
    xx <- thisValues2; xxx <- seq(1, length(thisValues2))
    if(length(xx)>10){
      aIndex<-seq(1,length(xx), by=ceiling(length(xx)/10))
      axis2Labs  <- xx[aIndex]; axis2At <- xxx[aIndex]
    } else{
      axis2Labs  <- xx; axis2At <- xxx
    }
    axis(at=axis1At, labels=axis1Labs, side=1)
    axis(at=axis2At, labels=axis2Labs, side=2)
    par(las=0)
    mtext(varDesc$Decr[varDesc$Variable==thisIntVars[1]], side=1, adj=0.5, line=2.5, cex=lab_cex); mtext(varDesc$Decr[varDesc$Variable==thisIntVars[2]], side=2, adj=0.5, line=2.5, cex=lab_cex)
    box()
  }
  
  predLegendValues <- pretty(seq(0,thisMax, length.out=5))
  predLegendValues <- c(1,2,3,5,8)
  legendValues <- predLegendValues # absoulte, and percentage
  
  
  legendCols <- unlist(lapply(predLegendValues, getColor, log=TRUE))
  # makeBlankPlot()
  # par(xpd=NA)
  # legend(legend=legendValues, col=legendCols, x="right", bty="n", pch=15, pt.cex=1.5, title="Absolute relative\ndifference (%)", inset=-0.5)
 
  
  makeBlankPlot()
  par(xpd=NA)
  legend(legend=legendValues, col=legendCols, x="center", bty="n", pch=15, pt.cex=1.5, title="Absolute relative\ndifference (%)", inset=-0.6, horiz=FALSE)
   # 
  # par(xpd=TRUE)
  # par(mar=c(0,0,0,0))
  # makeBlankPlot()
  # legend(legend=chaosRunDescr, col="white", x="center", bty="n")
  dev.off()
  
  
}

# use this model to give a prediction for each age-structured species group, and compare this to the actual range for each
thisRankings <- allTheRankings[allTheRankings$Code %in% asCodes,]
thisRankings$Predicted <- (predict.glm(this_model, newdata = thisRankings))^3
thisXmin <- 0; thisXmax <- max(as_data2fit$ratioDiffLast20years, na.rm=TRUE)

thisMax <- max(as_data2fit[,thisRespVar], na.rm=TRUE)
thisMin <- min(as_data2fit[,thisRespVar], na.rm=TRUE)

par(mfrow=c(3,3), mar=c(3,3,1,1), lend=1)
for(g in 1:ng){
  testCode <- as.character(groupsDF$Code[g])
  thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
       testData <- as_data2fit[as_data2fit$Code==testCode,]
      if(dim(testData)[1]>0){
        thisPredicted <- predict.glm(this_model, newdata = testData[1,])
        testPred <- exp(thisPredicted)
        thisHist <- hist(testData[,thisRespVar], plot=FALSE)
        # thisMax <- max(c(thisPredicted, thisHist$mids), na.rm=TRUE)
        # thisMin <- min(c(thisPredicted, thisHist$mids), na.rm=TRUE)
        # 
        plot(x=thisHist$mids, y=thisHist$density,xlim=c(thisMin, thisMax), type="h", lwd=10, col=myGrey)
        abline(v=thisPredicted, col="red", lwd=5)
        mtext(paste(testCode, sep=""), side=3, adj=0)
    }
    
  }
}
####################################################################
# mtext(expression(L[infinity]~"(cm)"), side=1, line=2)

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

########################################
## bp and as combined

this_k <- 0
bp_data2fit$logRatio <- unlist(lapply(100*bp_data2fit$ratioDiffLast20years, FUN=function(x){log(x+this_k, base=exp(1))}))
bp_data2fit$cubeRatio <- unlist(lapply(bp_data2fit$ratioDiffLast20years, FUN=function(x){(x)^(1/3)}))

bp_data2fit$logRelDiff <- unlist(lapply(bp_data2fit$relDiffLast20years, FUN=function(x){log(x+this_k, base=exp(1))}))
bp_data2fit$cubeRoot <- unlist(lapply(bp_data2fit$relDiffLast20years, FUN=function(x){(x)^(1/3)}))
bp_data2fit$sqrRoot <- unlist(lapply(bp_data2fit$relDiffLast20years, FUN=function(x){(x+0.5)^(1/this_k)}))



# thisVar <- "cubeRoot"
# thisVar <- "logRelDiff"
# thisVar <- "relDiffLast20years"
# thisVar <- "sqrRoot"
thisVar <- "logRatio"
# thisVar <- "cubeRatio"
## locate outlier
hist(bp_data2fit[,c(thisVar)])
# any outliers?
thisLimits <- boxplot(bp_data2fit[,c(thisVar)], plot=FALSE)$stats[c(1,5)]
outIndex <- bp_data2fit[,c(thisVar)]< min(thisLimits) | bp_data2fit[,c(thisVar)]> max(thisLimits)
table(outIndex)
thisBP_data2fit <- bp_data2fit[!outIndex,]

# check for outliers in all explanatory variables
single_bp_factors <- unique(unlist(str_split(gsub(":", " ", bp_factors)," ")))
nvars <- length(single_bp_factors)
for(v in 1:nvars){
  testVar <- single_bp_factors[v]
  thisLimits <- boxplot(bp_data2fit[,c(testVar)], plot=FALSE)$stats[c(1,5)]
  outIndex <- bp_data2fit[,c(testVar)]< min(thisLimits) | bp_data2fit[,c(testVar)]> max(thisLimits)
  x <- grep("TRUE", names(table(outIndex)))
  if(length(x)>0){
    numOutliers <- table(outIndex)[x]
    cat(paste("losing ", round(numOutliers/dim(bp_data2fit)[1],2), "from ", single_bp_factors[v]))
    thisBP_data2fit <- thisBP_data2fit[!outIndex,]
  }
}



this_df <- getBestLModel(thisVar=thisVar, thisFactors = bp_factors, cutoff=0.01, data2fit=thisBP_data2fit)
thisFormula <- as.formula(paste(thisVar, "~", paste(this_df$Factor, collapse=" + "), sep=""))
this_df

this_model <- glm(thisFormula, data=thisBP_data2fit)
# # check for any gaps in the explanatory vars - shouldn't need this
# xx <- unique(unlist(str_split(gsub(":"," ",thisFormula)," "))); thisExpVars<-xx[!(xx %in% c("~", "+"))]; nexvars <- length(thisExpVars)
# for(e in 1:nexvars){
#   index <- !is.na(thisBP_data2fit[,c(thisExpVars[e])])
#   if(length(as.double(table(index)[FALSE]))>0){
#     cat(paste("losing ", as.double(table(index)[FALSE]), "from ", thisExpVars[e]))
#     thisBP_data2fit <- thisBP_data2fit[index,]
#   }
# }


testFit <- predict.glm(this_model, newdata = thisBP_data2fit)
hist(testFit)
thisMin <- min(c(testFit, thisBP_data2fit[,c(thisVar)]), na.rm=TRUE)
thisMax <- max(c(testFit, thisBP_data2fit[,c(thisVar)]), na.rm=TRUE)
thisLine <- seq(thisMin, thisMax, length.out = 100)
# par(mfrow=c(3,2), mar=c(2,4,1,1))
plot(x=thisBP_data2fit[,c(thisVar)], y=testFit, pch=20, ylim=c(thisMin, thisMax), xlim=c(thisMin, thisMax))
points(x=thisLine, y=thisLine, type="l", col="red", lwd=2, lty=2)
mtext(this_k, side=3, adj=0, font=2)

plot(x=exp(thisBP_data2fit[,c(thisVar)]), y=exp(testFit), pch=20)

plot(this_model)



# #####################################################################
# # if don't use the fishing ones:
# getIndexfromC <- function(c){
#   index <- c()
#   for(r in 1:nruns){
#     for(g in 1:ng){
#       this_i <- ng * nruns * (c-1) + ng * (r-1) + g
#       index <- c(index, this_i)
#     }
#   }
#   return(index)
# }
# # c1index <- getIndexfromC(c=1); c3index <- getIndexfromC(c=3)
# chaosDescIndex <- c(1,3); chaosRunIndex <- unlist(lapply(chaosDescIndex, getIndexfromC))
# 
# this_data2fit <- allData2fit[chaosRunIndex,]
# index <- is.na(this_data2fit$relDiff)
# 
# this_data2fit<-this_data2fit[!index,]
