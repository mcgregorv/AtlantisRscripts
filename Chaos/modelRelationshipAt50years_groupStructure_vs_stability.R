# analyse stability using tracers rel to base runs - tracers were read in and stored in readInChaosRunsAndSaveForAnalyses.R
# use allRankings as possible explanatory variables - this table was created in network analyses
# have 2 main sets of runs so far - SampleA and SampleKEY - both have fished and non-fished version
# also including the all up or all down samples - so 3 sets of 2 (in total 6) with half including fishing
# Use DChaosBASE and DChaosFISH to compare
basePath <-paste(DIR$'Base',"AtlantisModels\\chaos\\",sep="")
thisChaosVersion <- "SampleA";
thisChaosVersion <- "FISHSampleA";
# thisChaosVersion <- "SampleKEY";
thisChaosVersion <- "FISHSampleKEY";
thisChaosVersion <- "SampleKEYPlusBP"
# nChaosRuns <- 35
source(paste(DIR$'General functions',"chaosFunctions.R", sep=""))
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
storeRatioResponse <- storeRelResponse;
storeAllTracersFromALlRuns <- storeRelResponse
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
      points(theseChaosTracers[r,], type="l", col=myGreen)
    }
  }
  storeRelResponse[c,1:this_nruns,,]<-allRelTracers
  storeRatioResponse[c,1:this_nruns,,]<-allRatioTracers
  storeAllTracersFromALlRuns[c,1:this_nruns,,] <- storeNarray
  
  rm(storeNarray, tracersArray)
}

# write these out
# save(list=c("storeRelResponse", "storeRatioResponse", "storeAllTracersFromALlRuns"), file=paste(basePath, "AllTracersFromAllCHAOSRuns", sep=""))



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

plotVar <- "half_ratioDiff"
thisMax <- max(allData2fit[c(plotVar)], na.rm=TRUE)
par(mfrow=c(2,2), mar=c(3,3,1,1))
for(g in 1:ng){
  thisCode <- thisGroupsDF$Code[g]
  hist(allData2fit[allData2fit$Code==thisCode,c(plotVar)], xlim=c(0, thisMax), ylab="", xlab="", main=""); mtext(thisCode, side=3, adj=0, line=-1)
  hist(allData2fit[allData2fit$Code==thisCode,c("ratioDiff")], add=TRUE, col=myOrange); 
  
}


# write out data2fit, so can call it in elsewhere
# take the X, X.1 columns out first
index <- grep("^X", colnames(allData2fit), invert = TRUE)
allData2fit <- allData2fit[,index]
# write.csv(allData2fit, paste(DIR$'Tables',"Chaos_allData2fit_inclHalfWay.csv", sep=""), row.names = FALSE)

test<- is.na(allData2fit$half_relDiff)
table(test) # some runs did not complete - take them out for this analyses as it  won't be clear which species group it was because of
test<- is.na(allData2fit$half_ratioDiff)
table(test)

data2fit <- allData2fit[!test,]
# turn B0 to tonnes . nope, already are
# data2fit$B0 <- data2fit$B0 * X_CN * mg_2_tonne
# change to 1000's though - easier later to plot, might as well do it here
data2fit$B0 <- data2fit$B0/1000

thisFactors <- c("Informance", "TL", "Keystone", "Response", "NumL1cons", "Lifespan", "propAdM", "propJuvM", "B0", "PropByTopPrey","ChaosRun", "Linf")
# take out Linf as cuts out whales and such
thisFactors <- c("Informance", "TL", "Keystone", "Response", "NumL1cons", "Lifespan", "propAdM", "propJuvM", "B0", "PropByTopPrey","ChaosRun")

varDesc <- data.frame(matrix(NA, nrow=length(thisFactors), ncol=2))
colnames(varDesc)<- c("Variable", "Decr")
varDesc$Variable <- thisFactors
# varDesc$Decr <- c("Trophic level", "Keystone ranking", "Responsive ranking", "Level 1 connections", "Lifespan (years)", 
#                   "Additional adult mortality", "Additional juvenile mortality", expression(B[0]~"('000 tonnes)"),
#                   "Diet proportion by top prey", expression(L[infinity]~"(cm)"))
thisDescriptions <- c("Informance", "Trophic level",  "Keystone ranking", "Responsive ranking", "Level 1 connections", "Lifespan (years)", 
                      "Additional adult mortality", "Additional juvenile mortality", "", "Diet proportion by top prey","Simulation set")

varDesc$Decr <- thisDescriptions
# varDesc$Decr[9]<- expression(B[0]~"('000 tonnes)")

# thisDescriptions <- c("Informance", "Trophic level",  "Keystone ranking", "Responsive ranking", "Level 1 connections", "Lifespan (years)", 
#   "Additional adult mortality", "Additional juvenile mortality", expression(B[0]~"('000 tonnes)"),
#   "Diet proportion by top prey","Simulation set")


## keep only factors that work for biomass pool groups as well
bp_ind_factors <- c("TL",     "NumL1cons",     "B0", "PropByTopPrey")


allInts <- getIntVars(bp_ind_factors); bp_factors <- c(bp_ind_factors, allInts)

data2fit$log_halfRelDiff <- unlist(lapply(data2fit$half_relDiff, FUN=function(x){log(abs(x)+1)}))
data2fit$cube_halfRelDiff <- unlist(lapply(data2fit$half_relDiff, FUN=function(x){(abs(x))^(1/6)}))

data2fit$cubeRel <- unlist(lapply(data2fit$relDiffLast20years, FUN=function(x){abs(x)^(1/6)}))
data2fit$logRel <- unlist(lapply(data2fit$relDiffLast20years, FUN=function(x){log((abs(x)+1))}))

justBPdata <- data2fit[!(data2fit$Code %in% asCodes),]

# first check regression tree
thisRespVar<-"cube_halfRelDiff"
singleFormula <- as.formula(paste(thisRespVar, "~", paste(bp_ind_factors, collapse=" + "), sep=""))
treeModel <- tree(formula=singleFormula, data=justBPdata)
plot(treeModel)
text(treeModel)
treeFactors <- as.character(summary(treeModel)$used)
testFactors <-c(treeFactors, getIntVars(treeFactors))
halfBPdf <- getBestLModel(thisVar=thisRespVar, thisFactors = bp_factors, cutoff=0.02, data2fit=justBPdata)
halfBPFormula <- as.formula(paste(thisRespVar, "~", paste(halfBPdf$Factor, collapse=" + "), sep=""))
halfBPdf

halfBPModel <- glm(halfBPFormula, data=justBPdata)

testFit <- predict.glm(halfBPModel, newdata = justBPdata)
hist(testFit); hist(justBPdata[,c(thisRespVar)])
thisMin <- min(c(testFit, justBPdata[,c(thisRespVar)]), na.rm=TRUE)
thisMax <- max(c(testFit, justBPdata[,c(thisRespVar)]), na.rm=TRUE)
thisLine <- seq(thisMin, thisMax, length.out = 100)
plot(x=justBPdata[,c(thisRespVar)], y=testFit, pch=20, ylim=c(thisMin, thisMax), xlim=c(thisMin, thisMax))
points(x=thisLine, y=thisLine, type="l", col="red", lwd=2, lty=2)
par(mfrow=c(2,2))

# at end of 150 year run
thisRespVar<-"cubeRel"
singleFormula <- as.formula(paste(thisRespVar, "~", paste(bp_ind_factors, collapse=" + "), sep=""))
treeModel <- tree(formula=singleFormula, data=justBPdata)
plot(treeModel)
text(treeModel)
treeFactors <- as.character(summary(treeModel)$used)
testFactors <-c(treeFactors, getIntVars(treeFactors))

testFormula <- as.formula(paste(thisRespVar,"~",paste(treeFactors,collapse="+"), sep=""))
testModel <- glm(testFormula, data=justBPdata)
thisRsq <- 1-(summary(testModel)$"deviance"/summary(testModel)$"null.deviance")

finalBPdf <- getBestLModel(thisVar=thisRespVar, thisFactors = bp_factors, cutoff=0.02, data2fit=justBPdata)
finalBPFormula <- as.formula(paste(thisRespVar, "~", paste(finalBPdf$Factor, collapse=" + "), sep=""))
finalBPdf

finalBPModel <- glm(finalBPFormula, data=justBPdata)
summary(finalBPModel)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       7.807e-01  7.006e-03  111.43   <2e-16 ***
#   TL:PropByTopPrey -7.510e-02  4.520e-03  -16.61   <2e-16 ***
#   NumL1cons:B0     -2.888e-07  1.643e-08  -17.58   <2e-16 ***
#   TL:B0             1.195e-06  9.339e-08   12.79   <2e-16 ***
  
testFit <- predict.glm(finalBPModel, newdata = justBPdata)
hist(testFit); hist(justBPdata[,c(thisRespVar)])
thisMin <- min(c(testFit, justBPdata[,c(thisRespVar)]), na.rm=TRUE)
thisMax <- max(c(testFit, justBPdata[,c(thisRespVar)]), na.rm=TRUE)
thisLine <- seq(thisMin, thisMax, length.out = 100)
plot(x=justBPdata[,c(thisRespVar)], y=testFit, pch=20, ylim=c(thisMin, thisMax), xlim=c(thisMin, thisMax))
points(x=thisLine, y=thisLine, type="l", col="red", lwd=2, lty=2)
par(mfrow=c(2,2))

## whole dataset, but only bp predictors
thisRespVar<-"cubeRel"
singleFormula <- as.formula(paste(thisRespVar, "~", paste(bp_ind_factors, collapse=" + "), sep=""))
treeModel <- tree(formula=singleFormula, data=data2fit)
plot(treeModel)
text(treeModel)
treeVars <- as.character(summary(treeModel)$used); treeIntVars <- getIntVars(treeVars)

  
finalBPAlldf <- getBestLModel(thisVar=thisRespVar, thisFactors = bp_factors, cutoff=0.02, data2fit=data2fit)
finalBPAllFormula <- as.formula(paste(thisRespVar, "~", paste(finalBPAlldf$Factor, collapse=" + "), sep=""))
finalBPAlldf

finalBPAllModel <- glm(finalBPAllFormula, data=data2fit)
testFit <- predict.glm(finalBPAllModel, newdata = data2fit)
hist(testFit); hist(data2fit[,c(thisRespVar)])
thisMin <- min(c(testFit, data2fit[,c(thisRespVar)]), na.rm=TRUE)
thisMax <- max(c(testFit, data2fit[,c(thisRespVar)]), na.rm=TRUE)
thisLine <- seq(thisMin, thisMax, length.out = 100)
plot(x=data2fit[,c(thisRespVar)], y=testFit, pch=20, ylim=c(thisMin, thisMax), xlim=c(thisMin, thisMax))
points(x=thisLine, y=thisLine, type="l", col="red", lwd=2, lty=2)
par(mfrow=c(2,2))
summary(finalBPAllModel)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   7.114e-01  4.061e-03  175.18   <2e-16 ***
#   TL           -7.287e-02  1.809e-03  -40.28   <2e-16 ***
#   TL:NumL1cons  1.501e-03  8.299e-05   18.08   <2e-16 ***
  
thisRespVar<-"cube_halfRelDiff"

halfBPAlldf <- getBestLModel(thisVar=thisRespVar, thisFactors = bp_factors, cutoff=0.02, data2fit=data2fit)
halfBPAllFormula <- as.formula(paste(thisRespVar, "~", paste(halfBPAlldf$Factor, collapse=" + "), sep=""))
halfBPAlldf

halfBPAllModel <- glm(halfBPAllFormula, data=data2fit)
testFit <- predict.glm(halfBPAllModel, newdata = data2fit)
hist(testFit); hist(data2fit[,c(thisRespVar)])
thisMin <- min(c(testFit, data2fit[,c(thisRespVar)]), na.rm=TRUE)
thisMax <- max(c(testFit, data2fit[,c(thisRespVar)]), na.rm=TRUE)
thisLine <- seq(thisMin, thisMax, length.out = 100)
plot(x=data2fit[,c(thisRespVar)], y=testFit, pch=20, ylim=c(thisMin, thisMax), xlim=c(thisMin, thisMax))
points(x=thisLine, y=thisLine, type="l", col="red", lwd=2, lty=2)
par(mfrow=c(2,2))
summary(halfBPAllModel)
## all terms are significant
# summary(halfBPAllModel)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   7.529e-01  3.693e-03  203.88   <2e-16 ***
#   TL           -6.569e-02  1.543e-03  -42.57   <2e-16 ***
#   NumL1cons:B0 -6.291e-08  3.360e-09  -18.72   <2e-16 ***
#   TL:NumL1cons  1.058e-03  6.875e-05   15.38   <2e-16 ***
#   ---
#
##################################### plot int effects

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

lab_cex <- 1
x <- grep(":", halfBPdf$Factor); intVars <- halfBPdf$Factor[x]; nivars <- length(intVars)
plotLetters<-c("A", "B", "C")
pdf(paste(plotPath,"InteractionEffects_BPonlydata50years.pdf", sep=""), height=3, width=12)
par(mfrow=c(1,3), mar=c(4,5.5,1.5,1), oma=c(0,0,0,7))
for(i in 1:nivars){
  thisIntVar <- intVars[i]; 
  plotIntEffect(thisIntVar = intVars[i], this_model=halfBPModel, thisMax=0.3, updateMax = FALSE)
  mtext(plotLetters[i], side=3, adj=0, line=0.1, font=2)
}
plotSingleEffect(thisVar="NumL1cons", this_model=halfBPModel)
mtext(plotLetters[3], side=3, adj=0, line=0.1, font=2)
predLegendValues <- pretty(seq(0,0.3, length.out=5))
legendValues <- predLegendValues # absoulte, and percentage


legendCols <- unlist(lapply(predLegendValues, getColor, log=FALSE))
# makeBlankPlot()
par(xpd=NA)
legend(legend=legendValues, col=legendCols, x="right", bty="n", pch=15, pt.cex=1.5, title="Absolute relative\ndifference (%)", inset=-0.35)

dev.off()

# do the end-of-run version
x <- grep(":", finalBPdf$Factor); intVars <- finalBPdf$Factor[x]; nivars <- length(intVars)
plotLetters<-c("A", "B", "C")
pdf(paste(plotPath,"InteractionEffects_BPonlydata150years.pdf", sep=""), height=3, width=12)
par(mfrow=c(1,3), mar=c(4,5.5,1.5,1), oma=c(0,0,0,7))
for(i in 1:nivars){
  thisIntVar <- intVars[i]; 
  plotIntEffect(thisIntVar = intVars[i], this_model=finalBPModel, thisMax=0.3, updateMax = FALSE)
  mtext(plotLetters[i], side=3, adj=0, line=0.1, font=2)
}
predLegendValues <- pretty(seq(0,0.3, length.out=5))
legendValues <- predLegendValues # absoulte, and percentage


legendCols <- unlist(lapply(predLegendValues, getColor, log=FALSE))
# makeBlankPlot()
par(xpd=NA)
legend(legend=legendValues, col=legendCols, x="right", bty="n", pch=15, pt.cex=1.5, title="Absolute relative\ndifference (%)", inset=-0.35)

dev.off()

#############################################
## same figures, using all data
##########################################
x <- grep(":", halfBPAlldf$Factor); intVars <- halfBPAlldf$Factor[x]; nivars <- length(intVars)
plotLetters<-c("A", "B", "C")
pdf(paste(plotPath,"InteractionEffects_ALLdata50years.pdf", sep=""), height=3, width=12)
par(mfrow=c(1,3), mar=c(4,5.5,1.5,1), oma=c(0,0,0,7))
plotSingleEffect(thisVar="TL", this_model=halfBPAllModel)
mtext(plotLetters[1], side=3, adj=0, line=0.1, font=2)
for(i in 1:nivars){
  thisIntVar <- intVars[(i)]; 
  plotIntEffect(thisIntVar = intVars[i], this_model=halfBPAllModel, thisMax=0.3, updateMax = FALSE)
  mtext(plotLetters[(i+1)], side=3, adj=0, line=0.1, font=2)
}
predLegendValues <- pretty(seq(0,0.3, length.out=5))
legendValues <- predLegendValues # absoulte, and percentage
legendCols <- unlist(lapply(predLegendValues, getColor, log=FALSE))
# makeBlankPlot()
par(xpd=NA)
legend(legend=legendValues, col=legendCols, x="right", bty="n", pch=15, pt.cex=1.5, title="Absolute relative\ndifference (%)", inset=-0.35)

dev.off()

# do the end-of-run version
x <- grep(":", finalBPAlldf$Factor); intVars <- finalBPAlldf$Factor[x]; nivars <- length(intVars)
plotLetters<-c("A", "B", "C")
pdf(paste(plotPath,"InteractionEffects_ALLdata150years.pdf", sep=""), height=4, width=12)
par(mfrow=c(1,2), mar=c(4,5.5,1.5,1), oma=c(0,0,0,10))
plotSingleEffect(thisVar="TL", this_model=finalBPAllModel)
mtext(plotLetters[1], side=3, adj=0, line=0.1, font=2)
for(i in 1:nivars){
  thisIntVar <- intVars[i]; 
  plotIntEffect(thisIntVar = intVars[i], this_model=finalBPAllModel, thisMax=0.3, updateMax = FALSE)
  mtext(plotLetters[(i+1)], side=3, adj=0, line=0.1, font=2)
}
predLegendValues <- pretty(seq(0,0.3, length.out=5))
legendValues <- predLegendValues # absoulte, and percentage
legendCols <- unlist(lapply(predLegendValues, getColor, log=FALSE))
# makeBlankPlot()
par(xpd=NA)
legend(legend=legendValues, col=legendCols, x="right", bty="n", pch=15, pt.cex=1.5, title="Absolute relative\ndifference (%)", inset=-0.5)

dev.off()

# ##############################################################################################################################################################
###############################################################################
###############################################################################

#############################################
## same figures, using age-structured data
##########################################
# prep model and data first
as_data2fit <- data2fit[data2fit$Code %in% asCodes,]
single_factors <- c("TL", "Keystone",  "Response", "Informance", "NumL1cons", "Lifespan", "propAdM",  "propJuvM",  "B0" ,   "PropByTopPrey", "Linf" ,"ChaosRun")
# no Linf as cuts out whales and such
# single_factors <- c("TL", "Response", "Keystone", "Informance", "NumL1cons", "Lifespan", "propAdM",  "propJuvM",  "B0" ,   "PropByTopPrey","ChaosRun")
as_intFactors <- getIntVars(single_factors); as_factors <- c(single_factors, as_intFactors)
as_factors<-c(single_factors, as_allInts)
as_data2fit$ChaosRun <- factor(as_data2fit$ChaosRun, levels=c("A", "B", "C", "D","E", "F"))
# do the runs look different..? yes at 50 years, no at 151 years
thisMax <- max(abs(as_data2fit$half_relDiff), na.rm=TRUE)*100
thisMax <- 30
plot(x=1:nchaosVersions, y=rep(0, nchaosVersions), type="n", ylim=c(0, thisMax))
for(c in 1:nchaosVersions){
  thisChaos <- levels(as_data2fit$ChaosRun)[c]
  endTestData <- 100 * abs(as_data2fit$half_relDiff[as_data2fit$ChaosRun==thisChaos])
  testData <- 100*abs(as_data2fit$relDiff[as_data2fit$ChaosRun==thisChaos])
  boxplot(endTestData, outline=FALSE, add=TRUE, at=c)
  boxplot(testData, outline=FALSE, add=TRUE, at=c, col=myOrange_trans, border=myOrange)
}
## fit the model - half first

thisRespVar<-"cube_halfRelDiff"

halfASdf <- getBestLModel(thisVar=thisRespVar, thisFactors = as_factors, cutoff=0.02, data2fit=as_data2fit)
halfASFormula <- as.formula(paste(thisRespVar, "~", paste(halfASdf$Factor, collapse=" + "), sep=""))
halfASdf

halfASModel <- glm(halfASFormula, data=as_data2fit)
testFit <- predict.glm(halfASModel, newdata = as_data2fit)
hist(testFit); hist(as_data2fit[,c(thisRespVar)])
thisMin <- min(c(testFit, as_data2fit[,c(thisRespVar)]), na.rm=TRUE)
thisMax <- max(c(testFit, as_data2fit[,c(thisRespVar)]), na.rm=TRUE)
thisLine <- seq(thisMin, thisMax, length.out = 100)
plot(x=as_data2fit[,c(thisRespVar)], y=testFit, pch=20, ylim=c(thisMin, thisMax), xlim=c(thisMin, thisMax))
points(x=thisLine, y=thisLine, type="l", col="red", lwd=2, lty=2)
par(mfrow=c(2,2))
summary(halfASModel)

## final - at 150 years
thisRespVar<-"cubeRel"

finalASdf <- getBestLModel(thisVar=thisRespVar, thisFactors = as_factors, cutoff=0.02, data2fit=as_data2fit)
finalASFormula <- as.formula(paste(thisRespVar, "~", paste(finalASdf$Factor, collapse=" + "), sep=""))
finalASdf

finalASModel <- glm(finalASFormula, data=as_data2fit)
testFit <- predict.glm(finalASModel, newdata = as_data2fit)
hist(testFit); hist(as_data2fit[,c(thisRespVar)])
thisMin <- min(c(testFit, as_data2fit[,c(thisRespVar)]), na.rm=TRUE)
thisMax <- max(c(testFit, as_data2fit[,c(thisRespVar)]), na.rm=TRUE)
thisLine <- seq(thisMin, thisMax, length.out = 100)
plot(x=as_data2fit[,c(thisRespVar)], y=testFit, pch=20, ylim=c(thisMin, thisMax), xlim=c(thisMin, thisMax))
points(x=thisLine, y=thisLine, type="l", col="red", lwd=2, lty=2)
par(mfrow=c(2,2))
summary(finalASModel)
## all significant
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          0.6136494  0.0110219  55.675  < 2e-16 ***
#   NumL1cons:ChaosRunA  0.0054007  0.0011189   4.827 1.43e-06 ***
#   NumL1cons:ChaosRunB  0.0052928  0.0011189   4.730 2.31e-06 ***
#   NumL1cons:ChaosRunC  0.0070742  0.0006293  11.242  < 2e-16 ***
#   NumL1cons:ChaosRunD  0.0073447  0.0005879  12.493  < 2e-16 ***
#   NumL1cons:ChaosRunE  0.0098800  0.0005790  17.065  < 2e-16 ***
#   NumL1cons:ChaosRunF  0.0098277  0.0005879  16.716  < 2e-16 ***
#   ChaosRunA:TL        -0.0403683  0.0049395  -8.173 3.88e-16 ***
#   ChaosRunB:TL        -0.0478087  0.0049395  -9.679  < 2e-16 ***
#   ChaosRunC:TL        -0.0351743  0.0030304 -11.607  < 2e-16 ***
#   ChaosRunD:TL        -0.0456226  0.0028789 -15.847  < 2e-16 ***
#   ChaosRunE:TL        -0.0547936  0.0028466 -19.249  < 2e-16 ***
#   ChaosRunF:TL        -0.0624926  0.0028789 -21.707  < 2e-16 ***
#   Keystone:propAdM    -0.0026063  0.0002136 -12.202  < 2e-16 ***

####################################################
#### now the effects plots (age-structured) ############################
##################################################
x <- grep(":", halfASdf$Factor); intVars <- halfASdf$Factor[x]; nivars <- length(intVars)
plotLetters<-c("A", "B", "C")
pdf(paste(plotPath,"InteractionEffects_ASdata50years.pdf", sep=""), height=3, width=12)
par(mfrow=c(1,3), mar=c(4,5.5,1.5,1), oma=c(0,0,0,7))
plotSingleEffect(thisVar="TL", this_model=halfASModel)
mtext(plotLetters[1], side=3, adj=0, line=0.1, font=2)
for(i in 1:nivars){
  thisIntVar <- intVars[(i)]; 
  plotIntEffect(thisIntVar = intVars[i], this_model=halfASModel, thisMax=0.3, updateMax = FALSE)
  mtext(plotLetters[(i+1)], side=3, adj=0, line=0.1, font=2)
}
predLegendValues <- pretty(seq(0,0.3, length.out=5))
legendValues <- predLegendValues # absoulte, and percentage
legendCols <- unlist(lapply(predLegendValues, getColor, log=FALSE))
# makeBlankPlot()
par(xpd=NA)
legend(legend=legendValues, col=legendCols, x="right", bty="n", pch=15, pt.cex=1.5, title="Absolute relative\ndifference (%)", inset=-0.35)

dev.off()

# do the end-of-run version
x <- grep(":", finalASdf$Factor); intVars <- finalASdf$Factor[x]; nivars <- length(intVars)
plotLetters<-c("A", "B", "C")
pdf(paste(plotPath,"InteractionEffects_ASdata150years.pdf", sep=""), height=4, width=12)
par(mfrow=c(1,2), mar=c(4,5.5,1.5,1), oma=c(0,0,0,10))
# plotSingleEffect(thisVar="TL", this_model=finalASModel)
# mtext(plotLetters[1], side=3, adj=0, line=0.1, font=2)
for(i in 1:nivars){
  thisIntVar <- intVars[i]; 
  plotIntEffect(thisIntVar = intVars[i], this_model=finalASModel, thisMax=0.3, updateMax = FALSE)
  mtext(plotLetters[(i+1)], side=3, adj=0, line=0.1, font=2)
}
predLegendValues <- pretty(seq(0,0.3, length.out=5))
legendValues <- predLegendValues # absoulte, and percentage
legendCols <- unlist(lapply(predLegendValues, getColor, log=FALSE))
# makeBlankPlot()
par(xpd=NA)
legend(legend=legendValues, col=legendCols, x="right", bty="n", pch=15, pt.cex=1.5, title="Absolute relative\ndifference (%)", inset=-0.5)

dev.off()
###############################################################################
##############################################################################################################################################################
 plot(this_model)
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



lab_cex <- 1
x <- grep(":", this_df$Factor); intVars <- this_df$Factor[x]; nivars <- length(intVars)
plotLetters<-c("A", "B", "C", "D")
pdf(paste(plotPath,"InteractionEffects_BPdata_allPointsPostBurnin.pdf", sep=""), height=3, width=12)
par(mfrow=c(1,3), mar=c(4,5.5,1.5,1), oma=c(0,0,0,7))
for(i in 1:nivars){
  thisIntVar <- intVars[i]; 
  thisMin <- 0
  plotIntEffect_noChaosRun(thisIntVar = thisIntVar, this_model=this_model, inputMax=0.1, updateMax = TRUE, thisData=testMelt[thisIndex,])
  mtext(plotLetters[i], side=3, adj=0, line=0.1, font=2)
}
predLegendValues <- pretty(seq(0,0.1, length.out=5))
legendValues <- predLegendValues # absoulte, and percentage
legendCols <- unlist(lapply(predLegendValues, getColor, log=FALSE, thisMax=0.1))
# makeBlankPlot()
par(xpd=NA)
legend(legend=legendValues*100, col=legendCols, x="right", bty="n", pch=15, pt.cex=1.5, title="Absolute relative\ndifference (%)", inset=-0.35, cex=1.1)
dev.off()
