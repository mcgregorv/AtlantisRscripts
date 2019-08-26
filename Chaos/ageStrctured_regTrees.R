source(paste(DIR$'General functions',"chaosFunctions.R", sep=""))
basePath <-paste(DIR$'Base',"AtlantisModels\\chaos\\",sep="")

#data were prepared in modelRelationshipAt50years_groupStructure_vs_stability.R
allData2fit<-read.csv( paste(DIR$'Tables',"Chaos_allData2fit_inclHalfWay.csv", sep=""))

groupsDF<- read.csv(paste(basePath,"..\\CRAM_groups.csv", sep="")); ng<-dim(groupsDF)[1]

as_codes <- as.character(groupsDF$Code[groupsDF$NumCohorts>1])
asData2fit <- allData2fit[allData2fit$Code %in% as_codes,]

hist(asData2fit$relDiffLast20years)

thisFactors <- c("Informance", "TL", "Keystone", "Response", "NumL1cons", "Lifespan", "propAdM", "propJuvM", "B0", "PropByTopPrey","ChaosRun", "Linf")
# note Linf  cuts out whales and such, might want to take it out
# thisFactors <- c("Informance", "TL", "Keystone", "Response", "NumL1cons", "Lifespan", "propAdM", "propJuvM", "B0", "PropByTopPrey","ChaosRun")

chaosVersions <- c("DChaosNtracers", "DChaosNtracersFISH", paste("ChaosNtracers",c("SampleInformanceAll", "FISHSampleInformanceAll", "SampleKEYPlusBP", "FISHSampleKEYPlusBP"), sep="")); nchaosVersions <- length(chaosVersions)
chaosRunDescr <- c("A: All up or down - no fishing", "B: All up or down - with fishing", "C: Weighted randon - no fishing", "D: Weighted random - with fishing", 
                   "E: Weighted Keystone - no fishing", "F: Weighted keystone - with fishing")


allInts <- getIntVars(thisFactors); all_factors <- c(thisFactors, allInts)

#turn categorical responses to factors
cat_factors <- c("Informance", "ChaosRun"); 
cat_factors <- c("ChaosRun")
ncat <- length(cat_factors)
for(c in 1:ncat){
  thisCat <- cat_factors[c]
  thisLevels <- sort(unique(asData2fit[,c(thisCat)]))
  asData2fit[,c(thisCat)] <- factor(asData2fit[,c(thisCat)], levels=thisLevels)
}

# take out any missing values
fitData <- asData2fit
for(f in 1:length(thisFactors)){
  index <- !is.na(fitData[,c(thisFactors[f])])
  xx <- rep(1,length(index))[!index]
  if(length(xx)>0){
    cat("drop ",sum(xx, na.rm=TRUE), "for ", thisFactors[f],"; proportion = ",signif(sum(xx, na.rm=TRUE)/length(index),2) ,"\n")
    fitData<-fitData[index,]
  }
}
fitData$logMaxDiff <- log(fitData$maxDiff)
fitData$logMedDiff <- log(fitData$medDiff)

thisRespVar<-"logMedDiff"

singleFormula <- as.formula(paste(thisRespVar, "~", paste(thisFactors, collapse=" + "), sep=""))
treeModel <- tree(formula=singleFormula, data=fitData)
plot(treeModel)
text(treeModel)
treeFactors <- as.character(summary(treeModel)$used)
testFactors <-c(treeFactors, getIntVars(treeFactors))

testGlm <- getBestLModel(thisVar=thisRespVar, thisFactors = testFactors, cutoff = 0.02, data2fit=fitData)

thisFormula <- as.formula(paste(thisRespVar, "~", paste(treeFactors, collapse=" + "), sep=""))
thisModel <- glm(thisFormula, data=fitData)
thisFitted <- predict.glm(thisModel, newdata = fitData)
plot(x=fitData[,c(thisRespVar)], y=thisFitted)


hist(fitData$logMaxDiff)



