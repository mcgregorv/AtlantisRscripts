source(paste(DIR$'General functions',"chaosFunctions.R", sep=""))
basePath <-paste(DIR$'Base',"AtlantisModels\\chaos\\",sep="")
plotPath <- paste(DIR$'Reports',"Chaos\\Figures\\", sep="") ## paper version
plotPath <- paste(DIR$'Figures',"Chaos\\", sep="")

#data were prepared in modelRelationshipAt50years_groupStructure_vs_stability.R
# allData2fit<-read.csv( paste(DIR$'Tables',"Chaos_allData2fit_inclHalfWay.csv", sep=""))
allData2fit <- read.csv(paste(DIR$'Tables',"Chaos_modelData_allPointsPostBurnin.csv", sep=""))

groupsDF<- read.csv(paste(basePath,"..\\CRAM_groups.csv", sep="")); ng<-dim(groupsDF)[1]

hist(allData2fit$RelBiomass)

thisFactors <- c("TL",  "NumL1cons", "Lifespan",  "B0", "PropByTopPrey")
# note Linf  cuts out whales and such, might want to take it out

allInts <- getIntVars(thisFactors); all_factors <- c(thisFactors, allInts)


# take out any missing values
fitData <- allData2fit
for(f in 1:length(thisFactors)){
  index <- !is.na(fitData[,c(thisFactors[f])])
  xx <- rep(1,length(index))[!index]
  if(length(xx)>0){
    cat("drop ",sum(xx, na.rm=TRUE), "for ", thisFactors[f],"; proportion = ",signif(sum(xx, na.rm=TRUE)/length(index),2) ,"\n")
    fitData<-fitData[index,]
  }
}

fitData$absRel <- abs(fitData$RelBiomass)
fitData$cubeRel <- (fitData$absRel)^(1/3)

fitData$root5thRel <- fitData$absRel^(1/5)
fitData$logRel <- log(fitData$RelBiomass +1)


# possibleRespVars <-c("logRatio", "logRel", "cubeRatio", "cubeRel", "root5thRel", "root6thRel", "logMed", "cubeMed", "root5thMed", "root6thMed")
# respVar_dfs <- NULL
# for(r in 1:(length(possibleRespVars))){
#   half_df <- getBestLModel(thisVar=possibleRespVars[r], thisFactors = all_factors, cutoff=0.02, data2fit=fitData)
#   respVar_dfs[[r]]<- half_df
# }
# names(respVar_dfs)<- possibleRespVars

thisRespVar<-"root5thRel"
this_df <- getBestLModel(thisVar = thisRespVar, thisFactors = all_factors, cutoff = 0.02, data2fit = fitData)


this_formula <- as.formula(paste(thisRespVar, "~", paste(this_df$Factor, collapse=" + "), sep=""))
this_model <- glm(this_formula, data=fitData)
testFit <- predict.glm(this_model, newdata = fitData)
par(mfrow=c(2,1))
hist_fit <- hist(testFit, plot=FALSE); hist_data <- hist(fitData[,c(thisRespVar)], plot=FALSE)
thisMin <- min(c(testFit, fitData[,c(thisRespVar)]), na.rm=TRUE)
thisMax <- max(c(testFit, fitData[,c(thisRespVar)]), na.rm=TRUE)
thisLine <- seq(thisMin, thisMax, length.out = 100)
par(mfrow=c(1,1))
plot(x=fitData[,c(thisRespVar)], y=testFit, pch=20, ylim=c(thisMin, thisMax), xlim=c(thisMin, thisMax))
points(x=thisLine, y=thisLine, type="l", col="red", lwd=2, lty=2)

# save the model out so can read it in else where
save(this_model, file=paste(DIR$'Data', "Chaos\\ALL_glmModel", sep=""))

# testResids <- fitData$root5thRel-testFit
# hist(fitData$root5thRel)
# hist(testFit)
# plot(x=fitData$root5thRel, y=testResids)

##################################### plot int effects
# thisFactors
# [1] "Informance"    "TL"            "Keystone"      "Response"      "NumL1cons"     "Lifespan"      "propAdM"       "propJuvM"      "B0"           
# [10] "PropByTopPrey" "ChaosRun"      "Linf"    

varDesc <- data.frame(matrix(NA, nrow=length(thisFactors), ncol=2))
colnames(varDesc)<- c("Variable", "Decr")
varDesc$Variable <- thisFactors
# varDesc$Decr <- c("Trophic level", "Keystone ranking", "Responsive ranking", "Level 1 connections", "Lifespan (years)", 
#                   "Additional adult mortality", "Additional juvenile mortality", expression(B[0]~"('000 tonnes)"),
#                   "Diet proportion by top prey", expression(L[infinity]~"(cm)"))
thisDescriptions <- c("Trophic level",  "Level 1 connections", "Lifespan (years)", "", "Diet proportion by top prey")

varDesc$Decr <- thisDescriptions


thisPredVars <- getPredVars(this_df); npvars <- length(thisPredVars)
thisPredVarLengths <- rep(NA, npvars)
thisPredVarValues <- NULL
for(v in 1:npvars){
  thisVar <- thisPredVars[v]
  tempVarValues<-getVarValues(thisVar, thisData=fitData)
  ## fix the odd things - don't bother with zero biomass, lifespan, Linf, 
  # and keystonennes and responsiveness should just have unique values
  if(thisVar %in% c("B0", "Lifespan", "Linf", "NumL1cons")){
    tempVarValues <- tempVarValues[tempVarValues!=0]
    # } else if(thisVar %in% c("Keystone", "Informance")){
  } else if(thisVar %in% c("Informance", "Keystone","Response")){
    thisUnique <- sort(unique(fitData[,c(thisVar)]))
    tempVarValues <- thisUnique
  }
  thisPredVarValues[[v]] <- tempVarValues
  thisPredVarLengths[v]<-length(thisPredVarValues[[v]])
}

lab_cex <- 1
x <- grep(":", this_df$Factor); intVars <- this_df$Factor[x]; nivars <- length(intVars)
plotLetters<-c("A", "B", "C", "D")
pdf(paste(plotPath,"InteractionEffects_ALLdata_allPointsPostBurnin.pdf", sep=""), height=6, width=12)
par(mfrow=c(2,3), mar=c(4,5.5,1.5,1), oma=c(0,0,0,0))
for(i in 1:nivars){
  thisIntVar <- intVars[i]; 
  thisMin <- 0
  plotIntEffect_noChaosRun(thisIntVar = thisIntVar, this_model=this_model, inputMax=0.1, updateMax = FALSE, thisData=fitData)
  mtext(plotLetters[i], side=3, adj=0, line=0.1, font=2)
}
plotSingleEffect(thisVar="TL", this_model=this_model)
mtext(plotLetters[4], side=3, adj=0, line=0.1, font=2)
predLegendValues <- pretty(seq(0,0.1, length.out=5))
legendValues <- predLegendValues # absoulte, and percentage
legendCols <- unlist(lapply(predLegendValues, getColor, log=FALSE, thisMax=0.1))
makeBlankPlot()
par(xpd=NA)
legend(legend=legendValues*100, col=legendCols, x="center", bty="n", pch=15, pt.cex=1.5, title="Absolute relative\ndifference (%)", inset=0, cex=1.5)
dev.off()

