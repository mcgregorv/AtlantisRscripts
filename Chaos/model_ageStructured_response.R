source(paste(DIR$'General functions',"chaosFunctions.R", sep=""))
basePath <-paste(DIR$'Base',"AtlantisModels\\chaos\\",sep="")
plotPath <- paste(DIR$'Reports',"Chaos\\Figures\\", sep="") ## paper version
plotPath <- paste(DIR$'Figures',"Chaos\\", sep="")

#data were prepared in modelRelationshipAt50years_groupStructure_vs_stability.R
allData2fit<-read.csv( paste(DIR$'Tables',"Chaos_allData2fit_inclHalfWay.csv", sep=""))
# allData2fit <- read.csv(paste(DIR$'Tables',"Chaos_modelData_allPointsPostBurnin.csv", sep=""))

groupsDF<- read.csv(paste(basePath,"..\\CRAM_groups.csv", sep="")); ng<-dim(groupsDF)[1]

as_codes <- as.character(groupsDF$Code[groupsDF$NumCohorts>1])
asData2fit <- allData2fit[allData2fit$Code %in% as_codes,]

hist(asData2fit$relDiffLast20years)

thisFactors <- c("Informance", "TL", "Keystone", "Response", "NumL1cons", "Lifespan", "propAdM", "propJuvM", "B0", "PropByTopPrey","ChaosRun", "Linf")
# note Linf  cuts out whales and such, might want to take it out
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

fitData$logRatio <- log(fitData$ratioDiffLast20years+1)
fitData$logRel <- log(fitData$relDiffLast20years + 1)
fitData$cubeRatio <- fitData$ratioDiffLast20years^(1/3)
fitData$cubeRel <- (fitData$relDiffLast20years)^(1/3)
fitData$root5thRel <- fitData$relDiffLast20years^(1/5)
fitData$root6thRel <- fitData$relDiffLast20years^(1/6)
fitData$root7thRel <- fitData$relDiffLast20years^(1/7)
fitData$logMed <- log(fitData$medDiff+1)
fitData$cubeMed <- fitData$medDiff^(1/3)
fitData$root5thMed <- fitData$medDiff^(1/5)
fitData$root6thMed <- fitData$medDiff^(1/6)
# 

thisRespVar<-"root5thRel"
this_dfs<-NULL
for(c in 1:nchaosVersions){
  thisChaosRun <- chaosVersions[c]; chaosRunLetter <- substr(thisChaosRun,start=1, stop=1)
  chaosRunIndex <-  fitData$ChaosRun==chaosRunLetter 
  half_df <- getBestLModel(thisVar=thisRespVar, thisFactors = all_factors[grep("Chaos", all_factors, ignore.case = TRUE, invert = TRUE)], cutoff=0.02, data2fit=fitData[chaosRunIndex,])
  this_dfs[[c]]<- half_df
}

possibleRespVars <-c("logRatio", "logRel", "cubeRatio", "cubeRel", "root5thRel", "root6thRel", "logMed", "cubeMed", "root5thMed", "root6thMed")
respVar_dfs <- NULL
for(r in 1:(length(possibleRespVars))){
  half_df <- getBestLModel(thisVar=possibleRespVars[r], thisFactors = all_factors, cutoff=0.02, data2fit=fitData)
  respVar_dfs[[r]]<- half_df
}
names(respVar_dfs)<- possibleRespVars

thisRespVar<-"root5thRel"
this_df <- getBestLModel(thisVar = thisRespVar, thisFactors = all_factors, cutoff = 0.02, data2fit = fitData)


this_formula <- as.formula(paste(thisRespVar, "~", paste(this_df$Factor, collapse=" + "), sep=""))
this_model <- glm(this_formula, data=fitData)
testFit <- predict.glm(this_model, newdata = fitData)
par(mfrow=c(2,1))
hist(testFit); hist(fitData[,c(thisRespVar)])
thisMin <- min(c(testFit, fitData[,c(thisRespVar)]), na.rm=TRUE)
thisMax <- max(c(testFit, fitData[,c(thisRespVar)]), na.rm=TRUE)
thisLine <- seq(thisMin, thisMax, length.out = 100)
par(mfrow=c(1,1))
plot(x=fitData[,c(thisRespVar)], y=testFit, pch=20, ylim=c(thisMin, thisMax), xlim=c(thisMin, thisMax))
points(x=thisLine, y=thisLine, type="l", col="red", lwd=2, lty=2)

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
thisDescriptions <- c("Informance", "Trophic level",  "Keystone ranking", "Responsive ranking", "Level 1 connections", "Lifespan (years)", 
                      "Additional adult mortality", "Additional juvenile mortality", "", "Diet proportion by top prey","Simulation set", "Linf")

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
thisChaosRuns<-c("A", "B", "C", "D", "E", "F")
lab_cex <- 1
x <- grep(":", this_df$Factor); intVars <- this_df$Factor[x]; nivars <- length(intVars)
plotLetters<-c("A", "B", "C")
pdf(paste(plotPath,"InteractionEffects_ASdata.pdf", sep=""), height=6, width=12)
par(mfrow=c(2,2), mar=c(4,5.5,1.5,1), oma=c(0,0,0,7))
for(i in 1:nivars){
  thisIntVar <- intVars[i]; 
  thisMin <- 0
  plotIntEffect(thisIntVar = thisIntVar, this_model=this_model, inputMax=0.3, updateMax = TRUE, thisData=fitData)
  # mtext(plotLetters[i], side=3, adj=0, line=0.1, font=2)
}

predLegendValues <- pretty(seq(0,0.3, length.out=5))
legendValues <- predLegendValues # absoulte, and percentage
legendCols <- unlist(lapply(predLegendValues, getColor, log=FALSE))
# makeBlankPlot()
par(xpd=NA)
legend(legend=legendValues, col=legendCols, x="right", bty="n", pch=15, pt.cex=1.5, title="Absolute relative\ndifference (%)", inset=-0.35)

dev.off()




# testing testing - not used 
# which groups see to be more or less changed by the end of the runs..?
# nag <- length(as_codes)
# plot(x=seq(1, nag), y=rep(0, nag), type="n")
# for(g in 1:nag){
#   thisCode <- as_codes[g]
#   thisData <- fitData[fitData$Code==thisCode,c("relDiffLast20years")]
#   boxplot(thisData, add=TRUE, at=g, col=myGreen_trans, border=myGreen, outline=FALSE)
# }
# # WHICH groups are within +/1 1?
# xx <- tapply(abs(asData2fit$relDiffLast20years), asData2fit$Code,FUN=function(x){getPercentile(x, k=0.50)[2]})
# maxAbsByGroup <- sort(xx[!is.na(xx)], decreasing = TRUE)
# 
# groupsByMax <- names(maxAbsByGroup)
# par(mfrow=c(2,1))
# thisMax<-1.2
# plot(x=seq(1, nag), y=rep(0, nag), type="n", ylim=c(0, thisMax), xaxt="n", xlab="")
# par(las=2)
# axis(at=1:nag, labels=groupsByMax, side=1)
# for(g in 1:nag){
#   thisCode <- groupsByMax[g]
#   thisData <- abs(fitData[fitData$Code==thisCode,c("relDiffLast20years")])
#   boxplot(thisData, add=TRUE, at=g, col=myBlue_trans, border=myBlue, outline=FALSE)
# }
# polygon(x=c(0.5,4.4,4.4,0.5), y=c(-0.02,-0.02,thisMax,thisMax), col=NA, border=myOrange)
# polygon(x=c(4.5,11.4,11.4,4.5), y=c(-0.02,-0.02,thisMax,thisMax), col=NA, border=myGreen)
# polygon(x=c(11.5,19.4,19.4,11.5), y=c(-0.02,-0.02,thisMax,thisMax), col=NA, border=myAqua)
# polygon(x=c(19.5,29.4,29.4,19.5), y=c(-0.02,-0.02,thisMax,thisMax), col=NA, border=myPurple)
# polygon(x=c(29.5,37.4,37.4,29.5), y=c(-0.02,-0.02,thisMax,thisMax), col=NA, border="red")
# 
# allTheRankings <- read.csv(paste(DIR$'Tables', "allTheRankings.csv", sep="")); nrankings <- dim(allTheRankings)[2]
# # add in informance
# rankDF<-read.csv(paste(DIR$'Tables',"interactionEffectsRANKINGs.csv",sep=""))
# 
# allTheRankings$Informance <- rankDF$informPerformRank[match(allTheRankings$Code,rankDF$Code)]
# 
# 
# summaryFit <- allTheRankings[allTheRankings$Code %in% as_codes,]
# summaryFit$ResponseCat <-"E"
# summaryFit$ResponseCat[summaryFit$Code %in% c("SPE", "GSH", "DPI", "IVS")]<-"A"
# summaryFit$ResponseCat[summaryFit$Code %in% c("EIS", "IVH", "BEE", "ELI", "CBO", "ASQ", "MJE")]<-"B"
# summaryFit$ResponseCat[summaryFit$Code %in% c( "SPD", "PFS", "PFM", "HOK", "BIS", "MAC", "LDO", "ELP")]<-"C"
# summaryFit$ResponseCat[summaryFit$Code %in% c("CET", "SB",  "LIN", "PIN", "ETB", "BAL", "ORH", "PFL", "HAK", "CEP", "JAV")]<-"D"
# 
# thisFactors
# responseLevels <- c("A", "B", "C", "D", "E"); nres <- length(responseLevels)
# par(las=1)
# thisExpl <- "Keystone"
# for(f in 1:length(thisFactors)){
#   thisExpl <- thisFactors[f]
#   if(!(thisExpl %in% cat_factors)){
#     thisYmax <- max(summaryFit[,c(thisExpl)], na.rm=TRUE); thisYmin <- min(summaryFit[,c(thisExpl)], na.rm=TRUE)
#     plot(x=1:nres, y=rep(0,nres), ylim=c(thisYmin, thisYmax), type="n", xaxt="n", xlab="")
#     axis(at=1:nres, labels=responseLevels, side=1)
#     for(r in 1:nres){
#       thisRes <- responseLevels[r]
#       thisData <- summaryFit[summaryFit$ResponseCat==thisRes, c(thisExpl)]
#       boxplot(thisData, add=TRUE, at=r, col=myRed_trans, border=myRed)
#     }
#     mtext(thisExpl, side=3, adj=0)
#   }
# }
#  ## at 50 years ..?
# xx <- tapply(abs(asData2fit$half_relDiff), asData2fit$Code,FUN=function(x){getPercentile(x, k=0.50)[2]})
# maxAbsByGroup <- sort(xx[!is.na(xx)], decreasing = TRUE)
# 
# plot(maxAbsByGroup, type="h", lwd=5, col=myBlue)
# 
# groupsByMax <- names(maxAbsByGroup)
# par(mfrow=c(2,1))
# thisMax<-0.6
# plot(x=seq(1, nag), y=rep(0, nag), type="n", ylim=c(0, thisMax), xaxt="n", xlab="")
# par(las=2)
# axis(at=1:nag, labels=groupsByMax, side=1)
# for(g in 1:nag){
#   thisCode <- groupsByMax[g]
#   thisData <- abs(asData2fit[asData2fit$Code==thisCode,c("half_relDiff")])
#   boxplot(thisData, add=TRUE, at=g, col=myBlue_trans, border=myBlue, outline=FALSE)
# }
# groups <- NULL
# groups[[1]]<-1:4; groups[[2]]<-5:9; groups[[3]]<-10:15; groups[[4]]<-16:24; groups[[5]]<-25:nag
# colByPoly <- colorRampPalette(colors=c(myGold, myGreen, myAqua, myPurple, "red"))(length(groups))
# for(p in 1:length(groups)){
#   x0<-min(groups[[p]])-0.5; x1 <- max(groups[[p]])+0.4
#   this_x <- c(x0, x1, x1, x0); this_y <- c(-0.02,-0.02,thisMax,thisMax)
#   polygon(x=this_x, y=this_y, col=NA, border=colByPoly[p])
# }
# 
# 
# summaryHalfFit <- allTheRankings[allTheRankings$Code %in% as_codes,]
# summaryHalfFit$ResponseCat <-"E"
# summaryHalfFit$ResponseCat[match(groupsByMax[groups[[1]]], summaryHalfFit$Code)]<-"A"
# summaryHalfFit$ResponseCat[match(groupsByMax[groups[[2]]], summaryHalfFit$Code)]<-"B"
# summaryHalfFit$ResponseCat[match(groupsByMax[groups[[3]]], summaryHalfFit$Code)]<-"C"
# summaryHalfFit$ResponseCat[match(groupsByMax[groups[[4]]], summaryHalfFit$Code)]<-"D"
# 
# thisFactors
# responseLevels <- c("A", "B", "C", "D", "E"); nres <- length(responseLevels)
# par(las=1)
# thisExpl <- "Keystone"
# for(f in 1:length(thisFactors)){
#   thisExpl <- thisFactors[f]
#   if(!(thisExpl %in% cat_factors)){
#     thisYmax <- max(summaryFit[,c(thisExpl)], na.rm=TRUE); thisYmin <- min(summaryFit[,c(thisExpl)], na.rm=TRUE)
#     plot(x=1:nres, y=rep(0,nres), ylim=c(thisYmin, thisYmax), type="n", xaxt="n", xlab="")
#     axis(at=1:nres, labels=responseLevels, side=1)
#     for(r in 1:nres){
#       thisRes <- responseLevels[r]
#       thisData <- summaryFit[summaryFit$ResponseCat==thisRes, c(thisExpl)]
#       boxplot(thisData, add=TRUE, at=r, col=myRed_trans, border=myRed)
#     }
#     mtext(thisExpl, side=3, adj=0)
#   }
# }
# 
# respvar <- "ResponseCat"
# getBestLModel(thisVar=respvar, thisFactors = c("TL", "Informance"), cutoff=0.01, data2fit = summaryHalfFit)
# 
# varDesc <- data.frame(matrix(NA, nrow=length(thisFactors), ncol=2))
# colnames(varDesc)<- c("Variable", "Decr")
# varDesc$Variable <- thisFactors
# thisDescriptions <- c("Trophic level", "Keystone ranking", "Responsive ranking", "Level 1 connections", "Lifespan (years)",
#                   "Additional adult mortality", "Additional juvenile mortality", expression(B[0]~"('000 tonnes)"),
#                   "Diet proportion by top prey", expression(L[infinity]~"(cm)"))
# # thisDescriptions <- c("Informance", "Trophic level",  "Keystone ranking", "Responsive ranking", "Level 1 connections", "Lifespan (years)", 
# #                       "Additional adult mortality", "Additional juvenile mortality", "", "Diet proportion by top prey","Simulation set")
# # 
# varDesc$Decr <- thisDescriptions
# 

