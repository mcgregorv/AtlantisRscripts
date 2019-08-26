nruns <- 10; sampleVersion<-"A"
runSetupFile <- paste(basePath, "Chaos\\Setup\\", sep="")
thisSampleFile <- paste(runSetupFile,"chaosSampleSetup",sampleVersion,"n",nruns,".csv", sep="")

storeSampledScalars<-read.csv(thisSampleFile)

## plot the resulting spreads, colored by informance
colRamp<-colorRampPalette(colors=c("red",myOrange,myAqua,myBlue))(4)

uniqueGradings<-sort(unique(rankTable$informPerformRank))
getColor<-function(x){
  # y<-(x-17 + 0.5)*2
  thisCol<-colRamp[match(x,uniqueGradings)]
  return(thisCol)
}
gradingsByVar <- rep(NA, nvars)
for(i in 1:nvars){
  thisVar <- vars[i]
  thisCohort <- get_first_number(thisVar); thisVar <-gsub(thisCohort, "", thisVar)
  thisName <- gsub("_N|_Nums", "", thisVar); thisName <- gsub("_", " ", thisName)
  thisCode <- rankTable$Code[grep(thisName, rankTable$Group)]
  thisRank<-rankTable$informPerformRank[grep(thisName, rankTable$Group)]
  if(!(thisRank %in% uniqueGradings)){ thisRank <- NA}
  gradingsByVar[i]<-thisRank
}
colByPer<-unlist(lapply(gradingsByVar, getColor))
colByPer[is.na(colByPer)]<-myGrey

addedSampledScalars <- 1+storeSampledScalars

thisMax <- max(addedSampledScalars, na.rm=TRUE); thisMin <- min(addedSampledScalars, na.rm=TRUE)
pdf(paste(plotPath, "scalesUsed",sampleVersion,"n",nruns,".pdf", sep=""))
plot(as.double(addedSampledScalars[1,]), type="n", ylim=c(thisMin, thisMax), xlab="Variable", ylab="Scalar")
for(r in 1:nruns){
  points(as.double(addedSampledScalars[r,]), pch=20, col=colByPer)
}
dev.off()

