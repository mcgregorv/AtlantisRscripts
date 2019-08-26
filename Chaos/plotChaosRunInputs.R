# read in the initial conditions used for the chaos runs and plot
basePath <-paste(DIR$'Base',"AtlantisModels\\chaos\\",sep="")

thisCex=1.65

plotPath <- paste(DIR$'Reports',"Chaos\\Figures\\", sep="") ## paper version
plotPath <- paste(DIR$'Figures',"Chaos\\", sep="")

rankTable <- read.csv(paste(DIR$'Table',"allTheRankings.csv", sep=""))
groupsDF <- read.csv(paste(basePath,"..\\CRAM_groups.csv", sep="")); ng <- dim(groupsDF)[1]
## vars were created in setUpSampleScalers_keystoneSpecies.R - which was used for the random and keystone versions
vars <- read.csv(paste(DIR$'Table',"VARSforInitialConditions.csv", sep=""))[,1]; nvars <- length(vars)

nruns <- 10; sampleKeyVersion<-"KEY" # this was the original 10-run version, which had only keystone (which were all age-structured) groups varied
sampleKeyVersion <-"KEYPlusBP"; nruns=35 ## has top 10 keystone and all biomass-pool gorups varied
runSetupFile <- paste(basePath, "Setup\\", sep="")
thisSampleFile <- paste(runSetupFile,"chaosSampleSetup",sampleKeyVersion,"n",nruns,".csv", sep="")
keystoneSamples <- read.csv(file=thisSampleFile)
informanceRankings <- read.csv(paste(DIR$'Tables',"interactionEffectsRANKINGs.csv", sep=""))
sampleVersion<-"A" # 10 runs, varied by informance
sampleVersion <-"InformanceAll" # same as A, but 35 runs

thisSampleFile <- paste(runSetupFile,"chaosSampleSetup",sampleVersion,"n",nruns,".csv", sep="")
randomSamples <- read.csv(file=thisSampleFile)

ksByVar <- rep(NA, nvars); informanceByVar <- rep(NA, nvars)
for(i in 1:nvars){
  thisVar <- vars[i]
  thisCohort <- get_first_number(thisVar); thisVar <-gsub(thisCohort, "", thisVar)
  thisName <- gsub("_N|_Nums", "", thisVar); 
  thisCode <- groupsDF$Code[grep(thisName, groupsDF$Name)]
  thisNumCohorts <- groupsDF$NumCohorts[groupsDF$Code==thisCode]
  thisKS<-rankTable$Keystone[rankTable$Code==as.character(thisCode)]
  ksByVar[i]<-thisKS
  if(is.na(thisCohort)){ksByVar[i] <- 0}
  informanceByVar[i]<-informanceRankings$informPerformRank[informanceRankings$Code==as.character(thisCode)]
  if(is.na(thisCohort)){informanceByVar[i] <- 16}
}

getColor<-function(x){
  thisCol<-myGrey
  if(!is.na(x)>0){
    if(x<=10){thisCol<-myOrange}
  }
  return(thisCol)
}
colByPer<-unlist(lapply(ksByVar, getColor))
colByPer[is.na(colByPer)]<-myGrey

addedSampledScalars <- 1+keystoneSamples

thisMax <- max(addedSampledScalars, na.rm=TRUE); thisMin <- min(addedSampledScalars, na.rm=TRUE)
# pdf(paste(plotPath, ""))
par(xaxs="r")
plot(as.double(addedSampledScalars[1,]), type="n", ylim=c(thisMin, thisMax), xlim=c(1,nvars), xlab="Variable", ylab="Scalar")
abline(h=seq(0.6,1.6,by=0.1), col=myGrey_trans)
abline(v=seq(0,361,by=25), col=myGrey_trans)
for(r in 1:nruns){
  points(as.double(addedSampledScalars[r,]), pch=20, col=colByPer)
}

randomAddedSampledScalars <- 1+randomSamples

thisMax <- max(randomAddedSampledScalars, na.rm=TRUE); thisMin <- min(randomAddedSampledScalars, na.rm=TRUE)
# pdf(paste(plotPath, ""))
plot(as.double(randomAddedSampledScalars[1,]), type="n", ylim=c(thisMin, thisMax), xlab="Variable", ylab="Scalar")
abline(h=seq(0.6,1.6,by=0.1), col=myGrey_trans)
abline(v=seq(0,361,by=25), col=myGrey_trans)
for(r in 1:nruns){
  points(as.double(randomAddedSampledScalars[r,]), pch=20, col=colByPer)
}
# index tracers that are biomass pools
biomPoolIndex <- grep("_Nums", vars, invert = TRUE); nbp <- length(biomPoolIndex)
for(b in 1:nbp){
  this_b <- biomPoolIndex[b]
  points(x=rep(this_b,nruns), y=as.double(randomAddedSampledScalars[,this_b]), pch=20, col=myGreen)
}

## how were the initial scalars defined..?
# for the keystone runs, those with high keystoneness were sampled from normal distributions with mean 0, sd 0.25 (95%CIs ~ +/- 0.5)
## for the random runs, sd was set based on informance ranking
nageSdev<-0.05
## there are 4 informance rankings, so assign 4 standard deviations
resultingScalars <- rev(c(0.05, 0.1, 0.2, 0.5))
ageDevs <- (resultingScalars/1.96)
par(mfrow=c(2,2))
for(a in 1:length(ageDevs)){
  thisX=seq(-1,1,length.out=100); thisY<-dnorm(thisX, sd=ageDevs[a])
  plot(x=thisX, y=thisY, type="l")
  abline(v=c(-1,1)*resultingScalars[a], col=myOrange, lwd=2)
}

thisMin<-0.4; thisMax <- 1.6
# pdf(paste(plotPath, "KeystoneScalarSampleDistributions.pdf", sep=""), height=3, width=5)
pdf(paste(plotPath, sampleKeyVersion,"_ScalarSampleDistributions.pdf", sep=""), height=3, width=5)
par(fig=c(0,0.8,0,1), mar=c(4,4,1,0))
plot(as.double(addedSampledScalars[1,]), type="n", ylim=c(thisMin, thisMax), xlim=c(1,nvars), xlab="Variable", ylab="Scalar")
abline(h=seq(thisMin, thisMax,by=0.1), col=myGrey_trans)
abline(v=seq(0,361,by=25), col=myGrey_trans)
for(r in 1:nruns){
  points(as.double(addedSampledScalars[r,]), pch=20, cex=0.5, col=colByPer)
  # points(as.double(addedSampledScalars[r,]))
}
par(fig=c(0.8,1,0,1), mar=c(4,0,1,0.1))
par(new=TRUE)
thisX=seq(-1,1,length.out=100); thisY<-dnorm(thisX, sd=0.25)
plot(y=thisX+1, x=thisY, type="l", col=myOrange, lwd=2, xlab="Probability", ylim=c(thisMin, thisMax), yaxt="n", ylab="")
# abline(h=c(-1,1)*(0.25*1.96), col=myGrey, lwd=2, lty=3)
dev.off()


## the random version
colByDiv <- c(myOrange, myRed, myBlue, myGreen)
infColIndex <- c(16,18,19,20)
getColorInformance<-function(x){
  thisCol<-myBlue
  if(!is.na(x)>0){
    thisCol<-colByDiv[match(x,infColIndex)]
  }
  return(thisCol)
}
colByInf<-unlist(lapply(informanceByVar, getColorInformance))
colByInf[is.na(colByInf)]<-myBlue ## this matches sd = 0.05, which is what was used for non-age structured (biomass pool) groups

pdf(paste(plotPath, sampleVersion, "_RandomScalarSampleDistributions.pdf", sep=""), height=3, width=5)
par(fig=c(0,0.8,0,1), mar=c(4,4,1,0))
plot(as.double(randomAddedSampledScalars[1,]), type="n", ylim=c(thisMin, thisMax), xlim=c(1,nvars), xlab="Variable", ylab="Scalar")
abline(h=seq(thisMin, thisMax,by=0.1), col=myGrey_trans)
abline(v=seq(0,361,by=25), col=myGrey_trans)
for(r in 1:nruns){
  points(as.double(randomAddedSampledScalars[r,]), pch=20, cex=0.5, col=colByInf)
}
par(fig=c(0.8,1,0,1), mar=c(4,0,1,0.1))
par(new=TRUE)
thisX=seq(-1,1,length.out=1000); thisY<-dnorm(thisX, sd=0.05)
plot(y=thisX+1, x=thisY, type="l", col=myGrey, lwd=2, xlab="Probability", ylim=c(thisMin, thisMax), yaxt="n", ylab="", xlim=c(0,16))
for(a in 1:length(ageDevs)){
  thisX=seq(-1,1,length.out=1000); thisY<-dnorm(thisX, sd=ageDevs[a])
  points(y=thisX+1, x=thisY, type="l", col=colByDiv[a], lwd=2)
  
}
dev.off()


