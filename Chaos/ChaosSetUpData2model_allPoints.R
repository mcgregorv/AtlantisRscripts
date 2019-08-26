basePath <-paste(DIR$'Base',"AtlantisModels\\chaos\\",sep="")

allTheRankings <- read.csv(paste(DIR$'Tables', "allTheRankings.csv", sep=""))
## add informance - how well defined the group is combined with how well they performed in the historical model
rankDF<-read.csv(paste(DIR$'Tables',"interactionEffectsRANKINGs.csv",sep=""))

chaosVersions <- c("DChaosNtracers", "DChaosNtracersFISH", paste("ChaosNtracers",c("SampleInformanceAll", "FISHSampleInformanceAll", "SampleKEYPlusBP", "FISHSampleKEYPlusBP"), sep="")); nchaosVersions <- length(chaosVersions)
chaosRunDescr <- c("A: All up or down - no fishing", "B: All up or down - with fishing", "C: Weighted randon - no fishing", "D: Weighted random - with fishing", 
                   "E: Weighted Keystone - no fishing", "F: Weighted keystone - with fishing")
nchaosVersions<-length(chaosVersions)
nruns <- 50
groupsDF <- read.csv(paste(basePath, "..\\CRAM_groups.csv",sep="")); ng <- dim(groupsDF)[1]
nts <- 151; nlayers <-6
timeAxis <- seq(1900,2015,by=10); timeAxisAt <- timeAxis-1900
# thisGroupsDF <- groupsDF[groupsDF$Code != "DC",]; ng <- dim(thisGroupsDF)[1]
# keepGindex <- grep("DC", groupsDF$Code, invert = TRUE)

# load baseArray - has tracers from Chaos BASE and FISH runs (c("outputDChaosFISH", "outputDChaosBASE"))
load(paste(basePath,"ChaosNtracersBASEandBaseFISH_baseAarray",sep="")) # first in dim 1 is fishing

allTracersArray<- array(NA, dim=c(nchaosVersions, nruns, ng, nts))
allRelTracersArray<-0*allTracersArray

for(c in 1:nchaosVersions){
  thisChaosVersion <- chaosVersions[c]
  load(paste(basePath,thisChaosVersion,sep=""))
  this_nruns <- dim(storeNarray)[1]
  allTracersArray[c,1:this_nruns,,]<- storeNarray
  if(length(grep("FISH", thisChaosVersion))>0){
    thisBaseData <- baseArray[1,,]
  } else{
    thisBaseData <- baseArray[2,,]
  }
  for(r in 1:this_nruns){
    allRelTracersArray[c,r,,]<- (storeNarray[r,,]-thisBaseData)/thisBaseData
  }
}
timeIndex <- 35:nts
pdf(paste(plotPath,"AllChaosRelative2BaseBiomas.pdf", sep=""), height=12, width=10)
par(mfrow=c(11,5), mar=c(4,4,1,1))
for(g in 1:ng){
  if(groupsDF$Code[g]!="DC"){
    thisName <- gsub("_", " ", groupsDF$LongName[g])
    thisData <- allRelTracersArray[,,g,timeIndex]; thisMax <- max(thisData, na.rm=TRUE)
    plot(thisData[1,1,], type="n", ylim=c(0, thisMax), ylab="Biomass (tonnes)", xlab="", xaxt="n")
    for(c in 1:nchaosVersions){
      thisChaosDescrip <- chaosRunDescr[c]
      if(length(grep("no fish", thisChaosDescrip))>0){
        thisCol<-myBlue_trans
      } else{
        thisCol=myOrange_trans
      }
      for(r in 1:nruns){
        points(thisData[c,r,], type="l", lwd=2, col=thisCol)
      }
    }
    axis(at=timeAxisAt, labels=timeAxis, side=1)
    mtext(thisName, side=3, adj=0)
  }
}
dev.off()

modelData <- data.frame(matrix(NA, nrow=0, ncol=dim(allTheRankings)[2]+2))
colnames(modelData)<- c(colnames(allTheRankings), "Informance", "RelBiomass")

for(g in 1:ng){
  thisCode <- as.character(groupsDF$Code[g]); thisNumCohorts <- groupsDF$NumCohorts[g]
  if(thisCode != "DC"){
    thisBiolData <- allTheRankings[allTheRankings$Code==thisCode,]
    
    if(thisNumCohorts==1){
      thisBiolData$Keystone <- 0
    }
    thisInf <- rankDF$informPerformRank[rankDF$Code==thisCode]
    if(length(thisInf)==0){thisInf <- 16}
    
    
    thisData <- allRelTracersArray[,,g,timeIndex]
    xx <- as.vector(thisData); thisRels <- xx[!is.na(xx)]; nr <- length(thisRels)
    this_df <- data.frame(matrix(NA, nrow=nr, ncol=dim(allTheRankings)[2]+2))
    colnames(this_df)<- c(colnames(allTheRankings), "Informance", "RelBiomass")
    this_df[,colnames(allTheRankings)]<- thisBiolData
    this_df$Informance <- rep(thisInf, nr)
    this_df$RelBiomass <- thisRels
    modelData <- rbind(modelData, this_df)
  }
}

# write.csv(modelData, paste(DIR$'Tables',"Chaos_modelData_allPointsPostBurnin.csv", sep=""), row.names = FALSE)



