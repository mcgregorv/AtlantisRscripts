basePath <-paste(DIR$'Base',"AtlantisModels\\chaos\\",sep="")

chaosVersions <- c("DChaosNtracers", "DChaosNtracersFISH", paste("ChaosNtracers",c("SampleInformanceAll", "FISHSampleInformanceAll", "SampleKEYPlusBP", "FISHSampleKEYPlusBP"), sep="")); nchaosVersions <- length(chaosVersions)
chaosRunDescr <- c("A: All up or down - no fishing", "B: All up or down - with fishing", "C: Weighted randon - no fishing", "D: Weighted random - with fishing", 
                   "E: Weighted Keystone - no fishing", "F: Weighted keystone - with fishing")
nchaosVersions<-length(chaosVersions)
nruns <- 50
groupsDF <- read.csv(paste(basePath, "..\\CRAM_groups.csv",sep="")); ng <- dim(groupsDF)[1]
nts <- 151; nlayers <-6
timeAxis <- seq(1870,2015,by=10); timeAxisAt <- timeAxis-1865
# thisGroupsDF <- groupsDF[groupsDF$Code != "DC",]; ng <- dim(thisGroupsDF)[1]
# keepGindex <- grep("DC", groupsDF$Code, invert = TRUE)

allTracersArray<- array(NA, dim=c(nchaosVersions, nruns, ng, nts))

for(c in 1:nchaosVersions){
  thisChaosVersion <- chaosVersions[c]
  # bring in storeNarray
  # load(paste(basePath,"ChaosNtracers",thisChaosVersion,sep=""))
  load(paste(basePath,thisChaosVersion,sep=""))
  this_nruns <- dim(storeNarray)[1]
  allTracersArray[c,1:this_nruns,,]<- storeNarray
}
pdf(paste(plotPath,"AllChaosBiomasTracers.pdf", sep=""), height=12, width=10)
par(mfrow=c(11,5), mar=c(4,4,1,1))
for(g in 1:ng){
  if(groupsDF$Code[g]!="DC"){
    thisName <- gsub("_", " ", groupsDF$LongName[g])
    thisData <- allTracersArray[,,g,]; thisMax <- max(thisData, na.rm=TRUE)
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


# load baseArray - has tracers from Chaos BASE and FISH runs (c("outputDChaosFISH", "outputDChaosBASE"))
load(paste(basePath,"ChaosNtracersBASEandBaseFISH_baseAarray",sep=""))


