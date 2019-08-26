# read in all the tracers, store in an array, and save this
# have 2 main sets of runs so far - SampleA and SampleKEY - both have fished and non-fished version
# Use DChaosBASE and DChaosFISH to compare
basePath <-paste(DIR$'Base',"AtlantisModels\\chaos\\",sep="")
thisChaosVersion <- "SampleA"; 
thisChaosVersion <- "FISHSampleA"; 
thisChaosVersion <- "SampleKEY"; 
thisChaosVersion <- "FISHSampleKEYPlusBP"; 
# thisChaosVersion <- "SampleKEYPlusBP"
thisChaosVersion <- "FISHSampleInformanceAll"
# thisChaosVersion <- "SampleInformanceAll"

nChaosRuns <- 35
# nChaosRuns <- 10
groupsDF <- read.csv(paste(basePath, "..\\CRAM_groups.csv",sep="")); ng <- dim(groupsDF)[1]
nts <- 151; nlayers <-6
storeNarray <- array(NA, dim=c(nChaosRuns, ng, nts))
store_nts <- rep(NA, nChaosRuns)

for(r in 1:nChaosRuns){
  thisNcFile <- paste(basePath,"outputChaos", thisChaosVersion,r,"\\output.nc", sep="")
  if(file.exists(thisNcFile)){
    ThisNC.nc <- nc_open(thisNcFile)
    thisVol<-ncvar_get(ThisNC.nc,"volume")
    this_nts <- dim(thisVol)[3]
    store_nts[r] <- this_nts
    # if(this_nts==nts){
      for(g in 1:ng){
        thisTracer <- paste(str_trim(groupsDF$Name[g],side="both"), "_N", sep="")
        thisData <- ncvar_get(ThisNC.nc, thisTracer)
        if(length(dim(thisData))==3){
          thisBiomass <- apply(thisData*thisVol, 3, sum) * mg_2_tonne * X_CN
        } else{
          thisBiomass <- apply(thisData*thisVol[nlayers,,],2,sum) *mg_2_tonne *X_CN
        }
        storeNarray[r,g,1:this_nts]<-thisBiomass
      }
    # }
  }
}

store_nts
# # save it
# save(list=c("storeNarray"), file=paste(basePath,"ChaosNtracers",thisChaosVersion,sep=""))

par(mfrow=c(2,2))
for(g in 1:ng){
  thisYmax<-max(storeNarray[,g,], na.rm=TRUE)
  plot(storeNarray[2,g,], type="l", ylim=c(0,thisYmax))
  for(r in 3:nChaosRuns){
    points(storeNarray[r,g,], type="l", col=myBlue_trans)
  }
  par(las=1)
  mtext(gsub("_"," ", groupsDF$Name[g]), side=3, adj=0)
}

## may as well grab the base runs and store them too
# baseRuns <- c("outputDChaosFISH", "outputDChaosBASE"); nBaseRuns <- length(baseRuns)
# baseArray <- array(NA, dim=c(nBaseRuns, ng, nts))
# for(r in 1:nBaseRuns){
#   ThisNC.nc <- nc_open(paste(basePath,baseRuns[r],"\\output.nc", sep=""))
#   thisVol<-ncvar_get(ThisNC.nc,"volume")
#   this_nts <- dim(thisVol)[3]
#   if(this_nts==nts){
#     for(g in 1:ng){
#       thisTracer <- paste(str_trim(groupsDF$Name[g],side="both"), "_N", sep="")
#       thisData <- ncvar_get(ThisNC.nc, thisTracer)
#       if(length(dim(thisData))==3){
#         thisBiomass <- apply(thisData*thisVol, 3, sum) * mg_2_tonne * X_CN
#       } else{
#         thisBiomass <- apply(thisData*thisVol[nlayers,,],2,sum) *mg_2_tonne *X_CN
#       }
#       baseArray[r,g,]<-thisBiomass
#     }
#   }
#   
# }

# save it
# save(list=c("baseArray"), file=paste(basePath,"ChaosNtracersBASEandBaseFISH_baseAarray",sep=""))



