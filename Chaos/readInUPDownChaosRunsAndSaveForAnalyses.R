# read in all the tracers, store in an array, and save this
# have 2 main sets of runs so far - SampleA and SampleKEY - both have fished and non-fished version
# Use DChaosBASE and DChaosFISH to compare
basePath <-paste(DIR$'Base',"AtlantisModels\\chaos\\",sep="")
thisChaosVersion <- "FISH"
thisChaosVersion <- ""
# upDown <- c("",paste(c("up","down"), sort(rep(c("01","02","005","05"),2)), sep=""))
upDown <- c(paste(c("up","down"), sort(rep(c("01","02","005","05"),2)), sep=""))

nChaosRuns <- length(upDown)
# nChaosRuns <- 10
groupsDF <- read.csv(paste(basePath, "..\\CRAM_groups.csv",sep="")); ng <- dim(groupsDF)[1]
nts <- 151; nlayers <-6
storeNarray <- array(NA, dim=c(nChaosRuns, ng, nts))
store_nts <- rep(NA, nChaosRuns)

for(r in 1:nChaosRuns){
  thisNcFile <- paste(basePath,"outputDChaos", thisChaosVersion,upDown[r],"\\output.nc", sep="")
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
  }
}

store_nts
# # save it
# save(list=c("storeNarray"), file=paste(basePath,"DChaosNtracers",thisChaosVersion,sep=""))
