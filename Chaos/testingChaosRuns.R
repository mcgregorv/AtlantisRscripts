#plot all tracers for a given box and layer

this_run<-"Chaos"
this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")

this_out <- paste("D",c("ChaosUp05", "ChaosUp02", "ChaosUp01", "ChaosUp005","ChaosBASE", "ChaosDown005", "ChaosDown01", "ChaosDown02", "ChaosDown05"), sep="");
plotDescrip<-"DChaosUpDownBASE"; baseRunIndex<-5

this_out <- c("DChaosBASE" ,paste("ChaosSampleA", 1:10, sep=""));
plotDescrip<-"ChaosSampleA"; baseRunIndex<-1


thisRunTemplate <- "ChaosFISH"
thisRunTemplate <- "Chaos|short"
thisRunTemplate <- "DChaos"
thisRunTemplate <- "ChaosSampleA"

thisDesc<-plotDescrip

nlayers<-6

groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

plotPath<-paste(this_path,"..\\Figures\\", plotDescrip,sep="")

nruns<-length(this_out)
runCols<-colorRampPalette(colors=c("midnightblue",myBlue,myAqua,myGold,  myOrange, "red"))(nruns)

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1865
fishingStartYear<-1865
modelStartYear<-1865

nc_list<-NULL; nts_list<-NULL; min_nts<-1e+12
for(r in 1:nruns){
  outPath<-paste(this_path,"output",this_out[r],"\\",sep="")
  nc_list[[r]]<-nc_open(paste(outPath,"output.nc",sep=""))
  thisVol<-ncvar_get(nc_list[[r]],"volume")
  thisDz<-ncvar_get(nc_list[[r]],"dz")
  nts_list[[r]]<-dim(thisVol)[3] #number of timesteps
  if(nts_list[[r]]<min_nts){min_nts<-nts_list[[r]]}
}

nts_list
max_nts<-max(nts_list, na.rm=TRUE)
burnin<-rep(0,nruns)
timeList<-NULL; timeMin <- 30000; timeMax <- 0
for(r in 1:nruns){
  this_nts<-nts_list[[r]]; this_burnin <- burnin[r]
  if(is.na(this_burnin)){this_burnin <-0}
  thisYear0<-1900 - this_burnin + 1
  this_time <- thisYear0 : (thisYear0 + this_nts -1)
  timeList[[r]]<-this_time
  if(max(this_time) > timeMax){timeMax<-max(this_time)}
  if(min(this_time) < timeMin){timeMin <- min(this_time)}
}
xLabsTemp<-seq(0,(max_nts*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0+burnin[1]

#get all tracer names
allTracers<-sort(names(nc_list[[r]]$var))
temp<-allTracers[grep("_N",allTracers)]; tracers2plot<-temp[grep("Nums",temp,invert = TRUE)]; 
tracers2plot<-tracers2plot[grep("Carrion_N", tracers2plot, invert=TRUE)] ## take out carrion as not using and temperature as forced
# tracers2plot<-c(tracers2plot,"Oxygen","Si", "NO3")
ntracers<-length(tracers2plot)

## grab just the age-structured so can look at those on their own
ageNames <- str_trim(groupsDF$Name[groupsDF$NumCohorts>1], side="both"); nageGroups <- length(ageNames)
ageTracers <- paste(ageNames, "_N", sep="")
ageIndex <- tracers2plot %in% ageTracers


## populate tracers array - i think at this stage, just summarise over space (so biomass as time series)
tracersArray <- array(NA, dim=c(nruns, ntracers, max_nts))
storeWeightsArray<-array(NA, dim=c(nruns, ntracers, 10, max_nts)); storeNumbersArray<-0*storeWeightsArray
for(r in 1:nruns){
  cat(r,"--")
  ThisNC.nc<-nc_list[[r]]
  thisVol<-ncvar_get(ThisNC.nc,"volume"); thisData<-ncvar_get(ThisNC.nc,"Pelagic_fish_sml_N")
  this_nts<-nts_list[[r]]
  if(!is.na(this_nts)){
     for(t in 1:ntracers){
      thisTracer<-tracers2plot[t]; thisName<-gsub("_N","",thisTracer); xx<-grep(thisName,groupsDF$Name)
      thisCode<-groupsDF$Code[xx]; thisNumCohorts<-groupsDF$NumCohorts[xx]
      #do biomass for all
      thisData<-ncvar_get(ThisNC.nc,thisTracer)
      if(length(dim(thisData))==3){
        xx<-apply(thisData*thisVol,3,sum)*mg_2_tonne*X_CN
      } else{
        xx<-apply(thisData * thisVol[nlayers,,], 2, sum) * mg_2_tonne * X_CN
      }
      tracersArray[r,t,1:min(this_nts, max_nts)]<-xx[1:min(this_nts, max_nts)]
      if(length(thisNumCohorts)>0){
        if(thisNumCohorts>1){
          for(c in 1:thisNumCohorts){
            thisTracer<-paste(thisName,c,"_Nums",sep=""); thisData<-ncvar_get(ThisNC.nc, thisTracer)
            xx<-apply(thisData,3,sum, na.rm=TRUE)
            storeNumbersArray[r,t,c,1:min(this_nts,max_nts)]<-xx[1:min(this_nts,max_nts)]
            thisTracer<-paste(thisName,c,"_ResN",sep=""); thisData<-ncvar_get(ThisNC.nc, thisTracer)
            xx<-apply(thisData,3,nonZeroMean)
            thisTracer<-paste(thisName,c,"_StructN",sep=""); thisData<-ncvar_get(ThisNC.nc, thisTracer)
            yy<-apply(thisData,3,nonZeroMean)
            storeWeightsArray[r,t,c,1:min(this_nts,max_nts)]<-((xx+yy)*mg_2_tonne * X_CN)[1:min(this_nts,max_nts)]
          }
        }
      }
    }
  }
}

## dump it out so can just read it back in to use

# tracers relative to the base run - so can check for convergence
relTracersArray <- 0*tracersArray
for(r in 1:nruns){
  relTracersArray[r,,] <- 100*(tracersArray[r,,]- tracersArray[baseRunIndex,,])/ tracersArray[baseRunIndex,,]
}

## output so can read in elsewhere
save(list=c("relTracersArray", "tracers2plot", "tracersArray", "this_out", "runCols", "timeList", "nts_list"), file=paste(this_path,thisDesc,"modelTracers",sep=""))
# load(paste(this_path,thisDesc,"modelTracers",sep=""));

par(mfrow=c(2,1))
plot(tracersArray[r,1,], type="l")
points(tracersArray[baseRunIndex, 1, ], type="l", lty=2)
plot(relTracersArray[r,1,], type="l")
abline(h=0,col="red")

meanChangeByTimeRun <- apply(relTracersArray, c(1,3), mean, na.rm=TRUE)
maxChangeByTimeRun <- apply(relTracersArray, c(1,3), max, na.rm=TRUE)
minChangeByTimeRun <- apply(relTracersArray, c(1,3), min, na.rm=TRUE)
plot(meanChangeByTimeRun[2,1:nts_list[[baseRunIndex]]], type="l", ylim=c(min(meanChangeByTimeRun[1:nts_list[[baseRunIndex]]]), max(meanChangeByTimeRun[1:nts_list[[baseRunIndex]]])))
for(r in 1:nruns){
  points(meanChangeByTimeRun[r,], type="l", col=runCols[r])
}
makeBlankPlot()
legend(legend=this_out, col=runCols, lty=1, x="center", bty="n", ncol=2)

pdf(paste(plotPath, thisDesc, "CompareRunBounds_ltd.pdf", sep=""))
par(mfrow=c(3,3))
for(r in 1:nruns){
  if(r !=baseRunIndex){
    thisRunName<-this_out[r]; 
    xx<-gsub(thisRunTemplate,"", thisRunName); 
    temp <- gsub("Up|Down", "", xx); thisDir<-gsub(temp, "", xx)
    if(nchar(temp)==2){
      thisNum <- as.double(temp)*10
    }else{
      thisNum <- as.double(temp)
    }
    thisText<-paste(thisNum, "% ", thisDir, sep="")
    # plot(meanChangeByTimeRun[r,], type="l", ylim=c(min(minChangeByTimeRun[r,]), max(maxChangeByTimeRun[r,])), col=runCols[r])
    plot(meanChangeByTimeRun[r,], type="l", ylim=c(-100,100), col=runCols[r], ylab="Percentage change from Base run", xlab="Year")
    thisPlot_nts <- min(nts_list[[baseRunIndex]], nts_list[[r]])
    thisX<-seq(1, thisPlot_nts)
    thisY1<-minChangeByTimeRun[r,1:thisPlot_nts]; thisY2<-rev(maxChangeByTimeRun[r,1:thisPlot_nts])
    polygon(x=c(thisX, rev(thisX)), y=c(thisY1, thisY2), col=paste(runCols[r], "88", sep=""), border=NA)
    for(t in 1:ntracers){
      points(relTracersArray[r,t,1:thisPlot_nts], type="l", col=paste(runCols[r], "88", sep=""))
    }
    mtext(thisText)
    abline(h=0, col="black"); abline(h=c(-10,10), col="black", lty=2)
  }
}
dev.off()

## just age-structured groups
as_meanChangeByTimeRun <- apply(relTracersArray[,ageIndex,], c(1,3), mean, na.rm=TRUE)
as_maxChangeByTimeRun <- apply(relTracersArray[,ageIndex,], c(1,3), max, na.rm=TRUE)
as_minChangeByTimeRun <- apply(relTracersArray[,ageIndex,], c(1,3), min, na.rm=TRUE)

# ts<-60; r=1
# test<-relTracersArray[r,ageIndex,ts]; as_minChangeByTimeRun[r,ts]
# min(test)

pdf(paste(plotPath, thisDesc, "CompareRunBounds_as_ltd.pdf", sep=""))
par(mfrow=c(3,3))
for(r in 1:nruns){
  if(r !=baseRunIndex){
    thisRunName<-this_out[r]; 
    xx<-gsub(thisRunTemplate,"", thisRunName); 
    temp <- gsub("Up|Down", "", xx); thisDir<-gsub(temp, "", xx)
    if(nchar(temp)==2){
      thisNum <- as.double(temp)*10
    }else{
      thisNum <- as.double(temp)
    }
    thisText<-paste(thisNum, "% ", thisDir, sep="")
    # plot(meanChangeByTimeRun[r,], type="l", ylim=c(min(minChangeByTimeRun[r,]), max(maxChangeByTimeRun[r,])), col=runCols[r])
    plot(as_meanChangeByTimeRun[r,], type="l", ylim=c(-200,200), col=runCols[r], ylab="Percentage change from Base run", xlab="Year")
    thisPlot_nts <- min(nts_list[[baseRunIndex]], nts_list[[r]])
    thisX<-seq(1, thisPlot_nts)
    thisY1<-as_minChangeByTimeRun[r,1:thisPlot_nts]; thisY2<-rev(as_maxChangeByTimeRun[r,1:thisPlot_nts])
    polygon(x=c(thisX, rev(thisX)), y=c(thisY1, thisY2), col=paste(runCols[r], "88", sep=""), border=NA)
    for(t in 1:nageGroups){
      this_t <- (1:ntracers)[ageIndex][t]
      points(relTracersArray[r,this_t,1:thisPlot_nts], type="l", col=paste(runCols[r], "88", sep=""))
    }
    mtext(thisText)
    abline(h=0, col="black"); abline(h=c(-10,10), col="black", lty=2)
  }
}
dev.off()

## just not-age-structured groups
nas_meanChangeByTimeRun <- apply(relTracersArray[,!ageIndex,], c(1,3), mean, na.rm=TRUE)
nas_maxChangeByTimeRun <- apply(relTracersArray[,!ageIndex,], c(1,3), max, na.rm=TRUE)
nas_minChangeByTimeRun <- apply(relTracersArray[,!ageIndex,], c(1,3), min, na.rm=TRUE)

ts<-1; r=1
test<-relTracersArray[r,!ageIndex,ts]; as_minChangeByTimeRun[r,ts]
min(test)
tracers2plot[!ageIndex][test<20]

pdf(paste(plotPath, thisDesc, "CompareRunBounds_not_as_ltd.pdf", sep=""))
par(mfrow=c(3,3))
for(r in 1:nruns){
  if(r !=baseRunIndex){
    thisRunName<-this_out[r]; 
    xx<-gsub(thisRunTemplate,"", thisRunName); 
    temp <- gsub("Up|Down", "", xx); thisDir<-gsub(temp, "", xx)
    if(nchar(temp)==2){
      thisNum <- as.double(temp)*10
    }else{
      thisNum <- as.double(temp)
    }
    thisText<-paste(thisNum, "% ", thisDir, sep="")
    thisMax <- 2*thisNum; thisMin <- (-1)*thisMax
    thisPlot_nts <- min(nts_list[[baseRunIndex]], nts_list[[r]])
    plot(nas_meanChangeByTimeRun[r,], type="l", ylim=c(thisMin, thisMax), col=runCols[r], ylab="Percentage change from Base run", xlab="Year")
    thisX<-seq(1, thisPlot_nts)
    thisY1<-nas_minChangeByTimeRun[r,1:thisPlot_nts]; thisY2<-rev(nas_maxChangeByTimeRun[r,1:thisPlot_nts])
    polygon(x=c(thisX, rev(thisX)), y=c(thisY1, thisY2), col=paste(runCols[r], "88", sep=""), border=NA)
    for(t in 1:nageGroups){
      this_t <- (1:ntracers)[!ageIndex][t]
      points(relTracersArray[r,this_t,1:thisPlot_nts], type="l", col=paste(runCols[r], "88", sep=""))
    }
    mtext(thisText)
    abline(h=0, col="black"); abline(h=c(-10,10), col="black", lty=2)
  }
}
dev.off()

####
## which group/s goes wildly out in the last run..?
absMaxByGroup <- apply(abs(relTracersArray[,,40:50]), c(1,2), max, na.rm=TRUE)
pdf(paste(plotPath, thisDesc, "RankTracers.pdf", sep=""))
par(mfrow=c(4,2))
for(r in 1:nruns){
  if(r != baseRunIndex){
    thisOrder <- order(absMaxByGroup[r,], decreasing = TRUE)
    ## only do those greater than 5%
    thisIndex <- absMaxByGroup[r,thisOrder]>10
    if(length(absMaxByGroup[thisIndex])>0){
      plot(absMaxByGroup[r,thisOrder][1:15], type="h", xaxt="n", xlab="")
      par(las=2)
      axis(at=(1:ntracers)[1:15], labels=tracers2plot[thisOrder][1:15], side=1)
      par(las=1)
      mtext(this_out[r])
    }
  }
}
dev.off()




