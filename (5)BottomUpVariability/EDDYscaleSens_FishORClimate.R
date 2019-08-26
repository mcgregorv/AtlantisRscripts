## catch histories were created in setUpScenario_ts_files_versionC1.R - they are full historic + 50 year future runs
## biol.prm files for eddy_scale are  created first

thisPath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

version<-"C"
runPath<-paste(thisPath,"eddySens\\",sep="")

groupsDF<-read.csv(paste(thisPath,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

plotPath <- paste(DIR$'Reports',"//(02)ClimateChange//Figures//InitialPP//", sep="") ## overwrite report figures
plotPath<-paste(thisPath,"Figures\\eddySens\\",version,"",sep="")

## read in trophic levels so can plot by these
trophicLevels<-read.csv(paste(thisPath, "base\\EWEbase\\CRAMGroupsTL.csv", sep=""))

load(paste(thisPath,"eddySens\\modelTracers_V",version,sep="")) ## brings in relTracers and storeTracers
## dims are run, tracer, layer, box, time for both
nruns<-dim(lookup_df)[1]

baseModelPath<-paste(thisPath,"base\\outputBase\\",sep="")
ThisNC.nc<-nc_open(paste(baseModelPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
nts<-dim(thisVol)[3]; nboxes<-dim(thisVol)[2]; nlayers<-dim(thisVol)[1]
allTracers<-names(ThisNC.nc$var); 
expTracers<-allTracers[grep("Nums|StructN|ResN", allTracers, invert = TRUE)] # take out numbers and individual weights
ntracers<-length(expTracers)

## relative tracers by timestep
baseRunIndex<-3
thisCex<-1.2

## what has the strongest effect - fishing or primary production..?
## check out key fisheries species groups first
comNames<-str_trim(groupsDF$Name[groupsDF$Code %in% c( "HOK", "HAK", "IVS", "LIN", "ORH", "PFM", "SSO")], side="both"); 


ncoms<-length(comNames)
comTracers<-paste(comNames,"_N", sep="")
comTracerIndex<-grep(paste(comTracers,collapse="|"), expTracers)
storeComBiom<-array(NA, dim=c(nruns, ncoms, dim(storeTracers)[5]))
for(c in 1:ncoms){
  thisBiomasses<-storeBiomass[,comTracerIndex[c],]
  thisBaseBiomass<-storeBiomass[baseRunIndex,comTracerIndex[c],]
  ## check the base biomass for very small levels (such as if the base run has crashed)
  ## and in this case, replace with NA so doesn't come into the analyses
  crashIndex<-(thisBaseBiomass/thisBaseBiomass[1])<1e-8
  thisBaseBiomass[crashIndex]<-NA
  if(sum(crashIndex)>0){cat(paste("check ",expTracers[comTracerIndex[c]], "\n"))}
  for(r in 1:nruns){
    storeComBiom[r,c,]<-( thisBiomasses[r,] - thisBaseBiomass)/thisBaseBiomass
  }
}

finalRel<-storeComBiom[,,202]
colnames(finalRel)<-expTracers[comTracerIndex]

catchSens<-unlist(lapply(lookup_df$Scenario, get_first_number)); ppSens<-c("PP base", "PP reduced")[lookup_df$Eddy]
runDescr<-paste("Catch ",catchSens,"%, ", ppSens, sep="")

rownames(finalRel)<-runDescr
runCatches<-unlist(lapply(gsub("%","",runDescr), get_first_number))
catchOrder<-order(runCatches)

thisCex<-1.5

# pdf(paste(plotPath,"TopFished_proportionChange.pdf", sep=""))
par(lend=1, las=1, mar=c(6,6,1,1), mfrow=c(1,1), mfrow=c(3,3))
for(c in 1:ncoms){
  thisCom<-comNames[c]
  setEPS()
  postscript(paste(plotPath,"BiomassResponse_",thisCom,".eps", sep=""), height=5)
  par(mar=c(8,7,2.5,0.5))
  plot(finalRel[catchOrder,c], type="h", lwd=5, col=c(myBlue,myOrange), cex.lab=thisCex, cex.axis=thisCex
       , xaxt="n", xlab="", ylab="Proportional change\nin biomass")
  par(las=2)
  axis(at=seq(1.5,nruns,by=2), labels = paste(unique(runCatches[catchOrder]),"%",sep=""), cex.lab=thisCex, cex.axis=thisCex, side=1)
  par(las=1)
  mtext(expTracers[comTracerIndex[c]], cex=thisCex, side=3, adj=0)
  mtext("Percentage of status quo catch", side=1, adj=0.5, line=6.5, cex=thisCex)
  dev.off()
}
setEPS()
postscript(paste(plotPath,"BiomassResponse_LEGEND.eps", sep=""), height=5)
par(mar=c(0,0,0,0), lend=1)
makeBlankPlot()
legend(legend=c("Status quo primary production", "Reduced primary production"), cex=thisCex, col=c(myBlue, myOrange), lwd=7, seg.len=3, x="center", bty="n")
dev.off()
# dev.off()

### bring in catches for timeseries plots
catchHistories<- read.csv(paste(DIR$'Tables',"ALLcatchHistories.csv",sep=""))
yearAveIndex<-catchHistories[,1] %in% 2010:2014
catchYears<-catchHistories[,1]

#look at the time series plots too
# par(mfrow=c(3,2), mar=c(4,4,1,1))
timeIndex<-120:202; x_labs<-timeIndex + 1865
biomassYears<-(1865:(1865 + dim(storeBiomass)[3]-1))[timeIndex]; biomassYears<-biomassYears[!is.na(biomassYears)]
futureCatchYears<-biomassYears[biomassYears>max(catchYears)]
for(c in 1:ncoms){
  thisCom<-comNames[c]
  thisCode<-groupsDF$Code[groupsDF$Name == thisCom]
  thisHistCatch<-catchHistories[,c(as.character(thisCode))]/1000; thisAveCatch<-mean(thisHistCatch[yearAveIndex])
  thisCatchAxis<-pretty(seq(0,max(thisHistCatch, na.rm=TRUE), length.out=5))
  futureCatch<-rep(thisAveCatch, length(futureCatchYears))
  thisBiomasses<-storeBiomass[,comTracerIndex[c],]
  plotBiomass<-thisBiomasses[1,timeIndex]; plotBiomass<-plotBiomass[!is.na(plotBiomass)]
  thisYmax<-max(thisBiomasses[,timeIndex], na.rm=TRUE)
  setEPS()
  postscript(paste(plotPath,"BiomassTimeSeries_",thisCom,".eps", sep=""), height=4.5)
  # png(paste(plotPath,"BiomassTimeSeries_",thisCom,".png", sep=""), height=400, width=500)
  par(mar=c(4.5,4.5,2,4.5), lend=1)
  plot(x=biomassYears, y=thisHistCatch[match(biomassYears, catchYears)], type="h", lwd=2, col=myGrey, xaxt="n", yaxt="n", ylab="", xlab="")
  points(x=futureCatchYears, y=futureCatch, type="h", lwd=2, col=myGrey)
  axis(at=thisCatchAxis, labels=thisCatchAxis, side=4, col.axis=myGrey, font=2, cex.axis=thisCex)
  mtext("Catch (tonnes)", side=4,adj=0.5, line=3, col=myGrey, font=2, cex=thisCex)
  par(new=TRUE)
  plot(x=biomassYears, y=plotBiomass, type="n", cex.lab=thisCex, cex.axis=thisCex, ylim=c(0,thisYmax), ylab="Biomass (tonnes)", xlab="Year")
  # axis(at=1:(length(timeIndex)), labels = x_labs, side=1)
  mtext(expTracers[comTracerIndex[c]], side=3, adj=0, cex=thisCex)
  for(r in 1:nruns){
    thisCol<-c(myBlue, myOrange)[lookup_df$Eddy[r]]
    points(x=biomassYears, y=thisBiomasses[r,timeIndex], type="l", lwd=2, col=thisCol)
  }
  dev.off()
}
setEPS()
postscript(paste(plotPath,"BiomassTimeSeries_LEGEND.eps", sep=""), height=5)
par(mar=c(0,0,0,0), lend=1)
makeBlankPlot()
legend(legend=c("Status quo primary production", "Reduced primary production","Catches"), cex=thisCex, col=c(myBlue, myOrange,myGrey), lwd=c(3,3,5), seg.len=3, x="center", bty="n")
dev.off()


## by TL do proportional change
## read in trophic levels so can plot by these
trophicLevels<-read.csv(paste(thisPath, "base\\EWEbase\\CRAMGroupsTL.csv", sep=""))
trophicLevelsIso<-read.csv(paste(thisPath, "base\\EWEbase\\CRAM_trophicLevels_isotopes.csv", sep=""))

## do status quo with and without PP reduced
nfruns<-nruns/2
timeIndex<-180:190
timeIndex<-192:202

catchLevels<-unique(unlist(lapply(lookup_df$Scenario, get_first_number)))

cleanNAs<-function(x){
  y<-x
  if(length(x)==0){y<-NA}
  else if(is.na(x)){y<-NA}
  else if(x==""){y<-NA}
  else if(x=="NaN"){y<-NA}
  else if(x=="Inf"){y<-NA}
  return(y)
}

thisTrophicLevels<-seq(0.5,5.5, by=0.5); 
thisTrophicLevels<-seq(1,5, by=1); 
nTLs<-length(thisTrophicLevels)
sumByTL<-array(NA, dim=c(nfruns, (nTLs+2))); minByTL<-sumByTL; meanByTL <- sumByTL; maxByTL<-sumByTL
medByTL<-sumByTL
for(f in 1:nfruns){
  run1<-2*f-1                  
  runIndex<-run1:(run1+1)
  thisBiomasses<-storeBiomass[runIndex,,timeIndex]
  ## total biomasses
  thisTotalBiomasses<-apply(thisBiomasses, c(1,3), sum, na.rm=TRUE)
  
  thisRelBiomass<-(thisTotalBiomasses[2,] - thisTotalBiomasses[1,])/thisTotalBiomasses[1,]
  thisRelBiomass<-unlist(lapply(thisRelBiomass, cleanNAs))
  sumByTL[f,7]<-mean(thisRelBiomass)
  
  ppcodes<-c("DF", "MA", "MB", "PL", "PS")
  primaryProducers<-c(paste(groupsDF$Name[groupsDF$Code %in% ppcodes],"_N", sep=""))
  testTracerIndex<-grep(paste(primaryProducers,collapse="|", sep=""), expTracers)
  expTracers[testTracerIndex]
  
  thisTotalBiomasses<-apply(thisBiomasses[,testTracerIndex,], c(1,3), sum, na.rm=TRUE)
  
  thisRelBiomass<-(thisTotalBiomasses[2,] - thisTotalBiomasses[1,])/thisTotalBiomasses[1,]
  thisRelBiomass<-unlist(lapply(thisRelBiomass, cleanNAs))
  
  sumByTL[f,1]<-mean(thisRelBiomass, na.rm=TRUE)
  
  thisIndBiomasses<-thisBiomasses[,testTracerIndex,]
  thisIndRel<- (thisIndBiomasses[2,,] - thisIndBiomasses[1,,])/thisIndBiomasses[1,,]
  thisMeanByGroup<- apply(thisIndRel, 1, mean, na.rm=TRUE)
  minByTL[f, (1)] <- min(thisMeanByGroup, na.rm=TRUE) ## the mean is wrt time (2045:2055)
  maxByTL[f, (1)] <- max(thisMeanByGroup, na.rm=TRUE)
  meanByTL[f, (1)] <- mean(thisMeanByGroup, na.rm=TRUE)
  medByTL[f, 1] <- median(thisMeanByGroup, na.rm=TRUE)

  for(l in 1:nTLs){
    thisTL<-thisTrophicLevels[l]
    levelNames<-paste(str_trim(trophicLevels$Name[trunc(trophicLevels$TL)==thisTL],side="both"),"_N", sep="")
    levelNames <- paste("^",levelNames[!(levelNames %in% primaryProducers)], sep="")
    testTracerIndex<-grep(paste(levelNames,collapse="|", sep=""), expTracers)
    expTracers[testTracerIndex]
    
    thisIndBiomasses<-thisBiomasses[,testTracerIndex,]
    thisIndRel<- (thisIndBiomasses[2,,] - thisIndBiomasses[1,,])/thisIndBiomasses[1,,]
    thisMeanByGroup<- apply(thisIndRel, 1, mean, na.rm=TRUE)
    
    thisTotalBiomasses<-apply(thisBiomasses[,testTracerIndex,], c(1,3), sum, na.rm=TRUE)
    
    thisRelBiomass<-(thisTotalBiomasses[2,] - thisTotalBiomasses[1,])/thisTotalBiomasses[1,]
    thisRelBiomass<-unlist(lapply(thisRelBiomass, cleanNAs))
    
    sumByTL[f, (l+1)] <- mean(thisRelBiomass, na.rm=TRUE)
    minByTL[f, (l+1)] <- min(thisMeanByGroup, na.rm=TRUE) ## the mean is wrt time (2045:2055)
    maxByTL[f, (l+1)] <- max(thisMeanByGroup, na.rm=TRUE)
    meanByTL[f, (l+1)] <- mean(thisMeanByGroup, na.rm=TRUE)
    medByTL[f, (l+1)] <- median(thisMeanByGroup, na.rm=TRUE)
    
  }
}

colByFish<-colorRampPalette(colors=c(myGreen, myAqua, myBlue, myRed))(nfruns)

orderFish<-order(catchLevels)

xbase<-seq(1,(nTLs +2))
thisYmin<-min(medByTL, na.rm=TRUE)
yaxis_at<-pretty(seq(0, thisYmin, length.out=5)); yaxis_lab<- 100 * abs(yaxis_at)
setEPS()
postscript(paste(plotPath,"PropDecreaseByTL.eps", sep=""), height=5)
par(lend=1, las=1)
plot(medByTL[1,1:(nTLs+1)], ylim=c(thisYmax, thisYmin), cex.lab=thisCex, cex.axis=thisCex, type="n",ylab="Percentage biomass decrease", xlab="Trophic level", xaxt="n", yaxt="n")
axis(at=1:(nTLs +1), labels=c("PP", thisTrophicLevels), side=1, cex.axis=thisCex, cex.lab=thisCex)
axis(at=yaxis_at, labels = yaxis_lab, side=2, cex.lab=thisCex, cex.axis=thisCex)
for(f in 1:nfruns){
  thisShift<-f*0.05 -0.1
  points(x=xbase + thisShift, y=medByTL[orderFish[f],], type="h", lwd=5, col=colByFish[f])
}
legend(legend=catchLevels[orderFish], col=colByFish, lwd=3, seg.len=3, x="topright", bty="n")
dev.off()

thisYmax<-max(maxByTL, na.rm=TRUE); thisYmin<-min(minByTL, na.rm=TRUE)
yaxis_at<-pretty(seq(0, thisYmin, length.out=5)); yaxis_lab<- 100 * abs(yaxis_at)
plot(meanByTL[1,1:(nTLs+1)], ylim=c(thisYmax, thisYmin), cex.lab=thisCex, cex.axis=thisCex, type="n",ylab="Percentage biomass decrease", xlab="Trophic level", xaxt="n", yaxt="n")
axis(at=1:(nTLs +1), labels=c("PP", thisTrophicLevels), side=1, cex.axis=thisCex, cex.lab=thisCex)
axis(at=yaxis_at, labels = yaxis_lab, side=2, cex.lab=thisCex, cex.axis=thisCex)
for(f in 1:nfruns){
  thisShift<-f*0.05 -0.1
  points(x=xbase + thisShift, y=meanByTL[orderFish[f],], col=colByFish[f], type="l", lwd=2)
  segments(x0=xbase + thisShift, y0=minByTL[orderFish[f],], x1=xbase + thisShift, y1=maxByTL[orderFish[f],], col=paste(colByFish[f],"88", sep=""), lwd=2.5)
}



## for each group, get rel biomass - then plot against respective trophic level
storeRelByGroup<-array(NA, dim=c(nfruns, ntracers))
for(f in 1:nfruns){
  run1<-2*f-1                  
  runIndex<-run1:(run1+1)
  thisBiomasses<-storeBiomass[runIndex,,timeIndex]
  for(t in 1:ntracers){
    thisTracer<-expTracers[t]
    if(length(grep("_N", thisTracer))>0){
      if(!(thisTracer %in% primaryProducers)){
        tempBiomass<-thisBiomasses[,t,]
        thisRelBiomass<- (tempBiomass[2,] - tempBiomass[1,]) / tempBiomass[2,]
        storeRelByGroup[f,t]<- mean(thisRelBiomass, na.rm=TRUE)
      }
    }
  }
}
## need trophic level in same order as expTracers
expTLs<-rep(NA, ntracers)
for(t in 1:ntracers){
  thisTracer <- expTracers[t]
  if(length(grep("_N", thisTracer))>0){
    thisName<-gsub("_N", "", thisTracer)
    if(!(thisTracer %in% primaryProducers)){
      thisCode<-groupsDF$Code[grep(thisName, groupsDF$Name)]
      if(thisName=="Zoo"){thisCode<-"ZM"}
      thisTL<-trophicLevels$TL[trophicLevels$Code==thisCode]
      expTLs[t]<-thisTL
    }
  }
}  
thisMin<-min(storeRelByGroup, na.rm=TRUE); thisMax<-max(storeRelByGroup, na.rm=TRUE)
plot(x=seq(0,6), y=rep(0,7), type="n", ylab="Percentage biomass decrease", xlab="Trophic level", ylim=c(thisMax, thisMin), xlim=c(0.5,5.5))
for(f in 1:6){
  points(x=expTLs, y=storeRelByGroup[f,], col=colByFish[f], pch=20, cex=1.4)
}
## which ones are TL 1?
testIndex<-round(expTLs)==1 & !is.na(expTLs); expTracers[testIndex]
tempBiomass<-thisBiomasses[,testIndex,]
totalBiomasses<-apply(tempBiomass, c(1,3), sum, na.rm=TRUE)
relTot<- (totalBiomasses[2,] - totalBiomasses[1,])/ totalBiomasses[2,]
relInd<- (tempBiomass[2,,] - tempBiomass[1,,])/ tempBiomass[2,,]

## total biomasses
for(f in 1:nfruns){
  run1<-2*f-1                  
  runIndex<-run1:(run1+1)
  thisBiomasses<-storeBiomass[runIndex,,timeIndex]
  
  thisTotalBiomasses<-apply(thisBiomasses, c(1,3), sum, na.rm=TRUE)
  
  thisRelBiomass<-(thisTotalBiomasses[2,] - thisTotalBiomasses[1,])/thisTotalBiomasses[1,]
  thisRelBiomass<-unlist(lapply(thisRelBiomass, cleanNAs))
  sumByTL[f,7]<-mean(thisRelBiomass)
  
  ppcodes<-c("DF", "MA", "MB", "PL", "PS")
  primaryProducers<-c(paste(groupsDF$Name[groupsDF$Code %in% ppcodes],"_N", sep=""))
  testTracerIndex<-grep(paste(primaryProducers,collapse="|", sep=""), expTracers)
  expTracers[testTracerIndex]
  
  thisTotalBiomasses<-apply(thisBiomasses[,testTracerIndex,], c(1,3), sum, na.rm=TRUE)
  
  thisRelBiomass<-(thisTotalBiomasses[2,] - thisTotalBiomasses[1,])/thisTotalBiomasses[1,]
  thisRelBiomass<-unlist(lapply(thisRelBiomass, cleanNAs))
  
  sumByTL[f,1]<-mean(thisRelBiomass, na.rm=TRUE)
  
  for(l in 1:5){
    levelNames<-paste(str_trim(trophicLevels$Name[trunc(trophicLevels$TL)==l],side="both"),"_N", sep="")
    levelNames <- paste("^",levelNames[!(levelNames %in% primaryProducers)], sep="")
    testTracerIndex<-grep(paste(levelNames,collapse="|", sep=""), expTracers)
    expTracers[testTracerIndex]
    
    thisTotalBiomasses<-apply(thisBiomasses[,testTracerIndex,], c(1,3), sum, na.rm=TRUE)
    
    thisRelBiomass<-(thisTotalBiomasses[2,] - thisTotalBiomasses[1,])/thisTotalBiomasses[1,]
    thisRelBiomass<-unlist(lapply(thisRelBiomass, cleanNAs))
    
    sumByTL[f, (l+1)] <- mean(thisRelBiomass, na.rm=TRUE)
    
  }
}

yaxis_at<-pretty(seq(0, min(sumByTL), length.out=5)); yaxis_lab<- 100 * abs(yaxis_at)
# 
# thisCex<-1.2
# setEPS()
# postscript(paste(plotPath,"PropDecreaseByTL.eps", sep=""), height=5)
# par(lend=1, las=1)
# plot(x=seq(1,(nTLs +2)), y=rep(0,(nTLs +2)), type="n", ylim=c(max(sumByTL, na.rm=TRUE), min(sumByTL, na.rm=TRUE)), cex.axis=thisCex, yaxt="n", xaxt="n", xlab="Trophic level", ylab="Percentage decrease", cex.lab=thisCex)
# axis(at=1:(nTLs +2), labels=c("PP", thisTrophicLevels, "Combined"), side=1, cex.axis=thisCex, cex.lab=thisCex)
# axis(at=yaxis_at, labels = yaxis_lab, side=2, cex.lab=thisCex, cex.axis=thisCex)
# for(f in 1:nfruns){
#   thisShift<-0.05*f-0.15
#   points(x=seq(1,(nTLs +2))+thisShift, y=sumByTL[orderFish[f],], type="h", lwd=4, col=colByFish[f])
# }
# dev.off()

makeBlankPlot()
legend(legend=paste(catchLevels[orderFish]," %", sep=""), col=colByFish, lwd=3, x="topright", bty="n", title="Catch levels")

trophicLevels$TL[trophicLevels$Code=="DC"]<-NA
timeIndex<-150:202
runIndex<-1:2; thisTitle<-"Zero catches" ## set the first one to be the base
# runIndex<-3:4; thisTitle<-"Status quo catches" ## set the first one to be the base
# runIndex<-5:6; thisTitle<-"Catches increased by 20%" ## set the first one to be the base
# runIndex<-7:8; thisTitle<-"Catches increased by 50%"
# runIndex<-9:10; thisTitle<-"Catches decreased by 50%"
# runIndex<-11:12; thisTitle<-"Catches decreased by 20%"
relByTL<-rep(0,6); relByTLtime<-array(NA, dim=c(6,length(timeIndex)))
for(l in 0:5){
  ## if l==0, grab primary producers - and otherwise, make sure don't grab primary producers
  if(l==0){
    thisCodes<-c("DF", "MA", "MB", "PL", "PS")
    primaryProducers<-c(paste(groupsDF$Name[groupsDF$Code %in% thisCodes],"_N", sep=""))
    thisTracerIndex<-grep(paste(primaryProducers,collapse="|", sep=""), expTracers)
  }else{
    thisCodes<-trophicLevels$Code[trunc(trophicLevels$TL)==l]
    thisCodes <- thisCodes[!(thisCodes %in% ppcodes)]
    thisTracerIndex<-grep(paste(groupsDF$Name[groupsDF$Code %in% thisCodes], collapse="|"), expTracers)
    
  }
  thisTracers<-expTracers[thisTracerIndex]
  thisData<-storeBiomass[runIndex,thisTracerIndex,timeIndex]; ntt<-length(thisTracers)
  # get total biomass for this trophic level in each run
  totalBiom<- apply(thisData, c(1,3), sum, na.rm=TRUE)
  thisRel<- (totalBiom[1,] - totalBiom[2, ])/ totalBiom[1,]
  relByTLtime[l+1,]<-thisRel
}

colnames(relByTLtime)<-timeIndex + 1865
toPlot<-melt(relByTLtime); colnames(toPlot)<-c("TL", "Year","Value")
toPlot$TL<-toPlot$TL-1
bp<-ggplot(data = toPlot, aes(x = Year, fill = TL, y = Value)) + 
  geom_bar(stat = 'identity')

pdf(paste(plotPath,"comparePPruns_", gsub(" |,|%","",thisTitle), ".pdf", sep=""), height=4, width=5)
par(mar=c(4,3.5,1,1))
bp + labs(y="Proportional decrease", x="Year", title=thisTitle) + theme_igray() +
  theme(axis.text=element_text(size=12, angle=0),axis.title=element_text(size=12)) + 
  guides(fill=guide_legend(title="Trophic Level"))  
dev.off() 

## how does the proportional decerase by trophic level compare to that for PP total?
nfishRuns<-nruns/2
timeIndex<-150:202
nfruns<-4

propReduced<-array(NA, dim=c(nfruns,5,length(timeIndex)))
for(f in 1:nfruns){
  run1<-2*f-1
  thisBaseBiomass <- storeBiomass[run1:(run1+1),,]
  biomassByTL<-array(NA, dim=c(2, 6, length(timeIndex)))
  for(l in 0:5){
    ## if l==0, grab primary producers - and otherwise, make sure don't grab primary producers
    if(l==0){
      thisCodes<-c("DF", "MA", "MB", "PL", "PS")
      primaryProducers<-c(paste(groupsDF$Name[groupsDF$Code %in% thisCodes],"_N", sep=""))
      thisTracerIndex<-grep(paste(primaryProducers,collapse="|", sep=""), expTracers)
    }else{
      thisCodes<-trophicLevels$Code[trunc(trophicLevels$TL)==l]
      thisCodes <- thisCodes[!(thisCodes %in% ppcodes)]
      thisTracerIndex<-grep(paste(groupsDF$Name[groupsDF$Code %in% thisCodes], collapse="|"), expTracers)
      
    }
    thisTracers<-expTracers[thisTracerIndex]
    thisData<-thisBaseBiomass[,thisTracerIndex,timeIndex]; ntt<-length(thisTracers)
    totalBiom<- apply(thisData, c(1,3), sum, na.rm=TRUE)
    biomassByTL[,(l+1),]<-totalBiom
    if(l > 0){
      propReduced[f,l,]<-(totalBiom[1,] - totalBiom[2,])/ totalBiom[1,]
    }
  }
}
par(mfrow=c(3,2))
for(l in 1:5){
  thisMax<-max(propReduced[,l,], na.rm=TRUE); thisMin<-min(propReduced, na.rm=TRUE)
  plot(propReduced[1,l,], type="n", ylim=c(thisMin,thisMax))
  for(r in 1:nfruns){
    points(propReduced[r,l,], type="l",lty=r)
  }
}

colByTL<-colorRampPalette(colors=c(myGold, "red", myPurple,myBlue, myAqua))(6)

# which group/s in TL 1 drop in the reduced PP runs 40 years in?
l=1
runIndex<-1:2
thisCodes<-trophicLevels$Code[trunc(trophicLevels$TL)==l]
thisCodes <- thisCodes[!(thisCodes %in% ppcodes)]
thisTracerIndex<-grep(paste(paste(groupsDF$Name[groupsDF$Code %in% thisCodes], "_N",sep=""), collapse="|"), expTracers)
thisTracers<-expTracers[thisTracerIndex]
thisData<-storeBiomass[runIndex,thisTracerIndex,timeIndex]; ntt<-length(thisTracers)

thisTotals<-apply(thisData, c(1,3), sum, na.rm=TRUE)
thisTotalRel<- (thisTotals[1,] - thisTotals[2, ])/ thisTotals[1,]
plot(thisTotalRel, type="l")
thisMax<-max(thisTotals)
plot(thisTotals[1,], type="l", col=myBlue, ylim=c(0,thisMax))
points(thisTotals[2,], type="l", col=myOrange)
legend(legend=c("Base","Reduced PP"), col=c(myBlue, myOrange), x="bottomleft", bty="n", lty=1)

for(t in 1:length(thisTracers)){
  thisMax<-max(thisData[,t,], na.rm=TRUE)
  plot(thisData[1,t,], type="l", col=myBlue, ylim=c(0,thisMax)); points(thisData[2,t,], type="l", col=myOrange)
  legend(legend=c("Base","Reduced PP"), col=c(myBlue, myOrange), x="bottomleft", bty="n", lty=1)
  mtext(thisTracers[t], side=3, adj=0)
}

test<-(thisData[1,,] - thisData[2,,])/thisData[1,,]
testCol<-colorRampPalette(colors=c(myGold,myGreen, myAqua, myBlue,myPurple))(length(thisTracers))
thisMax<-max(test, na.rm=TRUE); thisMin<-min(test, na.rm=TRUE)
par(mfrow=c(2,1))
plot(test[1,], type="n", ylim=c(thisMin, thisMax))
for(t in 1:dim(test)[1]){
  points(test[t,], type="l", lty=t, col=testCol[t])
}
makeBlankPlot()
legend(legend=thisTracers, lty=seq(1,length(thisTracers)), x="center", bty="n", col=testCol)

####################
# ## output biomass (B0) and perhaps at time, for Monique
# tracerIndex<-grep("_N", expTracers)
# outBiom<-storeBiomass[3,tracerIndex,]
# rownames(outBiom)<-expTracers[tracerIndex]
# colnames(outBiom)<-1865:(1865+202-1)
# write.csv(outBiom, paste(thisPath, "base\\EWEbase\\biomassAtlantisOut.csv", sep=""))
# ## and area
# areaDyn<-sum(thisVol[6,2:25,1])

