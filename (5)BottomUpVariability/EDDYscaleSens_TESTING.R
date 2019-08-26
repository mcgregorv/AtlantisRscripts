## catch histories were created in setUpScenario_ts_files_versionC1.R - they are full historic + 50 year future runs
## biol.prm files for eddy_scale are  created first

thisPath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

version<-"C"
runPath<-paste(thisPath,"eddySens\\",sep="")

groupsDF<-read.csv(paste(thisPath,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

plotPath<-paste(thisPath,"Figures", "\\eddySens\\PP_", sep="")

## set up DF with runs specified
## set up DF with runs specified
## set up DF with runs specified
scenarios<-c("All0catch") ## C
scenarios<-c("All0catch", "All50catch", "All100catch",  "All150catch") ## C
nscenarios<-length(scenarios)
eddysens<-c("EddyNOchange", "Eddy")
neddys<-length(eddysens); eddyIndex<-seq(1,neddys)
nruns<-nscenarios * neddys

lookup_df<-data.frame(array(NA, dim=c(nruns,3)))
colnames(lookup_df)<-c("Run", "Scenario", "Eddy")
lookup_df$Run<-seq(1,nruns)
lookup_df$Scenario<-rep( sort(rep(scenarios, neddys)))
lookup_df$Eddy<-rep(eddyIndex, (nscenarios))

# lookup_df<-lookup_df[!(lookup_df$Scenario == "All150catch"),]

nruns=dim(lookup_df)[1]

baseModelPath<-paste(thisPath,"base\\outputBase\\",sep="")
ThisNC.nc<-nc_open(paste(baseModelPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
nts<-dim(thisVol)[3]; nboxes<-dim(thisVol)[2]; nlayers<-dim(thisVol)[1]

nts_list<-rep(NA, nruns)
for(r in 1:nruns){
  cat(r,"--")
  thisRunPath<-paste(runPath,"outputEDDY",version,"_",r,"\\",sep="")
  runNC.nc<-nc_open(paste(thisRunPath, "output.nc",sep=""))
  runVol<-ncvar_get(runNC.nc,"volume"); run_nts<-dim(runVol)[3]
  nts_list[r]<-run_nts
}
nts_list
nts<-max(nts_list)
allTracers<-names(ThisNC.nc$var); 
expTracers<-allTracers[grep("Nums|StructN|ResN", allTracers, invert = TRUE)] # take out numbers and individual weights
ntracers<-length(expTracers)

storeTracers<-array(NA, dim=c(nruns, ntracers,  nlayers, nboxes, nts))

for(r in 1:nruns){
  cat(r,"--")
  thisRunPath<-paste(runPath,"outputEDDY",version,"_",r,"\\",sep="")
  runNC.nc<-nc_open(paste(thisRunPath, "output.nc",sep=""))
  runVol<-ncvar_get(runNC.nc,"volume"); run_nts<-dim(runVol)[3]
  for(t in 1:ntracers){
    thisTracer<-expTracers[t]
    thisData<-ncvar_get(runNC.nc,thisTracer)
    test<-grep("_N", thisTracer)
    if(length(dim(thisData))==3){
      storeTracers[r,t,,,1:min(run_nts, nts)]<-thisData[,,1:(min(run_nts, nts))]
    } else if(length(dim(thisData))==2 & dim(thisData)[1]==nboxes){
      storeTracers[r,t,nlayers,,1:min(run_nts,nts)]<-thisData[,1:(min(run_nts, nts))]
    }
  }
}

allBiomassTracers<-array(NA, dim=dim(storeTracers)[c(1,2,5)])
for(r in 1:nruns){
  runVol<-storeTracers[r,expTracers=="volume",,,]
  for(t in 1:ntracers){
    test1<-grep("_N", expTracers[t])
    if(length(test1)>0){
      thisData<-apply(storeTracers[r,t,,,] * runVol, 3, sum, na.rm=TRUE) * mg_2_tonne * X_CN
      allBiomassTracers[r,t,]<-thisData
    }
  }
}
zero2NA<-function(x){
  y<-x
  if(length(x)>0){
    if(!is.na(x)){
      if(x==0){
        y<-NA
      }
    }
  }
  return(y)
}
# 
# library(plot3D)
## 

colByRun<-colorRampPalette(colors=c(myGold,myGreen, myBlue,"midnightblue",myRed,"red"))(nruns)
ltyByRun<-rep(c(1,2), ceiling(nruns/2))
runDesc<-c("Zero catch", "Zero catch and reduce PP", "Base", "Reduce PP", "Increase catch", "Increase catch and reduce PP", "Reduce catch", "Reduce catch and PP")
# runDesc<-c("Zero catch", "Zero catch and reduce PP", "Base", "Reduce PP",  "Reduce catch", "Reduce catch and PP")

#overall, seems zero catch has more impact than half PP - how do the key fisheries species look..?
topFishCodes<-c("HOK", "ORH", "SSO", "LIN", "BOE", "PFM")
topFishTracers<-paste(str_trim(groupsDF$Name[groupsDF$Code %in% topFishCodes]), "_N", sep="")
testTracerIndex<-grep(paste(topFishTracers,collapse = "|"), expTracers)
ntt<-length(topFishTracers)
topBiomassData<-allBiomassTracers[,testTracerIndex,]

timeIndex<-130:202
year0<-1865; allYears<-1865:(1865+nts-1); thisYears<-allYears[timeIndex]

##just plot time series first
pdf(paste(plotPath,"Topfisheries_tracers.pdf", sep=""), height=7, width=5)
par(mfrow=c(4,2), mar=c(4,4,1,1), las=1)
for(t in 1:ntt){
  thisTracer<-expTracers[testTracerIndex][t]; thisTdata<-topBiomassData[,t,timeIndex]
  thisYmax<-max(thisTdata, na.rm=TRUE); thisYmin<-min(thisTdata, na.rm=TRUE)
  plot(x=thisYears, y=thisTdata[1,], type="n", ylim=c(0, thisYmax), ylab="Biomass (tonnes)", xlab="'Years"); mtext(thisTracer, side=3, adj=0)
  abline(v=2015,col=myGrey, lwd=2, lty=3)
  for(r in 1:nruns){
    points(x=thisYears, y=thisTdata[r,], type="l", lwd=2, col=colByRun[r], lty=ltyByRun[r])
  }
  
}
# par(mfrow=c(1,1))
makeBlankPlot()
legend(legend=runDesc, col=colByRun, lwd=2, x="center", bty="n")
dev.off()


## check fishing
r=1; thisRunPath<-paste(runPath,"outputEDDY",version,"_",r,"\\",sep="")
# 
# thisPath<-paste(paste(basePath,"SSRsens\\SSR_",Version,"\\outputSSR",Version,"1_",r,"\\", sep=""))
# thisOutFile<-paste(thisPath,"outputCATCH.nc", sep=""); catchNC.nc<-nc_open(thisOutFile); 


catchNC.nc<-nc_open(paste(thisRunPath, "outputCATCH.nc",sep="")); 
catchVars<-unique(sort(names(catchNC.nc$var)))
x<-grep("Orange", catchVars)
hokiCatchByCohort<-array(NA, dim=c(10, nts))
for(c in 1:10){
  thisVar<-paste("Hoki", c, "_Catch", sep=""); thisCatchData<-ncvar_get(catchNC.nc, thisVar)
  hokiCatchByCohort[c,]<-apply(thisCatchData, 2, sum, na.rm=TRUE)
}
test<-apply(hokiCatchByCohort, 2, sum, na.rm=TRUE)
plot(x=allYears, y=test)
abline(v=2015)


#at a given ts, how do they compare..?
this_tsIndex<-170; start_tsIndex<-150 # 170 is 20 years into future (2015+20)
compBiomData<-topBiomassData[,,this_tsIndex]; baseBiomData<-topBiomassData[,,start_tsIndex]
## or base is taken to be status quo
tempBase<-topBiomassData[3,,this_tsIndex]
baseBiomData <- matrix(rep(tempBase, nrow=nruns), ncol=length(testTracerIndex), nrow=nruns, byrow = TRUE)
relCompData<-(compBiomData - baseBiomData ) / compBiomData
colnames(relCompData)<-expTracers[testTracerIndex]
rownames(relCompData)<-runDesc
signif(relCompData, 2)

ntt<-length(testTracerIndex)

## plot 
shift<-seq(-0.1,0.1,length.out=nruns)
thisYmin<-min(relCompData, na.rm=TRUE); thisYmax<-max(relCompData, na.rm=TRUE)
par(mar=c(8,4,1,1))
plot(x=seq(1,ntt), y=rep(0,ntt), type="n", xlab="", ylab="Proportional change in biomass (tonnes)", xaxt="n", ylim=c(-10 , 2))
axis(at=seq(1,ntt), labels=gsub("_|N", " ", topFishTracers), side=1)
for(t in 1:ntt){
  thisX<-t+shift
  thisY<-relCompData[,t]
  points(x=thisX, y=thisY, type="h", lwd=5,col=colByRun)
}
legend(legend=runDesc, col=colByRun, lwd=3, x="left",bty="n")

timeIndex<-145:min(nts_list)
#take the base run to be that with status quo fishing and no pp change (run 3). compare the others to this
baseRunIndex<-3; baseRunBiomassTracers<-allBiomassTracers[baseRunIndex,,]
## turn zeros to NAs in base biomass
baseRunBiomassTracers<-apply(baseRunBiomassTracers, c(1,2), zero2NA)
relativeBiomass2base<-0*allBiomassTracers[,,timeIndex]
for(r in 1:nruns){
  thisData<-allBiomassTracers[r,,timeIndex]; thisRel<-thisData/baseRunBiomassTracers
  relativeBiomass2base[r,,]<-thisRel
}


#overall, seems zero catch has more impact than half PP - how do the key fisheries species look..?
topFishCodes<-c("HOK", "ORH", "SSO", "LIN", "BOE", "PFM")
topFishTracers<-paste(str_trim(groupsDF$Name[groupsDF$Code %in% topFishCodes]), "_N", sep="")
testTracerIndex<-grep(paste(topFishTracers,collapse = "|"), expTracers)
ntt<-length(topFishTracers)
plot(seq(1,ntt), type="n", ylim=c(-2.5,0.2), xaxt="n", xlab="")
axis(at=seq(1,ntt), labels=gsub("_|_N"," ", topFishTracers), side=1)

baseX<-seq(1,dim(relativeBiomass2base)[3])
for(t in 1:ntt){
  thisTracer<-paste("^",topFishTracers[t],sep=""); 
  thisData<-allBiomassTracers[,grep(thisTracer, expTracers),]; 
  thisBaseData<-baseRunBiomassTracers[grep(thisTracer, expTracers), ]
  for(r in 1:nruns){
    if(r != baseRunIndex){
      thisRel<-(thisBaseData - thisData[r,]) / thisBaseData
      thisMeanRel<-nonZeroMean(thisRel); thisMax<-max(thisRel, na.rm=TRUE); thisMin<-min(thisRel, na.rm=TRUE)
      thisX<-t * r/baseRunIndex
      points(x=thisX, y=thisMeanRel, type="h", col=colByRun[r], lwd=8)
    }
  }
}

test<-allBiomassTracers[3,grep("Hoki", expTracers),]

otherIsFishedCodes<-groupsDF$Code[groupsDF$IsFished==1 & !(groupsDF$Code %in% topFishCodes)]
topFishTracers<-paste(str_trim(groupsDF$Name[groupsDF$Code %in% otherIsFishedCodes]), "_N", sep="")
testTracerIndex<-grep(paste(topFishTracers,collapse = "|"), expTracers)
ntt<-length(topFishTracers)
plot(seq(1,ntt), type="n", ylim=c(-2.5,1), xaxt="n", xlab="")
axis(at=seq(1,ntt), labels=gsub("_|_N"," ", topFishTracers), side=1)
for(t in 1:ntt){
  thisTracer<-paste("^",topFishTracers[t],sep=""); 
  thisData<-allBiomassTracers[1:3,grep(thisTracer, expTracers),]; thisRatio<-(thisData[1,]-thisData[2])/thisData[1,]; thisMean<-mean(thisRatio, na.rm=TRUE)
  points(x=t, y=thisMean, type="h", lwd=3)
}







testTracer<-"PicoPhytopl_N"; 
# testTracer<-"Diatom_N";
# testTracer<-"Macroalgae_N"
t=grep(testTracer, expTracers)[1]
r=2
thisData<-apply(storeTracers[r,t,,,] * runVol, 3, sum, na.rm=TRUE) * mg_2_tonne * X_CN
TEST2<-thisData
r=1
thisData<-apply(storeTracers[r,t,,,] * runVol, 3, sum, na.rm=TRUE) * mg_2_tonne * X_CN
TEST1<-thisData
TEST2-TEST1
r=3
thisData<-apply(storeTracers[r,t,,,] * runVol, 3, sum, na.rm=TRUE) * mg_2_tonne * X_CN
TEST3<-thisData
TEST3-TEST1


timeIndex<-145:min(nts_list)

colByEddy<-c(myOrange, myBlue,myGreen)

testTracerIndex<-grep("Zoo_N",expTracers); zooTracers<- expTracers[testTracerIndex]
testTracerIndex<-grep("Bac", expTracers); bacteriaTracers<-expTracers[testTracerIndex]
testTracerIndex<-grep("Det_N",expTracers); detTracers<-expTracers[testTracerIndex]

thisTracers<-zooTracers
thisTracers<-bacteriaTracers
thisTracers<-detTracers

par(mfrow=c(2,2))
for(testTracer in thisTracers){
  
  testTracerIndex<-expTracers==testTracer
  testBiomass<-allBiomassTracers[,testTracerIndex,timeIndex]
  thisYmax<-max(testBiomass, na.rm=TRUE); thisYmin<-min(testBiomass, na.rm=TRUE)
  plot(testBiomass[1,], type="n", ylim=c(thisYmin,thisYmax))
  for(r in 1:2){
    thisCol<-colByEddy[r]; thisLty=r
    points(testBiomass[r,], type="l", lwd=3, col=thisCol, lty=thisLty)
  }
  legend(legend=c("Base", "Reduced PP"), col=colByEddy,lty=1:2, lwd=2, x="bottomright")
  mtext(testTracer, side=3)
}


ppcodes<-c("DF", "MA", "MB", "PL", "PS")
primaryProducers<-c(paste(groupsDF$Name[groupsDF$Code %in% ppcodes],"_N", sep=""))
testTracerIndex<-grep(paste(primaryProducers,collapse="|", sep=""), expTracers)
expTracers[testTracerIndex]
par(mfrow=c(3,2))
for(testTracer in primaryProducers){
  
  testTracerIndex<-expTracers==testTracer
  testBiomass<-allBiomassTracers[,testTracerIndex,timeIndex]
  thisYmax<-max(testBiomass, na.rm=TRUE); thisYmin<-min(testBiomass, na.rm=TRUE)
  plot(testBiomass[1,], type="n", ylim=c(thisYmin,thisYmax))
  for(r in 1:3){
    thisCol<-colByEddy[r]; thisLty=r
    points(testBiomass[r,], type="l", lwd=3, col=thisCol, lty=thisLty)
  }
  mtext(testTracer, side=3)
}
makeBlankPlot();   legend(legend=c("Base", "Reduced PP"), col=colByEddy,lty=1:2, lwd=2, x="center")

#read in trophic levels and plot change by trophic level (perhaps med, LQ, UQ)
trophicLevels<-read.csv(paste0(DIR$'Base',"ATLANTISmodels\\base\\EWEbase\\","CRAMGroupsTL.csv"))
#take out DC as not used in no-fishing runs
trophicLevels$TL[trophicLevels$Code=="DC"]<-NA

runIndex<-c(1,2); thisTitle<-"Reduced PP wrt base PP, both no fishing " ## set the first one to be the base
# runIndex<-c(3,4); thisTitle<-"Reduced PP wrt base PP, both status quo fishing " ## set the first one to be the base
# runIndex<-c(2,4); thisTitle<-"Status quo wrt base fishing, both Reduced PP " ## set the first one to be the base
# runIndex<-c(1,3); thisTitle<-"Status quo wrt no fishing"
relByTL<-rep(0,6); relByTLtime<-array(NA, dim=c(6,length(timeIndex)))
for(l in 0:5){
  thisCodes<-trophicLevels$Code[trunc(trophicLevels$TL)==l]
  thisTracerIndex<-grep(paste(groupsDF$Name[groupsDF$Code %in% thisCodes], collapse="|"), expTracers)
  thisTracers<-expTracers[thisTracerIndex]
  thisData<-allBiomassTracers[runIndex,thisTracerIndex,timeIndex]; ntt<-length(thisTracers)
  storeRel<-rep(0,ntt); storeRelTime<-array(NA, dim=c(ntt, dim(thisData)[3]))
  for(t in 1:ntt){
    thisRel<-(thisData[1,t,] - thisData[2,t,])/thisData[1,t,]
    storeRel[t]<-nonZeroMean(thisRel); storeRelTime[t,]<-thisRel
  }
  relByTL[l+1]<-median(storeRel, na.rm=TRUE)
  relByTLtime[l+1,]<-apply(storeRelTime,2,median, na.rm=TRUE)
}
toPlot<-melt(relByTLtime); colnames(toPlot)<-c("TL", "Year","Value")
toPlot$TL<-toPlot$TL-1
bp<-ggplot(data = toPlot, aes(x = Year, fill = TL, y = Value)) + 
  geom_bar(stat = 'identity')


pdf(paste(plotPath,"comparePPruns_", gsub(" |,","",thisTitle), ".pdf", sep=""), height=4, width=5)
par(mar=c(4,3.5,1,1))
bp + labs(y="Proportional decrease", x="Year", title=thisTitle) + theme_igray() +
  theme(axis.text=element_text(size=12, angle=0),axis.title=element_text(size=12)) + 
  guides(fill=guide_legend(title="Trophic Level"))  
 dev.off() 



primaryProducers<-c("DinoFlag_N", "Macroalgae_N", "MicroPB_N", "Diatom_N", "PicoPhytopl_N")
thisTracerIndex<-grep(paste(primaryProducers, collapse="|"), expTracers)
## total primary productivity change
combinedPPbiomass<-apply(allBiomassTracers[1:2,thisTracerIndex,],c(1,3), sum, na.rm=TRUE)
combinedPPratio<-(combinedPPbiomass[1,] - combinedPPbiomass[2,]) / combinedPPbiomass[1,]
plot(combinedPPratio, type="l")

toPlot2<-data.frame(melt(combinedPPratio[1:118]))
toPlot2$x<-seq(1,dim(toPlot2)[1]); colnames(toPlot2)<-c("y","x")

lp<-ggplot(data=toPlot2, aes(y=y, x=x)) + geom_line(size=1.2) 
lp + labs(y="Proportional decrease", x="Year") + theme_igray() +
  theme(axis.text=element_text(size=12, angle=0),axis.title=element_text(size=12)) + 
  guides(fill=guide_legend(title="Trophic Level"))  


# dev.off()
## within each trophic level, what actually happend- some possibly went up and some down, although the overal seems to have been
## always a reduction
## take carrion out of trophic levels
trophicLevels$TL[trophicLevels$Code=="DC"]<-NA
thisTL=5
thistrophicTracers<-paste(str_trim(groupsDF$Name[groupsDF$Code %in% trophicLevels$Code[round(trophicLevels$TL)==thisTL]]), "_N", sep="")
thisTracerIndex<-grep(paste(thistrophicTracers, collapse="|"), expTracers)
thisBiomData<-allBiomassTracers[1:2, thisTracerIndex, 1:min(nts_list)]
thisRelBiom<-0*thisBiomData[1,,]
for(t in 1:length(thistrophicTracers)){
  thisRelBiom[t,]<- (thisBiomData[1,t,] - thisBiomData[2,t, ]) / thisBiomData[1,t,]
  plot(thisRelBiom[t,], type="l"); mtext(expTracers[thisTracerIndex][t])
}


## do all trophic levels, 


colByGroup<-colorRampPalette(colors=rev(c("midnightblue", myBlue,myAqua)))(length(thisTracerIndex))

toPlot<-melt(thisRelBiom)
colnames(toPlot)<-c("Species","Year","Value")
toPlot$Species<-expTracers[thisTracerIndex][match(toPlot$Species,1:length(thisTracerIndex))]
bp<-ggplot(data = toPlot, aes(x = Year, fill = Species, y = Value)) + 
  geom_bar(stat = 'identity')

par(mar=c(4,3.5,1,1))
bp + labs(y="Proportional decrease", x="Year") + theme_igray() +
  theme(axis.text=element_text(size=12, angle=0),axis.title=element_text(size=12)) + 
  guides(fill=guide_legend(title="Species group"))   + scale_fill_manual(values=colByGroup)  






#what is the total (or median..?) pp change?
timeIndex<-seq(1,114)
thisTracers<-expTracers[thisTracerIndex]
thisData<-allBiomassTracers[1:2,thisTracerIndex,timeIndex]; ntt<-length(thisTracers)
storeRel<-rep(0,ntt); storeRelTime<-array(NA, dim=c(ntt, length(timeIndex)))
for(t in 1:ntt){
  thisRel<-(thisData[1,t,] - thisData[2,t,])/thisData[1,t,]
  storeRel[t]<-nonZeroMean(thisRel); storeRelTime[t,]<-thisRel
}
ppSumByTime<-apply(storeRelTime,2,sum, na.rm=TRUE)
ppMedByTime<-apply(storeRelTime,2, median, na.rm=TRUE)

colByPP<-colorRampPalette(colors=c(myGold,"red",myBlue,myAqua,myGreen))(5)



toPlot<-melt(storeRelTime)
colnames(toPlot)<-c("PrimaryProducer","Year","Value")
toPlot$`PrimaryProducer`<-gsub("_N|_", "", expTracers[thisTracerIndex][match(toPlot$`PrimaryProducer`, seq(1,length(primaryProducers)))])
bp<-ggplot(data = toPlot, aes(x = Year, fill = PrimaryProducer, y = Value)) + scale_fill_manual(values=colByPP)  + 
  geom_bar(stat = 'identity')

par(mar=c(4,3.5,1,1))
bp + labs(y="Proportional decrease", x="Year") + theme_igray() +
  theme(axis.text=element_text(size=12, angle=0),axis.title=element_text(size=12)) + guides(fill=guide_legend(title="Primary producer")) 


testData<-ratioPPtracers[thisRuns,testTracerIndex,timeIndex]

testMean<-apply(testData,c(1,2), mean, na.rm=TRUE); ymax<-max(testMean, na.rm=TRUE)
colByZ<-colorRampPalette(colors=c(myGold, myGreen,myAqua,myBlue,myRed))(length(testTracerIndex))
plot(testMean[,1], ylim=c(0.8,ymax*1.1), type="n", xlab="Run", ylab="Biomass relative to base")
for(z in 1:length(testTracerIndex)){
  points(testMean[,z], pch=20, col=colByZ[z],cex=1.5)
}
legend(legend=expTracers[testTracerIndex], pch=20, pt.cex=1.5, col=colByZ, x="topleft", bty="n", ncol=3)

testTracerIndex<-grep("_N", expTracers)
thisTracers<-expTracers[testTracerIndex]
ntt<-length(thisTracers)
par(lend=1, las=2, mfrow=c(1,1), mar=c(10,4,1,1))
plot(seq(1,ntt), type="n", ylim=c(-2.5,1), xaxt="n", xlab="")
axis(at=seq(1,ntt), labels=gsub("_|_N"," ", thisTracers), side=1)
for(t in 1:ntt){
  thisTracer<-paste("^",thisTracers[t],sep=""); 
  thisData<-allBiomassTracers[3:2,grep(thisTracer, expTracers),]; thisRatio<-(thisData[1,]-thisData[2])/thisData[1,]; thisMean<-mean(thisRatio, na.rm=TRUE)
  points(x=t, y=thisMean, type="h", lwd=3)
}

topFishCodes<-c("HOK", "ORH", "SSO", "LIN", "BOE", "PFM")
topFishTracers<-paste(str_trim(groupsDF$Name[groupsDF$Code %in% topFishCodes]), "_N", sep="")
testTracerIndex<-grep(paste(topFishTracers,collapse = "|"), expTracers)
ntt<-length(topFishTracers)
plot(seq(1,ntt), type="n", ylim=c(-2.5,1), xaxt="n", xlab="")
axis(at=seq(1,ntt), labels=gsub("_|_N"," ", topFishTracers), side=1)
for(t in 1:ntt){
  thisTracer<-paste("^",topFishTracers[t],sep=""); 
  thisData<-allBiomassTracers[2:3,grep(thisTracer, expTracers),]; thisRatio<-(thisData[1,]-thisData[2])/thisData[1,]; thisMean<-mean(thisRatio, na.rm=TRUE)
  points(x=t, y=thisMean, type="h", lwd=3)
}

otherIsFishedCodes<-groupsDF$Code[groupsDF$IsFished==1 & !(groupsDF$Code %in% topFishCodes)]
topFishTracers<-paste(str_trim(groupsDF$Name[groupsDF$Code %in% otherIsFishedCodes]), "_N", sep="")
testTracerIndex<-grep(paste(topFishTracers,collapse = "|"), expTracers)
ntt<-length(topFishTracers)
plot(seq(1,ntt), type="n", ylim=c(-2.5,1), xaxt="n", xlab="")
axis(at=seq(1,ntt), labels=gsub("_|_N"," ", topFishTracers), side=1)
for(t in 1:ntt){
  thisTracer<-paste("^",topFishTracers[t],sep=""); 
  thisData<-allBiomassTracers[1:3,grep(thisTracer, expTracers),]; thisRatio<-(thisData[1,]-thisData[2])/thisData[1,]; thisMean<-mean(thisRatio, na.rm=TRUE)
  points(x=t, y=thisMean, type="h", lwd=3)
}




