## catch histories were created in setUpScenario_ts_files_versionC1.R - they are full historic + 50 year future runs
## biol.prm files for eddy_scale are  created first

thisPath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

version<-"B"
runPath<-paste(thisPath,"eddySens\\",sep="")

groupsDF<-read.csv(paste(thisPath,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

plotPath<-paste(thisPath,"Figures\\eddySens\\",version,"",sep="")

## read in trophic levels so can plot by these
trophicLevels<-read.csv(paste(thisPath, "base\\EWEbase\\CRAMGroupsTL.csv", sep=""))

# eddysens<-c(0.025,0.05,0.075,0.1,0.25,0.5,0.75,1)

# ## set up DF with runs specified
# scenarios<-c("All0catch", "All50catch","All100catch", "All150catch", "Hoki50catch", "Hoki150catch"); nscenarios<-length(scenarios)
# eddysens<-c(0.01,0.05,0.1,0.5,1,5,10,50,100); neddys<-length(eddysens); eddyIndex<-seq(1,neddys)
# nruns<-nscenarios * neddys
# 
# lookup_df<-data.frame(array(NA, dim=c(nruns,3)))
# colnames(lookup_df)<-c("Run", "Scenario", "Eddy")
# lookup_df$Run<-seq(1,nruns)
# lookup_df$Scenario<-rep( sort(rep(scenarios, neddys)))
# lookup_df$Eddy<-rep(eddyIndex, (nscenarios))


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
relRunByTS<-array(NA, dim=c(nruns, ntracers, dim(relTracers)[5]))
baseSumData<-apply(storeTracers[nruns,,,,], c(1,4), nonZeroMean)
for(r in 1:nruns){
  thisData<-storeTracers[r,,,,]
  thisSumData<-apply(thisData,c(1,4), nonZeroMean)
  relRunByTS[r,,]<-thisSumData/baseSumData
}
cleanNAs<-function(x){
  y<-x
  if(length(x)==0){y<-NA}
  else if(is.na(x)){y<-NA}
  else if(x==""){y<-NA}
  else if(x=="NaN"){y<-NA}
  else if(x=="Inf"){y<-NA}
  return(y)
}
relRunByTS<-apply(relRunByTS,seq(1,3), cleanNAs)

## want array of relative biomasses just by group and time, not by cell, and just for run 1
par(mfrow=c(3,1))
for(r in 1:nruns){
  runTracers<-apply(storeTracers[r,,,,], c(1,4), nonZeroMean)
  baseTracers<-apply(storeTracers[nruns,,,,], c(1,4), nonZeroMean)
  runRel<-runTracers/baseTracers
  
  timeIndex<-seq(35,100) #this is just as far as they have run so far
  
  #check out how much the primary producers changed - perhaps an average timeseries
  ppcodes<-c("DF", "MA", "MB", "PL", "PS")
  primaryProducers<-c(paste(groupsDF$Name[groupsDF$Code %in% ppcodes],"_N", sep=""), paste(groupsDF$Name[groupsDF$Code %in% ppcodes],"_S", sep=""))
  testTracerIndex<-grep(paste(primaryProducers,collapse="|", sep=""), expTracers)
  expTracers[testTracerIndex]
  
  testData<-runRel[testTracerIndex,timeIndex]
  
  testMean<-apply(testData,1,mean, na.rm=TRUE)
  ppRange<-c(min(testMean), max(testMean))
  
  rangeByTL<-data.frame(matrix(NA, ncol=3, nrow=6)); colnames(rangeByTL)<-c("TL","min", "max"); rangeByTL$TL<-seq(0,5)
  rangeByTL[1,]<-c("PP",ppRange)
  for(l in 1:5){
    levelNames<-trophicLevels$Name[trunc(trophicLevels$TL)==l]
    #take out any PPs so don't have them twice
    # levelNames<-levelNames[!(levelNames %in% gsub("_N|_S","",c(primaryProducers)))]
    testTracerIndex<-grep(paste(levelNames,collapse="|", sep=""), expTracers)
    testNames<-expTracers[testTracerIndex]
    testTracerIndex<-testTracerIndex[grep("Cover", expTracers[testTracerIndex], invert=TRUE)]
    expTracers[testTracerIndex]
  
    testData<-runRel[testTracerIndex,timeIndex]
    
    testMean<-apply(testData,1,mean, na.rm=TRUE)
    
    # testMean<-apply(testData,c(1,4), mean, na.rm=TRUE)
    thisRange<-c(min(testMean), max(testMean))
    rangeByTL[l+1,c("min", "max")]<-thisRange
  }
  pdf(paste(plotPath,"meanRelBiomByTL_run",r,".pdf", sep=""), height=4,width=5)
  plot(x=seq(0,5), y=rep(0,6), ylim=c(0,max(as.double(rangeByTL$max))), type="n", xaxt="n", xlab="", ylab="Biomass relative to base model")
  axis(at=seq(0,5),labels=c("PP",seq(1,5)), side=1)
  for(i in 1:dim(rangeByTL)[1]){
    y1<-rangeByTL$min[i]; y2<-rangeByTL$max[i]
    polygon(x=c(i-1.1,i-0.9,i-0.9,i-1.1), y=c(y1, y1, y2, y2), col=myBlue, border=NA)
  }
  mtext(r,side=3,adj=0)
  dev.off()
}

######
timeIndex<-35:135
thisRuns<-1:8
testTracerIndex<-grep("_N", expTracers); testTracers<-expTracers[testTracerIndex]

## take Z's
testTracerIndex<-grep("Zoo_N",expTracers); expTracers[testTracerIndex]

testData<-relRunByTS[thisRuns,testTracerIndex,timeIndex]

testMean<-apply(testData,c(1,2), mean, na.rm=TRUE); ymax<-max(testMean, na.rm=TRUE)
colByZ<-c(myGold, myGreen,myAqua,myBlue)
plot(testMean[,1], ylim=c(0,ymax), type="n", xlab="Run", ylab="Biomass relative to base")
for(z in 1:length(testTracerIndex)){
  points(testMean[,z], pch=20, col=colByZ[z],cex=1.5)
}
legend(legend=expTracers[testTracerIndex], pch=20, pt.cex=1.5, col=colByZ, x="bottomright", bty="n")


ppcodes<-c("DF", "MA", "MB", "PL", "PS")
primaryProducers<-c(paste(groupsDF$Name[groupsDF$Code %in% ppcodes],"_N", sep=""), paste(groupsDF$Name[groupsDF$Code %in% ppcodes],"_S", sep=""))
testTracerIndex<-grep(paste(primaryProducers,collapse="|", sep=""), expTracers)
expTracers[testTracerIndex]


testData<-relRunByTS[thisRuns,testTracerIndex,timeIndex]

testMean<-apply(testData,c(1,2), mean, na.rm=TRUE); ymax<-max(testMean, na.rm=TRUE)
colByZ<-colorRampPalette(colors=c(myGold, myGreen,myAqua,myBlue,myRed))(length(testTracerIndex))
plot(testMean[,1], ylim=c(0,ymax*1.2), type="n", xlab="Run", ylab="Biomass relative to base")
for(z in 1:length(testTracerIndex)){
  points(testMean[,z], pch=20, col=colByZ[z],cex=1.5)
}
legend(legend=expTracers[testTracerIndex], pch=20, pt.cex=1.5, col=colByZ, x="topleft", bty="n")

#########################
## rel change by polygon
dynBoxes<-seq(2,25); ndboxes<-length(dynBoxes)
relChangeByBoxTracer<-array(NA, dim=c(nruns, ntracers,nboxes)); relChangeByBox<-array(NA, dim=c(nruns, nboxes))
baseByBoxTracer<-apply(storeTracers[nruns,,,,],c(1,3), nonZeroMean)
baseByBox<-apply(storeTracers[nruns,,,,],3, nonZeroMean)
## set v. smalls to NA - otherwise get very big relatives when only small amounts
setSmall2zero<-function(x,z=1e-6){
  y<-x
  if(!is.na(y)){
    if(x<z){y<-NA}
  }
  return(y)
}
baseByBoxTracer<-apply(baseByBoxTracer,c(1,2), setSmall2zero)
baseByBox<-apply()
for(r in 1:nruns){
  thisByBox<-apply(storeTracers[r,,,,], c(1,3), nonZeroMean)
  relChangeByBoxTracer[r,,]<-thisByBox/baseByBoxTracer
  thisByBox<-apply(storeTracers[r,,,,], 3, nonZeroMean)
  relChangeByBox[r,]<-thisByBox/baseByBox
}
meanByBox<-apply(relChangeByBoxTracer,3, mean, na.rm=TRUE)

#read in shape file
# pdf(paste(plotPath,"test.pdf", sep=""), height=10)
shapeFile<-paste(DIR$'Base',"ATLANTISmodels\\inputs\\bgm\\CHAT30_LL",sep="")
sdata<-read.shapefile(shapeFile)
shape<-formatShape(shapeFile=shapeFile)
ns<-length(shape)
SpDF <- SpatialPolygonsDataFrame(shape,data.frame( z=1:ns, row.names=paste("P",seq(1,(ns)),sep="")))
labels<-seq(1,(ns))
plot(shape)
LABELpOS<-polygonsLabel(shape, labels = labels, cex=.1,doPlot=FALSE)
labeldf<-data.frame(cbind("x"=LABELpOS[1:ns],"y"=LABELpOS[(ns+1):(2*ns)]))
# dev.off()

getColor<-function(x,thisMax){
  thisCol<-"white"
  if(!is.na(x) & x>0){
    if(x<thisMin){x<-thisMin}
    if(x>thisMax){x<-thisMax}
    y<-round((x-thisMin)/(thisMax-thisMin),2)*100+1
    thisCol<-thisColRamp[y]
  }
  return(thisCol)
}
thisColRamp<-colorRampPalette(colors=c(myLightAqua,myBlue,"midnightblue"))(101)

thisMax<-max(meanByBox,na.rm=TRUE); 
thisMax<-1; thisMin<-0.75
boxColors<-unlist(lapply(meanByBox, getColor, thisMax=thisMax))


dynBoxes<-seq(2,25)
# pdf(paste(plotPath,"pprelResponse.pdf", sep=""), height=3,width=4)
plot(shape)
map('nzHires',add=TRUE,col="black",lwd=2)
map.axes()
for(b in 1:dim(labeldf)[1]){
  thisCol<-"white"
  if(b %in% dynBoxes){thisCol<-boxColors[b]}
  polygon(sdata$shp$shp[[b]]$points,col=thisCol)
}
# dev.off()

# 
# xx<-pretty(seq(thisMin,thisMax,length.out = 5)); 
# legendText<-xx[xx<thisMax]; legendCols<-unlist(lapply(legendText,getColor, thisMax=thisMax))
# pdf(paste(plotPath,"ROMSCVboxLEGEND.pdf",sep=""),height=3,width=2.5)
# par(mar=c(0,0,0,0))
# makeBlankPlot(); legend(legend=as.character(formatC(legendText,format="f",digits=4)), col=legendCols, x="center",bty="n",pch=15, title="CV",pt.cex=1.5)
# dev.off()


## which commercial species are affected..?
#check out the top 5 first
comNames<-str_trim(groupsDF$Name[groupsDF$Code %in% c("BOE", "HOK", "LIN", "ORH", "PFM", "SSO")], side="both"); ncoms<-length(comNames)
comTracers<-paste(comNames,"_N", sep="")
comTracerIndex<-grep(paste(comTracers,collapse="|"), expTracers)
storeComBiom<-array(NA, dim=c(nruns, ncoms, dim(storeTracers)[5]))
for(c in 1:ncoms){
  thisData<-storeTracers[,comTracerIndex[c],,,]
  thisBiomasses<-array(NA, dim=dim(storeTracers)[c(1,5)])
  thisBaseData<-storeTracers[nruns,comTracerIndex[c],,,]; thisBaseVol<-storeTracers[nruns,grep("vol", expTracers),,,]
  thisBaseBiomass<-apply(thisBaseData * thisBaseVol, 3, sum) *mg_2_tonne * X_CN
  for(r in 1:nruns){
    thisRunData<-storeTracers[r,comTracerIndex[c],,,]; thisRunVol<-storeTracers[r,grep("vol", expTracers),,,]
    thisB<-apply(thisRunData * thisRunVol, 3, sum) *mg_2_tonne * X_CN
    storeComBiom[r,c,]<-thisB/thisBaseBiomass
  }
}

c=1
thisData<-storeComBiom[,c,]

test<-apply(storeComBiom, 2, nonZeroMean)

par(mfrow=c(3,2))
for(c in 1:ncoms){
  plot(test[,c],pch=20, ylim=c(0,1.3)); mtext(expTracers[comTracerIndex][c], side=3, adj=0)
}
