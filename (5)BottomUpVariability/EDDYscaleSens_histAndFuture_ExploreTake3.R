## catch histories were created in setUpScenario_ts_files_versionC1.R - they are full historic + 50 year future runs
## biol.prm files for eddy_scale are  created first

thisPath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

version<-"C"
runPath<-paste(thisPath,"eddySens\\",sep="")

groupsDF<-read.csv(paste(thisPath,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

plotPath<-paste(thisPath,"Figures\\eddySens\\",version,"",sep="")

## read in trophic levels so can plot by these
trophicLevels<-read.csv(paste(thisPath, "base\\EWEbase\\CRAMGroupsTL.csv", sep=""))

load(paste(thisPath,"eddySens\\modelTracers_Vsave",version,sep="")) ## brings in relTracers and storeTracers
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
relRunByTS<-array(NA, dim=c(nruns, ntracers, dim(storeBiomass)[3]))
baseBiomassData<-storeBiomass[baseRunIndex,,]
for(r in 1:nruns){
  thisData<-storeBiomass[r,,]
  relRunByTS[r,,]<-thisData / baseBiomassData
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

## want to check if those with reduced pp growth are different to the others
redPPIndex<-lookup_df$Eddy==2
redPPtraceres<-relRunByTS[redPPIndex,,]; fullPPtracers<-relRunByTS[!redPPIndex,,]

ratioPPtracers<-redPPtraceres/fullPPtracers

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

testTracer<-"PicoPhytopl_N"; 
testTracer<-"Diatom_N"; 

timeIndex<-150:202

colByEddy<-c(myOrange, myBlue)
ppcodes<-c("DF", "MA", "MB", "PL", "PS")
primaryProducers<-c(paste(groupsDF$Name[groupsDF$Code %in% ppcodes],"_N", sep=""))

par(mfrow=c(3,2), mar=c(4,4,1,1))
for(testTracer in primaryProducers){
  
  testTracerIndex<-expTracers==testTracer
  testBiomass<-allBiomassTracers[,testTracerIndex,timeIndex]
  thisYmax<-max(testBiomass, na.rm=TRUE); thisYmin<-min(testBiomass, na.rm=TRUE)
  plot(testBiomass[1,], type="n", ylim=c(thisYmin,thisYmax))
  # for(r in 1:nruns){
  for(r in 3:4){
    thisCol<-colByEddy[lookup_df$Eddy[r]]; thisLty=lookup_df$Eddy[r]
    points(testBiomass[r,], type="l", lwd=3, col=thisCol, lty=thisLty)
  }
  mtext(testTracer, side=3)
}
makeBlankPlot()
legend(legend=c( "Base", "Reduced PP"), col=colByEddy, lty=c(1,2), x="center", bty="n", ce==1.5, lwd=2, seg.len=3)


testRatioB<-testBiomass[redPPIndex,]/ testBiomass[!redPPIndex,]

######
thisRuns<-c(3:8,11:12) ## 
testTracerIndex<-grep("_N", expTracers); testTracers<-expTracers[testTracerIndex]

## take Z's
testTracerIndex<-grep("Zoo_N",expTracers); expTracers[testTracerIndex]

testData<-ratioPPtracers[thisRuns,testTracerIndex,timeIndex]

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


testData<-ratioPPtracers[thisRuns,testTracerIndex,timeIndex]

testMean<-apply(testData,c(1,2), mean, na.rm=TRUE); ymax<-max(testMean, na.rm=TRUE)
colByZ<-colorRampPalette(colors=c(myGold, myGreen,myAqua,myBlue,myRed))(length(testTracerIndex))
plot(testMean[,1], ylim=c(0.8,ymax*1.1), type="n", xlab="Run", ylab="Biomass relative to base")
for(z in 1:length(testTracerIndex)){
  points(testMean[,z], pch=20, col=colByZ[z],cex=1.5)
}
legend(legend=expTracers[testTracerIndex], pch=20, pt.cex=1.5, col=colByZ, x="topleft", bty="n", ncol=3)

#########################
## rel change by polygon
dynBoxes<-seq(2,25); ndboxes<-length(dynBoxes)
relChangeByBoxTracer<-array(NA, dim=c(nruns, ntracers,nboxes)); relChangeByBox<-array(NA, dim=c(nruns, nboxes))
baseByBoxTracer<-apply(storeTracers[baseRunIndex,,,,],c(1,3), nonZeroMean)
baseByBox<-apply(storeTracers[baseRunIndex,,,,],3, nonZeroMean)
## set v. smalls to NA - otherwise get very big relatives when only small amounts
setSmall2zero<-function(x,z=1e-6){
  y<-x
  if(!is.na(y)){
    if(x<z){y<-NA}
  }
  return(y)
}
baseByBoxTracer<-apply(baseByBoxTracer,c(1,2), setSmall2zero)
# baseByBox<-apply()
for(r in 1:nruns){
  thisByBox<-apply(storeTracers[r,,,,], c(1,3), nonZeroMean)
  relChangeByBoxTracer[r,,]<-(baseByBoxTracer - thisByBox)/baseByBoxTracer
  thisByBox<-apply(storeTracers[r,,,,], 3, nonZeroMean)
  relChangeByBox[r,]<-(baseByBox - thisByBox)/baseByBox
}
meanByBox<-apply(relChangeByBoxTracer,3, mean, na.rm=TRUE)
maxByBox<-apply(relChangeByBoxTracer,3, max, na.rm=TRUE)


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

thisMax<-max(meanByBox,na.rm=TRUE);  thisMin<-min(meanByBox, na.rm=TRUE)
boxColors<-unlist(lapply(meanByBox, getColor, thisMax=thisMax))

thisMax<-max(maxByBox,na.rm=TRUE);  thisMin<-min(maxByBox, na.rm=TRUE)
boxColors<-unlist(lapply(maxByBox, getColor, thisMax=thisMax))



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
baseRunIndex<-3
comNames<-str_trim(groupsDF$Name[groupsDF$Code %in% c("BOE", "HOK", "LIN", "ORH", "PFM", "SSO")], side="both"); ncoms<-length(comNames)
comTracers<-paste(comNames,"_N", sep="")
comTracerIndex<-grep(paste(comTracers,collapse="|"), expTracers)
storeComBiom<-array(NA, dim=c(nruns, ncoms, dim(storeTracers)[5]))
for(c in 1:ncoms){
  thisData<-storeTracers[,comTracerIndex[c],,,]
  thisBiomasses<-array(NA, dim=dim(storeTracers)[c(1,5)])
  thisBaseData<-storeTracers[baseRunIndex,comTracerIndex[c],,,]; thisBaseVol<-storeTracers[baseRunIndex,grep("vol", expTracers),,,]
  thisBaseBiomass<-apply(thisBaseData * thisBaseVol, 3, sum) *mg_2_tonne * X_CN
  for(r in 1:nruns){
    thisRunData<-storeTracers[r,comTracerIndex[c],,,]; thisRunVol<-storeTracers[r,grep("vol", expTracers),,,]
    thisB<-apply(thisRunData * thisRunVol, 3, sum) *mg_2_tonne * X_CN
    storeComBiom[r,c,]<-thisB/thisBaseBiomass
  }
}

c=1
thisData<-storeComBiom[,c,]

test<-apply(storeComBiom, c(1,2), nonZeroMean)

par(mfrow=c(3,2))
for(c in 1:ncoms){
  plot(test[,c],pch=20, ylim=c(0,1.3)); mtext(expTracers[comTracerIndex][c], side=3, adj=0)
}
