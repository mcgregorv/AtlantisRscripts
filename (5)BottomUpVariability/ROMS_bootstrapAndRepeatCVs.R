# tracers were read in and written out in exploringROMS_bootstrapCVs.R
basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
plotPath<-paste(DIR$'Figures',"ROMS\\BootVsRepeat_",sep="")
## repeat ROMS is Version D; bootstrap ROMS is Version B
VersionRepeat<-"D"; versionBoot<-"B"
outPath<-paste(basePath,"base\\ouputROMS",VersionRepeat,"\\",sep="")
dataOutPath<-paste(basePath,"BootstrapROMS\\modelOut",VersionRepeat,"_", sep="")
load(paste(dataOutPath,"modelTracers",sep="")); ##to bring "storeTracers", "storeTracersByCell", "storeTemperature", "store_nts"
## storeTracers and storeTracersByCell are already converted to tonnes
repeatTracers<-storeTracers; repeat_ntsList<-store_nts
repeatTemperature<-storeTemperature
repeatTracersByCell<-storeTracersByCell
##
dataOutPath<-paste(basePath,"BootstrapROMS\\modelOut",versionBoot,"_", sep="")
load(paste(dataOutPath,"modelTracers",sep="")); ##to bring "storeTracers", "storeTracersByCell", "storeTemperature", "store_nts"
## storeTracers and storeTracersByCell are already converted to tonnes
bootTracers<-storeTracers; boot_ntsList<-store_nts
bootTemperature<-storeTemperature
bootTracersByCell<-storeTracersByCell

nBootRuns<-dim(bootTracers)[1]; nRepeatRuns<-dim(repeatTracers)[1]

#####################
## join them together
nAllRuns<-dim(bootTracers)[1] + dim(repeatTracers)[1]
storeTracersByCell<-array(NA, dim=c(nAllRuns, dim(bootTracersByCell)[-1]))
storeTracersByCell[1:nBootRuns,,,,]<-bootTracersByCell
storeTracersByCell[(nBootRuns+1):nAllRuns,,,,]<-repeatTracersByCell

storeTracers<-array(NA, dim=c(nAllRuns, dim(bootTracers)[-1]))
storeTracers[1:nBootRuns,,]<-bootTracers
storeTracers[(nBootRuns+1):nAllRuns,,]<-repeatTracers

nts<-dim(repeatTracers)[3]
#############################################################################################
## cvs by tracers 
calcCV<-function(x){
  thisMean<-mean(x,na.rm=TRUE); thisVar<-var(x,na.rm=TRUE)
  thisCV<-sqrt(thisVar)/thisMean
  return(thisCV)
}
getColor<-function(x,thisMax){
  thisCol<-"white"
  if(!is.na(x) & x>0){
    if(x<thisMin){x<-thisMin}
    y<-round((x-thisMin)/(thisMax-thisMin),2)*100+1
    thisCol<-thisColRamp[y]
  }
  return(thisCol)
}
plotGrid<-function(x,y){
  thisX<-c(x-0.5,x-0.5,x+0.5,x+0.5); thisY<-c(y-0.5,y+0.5,y+0.5,y-0.5)
  thisCol<-plotColour[x,y]
  if(length(thisCol)>0){
    polygon(x=thisX,y=thisY,col=thisCol,border=NA)
  }
  return(NULL)
}

thisColRamp<-colorRampPalette(colors=c(myLightAqua,myBlue,"midnightblue"))(101)

storeCVsbyTracer<-array(NA,dim=c(ntracers,nts))
for(t in 1:ntracers){
  thisData<-storeTracersByCell[,t,,,]
  #add up tracers for each run/timestep
  sumTracers<-apply(thisData,c(1,4), sum, na.rm=TRUE)
  thisCVwrtTime<-apply(sumTracers,2,calcCV)
  storeCVsbyTracer[t,]<-thisCVwrtTime
}
colByDepth<-colorRampPalette(colors=c(myAqua,myBlue,'midnightblue'))(nlayers)
thisMax<-myRounding(max(storeCVsbyTracer,na.rm=TRUE),0.05, direction="up"); thisMin<-myRounding(min(storeCVsbyTracer, na.rm=TRUE), 0.05, direction="down")

plotColour<-apply(t(storeCVsbyTracer),c(1,2),getColor, thisMax=thisMax)
plotData<-t(storeCVsbyTracer)
tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotData)[1]),dim(plotData)[2]),"y"=sort(rep(seq(1,dim(plotData)[2]),dim(plotData)[1]))))

pdf(paste(plotPath,"CVbyTracerTime.pdf",sep=""),height=8,width=8)
par(mar=c(4,9,1.5,1),las=1)
plot(x=seq(1,dim(plotData)[1]),y=rep(dim(plotData)[2],dim(plotData)[1]),type="n",ylab="",xlab="Timestep (years)",xaxt="n",yaxt="n",ylim=c(0.5,dim(plotData)[2]+0.5))
axis(at=seq(1,dim(plotData)[2]),labels = gsub("_N|_"," ",Ntracers),side=2,las=1)
axis(at=seq(1,dim(plotData)[1]),labels=seq(1,dim(plotData)[1]),side=1,las=1)
temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
dev.off()


###########################################################
dynBoxIndex<-seq(2,25)
storeCVsbyDepth<-array(NA,dim=c(nlayers,nts))
layerIndex<-c(seq((nlayers-1),1),nlayers)
for(l in 1:nlayers){
  thisLayer<-layerIndex[l]
  thisData<-storeTracersByCell[,,thisLayer,dynBoxIndex,]
  #add up tracers for each run/timestep
  sumTracers<-apply(thisData,c(1,4), sum, na.rm=TRUE)
  thisCVwrtTime<-apply(sumTracers,2,calcCV)
  storeCVsbyDepth[l,]<-thisCVwrtTime
}

thisColRamp<-colorRampPalette(colors=c(myLightAqua,myBlue,"midnightblue"))(101)
thisMax<-max(storeCVsbyDepth,na.rm=TRUE); thisMin<-myRounding(min(storeCVsbyDepth, na.rm=TRUE),fraction=0.05, direction="down")
plotColour<-apply(t(storeCVsbyDepth),c(1,2),getColor, thisMax=thisMax)
plotData<-t(storeCVsbyDepth)
tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotData)[1]),dim(plotData)[2]),"y"=sort(rep(seq(1,dim(plotData)[2]),dim(plotData)[1]))))

depthLabels<-c("0-200 m","200-400 m","400-600 m","600-800 m","800-1300+ m","Sediment")
pdf(paste(plotPath,"CVbyDepthTime.pdf",sep=""),height=4,width=7)
par(mar=c(4,7,1.5,1),las=1)
plot(x=seq(1,dim(plotData)[1]),y=rep(dim(plotData)[2],dim(plotData)[1]),type="n",ylab="",xlab="Timestep (years)",xaxt="n",yaxt="n",ylim=c(0.5,dim(plotData)[2]+0.5))
axis(at=seq(1,dim(plotData)[2]),labels = depthLabels,side=2,las=1)
axis(at=seq(1,dim(plotData)[1]),labels=seq(1,dim(plotData)[1]),side=1,las=1)
temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
dev.off()

xx<-pretty(seq(thisMin,thisMax,length.out = 5)); 
legendText<-xx[xx<thisMax]; legendCols<-unlist(lapply(legendText,getColor, thisMax=thisMax))
pdf(paste(plotPath,"CVbyDepthTimeLEGEND.pdf",sep=""),height=3,width=2.5)
par(mar=c(0,0,0,0))
makeBlankPlot(); legend(legend=as.character(formatC(legendText,format="f",digits=4)), col=legendCols, x="center",bty="n",pch=15, title="CV",pt.cex=1.5)
dev.off()

########################################################
## by box
nboxes<-30
dynBoxIndex<-seq(2,25); 
storeCVsbyBoxTime<-array(NA,dim=c(nboxes,nts)); storeCVsbyBox<-rep(NA, nboxes)
for(b in 1:nboxes){
  thisData<-storeTracersByCell[,,,b,]
  #add up tracers for each run/timestep
  sumTracers<-apply(thisData,c(1,4), sum, na.rm=TRUE)
  thisCVwrtTime<-apply(sumTracers,2,calcCV)
  storeCVsbyBoxTime[b,]<-thisCVwrtTime
  #just by box
  thisSum<-apply(thisData,1,sum, na.rm=TRUE); thisCVwrtBox<-calcCV(thisSum)
  storeCVsbyBox[b]<-thisCVwrtBox
}
#for spatial, plot max cv by box
thisMax<-myRounding(max(storeCVsbyBox, na.rm=TRUE), 0.05, direction="up"); thisMin<-myRounding(min(storeCVsbyBox, na.rm=TRUE), 0.05, direction="down")
boxColors<-unlist(lapply(storeCVsbyBox, getColor, thisMax=thisMax))

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

pdf(paste(plotPath,"CVsBoxMap.pdf", sep=""), height=3,width=4)
plot(shape)
map('nzHires',add=TRUE,col="black",lwd=2)
map.axes()
for(b in 1:dim(labeldf)[1]){
  thisCol<-"white"
  if(b %in% dynBoxIndex){thisCol<-boxColors[b]}
  polygon(sdata$shp$shp[[b]]$points,col=thisCol)
}
dev.off()


xx<-pretty(seq(thisMin,thisMax,length.out = 5)); 
legendText<-xx[xx<thisMax]; legendCols<-unlist(lapply(legendText,getColor, thisMax=thisMax))
pdf(paste(plotPath,"ROMSCVboxLEGEND.pdf",sep=""),height=3,width=2.5)
par(mar=c(0,0,0,0))
makeBlankPlot(); legend(legend=as.character(formatC(legendText,format="f",digits=4)), col=legendCols, x="center",bty="n",pch=15, title="CV",pt.cex=1.5)
dev.off()

##################################
## by trophic level
#read in trophic levels
trophicLevelDF<-read.csv(paste(basePath,"inputs\\biol_prm\\CRAM_trophicLevels_isotopes.csv",sep=""))
trophicLevelDF$roundedTL<-myRounding(trophicLevelDF$Isotope,0.5, direction="down")
trophicLevelDF$roundedTL[is.na(trophicLevelDF$Isotope)]<-myRounding(trophicLevelDF$TrophicLevel2[is.na(trophicLevelDF$Isotope)],0.5, direction="down")

trophicLevels<-sort(unique(trophicLevelDF$roundedTL)); ntls<-length(trophicLevels)

cvByTrophicTime<-array(NA,dim=c(ntls,nts))
for(t in 1:ntls){
  thisCodes<-trophicLevelDF$Code[trophicLevelDF$roundedTL==trophicLevels[t]]
  thisNames<-unlist(lapply(groupsDF$Name[groupsDF$Code %in% thisCodes],str_trim,side="both"))
  thisTracers<-paste(thisNames,"_N", sep=""); tracerIndex<-grep(paste(thisTracers,collapse="|"),Ntracers)
  thisData<-storeTracersByCell[,tracerIndex,,,]
  ##
  thisData<-storeTracers[,tracerIndex,]
  timeDim<-length(dim(thisData))
  xx<-apply(thisData, c(1,timeDim), sum, na.rm=TRUE)
  cvByTime<-apply(xx, 2, calcCV)
  cvByTrophicTime[t,]<-cvByTime
  
  
  ## sum by run and time
  timeDim<-length(dim(thisData))
  xx<-apply(thisData, c(1,timeDim), sum, na.rm=TRUE)
  cvByTime<-apply(xx, 2, calcCV)
  cvByTrophicTime[t,]<-cvByTime
}
thisMax<-myRounding(max(cvByTrophicTime,na.rm=TRUE), 0.05, direction="up"); thisMin<-myRounding(min(cvByTrophicTime, na.rm=TRUE), 0.05, direction="down")
plotColour<-apply(t(cvByTrophicTime),c(1,2),getColor, thisMax=thisMax)
plotData<-t(cvByTrophicTime)
tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotData)[1]),dim(plotData)[2]),"y"=sort(rep(seq(1,dim(plotData)[2]),dim(plotData)[1]))))

pdf(paste(plotPath,"CVbyTrphicTime.pdf",sep=""),height=4,width=6.8)
par(mar=c(4,5,1.5,1),las=1)
plot(x=seq(1,dim(plotData)[1]),y=rep(dim(plotData)[2],dim(plotData)[1]),type="n",ylab="",xlab="Timestep (years)",xaxt="n",yaxt="n",ylim=c(0.5,dim(plotData)[2]+0.5))
axis(at=seq(1,dim(plotData)[2]),labels = trophicLevels,side=2,las=1)
axis(at=seq(1,dim(plotData)[1]),labels=seq(1,nts),side=1,las=1)
temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
par(las=0)
mtext("Trophic level",side=2,adj=0.5,line=3)
dev.off()



xx<-pretty(seq(thisMin,thisMax,length.out = 5)); 
legendText<-xx[xx<thisMax]; legendCols<-unlist(lapply(legendText,getColor, thisMax=thisMax))
pdf(paste(plotPath,"ROMSCVTrophicLevelLEGEND.pdf",sep=""),height=3,width=2.5)
par(mar=c(0,0,0,0))
makeBlankPlot(); legend(legend=as.character(formatC(legendText,format="f",digits=4)), col=legendCols, x="center",bty="n",pch=15, title="CV",pt.cex=1.5)
dev.off()







