#the ROMS variables we have are for XX years
#this script do some initial plots.
library(ggmap)
library(rgdal)
library(gtable)
library(maps)
library(mapdata)
library(rgeos)
source(paste(DIR$'General functions',"\\get_first_number.R",sep=""))
source(paste(DIR$'General functions',"formatShape.R",sep=""))

basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

outPath<-paste(basePath,"base\\ROMS1yearRepeat\\",sep="")

outPath<-paste(basePath,"base\\ouputROMSB\\",sep="")

nruns<-50; nts<-51; burnin<-35

nl<-6; nb<-30

plotGrid<-function(x,y){
  thisX<-c(x-0.5,x-0.5,x+0.5,x+0.5); thisY<-c(y-0.5,y+0.5,y+0.5,y-0.5)
  thisCol<-plotColour[x,y]
  if(length(thisCol)>0){
    polygon(x=thisX,y=thisY,col=thisCol,border=NA)
  }
  return(NULL)
}

#get _N tracers
allTracers<-names(ThisNC.nc$var)
x<-grep("_N",allTracers); temp<-allTracers[x]
y<-grep("Nums",temp, invert = TRUE); Ntracers<-temp[y]; ntracers<-length(Ntracers)

#set up array to store outputs
storeTracers<-array(NA, dim=c(nruns, ntracers,nts)); storeTracersByCell<-array(NA, dim=c(nruns, ntracers,nl,nb,nts))
for(r in 1:nruns){
  cat(r," -- ")
  # thisOutPath<-paste(outPath,"outputBASE50yr",r,"\\",sep="")
  thisOutPath<-paste(outPath, "outputROMSBootstrapB",r,"\\",sep="")
  thisOutFile<-paste(thisOutPath,"output.nc",sep="")
  ThisNC.nc<-nc_open(thisOutFile)
  thisVol<-ncvar_get(ThisNC.nc,"volume")
  this_nts<-dim(thisVol)[3]
  for(t in 1:ntracers){
    thisTracer<-Ntracers[t]; thisData<-ncvar_get(ThisNC.nc, thisTracer)
    if(length(dim(thisData))==3){
      xx<-apply(thisData*thisVol,3,sum)*mg_2_tonne*X_CN
      yy<-thisData*thisVol*mg_2_tonne*X_CN
      storeTracersByCell[r,t,,,1:this_nts]<-yy
    } else{
      xx<-apply(thisData*thisVol[nl,,],2,sum)*mg_2_tonne*X_CN
      yy<-thisData*thisVol[nl,,]*mg_2_tonne*X_CN
      storeTracersByCell[r,t,nl,,1:this_nts]<-yy
    }
    storeTracers[r,t,1:this_nts]<-xx
  }
}
storeTemperature<-array(NA,dim=c(nruns,nl,nb,nts)); thisTracer<-"Temp"
store_nts<-rep(NA,nruns)
for(r in 1:nruns){
  cat(r," -- ")
  # thisOutPath<-paste(outPath,"outputBASE50yr",r,"\\",sep="")
  thisOutPath<-paste(outPath, "outputROMSBootstrapB",r,"\\",sep="")
  thisOutFile<-paste(thisOutPath,"output.nc",sep="")
  ThisNC.nc<-nc_open(thisOutFile)
  thisVol<-ncvar_get(ThisNC.nc,"volume")
  thisData<-ncvar_get(ThisNC.nc, thisTracer); this_nts<-dim(thisData)[3]
  storeTemperature[r,,,1:this_nts]<-thisData
  store_nts[r]<-this_nts
}

#############################################################################################

calcCI<-function(x,i){
  y<-sort(x[!is.na(x)])
  z<-max(round(length(y)*i), 1);
  thisOut<-y[z]
  return(thisOut)
}

calcCV<-function(x){
  thisMean<-mean(x,na.rm=TRUE); thisVar<-var(x,na.rm=TRUE)
  thisCV<-sqrt(thisVar)/thisMean
  return(thisCV)
}
storeCVsByTracer<-array(NA,dim=c(ntracers,nts))
for(t in 1:ntracers){
  thisData<-storeTracersByCell[,t,,,];
  tempCVs<-apply(thisData,c(2,3,4),calcCV); meanCV<-nonZeroMean(tempCVs)
  meanCV<-apply(tempCVs,3,nonZeroMean)
  storeCVsByTracer[t,]<-meanCV
}

######################
cvByTracer<-apply(storeCVsByTracer,1,nonZeroMean)

storeCVsByBox<-array(NA,dim=c(nboxes,nts)); dynBoxes<-seq(2,25)
for(b in 1:nboxes){
  if(b %in% dynBoxes){
    thisData<-storeTracersByCell[,,,b,];
    tempCVs<-apply(thisData,c(2,3,4),calcCV); meanCV<-nonZeroMean(tempCVs)
    thisCVs<-apply(tempCVs,3,nonZeroMean)
    storeCVsByBox[b,]<-thisCVs
  }
}
cvByBox<-apply(storeCVsByBox,1,nonZeroMean)

storeCVsbyDepth<-array(NA,dim=c(nlayers,nts))
layerIndex<-c(seq((nlayers-1),1),nlayers)
for(l in 1:nlayers){
  thisLayer<-layerIndex[l]
  # thisData<-apply(storeTracersByCell[,,thisLayer,,],c(1,4),sum,na.rm=TRUE)
  thisData<-storeTracersByCell[,,thisLayer,dynBoxIndex,];
  tempCVs<-apply(thisData,c(2,3,4),calcCV); meanCV<-nonZeroMean(tempCVs)
  thisCVs<-apply(tempCVs,3,nonZeroMean)
  storeCVsbyDepth[l,]<-thisCVs
}
cvByDepth<-apply(storeCVsbyDepth,1,nonZeroMean)
plot_nts<-49
timestepsBinned<-myRounding(seq(1,plot_nts)/10,fraction=0.2,direction="down")*10; nby<-length(unique(timestepsBinned))

cvByDepthTimeBin<-array(NA,dim=c(nlayers,nby))
for(l in 1:nlayers){
  xx<-tapply(storeCVsbyDepth[l,1:plot_nts],timestepsBinned,mean,na.rm=TRUE)
  cvByDepthTimeBin[l,]<-xx
}
thisColRamp<-colorRampPalette(colors=c(myLightAqua,myBlue,"midnightblue"))(101)
thisMax<-max(cvByDepthTimeBin,na.rm=TRUE)
plotColour<-apply(t(cvByDepthTimeBin),c(1,2),getCol)
plotData<-t(cvByDepthTimeBin)
tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotData)[1]),dim(plotData)[2]),"y"=sort(rep(seq(1,dim(plotData)[2]),dim(plotData)[1]))))

depthLabels<-c("0-200 m","200-400 m","400-600 m","600-800 m","800-1300+ m","Sediment")
pdf(paste(plotPath,"ROMSrepeat1yearCVbyDepthTime.pdf",sep=""),height=4,width=7)
par(mar=c(4,7,1.5,1),las=1)
plot(x=seq(1,dim(plotData)[1]),y=rep(dim(plotData)[2],dim(plotData)[1]),type="n",ylab="",xlab="Timestep (years)",xaxt="n",yaxt="n",ylim=c(0.5,dim(plotData)[2]+0.5))
axis(at=seq(1,dim(plotData)[2]),labels = depthLabels,side=2,las=1)
axis(at=seq(1,dim(plotData)[1]),labels=unique(timestepsBinned),side=1,las=1)
temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
dev.off()

xx<-pretty(seq(0,thisMax,length.out = 5)); 
legendText<-xx[xx<thisMax]; legendCols<-unlist(lapply(legendText,getCol))
pdf(paste(plotPath,"ROMSrepeat1yearCVbyDepthTimeLEGEND.pdf",sep=""),height=3,width=2.5)
par(mar=c(0,0,0,0))
makeBlankPlot(); legend(legend=as.character(formatC(legendText,format="f",digits=4)), col=legendCols, x="center",bty="n",pch=15, title="CV",pt.cex=1.5)
dev.off()

##################################
plot(cvByDepthTimeBin[6,])
##just sediment layer CVs by box and time

storeSedCVsByBox<-array(NA,dim=c(nboxes,nts)); dynBoxes<-seq(2,25)
for(b in 1:nboxes){
  if(b %in% dynBoxes){
    thisData<-apply(storeTracersByCell[,,nlayers,b,],c(1,3),sum,na.rm=TRUE)
    thisCVs<-apply(thisData,2,calcCV)
  }
  storeSedCVsByBox[b,]<-thisCVs
}
plot(apply(storeSedCVsByBox[,1:48],2,mean))
abline(v=30)
abline(v=c(25,35))
toPlot<-apply(storeSedCVsByBox[dynBoxIndex,25:35],1,mean); thisMax<-max(toPlot)
plotCol<-unlist(lapply(toPlot,getCol))

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

pdf(paste(plotPath,"CVsSedMap.pdf", sep=""), height=3,width=4)
      plot(shape)
      map('nzHires',add=TRUE,col="black",lwd=2)
      map.axes()
      for(b in 1:dim(labeldf)[1]){
        thisCol<-"white"
        if(b %in% dynBoxes){thisCol<-plotCol[b]}
        polygon(sdata$shp$shp[[b]]$points,col=thisCol)
      }
      mtext("Sediment CVs",side=3,adj=0)
dev.off()



#########################################
## by total box depth
depthByBox<-round(apply(thisDz[,,1],2,sum),-1)[dynBoxIndex] #the rounding just takes off sediment 1 m
boxDepths<-sort(unique(depthByBox)); nboxDepths<-length(boxDepths)
cvByBoxDepthTime<-array(NA,dim=c(nboxDepths,nby))
for(d in 1:nboxDepths){
  thisData<-storeCVsByBox[dynBoxIndex,][depthByBox==boxDepths[d],]
  if(length(dim(thisData))>1){
    xx<-apply(thisData,2,mean,na.rm=TRUE)
  }else{
    xx<-as.double(thisData)
  }
  yy<-tapply(xx[1:plot_nts],timestepsBinned,mean,na.rm=TRUE)
  cvByBoxDepthTime[d,]<-yy
}
thisMax<-max(cvByBoxDepthTime,na.rm=TRUE)
plotColour<-apply(t(cvByBoxDepthTime),c(1,2),getCol)
plotData<-t(cvByBoxDepthTime)
tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotData)[1]),dim(plotData)[2]),"y"=sort(rep(seq(1,dim(plotData)[2]),dim(plotData)[1]))))

depthLabels<-paste(boxDepths," m",sep="")
pdf(paste(plotPath,"ROMSrepeat1yearCVbyBoxDepthTime.pdf",sep=""),height=4,width=6.8)
par(mar=c(4,5,1.5,1),las=1)
plot(x=seq(1,dim(plotData)[1]),y=rep(dim(plotData)[2],dim(plotData)[1]),type="n",ylab="",xlab="Timestep (years)",xaxt="n",yaxt="n",ylim=c(0.5,dim(plotData)[2]+0.5))
axis(at=seq(1,dim(plotData)[2]),labels = depthLabels,side=2,las=1)
axis(at=seq(1,dim(plotData)[1]),labels=unique(timestepsBinned),side=1,las=1)
temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
par(las=0)
mtext("Total box depth",side=2,adj=0.5,line=4)
dev.off()

xx<-pretty(seq(0,thisMax,length.out = 5)); 
legendText<-xx[xx<thisMax]; legendCols<-unlist(lapply(legendText,getCol))
pdf(paste(plotPath,"ROMSrepeat1yearCVbyBoxDepthTimeLEGEND.pdf",sep=""),height=3,width=2.5)
par(mar=c(0,0,0,0))
makeBlankPlot(); legend(legend=as.character(formatC(legendText,format="f",digits=3)), col=legendCols, x="center",bty="n",pch=15, title="CV",pt.cex=1.5)
dev.off()

##################################
## by trophic level
#read in trophic levels
trophicLevelDF<-read.csv(paste(basePath,"inputs\\biol_prm\\CRAM_trophicLevels_isotopes.csv",sep=""))
trophicLevelDF$roundedTL<-myRounding(trophicLevelDF$Isotope,1)
trophicLevelDF$roundedTL[is.na(trophicLevelDF$Isotope)]<-myRounding(trophicLevelDF$TrophicLevel2[is.na(trophicLevelDF$Isotope)],1)

trophicLevels<-sort(unique(trophicLevelDF$roundedTL)); ntls<-length(trophicLevels)

cvByTrophicTime<-array(NA,dim=c(ntls,nby))
for(t in 1:ntls){
  thisCodes<-trophicLevelDF$Code[trophicLevelDF$roundedTL==trophicLevels[t]]
  thisNames<-unlist(lapply(groupsDF$Name[groupsDF$Code %in% thisCodes],str_trim,side="both"))
  thisTracers<-paste(thisNames,"_N", sep=""); tracerIndex<-grep(paste(thisTracers,collapse="|"),Ntracers)
  thisData<-storeCVsByTracer[tracerIndex,]
  if(length(dim(thisData))>1){
    xx<-apply(thisData,2,mean,na.rm=TRUE)
  }else{
    xx<-as.double(thisData)
  }
  yy<-tapply(xx[1:plot_nts],timestepsBinned,mean,na.rm=TRUE)
  cvByTrophicTime[t,]<-yy
}
thisMax<-max(cvByTrophicTime,na.rm=TRUE)
plotColour<-apply(t(cvByTrophicTime),c(1,2),getCol)
plotData<-t(cvByTrophicTime)
tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotData)[1]),dim(plotData)[2]),"y"=sort(rep(seq(1,dim(plotData)[2]),dim(plotData)[1]))))

pdf(paste(plotPath,"ROMSrepeat1yearCVbyTrphicTime_inclZero.pdf",sep=""),height=4,width=6.8)
par(mar=c(4,5,1.5,1),las=1)
plot(x=seq(1,dim(plotData)[1]),y=rep(dim(plotData)[2],dim(plotData)[1]),type="n",ylab="",xlab="Timestep (years)",xaxt="n",yaxt="n",ylim=c(0.5,dim(plotData)[2]+0.5))
axis(at=seq(1,dim(plotData)[2]),labels = trophicLevels,side=2,las=1)
axis(at=seq(1,dim(plotData)[1]),labels=unique(timestepsBinned),side=1,las=1)
temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
par(las=0)
mtext("Trophic level",side=2,adj=0.5,line=3)
dev.off()


pdf(paste(plotPath,"ROMSrepeat1yearCVbyTrphicTime.pdf",sep=""),height=3,width=6.8)
index<-tempDF$y>1
par(mar=c(4,5,1.5,1),las=1)
plot(x=seq(1,dim(plotData)[1]),y=rep(dim(plotData)[2],dim(plotData)[1]),type="n",ylab="",xlab="Timestep (years)",xaxt="n",yaxt="n",ylim=c(1.5,dim(plotData)[2]+0.5))
axis(at=seq(1,dim(plotData)[2]),labels = trophicLevels,side=2,las=1)
axis(at=seq(1,dim(plotData)[1]),labels=unique(timestepsBinned),side=1,las=1)
temp<-Map(plotGrid,x=tempDF$x[index],y=tempDF$y[index])
par(las=0)
mtext("Trophic level",side=2,adj=0.5,line=3)
dev.off()

xx<-pretty(seq(0,thisMax,length.out = 5)); 
legendText<-xx[xx<thisMax]; legendCols<-unlist(lapply(legendText,getCol))
pdf(paste(plotPath,"ROMSrepeat1yearCVbyTrophicTimeLEGEND.pdf",sep=""),height=3,width=2.2)
par(mar=c(0,0,0,0))
makeBlankPlot(); legend(legend=as.character(formatC(legendText,format="f",digits=3)), col=legendCols, x="center",bty="n",pch=15, title="CV",pt.cex=1.5)
dev.off()


#overall, what is the mean and variance of cvs between runs?
cvsBetweenRuns<-apply(storeTracersByCell,c(2,3,4,5),calcCV)

meanAllCvs<-calcCV(cvsBetweenRuns)



#get out the trophic level ones
t=1
thisCodes<-trophicLevelDF$Code[trophicLevelDF$roundedTL==trophicLevels[t]]
thisNames<-unlist(lapply(groupsDF$Name[groupsDF$Code %in% thisCodes],str_trim,side="both"))
thisTracers<-paste(thisNames,"_N", sep=""); tracerIndex<-grep(paste(thisTracers,collapse="|"),Ntracers)
orderedTracers<-Ntracers[tracerIndex]
thisData<-storeTracersByCell[,tracerIndex,,,]; nd<-dim(thisData)[2]
relTracers<-array(NA,dim=c(nd, nruns, dim(thisData)[5])); absTracers<-relTracers
for(d in 1:nd){
  absTracers[d,,]<-apply(thisData,c(2,5),sum,na.rm=TRUE)
}

for(d in 1:nd){
  thisMax<-max(absTracers[d,,], na.rm=TRUE)
  plot(absTracers[d,1,],type="n",ylim=c(0,thisMax))
  for(r in 1:nruns){
    points(absTracers[d,r,],type="l",lwd=2,col=myBlue_trans)
  }
  tempCVs<-apply(absTracers[d,,],2,calcCV); meanCV<-mean(tempCVs)
  tempCVs<-apply(thisData[,1,,,],c(2,3,4),calcCV); meanCV<-nonZeroMean(tempCVs)
}








test<-storeTracersByCell[,10,5,10,30]








