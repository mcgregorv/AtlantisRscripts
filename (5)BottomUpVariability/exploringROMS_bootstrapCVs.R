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
source(paste(DIR$'General functions',"myRounding.R",sep=""))

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

basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

Version<-"D"
Version<-"B"
# outPath<-paste(basePath,"base\\ROMS1yearRepeat\\",sep="")
# outPath<-paste(basePath,"base\\ouputROMS",Version,"\\",sep="")
outPath<-paste(basePath,"base\\outputROMS",Version,"\\",sep="")
dataOutPath<-paste(basePath,"BootstrapROMS\\modelOut",Version,"_", sep="")

plotPath<-paste(DIR$'Figures',"ROMS\\",Version,sep="")

nruns<-50; nts<-51; burnin<-35 ## this is for bootstrap ROMS
# nruns<-9; nts<-51; burnin<-35 ##this is for repeat each ROMS

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
BaseNC.nc<-nc_open(paste(outPath, "outputROMSBootstrap",Version,"1\\output.nc",sep=""))
allTracers<-names(BaseNC.nc$var)
thisVol<-ncvar_get(BaseNC.nc, "volume"); nlayers<-dim(thisVol)[1]; nts<-dim(thisVol)[3]
x<-grep("_N",allTracers); temp<-allTracers[x]
y<-grep("Nums",temp, invert = TRUE); Ntracers<-temp[y]; ntracers<-length(Ntracers)

#set up array to store outputs
storeTracers<-array(NA, dim=c(nruns, ntracers,nts)); storeTracersByCell<-array(NA, dim=c(nruns, ntracers,nl,nb,nts))
storeTemperature<-array(NA,dim=c(nruns,nl,nb,nts)); 
store_nts<-rep(NA,nruns)
for(r in 1:nruns){
  cat(r," -- ")
  # thisOutPath<-paste(outPath,"outputBASE50yr",r,"\\",sep="")
  thisOutPath<-paste(outPath, "outputROMSBootstrap",Version,r,"\\",sep="")
  thisOutFile<-paste(thisOutPath,"output.nc",sep="")
  ThisNC.nc<-nc_open(thisOutFile)
  thisVol<-ncvar_get(ThisNC.nc,"volume")
  thisTracer<-"Temp"
  thisData<-ncvar_get(ThisNC.nc, thisTracer); this_nts<-dim(thisData)[3]
  storeTemperature[r,,,1:min(nts,this_nts)]<-thisData[,,1:min(nts,this_nts)]
  store_nts[r]<-this_nts
  for(t in 1:ntracers){
    thisTracer<-Ntracers[t]; thisData<-ncvar_get(ThisNC.nc, thisTracer)
    if(length(dim(thisData))==3){
      xx<-apply(thisData*thisVol,3,sum)*mg_2_tonne*X_CN
      yy<-thisData*thisVol*mg_2_tonne*X_CN
      storeTracersByCell[r,t,,,1:min(nts,this_nts)]<-yy[,,1:min(nts,this_nts)]
    } else{
      xx<-apply(thisData*thisVol[nl,,],2,sum)*mg_2_tonne*X_CN
      yy<-thisData*thisVol[nl,,]*mg_2_tonne*X_CN
      storeTracersByCell[r,t,nl,,1:min(nts,this_nts)]<-yy[,1:min(nts,this_nts)]
    }
    storeTracers[r,t,1:min(nts,this_nts)]<-xx[1:min(nts,this_nts)]
  }
}


#write them out - commented out so only write out if intending to
# save(list=c("storeTracers", "storeTracersByCell", "storeTemperature", "store_nts"),file=paste(dataOutPath,"modelTracers",sep=""))

# 
# load(paste(dataOutPath,"modelTracers",sep="")); ##to bring "storeTracers", "storeTracersByCell", "storeTemperature", "store_nts"
# ## . dim = dim=c(nruns, ntracers, nyears), c(nruns, ntracers,10, nyears), c(nruns, ntracers,10, nyears), c(nruns, ntracers, nyears) resp.


## do straight _N tracers
pdf(paste(plotPath,"TracersFromROMSbootstraps.pdf", sep=""), height=10)
par(mfrow=c(7,1), mar=c(4,4,1,1))
for(t in 1:ntracers){
  thisTracer<-Ntracers[t]
  thisData<-storeTracers[,t,]
  thisMax<-max(thisData,na.rm=TRUE)*1.1
  plot(thisData[1,], type="n",ylim=c(0,thisMax),xlab="Timestep (years)",ylab="Tonnes")
  for(r in 1:nruns){
    points(thisData[r,],type="l",lwd=2,col=myGrey_trans)
  }
  mtext(thisTracer,side=3,adj=0)
}
dev.off()

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

#what is the variance difference when summarised by timestep?
timeTracerCVs<-apply(storeTracers,c(2,3),calcCV)
test<-apply(timeTracerCVs,1,max,na.rm=TRUE)

#############################################################################################
## cvs by tracers again
storeCVsbyTracer<-array(NA,dim=c(ntracers,nts))
for(t in 1:ntracers){
  thisData<-storeTracersByCell[,t,,,]
  #add up tracers for each run/timestep
  sumTracers<-apply(thisData,c(1,4), sum, na.rm=TRUE)
  thisCVwrtTime<-apply(sumTracers,2,calcCV)
  storeCVsbyTracer[t,]<-thisCVwrtTime
}
thisColRamp<-colorRampPalette(colors=c(myLightAqua,myBlue,"midnightblue"))(101)

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

xx<-pretty(seq(thisMin,thisMax,length.out = 5)); 
legendText<-xx[xx<thisMax]; legendCols<-unlist(lapply(legendText,getColor, thisMax=thisMax))
pdf(paste(plotPath,"CVbyTracerTimeLEGEND.pdf",sep=""),height=3,width=2.5)
par(mar=c(0,0,0,0))
makeBlankPlot(); legend(legend=as.character(formatC(legendText,format="f",digits=4)), col=legendCols, x="center",bty="n",pch=15, title="CV",pt.cex=1.5)
dev.off()

maxCVsbyTracer<-apply(storeCVsbyTracer, 1, max, na.rm=TRUE)

indexMax<-maxCVsbyTracer>0.1
maxTracers<-gsub("_N|_"," ",Ntracers[indexMax])
pdf(paste(plotPath,"CVmaxBytracer_above58percent.pdf", sep=""), height=4, width=5)
par(lend=1, mar=c(9,4,1,1))
plot(maxCVsbyTracer[indexMax], type="h", lwd=5, xaxt="n", xlab="", ylab="CV")
par(las=2)
axis(at=seq(1, length(maxTracers)), labels=maxTracers, side=1)
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
colByDepth<-colorRampPalette(colors=c(myAqua,myBlue,'midnightblue'))(nlayers)
thisMax<-max(storeCVsbyDepth,na.rm=TRUE)
plot(storeCVsbyDepth[1,],type="n",ylim=c(0,1.2*thisMax))
for(l in 1:nlayers){
  points(storeCVsbyDepth[l,],pch=20,col=colByDepth[l])
}

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

dynBoxIndex<-seq(2,25)
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
thisMax<-myRounding(max(storeCVsbyBox, na.rm=TRUE), 0.002, direction="up"); thisMin<-myRounding(min(storeCVsbyBox, na.rm=TRUE), 0.05, direction="down")
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

dynBoxes<-seq(2,25)
pdf(paste(plotPath,"CVsBoxMap.pdf", sep=""), height=3,width=4)
plot(shape)
map('nzHires',add=TRUE,col="black",lwd=2)
map.axes()
for(b in 1:dim(labeldf)[1]){
  thisCol<-"white"
  if(b %in% dynBoxes){thisCol<-boxColors[b]}
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
TLpath<-paste(DIR$'Base',"ATLANTISmodels\\base\\EWEbase\\", sep="")
trophicLevelDF<-read.csv(paste(TLpath,"CRAMGroupsTL.csv",sep=""))
trophicLevelDF$roundedTL<-myRounding(trophicLevelDF$TL,1,direction = "down")

trophicLevels<-sort(unique(trophicLevelDF$roundedTL)); ntls<-length(trophicLevels)

cvByTrophicTime<-array(NA,dim=c(ntls,nts))
for(t in 1:ntls){
  thisCodes<-trophicLevelDF$Code[trophicLevelDF$roundedTL==trophicLevels[t]]
  thisNames<-unlist(lapply(groupsDF$Name[groupsDF$Code %in% thisCodes],str_trim,side="both"))
  thisTracers<-paste(thisNames,"_N", sep=""); tracerIndex<-grep(paste(thisTracers,collapse="|"),Ntracers)
  thisData<-storeTracersByCell[,tracerIndex,,,]
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







