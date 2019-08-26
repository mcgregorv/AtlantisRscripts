#plot selected tracers by box (not layer) comparing multiple runs
library(ggmap)
library(rgdal)
library(gtable)
library(maps)
library(mapdata)
library(rgeos)
source(paste(DIR$'General functions',"\\get_first_number.R",sep=""))
source(paste(DIR$'General functions',"formatShape.R",sep=""))


# this_out<-"outputFishOFF"; 
this_out <-"outputTestsSCA4"
runFolder<-"TBGB_JP2"; 
# this_out<-"output"; runFolder<-"TBGBReportBase";
# this_out <-"MyRun_Fish1899_better1_codeupdate_rewriteDmatrix1"; runFolder<-"TBGBfish"
thisDesc <- paste(runFolder, this_out,sep="")

basePath<-  paste(DIR$'Base', "TBGB\\",runFolder,"\\",sep="")
outPath<-paste(basePath,"\\",this_out,"\\",sep="")
plotPath<-paste(basePath,"..\\Figures\\Testing\\Spatial\\",thisDesc, sep="")

year0<-1899

daysTimeStep<-73; timestepsPerYear <- 365/daysTimeStep

  ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
  thisVol<-ncvar_get(ThisNC.nc,"volume")
  thisDz<-ncvar_get(ThisNC.nc,"dz"); nlayers<-dim(thisVol)[1]
 
max_nts<-dim(thisVol)[3]; nts<-max_nts
max_ndays<-max_nts * daysTimeStep; 
nboxes<-dim(thisVol)[2]
xLabsTemp<-seq(0,(max_nts),by=timestepsPerYear)
xLabsAt<-xLabsTemp
xLabs<-xLabsTemp*daysTimeStep/365+year0

## get allTracers, then select ones to plot
allTracers<-sort(names(ThisNC.nc$var))
temp <- allTracers[grep("_N", allTracers)]; Ntracers<-temp[grep("Nums", temp, invert = TRUE)]

allPlotTracers<-c(Ntracers,"Oxygen", "NO3", "NH3")
nPT<-length(allPlotTracers)

#populate array of tracers
tracersArray<-array(NA, dim=c( nPT,nboxes, max_nts))

  for(t in 1:nPT){
    thisTracer<-allPlotTracers[t]
    thisData<-ncvar_get(ThisNC.nc, thisTracer)
    if(length(dim(thisData))==2){
      tracersArray[t,,1:max_nts]<-thisData * thisVol[nlayers,,] * X_CN * mg_2_tonne
    } else{
      tracersArray[t,,1:max_nts]<-apply(thisData * thisVol,c(2,3), sum, na.rm=TRUE) * mg_2_tonne * X_CN
    }
  }

getColor<-function(x,thisMax){
  thisCol<-"white"
  if(!is.na(x) & x>0){
    y<-round(x/thisMax,2)*100+1
    thisCol<-thisColRamp[y]
  }
  return(thisCol)
}
## log scale
getColor<-function(x,thisMax){
  thisCol<-"white"
  if(!is.na(x) & x>0){
    if(x<1){
      y<-1
    } else{
      y<-round(log10(x)/log10(thisMax),2)*100+1
    }
    thisCol<-thisColRamp[y]
  }
  return(thisCol)
}
thisColRamp <- colorRampPalette(colors=c(myLightBlue, myBlue, "midnightblue"))(101)


plotGrid<-function(x,y){
  thisX<-c(x-0.5,x-0.5,x+0.5,x+0.5); thisY<-c(y-0.5,y+0.5,y+0.5,y-0.5)
  thisCol<-plotColour[x,y]
  if(length(thisCol)>0){
    polygon(x=thisX,y=thisY,col=thisCol,border=NA)
  }
  return(NULL)
}
plotXAxis <- pretty(seq(min(xLabs), max(xLabs),length.out = 5))
xAxisIndex<-xLabs %in% plotXAxis

xx<-10^pretty(seq(log10(0.1), log10(thisMax), length.out=5)); legendText<-c(xx[xx<=thisMax],signif(thisMax,2))
legendCols<-unlist(lapply(legendText, getColor, thisMax=thisMax))

pdf(paste(plotPath,"SpatialLEGEND.pdf",sep=""), height=3, width=3)
par(mfrow=c(1,1), mar=c(0,0,0,0))
makeBlankPlot()
legend(legend=legendText, col=legendCols, pch=15, x="center", title="Biomass (tonnes)", bty="n")
dev.off()


## read the shape file in for snapshot comparisons
source(paste(DIR$'General functions',"formatShape.R",sep=""))
#
thisPath<-paste(DIR$'Base',"\\TBGB\\TBGBSpatialPolygons\\",sep="")
# 
# plotPath<-paste(DIR$'Base',"TBGB\\Figures\\",sep="")
shapeFile<-paste(thisPath,"P0to26_depth",sep="")

sdata<-read.shapefile(shapeFile)
shape<-formatShape(shapeFile=shapeFile)
ns<-length(shape)
SpDF <- SpatialPolygonsDataFrame(shape,data.frame( z=1:ns, row.names=paste("P",seq(1,(ns)),sep="")))
labels<-seq(1,(ns))
plot(shape)
LABELpOS<-polygonsLabel(shape, labels = labels, cex=.1,doPlot=FALSE)
labeldf<-data.frame(cbind("x"=LABELpOS[1:ns],"y"=LABELpOS[(ns+1):(2*ns)]))

thisTime<-1
thisTime<-112

timeSteps <- c(1, nts)

pdf(paste(plotPath,"SpatialMapSnapshots_allGroups.pdf", sep=""), width=10)
par(mfrow=c(3,4), mar=c(4,4,1,1))
for(t in 1:nPT){
  plotData<-t(tracersArray[t,,]); thisMax <- max(plotData, na.rm=TRUE)
  plotColour<-apply(plotData,c(1,2),getColor,thisMax)
  tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotData)[1]),dim(plotData)[2]),"y"=sort(rep(seq(1,dim(plotData)[2]),dim(plotData)[1]))))
  par(mar=c(6,4,1,1))
  plot(x=seq(1,dim(plotData)[1]),y=rep(dim(plotData)[2],dim(plotData)[1]),type="n",xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,dim(plotData)[2]))
  mtext(allPlotTracers[t],side=3, adj=0)
  axis(at=xLabsAt,labels = xLabs,side=1,las=2)
  axis(at=seq(1,dim(plotData)[2]),labels=colnames(plotData),side=2,las=1)
  temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
  box()
  mtext("Box",side=2,adj=0.5,line=3)
  
  xx<-apply(plotData, 1, sum, na.rm=TRUE)
  par(mar=c(10,4,1,1))
  plot(xx, type="l", ylab="Biomass (tonnes)", xlab="Year", xaxt="n", ylim=c(0, max(xx)), col="midnightblue", lwd=2); 
  axis(at=xLabsAt[xAxisIndex],labels = xLabs[xAxisIndex],side=1,las=2)
  abline(h=xx[1], col="red", lwd=2, lty=2)
  
  
   for(thisTime in timeSteps){
     if(thisTime > dim(plotColour)[1]){thisTime <- dim(plotColour)[1]-1}
    plot(shape)
    map('nzHires',add=TRUE,col="black",lwd=2)
    map.axes()
    for(plotB in 1:dim(labeldf)[1]){
      Blab <- plotB +1
      if(Blab==(dim(labeldf)[1]+1)){Blab<-1}
      thisCol<-plotColour[thisTime, Blab]
      test<-sdata$shp$shp[[plotB]]$points
      thisX<-mean(test$'X'); thisY<-mean(test$'Y')
      text(Blab, x=thisX, y=thisY)
      polygon(sdata$shp$shp[[plotB]]$points,col=thisCol)
    }
    mtext(paste(allPlotTracers[t],", timestep ",thisTime, sep=""), side=3, adj=0)
  }
}
dev.off()

test<- ncvar_get(ThisNC.nc, "MicroPB_N")

# 
# pdf(paste(plotPath,"SpatialHeatMaps_allGroups.pdf",sep=""))
# par(mfrow=c(3,2), mar=c(4,4,1,1))
# for(t in 1:nPT){
#   plotData<-t(tracersArray[t,,]); thisMax <- max(plotData, na.rm=TRUE)
#   plotColour<-apply(plotData,c(1,2),getColor,thisMax)
#   tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotData)[1]),dim(plotData)[2]),"y"=sort(rep(seq(1,dim(plotData)[2]),dim(plotData)[1]))))
#   par(mar=c(6,4,1,1))
#   plot(x=seq(1,dim(plotData)[1]),y=rep(dim(plotData)[2],dim(plotData)[1]),type="n",xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,dim(plotData)[2]))
#   mtext(allPlotTracers[t],side=3, adj=0)
#   axis(at=xLabsAt[xAxisIndex],labels = xLabs[xAxisIndex],side=1,las=2)
#   axis(at=seq(1,dim(plotData)[2]),labels=colnames(plotData),side=2,las=1)
#   temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
#   box()
#   mtext("Box",side=2,adj=0.5,line=3)
#   
#   xx<-apply(plotData, 1, sum, na.rm=TRUE)
#   par(mar=c(10,4,1,1))
#   plot(xx, type="l", ylab="Biomass (tonnes)", xlab="Year", xaxt="n", ylim=c(0, max(xx)), col="midnightblue", lwd=2); 
#   axis(at=xLabsAt[xAxisIndex],labels = xLabs[xAxisIndex],side=1,las=2)
#   abline(h=xx[1], col="red", lwd=2, lty=2)
# }
# dev.off()
