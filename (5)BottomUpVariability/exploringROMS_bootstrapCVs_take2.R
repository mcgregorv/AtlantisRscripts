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
Version<-"shortB"
plotPath<-paste(DIR$'Figures',"ROMS\\testing",Version,sep="")

# outPath<-paste(basePath,"base\\ROMS1yearRepeat\\",sep="")
outPath<-paste(basePath,"base\\ouputROMS",Version,"\\",sep="")
outPath<-paste(basePath,"base\\outputROMS",Version,"\\",sep="")
dataOutPath<-paste(basePath,"BootstrapROMS\\modelOut",Version,"_", sep="")

load(paste(dataOutPath,"modelTracers",sep="")); ##to bring "to bring "storeTracers", "storeTracersByCell", "storeTemperature", "store_nts".

#get _N tracers
BaseNC.nc<-nc_open(paste(outPath, "outputROMSBootstrap",Version,"1\\output.nc",sep=""))
allTracers<-names(BaseNC.nc$var)
thisVol<-ncvar_get(BaseNC.nc, "volume"); nlayers<-dim(thisVol)[1]; nts<-dim(thisVol)[3]
thisDZ<-ncvar_get(BaseNC.nc,"dz")
x<-grep("_N",allTracers); temp<-allTracers[x]
y<-grep("Nums",temp, invert = TRUE); Ntracers<-temp[y]; ntracers<-length(Ntracers)

storeCVsByCellAndTracer<-array(NA, dim=c(dim(storeTracersByCell)[-1]))
testGroup<-"DinoF"; t<-grep(testGroup,Ntracers)

for(t in 1:ntracers){
  testData<-storeTracersByCell[,t,,,]
  testCVs<-array(NA, dim=dim(testData)[-1])
  for(l in 1:dim(testCVs)[1]){
    for(b in 1:dim(testCVs)[2]){
      for(y in 1:dim(testCVs)[3]){
        this_x<-testData[,l,b,y]
        thisVar<-var(this_x); thisMean<-mean(this_x); thisCV<-max(0,sqrt(thisVar)/thisMean)
        storeCVsByCellAndTracer[t,l,b,y]<-thisCV
      }
    }
  }
}
skipCarrion<-seq(1,ntracers)==grep("Carrion",Ntracers)

maxByGroup<-apply(storeCVsByCellAndTracer,1,max, na.rm=TRUE)
nzMeanByGroup<-apply(storeCVsByCellAndTracer,1,nonZeroMean)

plotNames<-gsub("_|_N"," ",Ntracers[!skipCarrion][rev(order(nzMeanByGroup[!skipCarrion]))])

par(lend=1, mar=c(8,4,1,1))
plot(100*nzMeanByGroup[!skipCarrion][rev(order(nzMeanByGroup[!skipCarrion]))], type="h", lwd=5, xaxt="n", xlab="", ylab="CV")
par(las=2)
axis(at=seq(1,ntracers-1), labels = plotNames, side=1)

nzMeanByBox<-apply(storeCVsByCellAndTracer[,-6,,],3,nonZeroMean) #wihtout sediment

dynBoxes<-seq(2,25)
volByBox<-apply(thisVol[,,1],2,sum)
plot(x=volByBox[dynBoxes], y=nzMeanByBox[dynBoxes], col=myOrange, pch=20)

depthByBox<-apply(thisDZ[,,1],2,sum)
plot(x=depthByBox[dynBoxes], y=nzMeanByBox[dynBoxes], col=myOrange, pch=20)

# color boxes on map
#for spatial, plot max cv by box
thisColRamp<-colorRampPalette(colors=c(myLightAqua,myBlue,"midnightblue"))(101)
thisMax<-myRounding(max(nzMeanByBox, na.rm=TRUE), 0.05, direction="up"); thisMin<-myRounding(min(nzMeanByBox, na.rm=TRUE), 0.05, direction="down")
boxColors<-unlist(lapply(nzMeanByBox, getColor, thisMax=thisMax))

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

calcCV<-function(x){
  thisVar<-var(x, na.rm=TRUE); thisMean<-mean(x, na.rm=TRUE)
  thisCV<-sqrt(thisVar)/thisMean
  return(thisCV)
}

#####################################
## by depth
dynBoxIndex<-seq(2,25)
timeIndex<-seq(10,40);
timeIndex<-seq(1,nts);
nits<-length(timeIndex)
storeCVsbyDepth<-array(NA,dim=c(nlayers,nits))
layerIndex<-c(seq((nlayers-1),1),nlayers)
for(l in 1:nlayers){
  thisLayer<-layerIndex[l]
  # thisData<-storeCVsByCellAndTracer[,thisLayer,,]
  # thisNZmeanByTime<-apply(thisData,3,mean, na.rm=TRUE)
  # storeCVsbyDepth[l,]<-thisNZmeanByTime
  
  ## alt 
  thisData<-storeTracersByCell[,,thisLayer,,timeIndex]; thisX<-apply(thisData,c(1,4), sum, na.rm=TRUE)
  thisCbyTime<-apply(thisX,2,calcCV)
  
  storeCVsbyDepth[thisLayer,]<-thisCbyTime
}
colByDepth<-colorRampPalette(colors=c(myAqua,myBlue,'midnightblue'))(nlayers)
thisMax<-max(storeCVsbyDepth,na.rm=TRUE)
plot(storeCVsbyDepth[1,],type="n",ylim=c(0,1.2*thisMax))
for(l in 1:nlayers){
  points(storeCVsbyDepth[l,],pch=20,col=colByDepth[l])
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

## check out all layers by tracer
for(l in 1:nlayers){
  thisLayer<-layerIndex[l]
  pdf(paste(plotPath,"CVsbyTracer_layer",thisLayer,".pdf", sep=""))
  par(mfrow=c(5,2), mar=c(4,4,1,1))
  thisData<-storeTracersByCell[,,thisLayer,,]
  for(t in 1:ntracers){
    thisX<-apply(thisData[,t,,],c(1,3),sum, na.rm=TRUE)
    thisCV<-apply(thisX,2,calcCV)
    if(sum(thisCV, na.rm=TRUE)>0){
      plot(thisCV, pch=20, col=myBlue);mtext(Ntracers[t],side=3,adj=0)
    }
  }
  dev.off()
}

## check out all layers by tracer
for(l in 1:nlayers){
  thisLayer<-layerIndex[l]
  pdf(paste(plotPath,"Tracer_layer",thisLayer,".pdf", sep=""))
  par(mfrow=c(5,2), mar=c(4,4,1,1))
  thisData<-storeTracersByCell[,,thisLayer,,]
  for(t in 1:ntracers){
    thisX<-apply(thisData[,t,,],c(1,3),sum, na.rm=TRUE)
    thisMax<-max(thisX, na.rm=TRUE)
    if(thisMax>0){
      plot(thisX[1,], type="l",col=myGrey_trans,lwd=2);mtext(Ntracers[t],side=3,adj=0)
      for(r in 2:50){
        points(thisX[r,], type="l",col=myGrey_trans,lwd=2)
      }
    }
  }
  dev.off()
}

## loop through each box and plot boxplot of group cvs by layer
nboxes<-dim(storeTracersByCell)[4]
getCut<-function(x, cut=0.33){
  yy<-sort(x); z<-round(cut*length(x))
  cut<-yy[z]
  return(cut)
}

timeIndex<-seq(10,40)

getCut<-function(x, cut=0.98){
  total<-sum(x, na.rm=TRUE)
  yy<-rev(sort(x)); yyCum<-cumsum(yy); z<-min(yyCum[yyCum> cut*total]); cut<-yy[yyCum==z]
  return(cut)
}
pdf(paste(plotPath,"CVsByBoxLayerTracer.pdf", sep=""))
storeTopGroupScores<-rep(0,ntracers)
for(B in dynBoxes){
  # if(B %in% dynBoxes){
    par(mfrow=c(5,2))
    # jpeg(paste(plotPath,"CVsByBoxLayerTracer",B,"_map.jpg", sep=""))
    # plot(shape)
    # map('nzHires',add=TRUE,col="black",lwd=2)
    # map.axes()
    # polygon(sdata$shp$shp[[B]]$points,col=myBlue)
    # dev.off()
    thisData<-storeCVsByCellAndTracer[,,B,timeIndex]
    temp<-apply(thisData, c(1,2), nonZeroMean); test<-apply(temp, 1, sum, na.rm=TRUE); lowerCut<-getCut(test, cut=0.9)
    cutIndex<-test>=lowerCut; keep<-temp[cutIndex,]; thisGroupNames<-gsub("_", " ", Ntracers[cutIndex])
    storeTopGroupScores[cutIndex]<-storeTopGroupScores[cutIndex]+1
    toPlot<-melt(keep); toPlot$value[toPlot$value=="NaN"]<-NA; colnames(toPlot)<-c("Group", "Layer", "CV")
    toPlot$Group<-thisGroupNames[match(toPlot$Group, seq(1,ntracers))]; toPlot$Layer<-factor(toPlot$Layer, levels=seq(1,nlayers))
    bp <- ggplot(data = toPlot, aes(x = Layer, fill = Group, y = CV)) + 
      geom_bar(stat = 'identity')
    pdf(paste(plotPath,"CVsByBoxLayerTracer",B,"_boxPlot.pdf", sep=""))
   print( bp)
    dev.off()

  # }
}
# dev.off()
nzIndex<-storeTopGroupScores>0; plotNames<-gsub("_|N", " ", rev(Ntracers[nzIndex][order(storeTopGroupScores[nzIndex])]))
par(mar=c(8,4,1,1), mfrow=c(1,1), lend=1)
plot(rev(sort(storeTopGroupScores[nzIndex])), type="h", lwd=5, xaxt="n", xlab="")
par(las=2)
axis(at=seq(1,length(storeTopGroupScores[nzIndex])), labels=plotNames, side=1)

## plot maps for the top ones
pdf(paste(plotPath,"TopTracerCVsSpatial.pdf", sep=""))
par(mfrow=c(5,2), mar=c(4,4,1,1))
for(s in 1:ntracers){
  if(rev(sort(storeTopGroupScores))[s]>0){
    t<-rev(order(storeTopGroupScores))[s]
    thisData<-storeCVsByCellAndTracer[t,,,]
    
    nzMeanByBox<-apply(storeCVsByCellAndTracer[t,,,],3,nonZeroMean) #wihtout sediment
    
    thisMax<-myRounding(max(nzMeanByBox, na.rm=TRUE), 0.05, direction="up"); thisMin<-myRounding(min(nzMeanByBox, na.rm=TRUE), 0.05, direction="down")
    boxColors<-unlist(lapply(nzMeanByBox, getColor, thisMax=thisMax))
    plot(shape)
    map('nzHires',add=TRUE,col="black",lwd=2)
    map.axes()
    for(b in 1:dim(labeldf)[1]){
      thisCol<-"white"
      if(b %in% dynBoxes){thisCol<-boxColors[b]}
      polygon(sdata$shp$shp[[b]]$points,col=thisCol)
    }
    par(las=0)
    mtext(paste(gsub("_|N", " ", Ntracers[t]), " count:",storeTopGroupScores[t],sep=""), side=3, adj=1, line=0)
  }
}
dev.off()


Ntracers[storeTopGroupScores==0]

for(B in 1:nboxes){
  if(B %in% dynBoxes){
    par(mfrow=c(5,2))
    pdf(paste(plotPath,"CVsByBoxLayerTracer",B,"_map.pdf", sep=""))
    plot(shape)
    map('nzHires',add=TRUE,col="black",lwd=2)
    map.axes()
    polygon(sdata$shp$shp[[B]]$points,col=myBlue)
    dev.off()
  }
}
## add them to latex input file
writeFile<-paste(DIR$'Reports',"(p1)ROMSBootstrap\\testingCVs.txt", sep="")
cat("", file=writeFile, append=FALSE)
for(b in 1:nboxes){
  if(b %in% dynBoxes){
    thisText<-paste("	\\begin{figure}[H]
		\\centering
		\\includegraphics[height=7cm]{testingBCVsByBoxLayerTracer",b,"_map.pdf}
     \\includegraphics[height=8cm]{testingBCVsByBoxLayerTracer",b,"_boxPlot.pdf}
		\\caption{.}\\label{gif:box",b,"}
	\\end{figure}
                    ",sep="")
    cat(thisText, file=writeFile, append=TRUE)
  }
}


## box 24 has the highest cvs - check which layers and species this is
thisData<-storeCVsByCellAndTracer[,,25,]
nzMeanByLayerTracer<-apply(thisData, c(1,2), nonZeroMean)
test<-apply(thisData, 1, nonZeroMean)
thisSubSet<-nzMeanByLayerTracer[test>0.5,]
colByLayer<-colorRampPalette(colors=c(myAqua,myBlue,"midnightblue"))(nlayers)
par(mar=c(8,4,1,1))
plot(thisSubSet[,1],ylim=c(0,max(thisSubSet,na.rm = TRUE)), type="n", xaxt="n", xlab="")
for(l in 1:nlayers){
  points(thisSubSet[,l],pch=l, col=colByLayer[l], cex=1.5)
}
par(las=2)
axis(at=seq(1,dim(thisSubSet)[1]), labels = Ntracers[test>0.5], side=1)

t=49 # ZM
testZM<-storeTracersByCell[,t,,24,]
thisData<-storeCVsByCellAndTracer[t,,,]
meanByLayerBox<-apply(thisData,c(1,2), nonZeroMean)

## how are diatom, picophytopl and zoo over the other boxes?
t=grep("Diatom",Ntracers)
t=grep("PicoPhy", Ntracers)
t=grep("^Zoo_N",Ntracers)
thisData<-storeCVsByCellAndTracer[t,,,]
meanByBox<-apply(thisData, 2, nonZeroMean)

par(las=1)
plot(meanByBox)
mtext(Ntracers[t],side=3)
