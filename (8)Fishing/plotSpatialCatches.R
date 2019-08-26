# plot forced catch history spatially and temporally for a given species group

#summarise catch by fishing fleet, also which groups are caught in which fishery
library(xlsx)

this_run<-"TBGB"
this_out<-"TBGB_JP2"

this_path<-paste(DIR$'Base',"\\",this_run,"\\",sep="")
inputsPath<-paste(this_path,this_out,"\\",sep="")

plotPath<-paste(DIR$'Figures',"CatchHistories\\TBGB",sep="")


#read in box file
nbox<-25 #this is the number of dynamic boxes

#set catch history path
CH_path<-paste(inputsPath,"catch_history\\",sep="")

getColor<-function(x,thisMax){
  thisCol<-"white"
  if(!is.na(x) & x>0){
    y<-round(x/thisMax,2)*100+1
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
# plotXAxis <- pretty(seq(min(xLabs), max(xLabs),length.out = 5))
# xAxisIndex<-xLabs %in% plotXAxis


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

timeSteps <- c(100, 1000)

thisCode <- "SCA"

thisCatchData <- read.xlsx(paste(CH_path, thisCode,".xlsx", sep=""), sheet=3)
timesteps <- thisCatchData[,1]; plotData <- thisCatchData[,-1]
monthTimeSteps <- (timesteps-1)/30
yearTimeSteps <- (monthTimeSteps+1)/12; yearIndex <- round(yearTimeSteps)==yearTimeSteps
timeAxisAt <- (1:length(timesteps))[yearIndex]; yearAxisLabs <- yearTimeSteps[yearIndex]+1899
thisMax <- max(plotData, na.rm=TRUE)

catchByTime <- data.frame(cbind("catch"=apply(plotData,1,sum, na.rm=TRUE), "time"=round(yearTimeSteps)))
catchByYear <- tapply(catchByTime$catch, catchByTime$time, sum, na.rm=TRUE)
plot(x=names(catchByYear), y=catchByYear, type="h")

###########################
# bring in by fleet and month
fleets <- c(8,12)
thisCatchData <- read.xlsx(paste(CH_path, thisCode,".xlsx", sep=""), sheet=1)
years<-1900:2012; nyears <- length(years)
allData <- array(NA, dim=c(2,24,nyears, 12)) # fleet, box, year, month
for(f in 1:2){
  for(b in 1:24){
    for(y in 1:nyears){
      for(m in 1:12){
        thisRow <- thisCatchData[,1]==years[y] & thisCatchData[,2]==m
        thisCol <- (thisCatchData[1,]==fleets[f] & thisCatchData[2,]==b)
        thisDatum <- thisCatchData[thisRow, thisCol]
        if(length(thisDatum)>0){allData[f,b,y,m]<- thisDatum}
      }
    }
  }
}
catchByYear <- apply(allData, 3, sum, na.rm=TRUE)
catchByBoxYear <- apply(allData, c(2,3), sum, na.rm=TRUE)
plotData <- t(catchByBoxYear); thisMax <- max(plotData, na.rm=TRUE)
xx<-pretty(seq(0, (thisMax), length.out=5)); legendText<-c(xx[signif(xx,2)<thisMax])
legendCols<-unlist(lapply(legendText, getColor, thisMax=thisMax))

pdf(paste(plotPath,"SpatialLEGEND.pdf",sep=""), height=3, width=3)
par(mfrow=c(1,1), mar=c(0,0,0,0))
makeBlankPlot()
legend(legend=legendText/1000, col=legendCols, pch=15, x="center", title="Catch (tonnes)", bty="n")
dev.off()

pdf(paste(plotPath,thisCode, "CatchesBySpaceTime.pdf", sep=""), width=14, height=7)
par(mfrow=c(2,1), mar=c(4,4,1,1))
plotColour<-apply(plotData,c(1,2),getColor,thisMax)
tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotData)[1]),dim(plotData)[2]),"y"=sort(rep(seq(1,dim(plotData)[2]),dim(plotData)[1]))))
par(mar=c(6,4,1,1))
plot(x=seq(1,dim(plotData)[1]),y=rep(dim(plotData)[2],dim(plotData)[1]),type="n",xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,dim(plotData)[2]))
axis(at=(1:nyears)[pretty(1:nyears)],labels = years[pretty(1:nyears)],side=1,las=2)
axis(at=seq(1,dim(plotData)[2]),labels=colnames(plotData),side=2,las=1)
temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
box()
mtext("Box",side=2,adj=0.5,line=3)

plot(x=years, y=catchByYear, type="h", ylab="Catch (tonnes)", xlab="Year",  ylim=c(0, max(catchByYear)), col="midnightblue", lwd=2); 
dev.off()

pdf(paste(plotPath,thisCode,"SpatialCatches.pdf", sep=""))
par(mfrow=c(8,3), mar=c(0,0,0,0))
for(y in 1:nyears){
  plot(shape)
  map('nzHires',add=TRUE,col="black",lwd=2)
  map.axes()
  for(plotB in 1:dim(labeldf)[1]){
    Blab <- plotB 
    if(Blab==(dim(labeldf)[1]+1)){Blab<-1}# check;  don't think I need this, not including boundary box here..
    if(Blab<=dim(plotColour)[2]){
      thisCol<-plotColour[y, Blab]
      test<-sdata$shp$shp[[plotB]]$points
      thisX<-mean(test$'X'); thisY<-mean(test$'Y')
      text(Blab, x=thisX, y=thisY)
      polygon(sdata$shp$shp[[plotB]]$points,col=thisCol)
    }
  }
  mtext(paste(years[y], sep=""), side=3, adj=0, line=-1)
}
dev.off()





