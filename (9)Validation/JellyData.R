library(nzPlot)

#read in data
thisFile<-paste(DIR$'Data',"NZ_Medusa_jellyData.csv", sep="")
thisData<-read.csv(thisFile)
plotLon<-unlist(lapply(thisData$Lon,FUN=function(x){360+x}))

thisMax<-max(thisData[,c(3:dim(thisData)[2])]); thisMin<-min(thisData[,c(3:dim(thisData)[2])])
getCol<-function(x){
  thisCol<-"white"
  y=round(((x-thisMin)/(thisMax-thisMin)),2)*100+1
  if(sum(x,na.rm=TRUE)>0){
    thisCol<-thisColRamp[y]
  }
  # return(y)
  return(thisCol)
}
thisColRamp<-colorRampPalette(colors=c(myLightAqua,myBlue,"midnightblue"))(101)

thisDataCols<-apply(thisData[,3:dim(thisData)[2]], c(1,2), getCol)

nyears<-dim(thisDataCols)[2]; years<-gsub("X","",colnames(thisData)[3:dim(thisData)[2]])

pdf(paste(DIR$'Figures',"ChathamRiseNZPlotMap.pdf", sep=""), height=10, width=8)
par(mar=c(3,5,1,0),oma=c(0,0,0,0), mfrow=c(5,2))
for(y in 1:nyears){
  thisYear<-years[y]
  nz(xlim = c(165, 187), ylim = c(-50, -34),xlab=" ",ylab=" ",cex=thisCex,fill.col="DarkOliveGreen3")
  nz.depth(contour = 1000, col="lightblue"); nz.depth(contour = 500, col="lightblue"); nz.depth(contour = 200, col="lightblue");
  
  nz.points(x=plotLon, y=thisData$Lat, col=thisDataCols[,y], pch=20)
  # tidy up land
  nz.polygon(nz.coast$long, nz.coast$lat,  col = "DarkOliveGreen3",  border = T, xpd = F)
  nz.polygon(nz.islands$long, nz.islands$lat,  col = "DarkOliveGreen3",  border = T, xpd = F)
  nz.text("200 m, 500 m and\n1000 m contours",x=183.8, y=-34.8, col="lightblue", cex=0.5)
  # nz.text("Chatham Rise", x=179, y=-43.4, col="black", font=2)
  mtext(thisYear, side=3,adj=0)
}
dev.off()

xx<-pretty(seq(thisMin, thisMax, length.out=10))
legendText<-xx[xx<thisMax]; legendCols<-unlist(lapply(legendText, getCol))

pdf(paste(DIR$'Figures',"ChathamRiseNZPlotMapLEGEND.pdf", sep=""),height=2,width=2)
makeBlankPlot()
legend(legend=legendText,col=legendCols,pch=20,x="center", bty="n",pt.cex=1.5)
dev.off()

## median over all years
medBySpace<-apply(thisData[,c(3:dim(thisData)[2])],1,median,na.rm=TRUE); medCols<-unlist(lapply(medBySpace, getCol))
medByYear<-apply(thisData[,c(3:dim(thisData)[2])], 2, median, na.rm=TRUE)


pdf(paste(DIR$'Figures',"MedusaByYearMedian.pdf", sep=""), height=5, width=7)
par(mar=c(4,5,1,1))
plot(x=years,y=medByYear, type="l",lwd=2,col=myBlue, xlab="Year", ylab=expression("Medusa (mg C " ~ m^3 ~")"))
dev.off()

pdf(paste(DIR$'Figures',"ChathamRiseNZPlotMapMedian.pdf", sep=""), height=5, width=5)
par(mar=c(3,5,1,0),oma=c(0,0,0,0), mfrow=c(1,1))
  nz(xlim = c(165, 187), ylim = c(-50, -34),xlab=" ",ylab=" ",cex=thisCex,fill.col="DarkOliveGreen3")
  nz.depth(contour = 1000, col="lightblue"); nz.depth(contour = 500, col="lightblue"); nz.depth(contour = 200, col="lightblue");
  
  nz.points(x=plotLon, y=thisData$Lat, col=medCols, pch=20)
  # tidy up land
  nz.polygon(nz.coast$long, nz.coast$lat,  col = "DarkOliveGreen3",  border = T, xpd = F)
  nz.polygon(nz.islands$long, nz.islands$lat,  col = "DarkOliveGreen3",  border = T, xpd = F)
  nz.text("200 m, 500 m and\n1000 m contours",x=183.8, y=-34.8, col="lightblue", cex=0.5)
dev.off()


lonLim<-c(165,195); latLim=c(-58,-30)
par(mar=c(4,4,1,1),oma=c(0,2,0,0))
map('nzHires',add=FALSE,col="black",lwd=2,xlim = lonLim,ylim=latLim)
par(las=1)
map.axes()
points(x=plotLon, y=thisData$Lat, col=medCols, pch=20)



