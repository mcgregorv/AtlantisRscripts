library(ggmap)
library(rgdal)
library(gtable)
library(maps)
library(mapdata)
library(rgeos)
source(paste(DIR$'General functions',"\\get_first_number.R",sep=""))
source(paste(DIR$'General functions',"formatShape.R",sep=""))

################################
#doing for BAL first
thisGroup<-"Chl_a"
################################

this_out<-paste("test",sep="")
this_path<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
outPath<-paste(this_path,"\\",this_run,"\\output",this_out,"\\",sep="")
plotPath<-paste(this_path,"\\Figures\\",this_run,"\\",this_out,"",sep="")

this_nbox<-length(sdata$shp$shp)
boxes<-seq(1,this_nbox)

groupsDF<-read.csv(paste(this_path,"CRAM_groups.csv",sep=""))

if(thisGroup %in% groupsDF$Code){
  thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisGroup],side="both")
} else{
  thisName<-thisGroup; thisVar<-thisGroup
}

#read in nc file
ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")

nts<-dim(thisVol)[3] #number of timesteps

if(thisVar !=thisName){thisVar<-paste(thisName,"_N",sep="")}
temp<-ncvar_get(ThisNC.nc,thisVar)
thisData<-apply(temp*thisVol,c(2,3),sum)
thisMax<-max(thisData)*1.01
nonZeroMin<-thisMax
for(t in 1:nts){
  x<-thisData[,t]
  xx<-x[x>0]
  thisMin<-min(xx)*0.99
  if(thisMin<nonZeroMin){nonZeroMin<-thisMin}
}
tracerCols<-c("white",colorRampPalette(colors=c(myGreen,myAqua,myBlue,"midnightblue"))(10))
chl_colors<-c(colorRampPalette(colors=c("white",myAqua,myBlue,"midnightblue"))(10))

getCol<-function(x,log=FALSE){
  if(is.na(x)){
    ycol="white"
  }else if(is.null(x)){
    ycol="white"
  } else{
    if(log){
      y<-round(((log(x)-log(nonZeroMin))/(log(thisMax)-log(nonZeroMin))),1)*10
    }else{
      y<-ceiling(10*(x-nonZeroMin)/(thisMax-nonZeroMin))+1
    }
    ycol<-chl_colors[y]
  }
  return(ycol)
}
#read in shape file
shapeFile<-paste(DIR$'Base',"ATLANTISmodels\\inputs\\bgm\\CHAT30_LL",sep="")
sdata<-read.shapefile(shapeFile)
shape<-formatShape(shapeFile=shapeFile)
ns<-length(shape)
SpDF <- SpatialPolygonsDataFrame(shape,data.frame( z=1:ns, row.names=paste("P",seq(1,(ns)),sep="")))
labels<-seq(1,(ns))
pdf("test.pdf")
plot(shape)
LABELpOS<-polygonsLabel(shape, labels = labels, cex=.1,doPlot=FALSE)
dev.off()
labeldf<-data.frame(cbind("x"=LABELpOS[1:ns],"y"=LABELpOS[(ns+1):(2*ns)]))

t=1
pdf(paste(plotPath,"spatialTracer_",thisGroup,"_IC.pdf",sep=""),height=5)
timeData<-thisData[,t]
timeColors<-unlist(lapply(timeData,getCol))
plot(shape)
map('nzHires',add=TRUE,col="black",lwd=2)
map.axes()
for(plotB in 1:dim(labeldf)[1]){
  polygon(sdata$shp$shp[[plotB]]$points,col=timeColors[plotB])
}
dev.off()

t=1
pdf(paste(plotPath,"spatialTracer_",thisGroup,"_logIC.pdf",sep=""),height=5)
timeData<-thisData[,t]
timeColors<-unlist(lapply(timeData,getCol,log=TRUE))
plot(shape)
map('nzHires',add=TRUE,col="black",lwd=2)
map.axes()
for(plotB in 1:dim(labeldf)[1]){
  thisCol<-getCol(timeData[plotB],log=TRUE)
  # cat(plotB,", ",thisCol,"\n")
  polygon(sdata$shp$shp[[plotB]]$points,col=thisCol)
}
dev.off()

pdf(paste(plotPath,"spatialTracer_",thisGroup,".pdf",sep=""))
par(mfrow=c(3,2),mar=c(0,0,0,0),oma=c(0,0,0,0))
for(t in 1:nts){
  timeData<-thisData[,t]
  timeColors<-unlist(lapply(timeData,getCol))
  plot(shape)
  map('nzHires',add=TRUE,col="black",lwd=2)
  map.axes()
  for(plotB in 1:dim(labeldf)[1]){
    polygon(sdata$shp$shp[[plotB]]$points,col=timeColors[plotB])
  }
  mtext(paste(thisGroup,", ts",t,sep=""),side=3,adj=0,line=-1)
}
dev.off()


