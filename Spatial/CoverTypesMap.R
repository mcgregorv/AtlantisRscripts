## read in initial conditions and plot spatial distributions. The read in tracers and plot spatial distribution to compare
# a bit irrelevant for 'verts' as they are not doing free movement, but comparison important for zoo's and below
library(ggmap)
library(rgdal)
library(gtable)
library(maps)
library(mapdata)
library(rgeos)
source(paste(DIR$'General functions',"\\get_first_number.R",sep=""))
source(paste(DIR$'General functions',"formatShape.R",sep=""))
this_run<-"base"
this_out<-paste("BASE",sep="")

basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
plotPath<-paste(DIR$'Figures',"Spatial\\", sep="")

covers<-c("Filter_Other_Cover", "Macroalgae_Cover",   "MicroPB_Cover","reef", "flat", "soft", "canyon"); ncovers<-length(covers)

test_nts<-50

ThisNC.nc<-nc_open(paste(basePath,"base\\output",this_out,"\\output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc, "volume"); nboxes<-dim(thisVol)[2]; nts<-dim(thisVol)[3]; nlayers<-dim(thisVol)[1]

storeCovers<-array(NA,dim=c(ncovers, nboxes, nts))
for(c in 1:ncovers){
  thisCover<-covers[c]
  thisData<-ncvar_get(ThisNC.nc, thisCover)
  storeCovers[c,,]<-thisData
}

thisColRamp<-colorRampPalette(colors=c(myLightAqua, myBlue,"midnightblue"))(101)
getCol<-function(x){
  thisCol<-"white"
  y=round((x/thisMax),2)*100+1
  if(sum(x,na.rm=TRUE)>0){
    thisCol<-thisColRamp[y]
  }
  return(thisCol)
}

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


# pdf(paste(plotPath,"DensMapsTSnoAge_",this_nts,".pdf", sep=""), height=11)
par(lend=1, mfrow=c(4,2), mar=c(2,0.5,1.5,1.5))
for(c in 1:ncovers){
  thisData<-storeCovers[c,,test_nts]
  thisData[!dynBoxIndex]<-0; thisMax<-max(thisData, na.rm=TRUE)
  thisColors<-unlist(lapply(thisData, getCol))
  plot(shape)
  # map('nzHires',add=TRUE,col="black",lwd=2)
  # map.axes()
  for(b in 1:dim(labeldf)[1]){
    polygon(sdata$shp$shp[[b]]$points,col=thisColors[b])
  }
  mtext(covers[c],side=3,adj=0)
}
# dev.off()

##read in initial conditions and check spatial distribution
ThisIC.nc<-nc_open(paste(basePath,"CRAM_input_short.nc", sep=""))
par(lend=1, mfrow=c(4,2), mar=c(2,0.5,1.5,1.5))
for(c in 1:ncovers){
  thisData<-ncvar_get(ThisIC.nc, covers[c])[,1]
  thisData[!dynBoxIndex]<-0; thisMax<-max(thisData, na.rm=TRUE)
  thisColors<-unlist(lapply(thisData, getCol))
  plot(shape)
  # map('nzHires',add=TRUE,col="black",lwd=2)
  # map.axes()
  for(b in 1:dim(labeldf)[1]){
    polygon(sdata$shp$shp[[b]]$points,col=thisColors[b])
  }
  mtext(covers[c],side=3,adj=0)
}

par(lend=1, mfrow=c(2,2), mar=c(2,0.5,1.5,1.5))
for(c in 1:3){
  thisTracer<-gsub("Cover","N",covers[c])
  thisData<-ncvar_get(ThisIC.nc, thisTracer)
  if(length(dim(thisData))==2){
    thisData<-thisData[,1]
  } else{
    thisData<-apply(thisData[,,1],2,sum,na.rm=TRUE)
  }
  thisData[!dynBoxIndex]<-0; thisMax<-max(thisData, na.rm=TRUE)
  thisColors<-unlist(lapply(thisData, getCol))
  plot(shape)
  # map('nzHires',add=TRUE,col="black",lwd=2)
  # map.axes()
  for(b in 1:dim(labeldf)[1]){
    polygon(sdata$shp$shp[[b]]$points,col=thisColors[b])
  }
  mtext(thisTracer,side=3,adj=0)
}


par(lend=1, mfrow=c(2,2), mar=c(2,0.5,1.5,1.5))
for(c in 1:3){
  thisTracer<-gsub("Cover","N",covers[c])
  thisData<-ncvar_get(ThisNC.nc, thisTracer)
  if(length(dim(thisData))==2){
    thisData<-thisData[,1]
  } else{
    thisData<-apply(thisData[,,1],2,sum,na.rm=TRUE)
  }
  thisData[!dynBoxIndex]<-0; thisMax<-max(thisData, na.rm=TRUE)
  thisColors<-unlist(lapply(thisData, getCol))
  plot(shape)
  # map('nzHires',add=TRUE,col="black",lwd=2)
  # map.axes()
  for(b in 1:dim(labeldf)[1]){
    polygon(sdata$shp$shp[[b]]$points,col=thisColors[b])
  }
  mtext(thisTracer,side=3,adj=0)
}



