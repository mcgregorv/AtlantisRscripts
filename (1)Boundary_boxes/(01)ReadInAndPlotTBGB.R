library(ggmap)
library(rgdal)
library(gtable)
library(maps)
library(mapdata)
library(rgeos)
source(paste(DIR$'General functions',"\\get_first_number.R",sep=""))
source(paste(DIR$'General functions',"formatShape.R",sep=""))

thisPath<-paste(DIR$'Base',"\\TBGB\\TBGBSpatialPolygons\\",sep="")

plotPath<-paste(DIR$'Base',"TBGB\\Figures\\",sep="")
shapeFile<-paste(thisPath,"P0to26_depth",sep="")

sdata<-read.shapefile(shapeFile)
shape<-formatShape(shapeFile=shapeFile)
ns<-length(shape)
SpDF <- SpatialPolygonsDataFrame(shape,data.frame( z=1:ns, row.names=paste("P",seq(1,(ns)),sep="")))
labels<-seq(1,(ns))
plot(shape)
LABELpOS<-polygonsLabel(shape, labels = labels, cex=.1,doPlot=FALSE)
labeldf<-data.frame(cbind("x"=LABELpOS[1:ns],"y"=LABELpOS[(ns+1):(2*ns)]))


plot(shape)
map('nzHires',add=TRUE,col="black",lwd=2)
map.axes()
for(plotB in 1:dim(labeldf)[1]){
  Blab <- plotB
  if(Blab==dim(labeldf)[1]){Blab<-0}
  test<-sdata$shp$shp[[plotB]]$points
  thisX<-mean(test$'X'); thisY<-mean(test$'Y')
  text(Blab, x=thisX, y=thisY)
  # polygon(sdata$shp$shp[[plotB]]$points,col=thisCol)
}

