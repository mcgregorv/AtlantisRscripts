## Vidette McGregor NIWA ##
formatShape<-function(shapeFile,projection=NULL){
  shape<-read.shapefile(shapeFile)
  
  shapeList<-list()
  nshape<-length(shape$shp$shp)
  for(s in 1:nshape){
    thisP<-Polygon(shape$shp$shp[[s]]$points)
    shapeList[[s]]<-Polygons(list(thisP),paste("P",s,sep=""))
  }
  shapeP<-SpatialPolygons(shapeList)
  # attr<-data.frame(attr1=1:nshape,row.names=paste("P",seq(1,nshape),sep=""))
  # spatialPDF<-SpatialPolygonsDataFrame(shapeP,attr)
  # if(is.null(projection)){
  #   proj4string(spatialPDF) <- CRS("+init=EPSG:4326")
  # }else{
  #   proj4string(spatialPDF) <- CRS(projection)
  # }
  return(shapeP)
}

