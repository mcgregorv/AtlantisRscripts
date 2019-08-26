getOxygen_WOA<-function(dataPath){
  # dataPath<-paste(DIR$'Data',"ocean_physics\\",sep="")

  source(paste(DIR$'General functions',"fixNegLongs.R",sep=""))

  oceanNC<-nc_open(paste(dataPath,"woa13_all_o00_01.nc",sep=""))
  oceanVars<-sort(unique(names(oceanNC$var)))
  
  plotPath<-paste(DIR$'Figures',"Nutrients\\",sep="")
  
  oceanData<-NULL
  for(v in 1:length(oceanVars)){
    thisVar<-oceanVars[v]
    oceanData[[thisVar]]<-ncvar_get(oceanNC,thisVar)
  }

  ndepths<-length(oceanData$depth_bnds[2,])

  latIndex<-oceanData$lat_bnds[2,]<=(-42) & oceanData$lat_bnds[2,]>(-48)

  ##need to edit lons to not have negatives (1-360)
  allLons<-unlist(lapply(oceanData$lon_bnds[2,],fixNegLongs))
  lonIndex<-allLons>=172 & allLons<=190
  
  lats<-sort(unique(oceanData$lat_bnds[2,latIndex]))
  lons<-sort(allLons[lonIndex])
  
  oceanDF<-data.frame(cbind(sort(rep(lats,length(lons))),rep(sort(lons),length(lats))))
  colnames(oceanDF)<-c("lat","lon")
  oceanOxygenByDepth<-data.frame(matrix(NA,nrow=dim(oceanDF)[1],ncol=ndepths)); colnames(oceanOxygenByDepth)<-oceanData$depth_bnds[2,]

  for(i in 1:dim(oceanDF)[1]){
    latIndex<-oceanData$lat_bnds[2,]==oceanDF$lat[i]
    lonIndex<-allLons==oceanDF$lon[i]
    
    thisData<-oceanData$o_an[lonIndex,latIndex,] 
    if(length(thisData)>0){
      oceanOxygenByDepth[i,]<-thisData
    }
    
  }
  oceanOxygenByDepth$lat<-oceanDF$lat; oceanOxygenByDepth$lon<-oceanDF$lon
  return(oceanOxygenByDepth)
}