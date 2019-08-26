getNitrate_WOA<-function(dataPath){
  # dataPath<-paste(DIR$'Data',"ocean_physics\\",sep="")
  # plotPath<-paste(DIR$'Figures',"Nutrients\\",sep="")
  # 
  ## this is nitrate 
  oceanNC<-nc_open(paste(dataPath,"woa13_all_n01_01.nc",sep=""))
  oceanVars<-sort(unique(names(oceanNC$var)))

  source(paste(DIR$'General functions',"fixNegLongs.R",sep=""))

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
  oceanNitrateByDepth<-data.frame(matrix(NA,nrow=dim(oceanDF)[1],ncol=ndepths)); colnames(oceanNitrateByDepth)<-oceanData$depth_bnds[2,]

  for(i in 1:dim(oceanDF)[1]){
    latIndex<-oceanData$lat_bnds[2,]==oceanDF$lat[i]
    lonIndex<-allLons==oceanDF$lon[i]
    
    thisData<-oceanData$n_an[lonIndex,latIndex,] 
    if(length(thisData)>0){
      oceanNitrateByDepth[i,]<-thisData
    }
    
  }

  oceanNitrateByDepth$lat<-oceanDF$lat; oceanNitrateByDepth$lon<-oceanDF$lon
  
  return(oceanNitrateByDepth)
}