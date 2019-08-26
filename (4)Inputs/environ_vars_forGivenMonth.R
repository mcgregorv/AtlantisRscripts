##read in the .nc file, pull out data for given month and put it into a .csv file
##temp
basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

#time is in seconds from 1996. There is one value per 12 hour timestep
tempNC.nc<-nc_open(paste(basePath,"inputs\\ROMS\\Chatham30_tempAll.nc",sep=""))

tempVars<-sort(names(tempNC.nc$var))
tempDF<-ncvar_get(tempNC.nc,"temperature")

#just want the SST, so take max within each box
SSTdf<-apply(tempDF,c(2,3),max,na.rm=TRUE)

timeSteps<-seq(0,ceiling(dim(SSTdf)[2]/2),by=0.5)[1:dim(SSTdf)[2]]

##which ones are in January? 
timeYears<-unlist(lapply(timeSteps,FUN=function(x){trunc(x/365)}))+1996
thisYears<-sort(unique(timeYears))

meanJanSST<-data.frame(matrix(NA,ncol=length(thisYears),nrow=30)); colnames(meanJanSST)<-thisYears

for(y in 1:length(thisYears)){
  thisYear<-thisYears[y]
  thisData<-SSTdf[,timeYears==thisYear][,1:62]
  meanByBox<-apply(thisData,1,mean,na.rm=TRUE)
  meanJanSST[,y]<-meanByBox
}

write.csv(meanJanSST,file=paste(DIR$'Tables',"MeanSSTJanuaryByBoxNumber.csv",sep=""),row.names = FALSE)
