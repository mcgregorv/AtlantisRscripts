thisPath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

version<-"C"
runPath<-paste(thisPath,"eddySens\\",sep="")

runs<-paste("outputEDDYCTESTING_",1:2, sep="")
nruns<-length(runs)

tracersOfInterest<-paste(groupsDF$Name[groupsDF$Code %in% c("PL", "PS", "MA", "MB", "DF")], "_N", sep="")
ntt<-length(tracersOfInterest)

dynBoxes<-2:25

storeTOI<-array(NA, dim=c(nruns, ntt, 151))
for(r in 1:nruns){
  thisRunPath<-paste(runPath,runs[r],"\\", sep="")
  runNC.nc<-nc_open(paste(thisRunPath,"output.nc", sep=""))
  runVol<-ncvar_get(runNC.nc, "volume")
  for(t in 1:ntt){
    thisTracer<-tracersOfInterest[t]; thisData<-ncvar_get(runNC.nc, thisTracer)
    if(length(dim(thisData))==3){
      xx<-apply(thisData[-6,dynBoxes,] * runVol[-6,dynBoxes,], 3, sum, na.rm=TRUE) * mg_2_tonne * X_CN
    } else{
      xx<-apply(thisData[dynBoxes,] * runVol[nlayers,dynBoxes,], 2, sum, na.rm=TRUE) * mg_2_tonne * X_CN
    }
    storeTOI[r,t,1:(length(xx))]<-xx
  }
}
tracersOfInterest

plot(storeTOI[1,5,])
points(storeTOI[2,5,],type="l")

t=5; thisTracer<-tracersOfInterest[t]
r=1
thisRunPath<-paste(runPath,runs[r],"\\", sep="")
runNC.nc<-nc_open(paste(thisRunPath,"output.nc", sep=""))
runVol1<-ncvar_get(runNC.nc, "volume")
data1<-ncvar_get(runNC.nc, thisTracer)
xx1<-apply(data1[-6,dynBoxes,] * runVol1[-6,dynBoxes,], 3, sum, na.rm=TRUE) * mg_2_tonne * X_CN

r=2
thisRunPath<-paste(runPath,runs[r],"\\", sep="")
runNC.nc<-nc_open(paste(thisRunPath,"output.nc", sep=""))
runVol2<-ncvar_get(runNC.nc, "volume")
data2<-ncvar_get(runNC.nc, thisTracer)
xx2<-apply(data2[-6,dynBoxes,] * runVol2[-6,dynBoxes,], 3, sum, na.rm=TRUE) * mg_2_tonne * X_CN


plot(xx1, type="l")
points(xx2, pch=20)
xx2/xx1[1:length(xx2)]

xx2<-apply(data2 * runVol2, 3, sum, na.rm=TRUE) * mg_2_tonne * X_CN
xx1<-apply(data1 * runVol1, 3, sum, na.rm=TRUE) * mg_2_tonne * X_CN
(xx1 - xx2)/xx1

storeTOI[,5,]

logLines<-readLines(paste(runPath,"outputEDDYCTESTING_2\\log.txt", sep=""))
x<-grep("PS sp_grow:", logLines); thisLines<-logLines[x]; timeLines<-logLines[(x-1)]
getTime<-function(l){
  time<-as.double(unlist(str_split(l," "))[2])
  return(time)
}
thisTimes<-unlist(lapply(timeLines, getTime))
getBiom<-function(l){
  biom<-as.double(unlist(str_split(gsub(",","",l)," "))[5])
  return(biom)
}
thisBiols<-unlist(lapply(thisLines, getBiom))

r=2
runNC.nc<-nc_open(paste(thisRunPath,"output.nc", sep=""))
runVol<-ncvar_get(runNC.nc, "volume")
storeTOI<-array(NA, dim=c(ntt,dim(runVol)[3]))
for(t in 1:ntt){
  thisTracer<-tracersOfInterest[t]; thisData<-ncvar_get(runNC.nc, thisTracer)
  if(length(dim(thisData))==3){
    xx<-apply(thisData * runVol, 3, sum, na.rm=TRUE) * mg_2_tonne * X_CN
  } else{
    xx<-apply(thisData * runVol[nlayers,,], 2, sum, na.rm=TRUE) * mg_2_tonne * X_CN
  }
  storeTOI[t,1:(length(xx))]<-xx
}
plot(storeTOI[5,])

dynBoxes<-seq(2,25)

test<-apply(thisData, 3, nonZeroMean)

test1<-apply(thisData[,dynBoxes,], 3, nonZeroMean)

















