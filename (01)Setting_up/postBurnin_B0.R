#read in nc file, and grab the biomass' for all groups after burnin

#for those groups that have B0 estimated prior to model run (with some confidence) read these in and compare

#plot all tracers for a given box and layer
this_run<-"base"
this_out<-"BASEH3"

## all boxes
nboxes<-30
boxIndex<-seq(1,nboxes)

mg_2_tonne<-2e-8; X_CN<-5.7

burnin<-35 #number of years to skip 

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")

#read in tracers
ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume") #used in conversion mg per m^3 to biomass

#read in B0's 
thisB0df<-read.csv(paste(this_path,"..\\CRAM_B0.csv",sep=""))

groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

storeB0postBurnin<-rep(NA, ng)

for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[g], side="both")
  thisTracer<-paste(thisName,"_N", sep=""); thisData<-ncvar_get(ThisNC.nc, thisTracer)
  if(length(dim(thisData))==3){
    #then it is per m^3, so use volume
    xx<-apply(thisData*thisVol, 3, sum) * mg_2_tonne *X_CN ## convert to tonnes
  } else{
    # then it is per m^2, so use area
    xx<-apply(thisData * thisVol[nlayers,,], 2, sum) * mg_2_tonne *X_CN ## convert to tonnes
  }
  storeB0postBurnin[g]<-xx[burnin]
}

outData<-data.frame(cbind("Code"=as.character(groupsDF$Code), "Bburnin"=unlist(storeB0postBurnin)))
outData$B0<-thisB0df$B0
outData$Bburnin<-as.double(as.character(outData$Bburnin))
write.csv(outData, paste(DIR$'Tables',"B0vsBburnin.csv", sep=""), row.names = FALSE)

plot(x=outData$B0, y=outData$Bburnin)
index<-!is.na(outData$B0)
plot(outData$Bburnin[index]/outData$B0[index], pch=20, xaxt="n", ylab="Bburnin/B0", xlab="")
abline(h=1, col="red", lty=2)
par(las=2)
axis(at=seq(1,length(outData$B0[index])), labels=outData$Code[index], side=1)
abline(v=c(10,25), col="grey")


thisArea<-sum(thisVol[6,2:25,1])
