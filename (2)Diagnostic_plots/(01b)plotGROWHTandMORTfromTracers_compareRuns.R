#plot growth and mortality curves based on tracers for a given group
this_run<-"base"

plotPath<-paste(this_path,"..\\Figures\\BALnumbers_",sep="")

this_out<-c("outputXXX_mQA1", "outputXXX_mQA2", "outputXXX_mQA3","outputXXX_mQA4")
# this_out<-c("ShortLonger","ShortLongerCompare"); 
this_out<-c("TEST150yrfish_HOK", "FISHXXX_mLFC1", "FISHXXX_mLFC2", "FISHXXX_mLFC3", "FISHXXX_mLFC4", "FISHXXX_mLFC5");


nruns<-length(this_out)
burnin<-1 #number of years to skip in plot

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1900
fishingStartYear<-1900
modelStartYear<-1900

nc_list<-NULL; nts_list<-NULL; min_nts<-1e+12; max_nts<-0
for(r in 1:nruns){
  outPath<-paste(this_path,"output",this_out[r],"\\",sep="")
  
  nc_list[[r]]<-nc_open(paste(outPath,"output.nc",sep=""))
  thisVol<-ncvar_get(nc_list[[r]],"volume")
  thisDz<-ncvar_get(nc_list[[r]],"dz")
  nts_list[[r]]<-dim(thisVol)[3]-burnin #number of timesteps
  if(nts_list[[r]]<min_nts){min_nts<-nts_list[[r]]}
  if(nts_list[[r]]>max_nts){max_nts<-nts_list[[r]]}
}

xLabsTemp<-seq(0,(min_nts*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0

#get all tracer names
allTracers<-sort(names(nc_list[[r]]$var))
tracers2plot<-allTracers[grep("Nums",allTracers,invert = FALSE)]; ntracers<-length(tracers2plot)

runCols<-colorRampPalette(colors=c(myGreen,myAqua,myBlue,"midnightblue"))(nruns)

##loop through runs and grab out weights of individuals by cohort and timestep and numbers by cohort and timestep and store
thisCode<-"BAL"; thisNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisCode]; thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode],side="both")
numbers_array<-array(NA,dim=c(nruns, thisNumCohorts, max_nts)); weights_array<-numbers_array

for(r in 1:nruns){
  for(c in 1:thisNumCohorts){
    thisVar<-paste(thisName,c,"_ResN",sep=""); thisTemp<-ncvar_get(nc_list[[r]],thisVar); thisRN<-apply(thisTemp,3,nonZeroMean)
    thisVar<-paste(thisName,c,"_StructN",sep=""); thisTemp<-ncvar_get(nc_list[[r]],thisVar); thisSN<-apply(thisTemp,3,nonZeroMean)
    weights_array[r,c,1:nts_list[[r]]]<-(thisRN+thisSN)[1:nts_list[[r]]]
    ##numbers
    thisVar<-paste(thisName,c,"_Nums",sep=""); thisTemp<-ncvar_get(nc_list[[r]], thisVar); thisNums<-apply(thisTemp,3,sum)
    numbers_array[r,c,1:nts_list[[r]]]<-thisNums[1:nts_list[[r]]]
  }
}

weightMax<-max(weights_array,na.rm=TRUE); numMax<-max(numbers_array,na.rm=TRUE)


legendText<-paste(this_out,unlist(nts_list))


## now plot the curves at timestep 1, min_nts and run max nts
plotsFile<-paste(plotPath,"GrowthAndMortCurves_",thisCode,".jpg",sep="")
jpeg(plotsFile,quality=300, width=800, height = 400)
par(mfrow=c(1,2), mar=c(4,4,1,1), oma=c(1,1,1,1))
plot(x=seq(1,thisNumCohorts),y=seq(0,weightMax,length.out=thisNumCohorts),type="n",ylab="Weight per individual", xlab="Age class")
for(r in 1:nruns){
  for(c in 1:thisNumCohorts){
    points(x=rep(c,length(weights_array[r,c,])), y=weights_array[r,c,],col=runCols[r],pch=r,lwd=2)
  }
  points(x=seq(1,thisNumCohorts), y=weights_array[r,,nts_list[[r]]], type="l", lty=1, col=runCols[r])
}
points(x=seq(1,thisNumCohorts),y=weights_array[1,,1],type="l",col="red",lwd=2)

plot(x=seq(1,thisNumCohorts),y=seq(0,numMax,length.out=thisNumCohorts),type="n",ylab="Numbers", xlab="Age class")
for(r in 1:nruns){
  for(c in 1:thisNumCohorts){
    points(x=rep(c,length(numbers_array[r,c,])), y=numbers_array[r,c,],col=runCols[r],pch=r,lwd=2)
  }
  points(x=seq(1,thisNumCohorts), y=numbers_array[r,,nts_list[[r]]], type="l", lty=1, col=runCols[r])
}
points(x=seq(1,thisNumCohorts),y=numbers_array[1,,1],type="l",col="red",lwd=2)
legend(legend=legendText, col=runCols, pch=seq(1,nruns), x="topright", bty="n")
dev.off()


