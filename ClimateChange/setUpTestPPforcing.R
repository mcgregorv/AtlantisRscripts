nFutureYears<-50

testPPtimeseriesScalars <- seq(0,10,length.out=nFutureYears)

startYear<-1890
scaleYears<-startYear:(startYear+nFutureYears-1)

nts<-202
allYears<-1865:(1865+nts-1)

## read in base model as will use this to scale PP
# for the real version, will be using timeseries from ESM models
baseOut<-"BaseNEW"
basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
ThisNC.nc<-nc_open(paste(basePath,"base\\output",this_out,"\\output.nc", sep=""))
BaseNC.nc<-nc_open(paste(basePath,"base\\",baseOut,"\\output.nc", sep=""))

baseVol<-ncvar_get(BaseNC.nc, "volume"); nlayers<-dim(baseVol)[1]
dynBoxes<-2:25

thisTracer<-"Diatom_N"

diatomData<-ncvar_get(BaseNC.nc, thisTracer)
diatomBiomass<-apply(diatomData[-nlayers,dynBoxes,] * baseVol[-nlayers, dynBoxes,], 3, sum) * mg_2_tonne * X_CN
diatomMean<-apply(diatomData[-nlayers,dynBoxes,], 3, nonZeroMean)
volByTime<-apply(baseVol[-nlayers,dynBoxes,], 3, sum)
test<-diatomMean * volByTime *mg_2_tonne * X_CN

# want to rescale the non zero mean, then re-aporiton the new mean by the previous distrution wrt depth layer and box
scaleYearIndex <- allYears %in% scaleYears

newMean <- diatomMean[scaleYearIndex] * (1- testPPtimeseriesScalars/100)

sumByScaleYear<-apply(diatomData[-nlayers, dynBoxes, ], 3, sum)[scaleYearIndex]
scaledDiatomData<-0 * diatomData[,,scaleYearIndex]
for(y in 1:dim(scaledDiatomData)[3]){
  temp<-diatomData[,,scaleYearIndex][,,y]/diatomMean[scaleYearIndex][y]
  scaledDiatomData[,,y]<- temp * newMean[y]
}

test<-apply(scaledDiatomData[-nlayers, dynBoxes,] * baseVol[-nlayers, dynBoxes,scaleYearIndex], 3, sum)*mg_2_tonne *X_CN
testBase<-diatomBiomass[scaleYearIndex]

plot(testBase,type="l", ylim=c(0,max(testBase)))
points(test,type="l", col="red")

## now add the scaled PP data to the force file
thisForceFile<-paste(basePath, "\\inputs\\PPforcing_",thisTracer,".txt", sep="")
thisHeader<-paste("netcdf Chatham_PPtest {
dimensions:
	t = UNLIMITED ; // (4 currently)
	b = 30 ;
	z = 6 ;
variables:
	double t(t) ;
		t:units = \"seconds since ",startYear,"-01-01 00:00:00 +10\" ;
		t:dt = 43200. ;
	double Diatom_N(t, b, z) ;
		Diatom_N:_FillValue = -999.;
		Diatom_N:missing_value = -999.;
		Diatom_N:valid_min = 0.;
		Diatom_N:valid_max = 5000.;
		Diatom_N:units = \"mg N m-3\" ;

// global attributes:
		:title = \"trivial\" ;
		:geometry = \"CHAT30_aea.bgm\" ;
		:parameters = \"\" ;
data:\n\n",sep="")
cat(thisHeader, file=thisForceFile, append=FALSE)

#need the timesteps in seconds
scaleSeconds <- (scaleYears - startYear) * 60*60*24*365

cat("t = ",paste(scaleSeconds, collapse=", ", sep=""), ";\n", file=thisForceFile, append=TRUE)

## now add the forced data
cat(paste(thisTracer, " = \n", sep=""), file=thisForceFile, append=TRUE)
for(y in 1:nFutureYears){
  thisData<-scaledDiatomData[,,y]
  for(b in 1:dim(thisData)[2]){
    xx<-signif(thisData[,b],4)
    ## check no negatives
    xx[xx<0]<-1e-8
    cat(paste(xx,collapse=", "), file=thisForceFile, append=TRUE); 
    if(b==dim(thisData)[2] & y==nFutureYears){
      cat(";\n", file=thisForceFile, append=TRUE)
    } else{
      cat(",\n", file=thisForceFile, append=TRUE)
    }
    
    
  }
}
cat("}", file=thisForceFile, append=TRUE)


# 
# testNC.nc<-nc_open(paste(basePath,"TBGBMM\\output\\output.nc", sep=""))
# 
# testDino<-ncvar_get(testNC.nc, "Diatom_N")
# testVol<-ncvar_get(testNC.nc, "volume")
# testNlayers<-dim(testVol)[1]; testDynBoxes<-2:dim(testVol)[2]
# 
# testDinoTonnes <- apply(testDino[-testNlayers, testDynBoxes,] * testVol[-testNlayers, testDynBoxes,], 3, sum) * mg_2_tonne * X_CN
# 
# modelArea<-sum(testVol[nlayers,,1])



