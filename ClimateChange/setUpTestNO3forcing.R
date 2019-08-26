## est change of -1 micromol per m^3

## read in outputs from the base run to check what No3 we have in the base model
decreasePerc <- 0.18 # based on 1 mmol decrease of NO3 per m^3 from average of 5.31 mmol m-3
startYear<-1900

nboxes<-30

mg_2_tonne<-2e-8; X_CN<-5.7

this_run<-"base"
this_run<-"ClimateChange"

base_out<-"TestOutIN5yrNewInput"
# this_out<-"FISH4"


this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep=""); basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="");
baseOutPath<-paste(this_path,"output",base_out,"\\", sep="")


#read in groups csv file to get codes for fished groups
groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

daysTimeStep<-1
numStepsPerYear<-365/daysTimeStep
year0<-1900
fishingStartYear<-1900
modelStartYear<-1900

BaseNC.nc<-nc_open(paste(baseOutPath,"output.nc",sep=""))
baseVol<-ncvar_get(BaseNC.nc, "volume")
thisVol<-ncvar_get(BaseNC.nc,"volume")
thisDz<-ncvar_get(BaseNC.nc,"dz")
ntsteps<-dim(thisVol)[3]

nts<-dim(baseVol)[3] #number of timesteps
cat(paste(nts,"\n"))
xLabsTemp<-seq(0,(nts*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0

#get all tracer names
allTracers<-sort(names(BaseNC.nc$var))

## base NO3 data - use both NO3 and NH3 (which I think is actually NH4)
# burnin<-75
dynBoxes<-2:26
baseNO3<-ncvar_get(BaseNC.nc,"NO3") ; baseNH3<- ncvar_get(BaseNC.nc, "NH3");
baseN <- baseNO3 + baseNH3
meanN <- apply(baseN[,dynBoxes,], 3, nonZeroMean)
thisMax<-max(c(meanN), na.rm=TRUE); thisMin<-min(c(meanN), na.rm=TRUE)
par(mfrow=c(2,1), mar=c(4,4,1,1))
plot(meanN,type="l", ylim=c(thisMin,thisMax), xaxt="n"); 
axis(at=xLabsAt, labels = xLabs, side=1)




# grab a snapshot year - the second year after burnin seems fine (73 timesteps per year, so 73*2 : 73*3-1)
ssTimeSteps <- 365:729
NO3Snapshot <- baseNO3[,,ssTimeSteps]; NH3Snapshot <- baseNH3[,,ssTimeSteps]
scaledNO3 <- NO3Snapshot * (1-decreasePerc); scaledNH3 <- NH3Snapshot * (1-decreasePerc)
yy<-apply(scaledNO3[,dynBoxes,],3,nonZeroMean) + apply(scaledNH3[,dynBoxes,], 3, nonZeroMean)

par(new=TRUE)
plot(x=ssTimeSteps, y=yy, xlim=c(0,dim(baseNO3)[3]), type="l", col=myLightBlue, yaxt="n")
thisMax<-max(c(yy,meanNO3[ssTimeSteps])); thisMin<-min(c(yy, meanNO3))
plot(yy, type="l", col=myLightBlue, ylim=c(thisMin,thisMax))
points(meanNO3[ssTimeSteps], type="l", col=myGreen)

plot(yy/meanN[ssTimeSteps], type="l")


#####################################
## check out silicate too
siDecrease <- 0.056
baseSi<-ncvar_get(BaseNC.nc,"Si");
meanSi <- apply(baseSi[,dynBoxes,], 3, nonZeroMean)
# use same timesteps
SiSnapshot <- baseSi[,,ssTimeSteps]
scaledSi <- SiSnapshot * (1-siDecrease)
yy<-apply(SiSnapshot[,dynBoxes,],3,nonZeroMean); yyScaled <- apply(scaledSi[,dynBoxes,],3,nonZeroMean);
thisMax<-max(c(meanSi,yy,yyScaled), na.rm=TRUE); thisMin<-min(c(meanSi, yy, yyScaled), na.rm=TRUE)
plot(meanSi,type="l", ylim=c(thisMin,thisMax), xaxt="n");
axis(at=xLabsAt, labels = xLabs, side=1)

par(new=TRUE)
plot(x=ssTimeSteps, y=yy, xlim=c(0,dim(baseSi)[3]), type="l", col=myLightBlue, ylim=c(thisMin, thisMax))
points(x=ssTimeSteps, y=yyScaled, type="l", col=myOrange, lty=2)

thisMin<-min(c(yy,yyScaled), na.rm=TRUE); thisMax<-max(c(yy,yyScaled), na.rm=TRUE)
plot(rep(yy,3), type="l", col=myLightBlue, ylim=c(thisMin, thisMax)); points(rep(yyScaled,3), col=myOrange, lty=2, type="l")

#####################################

thisTracer<-"NO3"; tracersScaledData<-scaledNO3

# thisTracer<-"NH3"; tracersScaledData<-scaledNH3

# thisTracer <- "Si"; tracersScaledData<-scaledSi



# time_int<-12*60*60
# scaleSeconds <- seq(0,364.5* 60*60*24, by=time_int) 

time_int<-24*60*60
scaleSeconds <- seq(time_int*1,365*10* time_int, by=time_int) 

## now add the scaled PP data to the force file
thisForceFile<-paste(basePath, "\\inputs\\PPforcing_",thisTracer,".txt", sep="")
# thisForceFile<-paste(basePath, "\\inputs\\PPforcing_",thisTracer,"TESTtime.txt", sep="")
thisHeader<-paste("netcdf PPforcing_",thisTracer," {
                  dimensions:
                  t = UNLIMITED ; // (4 currently)
                  b = 30 ;
                  z = 6 ;
                  variables:
                  double t(t) ;
                  t:units = \"seconds since ",startYear,"-01-01 00:00:00 +10\" ;
                  t:dt = ",time_int,". ;
                  double ",thisTracer,"(t, b, z) ;
                  ",thisTracer,":_FillValue = -999.;
                  ",thisTracer,":missing_value = -999.;
                  ",thisTracer,":valid_min = 0.;
                  ",thisTracer,":valid_max = 5000.;
                  ",thisTracer,":units = \"mg N m-3\" ;
                  
                  // global attributes:
                  :title = \"trivial\" ;
                  :geometry = \"CHAT30_aea.bgm\" ;
                  :parameters = \"\" ;
                  data:\n\n",sep="")
cat(thisHeader, file=thisForceFile, append=FALSE)

#need the timesteps in seconds
## just do one year
cat("t = ",paste(scaleSeconds, collapse=", ", sep=""), ";\n", file=thisForceFile, append=TRUE)

## now add the forced data
cat(paste(thisTracer, " = ", sep=""), file=thisForceFile, append=TRUE)
for(t in 1:length(scaleSeconds)){
  #need to index day of year - as reusing the scaled tracer for one year
  dayOfYear <- (scaleSeconds[t]/time_int) %% 365 +1
  # thisData<-tracersScaledData[,,t]
  thisData<-tracersScaledData[,,dayOfYear]
  # if(scaleSeconds[t]<startLow){thisData <- apply(thisData, c(1,2), FUN=function(x){gsub("0", "_", x)})}
  for(b in 1:dim(thisData)[2]){
    if(is.character(thisData[,b])==TRUE){
      xx<-thisData[,b]
    } else{
      xx<-signif(thisData[,b],4)
      ## check no negatives
      xx[xx<0]<-1e-8
    }
    ## need to swap the zeros and the WC non-zeros
    wcvars<-xx[1:(nlayers-1)]; indexZeros<-wcvars==0
    newWCvars <- c(wcvars[!indexZeros], wcvars[indexZeros])
    xx<-c(newWCvars, xx[nlayers])
    cat(paste(xx,collapse=", "), file=thisForceFile, append=TRUE); 
    if(b==dim(thisData)[2] & t==length(scaleSeconds)){
      cat(";\n", file=thisForceFile, append=TRUE)
    } else{
      cat(",\n", file=thisForceFile, append=TRUE)
    }
    
    
  }
}
cat("}", file=thisForceFile, append=TRUE)



