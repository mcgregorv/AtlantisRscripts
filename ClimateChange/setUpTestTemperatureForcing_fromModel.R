## model was fitted in setUpTestTemperatureForcing_fitModels.R
basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

load(paste(basePath,"ClimateChange\\modelDepth", sep="")) ## brings in modelDepth
## create 2 versions of forcing files - one repeats 1.2 SST, the other 2.6 SST - 5 year burnin + 10 year run
# read in temperature from atlantis run (or could read from input files)
nFutureYears<-10
 #increase all by 2.6%
testTempIncreases<-rep(2.6,nFutureYears); modelInt<-0.48 ## 0.48 is how much the 1999 mean was above the ROMS mean (deg celsius). Used because 1999 was used to fit modelDepth

startYear<-1870
scaleYears<-startYear:(startYear+nFutureYears-1)
dynBoxes<-2:25
nboxes<-30; nlayers <-6

thisTracer<-"temperature"

testDepths<-data.frame("Depth"=seq(100,1000,by=100)); predTempChanges<-(predict(modelDepth, newdata=testDepths)/modelInt ) * testTempIncreases[1]
par(lend=1)
plot(x=testDepths$Depth, y=predTempChanges, type="h", lwd=5, col=myBlue, ylab="Temperature increase", xlab="Depth (m)")

## read in temperature initial conditions for the year with highest temperature
# for the real version, will be using timeseries from ESM models, and in between that may need to sample all base years 1996:2004
inputROMSPath<-paste(basePath,"inputs\\ROMS\\",sep=""); 
baseYear<-1999
tempNCfile<-paste(inputROMSPath,"Chatham_temp_year",baseYear,".nc",sep="")
TempIC<-nc_open(tempNCfile)
temperatureIC <- ncvar_get(TempIC, thisTracer)

## have model for fitting relationship between delta SST / delta temperature using depth and month
## delta SST is going to be testPPtimeseriesScalars for this year
# need depth by cell
load(paste(basePath,"ClimateChange\\inputDepthsByCell", sep=""))
# need months
inputYears <- sort(rep(scaleYears, 365*2)); 
monthDays<-c(31,28,31,30,31,30,31,31,30,31,30,31)
monthYears <- c(rep(1, monthDays[1]), rep(2, monthDays[2]), rep(3, monthDays[3]), rep(4, monthDays[4]), rep(5, monthDays[5]), rep(6, monthDays[6])
                , rep(7, monthDays[7]), rep(8, monthDays[8]), rep(9, monthDays[9]), rep(10, monthDays[10]), rep(11, monthDays[11]), rep(12, monthDays[12]))
inputMonths <- rep(sort(rep(monthYears, 2)), nFutureYears)

# set up an array to store scaled temperatures - need dimensions time, box, layer
# time needs to have day and night, for all years wer're forcing (730* nFutureYears)
scaledTempDF<-array(NA, dim=c(nlayers, nboxes, 730*nFutureYears))
## the month and depth parts are the same for each year - so can define here
tempMonths <- sort(rep(monthYears,2)) # just times by 2 to get night and day 
## delta SST is just testTempScalars

baseSST <- apply(temperatureIC, c(2,3), max, na.rm=TRUE)
# have a look at model increase temp / base temp at depth
ratioTemp<-temperatureIC * 0
for(b in dynBoxes){
  for(l in 1:nlayers){
    thisDepth <- inputDepthsByCell[l,b]
    if(!is.na(thisDepth)){
      tempdf<-data.frame("Depth"=thisDepth); tempPred
    }
  }
}


for(y in 1:nFutureYears){
  thisYear <- scaleYears[y]
  cat(thisYear,"--")
  # index which timesteps are for this year
  thisYearIndex <- inputYears == thisYear
  this_ts<-seq(1,length(thisYearIndex))[thisYearIndex]
  thisNewSST <- baseSST + testTempIncreases[y]
  # storeThisYearsTemps <- 0 * temperatureIC
  # storeDeltaTemp <- array(NA, dim=c(nlayers, nboxes, 730)); storeRelTemp <- storeDeltaTemp
  for(t in this_ts){
    base_t_Index <- (t %% 730)+1
    thisMonth <- tempMonths[base_t_Index]
    for(b in 1:nboxes){
      if(b %in% dynBoxes){
        thisSurfaceDepth <- min(inputDepthsByCell[,b], na.rm=TRUE)
        for(l in 1:(nlayers-1)){ ## leave sediment as is
          # get the base temperature for this cell and multiply by deltaTemp to get newTemp
          thisBaseTemp<-temperatureIC[l,b,base_t_Index]
          if(!is.na(thisBaseTemp)){
            thisDepth <- inputDepthsByCell[l,b]
            thisInt <- (predict(modelDepth, newdata = data.frame("Depth"=thisDepth))/modelInt) * testTempIncreases[y]
            ## if this is the most surface layer for this box, just scale the temperature by thisScalar - don't need to fit model
            if(thisDepth == thisSurfaceDepth){
              thisNewTemp <- thisBaseTemp + testTempIncreases[y]
            } else{
              thisNewTemp<-thisBaseTemp  + thisInt
            }
            scaledTempDF[l,b,t]<-thisNewTemp
            # storeThisYearsTemps[l,b,t]<-thisNewTemp
          }

        }
      } else{
        # just set at 15 for boundary
        thisNewTemp <-15
        scaledTempDF[,b,t]<-thisNewTemp
        # storeThisYearsTemps[,b,t]<-thisNewTemp
      }
      # storeThisYearsTemps[nlayers,b,t]<-temperatureIC[nlayers,b,t]
    }
  }
}
testMeanByTS <- apply(scaledTempDF[,dynBoxes,], 3, nonZeroMean)
baseMeanByTS <- apply(temperatureIC[,dynBoxes,], 3, nonZeroMean)
thisYmax<-max(c(max(testMeanByTS, na.rm=TRUE), max(baseMeanByTS, na.rm=TRUE)))
plot(testMeanByTS, type="l", ylim=c(8,thisYmax))
points(baseMeanByTS, type="l", col="red")



## now add the scaled temperature data to the force file
thisForceFile<-paste(basePath, "\\inputs\\ClimateForcing_",thisTracer,"TEST26.txt", sep="")
thisHeader<-paste("netcdf Chatham_Temptest {
                  dimensions:
                  t = UNLIMITED ; // (4 currently)
                  b = 30 ;
                  z = 6 ;
                  variables:
                  double t(t) ;
                  t:units = \"seconds since ",startYear,"-01-01 00:00:00 +10\" ;
                  t:dt = 43200. ;
                  double temperature(t, b, z) ;
                  temperature:_FillValue = -999.;
                  temperature:missing_value = -999.;
                  temperature:valid_min = 0.;
                  temperature:valid_max = 300.;
                  temperature:units = \"degrees celsius\" ;
                  
                  // global attributes:
                  :title = \"trivial\" ;
                  :geometry = \"CHAT30_aea.bgm\" ;
                  :parameters = \"\" ;
                  data:\n\n",sep="")
cat(thisHeader, file=thisForceFile, append=FALSE)

#need the timesteps in seconds
ts_int_seconds <- 60*60*12; startSeconds <- startYear * 60*60*24*365; endSeconds <- (max(scaleYears)) * 60*60*24*365
scaleSeconds <- seq(startSeconds, endSeconds, by=ts_int_seconds)

cat("t = ",paste(scaleSeconds, collapse=", ", sep=""), ";\n", file=thisForceFile, append=TRUE)

## now add the forced data
cat(paste(thisTracer, " = \n", sep=""), file=thisForceFile, append=TRUE)
for(t in 1:length(scaleSeconds)){
  thisData<-scaledTempDF[,,t]
  for(b in 1:dim(thisData)[2]){
    xx<-signif(thisData[,b],4)
    ## check no negatives
    xx[xx<0]<-1e-8
    xx[is.na(xx)]<-0
    xx[nlayers]<-min(xx[xx>0])
    cat(paste(xx,collapse=", "), file=thisForceFile, append=TRUE); 
    if(b==dim(thisData)[2] & t==length(scaleSeconds)){
      cat(";\n", file=thisForceFile, append=TRUE)
    } else{
      cat(",\n", file=thisForceFile, append=TRUE)
    }
    
    
  }
}
cat("}", file=thisForceFile, append=TRUE)





# cat(paste(thisTracer, " = \n", sep=""), file=thisForceFile, append=TRUE)
# for(y in 1:nFutureYears){
#   thisData<-scaledTempDF[,,y]
#   for(b in 1:dim(thisData)[2]){
#     xx<-signif(thisData[,b],4)
#     ## check no negatives
#     xx[xx<0]<-1e-8
#     xx[is.na(xx)]<-0
#     xx[nlayers]<-min(xx[xx>0])
#     cat(paste(xx,collapse=", "), file=thisForceFile, append=TRUE); 
#     if(b==dim(thisData)[2] & y==nFutureYears){
#       cat(";\n", file=thisForceFile, append=TRUE)
#     } else{
#       cat(",\n", file=thisForceFile, append=TRUE)
#     }
#     
#     
#   }
# }
# cat("}", file=thisForceFile, append=TRUE)
# 
# 
# 
