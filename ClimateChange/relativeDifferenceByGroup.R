#read in all climate change runs and plot differences in response to N, Si, temperature
# by group - PPs first, then zooplankton, then inverts, then fish,... ?

this_run<-"ClimateChange"
this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")

this_out<-c("TestOutIN50yr", "TestNFromBurnin","TestSiFromBurnin", "TestTemperatureFromBurnin", "TestClimateFromBurnin") # have the climate run going through now (N, Si, temperature)
# this_out<-c("TestOutIN50yr", "TestNFromBurnin","TestSiFromBurnin", "TestTemperatureFromBurnin") 
plotDescrip<-"TestChangesFromBurnin"; thisRuns<-""

groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

plotPath<-paste(this_path,"..\\Figures\\", plotDescrip,sep="")

nruns<-length(this_out)
burnin<-rep(1,nruns) #number of years to skip in plot

runCols<-colorRampPalette(colors=c("midnightblue",myBlue,myAqua,myGold,  myOrange, "red"))(nruns)

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1865
fishingStartYear<-1865
modelStartYear<-1865

nc_list<-NULL; nts_list<-NULL; min_nts<-1e+12
for(r in 1:nruns){
  outPath<-paste(this_path,"output",this_out[r],"\\",sep="")
  nc_list[[r]]<-nc_open(paste(outPath,"output.nc",sep=""))
  thisVol<-ncvar_get(nc_list[[r]],"volume")
  thisDz<-ncvar_get(nc_list[[r]],"dz")
  nts_list[[r]]<-dim(thisVol)[3]-burnin[r] #number of timesteps
  if(nts_list[[r]]<min_nts){min_nts<-nts_list[[r]]}
}
nts_list
max_nts<-max(nts_list, na.rm=TRUE)
daysTimeStep=365; numStepsPerYear=365/daysTimeStep
xLabsTemp<-seq(0,(max_nts*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0+burnin[1]

#get all tracer names
allTracers<-sort(names(nc_list[[r]]$var))
temp<-allTracers[grep("_N",allTracers)]; tracers2plot<-temp[grep("Nums",temp,invert = TRUE)]; 
tracers2plot<-c(tracers2plot,"Oxygen","Temp","Si", "NO3")
ntracers<-length(tracers2plot)

dynBoxes<-2:26

ppCodes <- c("DF", "MA", "MB", "PL", "PS"); nPPs <- length(ppCodes)
ppNames <- str_trim(groupsDF$Name[groupsDF$Code %in% ppCodes], side="both")
ppTracers <- paste(ppNames,"_N", sep="")

## first store all tracers so can get them easily - just the dynamic boxes combined

## get differences for all PPs, as a time-series
# this_out<-c("TestOutIN50yr", "TestNFromBurnin","TestSiFromBurnin", "TestTemperatureFromBurnin", "TestClimateFromBurnin")
# noChange, N + temp, Si + temp, Temp, N + temp + Si


storePPrelchanges <- array(NA, dim=c(nruns, nPPs, nts))
for(r in 1:nruns){
  thisNC.nc <- nc_list[[r]]
  for(p in 1:nPP){
    thisTracer<-ppTracers[p]
    thisData <- ncvar_get(thisNC.nc, thisTracer)
  }
}







