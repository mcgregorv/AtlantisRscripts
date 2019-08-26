#bring in the run, the observed data, the catch history and plot together
source(paste(DIR$'General functions',"getCIfromCV.R", sep=""))
source(paste(DIR$'General functions',"getCIfromCVs.R", sep=""))
nboxes<-30

mg_2_tonne<-2e-8; X_CN<-5.7

catchYears<-seq(1900,2014); ny<-length(catchYears)

#################################### TRACERS
#plot all tracers for a given box and layer
this_run<-"base"

# base_out<-"ORHmum4base"
# this_out<-"ORHmum4fish"

base_out<-"BASE4"
this_out<-"FISH4"
# 
# this_out<-"TEST150yrfish"
# base_out<-"TEST150yrbase"

burnin<-35 #number of years to skip in plot

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")
baseOutPath<-paste(this_path,"output",base_out,"\\", sep="")
r<-""


#read in groups csv file to get codes for fished groups
groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
catchCodes<-groupsDF$Code[groupsDF$IsFished==1]; ncg<-length(catchCodes)

catchPath<-paste(this_path,"inputs\\catch_history\\",sep="")
#read in catch ARRAY
catch_array<-read.csv(paste(DIR$'Tables',"ALLcatchHistories.csv",sep=""))

# plotPath<-paste(DIR$'Figures',"Validation\\",sep="") ## use this to overwrite the figure used in the report
plotPath<-paste(DIR$'Figures',"\\Validation\\testing\\",this_out,sep="")

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1900
fishingStartYear<-1900
modelStartYear<-1900

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
ntsteps<-dim(thisVol)[3]

BaseNC.nc<-nc_open(paste(baseOutPath,"output.nc",sep=""))
baseVol<-ncvar_get(BaseNC.nc, "volume")

nts<-dim(thisVol)[3] #number of timesteps
cat(paste(nts,"\n"))
xLabsTemp<-seq(0,(nts*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0

#get all tracer names
allTracers<-sort(names(ThisNC.nc$var))

B0data<-read.csv(paste(this_path,"..\\CRAM_B0.csv",sep=""))

###########################################
##plot them all seperately too
###################################
topCatchYears<-seq(1950,2014); biomassYears<-topCatchYears
topCatchYears<-seq(1900,2014); biomassYears<-topCatchYears

modelArea<-sum(baseVol[nlayers,,1])

storeTracersArray <- array(NA, dim=c(ntracers, nts))

for(t in 1:ntracers){
  thisTracer<-Ntracers[t]; thisName<-gsub("_N","",thisTracer); xx<-grep(thisName,groupsDF$Name)
  thisCode<-groupsDF$Code[xx]; thisNumCohorts<-groupsDF$NumCohorts[xx]
  #do biomass for all
  thisData<-ncvar_get(ThisNC.nc,thisTracer)
  if(length(dim(thisData))==3){
    xx<-apply(thisData*thisVol,3,sum)*mg_2_tonne*X_CN
  } else{
    xx<-apply(thisData * thisVol[nlayers,,], 2, sum) * mg_2_tonne * X_CN
  }
  storeTracersArray[t,]<-xx[1:nts]
 
}

# save(list=c("storeTracersArray", "Ntracers"), file=paste())