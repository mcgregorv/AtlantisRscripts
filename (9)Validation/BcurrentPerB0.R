#read in fishing run and create table with Bcurrent (%B0) so can compare with stock assessment for those that have one

#bring in the run, the observed data, the catch history and plot together
source(paste(DIR$'General functions',"getCIfromCVs.R", sep=""))
nboxes<-30

mg_2_tonne<-2e-8; X_CN<-5.7

catchYears<-seq(1900,2014); ny<-length(catchYears)

#################################### TRACERS
#plot all tracers for a given box and layer
this_run<-"base"
this_out<-"TEST150yrfish_HOK"
# 
# this_out<-"TEST150yrfish"
# base_out<-"TEST150yrbase"

burnin<-35 #number of years to skip in plot

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")

#read in groups csv file to get codes for fished groups
groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
catchCodes<-groupsDF$Code[groupsDF$IsFished==1]; ncg<-length(catchCodes)

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1900
fishingStartYear<-1900
modelStartYear<-1900

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
ntsteps<-dim(thisVol)[3]

thisProp<-0.9
storeBcursPercB0<-rep(NA,ng)
for(g in 1:ng){
    thisCode<-groupsDF$Code[g]; isFished<-groupsDF$IsFished[g]; thisName<-str_trim(groupsDF$Name[g], side="both")
  if(isFished==1){
    thisVar<-paste(thisName,"_N", sep="")
    temp<-ncvar_get(ThisNC.nc, thisVar)
    xx<-apply(temp*thisVol,3,sum)*mg_2_tonne * X_CN
    thisB0<-xx[burnin]
    storeBcursPercB0[g]<-round(xx[length(xx)]/thisB0,2)*100
  }
}

