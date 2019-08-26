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

base_out<-"BASE"
this_out<-"FISH"
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

dynamicBoxes<-2:25
sum(baseVol[nlayers, , 1])

nts<-dim(thisVol)[3]-burnin 

storeBiomasses<-array(NA, dim=c(ng, nts))
for(g in 1:ng){
  thisName<-str_trim(groupsDF$Name[g], side="both")
  thisTracer<-paste(thisName, "_N", sep="")
  thisData<-ncvar_get(ThisNC.nc, thisTracer)
  if(length(dim(thisData))== 3){
    yy<-apply(thisData * thisVol, 3, sum) * mg_2_tonne * X_CN
  } else{
    yy <- apply(thisData * thisVol[nlayers, ,], 2, sum) * mg_2_tonne * X_CN
  }
  storeBiomasses[g,]<-yy[burnin:(nts+burnin-1)]
}
save("storeBiomasses", file=paste(DIR$'Base',"\\ATLANTISmodels\\base\\EWEbase\\CRAM_biomasses", sep=""))

# load(file=paste(DIR$'Base',"\\ATLANTISmodels\\base\\EWEbase\\CRAM_biomasses", sep=""))
#    