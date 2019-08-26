#plot all tracers for a given box and layer
this_run<-"base"

this_out<-"EatMoreAS5" #I've only used a current run to get the number of boxes and layers (nlayers and nboxes below), 
# but you could enter these manually

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")

#this is the file you are about to create
mortForcingFile<-paste(this_path,"mortForcingTEST.txt",sep="")

#you need to edit this - use the same one that is being used in the initial conditions
thisBgmFile<-"CHAT30_aea.bgm"

# not sure if it has to be, but makes sense to set this to the same as model intitial conditions
modelStartYear<-1900

#define timesteps for mortality changes
timeSteps<-c(5,10,50,150,300); nts<-length(timeSteps) #time steps are in days. It uses the first mort scalar up until and including the first timestep
# then interpolates between each of the next pairs of timesteps, then if you set externalBiologyForcingFile_rewind 1 in your force.prm file, it 
# will reuse the scalings in order, and if set to 0 it remains at the final value.

#set scalars for mortality. most likely will want first value to be 1 since 
mortScalars<-c(1,5,0.2,1,10)  ##needs to be same length as timesteps - will be applied in same order as timesteps

#just getting numLayers and numBoxes - can replace and do manually if prefer
ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
numLayers<-dim(thisVol)[1]; numBoxes<-dim(thisVol)[2]

##number of seconds in a day - have set up the file in seconds, so need this to convert timestep days to second
secPerDay<-24*60*60

#replace this with "SCA"
thisGroupCode<-"HOK"

#this is to hold arrays of mortality scalars in required dimensions. 
# I'm just using same value for each cell (which is what you'd probably want to do to start with), 
# but you could use different values for different boxes/layers
baseMortArray<-array(1,dim=c(numLayers,numBoxes))
#populate the list of arrays (one array for each timestep for which you are changing the mortality scalar)
mortList<-NULL
for(t in 1:nts){
  mortList[[t]]<-baseMortArray*mortScalars[t]
}

##start the moratlity forcing file
cat(paste("netcdf externalMort {
dimensions:
	t = UNLIMITED ; // (3 currently)
	b = ",numBoxes," ;
	z = ",numLayers," ;
variables:
	double t(t) ;
		t:units = \"seconds since ",modelStartYear,"-01-01 00:00:00 +10\" ;
		t:dt = 86400. ;
	double ",thisGroupCode,"mort_ff(t, b, z) ;
		",thisGroupCode,"mort_ff:_FillValue = 1. ;\n",sep=""),file=mortForcingFile,append=FALSE)


##add this bit
cat(paste("// global attributes:
		:title = \"SpatialForcing\" ;
		:geometry = \"",thisBgmFile,"\" ;
		:parameters = \"\" ;\n",sep=""),file=mortForcingFile,append=TRUE)

##start the data section and put in timesteps
cat("data:
  
  t = ",paste(timeSteps*secPerDay,collapse=", ")," ;\n", file=mortForcingFile,append=TRUE)

cat(paste(thisGroupCode,"mort_ff = \n",sep=""),file=mortForcingFile,append=TRUE)
for(t in 1:nts){
  for(b in 1:nboxes){
    for(l in 1:nlayers){
      if(l==numLayers & b==numBoxes & t==nts){ #if it is the last row, need to end with ;
        # and since only doing one group, it is also the end of the file, so end with }
        cat(paste(mortList[[t]][l,b],";}",sep=""),file=mortForcingFile,append=TRUE)
      } else{
        cat(paste(mortList[[t]][l,b],", ",sep=""),file=mortForcingFile,append=TRUE)
      }
      
    }
    cat("\n",file=mortForcingFile,append=TRUE)
  }
  
}

##then turn it into a .nc file using command prompt and 
# ncgen -o mortForcingTEST.nc mortForcingTEST.txt

##now you need to edit the force.prm file
# see here for more details https://confluence.csiro.au/pages/viewpage.action?spaceKey=Atlantis&title=External+Mortality%2C+Growth+and+Recruitment+Scaling
# you need these lines

##############################
## mortality forcing file
# use_external_scaling 1
# scale_all_mortality 1
# externalBiologyForcingFile mortForcingTEST.nc
# externalBiologyForcingFile_rewind 1



