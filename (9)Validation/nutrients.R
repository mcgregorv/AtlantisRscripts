#################################### TRACERS
#plot all tracers for a given box and layer
this_run<-"base"
this_out<-"TEST150yrfish"

burnin<-35 #number of years to skip in plot

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")
baseOutPath<-paste(this_path,"output",base_out,"\\", sep="")

#read in oxygen data
O2data<-read.csv(paste(DIR$'Data', "nutrients\\WOA_oxygenReady2Compare.csv", sep=""))
  
this_dz<-ncvar_get(ThisNC.nc, "dz")               
this_vol<-ncvar_get(ThisNC.nc,"volume")

