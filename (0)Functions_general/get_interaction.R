get_interaction<-function(ThisNC.nc,predatorCode,preyCode,groupsDF){
  #gets spatial and temporal interaction level for given predator and prey (how often are they in the same place at the same time?)
  # predatorCode<-"IVS"; preyCode<-"MB"
  
  #first plot biomass by trophic level
  this_run<-"base"
  this_out<-paste("test",sep="")

  this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
  outPath<-paste(this_path,"output",this_out,"\\",sep="")
  
  plotPath<-paste(this_path,"..\\Figures\\",this_run,"\\",this_out,"",sep="")
  
  #read in CRAM groups file
  groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep=""))
  
  #read in tracers from model run
  ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
  thisVol<-ncvar_get(ThisNC.nc,"volume")
  thisDz<-ncvar_get(ThisNC.nc,"dz")
  nts<-dim(thisVol)[3]; nlayers<-dim(thisVol)[1]; nboxes<-dim(thisVol)[2]
  
  name_pred<-paste(str_trim(groupsDF$Name[groupsDF$Code==predatorCode]),"_N",sep="")
  name_prey<-paste(str_trim(groupsDF$Name[groupsDF$Code==preyCode]),"_N",sep="")
  data_pred<-ncvar_get(ThisNC.nc,name_pred); ld_pred<-length(dim(data_pred))
  data_prey<-ncvar_get(ThisNC.nc,name_prey); ld_prey<-length(dim(data_prey))  
  if(ld_prey>ld_pred){
    data_prey<-data_prey[nlayers,,]
  } else if(ld_pred>ld_prey){
    data_pred<-data_pred[nlayers,,]
  }
  #for each datapoint, if both >0 give score 1 and 0 otherwise
  overlap<-0*data_pred
  for(t in 1:nts){
    for(b in 1:nboxes){
      if(length(dim(data_pred))==3){
        for(l in 1:nlayers){
          if(data_prey[l,b,t]>0 & data_pred[l,b,t]){overlap[l,b,t]<-1}
        }
      } else{
        if(data_prey[b,t]>0 & data_pred[b,t]){overlap[b,t]<-1}
      }
    }
  }
  overlap_sum<-apply(overlap,c(1,2),sum)
  #turn it into proportions of total posible
  # thisMax<-max(overlap_sum)
  # overlap_prop<-overlap_sum/thisMax
  overlap_prop<-overlap_sum/nts
  return(overlap_prop)
}