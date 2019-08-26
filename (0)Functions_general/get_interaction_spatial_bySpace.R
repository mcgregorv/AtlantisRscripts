get_interaction_spatial_bySpace<-function(predatorCode,preyCode,predatorAge,preyAge,groupsDF,ThisIC.nc,plotPath,biolLines){
  #gets spatial and temporal interaction level for given predator and prey (how often are they in the same place at the same time?)
  # predatorCode<-"HOK"; preyCode<-"ASQ"; predatorAge<-"ad"; preyAge<-"juv"
  
  thisVol<-ncvar_get(ThisIC.nc,"volume")
  thisDz<-ncvar_get(ThisIC.nc,"dz")
  if(length(dim(thisVol))==3){thisVol<-thisVol[,,1]}
  if(length(dim(thisDz))==3){thisDz<-thisDz[,,1]}
  nlayers<-dim(thisVol)[1]; nboxes<-dim(thisVol)[2]
  
  pred_numCohorts<-groupsDF$NumCohorts[groupsDF$Code==predatorCode]; prey_numCohorts<-groupsDF$NumCohorts[groupsDF$Code==preyCode]
  predName<-str_trim(groupsDF$Name[groupsDF$Code==predatorCode],side="both"); preyName<-str_trim(groupsDF$Name[groupsDF$Code==preyCode],side="both")
  
  if(predatorAge %in% c("juv","ad")){
    tempBiomassArray<-array(NA,dim=c(pred_numCohorts,dim(thisVol)))
    for(c in 1:pred_numCohorts){
      predSN<-ncvar_get(ThisIC.nc,paste(str_trim(predName,side="both"),c,"_StructN",sep="")); if(length(dim(predSN))==3){predSN<-predSN[,,1]}
      predRN<-ncvar_get(ThisIC.nc,paste(str_trim(predName,side="both"),c,"_ResN",sep="")); if(length(dim(predRN))==3){predRN<-predRN[,,1]}
      predNums<-ncvar_get(ThisIC.nc,paste(str_trim(predName,side="both"),c,"_Nums",sep="")); if(length(dim(predNums))==3){predNums<-predNums[,,1]}
      tempBiomassArray[c,,]<-predNums*(predSN+predRN)
    }
    predAgeMature<-get_first_number(biolLines[grep(paste(predatorCode,"_age_mat",sep=""),biolLines)]) 
    #only want the predatorAge and preyAge appropriate ones
    if(predatorAge=="juv"){
      temp<-tempBiomassArray[1:predAgeMature,,]
      if(length(dim(temp))==2){
        predBiomass<-temp
      } else{
        predBiomass<-apply(temp,c(2,3),sum,na.rm=TRUE)
      }
    } else{
      temp<-tempBiomassArray[(predAgeMature+1):pred_numCohorts,,]
      if(length(dim(temp))==2){
        predBiomass<-temp
      } else{
        predBiomass<-apply(temp,c(2,3),sum,na.rm=TRUE)
      }
    }
    data_pred<-predBiomass
  }else{
    name_pred<-paste(str_trim(groupsDF$Name[groupsDF$Code==predatorCode]),"_N",sep="")
    data_pred<-ncvar_get(ThisIC.nc,name_pred); 
    if(length(dim(data_pred))==3){data_pred<-data_pred[,,1]}
  }
  if(preyAge %in% c("juv","ad")){
    tempBiomassArray<-array(NA,dim=c(prey_numCohorts,dim(thisVol)))
    for(c in 1:prey_numCohorts){
      preySN<-ncvar_get(ThisIC.nc,paste(str_trim(preyName,side="both"),c,"_StructN",sep="")); if(length(dim(preySN))==3){preySN<-preySN[,,1]}
      preyRN<-ncvar_get(ThisIC.nc,paste(str_trim(preyName,side="both"),c,"_ResN",sep="")); if(length(dim(preyRN))==3){preyRN<-preyRN[,,1]}
      preyNums<-ncvar_get(ThisIC.nc,paste(str_trim(preyName,side="both"),c,"_Nums",sep="")); if(length(dim(preyNums))==3){preyNums<-preyNums[,,1]}
      tempBiomassArray[c,,]<-preyNums*(preySN+preyRN)
    }
    preyAgeMature<-get_first_number(biolLines[grep(paste(preyCode,"_age_mat",sep=""),biolLines)]) 
    #only want the predatorAge and preyAge appropriate ones
    if(preyAge=="juv"){
      temp<-tempBiomassArray[1:preyAgeMature,,]
      if(length(dim(temp))==2){
        preyBiomass<-temp
      }else{
        preyBiomass<-apply(temp,c(2,3),sum,na.rm=TRUE)
      }
    } else{
      temp<-tempBiomassArray[(preyAgeMature+1):prey_numCohorts,,]
      if(length(dim(temp))==2){
        preyBiomass<-temp
      }else{
        preyBiomass<-apply(temp,c(2,3),sum,na.rm=TRUE)
      }
    }
    data_prey<-preyBiomass
  }else{
    name_prey<-paste(str_trim(groupsDF$Name[groupsDF$Code==preyCode]),"_N",sep="")
    data_prey<-ncvar_get(ThisIC.nc,name_prey); 
    if(length(dim(data_prey))==3){data_prey<-data_prey[,,1]}
  }
  #make sure they are same dimensions
  ld_pred<-length(dim(data_pred)); ld_prey<-length(dim(data_prey))  
  if(ld_prey>ld_pred){
    data_prey<-data_prey[nlayers,,]
  } else if(ld_pred>ld_prey){
    data_pred<-data_pred[nlayers,,]
  }
  
  #if sediment of predator is zero, but it is in the bottom layer, it may be able to access sediment
  #should read this in from biol.prm file, but for now I'm just going to set sediment biomass to 0.1 times that in lowest wc layer
  if(length(dim(data_prey))==0 | sum(data_prey[nlayers,],na.rm=TRUE)>0){
    if(length(dim(data_pred))>1){
      if(sum(data_pred[nlayers,],na.rm=TRUE)==0 & sum(data_pred[1,],na.rm=TRUE)>0){
        data_pred[nlayers,]<-0.1*data_pred[1,]
      }
    }
  }
  #for each datapoint, if both >0 give score 1 and 0 otherwise
  overlap<-0*data_pred
  #als get the number of cells that have either the prey or the predator in them
  potentialOverlapSpace<-0*overlap
  for(b in 1:nboxes){
    if(length(dim(data_pred))==3){
      for(l in 1:nlayers){
        if(data_prey[l,b]>0 & data_pred[l,b]){overlap[l,b]<-1}
        if(data_prey[l,b]>0 | data_pred[l,b]){potentialOverlapSpace[l,b]<-1}
      }
    } else{
      if(data_prey[b]>0 & data_pred[b]){overlap[b]<-1}
      if(data_prey[b]>0 | data_pred[b]){potentialOverlapSpace[b]<-1}
    }
  }
  # overlap_sum<-sum(overlap,na.rm=TRUE)
  # potentialSum<-sum(potentialOverlapSpace,na.rm=TRUE)
  # #turn it into proportions of total posible
  # # thisMax<-max(overlap_sum)
  # # overlap_prop<-overlap_sum/thisMax
  # overlap_prop<-0
  # if(potentialSum>0){
  #   overlap_prop<-overlap_sum/potentialSum
  # }
  # return(overlap_prop)
  return(potentialOverlapSpace)
}