get_interaction_gapeSize<-function(predatorCode,preyCode,predatorAge,preyAge,groupsDF,ThisIC.nc,plotPath,biolLines){
  #gets spatial and temporal interaction level for given predator and prey (how often are they in the same place at the same time?)
  # predatorCode<-"HOK"; preyCode<-"ASQ"
  ICnames<-sort(names(ThisIC.nc$var))
  thisVol<-ncvar_get(ThisIC.nc,"volume")
  thisDz<-ncvar_get(ThisIC.nc,"dz")
  if(length(dim(thisVol))==3){thisVol<-thisVol[,,1]}
  if(length(dim(thisDz))==3){thisDz<-thisDz[,,1]}
  nlayers<-dim(thisVol)[1]; nboxes<-dim(thisVol)[2]
  
  pred_numCohorts<-groupsDF$NumCohorts[groupsDF$Code==predatorCode]; prey_numCohorts<-groupsDF$NumCohorts[groupsDF$Code==preyCode]
  predName<-str_trim(groupsDF$Name[groupsDF$Code==predatorCode],side="both"); preyName<-str_trim(groupsDF$Name[groupsDF$Code==preyCode],side="both")
    
  #default to probablity = 1, which is what it will stay at if either groups are not age-structured. 
  #need to check out how Atlantis does this, but this will do for now
  limitProp<-1
  
  ##check if both predator and prey are age-structured as that is the easiest and we'll do that first
  if(pred_numCohorts>1 & prey_numCohorts>1){
    #get the SN of both, as this is what gape size is used on
    preySN<-unlist(lapply(seq(1,prey_numCohorts),FUN=getSNRN,name=preyName,whichN="SN",ThisIC.nc=ThisIC.nc))
    predSN<-unlist(lapply(seq(1,pred_numCohorts),FUN=getSNRN,name=predName,whichN="SN",ThisIC.nc=ThisIC.nc))
    #get predator gapesize from biol.prm file
    thisVar<-paste("KLP_",predatorCode,sep=""); thisLine<-biolLines[grep(thisVar,biolLines)]
    thisKLP<-get_first_number(thisLine)
    thisVar<-paste("KUP_",predatorCode,sep=""); thisLine<-biolLines[grep(thisVar,biolLines)]
    thisKUP<-get_first_number(thisLine)
    lowerSNlimits<-thisKLP*predSN; upperSNlimits<-thisKUP*predSN
    #need to split predator and prey into adult and juv
    #this is actually the last cohort that is immature
    predAgeMature<-get_first_number(biolLines[grep(paste(predatorCode,"_age_mat",sep=""),biolLines)]) 
    preyAgeMature<-get_first_number(biolLines[grep(paste(preyCode,"_age_mat",sep=""),biolLines)])
    #only want the predatorAge and preyAge appropriate ones
    if(predatorAge=="juv"){
      lowerSNlimits<-lowerSNlimits[1:predAgeMature]; upperSNlimits<-upperSNlimits[1:predAgeMature];
    } else{
      lowerSNlimits<-lowerSNlimits[(predAgeMature+1):pred_numCohorts]; upperSNlimits<-upperSNlimits[(predAgeMature+1):pred_numCohorts];
    }
    if(preyAge=="juv"){
      preySN<-preySN[1:preyAgeMature]; 
    } else{
      preySN<-preySN[(preyAgeMature+1):prey_numCohorts];
    }
    #probability of prey being within predators gape size
    numCombinations<-length(preySN)*length(lowerSNlimits); limitCount<-0
    for(i in 1:length(preySN)){
      for(j in 1:length(lowerSNlimits)){
        if(preySN[i]>=lowerSNlimits[j] & preySN[i]<=upperSNlimits[j]){ limitCount<-limitCount+1}
      }
    }
    limitProp<-limitCount/numCombinations
  
  }
  return(limitProp)
}