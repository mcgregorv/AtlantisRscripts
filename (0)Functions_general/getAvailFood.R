getAvailFood<-function(groupCode,predAge,groupsDF,sampledAvails,biolLines,volume,thisLayer,thisBox){
  # groupCode<-"HOK"; predAge<-"ad"; thisLayer<-1; thisBox<-2; 
  # groupCode<-"ASQ"; predAge<-"juv"
  # groupCode<-"BC"; predAge<-""; thisLayer=6; thisBox=2;
  #for each predator cohort, loops through prey
  #for each prey/cohort, gets avail food first from matrix, then checks gape size, max deaths per day, and stores
  #ACTUALLY, SHOULD ALLOW TO EAT THEM ALL, AS THIS IS WHAT HAPPENS IN ATLANTIS. Eat down to zero, but penalise if eat more than M
  # cat(dim(volume))
  thisCellVol<-volume[thisLayer,thisBox]
  nlayers<-dim(volume)[1]
  thisRowName<-paste(groupCode,predAge,sep="")
  thisPreyAvail<-sampledAvails[grep(thisRowName,rownames(sampledAvails)),]
  thisNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==groupCode]
  thisName<-str_trim(groupsDF$Name[groupsDF$Code==groupCode])
  nonZeroPreys<-colnames(thisPreyAvail)[as.double(thisPreyAvail)>0]; npreys<-length(nonZeroPreys)
  if(thisNumCohorts>1){
    #get KLP and KUP for predator
    thisVar<-paste("KLP_",groupCode,sep=""); thisLine<-biolLines[grep(thisVar,biolLines)]
    thisKLP<-get_first_number(thisLine)
    thisVar<-paste("KUP_",groupCode,sep=""); thisLine<-biolLines[grep(thisVar,biolLines)]
    thisKUP<-get_first_number(thisLine)
    predSN<-unlist(lapply(seq(1,thisNumCohorts),FUN=getSNRN,name=thisName,whichN="SN",ThisIC.nc=ThisIC.nc))
    lowerSNlimits<-thisKLP*predSN; upperSNlimits<-thisKUP*predSN
    
    
    #get age mature for predator so can loop through included cohorts
    predAgeMature<-get_first_number(biolLines[grep(paste(groupCode,"_age_mat",sep=""),biolLines)]) 
    
    if(predAge=="juv"){
      firstCohort<-1; lastCohort<-predAgeMature
    }else{
      firstCohort<-predAgeMature+1; lastCohort<-thisNumCohorts
    }
    sumPreyBiomass<-rep(0,(lastCohort-firstCohort+1))
    for(c in firstCohort:lastCohort){
      cINdex<-c-firstCohort+1
      thisGapeL<-lowerSNlimits[c]; thisGapeU<-upperSNlimits[c]
      for(prey in 1:npreys){
        #loop through preys, for each calculate how much (mg N) within gape size in this cell, multiply by sampledAvail
        thisPreyVar<-nonZeroPreys[prey]
        # cat("prey ",thisPreyVar,",,")
        thisPreyCode<-gsub('juv|ad',"",thisPreyVar)
        thisPreyAge<-gsub(thisPreyCode,"",thisPreyVar)
        thisPreyNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisPreyCode]
        thisPreyName<-str_trim(groupsDF$Name[groupsDF$Code==thisPreyCode],side="both")
        ##can we pre-populate this? it's really slow in here
        thisPreyBiomass<-getPreyBiomassAvail(thisPreyCode,thisPreyName,thisPreyNumCohorts,thisPreyAge,ThisIC.nc,volume,thisLayer,thisBox,biolLines,thisGapeL,thisGapeU)
        if(!is.na(thisPreyBiomass)){   
          sumPreyBiomass[cINdex]<-sumPreyBiomass[cINdex]+thisPreyBiomass
        }
      }
    }
  } else{
    sumPreyBiomass<-0
    for(p in 1:npreys){
      thisPreyVar<-nonZeroPreys[p]
      thisPreyCode<-gsub('juv|ad',"",thisPreyVar)
      thisPreyAge<-gsub(thisPreyCode,"",thisPreyVar)
      thisPreyNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisPreyCode]
      thisPreyName<-str_trim(groupsDF$Name[groupsDF$Code==thisPreyCode],side="both")
      thisPreyBiomass<-getPreyBiomassAvail(thisPreyCode,thisPreyName,thisPreyNumCohorts,thisPreyAge,ThisIC.nc,volume,thisLayer,thisBox,biolLines,thisGapeL = NA, thisGapeU = NA)
      
      sumPreyBiomass<-sumPreyBiomass+thisPreyBiomass
    }
  }
  #return avail biomass
  return(sumPreyBiomass)
}