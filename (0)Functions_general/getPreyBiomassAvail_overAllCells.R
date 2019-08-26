getPreyBiomassAvail<-function(thisPreyCode,thisPreyName,thisPreyNumCohorts,thisPreyAge,ThisIC.nc,volume,biolLines,thisGapeL=NA,thisGapeU=NA){
  #there is a version of this that returns one value for a given cell. THis returns an array of dimensions matching 'volume'
  thisPreyBiomass<-0
  ##called in getAvailFood
  if(thisPreyNumCohorts>1){
    tempBiomassArray<-array(0,dim=c(thisPreyNumCohorts,dim(volume)))
    for(c in 1:thisPreyNumCohorts){
      preySN<-getSNRN(name=thisPreyName,cohort=c,whichN="SN",ThisIC.nc=ThisIC.nc)
      preyRN<-getSNRN(name=thisPreyName,cohort=c,whichN="RN",ThisIC.nc=ThisIC.nc)
      preyNums<-ncvar_get(ThisIC.nc,paste(str_trim(thisPreyName,side="both"),c,"_Nums",sep="")); if(length(dim(preyNums))==3){preyNums<-preyNums[,,1]}
      if(!is.na(thisGapeL) & !is.na(thisGapeU)){
        if(preySN>=thisGapeL & preySN<=thisGapeU){
          tempBiomassArray[c,,]<-preyNums*(preySN+preyRN)
        }
      }else{
        tempBiomassArray[c,,]<-preyNums*(preySN+preyRN)
      }
    }
    #get age mature, so know which cohorts to include
    preyAgeMature<-get_first_number(biolLines[grep(paste(thisPreyCode,"_age_mat",sep=""),biolLines)]) 
    if(thisPreyAge=="juv"){
      preyBiomass<-tempBiomassArray[1:preyAgeMature]
    } else{
      preyBiomass<-tempBiomassArray[(preyAgeMature+1):thisPreyNumCohorts]
    }
    thisPreyBiomass<-sum(preyBiomass,na.rm=TRUE)
  }else{
    #no cohorts, and at this stage no gape size limitations, so just get mg N then multiply by prop avail
    thisPreyICvar<-paste(thisPreyName,"_N",sep="")
    temp<-ncvar_get(ThisIC.nc,thisPreyICvar)
    if(length(dim(temp))==3){temp<-temp[,,1]}
    if(length(dim(temp))==2){
      if(temp[1,2]>9e+36){
        temp<-temp[,1]
      }
    }
    if(length(dim(temp))==0){
      #if prey is in the sed layer, it's only accessed if pred is there or if predator is in 
      # bottom water column layer. In the latter case, just make 0.01 * prey avail
      if(thisLayer==nlayers){
        thisPreyBiomass<-temp[thisBox]*thisCellVol
      }else if(thisLayer==1){
        thisPreyBiomass<-0.01*temp[thisBox]*thisCellVol
      }
    } else{
      thisPreyBiomass<-temp[thisLayer,thisBox]*thisCellVol
    }
  }
  return(thisPreyBiomass)
}
