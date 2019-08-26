getMortalityRequiredPerDay<-function(groupCode,groupsDF,biolLines,ThisIC.nc){
  # groupCode<-"HOK"
  #runs getAssumedNaturalMortality() to get M by day,
  #then multiplies by numbers spatially to get desired deaths per day
  thisDeathRequired<-NA
  # groupCode<-"HOK"
  #this uses mum in biol.prm, so assumes that has been calculated appropriately 
  #in my case, see Estimate_mum_byCOhort.R and Overwrite_mum.R
  thisNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==groupCode]
  thisName<-str_trim(groupsDF$Name[groupsDF$Code==groupCode],side="both")
  if(thisNumCohorts>1){
    #get mortality by day
    thisMday<-getAssumedNaturalMortality(groupCode,ThisIC.nc=ThisIC.nc,groupsDF=groupsDF,biolLines=biolLines)
    #this is instantaneous. to get proprtion
    thisMpropday<-1-exp(-thisMday)
    #get numbers by cell to covert to absolute food required
    #first get volume to get dimensions right
    thisVol<-ncvar_get(ThisIC.nc,"volume"); if(length(dim(thisVol))==3){thisVol<-thisVol[,,1]}
    deathRequired<-array(NA,dim=c(thisNumCohorts,dim(thisVol)))
    for(c in 1:thisNumCohorts){
      thisVar<-paste(thisName,c,"_Nums",sep="")
      xx<-ncvar_get(ThisIC.nc,thisVar); if(length(dim(xx))==3){xx<-xx[,,1]}
      deathRequired[c,,]<-xx
    }
  }
  return(thisFoodRequired)
}