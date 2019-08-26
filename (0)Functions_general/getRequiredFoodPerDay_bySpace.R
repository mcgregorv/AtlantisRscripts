getRequiredFoodPerDay_bySpace<-function(groupCode,groupsDF,biolLines,ThisIC.nc){
  # groupCode<-"HOK"
  #this uses mum in biol.prm, so assumes that has been calculated appropriately 
  #in my case, see Estimate_mum_byCOhort.R and Overwrite_mum.R
  thisVol<-ncvar_get(ThisIC.nc,"volume"); if(length(dim(thisVol))==3){thisVol<-thisVol[,,1]}
  thisFoodRequired<-array(NA,dim=c(10,dim(thisVol)))
  thisNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==groupCode]
  thisName<-str_trim(groupsDF$Name[groupsDF$Code==groupCode],side="both")
  if(thisNumCohorts>1){
    #get mum from biolLines
    thisVar<-paste("mum_",groupCode,sep="")
    thisMums<-get_first_number(biolLines[grep(thisVar,biolLines)+1],n="all")
    #get efficiency
    thisVar<-paste("E_",groupCode,sep="")
    thisE<-get_first_number(biolLines[grep(thisVar,biolLines)])
    #this is per individual
    indFoodRequired<-thisMums/thisE
    #get numbers by cell to covert to absolute food required
    #first get volume to get dimensions right
    
    foodRequired<-array(NA,dim=c(thisNumCohorts,dim(thisVol)))
    for(c in 1:thisNumCohorts){
      thisVar<-paste(thisName,c,"_Nums",sep="")
      xx<-ncvar_get(ThisIC.nc,thisVar); if(length(dim(xx))==3){xx<-xx[,,1]}
      foodRequired[c,,]<-xx*indFoodRequired[c]
    }
    thisFoodRequired[1:thisNumCohorts,,]<-foodRequired
  } else{
    #not sure what we want here - need to check Atlantis code for invert growth wrt food consumption
  }
  return(thisFoodRequired)
}