get_interaction_trophicLevel<-function(predatorCode,preyCode,predatorAge,preyAge,groupsTL){
  #gets trophic level interaction level for given predator and prey 
  # predatorCode<-"HOK"; preyCode<-"ASQ"
  
  limitProp<-0
  
  #this has predators eating up to round(TL,0.5,"up) and down to round(TL-1,0.5,"down)
  if(predatorAge=="juv"){
    predTL<-groupsTL$TrophicLevel1[groupsTL$Code==predatorCode]
  } else{
    predTL<-groupsTL$Isotope[groupsTL$Code==predatorCode]
    if(is.na(predTL)){predTL<-groupsTL$TrophicLevel2[groupsTL$Code==predatorCode]}
  }
  minTL<-myRounding(predTL-1,0.5,"down")
  maxTL<-myRounding(predTL,0.5,"up")
  if(preyAge=="juv"){
    preyTL<-groupsTL$TrophicLevel1[groupsTL$Code==preyCode]
  } else {
    preyTL<-groupsTL$TrophicLevel2[groupsTL$Code==preyCode]
  }
  
  if(preyTL>=minTL & preyTL<=maxTL){limitProp<-1}
  
  return(limitProp)
}