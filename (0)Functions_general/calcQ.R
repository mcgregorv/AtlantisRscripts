calcQ<-function(biomassByGroup){
  ## calculate kempton's Q diversity index
  ## ainsworth suggests dropping species groups that drop below ** of B0 before calculating Q
  #biomassByGroup is a vector of biomass' for each species group
  getMyBin<-function(x, bins=biomassBins){
    index<-bins<=x
    thisBiomassBin<-max(bins[index], na.rm=TRUE)
    bbIndex<-seq(1,length(bins))[bins==thisBiomassBin][1]
    if(length(bbIndex)==0){bbIndex<-length(bins)+1}
    return(bbIndex)
  }
  thisS<-length(biomassByGroup)
  ## define biomass bins
  xx<-sort(biomassByGroup); lx<-max(1e-6, xx[round(0.05*length(xx))]); ux<-xx[round(0.95*length(xx))]
  biomassBins<-exp(seq(log(lx), log(ux), length.out=20))
  ## asign biomass bin to each species group
  groupBiomassBins<-unlist(lapply(biomassByGroup, getMyBin))
  ## get cummulative number of species groups by biomass bin
  numGroupsByBin<-table(groupBiomassBins)
  cumNumbersByBin<-rep(0,length(biomassBins))
  for(b in 1:length(biomassBins)){
    x<-numGroupsByBin[as.double(names(numGroupsByBin))<=b]
    cumNumbersByBin[b]<-sum(x)
  }
  ## get biomassBin for 10th and 90th species percentiles
  bin10perc<-biomassBins[getMyBin(x=0.1*thisS, bins=cumNumbersByBin)]
  bin90perc<-biomassBins[getMyBin(x=0.9*thisS, bins=cumNumbersByBin)]
  
  ## calculate Q
  Q<-(0.8*thisS)/(log(bin90perc/bin10perc))
  return(Q)
}
