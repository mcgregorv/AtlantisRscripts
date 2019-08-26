getAssumedNaturalMortality<-function(groupCode,ThisIC.nc,groupsDF,biolLines){
  # groupCode<-"HOK"
  thisMday<-NA
  #we set up the initial conditions using an assumed natural mortality, such that numbers in cohort C = numbers in cohort c-1 *(exp(-M))
  calcMday<-function(N1,N2,ageClassSize){
    Mday<-(-1)*(1/(ageClassSize*365))*log((N2/N1),base=exp(1))
  }
  thisName<-str_trim(groupsDF$Name[groupsDF$Code==groupCode],side='both')
  thisNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==groupCode]
  if(thisNumCohorts>1){
    #get total numbers by cohort
    numsByCohort<-rep(0,thisNumCohorts)
    for(c in 1:thisNumCohorts){
      numsByCohort[c]<-getICNumbers(name=thisName,cohort=c,ThisIC.nc)
    }
    #get age class size
    thisVar<-paste(groupCode,"_AgeClassSize",sep="")
    thisAgeSize<-get_first_number(biolLines[grep(thisVar,biolLines)])
    thisMdays<-mapply(calcMday,numsByCohort[1:(thisNumCohorts-1)],numsByCohort[2:thisNumCohorts],ageClassSize=thisAgeSize)
    thisMday<-mean(thisMdays,na.rm=TRUE)
  }
  return(thisMday)
}
  