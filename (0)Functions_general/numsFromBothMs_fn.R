numsFromBothMs_fn<-function(Mj, Ma, nyears=250, BHpars, ageClassSize){
  #calc BH alpha for given R0
  simNums_array<-array(NA,dim=c(thisNumCohorts, nyears)); simNums_array[,1]<-numbers_array[,1]
  for(y in 2:nyears) {
    for(c in 1:thisNumCohorts){
      thisM<-Ma; if(c<3){thisM<-Mj}
      if(c==1){
        #calc recruits
        thisSSB<-sum(simNums_array[, y-1] * weight_array[,1]) / (mg_2_tonne * X_CN)
        thisRecruits<-calcR(S=thisSSB, a=BHpars[1], b=BHpars[2])
      } else{
        thisRecruits<-simNums_array[c-1, y-1]*((ageClassSize-1)/ageClassSize) #these are those that aged up from previous cohort
      }
      # add the recruits, take out the aged, and take out the mortality
      simNums_array[c,y] <- (simNums_array[c,y-1]*((ageClassSize-1)/ageClassSize) + thisRecruits)*exp((-1)*(thisM*ageClassSize))
      #if it's the last age class, kill any extra
      if(c==thisNumCohorts){
        simNums_array[c,y]<-simNums_array[c,1]
      }
    }
  }
  ##output the last timestep of numbers to compare with the initial conditions
  thisNumFromM<-simNums_array[,nyears]
  return(thisNumFromM)
}