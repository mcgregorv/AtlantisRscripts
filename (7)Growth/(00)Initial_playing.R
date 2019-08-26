simpleHII<-function(C,Gmax,E,N){
  # C is clearance rate (how many m^3 water searched per day)
  # Gmax is maximum growth rate - if one dominant prey, handling time can be swapped out by efficiency/max_growth
  # E is efficiency - the fraction of nitrogen in prey assimilated
  # N is density of prey (numbers per m^3) - although sometimes in Atlantis it is biomass, so we should test examples of this too
  A<-(C*N)/(1+(C/Gmax)*E*N) #this is the attack rate - multiply it by the density of predator to get total prey attacked per m^3 per day
  return(A)
}

dynamicHIIsinglePrey<-function(C,ht,N){
  ##ACTUALLY SAME AS SIMPLE IF ONE PREY TYPE
  # C is clearance rate (how many m^3 water searched per day)
  # ht is handling time, in days, how much time one predator takes to 'handle' one prey 
  #   (covers all time spent once they attack the prey through to starting search again)
  # N is density of prey (numbers per m^3) - although sometimes in Atlantis it is biomass, so we should test examples of this too
  A<-(C*N)/(1+(C*ht*N)) #this is the attack rate - multiply it by the density of predator to get total prey attacked per m^3 per day
  return(A)
}

dynamicHII<-function(C,ht,preys,Es){
  # C is clearance rate (how many m^3 water searched per day)
  # ht is handling time, in days, how much time one predator takes to 'handle' one prey 
  #   (covers all time spent once they attack the prey through to starting search again)
  # N is density of prey (numbers per m^3) - although sometimes in Atlantis it is biomass, so we should test examples of this too
  npreys<-length(preys)
  topSum<-0; bottomSum<-0
  for(n in 1:npreys){
    topSum<-topSum+Es[n]*preys[n]
    bottomSum<-bottomSum+ht*preys[n]
  }
  A<-(C*topSum)/(1+(C*bottomSum)) #this is the attack rate 
  return(A)
}


testC<-10; testGmax<-5; testE<-0.8; testN<-100; testHt<-0.1
testGmax<-testE/testHt

testA<-simpleHII(C=testC,Gmax = testGmax, E= testE, N=testN)
testDynA<-dynamicHII(C=testC,ht=testHt, N=testN)


vectorN<-seq(0,1000); vectorAttack<-unlist(lapply(vectorN,simpleHII,C=testC,Gmax = testGmax, E= testE))
plot(x=vectorN,y=vectorAttack,type="l",col=myBlue)

vectorDynAttack<-unlist(lapply(vectorN,dynamicHII,C=testC,ht=testHt))
points(x=vectorN,y=vectorDynAttack,type="l",col=myOrange,lty=2)

##what if there are 2 preys equal amounts, and the predator has different efficencies for both?
testPrey1<-seq(0,500); testPrey2<-testPrey1; testE1<-0.5; testE2<-0.2
testPreyDF<-data.frame(rbind(testPrey1,testPrey2))
testPredationDyn<-apply(testPreyDF,2,dynamicHII,Es=c(testE1,testE2),C=testC,ht=testHt)

testTotalPrey<-testPrey1+testPrey2

plot(testSimplePred/testPredationDyn,type="l",col=myRed,lwd=2)


plot(x=testTotalPrey,y=testPredationDyn,type="l",col=myOrange,lwd=2)

testSimplePred<-unlist(lapply(testTotalPrey,simpleHII,C=testC,Gmax = testGmax, E= testE1))
plot(x=testTotalPrey,y=testSimplePred,type="l",col=myBlue,lwd=2,lty=2)

points(x=testTotalPrey,y=testPredationDyn,type="l",col=myOrange,lwd=2)


#and if the preys have different proporitons - first prey 1 with E =0.5 is bigger
testPrey1<-seq(0,500); testPrey2<-0.5*testPrey1; testE1<-0.5; testE2<-0.2
testPreyDF<-data.frame(rbind(testPrey1,testPrey2))
testPredationDyn<-apply(testPreyDF,2,dynamicHII,Es=c(testE1,testE2),C=testC,ht=testHt)

testTotalPrey<-testPrey1+testPrey2

plot(x=testTotalPrey,y=testPredationDyn,type="l",col=myOrange,lwd=2)

testSimplePred<-unlist(lapply(testTotalPrey,simpleHII,C=testC,Gmax = testGmax, E= testE1))
plot(x=testTotalPrey,y=testSimplePred,type="l",col=myBlue,lwd=2,lty=2)

points(x=testTotalPrey,y=testPredationDyn,type="l",col=myOrange,lwd=2)

#difference is not quite as big - test simple is closer to test dynamic
plot(testSimplePred/testPredationDyn,type="l",col=myRed,lwd=2)



#and if the preys have different proporitons - first prey 1 with E =0.5 is bigger
testPrey1<-seq(0,500); testPrey2<-2*testPrey1; testE1<-0.5; testE2<-0.2
testPreyDF<-data.frame(rbind(testPrey1,testPrey2))
testPredationDyn<-apply(testPreyDF,2,dynamicHII,Es=c(testE1,testE2),C=testC,ht=testHt)

testTotalPrey<-testPrey1+testPrey2

plot(x=testTotalPrey,y=testPredationDyn,type="l",col=myOrange,lwd=2)

testSimplePred<-unlist(lapply(testTotalPrey,simpleHII,C=testC,Gmax = testGmax, E= testE1))
plot(x=testTotalPrey,y=testSimplePred,type="l",col=myBlue,lwd=2,lty=2)

points(x=testTotalPrey,y=testPredationDyn,type="l",col=myOrange,lwd=2)

#difference is bigger - test simple is further from test dynamic
plot(testSimplePred/testPredationDyn,type="l",col=myRed,lwd=2)


#if you don't have prey on the numerator, what does this do??
parslowHII<-function(C,Gmax,E,N){
  A<-(C)/(1+(C/Gmax)*E*N) #this is the attack rate - multiply it by the density of predator to get total prey attacked per m^3 per day
  return(A)
}
testParslow<-unlist(lapply(testTotalPrey,parslowHII,C=testC,Gmax = testGmax, E= testE1))
plot(x=testTotalPrey,y=testParslow,type="l",col=myBlue,lwd=2,lty=2,xlab="Prey",ylab="CLEAR")

#but if then later you multiply this rate by the prey..? 
testParslowPrey<-testParslow*testTotalPrey
plot(x=testTotalPrey,y=testParslowPrey,type="l",col=myBlue,lwd=2,lty=2,xlab="Prey",ylab="CLEAR*Prey")

points(x=testTotalPrey,y=testSimplePred,type="l",col=myRed,lwd=2,lty=3)
##this is the same

##but it affects when the mum limit kicks in 
