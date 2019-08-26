#first take an example species with 2 cohorts. Suppose they need to grow enough in cohort 1
predWeights<-c(1603,	4342,	7534, 10592,	13233,	15381,	17061,	18343,	19304,	20017) + 
  c(4246,	11504,	19963, 28064,	35063,	40753,	45205,	48601,	51149,	53037)
predWeights<-signif(predWeights,2)

#set proportions of body weight consumed per day for cohort 1 and cohort n
lambda0<-0.3; lambda1<-0.1; thisE<-0.5

lambda1<-0.16 #this gives a1~0.75 for hoki

maxWeight<-max(predWeights); minWeight<-min(predWeights)

calcGrowth<-function(x,out_a1=FALSE){
  a1<-1-((log(lambda0,base=exp(1))-log(lambda1,base=exp(1)))/(log(maxWeight,base=exp(1))-log(minWeight,base=exp(1))))
  a0<-thisE*lambda0*x^(1-a1)
  mum<-a0*x^a1
  if(out_a1==TRUE){
    return(a1)
  }else{
    return(mum)
  }
}

##this is mum
growthByCohort<-unlist(lapply(predWeights,calcGrowth))
a1_byCohort<-unlist(lapply(predWeights,calcGrowth,out_a1=TRUE))


plot(x=seq(1,length(predWeights)),y=predWeights,xlab="Cohort",ylab="mg N",type="l",lwd=2,col=myPurple)
points(x=seq(1,length(predWeights)),y=predWeights,pch=20,col=myPurple)

plot(x=seq(1,length(predWeights)),y=growthByCohort,xlab="Cohort",ylab="Maximum growth (mg N assimilated per day)",type="l",lwd=2,col=myPurple)
points(x=seq(1,length(predWeights)),y=growthByCohort,pch=20,col=myPurple)

#if each day they grow mum per day, how big do they get??
#do cohort 1 first
nt<-365
testGrowth<-rep(0,nt)
testGrowth[1]<-predWeights[1]
for(t in 2:nt){
  testGrowth[t]<-testGrowth[(t-1)]+growthByCohort[1]
}
plot(testGrowth,type="l")
abline(h=predWeights,col=myBlue_trans,lwd=2)

##################### another way..?
calcGrowth<-function(x){
  nx<-length(x)
  y<-rep(0,nx)
  for(i in 1:(nx-1)){
    y[i]<-(mean(x[i:(i+1)])-x[i])/365
  }
  y[nx]<-2*(x[nx]-mean(x[(nx-1):nx]))/365
  return(y)
}

calcGrowthMax<-function(x){
  nx<-length(x)
  y<-rep(0,nx)
  for(i in 1:(nx-1)){
    y[i]<-(x[(i+1)]-x[i])/365
  }
  y[nx]<-y[nx-1]/2
  return(y)
}

testMum<-calcGrowthMax(predWeights)
testGrowth<-rep(0,nt)
testGrowth[1]<-predWeights[1]
for(t in 2:nt){
  testGrowth[t]<-testGrowth[(t-1)]+testMum[1]
}
plot(testGrowth,type="l")
abline(h=predWeights,col=myBlue_trans,lwd=2)

