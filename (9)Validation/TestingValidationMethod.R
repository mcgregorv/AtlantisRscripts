#the two series need at least 3 matching data points to proceed with fit and plot
source(paste(DIR$'General functions',"makeBlankPlot.R",sep=""))
##test 1
testEst<-runif(11)*15

testObs<-seq(20,10)

#raw plot
plot(testEst,type="l",xlab="Timestep",ylab="Abundance",ylim=c(0,20))
points(testObs,pch=20)

if (sum(!is.na(testEst))>3){
  lmfit=nls(testEst~b+a*testObs,algorithm="port",start=c(a=1,b=0),lower=c(a=0,b=-Inf),upper=c(a=Inf,b=Inf))
  SS = sum(residuals(lmfit)^2)
  ScaledReal = abs(coef(lmfit)[1]) * testObs[!is.na(testObs)] + coef(lmfit)[2]
  resids<-testObs-testEst; meanResids<-mean(resids)
  newScaledO<-testObs-meanResids
  
  points(newScaledO,col=myBlue,pch=3)
  
  points(ScaledReal,pch=8,col=myOrange)
}
legend(legend=c("Observed","Scaled observed","Estimated"),lwd=c(NA,NA,1),col=c("black",myOrange,"black"),pch=c(20,8,NA),x="topright",bty="n")

resids<-testObs-testEst
residsData<-data.frame(cbind(seq(1,length(resids)),resids)); colnames(residsData)<-c("x","resids")
residslmfit=lm(x~resids,data=residsData)
residslm_a<-summary(residslmfit)$'coefficients'[[1]]; residslm_b<-summary(residslmfit)$'coefficients'[[2]]; 

##test 2
testObs<-seq(10,20)
#raw plot
plot(testEst,type="l",xlab="Timestep",ylab="Abundance",ylim=c(0,20))
points(testObs,pch=20)

if (sum(!is.na(testEst))>3){
  lmfit=nls(testEst~b+a*testObs,algorithm="port",start=c(a=1,b=0),lower=c(a=0,b=-Inf),upper=c(a=Inf,b=Inf))
  SS = sum(residuals(lmfit)^2)
  ScaledReal = abs(coef(lmfit)[1]) * testObs[!is.na(testObs)] + coef(lmfit)[2]
  
  
  thisDF<-data.frame
  lmfit<-lm(testEst~testObs,data=thisDF)
  lmfit=nls(testEst~b+testObs,algorithm="port",start=c(a=1,b=0),lower=c(a=0,b=-Inf),upper=c(a=Inf,b=Inf))
  ScaledReal = testObs- coef(lmfit)[2]
  

  points(ScaledReal,pch=8,col=myOrange)
}
legend(legend=c("Observed","Scaled observed","Estimated"),lwd=c(NA,NA,1),col=c("black",myOrange,"black"),pch=c(20,8,NA),x="topright",bty="n")



pdf(paste(DIR$'Figures',"testingScalingObservationsMethod.pdf",sep=""),height=7)
par(mar=c(3,3.5,1,1),oma=c(1,1,1,1),mfrow=c(5,2))
for(i in 1:9){
  testEst<-runif(11)*15
  
  testObs<-runif(11)*25
  
  #raw plot
  plot(testEst,type="l",xlab="Timestep",ylab="Abundance",ylim=c(0,30),lwd=2.5)
  points(testObs,pch=16,col=myGrey,cex=1)
  points(testObs,type="l",lwd=2,col=myGrey_trans)
  
  
  if (sum(!is.na(testEst))>3){
    lmfit=nls(testEst~b+a*testObs,algorithm="port",start=c(a=1,b=0),lower=c(a=0,b=-Inf),upper=c(a=Inf,b=Inf))
    SS = sum(residuals(lmfit)^2)
    ScaledReal = abs(coef(lmfit)[1]) * testObs[!is.na(testObs)] + coef(lmfit)[2]
    # compare magnitude of change in obs and est
    maxObsDiff<-max(testObs)-min(testObs); maxEstDiff<-max(testEst)-min(testEst)
    thisScalar<-maxObsDiff/maxEstDiff
    resids<-testObs/thisScalar-testEst; meanResids<-mean(resids)
    newScaledO<-(testObs/thisScalar-meanResids)
    
    points(newScaledO,col=myBlue_trans,pch=16,cex=1)
    points(newScaledO,col=myBlue_trans,type="l",lwd=2,cex=1.5)
    
    points(ScaledReal,pch=8,col=myOrange,cex=1.5)
  }
}
makeBlankPlot()
legend(legend=c("Estimated","Observed","Observed (scaling 1)","Observed (scaling alt.)"),col=c("black",myGrey_trans,myOrange,myBlue_trans),
       lwd=c(2,2,NA,3),pch=c(NA,16,8,16),x="center",bty="n")
dev.off()


##take 2 - so doesn't go below zero
pdf(paste(DIR$'Figures',"testingScalingObservationsMethod_take2.pdf",sep=""),height=7)
par(mar=c(3,4,1,1),oma=c(1,1,1,1),mfrow=c(5,2))
for(i in 1:4){
  testEst<-runif(11)*15
  testObs<-runif(11)*40
  
  #raw plot
  plot(testEst,type="l",xlab="Timestep",ylab="Abundance",ylim=c(0,50),lwd=2.5)
  points(testObs,pch=16,col=myGrey,cex=1)
  points(testObs,type="l",lwd=2,col=myGrey_trans)

  lmfit=nls(testEst~b+a*testObs,algorithm="port",start=c(a=1,b=0),lower=c(a=0,b=-Inf),upper=c(a=Inf,b=Inf))
  SS = sum(residuals(lmfit)^2)
  ScaledReal = abs(coef(lmfit)[1]) * testObs[!is.na(testObs)] + coef(lmfit)[2]
  
  points(ScaledReal,pch=8,col=myOrange,cex=1.5)
  points(ScaledReal,type="l",col=myOrange_trans,lwd=2)
  
  ##alt method
  #first scale by ratio of means
  meanObs<-mean(testObs); meanEst<-mean(testEst)
  step1Obs<-(testObs/meanObs)*meanEst ##now the means are the same
  # compare magnitude of change in obs and est
  maxObsDiff<-max(step1Obs)-min(step1Obs); maxEstDiff<-max(testEst)-min(testEst)
  thisScalar<-maxObsDiff/maxEstDiff
  
  newScaledO<-(step1Obs/thisScalar)
 
  plot(testEst,type="l",xlab="Timestep",ylab="Abundance",ylim=c(0,50),lwd=2.5)
  points(testObs,pch=16,col=myGrey,cex=1)
  points(testObs,type="l",lwd=2,col=myGrey_trans)
  
  points(newScaledO,col=myBlue_trans,pch=16,cex=1)
  points(newScaledO,col=myBlue_trans,type="l",lwd=2,cex=1.5)
  
}
makeBlankPlot()
legend(legend=c("Estimated","Observed","Observed (scaling 1)","Observed (scaling alt.)"),col=c("black",myGrey_trans,myOrange,myBlue_trans),
       lwd=c(2,2,2,2),pch=c(NA,16,8,16),x="center",bty="n")
dev.off()


##take 2 - so doesn't go below zero
pdf(paste(DIR$'Figures',"testingScalingObservationsMethod_samePlot.pdf",sep=""),height=7)
par(mar=c(3,4,1,1),oma=c(1,1,1,1),mfrow=c(5,2))
for(i in 1:4){
  testEst<-runif(11)*15
  testObs<-runif(11)*40
  
  #raw plot
  plot(testEst,type="l",xlab="Timestep",ylab="Abundance",ylim=c(0,45),lwd=2.5)
  points(testObs,pch=16,col=myGrey,cex=1)
  points(testObs,type="l",lwd=2,col=myGrey_trans)
  
  lmfit=nls(testEst~b+a*testObs,algorithm="port",start=c(a=1,b=0),lower=c(a=0,b=-Inf),upper=c(a=Inf,b=Inf))
  SS = sum(residuals(lmfit)^2)
  ScaledReal = abs(coef(lmfit)[1]) * testObs[!is.na(testObs)] + coef(lmfit)[2]
  
  points(ScaledReal,pch=8,col=myOrange,cex=1.5)
  points(ScaledReal,type="l",col=myOrange_trans,lwd=2)
  
  ##alt method
  #first scale by ratio of means
  meanObs<-mean(testObs); meanEst<-mean(testEst)
  step1Obs<-(testObs/meanObs)*meanEst ##now the means are the same
  # compare magnitude of change in obs and est
  maxObsDiff<-max(step1Obs)-min(step1Obs); maxEstDiff<-max(testEst)-min(testEst)
  thisScalar<-maxObsDiff/maxEstDiff
  
  newScaledO<-(step1Obs/thisScalar)

  points(newScaledO,col=myBlue_trans,pch=16,cex=1)
  points(newScaledO,col=myBlue_trans,type="l",lwd=2,cex=1.5)
  
  ##plot same thing zoomed in
  plot(testEst,type="l",xlab="Timestep",ylab="Abundance",ylim=c(0,25),lwd=2.5)

  
  points(ScaledReal,pch=8,col=myOrange,cex=1.5)
  points(ScaledReal,type="l",col=myOrange_trans,lwd=2)
  
  
  points(newScaledO,col=myBlue_trans,pch=16,cex=1)
  points(newScaledO,col=myBlue_trans,type="l",lwd=2,cex=1.5)
  
  
  
}
makeBlankPlot()
legend(legend=c("Estimated","Observed","Observed (scaling 1)","Observed (scaling alt.)"),col=c("black",myGrey_trans,myOrange,myBlue_trans),
       lwd=c(2,2,2,2),pch=c(NA,16,8,16),x="center",bty="n")
dev.off()


##take 2 - so doesn't go below zero
pdf(paste(DIR$'Figures',"testingScalingObservationsMethod_alt3.pdf",sep=""),height=7)
par(mar=c(3,4,1,1),oma=c(1,1,1,1),mfrow=c(5,2))
for(i in 1:4){
  testEst<-runif(11)*15
  testObs<-runif(11)*40
  
  #raw plot
  plot(testEst,type="l",xlab="Timestep",ylab="Abundance",ylim=c(0,45),lwd=2.5)
  points(testObs,pch=16,col=myGrey,cex=1)
  points(testObs,type="l",lwd=2,col=myGrey_trans)
  
  lmfit=nls(testEst~b+a*testObs,algorithm="port",start=c(a=1,b=0),lower=c(a=0,b=-Inf),upper=c(a=Inf,b=Inf))
  ScaledReal = abs(coef(lmfit)[1]) * testObs[!is.na(testObs)] + coef(lmfit)[2]
  
  points(ScaledReal,pch=8,col=myOrange,cex=1.5)
  points(ScaledReal,type="l",col=myOrange_trans,lwd=2)
  
  ##alt method
  #first scale by ratio of means
  meanObs<-mean(testObs); meanEst<-mean(testEst)
  step1Obs<-(testObs/meanObs)*meanEst ##now the means are the same
  # compare magnitude of change in obs and est
  maxObsDiff<-max(step1Obs)-min(step1Obs); maxEstDiff<-max(testEst)-min(testEst)
  thisScalar<-maxObsDiff/maxEstDiff
  
  newScaledO<-(step1Obs/thisScalar)
  
  points(newScaledO,col=myBlue_trans,pch=16,cex=1)
  points(newScaledO,col=myBlue_trans,type="l",lwd=2,cex=1.5)
  
  #alt 3
  lmfit=nls(testEst~b+a*testObs,algorithm="port",start=c(a=1,b=0),lower=c(a=0,b=0),upper=c(a=Inf,b=0))
  ScaledReal2 = abs(coef(lmfit)[1]) * testObs[!is.na(testObs)] + coef(lmfit)[2]
  
  points(ScaledReal2,pch=16,col=myRed,cex=1.5)
  points(ScaledReal2,type="l",col=myRed_trans,lwd=2)
  
  
  ##plot same thing zoomed in
  plot(testEst,type="l",xlab="Timestep",ylab="Abundance",ylim=c(0,max(ScaledReal,ScaledReal2,newScaledO,testEst)),lwd=2.5)
  
  
  points(ScaledReal,pch=8,col=myOrange,cex=1.5)
  points(ScaledReal,type="l",col=myOrange_trans,lwd=2)
  
  
  points(newScaledO,col=myBlue_trans,pch=16,cex=1)
  points(newScaledO,col=myBlue_trans,type="l",lwd=2,cex=1.5)
  
  points(ScaledReal2,pch=16,col=myRed,cex=1.5)
  points(ScaledReal2,type="l",col=myRed_trans,lwd=2)
  
  
}
makeBlankPlot()
legend(legend=c("Estimated","Observed","Observed (scaling 1)","Observed (scaling mult. only)","Observed (scaling alt.)"),col=c("black",myGrey_trans,myOrange_trans,myRed_trans,myBlue_trans),
       lwd=c(2,2,2,2,2),pch=c(NA,16,8,16,16),x="center",bty="n")
dev.off()

