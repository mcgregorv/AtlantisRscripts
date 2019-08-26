#How much do the ROMS vars vary between years?
#This one looks at temperature
basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
plotPath<-paste(basePath,"figures\\base\\ROMSBootstrap\\ROMSvariables\\",sep="")
nlayer<-6; nboxes<-30

year1<-1996

#read in the full ROMS files
ROMSpath<-paste(basePath,"inputs\\ROMS\\",sep="")

tempROMSFile<-paste(ROMSpath,"Chatham30_tempAll.nc",sep="")

fullTemp<-nc_open(tempROMSFile)

vnames_temp<-names(fullTemp$var); tempData<-ncvar_get(fullTemp,"temperature")
nts<-dim(tempData)[3] #they are 12 hour time steps

nyears<-nts/(365*2) #check it is a round number!
tsPerYear<-2*365
yearsVec<-seq(year1,(year1+nyears-1))
#split the data up by year
tempByYear<-NULL; 
for(y in 1:nyears){
  startY<-(y-1)*tsPerYear+1; endY<-startY+tsPerYear-1
  tempByYear[[y]]<-tempData[,,startY:endY]
}

transColors<-colorRampPalette(colors=c(myOrange_trans,myRed_trans,myPurple_trans,myBlue_trans,myGreen_trans))(nyears)

#for each time step, get the non-zero temperature average over space
nonZeroMean<-function(x){
  x[x==0]<-NA
  thisMean<-mean(x,na.rm=TRUE)
  return(thisMean)
}
xPointMovingAverage<-function(x,n){
  y<-0*x
  for(i in 1:length(x)){
    thisStart<-max(1,(i-n))
    thisVec<-x[thisStart:i]
    y[i]<-mean(thisVec,na.rm=TRUE)
  }
  return(y)
}
xMidPointMovingAverage<-function(x,n){
  y<-0*x
  for(i in 1:length(x)){
    thisStart<-max(1,(i-(trunc(n/2))))
    thisStop<-min(length(x),(i+trunc(n/2)))
    thisVec<-x[thisStart:thisStop]
    y[i]<-mean(thisVec,na.rm=TRUE)
  }
  return(y)
}

axistext<-pretty(seq(5,365,length.out=10))
axisat<-axistext*2

meanTempByDayYear<-array(NA,dim=c(730,nyears))
for(y in 1:nyears){
  thisYear<-year1+y-1
  thisTempData<-tempByYear[[y]]
  thisMeanTemp<-apply(thisTempData,3,nonZeroMean)
  meanTempByDayYear[,y]<-thisMeanTemp
  #do 10 ts moving average (5 day)
  pdf(paste(plotPath,"TemperatureMovingAverage",thisYear,".pdf",sep=""),height=4,width=4)
  par(mar=c(5,5,2,1))
  thisYlab <- 'Temperature degree~C'
  label1 <- 'degree~C'
  thisYlab<-parse(text = paste("Temperature ", "*degree ~ C", sep = ""))
  plot(thisMeanTemp,pch=20,cex=0.8,ylab=thisYlab,cex.axis=1.5,cex.lab=1.5,xlab="Day of year",xaxt="n")
  axis(at=axisat,labels=axistext,side=1,cex.axis=1.5)
  thisMovingMean<-xMidPointMovingAverage(x=thisMeanTemp,n=60)
  points(thisMovingMean,type="l",col=myBlue,lwd=2,lty=1)
  mtext(thisYear,side=3,adj=0,cex=1.5)
  dev.off()
}
#pop them all on one plot
jpeg(paste(plotPath,"TemperatureMovingAverage_allYears.jpg",sep=""))
par(mar=c(10,5,0.5,1))
plot(meanTempByDayYear[,1],pch=20,cex=0.8,ylab=thisYlab,cex.axis=1.5,cex.lab=1.5,xlab="Day of year",col=myGrey_trans,ylim=c(8.2,10.6),xaxt="n")
axis(at=axisat,labels=axistext,side=1,cex.axis=1.5)
for(y in 2:nyears){
  points(meanTempByDayYear[,y],pch=20,col=myGrey_trans)
}
for(y in 1:nyears){
  thisMovingMean<-xMidPointMovingAverage(x=meanTempByDayYear[,y],n=60)
  points(thisMovingMean,type="l",col=transColors[y],lwd=2)
}
par(xpd=TRUE)
legend(legend=seq(year1,(year1+nyears-1)),col=transColors,lwd=3,x="bottom",horiz=FALSE,inset=-0.4,ncol=round(nyears/3),bty="n")
dev.off()

##with overall trend line

tempMeanOverYears<-apply(meanTempByDayYear,1,nonZeroMean)
for(y in 1:nyears){
  thisYear<-year1+y-1
  pdf(paste(plotPath,"TemperatureMovingAverage_allYearsMean_",thisYear,".pdf",sep=""), width=4,height=4)
  par(mar=c(5,5,1.5,1))
  plot(meanTempByDayYear[,y],lwd=2.5,type="l",ylab=thisYlab,cex.axis=2,cex.lab=2,xlab="Day of year",col="midnightblue",ylim=c(8.2,10.6),xaxt="n")
  axis(at=axisat,labels=axistext,side=1,cex.axis=2)
  points(tempMeanOverYears,type="l",col=myGrey,lwd=2.5)
  mtext(thisYear,side=3,adj=0,cex=2)
  dev.off()
}
#legend
pdf(paste(plotPath,"TemperatureMovingAverage_allYearsMean_LEGEND.pdf",sep=""), width=5,height=1.5)
par(mar=c(0,0,0,0), lend=1)
makeBlankPlot()
legend(legend=c("Mean sea temperature by day over\n years 1996-2004","Mean sea temperature by day for each year"),
       col=c(myGrey,"midnightblue"), lty=c(1,1), lwd=c(3.5,3.5), seg.len=5, bty="n", x="center")
dev.off()


pdf(paste(plotPath,"TemperatureMovingAverage_allYearsMean_byYear.pdf",sep=""), height=9,width=7)
par(mar=c(5,5,1.5,1), mfrow=c(5,2))
tempMeanOverYears<-apply(meanTempByDayYear,1,nonZeroMean)
for(y in 1:nyears){
  thisYear<-year1+y-1
  plot(meanTempByDayYear[,y],lwd=2.5,type="l",ylab=thisYlab,cex.axis=1.5,cex.lab=1.5,xlab="Day of year",col="midnightblue",ylim=c(8.2,10.6),xaxt="n")
  axis(at=axisat,labels=axistext,side=1,cex.axis=1.5)
  points(tempMeanOverYears,type="l",col=myGrey,lwd=2.5)
  mtext(thisYear,side=3,adj=0,cex=1.5)

}
makeBlankPlot()
legend(legend=c("Mean sea temperature by day over\n years 1996-2004","Mean sea temperature by day and year"),
       col=c(myGrey,"midnightblue"), lty=c(1,1), lwd=c(2.5,2.5), bty="n", x="center")
dev.off()


#one plot
jpeg(paste(plotPath,"TemperatureMovingAverage_allYearsMean.jpg",sep=""))
tempMeanOverYears<-apply(meanTempByDayYear,1,nonZeroMean)
par(mar=c(10,5,0.5,1))
plot(meanTempByDayYear[,1],pch=20,cex=0.8,ylab=thisYlab,cex.axis=1.5,cex.lab=1.5,xlab="Day of year",col=myGrey_trans,ylim=c(8.2,10.6),xaxt="n")
axis(at=axisat,labels=axistext,side=1,cex.axis=1.5)
for(y in 2:nyears){
  points(meanTempByDayYear[,y],pch=20,col=myGrey_trans)
}
points(tempMeanOverYears,pch=20,col="midnightblue",lwd=0.5,cex=1)
thisMovingMean<-xMidPointMovingAverage(x=tempMeanOverYears,n=60)
points(thisMovingMean,type="l",col=myAqua,lwd=2)

par(xpd=TRUE)
legend(legend=seq(year1,(year1+nyears-1)),col=transColors,lwd=3,x="bottom",horiz=FALSE,inset=-0.4,ncol=round(nyears/3),bty="n")
dev.off()

#what does it look like when we don't have a dominant trend?
tsIndex<-seq(1,730)[seq(1,730)<200]; 
tsIndex<-seq(1,730)[seq(1,730)<=550 & seq(1,730)>=450]; 
ts1<-min(seq(1,730)[tsIndex]); ts2<-max(seq(1,730)[tsIndex])
jpeg(paste(plotPath,"TemperatureMovingAverage_Notrend",ts1,"-",ts2,".jpg",sep=""))
par(mar=c(10,5,0.5,1))
plot(meanTempByDayYear[tsIndex,1],pch=20,cex=0.8,ylab=thisYlab,cex.axis=1.5,cex.lab=1.5,xlab="Day of year",col=myGrey_trans,ylim=c(8.2,10.6),xaxt="n")
axis(at=axisat,labels=axistext,side=1,cex.axis=1.5)
for(y in 2:nyears){
  points(meanTempByDayYear[tsIndex,y],pch=20,col=myGrey_trans)
}
for(y in 1:nyears){
  thisMovingMean<-xMidPointMovingAverage(x=meanTempByDayYear[tsIndex,y],n=60)
  points(thisMovingMean,type="l",col=transColors[y],lwd=2)
}
par(xpd=TRUE)
legend(legend=seq(year1,(year1+nyears-1)),col=transColors,lwd=3,x="bottom",horiz=FALSE,inset=-0.4,ncol=round(nyears/3),bty="n")
dev.off()


jpeg(paste(plotPath,"TemperatureMovingAverage_allYearsMean_Nottrend",ts1,"-",ts2,".jpg",sep=""))
tempMeanOverYears<-apply(meanTempByDayYear,1,nonZeroMean)
par(mar=c(10,5,0.5,1))
plot(meanTempByDayYear[tsIndex,1],pch=20,cex=0.8,ylab=thisYlab,cex.axis=1.5,cex.lab=1.5,xlab="Day of year",col=myGrey_trans,ylim=c(8.2,10.6),xaxt="n")
axis(at=axisat,labels=axistext,side=1,cex.axis=1.5)
for(y in 2:nyears){
  points(meanTempByDayYear[tsIndex,y],pch=20,col=myGrey_trans)
}
points(tempMeanOverYears[tsIndex],pch=20,col="midnightblue",lwd=0.5,cex=1)
thisMovingMean<-xMidPointMovingAverage(x=tempMeanOverYears[tsIndex],n=60)
points(thisMovingMean,type="l",col=myAqua,lwd=2)

par(xpd=TRUE)
legend(legend=seq(year1,(year1+nyears-1)),col=transColors,lwd=3,x="bottom",horiz=FALSE,inset=-0.4,ncol=round(nyears/3),bty="n")
dev.off()

#how does it look when we bootstrap it?

#for each point, get the median, UQ and LQ
getUQ<-function(x){
  temp<-summary(x,na.rm=TRUE)
  return(temp[5])
}
getLQ<-function(x){
  temp<-summary(x,na.rm=TRUE)
  return(temp[2])
}
getUpperCI<-function(x){
  temp<-sort(x)
  index<-round(length(x)*0.975)
  return(temp[index])
}
getLowerCI<-function(x){
  temp<-sort(x)
  index<-round(length(x)*0.025)
  return(temp[index])
}

nsamples<-20; 
storeSampledYears<-array(NA,dim=c(730,nsamples))
for(s in 1:nsamples){
  thisSampleYear<-sample(seq(1,nyears),size=1)
  storeSampledYears[,s]<-meanTempByDayYear[,thisSampleYear]
}

allMedians<-apply(storeSampledYears,1,median,na.rm=TRUE)
allUQ<-apply(storeSampledYears,1,getUQ)
allLQ<-apply(storeSampledYears,1,getLQ)
allUpperCI<-apply(storeSampledYears,1,getUpperCI)
allLowerCI<-apply(storeSampledYears,1,getLowerCI)


jpeg(paste(plotPath,"TemperatureMovingAverage_allYearsMean_Bootstrap_n",nsamples,".jpg",sep=""))
par(mar=c(5,5,0.5,1))
plot(allMedians,type="l",ylab=thisYlab,cex.axis=1.5,cex.lab=1.5,xlab="Day of year",col=myGrey,ylim=c(8.2,10.6),xaxt="n",lwd=2.5)
axis(at=axisat,labels=axistext,side=1,cex.axis=1.5)
points(allUpperCI,type="l",col=myGrey,lwd=2,lty=2)
points(allLowerCI,type="l",col=myGrey,lwd=2,lty=2)
dev.off()
