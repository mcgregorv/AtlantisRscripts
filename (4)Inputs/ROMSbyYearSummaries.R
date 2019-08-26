#read in multi year ROMS and plot summaries for temp, salt, and poss movement
basePath<-paste(DIR$'Base',"ATLANTISmodels\\inputs\\ROMS\\",sep="")

tempRoms<-nc_open(paste(basePath,"Chatham30_tempAll.nc",sep=""))
temperature<-ncvar_get(tempRoms,"temperature")

getCI<-function(x,i){
  z<-x[!is.na(x)]
  y<-round(length(z)*i)
  if(y<1){y<-1}
  yy<-sort(z)[y]
  return(yy)
}

lCIbyTS<-apply(temperature,3,getCI, i=0.025)
uCIbyTS<-apply(temperature,3,getCI, i=0.975)
mCIbyTS<-apply(temperature,3,getCI, i=0.5)

xx<-as.vector(temperature[,,1])
ntsPerYear<-length(xx[!is.na(xx)])

yMax<-max(uCIbyTS); yMin<-min(lCIbyTS)

xaxisAt<-c(seq(1,length(mCIbyTS),by=365*2),length(mCIbyTS)); xaxisLab<-round(xaxisAt/(2*365))+1996
modelYears<-xaxisLab; modelTimeSteps<-seq(1,length(mCIbyTS),by=1)/(2*365)+1996

plot(mCIbyTS,type="l",lwd=2,col=myBlue,ylim=c(0,yMax*1.1),xaxt="n")
axis(at=xaxisAt,labels=xaxisLab,side=1)
polygon(x=c(seq(1,length(mCIbyTS)),seq(length(mCIbyTS),1)), y=c(lCIbyTS,rev(uCIbyTS)),col=myBlue_trans,border=NA)
abline(v=xaxisAt,col=myGrey_trans)

# https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html


#read in SST obs to compare
sstPath<-paste(DIR$'Base',"data\\nutrients\\temperature\\",sep="")
sstData<-nc_open(paste(sstPath, "sstmnmean.nc",sep=""))
tempData<-ncvar_get(sstData,"sst")
# 1.0 degree latitude x 1.0 degree longitude global grid (180x360).
# 89.5N - 89.5S, 0.5E - 359.5E.
degLat<-rev(seq(-89.5,89.5,by=1)); degLon<-seq(0.5,359.5,by=1)
degMonths<-seq(1,dim(tempData)[3])
degYears<-unlist(lapply(degMonths,FUN=function(x){ceiling((x-1)/12)})) + 1981
tempYears<-unique(degYears); nTempYears<-length(tempYears)
xaxsAt<-degMonths[seq(1,length(degMonths),by=12)]; xaxsLab<-unlist(lapply(xaxsAt,FUN=function(x){ceiling((x-1)/12)})) + 1981
sstTimesteps<-degMonths/12 + 1981+10/12

indexLat<-degLat<=(-43) & degLat>(-45); indexLon<-degLon>173 & degLon<=186
thisTempData<-tempData[indexLon,indexLat,]

toPlot<-apply(thisTempData,3,mean,na.rm=TRUE)
sstMax<-max(toPlot,na.rm=TRUE)

##long term averages
longTermSSTdata<-nc_open(paste(sstPath,"sstltm1961_1990.nc",sep=""))

ltsst<-ncvar_get(longTermSSTdata,"sst")
this_ltsst<-ltsst[indexLon,indexLat,]
histMean<-mean(this_ltsst)
histMeanByMonth<-apply(this_ltsst,3,mean)
xx<-rep(histMeanByMonth,(length(sstTimesteps)-1)/12); nx<-length(xx); nadd<-length(sstTimesteps)-1-nx
repHistMonthMeans<-c(histMeanByMonth[1],xx,histMeanByMonth[1:nadd])

pdf(paste(plotPath,"SST_1981-2017.pdf",sep=""),height=4,width=6)
par(mar=c(3.5,5,1,1))
plot(x=sstTimesteps, y=toPlot,type="l", xaxt="n",lwd=2, ylim=c(0,max(yMax, sstMax)),ylab=expression("SST ("~''^o~"Celsius)"),xlab="",cex.axis=1.5,cex.lab=1.5)
axis(at=xaxsLab, labels = xaxsLab, side=1,cex.axis=1.5,cex.lab=1.5)
polygon(x=c(modelTimeSteps, rev(modelTimeSteps)), y=c(rep(-1,length(modelTimeSteps)),rep(20,length(modelTimeSteps))),col=myBlue_trans,border=NA)
dev.off()

toOut<-data.frame(cbind(sstTimesteps,toPlot)); colnames(toOut)<-c("Time","SST")
write.csv(toOut,paste(DIR$'Tables',"SSTmeansTimeSeries.csv",sep=""),row.names = FALSE)

points(x=modelTimeSteps,y=mCIbyTS,type="l",lwd=2,col=myBlue)
polygon(x=c(modelTimeSteps, rev(modelTimeSteps)), y=c(lCIbyTS,rev(uCIbyTS)),col=myBlue_trans,border=NA)

points(x=sstTimesteps, y=repHistMonthMeans,lwd=2,col=myOrange, type="l")


plot(x=repHistMonthMeans, y=toPlot,pch=20, col=myBlue)
points(x=seq(0,sstMax), y=seq(0,sstMax), col=myOrange,lty=2,type="l", lwd=2)

## get mean by month for 1982 to 2016
sstMeanByMonth<-rep(NA,12)

getMonth<-function(x){
  y<-trunc(x); xx<-round((x-y)*12) + 1
  return(xx)
}
getYear<-function(x){
  y<-trunc(x)
  return(y)
}
sstTimestepMonths<-unlist(lapply(sstTimesteps, getMonth))  ; sstTimestepYears<-unlist(lapply(sstTimesteps, getYear))  
yearColRamp<-colorRampPalette(colors=c(myYellow,myGreen,myAqua,myBlue))(nTempYears)
getYearCol<-function(x){
  yearIndex<-tempYears==x
  thisCol<-yearColRamp[yearIndex]
  return(thisCol)
}
sstTimeYearCols<-unlist(lapply(sstTimestepYears, getYearCol))


sstMeanByMonth<-tapply(toPlot,sstTimestepMonths,mean,na.rm=TRUE)
minSST<-min(toPlot); maxSST<-max(toPlot)


##index years 1996-2004
indexYears<-sstTimesteps>=1996 & sstTimesteps<=2004
romsSSTbyMonth<-tapply(toPlot[indexYears],sstTimestepMonths[indexYears],mean,na.rm=TRUE)
pdf(paste(plotPath,"SSTmeanByMonthCompare.pdf",sep=""),height=4,width=5)
par(mar=c(5,5,1,1))
plot(histMeanByMonth,type="l",lwd=3,col=myOrange,ylim=c(0,17), ylab= expression("SST ("~''^o~"Celsius)"),xlab="Month",cex.axis=1.5,cex.lab=1.5)# histMeanByMonth is mean by month from 1961-1990
points(sstMeanByMonth,type="l",lwd=2.5,col="black",lty=2) #sstMeanByMonth is mean by month for 1981-2017
points(romsSSTbyMonth,type="l",lwd=2.5,col=myBlue)
par(lend=1)
legend(legend=c("1961-1990", "1981-2017", "1996-2004"),col=c(myOrange,"black",myBlue),lty=c(1,2,1), lwd=c(3,2.5,2.5),x="bottomleft",bty="n",seg.len=3,cex=1.5)
dev.off()



plot(x=seq(1,12), y=sstMeanByMonth,type="l", ylim=c(minSST, maxSST))
points(x=sstTimestepMonths,y=toPlot,pch=20,col=sstTimeYearCols)
points(x=seq(1,12), y=sstMeanByMonth,type="l",lwd=2)
points(x=seq(1,12), y=histMeanByMonth,type="l",lty=2,col=myOrange, lwd=2)

diffFromMean<-0*toPlot
for(m in 1:12){
  diffFromMean[sstTimestepMonths==m]<-toPlot[sstTimestepMonths==m]-sstMeanByMonth[m]
}

plot(x=sstTimesteps, y=diffFromMean, pch=20,col=myBlue_trans, cex=1.5)
abline(h=0,col=myOrange,lwd=3,lty=2)
abline(v=c(1996,2005), col=myGrey_trans,lwd=2)




