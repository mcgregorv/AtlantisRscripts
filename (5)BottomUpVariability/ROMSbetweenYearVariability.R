#How much do the ROMS vars vary between years?
#This one looks at temperature
basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
nlayer<-6; nboxes<-30

plotPath<-paste(DIR$'Figures',"ROMS\\",sep="")

year1<-1996

#read in the full ROMS files
ROMSpath<-paste(basePath,"inputs\\ROMS\\",sep="")

tempROMSFile<-paste(ROMSpath,"Chatham30_tempAll.nc",sep="")
hydroROMSFile<-paste(ROMSpath,"Chatham30_hydroAll.nc",sep="")
saltROMSFile<-paste(ROMSpath,"Chatham30_saltAll.nc",sep="")

fullHydro<-nc_open(hydroROMSFile)
fullTemp<-nc_open(tempROMSFile)
fullSalt<-nc_open(saltROMSFile)

vnames_salt<-names(fullSalt$var); vnames_temp<-names(fullTemp$var); vnames_hydro<-names(fullHydro$var)
saltData<-ncvar_get(fullSalt,"salinity"); tempData<-ncvar_get(fullTemp,"temperature")
nts<-dim(saltData)[3] #they are 12 hour time steps
hydroExchangeData<-ncvar_get(fullHydro,"exchange"); hydroDestBData<-ncvar_get(fullHydro,"dest_b"); hydroDestKData<-ncvar_get(fullHydro,"dest_k")

nyears<-nts/(365*2) #check it is a round number!
tsPerYear<-2*365
yearsVec<-seq(year1,(year1+nyears-1))
#split the data up by year
saltByYear<-NULL; tempByYear<-NULL; exchangeByYear<-NULL; destBByYear<-NULL;  destKByYear<-NULL
for(y in 1:nyears){
  startY<-(y-1)*tsPerYear+1; endY<-startY+tsPerYear-1
  # cat(paste(y, " has start year",startY,", endY ", endY," which in years are ",startY/(365*2)," ",endY/(365*2),"\n",sep=""))
  saltByYear[[y]]<-saltData[,,startY:endY]; tempByYear[[y]]<-tempData[,,startY:endY]; 
  exchangeByYear[[y]]<-hydroExchangeData[,,,startY:endY]
  destBByYear[[y]]<-hydroDestBData[,,,startY:endY]; destKByYear[[y]]<-hydroDestKData[,,,startY:endY]
}

transColors<-colorRampPalette(colors=c(myOrange_trans,myRed_trans,myPurple_trans,myBlue_trans,myGreen_trans))(nyears)

#############################
## TEMPERATURE
#############################
maxTemp<-max(tempData,na.rm=TRUE); minTemp<-min(tempData,na.rm=TRUE)
dummyVec<-0*as.vector(tempByYear[[1]])
plot(dummyVec,ylim=c(minTemp,maxTemp),type="n",ylab="Temperature ",xlab="")
jpeg(paste(plotPath,"tempCVByYear.jpg",sep=""))
par(mfrow=c(3,3))
for(y in 1:nyears){
  thisvector<-as.vector(tempByYear[[y]])
  thisCol<-transColors[y]
  plot(thisvector,type="l",lwd=3,col=thisCol,ylim=c(minTemp,maxTemp),ylab="Temperature",xlab="")
  mtext(yearsVec[y],side=3,adj=0)
}
dev.off()

tempVectors<-data.frame(matrix(NA,ncol=length(dummyVec),nrow=nyears))
for(y in 1:nyears){
  tempVectors[y,]<-as.vector(tempByYear[[y]])
}
varianceByPoint<-as.double(apply(tempVectors,2,var,na.rm=TRUE))

meanByPoint<-as.double(apply(tempVectors,2,mean,na.rm=TRUE))

cvByPoint<-sqrt(varianceByPoint)/meanByPoint

plot(varianceByPoint,col=myBlue_trans)

npoints<-length(cvByPoint)

seqIndex<-seq(1,npoints)
tsIndex<-ceiling(seqIndex/(nlayer*nboxes-1))
dayIndex<-round(tsIndex/2)==tsIndex/2

jpeg(paste(plotPath,"tempCVSummary.jpg",sep=""))
plot(cvByPoint,col=myBlue_trans,pch=20,xaxt="n",xlab="Day",ylab="CV")
xlabs<-unique(tsIndex/2)
xats<-xlabs*(nlayer*nboxes*2)
axisIndex<-pretty(seq(0,length(xats),length.out=3))
axis(at=xats[axisIndex],labels = xlabs[axisIndex],side=1)
mtext("Temperature",side=3,adj=0)
dev.off()

cvAsArray<-array(cvByPoint,dim=c(nlayer,nboxes,730))

testArray<-array(as.double(tempVectors[1,]),dim=c(nlayer,nboxes,730))

###########################
##SALINITY
################################
maxSalt<-max(saltData,na.rm=TRUE); minSalt<-min(saltData,na.rm=TRUE)
dummyVec<-0*as.vector(saltByYear[[1]])
plot(dummyVec,ylim=c(minSalt,maxSalt),type="n",ylab="Salinity ",xlab="")
jpeg(paste(plotPath,"saltCVByYear.jpg",sep=""))
par(mfrow=c(3,3))
for(y in 1:nyears){
  thisvector<-as.vector(saltByYear[[y]])
  thisvector[thisvector==0]<-NA
  thisCol<-transColors[y]
  plot(thisvector,type="l",lwd=3,col=thisCol,ylab="Salinity",xlab="")
  mtext(yearsVec[y],side=3,adj=0)
}
dev.off()

saltVectors<-data.frame(matrix(NA,ncol=length(dummyVec),nrow=nyears))
for(y in 1:nyears){
  saltVectors[y,]<-as.vector(saltByYear[[y]])
}
saltVarianceByPoint<-as.double(apply(saltVectors,2,var,na.rm=TRUE))

saltMeanByPoint<-as.double(apply(saltVectors,2,mean,na.rm=TRUE))

saltCvByPoint<-sqrt(saltVarianceByPoint)/saltMeanByPoint

jpeg(paste(plotPath,"SaltCVSummary.jpg",sep=""))
par(mfrow=c(1,1))
plot(saltCvByPoint,col=myBlue_trans,pch=20,xaxt="n",xlab="Day",ylab="CV")
xlabs<-unique(tsIndex/2)
xats<-xlabs*(nlayer*nboxes*2)
axisIndex<-pretty(seq(0,length(xats),length.out=3))
axis(at=xats[axisIndex],labels = xlabs[axisIndex],side=1)
mtext("Salinity",side=3,adj=0)
dev.off()

###########################
##EXCHANGE
################################
maxExchange<-max(hydroExchangeData,na.rm=TRUE); minExchange<-min(hydroExchangeData,na.rm=TRUE)
dummyVec<-0*as.vector(exchangeByYear[[1]])
plot(dummyVec,ylim=c(minExchange,maxExchange),type="n",ylab="Exchange ",xlab="")

jpeg(paste(plotPath,"ExchangeByYear.jpg",sep=""))
par(mfrow=c(3,3))
for(y in 1:nyears){
  thisvector<-as.vector(exchangeByYear[[y]])
  thisvector[thisvector==0]<-NA
  thisCol<-transColors[y]
  plot(thisvector,type="l",lwd=3,col=thisCol,ylab="Exchange",xlab="")
  mtext(yearsVec[y],side=3,adj=0)
}
dev.off()

exchangeVectors<-data.frame(matrix(NA,ncol=length(dummyVec),nrow=nyears))
for(y in 1:nyears){
  exchangeVectors[y,]<-as.vector(exchangeByYear[[y]])
}
exchangeVarianceByPoint<-as.double(apply(exchangeVectors,2,var,na.rm=TRUE))

exchangeMeanByPoint<-as.double(apply(exchangeVectors,2,mean,na.rm=TRUE))

exchangeCvByPoint<-sqrt(exchangeVarianceByPoint)/exchangeMeanByPoint

# pdf(paste(plotPath,"ExchangeCVSummary.pdf",sep=""),height=2,width=3)

jpeg(paste(plotPath,"ExchangeCVSummary.jpg",sep=""))
par(mfrow=c(1,1))
plot(exchangeCvByPoint,col=myBlue_trans,pch=20,xaxt="n",xlab="Day",ylab="CV",ylim=c(-500,500))
xlabs<-unique(tsIndex/2)
xats<-xlabs*(nlayer*nboxes*2)
axisIndex<-pretty(seq(0,length(xats),length.out=3))
# axis(at=xats[axisIndex],labels = xlabs[axisIndex],side=1)
mtext("Exchange",side=3,adj=0)
dev.off()

exchangeCVsArray<-array(exchangeCvByPoint,dim=c(21,nlayer,  nboxes, 730))












