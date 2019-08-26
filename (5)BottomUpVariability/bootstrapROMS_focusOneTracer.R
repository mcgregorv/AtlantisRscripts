#the ROMS variables we have are for XX years
#this script do some initial plots.
library(ggmap)
library(rgdal)
library(gtable)
library(maps)
library(mapdata)
library(rgeos)
source(paste(DIR$'General functions',"\\get_first_number.R",sep=""))
source(paste(DIR$'General functions',"formatShape.R",sep=""))

basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

thisTracer<-"Meiobenth_N"
thisTracer<-"Diatom_N"


version<-"A"

#read in the full ROMS files
ROMSpath<-paste(basePath,"inputs\\ROMS\\",sep="")

outPath<-paste(basePath,"base\\ROMSBootstrapOut\\",sep="")

nruns<-20

nl<-6; nb<-30

boundaryBoxes<-c(1,seq(26,30))

r=1; baseNC.nc<-nc_open(paste(outPath,"outputROMSBootstrap",version,r,"\\output.nc",sep=""))
baseVol<-ncvar_get(baseNC.nc,"volume"); base_nts<-dim(baseVol)[3]

storeThisData<-array(NA,dim=c(nl,nb,base_nts,nruns))
storeBiomassData<-0*storeThisData; storeDepthData<-0*storeThisData

for(r in 1:nruns){
  thisOutPath<-paste(outPath,"outputROMSBootstrap",version,r,"\\",sep="")
  thisOutFile<-paste(thisOutPath,"output.nc",sep="")
  ThisNC.nc<-nc_open(thisOutFile)
  thisData<-ncvar_get(ThisNC.nc,thisTracer); thisVol<-ncvar_get(ThisNC.nc,"volume")
  if(length(dim(thisData))==3){
    storeThisData[,,,r]<-thisData
    storeBiomassData[,,,r]<-thisData*thisVol
  }else{
    storeThisData[nl,,,r]<-thisData
    storeBiomassData[nl,,,r]<-thisData*(thisVol[nl,,])
  }
  thisDepthData<-ncvar_get(ThisNC.nc,"dz")
  storeDepthData[,,,r]<-thisDepthData
}
testDepth<-thisDepthData[,,1]
uniqueDepths<-sort(unique(colSums(testDepth)))-1
colorByDepth<-colorRampPalette(colors=c(myOrange,myRed,myPurple,myBlue,myAqua,myGreen))(length(uniqueDepths))
depthByCell<-0*thisDepthData[,,1]

for(b in 1:nb){
  for(l in 1:nl){
    if(l==nl){
      thisDepth<-sum(thisDepthData[,b,1])
    }else{
      thisDepth<-sum(thisDepthData[l:(nl-1),b,1])
    }
    depthByCell[l,b]<-thisDepth
  }
}
xALL<-as.vector(storeThisData)
x<-xALL[!is.na(xALL)]
h<-hist(x,col="white",border=NA,freq=FALSE)

thisTracerMax<-5e+3
jpg(paste(plotPath,"TestingBiomassByBoxLayerTime_",testVar,version,".jpg",sep=""))
plot(x=0,y=0,ylim=c(0,0.5),xlim=c(0,thisTracerMax),type="l",ylab="",xlab="")
for(b in 1:nb){
  for(l in 1:nl){
    if(l==nl){
      thisDepth<-depthByCell[l,b]-1;
    }else{
      thisDepth<-depthByCell[l,b]; 
    }
    thisCol<-colorByDepth[uniqueDepths==thisDepth]
    for(t in 1:base_nts){
      x<-storeThisData[l,b,t,]

        xvariance<-var(x,na.rm=TRUE)
        if(!is.na(xvariance>0)){
          if(xvariance>0){
            if(max(x,na.rm=TRUE)<5e+3){
            # h<-hist(x,freq=TRUE,col=thisCol)
            # mtext(paste("l ",l,", b ",b),side=3,adj=0)
              # cat(paste("l ",l,", b ",b," t ",t, "\n"),side=3,adj=0)
            h<-hist(x,plot=FALSE)
            xfit<-seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length=40)
            yfit<-dnorm(xfit,mean=mean(x,na.rm=TRUE),sd=sd(x,na.rm=TRUE))
            # yfit <- yfit*diff(h$mids[1:2])*length(x)
            # par(new=TRUE)
            points(xfit, yfit, col=thisCol, lwd=2,type="l")
          }
        }
      }
  
    }
  }
}
legend(legend=round(uniqueDepths),col=colorByDepth,lty=1,lwd=2,x="topright")
dev.off()

pdf(paste(plotPath,"TestingBiomassByBoxLayerTime_",testVar,version,".pdf",sep=""))
par(mfrow=c(3,2))
for(b in 1:nb){
  for(l in 1:nl){
    if(l==nl){
      thisDepth<-depthByCell[l,b]-1;
    }else{
      thisDepth<-depthByCell[l,b]; 
    }
    thisCol<-colorByDepth[uniqueDepths==thisDepth]
    for(t in 1:base_nts){
      x<-storeThisData[l,b,t,]
      xvariance<-var(x,na.rm=TRUE)
      if(!is.na(xvariance>0)){
        if(xvariance>0){
          h<-hist(x,freq=TRUE,col=thisCol)
          mtext(paste("l ",l,", b ",b),side=3,adj=0)
          h<-hist(x,plot=FALSE)
          xfit<-seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length=40)
          yfit<-dnorm(xfit,mean=mean(x,na.rm=TRUE),sd=sd(x,na.rm=TRUE))
          yfit <- yfit*diff(h$mids[1:2])*length(x)
          par(new=TRUE)
          plot(xfit, yfit, col=thisCol, lwd=2,type="n",ylab="",xlab="")

          # # points(xfit, yfit, col=thisCol, lwd=2,type="l")
        }
      }
    }
  }
}
dev.off()

meanBetweenRuns<-apply(storeBiomassData,c(1,2,3),mean,na.rm=TRUE)
varianceBetweenRuns<-apply(storeBiomassData,c(1,2,3),var,na.rm=TRUE)
cvBetweenRuns<-sqrt(varianceBetweenRuns)/meanBetweenRuns

cvMax<-max(cvBetweenRuns,na.rm=TRUE)
cvByLayer<-apply(cvBetweenRuns,1,mean,na.rm=TRUE)
cvByBoxLayer<-signif(apply(cvBetweenRuns,c(1,2),max,na.rm=TRUE),2)
meanByBoxLayer<-apply(meanBetweenRuns,c(1,2),mean,na.rm=TRUE)
depthByBox<-colSums(storeDepthData[,,1,1],na.rm=TRUE)
areaByBox<-ncvar_get(baseNC.nc,"volume")[nl,,1]

boundaryBoxIndex<-seq(1,nb)%in% boundaryBoxes
plot(x=depthByBox[!boundaryBoxIndex],y=cvByBoxLayer[6,!boundaryBoxIndex])
plot(x=areaByBox[!boundaryBoxIndex],y=cvByBoxLayer[6,!boundaryBoxIndex])

jpeg(paste(plotPath,"ROMSBootstrapCVByBoxDepth_",testVar,version,".jpg",sep=""))
plot(x=depthByBox[!boundaryBoxIndex],y=cvByBoxLayer[6,!boundaryBoxIndex],col=myBlue,pch=2,lwd=2,xlab="Box depth (m)",ylab="CV",cex=1.5,cex.lab=1.5,cex.axis=1.5)
dev.off()

#plot cvs between runs for each var
 thisVec<-as.vector(cvBetweenRuns[,!boundaryBoxIndex,])
  
  npoints<-length(thisVec)
  seqIndex<-seq(1,npoints)
  tsIndex<-ceiling(seqIndex/(nlayer*ndynBoxes-1))
  
  xlabs<-unique(tsIndex)
  xats<-xlabs*(nlayer*ndynBoxes)
  axisIndex<-pretty(seq(0,length(xats),length.out=3))
  
  jpeg(paste(plotPath,"ROMS_CV_",version,"_",thisVar,".jpg",sep=""))
  plot(thisVec,col=myBlue_trans,pch=20,xaxt="n",xlab="Year",ylab="CV",ylim=c(0,max(thisVec,na.rm=TRUE)),cex.axis=1.5,cex.lab=1.5,lwd=2,cex=1.5)
  axis(at=xats[axisIndex],labels = xlabs[axisIndex],side=1,cex.axis=1.5)
  mtext(thisVar,side=3,adj=0,cex=1.5)
  dev.off()
  
#get biomass for each run and see if any appear stable just looking at that
lastPoint<-12
r=1
thisData<-storeBiomassData[,,,r]
thisY<-apply(thisData*(baseVol[,,]),3,sum,na.rm=TRUE)
plot(thisY[1:lastPoint],type="l",col=rainbow(nruns)[r],ylim=c(0,2.2e+22))
for(r in 1:nruns){
  thisData<-storeBiomassData[,,,r]
  thisY<-apply((thisData[,,])*(baseVol[,,]),3,sum,na.rm=TRUE)
  cat(max(thisY[1:lastPoint],na.rm=TRUE),"\n")
  points(thisY[1:lastPoint],type="l",col=rainbow(nruns)[r])
} 

box19<-thisData[,19,]

basePL<-ncvar_get(baseNC.nc,"Diatom_N")
for(b in 1:nb){
  test<-basePL[,b,]
  testMin<-min(test,na.rm=TRUE)
  if(testMin<0){
    cat(paste("box ",b ,"\n"))
  }
}