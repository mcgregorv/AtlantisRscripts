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

version<-"B"

daysTimeStep<-5
numStepsPerYear<-365/daysTimeStep

#read in the full ROMS files
ROMSpath<-paste(basePath,"inputs\\ROMS\\",sep="")

outPath<-paste(basePath,"base\\ROMSBootstrapOut\\",sep="")

plotPath<-paste(basePath,"figures\\base\\ROMSBootstrap\\",version,"\\",sep="")

nruns<-50

nl<-6; nb<-30

r=1
thisOutPath<-paste(outPath,"outputROMSBootstrap",version,r,"\\",sep="")
thisOutFile<-paste(thisOutPath,"output.nc",sep="")
ThisNC.nc<-nc_open(thisOutFile)
thisVol<-ncvar_get(ThisNC.nc,"volume")
nts<-dim(thisVol)[3]

allTracers<-sort(names(ThisNC.nc$var))
skip<-c(grep("_Nums",allTracers),grep("_ResN",allTracers),grep("_StructN",allTracers))
allTracers[!(seq(1,length(allTracers)) %in% skip)]


NutrientTracers<-c("Chl_a", "Denitrifiction", "DON","NH3", "Nitrification", "NO3", "Oxygen", "Si", "salt","Temp", "vflux")
BacteraDetritusTracers<-c("Det_Si","Lab_Det_N", "Pelag_Bact_N", "Ref_Det_N", "Sed_Bact_N")
PrimaryProducerTracers<-c("Diatom_N","Diatom_S", "DinoFlag_N", "Macroalgae_N", "MicroPB_N", "MicroPB_S", "PicoPhytopl_N")
PrimaryConsumerTracers<-c("Benthic_Carniv_N", "Carniv_Zoo_N", "Deposit_Feeder_N","Meiobenth_N",  "Filter_Other_N", "Gelat_Zoo_N", "Invert_comm_Herb_N", "MicroZoo_N", "Zoo_N")
testTracers<-sort(unique(c(NutrientTracers, BacteraDetritusTracers, PrimaryConsumerTracers, PrimaryProducerTracers)))

nvars<-length(testTracers)

storeThisData<-array(NA,dim=c(nl,nb,nts,nruns,nvars))

for(r in 1:nruns){
  for(v in 1:nvars){
    testVar<-testTracers[v]
    cat(testVar,",,")
    thisOutPath<-paste(outPath,"outputROMSBootstrap",version,r,"\\",sep="")
    thisOutFile<-paste(thisOutPath,"output.nc",sep="")
    ThisNC.nc<-nc_open(thisOutFile)
    thisData<-ncvar_get(ThisNC.nc,testVar)
    if(length(dim(thisData))==3){
      storeThisData[,,,r,v]<-thisData
    }else{
      storeThisData[nl,,,r,v]<-thisData
    }
  }
  nc_close(ThisNC.nc)
}

meanBetweenRuns<-apply(storeThisData,c(1,2,3,5),mean,na.rm=TRUE)
varianceBetweenRuns<-apply(storeThisData,c(1,2,3,5),var,na.rm=TRUE)
cvBetweenRuns<-sqrt(varianceBetweenRuns)/meanBetweenRuns

boundaryBoxes<-c(1,seq(26,30))
boundaryBoxIndex<-seq(1,nb)%in% boundaryBoxes
ndynBoxes<-nb-length(boundaryBoxes)

thisMaxCV<-3.2

colorByBox<-c(myGrey,colorRampPalette(colors=c(myOrange,myRed,myPurple,myBlue,myAqua,myGreen))(ndynBoxes),rep(myGrey,5))

#using color By Box
for(v in 1:nvars){
  thisVar<-testTracers[v]
  thisCVs<-cvBetweenRuns[,,,v]
  thisNts<-dim(thisCVs)[3]
  thisMaxCV<-max(thisCVs,na.rm=TRUE)
  # pdf(paste(plotPath,"CVcoloredByBox_",thisVar,".pdf",sep=""))
  jpeg(paste(plotPath,"CVcoloredByBox_",thisVar,".jpg",sep=""))
  par(mar=c(5,5,1,1))
  plot(seq(1,thisNts),type="n",ylim=c(0,thisMaxCV),xlab="Day",ylab="CV",xaxt="n",cex.axis=1.5,cex.lab=1.5)
  axisAt<-pretty(seq(1,thisNts,length.out=10)); axisLab<-daysTimeStep*axisAt
  axis(at=axisAt,labels=axisLab,side=1,cex.axis=1.5)
  for(b in 1:nb){
    if(!b %in% boundaryBoxes){
      thisCol<-colorByBox[b]
      for(l in 1:nl){
        boxCVs<-thisCVs[l,b,]
        points(boxCVs[boxCVs>0],col=thisCol,pch=20)
      }
    }
  }
  dev.off()
}

#Plot map with legend colors
# jpeg(paste(plotPath,"mapWithColorByBox.jpg",sep=""))
pdf(paste(plotPath,"mapWithColorByBox.pdf",sep=""))
plot(shape)
map('nzHires',add=TRUE,col="black",lwd=2)
map.axes()
for(plotB in 1:dim(labeldf)[1]){
  thisCol<-colorByBox[plotB]
  polygon(sdata$shp$shp[[plotB]]$points,col=thisCol)
}
dev.off()

##now using color By Depth
r=1
thisOutPath<-paste(outPath,"outputROMSBootstrap",version,r,"\\",sep="")
thisOutFile<-paste(thisOutPath,"output.nc",sep="")
ThisNC.nc<-nc_open(thisOutFile)
thisDepthData<-ncvar_get(ThisNC.nc,"dz")
testDepth<-thisDepthData[,,1]
uniqueDepths<-sort(unique(round(colSums(testDepth))))-1
colorByDepth<-c(colorRampPalette(colors=c(myOrange,myRed,myPurple,myBlue,myAqua,myGreen))(length(uniqueDepths)-1),myGrey)
depthByCell<-0*thisDepthData[,,1]

for(b in 1:nb){
  for(l in 1:nl){
    if(l==nl){
      thisDepth<-round(sum(thisDepthData[,b,1]))
    }else{
      thisDepth<-round(sum(thisDepthData[l:(nl-1),b,1]))
    }
    depthByCell[l,b]<-thisDepth
  }
}
uniqueDepthRange<-cbind(c(0,uniqueDepths[1:(length(uniqueDepths)-1)]),uniqueDepths)

pdf(paste(plotPath,"mapWithColorByMaxDepth.pdf",sep=""))
plot(shape)
map('nzHires',add=TRUE,col="black",lwd=2)
map.axes()
for(plotB in 1:dim(labeldf)[1]){
  thisDepth<-round(max(depthByCell[,plotB],na.rm=TRUE))-1
  depthIndex<-uniqueDepthRange[,1]<=thisDepth & uniqueDepthRange[,2]>=thisDepth
  thisCol<-colorByDepth[depthIndex]
  polygon(sdata$shp$shp[[plotB]]$points,col=thisCol)
}
legend(legend=uniqueDepths,col=colorByDepth,lwd=5,x="bottom",ncol=4,bty="n",seg.len=1,title="Depth (m)")
dev.off()

#using color By Depth
for(v in 1:nvars){
  thisVar<-testTracers[v]
  thisCVs<-cvBetweenRuns[,,,v]
  thisNts<-dim(thisCVs)[3]
  thisMaxCV<-max(thisCVs,na.rm=TRUE)
  # pdf(paste(plotPath,"CVcoloredByBox_",thisVar,".pdf",sep=""))
  jpeg(paste(plotPath,"CVcoloredByDepth_",thisVar,".jpg",sep=""))
  par(mar=c(5,5,1,1))
  plot(seq(1,thisNts),type="n",ylim=c(0,thisMaxCV),xlab="Day",ylab="CV",xaxt="n",cex.axis=1.5,cex.lab=1.5)
  axisAt<-pretty(seq(1,thisNts,length.out=10)); axisLab<-daysTimeStep*axisAt
  axis(at=axisAt,labels=axisLab,side=1,cex.axis=1.5)
  for(b in 1:nb){
    if(!b %in% boundaryBoxes){
      for(l in 1:nl){
        thisDepth<-round(depthByCell[l,b])
        thisCol<-colorByDepth[match(thisDepth,uniqueDepths)]
        boxCVs<-thisCVs[l,b,]
        points(boxCVs[boxCVs>0],col=thisCol,pch=20)
      }
    }
  }
  dev.off()
}

#depth layers by color (legend)
pdf(paste(plotPath,"ColorByDepthLegend.pdf",sep=""),height=2,width=1)
par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(0,type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
par(xpd=TRUE,lend=1)
legend(legend=uniqueDepths, col=colorByDepth,lwd=7,x="center",bty="n",title="Depth (m)")
dev.off()
################################################################
##################################################################
##################################################################
