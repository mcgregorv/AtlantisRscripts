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
temp<-allTracers[!(seq(1,length(allTracers)) %in% skip)]
testTracers<-temp[grep("_N",temp)]

nvars<-length(testTracers)

storeNData<-array(NA,dim=c(nl,nb,nts,nruns,nvars))

for(r in 1:nruns){
  cat("\n",r,"--")
  thisOutPath<-paste(outPath,"outputROMSBootstrap",version,r,"\\",sep="")
  thisOutFile<-paste(thisOutPath,"output.nc",sep="")
  ThisNC.nc<-nc_open(thisOutFile)
  for(v in 1:nvars){
    testVar<-testTracers[v]
    cat(testVar,",,")
    thisData<-ncvar_get(ThisNC.nc,testVar)
    if(length(dim(thisData))==3){
      storeNData[,,,r,v]<-thisData*thisVol
    }else{
      storeNData[nl,,,r,v]<-thisData*thisVol[nl,,]
    }
  }
  nc_close(ThisNC.nc)
}

totalNByRun<-apply(storeNData,c(1,2,3,4),sum,na.rm=TRUE)

meanNBetweenRuns<-apply(totalNByRun,c(1,2,3),mean,na.rm=TRUE)
varianceNBetweenRuns<-apply(totalNByRun,c(1,2,3),var,na.rm=TRUE)
cvNBetweenRuns<-sqrt(varianceNBetweenRuns)/meanNBetweenRuns

##set up for color by depth
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
# shapeFile<-paste(DIR$'Base',"ATLANTISmodels\\inputs\\bgm\\CHAT30_LL",sep="")
# sdata<-read.shapefile(shapeFile)
# shape<-formatShape(shapeFile=shapeFile)
# ns<-length(shape)
# SpDF <- SpatialPolygonsDataFrame(shape,data.frame( z=1:ns, row.names=paste("P",seq(1,(ns)),sep="")))
# labels<-seq(1,(ns))
# plot(shape)
# LABELpOS<-polygonsLabel(shape, labels = labels, cex=.1,doPlot=FALSE)
# labeldf<-data.frame(cbind("x"=LABELpOS[1:ns],"y"=LABELpOS[(ns+1):(2*ns)]))

# pdf(paste(plotPath,"mapWithColorByMaxDepth.pdf",sep=""))
# plot(shape)
# map('nzHires',add=TRUE,col="black",lwd=2)
# map.axes()
# for(plotB in 1:dim(labeldf)[1]){
#   thisDepth<-round(max(depthByCell[,plotB],na.rm=TRUE))-1
#   depthIndex<-uniqueDepthRange[,1]<=thisDepth & uniqueDepthRange[,2]>=thisDepth
#   thisCol<-colorByDepth[depthIndex]
#   polygon(sdata$shp$shp[[plotB]]$points,col=thisCol)
# }
# legend(legend=uniqueDepths,col=colorByDepth,lwd=5,x="bottom",ncol=4,bty="n",seg.len=1,title="Depth (m)")
# dev.off()

#using color By Depth
boundaryBoxes<-c(1,26:30)

thisCVs<-cvNBetweenRuns
thisNts<-dim(thisCVs)[3]
thisMaxCV<-max(thisCVs,na.rm=TRUE)
jpeg(paste(plotPath,"CVcoloredByDepth_Ntotal.jpg",sep=""))
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


# #depth layers by color (legend)
# pdf(paste(plotPath,"ColorByDepthLegend.pdf",sep=""),height=2,width=1)
# par(mar=c(0,0,0,0),oma=c(0,0,0,0))
# plot(0,type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
# par(xpd=TRUE,lend=1)
# legend(legend=uniqueDepths, col=colorByDepth,lwd=7,x="center",bty="n",title="Depth (m)")
# dev.off()
  
###################################################
##its a bit hard to see the lower depths. How do they look on their own?
thisYmax<-max(cvNBetweenRuns,na.rm=TRUE)
##where is this max??
maxIndex<-cvNBetweenRuns>0.5; 
for(t in 1:dim(maxIndex)[3]){
  thisMax<-max(cvNBetweenRuns[,,t],na.rm=TRUE)
  if(thisMax>0.5){
    cat(t,"\n")
  }
}
boxNumberLookup<-t(array(rep(seq(1,nb),nl),dim=c(nb,nl)))

plotDepths<-c(uniqueDepths,(uniqueDepths+1))
sameScale=FALSE
cvMaxLimit<-0.1
for(thisDepth in plotDepths){
  depthIndex<-apply(depthByCell,2,FUN=function(x){round(x)==thisDepth})
  thisBoxes<-boxNumberLookup[depthIndex]
  test1<-cvNBetweenRuns[,,1][depthIndex]
  thisNCVs<-array(NA,dim=c(length(test1),dim(cvNBetweenRuns)[3]))
  for(t in 1:dim(cvNBetweenRuns)[3]){
    thisNCVs[,t]<-cvNBetweenRuns[,,t][depthIndex]
  }
  thisColorRamp<-colorRampPalette(colors=c(myGold,myRed,myPurple,myBlue,myAqua,myGreen))(length(test1))

  iIndex<-rowSums(thisNCVs,na.rm=TRUE)>0
  mfrowCol<-3
  ni<-length(seq(1,length(test1))[iIndex]); mfrowRow<-ceiling(ni/mfrowCol)
  if(mfrowRow==1){
    mfrowCol<-ni
  }
  if(mfrowRow>0){
    # jpeg(paste(plotPath,"CV_NTotal_byBoxForGivenDepth_",thisDepth,"_sameScale.jpg",sep=""))
    # par(mfrow=c(mfrowRow,mfrowCol),mar=c(4,4,0.5,0.5))
    for(i in 1:length(test1)){
      thisCol<-thisColorRamp[i]
      if(sum(thisNCVs[i,],na.rm=TRUE)>0){
        thisiYmax<-max(thisNCVs[i,])
        if(thisiYmax>=cvMaxLimit){
          if(!sameScale){
            thisYmax<-thisiYmax
            jpeg(paste(plotPath,"CV_NTotal_byBoxForGivenDepth_",thisDepth,"_box",thisBoxes[i],".jpg",sep=""))
          }else{
            jpeg(paste(plotPath,"CV_NTotal_byBoxForGivenDepth_",thisDepth,"_sameScale_box",thisBoxes[i],".jpg",sep=""))
          }
          plot(thisNCVs[i,],col=thisCol,xlab="",ylab="",ylim=c(0,thisYmax)) 
          mtext(paste("Box:",thisBoxes[i],sep=""),side=3,line=-1.5,adj=0.1)
          dev.off()
        }
      }
    }
    # dev.off()
  }
}

##do a map and color the polygons by max cv
colorByCV<-c(myGrey,colorRampPalette(colors=c(myGold,myGreen,myAqua,myBlue,myPurple,myRed,"red"))(10))
getCVColor<-function(cv,thisMax=1){
  x<-round(cv/thisMax,1)*10
  if(x>10){x<-10}
  if(x<1){x<-1}
  return(colorByCV[x])
}

legendCVMax<-0.53
pdf(paste(plotPath,"mapWithColorByMaxCV_N.pdf",sep=""))
plot(shape)
map('nzHires',add=TRUE,col="black",lwd=2)
map.axes()
for(plotB in 1:dim(labeldf)[1]){
  thisCVdata<-cvNBetweenRuns[,plotB,]
  thisMaxCV<-max(thisCVdata,na.rm=TRUE)
  thisCol<-getCVColor(thisMaxCV,thisMax=legendCVMax)
  polygon(sdata$shp$shp[[plotB]]$points,col=thisCol)
}
legendNums<-pretty(seq(0,legendCVMax,length.out=10)); legendCols<-unlist(lapply(legendNums,getCVColor,thisMax=legendCVMax))
legendText<-c(paste("<",legendNums[2],sep=""),legendNums[2:(length(legendNums)-1)],paste("",legendNums[length(legendNums)],sep=""))
legend(legend=legendText,col=legendCols,lwd=5,x="bottom",ncol=4,seg.len=1,title="Maximum CV")
dev.off()

##need it on a log scale
getCVLogColor<-function(cv,thisMax=1){
  x<-round(log(thisMax)/log(cv),1)*10
  if(x>10){x<-10}
  if(x<1){x<-1}
  return(colorByCV[x])
}
legendCVMax<-0.53
pdf(paste(plotPath,"mapWithColorByLogMaxCV_N.pdf",sep=""))
plot(shape)
map('nzHires',add=TRUE,col="black",lwd=2)
map.axes()
for(plotB in 1:dim(labeldf)[1]){
  thisCVdata<-cvNBetweenRuns[,plotB,]
  thisMaxCV<-max(thisCVdata,na.rm=TRUE)
  thisCol<-getCVLogColor(thisMaxCV,thisMax=legendCVMax)
  polygon(sdata$shp$shp[[plotB]]$points,col=thisCol)
}
test<-seq(0.1,legendCVMax,length.out=10)
exp(pretty(log(test)))

temp<-seq(log(0.05),log(0.53),length.out=10)
legendNums<-exp(temp); legendCols<-unlist(lapply(legendNums,getCVLogColor,thisMax=legendCVMax))
legendText<-signif(legendNums,2)
legend(legend=legendText,col=legendCols,lwd=5,x="bottom",ncol=4,seg.len=1,title="Maximum CV")

dev.off()

##plot cv in Ntotal against CV in temperature for each depth

##plot cv by trophic level




