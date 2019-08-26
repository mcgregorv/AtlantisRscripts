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
source(paste(DIR$'General functions',"myRounding.R",sep=""))

getColor<-function(x,thisMax){
  thisCol<-"white"
  if(!is.na(x) & x>0){
    if(x<thisMin){x<-thisMin}
    y<-round((x-thisMin)/(thisMax-thisMin),2)*100+1
    thisCol<-thisColRamp[y]
  }
  return(thisCol)
}
plotGrid<-function(x,y){
  thisX<-c(x-0.5,x-0.5,x+0.5,x+0.5); thisY<-c(y-0.5,y+0.5,y+0.5,y-0.5)
  thisCol<-plotColour[x,y]
  if(length(thisCol)>0){
    polygon(x=thisX,y=thisY,col=thisCol,border=NA)
  }
  return(NULL)
}

basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

Version<-"D"
Version<-"B"
Version<-"shortB"
plotPath<-paste(DIR$'Figures',"ROMS\\testing",Version,sep="")

# outPath<-paste(basePath,"base\\ROMS1yearRepeat\\",sep="")
outPath<-paste(basePath,"base\\ouputROMS",Version,"\\",sep="")
outPath<-paste(basePath,"base\\outputROMS",Version,"\\",sep="")
dataOutPath<-paste(basePath,"BootstrapROMS\\modelOut",Version,"_", sep="")

load(paste(dataOutPath,"modelTracers",sep="")); ##to bring "to bring "storeTracers", "storeTracersByCell", "storeTemperature", "store_nts".

#get _N tracers
BaseNC.nc<-nc_open(paste(outPath, "outputROMSBootstrap",Version,"1\\output.nc",sep=""))
allTracers<-names(BaseNC.nc$var)
thisVol<-ncvar_get(BaseNC.nc, "volume"); nlayers<-dim(thisVol)[1]; nts<-dim(thisVol)[3]
thisDZ<-ncvar_get(BaseNC.nc,"dz")
x<-grep("_N",allTracers); temp<-allTracers[x]
y<-grep("Nums",temp, invert = TRUE); Ntracers<-temp[y]; ntracers<-length(Ntracers)

## check out temperature - how much it differs between boxes and layers and time, and how much it varies between runs
## may as well grab the rest too
allNutTracers<-sort(c("Temp", "Light", "Si", "Nitrification", "NH3", "NO3", "Oxygen", "Chl_a", "eflux", "vflux", "salt"))
allNutTracers<-sort(allTracers[grep("_N|_StructN|_ResN",allTracers,invert = TRUE)])

storeNutTracers<-array(NA, dim=c(dim(storeTracersByCell)[1], length(allNutTracers), dim(storeTracersByCell)[3:5]))

nruns<-dim(storeNutTracers)[1]
ntracers<-length(allNutTracers)
nlayers<-dim(thisVol)[1]; nboxes<-dim(thisVol)[2]

for(r in 1:nruns){
  thisOutPath<-paste(outPath, "outputROMSBootstrap",Version,r,"\\",sep="")
  thisOutFile<-paste(thisOutPath,"output.nc",sep="")
  ThisNC.nc<-nc_open(thisOutFile)
  thisVol<-ncvar_get(ThisNC.nc,"volume")
  for(t in 1:length(allNutTracers)){
    thisTracer<-allNutTracers[t]
    thisData<-ncvar_get(ThisNC.nc, thisTracer)
    if(length(dim(thisData))==3){
      storeNutTracers[r,t,,,]<-thisData
    } else{
      storeNutTracers[r,t,nlayers,,]<-thisData
    }
  }
}

storeCVsByCellAndTracer<-array(NA, dim=c(dim(storeNutTracers)[-1]))

for(t in 1:ntracers){
  testData<-storeTracersByCell[,t,,,]
  testCVs<-array(NA, dim=dim(testData)[-1])
  for(l in 1:dim(testCVs)[1]){
    for(b in 1:dim(testCVs)[2]){
      for(y in 1:dim(testCVs)[3]){
        this_x<-testData[,l,b,y]
        thisVar<-var(this_x); thisMean<-mean(this_x); thisCV<-max(0,sqrt(thisVar)/thisMean)
        storeCVsByCellAndTracer[t,l,b,y]<-thisCV
      }
    }
  }
}

t=grep("Temp", allNutTracers)

thisData<-storeCVsByCellAndTracer[t,,,]
thisMeanCVByBox<-apply(thisData,2,nonZeroMean)

nonZeroMin<-function(x){
  x[x==0]<-NA
  thisMin<-min(x,na.rm=TRUE)
  return(thisMin)
}
dynBoxes<-seq(2,25)

pdf(paste(plotPath, "nutrientCVsByLayer.pdf",sep=""))
par(mfrow=c(3,2), mar=c(4,4,1,1), oma=c(1,1,4,1))
for(t in 1:ntracers){
   for(l in 1:nlayers){
    thisData<-storeNutTracers[,t,l,dynBoxes,]
    thisMeanByBox<-apply(thisData,2,nonZeroMean)
    thisMinByBox<-apply(thisData, 2, nonZeroMin)
    thisMaxByBox<-apply(thisData, 2, max, na.rm=TRUE)
    thisMin<-min(thisMinByBox, na.rm=TRUE); thisMax<-max(thisMaxByBox, na.rm=TRUE)
    if(thisMax>0){
      plot(thisMeanByBox, ylim=c(thisMin, thisMax), pch=20)
      for(b in 1:length(thisMinByBox)){
        segments(x0=b, y0=thisMinByBox[b], x1=b, y1=thisMaxByBox[b], col=myBlue)
      }
    }else{
      makeBlankPlot()
    }
    mtext(paste("Layer ",l,sep=""), side=3,adj=0)
  }
  mtext(allNutTracers[t], side=3, outer=TRUE, line=0)
}
dev.off()

## is it all runs that have high det-si in wc..?
daysPerTs<-5; timeStepYears<-seq(0,(nts-1))*5/365
axisLabs<-pretty(seq(1,max(timeStepYears), length.out=5))
asixAt<-seq(1,nts)[timeStepYears %in% axisLabs]
t=grep("Det",allNutTracers)
l<-5
pdf(paste(plotPath,"DetSi_surfaceLayerByBox.pdf", sep=""))
par(mfrow=c(3,2))
for(thisB in 2:25){
  thisData<-storeNutTracers[,t,l,thisB,]
  thisMax<-max(thisData, na.rm=TRUE)
  plot(thisData[1,], pch=20, col=myGrey_trans, ylim=c(0, thisMax), xlab="Year", ylab="DetSi", xaxt="n")
  for(r in 1:nruns){
    points(thisData[r,], pch=20, col=myGrey_trans)
  }
  axis(at=axisAt, labels = axisLabs, side=1)
  mtext(thisB,side=3)
}
dev.off()

l<-6
pdf(paste(plotPath,"DetSi_SedimentByBox.pdf", sep=""))
par(mfrow=c(3,2))
for(thisB in 2:25){
  thisData<-storeNutTracers[,t,l,thisB,]
  thisMax<-max(thisData, na.rm=TRUE)
  plot(thisData[1,], pch=20, col=myGrey_trans, ylim=c(0, thisMax), xlab="timestep", ylab="DetSi")
  for(r in 1:nruns){
    points(thisData[r,], pch=20, col=myGrey_trans)
  }
  mtext(thisB,side=3)
}
dev.off()

## salt
t=grep("salt",allNutTracers)
for(l in 1:6){
thisData<-storeNutTracers[,t,l,,]
test<-apply(thisData, c(2,3), nonZeroMean)
par(mfrow=c(1,1))
plot(test[25,], xaxt="n", ylim=c(min(test, na.rm=TRUE), max(test, na.rm=TRUE))); axis(at=axisAt, labels=axisLabs, side=1); mtext(l,side=3)
for(b in 1:30){
  thisCol<-rainbow(30)[b]
  if(b %in% c(1,26:30)){thisCol=myGrey_trans}
  points(test[b,], col=thisCol)
}
points(test[25,], col="black")
points(test[24,], col="midnightblue")
}

plot(these, xaxt="n", ylim=c(6000, 10400))
axis(at=axisAt, labels=axisLabs, side=1)
points(test[24,],col="red")
points(theothers,col=myBlue)

maxByRun<-apply(thisData, 1, nonZeroMean)
