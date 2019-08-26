library(ggmap)
library(rgdal)
library(gtable)
library(maps)
library(mapdata)
library(rgeos)
source(paste(DIR$'General functions',"\\get_first_number.R",sep=""))
source(paste(DIR$'General functions',"formatShape.R",sep=""))

this_run<-"base"
# this_run<-"simple"
# this_out<-"SimpleShort"
this_out<-paste("Base",sep="")
# this_out<-paste("Test",sep="")
# 
basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
outPath<-paste(basePath,this_run,"\\","output",this_out,"\\",sep="")

plotPath<-paste(basePath,"\\Figures\\",this_run,"\\",sep="")

skip_tsteps<-0 #min burn-in

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1900
fishingStartYear<-1900
modelStartYear<-1900

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
allTracers<-sort(names(ThisNC.nc$var))

this_nts<-dim(thisVol)[3] #number of timesteps
nbox<-dim(thisVol)[2]; nlayer<-dim(thisVol)[1]

# #read in shape file
shapeFile<-paste(DIR$'Base',"ATLANTISmodels\\inputs\\bgm\\CHAT30_LL",sep="")
sdata<-read.shapefile(shapeFile)
shape<-formatShape(shapeFile=shapeFile)
ns<-length(shape)
SpDF <- SpatialPolygonsDataFrame(shape,data.frame( z=1:ns, row.names=paste("P",seq(1,(ns)),sep="")))
labels<-seq(1,(ns))
plot(shape)
LABELpOS<-polygonsLabel(shape, labels = labels, cex=.1,doPlot=FALSE)
labeldf<-data.frame(cbind("x"=LABELpOS[1:ns],"y"=LABELpOS[(ns+1):(2*ns)]))

wcColors<-colorRampPalette(colors=c(myAqua,myBlue,"midnightblue"))((nlayer-1))
waterColumnColors<-data.frame(cbind("layer"=seq(0,(nlayer-2)),"color"=wcColors))
sedCol<-myGreen
#################
## Tracer to plot ##
###################
allTracers<-unique(sort(names(ThisNC.nc$var)))
numTracers<-allTracers[grep("Num",allTracers)]; nt<-length(numTracers)

NumbersArray<-array(0,dim=dim(thisVol))
DensityArray<-array(0,dim=dim(thisVol))

for(t in 1:nt){
  thisTracer<-numTracers[t]
  thisData<-ncvar_get(ThisNC.nc,thisTracer)
  NumbersArray<-NumbersArray+thisData
}
thisVol[thisVol==0]<-NA
DensityArray<-NumbersArray/thisVol

nts<-this_nts-skip_tsteps+1

  pdf(paste(plotPath,"Tracer3DwithMap_DENSITY.pdf",sep=""))
  par(mfrow=c(2,2))
  for(b in 1:nbox){
    temp<-DensityArray[,b,c(skip_tsteps:this_nts)]
    # temp<-NumbersArray[,b,]
    #get depth for this box
    thisDepth<-round(sum(thisDz[,b,1]-1))
    thisMax<-max(temp,na.rm=TRUE)
    if(thisMax>0){
      plot(x=0*as.double(temp[1,]),type="n",ylim=c(0,thisMax),ylab="Tracer abundance",xlab="Timestep")
      layerIndex<-rowSums(temp)>0;  layerIndex[nlayer]<-TRUE; layerIndex[is.na(layerIndex)]<-FALSE
      temp<-temp[layerIndex,]
      #last layer is sediment. Other layers are wter column and need to be reversed
      thisWCLayers<-rev(seq((nlayer-2),0)[layerIndex[1:(nlayer-1)]])
      thisNL<-length(thisWCLayers); thisWCcolors<-(waterColumnColors$color[match(thisWCLayers,waterColumnColors$layer)])
      #get last value
      if(length(dim(temp))==0){
        points(temp,col=sedCol,lwd=2,type="l",lty=2)
      }else{
        for(l in 1:thisNL){
          thisCol<-as.character(thisWCcolors[l])
          if(is.na(thisCol)){thisCol<-myGrey}
          points(temp[l,],col=thisCol,lwd=2,type="l")
        }
        points(temp[(thisNL+1),],col=sedCol,lwd=2,type="l",lty=2)
      }
      mtext(paste("Box ",b,sep=""),side=3,adj=0.5,outer=FALSE)
      # axis(at=fvalue,labels=fvalue,side=4)
      plot(shape)
      map('nzHires',add=TRUE,col="black",lwd=2)
      map.axes()
      for(plotB in 1:dim(labeldf)[1]){
        if(b== plotB){thisCol=myGreen}else{thisCol="white"}
        polygon(sdata$shp$shp[[plotB]]$points,col=thisCol)
      }
      mtext(paste("Depth= ",thisDepth," m",sep=""),side=3,adj=1,outer=FALSE)
    }
  }
  dev.off()

