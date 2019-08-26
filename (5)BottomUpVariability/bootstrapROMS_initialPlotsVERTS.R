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

version<-"A"

#read in the full ROMS files
ROMSpath<-paste(basePath,"inputs\\ROMS\\",sep="")

outPath<-paste(basePath,"base\\ROMSBootstrapOut\\",sep="")

plotPath<-paste(basePath,"Figures\\base\\ROMSBootstrap\\",version,"\\",sep="")


nruns<-20

nl<-6; nb<-30

#read in groups file
groupsDF<-read.csv(paste(basePath,"CRAM_Groups.csv",sep=""))
vertsDF<-groupsDF[groupsDF$NumCohorts>1,]

vertBiomassTracers<-unlist(lapply(vertsDF$Name,FUN=function(x){paste(str_trim(x,side="both"),"_N",sep="")}))
nvars<-length(vertBiomassTracers)
r=1
thisOutPath<-paste(outPath,"outputROMSBootstrap",version,r,"\\",sep="")
thisOutFile<-paste(thisOutPath,"output.nc",sep="")
ThisNC.nc<-nc_open(thisOutFile)
thisData<-ncvar_get(ThisNC.nc,testVar)
thisVol<-ncvar_get(ThisNC.nc,"volume")
nts<-dim(thisData)[3]

tsToPlot<-seq(10,50,by=10); ntsPlot<-length(tsToPlot)
#want the distribution of each, and the distribution of the cvs of each
storeBiomass<-array(NA,dim=c(nl,nb,ntsPlot,nvars,nruns))
storeCV<-array(NA,dim=c(nl,nb,ntsPlot,nvars))
##biomass first ##
for(r in 1:nruns){
  cat("\n",r,"--")
  thisNC.nc<-paste(outPath,"outputROMSBootstrap",version,r,"\\",sep="")
  thisVol<-ncvar_get(ThisNC.nc,"volume")
  this_nts<-dim(thisVol)[3]
  this_tsToPlot<-tsToPlot[tsToPlot<=this_nts]
  ts_index<-tsToPlot %in% this_tsToPlot
  for(v in 1:nvars){
    thisTracer<-vertBiomassTracers[v]
    cat(thisTracer," ")
    thisData<-ncvar_get(ThisNC.nc,thisTracer)
    if(length(dim(thisData))==3){
      temp<-(thisData*thisVol)[,,this_tsToPlot]
      storeBiomass[,,ts_index,v,r]<-temp
    }else{
      temp<-(thisData*thisVol[nl,,])[,this_tsToPlot]
      storeBiomass[nl,,ts_index,v,r]<-temp
    }
  }
}
## now CV ##
# storeBiomassVectorData<-array(NA,dim=c(nts,nruns,nvars))
storeBetweenRunCV<-array(NA,dim=c(nl,nb,ntsPlot,nvars)); 
for(v in 1:nvars){
  for(b in 1:nb){
    for(l in 1:nl){
      for(t in 1:ntsPlot){
        # b=2;l=5;t=2; v=1
        testVar<-vertBiomassTracers[v]
        # cat(testVar,",,")
        ##BETWEEN YEAR
        thisVector<-storeBiomass[l,b,t,v,]
        if(sum(thisVector,na.rm=TRUE)>0){
          #turn 0's into NA
          index<-thisVector==0; thisVector[index]<-NA
          thisvaraiance<-var(thisVector); thisMean<-mean(thisVector); thiscv<-sqrt(thisvaraiance)/thisMean
          storeBetweenRunCV[l,b,t,v]<-thiscv 
        }
      }
    } 
  }
}
##none of these have any variance??
##NUMBERS
numberTracers<-c()
for(v in 1:dim(vertsDF)[1]){
  thisName<-str_trim(vertsDF$Name[v],side="both")
  thisNumCohorts<-vertsDF$NumCohorts[v]
  thisNumTracers<-paste(thisName,seq(1,thisNumCohorts),"_Nums",sep="")
  numberTracers<-c(numberTracers,thisNumTracers)
}
numberTracers<-sort(unique(numberTracers))

testTracers<-numberTracers; nvars<-length(testTracers)

storeThisData<-array(NA,dim=c(nl,nb,ntsPlot,nruns,nvars))

for(r in 1:nruns){
  thisOutPath<-paste(outPath,"outputROMSBootstrap",version,r,"\\",sep="")
  thisOutFile<-paste(thisOutPath,"output.nc",sep="")
  ThisNC.nc<-nc_open(thisOutFile)
  thisVol<-ncvar_get(ThisNC.nc,"volume")
  this_nts<-dim(thisVol)[3]
  this_tsToPlot<-tsToPlot[tsToPlot<=this_nts]
  ts_index<-tsToPlot %in% this_tsToPlot; ts_otherIndex<-seq(1,this_nts) %in% this_tsToPlot
  for(v in 1:nvars){
    testVar<-testTracers[v]
    cat(testVar,",,")
    thisData<-ncvar_get(ThisNC.nc,testVar)
    if(length(dim(thisData))==3){
      storeThisData[,,ts_index,r,v]<-thisData[,,ts_otherIndex]
    } else{
      storeThisData[nl,,ts_index,r,v]<-thisData[,ts_otherIndex]
    }
    
  }
}

storeBetweenRunCV<-array(NA,dim=c(nl,nb,ntsPlot,nvars)); 
for(v in 1:nvars){
  for(b in 1:nb){
    for(l in 1:nl){
      for(t in 1:ntsPlot){
        # b=2;l=5;t=2; v=1
        testVar<-testTracers[v]
        # cat(testVar,",,")
        ##BETWEEN YEAR
        thisVector<-storeThisData[l,b,t,,v]
        #turn 0's into NA
        index<-thisVector==0; thisVector[index]<-NA
        thisvaraiance<-var(thisVector); thisMean<-mean(thisVector); thiscv<-sqrt(thisvaraiance)/thisMean
        storeBetweenRunCV[l,b,t,v]<-thiscv 
      }
    }
  } 
}

thisPlotPath<-paste(plotPath,"VERTS\\",sep="")
for(v in 1:nvars){
  thisVar<-testTracers[v]
  thisVec<-as.vector(storeBetweenRunCV[,,,v])
  thisMax<-max(thisVec,na.rm=TRUE)
  if(thisMax>0.01){
    cat(thisVar,"--")
    npoints<-length(thisVec)
    seqIndex<-seq(1,npoints)
    tsIndex<-ceiling(seqIndex/(nl*nb-1))
    
    xlabs<-unique(tsIndex)
    xats<-xlabs*(nl*nb)
    axisIndex<-pretty(seq(0,length(xats),length.out=3))
    
    jpeg(paste(thisPlotPath,"ROMS_CV_",thisVar,".jpg",sep=""))
    plot(thisVec,col=myBlue_trans,pch=20,xaxt="n",xlab="Year",ylab="CV",ylim=c(0,max(thisVec,na.rm=TRUE)))
    axis(at=xats[axisIndex],labels = xlabs[axisIndex],side=1)
    mtext(thisVar,side=3,adj=0)
    dev.off()
  }
  

}
