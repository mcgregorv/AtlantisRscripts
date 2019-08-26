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

version<-"A"

basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

#for reading in ROMS year permutations
ROMSpath<-paste(basePath,"inputs\\ROMS\\",sep="")

outPath<-paste(basePath,"base\\ROMS",version,"_output\\",sep="")

nruns<-50; nts<-21

nl<-6; nb<-30

#get _N tracers
allTracers<-names(ThisNC.nc$var)
x<-grep("_N",allTracers); temp<-allTracers[x]
y<-grep("Nums",temp, invert = TRUE); Ntracers<-temp[y]; ntracers<-length(Ntracers)

#set up array to store outputs
storeTracers<-array(NA, dim=c(nruns, ntracers,nts)); storeTracersByCell<-array(NA, dim=c(nruns, ntracers,nl,nb,nts))
for(r in 1:nruns){
  cat(r," -- ")
  thisOutPath<-paste(outPath,"outputROMSBootstrap",version,r,"\\",sep="")
  thisOutFile<-paste(thisOutPath,"output.nc",sep="")
  ThisNC.nc<-nc_open(thisOutFile)
  thisVol<-ncvar_get(ThisNC.nc,"volume")
  this_nts<-dim(thisData)[3]
  for(t in 1:ntracers){
    thisTracer<-Ntracers[t]; thisData<-ncvar_get(ThisNC.nc, thisTracer)
    if(length(dim(thisData))==3){
      xx<-apply(thisData*thisVol,3,sum)*mg_2_tonne*X_CN
      yy<-thisData*thisVol*mg_2_tonne*X_CN
      storeTracersByCell[r,t,,,]<-yy
    } else{
      xx<-apply(thisData*thisVol[nl,,],2,sum)*mg_2_tonne*X_CN
      yy<-thisData*thisVol[nl,,]*mg_2_tonne*X_CN
      storeTracersByCell[r,t,nl,,]<-yy
    }
    storeTracers[r,t,]<-xx
  }
}
storeTemperature<-array(NA,dim=c(nruns,nl,nb,nts)); thisTracer<-"Temp"
for(r in 1:nruns){
  cat(r," -- ")
  thisOutPath<-paste(outPath,"outputROMSBootstrap",version,r,"\\",sep="")
  thisOutFile<-paste(thisOutPath,"output.nc",sep="")
  ThisNC.nc<-nc_open(thisOutFile)
  thisVol<-ncvar_get(ThisNC.nc,"volume")
  thisData<-ncvar_get(ThisNC.nc, thisTracer)
  storeTemperature[r,,,]<-thisData
}

test<-apply(storeTracers,1,sum)
#colorByRun
runColors<-colorRampPalette(colors=c(myGold,"red",myRed,myPurple,myBlue,myAqua,myGreen,myDarkGreen))(nruns)

for(t in 1:ntracers){
  thisTracer<-Ntracers[t]; plotData<-storeTracers[,t,]; thisMax<-max(plotData)
  plot(plotData[1,],type="l", lwd=2,col=runColors[1], ylim=c(0,thisMax))
  for(r in 2:nruns){
    points(plotData[r,], type="l", lwd=2, col=runColors[r])
  }
  mtext(thisTracer,side=3,adj=0)
}

storeThisData<-array(NA,dim=c(nl,nb,nts,nruns,nvars))

storeBiomassVectorData<-array(NA,dim=c(nts,nruns,nvars))

for(v in 1:nvars){
  testVar<-testTracers[v]
  cat(testVar,",,")
  for(r in 1:nruns){
    thisOutPath<-paste(outPath,"outputROMSBootstrap",version,r,"\\",sep="")
    thisOutFile<-paste(thisOutPath,"output.nc",sep="")
    ThisNC.nc<-nc_open(thisOutFile)
    thisData<-ncvar_get(ThisNC.nc,testVar)
    if(length(dim(thisData))==3){
      storeThisData[,,,r,v]<-thisData
    }else{
      storeThisData[nl,,,r,v]<-thisData
    }
    if(length(dim(thisData))==3){
      storeBiomassVectorData[,r,v]<-apply(thisData*thisVol,3,sum,na.rm=TRUE)
    }
  }
}

storeBetweenRunVariance<-array(NA,dim=c(nl,nb,nts,nvars)); storeBetweenRunCV<-0*storeBetweenRunVariance
storeWithinRunVariance<-array(NA,dim=c(nl,nb,nruns,nvars)); storeWithinRunCV<-0*storeWithinRunVariance
for(v in 1:nvars){
  for(b in 1:nb){
    for(l in 1:nl){
      for(t in 1:nts){
        # b=2;l=5;t=2; v=1
        testVar<-testTracers[v]
        # cat(testVar,",,")
        ##BETWEEN YEAR
        thisVector<-storeThisData[l,b,t,,v]
        #turn 0's into NA
        index<-thisVector==0; thisVector[index]<-NA
        thisvaraiance<-var(thisVector); thisMean<-mean(thisVector); thiscv<-sqrt(thisvaraiance)/thisMean
        storeBetweenRunVariance[l,b,t,v]<-thisvaraiance
        storeBetweenRunCV[l,b,t,v]<-thiscv 
        if(!is.na(thiscv)){
          cat("l",l,", b",b," t",t,"\n")
        }

        
        # plot(thisVector,xlab="Run",ylab=testVar)
        # mtext(paste("CV=",signif(thiscv,2),", box ",b,":",l),side=3,adj=0)
      }
      # for(r in 1:nruns){
      #   ##WITHIN YEAR
      #   thisVector<-storeThisData[l,b,,r,v]
      #   #turn 0's into NA
      #   index<-thisVector==0; thisVector[index]<-NA
      #   thisvaraiance<-var(thisVector); thisMean<-mean(thisVector); thiscv<-sqrt(thisvaraiance)/thisMean
      #   storeWithinRunVariance[l,b,r,v]<-thisvaraiance
      #   storeWithinRunCV[l,b,r,v]<-thiscv 
      # }
    }
  } 
}

boundaryBoxes<-c(1,seq(26,30))
boundaryBoxIndex<-seq(1,nb)%in% boundaryBoxes
ndynBoxes<-nb-length(boundaryBoxes)
#plot cvs between runs for each var
lastYearToPlot<-10; lastPointToPlot<-lastYearToPlot*nlayer*ndynBoxes
for(v in 1:nvars){
  thisVar<-testTracers[v]
  thisVec<-as.vector(storeBetweenRunCV[,!boundaryBoxIndex,,v])
  
  npoints<-length(thisVec)
  seqIndex<-seq(1,npoints)
  tsIndex<-ceiling(seqIndex/(nlayer*ndynBoxes-1))
  
  xlabs<-unique(tsIndex)
  xats<-xlabs*(nlayer*ndynBoxes)
  axisIndex<-pretty(seq(0,length(xats),length.out=3))
  
  jpeg(paste(plotPath,"ROMS_CV_",version,"_",thisVar,".jpg",sep=""))
  plot(thisVec[1:lastPointToPlot],col=myBlue_trans,pch=20,xaxt="n",xlab="Year",ylab="CV",ylim=c(0,max(thisVec[1:lastPointToPlot],na.rm=TRUE)))
  # plot(thisVec,col=myBlue_trans,pch=20,xaxt="n",xlab="Year",ylab="CV",ylim=c(0,max(thisVec,na.rm=TRUE)))
  axis(at=xats[axisIndex],labels = xlabs[axisIndex],side=1)
  mtext(thisVar,side=3,adj=0)
  dev.off()
}
###############################################
# #TESTING TESTING
# test<-storeBetweenRunCV[,,1,v]
# test
# [1] 4.605451e-05
# > l=2; b=2; t=1; v=1
# > dim(storeThisData)
# [1]  6 30 12 20  1
# > storeThisData[l,b,t,,1]
# [1] 0.03005886 0.03006054 0.03005755 0.03005584 0.03005719 0.03005719 0.03005719 0.03005719 0.03006018 0.03005704 0.03006054 0.03005755 0.03005886 0.03006018
# [15] 0.03005755 0.03005974 0.03005886 0.03005886 0.03005755 0.03005719
# > x<-storeThisData[l,b,t,,1]
# > var(x)
# [1] 1.91634e-12
# > mean(x)
# [1] 0.03005828
###############################################

#do CVs as distributions
for(v in 1:nvars){
  thisVar<-testTracers[v]
  x<-as.vector(storeBetweenRunCV[,,,v])
  if(sum(x,na.rm=TRUE)>0){

    jpeg(paste(plotPath,"ROMS_CV_Density_",version,"_",thisVar,".jpg",sep=""))
    
    h<-hist(x,main="",freq=FALSE,ylab="Density",xlab="CV")
    
    xfit<-seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length=40) 
    yfit<-dnorm(xfit,mean=mean(x,na.rm=TRUE),sd=sd(x,na.rm=TRUE)) 
    yfit <- yfit*diff(h$mids[1:2])*length(x) 
    par(new=TRUE)
    plot(xfit, yfit, col=myBlue, lwd=2,type="l",xaxt="n",yaxt="n",ylab="",xlab="")
    mtext(thisVar,side=3,adj=0)
    
    dev.off()
  }

}
##plot them by box
colByBox<-colorRampPalette(colors=c(myOrange,myRed,myPurple,myBlue,myAqua,myGreen))(nboxes)
boundaryBoxes<-c(1,seq(26,30))
for(v in 1:nvars){
  thisVar<-testTracers[v]
  jpeg(paste(plotPath,"ROMS_CV_DensityByBOX_",version,"_",thisVar,".jpg",sep=""))
  firstPlot=TRUE
  
  for(b in 1:nboxes){
    if(!b %in% boundaryBoxes){
      thisCol<-colByBox[b]
      x<-as.vector(storeBetweenRunCV[,b,,v])
      if(sum(x,na.rm=TRUE)>0){
        
        h<-hist(x,plot=FALSE)
        
        xfit<-seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length=40) 
        yfit<-dnorm(xfit,mean=mean(x,na.rm=TRUE),sd=sd(x,na.rm=TRUE)) 
        yfit <- yfit*diff(h$mids[1:2])*length(x) 
        if(firstPlot){
          plot(xfit, yfit, col=myBlue, lwd=2,type="l",xaxt="n",yaxt="n",ylab="",xlab="")
          mtext(thisVar,side=3,adj=0)
          firstPlot<-FALSE
        }else{
          par(new=TRUE)
          plot(xfit, yfit, col=thisCol, lwd=2,type="l",xaxt="n",yaxt="n",ylab="",xlab="")
        }
      }
    }
  }
  dev.off()
}

##plot with a map for those that look interesting
#read in shape file
shapeFile<-paste(DIR$'Base',"ATLANTISmodels\\inputs\\bgm\\CHAT30_LL",sep="")
sdata<-read.shapefile(shapeFile)
shape<-formatShape(shapeFile=shapeFile)
ns<-length(shape)
SpDF <- SpatialPolygonsDataFrame(shape,data.frame( z=1:ns, row.names=paste("P",seq(1,(ns)),sep="")))
labels<-seq(1,(ns))
plot(shape)
LABELpOS<-polygonsLabel(shape, labels = labels, cex=.1,doPlot=FALSE)
labeldf<-data.frame(cbind("x"=LABELpOS[1:ns],"y"=LABELpOS[(ns+1):(2*ns)]))

mapTracers<-c("Diatom_N","Pelagic_fish_lge1_Nums")
mapTracers<-c("Meiobenth_N")
for(v in 1:nvars){
  thisVar<-testTracers[v]
  if(thisVar %in% mapTracers){
    pdf(paste(plotPath,"ROMS_CV_DensityByBOXWithMAP_",version,"_",thisVar,".pdf",sep=""))
    par(mfrow=c(2,2))
    for(b in 1:nboxes){
      if(!b %in% boundaryBoxes){
        #first plot all boxes, so have them in the background
        #do the background one with all boxes to get full distribution range
        x<-as.vector(storeBetweenRunCV[,,,v])
        h<-hist(x,plot=FALSE)
        xfit<-seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length=40) 
        yfit<-dnorm(xfit,mean=mean(x,na.rm=TRUE),sd=sd(x,na.rm=TRUE)) 
        yfit <- yfit*diff(h$mids[1:2])*length(x) 
        plot(xfit, yfit, col=myGrey, lwd=2,type="n",ylab="",xlab="")
        
        for(ball in 1:nboxes){
          if(!ball %in% boundaryBoxes){
            x<-as.vector(storeBetweenRunCV[,ball,,v])
            if(sum(x,na.rm=TRUE)>0){
              h<-hist(x,plot=FALSE)
              thisCol<-colByBox[ball]
              
              xfit<-seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length=40) 
              yfit<-dnorm(xfit,mean=mean(x,na.rm=TRUE),sd=sd(x,na.rm=TRUE)) 
              yfit <- yfit*diff(h$mids[1:2])*length(x) 
              par(new=TRUE)
              plot(xfit, yfit, col=thisCol, lwd=2,type="l",xaxt="n",yaxt="n",ylab="",xlab="")
            }
          }
        }
        ##
        thisCol<-"black"
        x<-as.vector(storeBetweenRunCV[,b,,v])
        if(sum(x,na.rm=TRUE)>0){
            h<-hist(x,plot=FALSE)
            
            xfit<-seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length=40) 
            yfit<-dnorm(xfit,mean=mean(x,na.rm=TRUE),sd=sd(x,na.rm=TRUE)) 
            yfit <- yfit*diff(h$mids[1:2])*length(x) 
            par(new=TRUE)
            plot(xfit, yfit, col=thisCol, lwd=3,type="l",ylab="",xlab="",xaxt="n",yaxt="n")
            mtext(thisVar,side=3,adj=0)
            
            #add the map
            plot(shape)
            map('nzHires',add=TRUE,col="black",lwd=2)
            map.axes()
            for(plotB in 1:dim(labeldf)[1]){
              if(b== plotB){thisCol=myGreen}else{thisCol="white"}
              polygon(sdata$shp$shp[[plotB]]$points,col=thisCol)
            }
            mtext(paste("Box ",b,sep=""),side=3,adj=1)
            # mtext(paste("Depth= ",thisDepth," m",sep=""),side=3,adj=1,outer=FALSE)
        }
      }
    }
    dev.off()
  }

  
}


#for each variable, check out the variance between runs relative to the variance within runs. For each given box and layer, then summarise over the model region


