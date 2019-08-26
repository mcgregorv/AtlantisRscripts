## read in initial conditions and plot spatial distributions. The read in tracers and plot spatial distribution to compare
# a bit irrelevant for 'verts' as they are not doing free movement, but comparison important for zoo's and below
library(ggmap)
library(rgdal)
library(gtable)
library(maps)
library(mapdata)
library(rgeos)
source(paste(DIR$'General functions',"\\get_first_number.R",sep=""))
source(paste(DIR$'General functions',"formatShape.R",sep=""))
this_run<-"base"
this_out<-paste("BASE50yr9",sep="")

basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
plotPath<-paste(DIR$'Figures',"Spatial\\",this_out, sep="")
## read in groups file and biol.prm lines
groupsDF<-read.csv(paste(basePath, "CRAM_groups.csv", sep="")); ng<-dim(groupsDF)[1]

# biol prm lines
thisBiolFile<-paste(basePath,"\\CRAM_BH_hybrid_biol.prm",sep="")
biolLines<-readLines(thisBiolFile)

this_file<-paste(basePath,"inputs\\biol_initial\\CRAM_Initial_biomass_distributions.xlsx",sep="")
wb <- loadWorkbook(this_file)
sheetNames<-names(wb)

test_nts<-50

mg_2_tonne<-2e-8; X_CN<-5.7

ThisNC.nc<-nc_open(paste(basePath,"base\\output",this_out,"\\output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc, "volume"); nboxes<-dim(thisVol)[2]; nts<-dim(thisVol)[3]; nlayers<-dim(thisVol)[1]
dynBoxIndex<-seq(1,nboxes) %in% seq(2,25)
volByBox<-apply(thisVol[,,1],2,sum); areaByBox<-thisVol[nlayers,,1]
storeInitialSpatial<-array(NA, dim=c(ng, nboxes))
storeInitialByArea<-array(NA, dim=c(ng,nboxes))
skip<-c("PL", "PS")
# skip<-c()
for(g in 1:ng){
  thisCode<-as.character(groupsDF$Code[g]); thisNumCohorts<-groupsDF$NumCohorts[g]
  if(!(thisCode %in% skip)){
    thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode], side="both");  thisTracer<-paste(thisName,"_N", sep=""); 
    ##use xlsx
    if(thisNumCohorts>1 | thisCode %in% c("PL", "PS")){
      sheetIndex<-grep(thisCode,sheetNames)
    } else{
      sheetIndex<-seq(1,length(sheetNames))[sheetNames==thisCode]
    }
    df2 <- read.xlsx(xlsxFile = this_file, sheet = sheetIndex, skipEmptyRows = TRUE)
    if(thisNumCohorts>1){
      icTonnesByBox<-as.double(df2[5:28,8])
      icTonnesPm3ByBox<-icTonnesByBox/volByBox[dynBoxIndex]
      icTonnesPm2ByBox<-icTonnesByBox/areaByBox[dynBoxIndex]
      storeInitialByArea[g,dynBoxIndex]<-icTonnesPm2ByBox
    } else {
      thisDim<-as.double(df2[1,2])
      thisNByBox<-df2[4:27,2]
      icTonnesPm3ByBox<-as.double(thisNByBox)*mg_2_tonne*X_CN
    }
    storeInitialSpatial[g,dynBoxIndex]<-icTonnesPm3ByBox
  }
}


totalByBox<-apply(storeInitialSpatial, 2, sum, na.rm=TRUE)

pdf(paste(plotPath,"InitialDens.pdf", sep=""))
par(lend=1, mfrow=c(5,2), mar=c(4,4.5,1.5,1.5))
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]
  if(sum(storeInitialSpatial[g,], na.rm=TRUE)>0){
    plot(storeInitialSpatial[g,], type="h", lwd=5, col=myBlue,xlab="Box",ylab=expression("Tonnes per "~m^3))
    mtext(thisCode,side=3,adj=0)
  }
}
dev.off()


pdf(paste(plotPath,"InitialDensByArea.pdf", sep=""))
par(lend=1, mfrow=c(5,2), mar=c(4,4.5,1.5,1.5))
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    if(sum(storeInitialByArea[g,], na.rm=TRUE)>0){
      plot(storeInitialByArea[g,], type="h", lwd=5, col=myBlue,xlab="Box",ylab=expression("Tonnes per "~m^2))
      mtext(thisCode,side=3,adj=0)
    }
  }
}
dev.off()

thisColRamp<-colorRampPalette(colors=c(myLightAqua, myBlue,"midnightblue"))(101)
getCol<-function(x){
  thisCol<-"white"
  y=round((x/thisMax),2)*100+1
  if(sum(x,na.rm=TRUE)>0){
    thisCol<-thisColRamp[y]
  }
  return(thisCol)
}

#read in shape file
# pdf(paste(plotPath,"test.pdf", sep=""), height=10)
shapeFile<-paste(DIR$'Base',"ATLANTISmodels\\inputs\\bgm\\CHAT30_LL",sep="")
sdata<-read.shapefile(shapeFile)
shape<-formatShape(shapeFile=shapeFile)
ns<-length(shape)
SpDF <- SpatialPolygonsDataFrame(shape,data.frame( z=1:ns, row.names=paste("P",seq(1,(ns)),sep="")))
labels<-seq(1,(ns))
plot(shape)
LABELpOS<-polygonsLabel(shape, labels = labels, cex=.1,doPlot=FALSE)
labeldf<-data.frame(cbind("x"=LABELpOS[1:ns],"y"=LABELpOS[(ns+1):(2*ns)]))
# dev.off()

pdf(paste(plotPath,"InitialDensByAreaMaps.pdf", sep=""), height=10)
par(lend=1, mfrow=c(8,2), mar=c(2,0.5,1.5,1.5))
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    if(sum(storeInitialByArea[g,], na.rm=TRUE)>0){
      thisData<-storeInitialByArea[g,]; thisMax<-max(thisData, na.rm=TRUE)
      thisColors<-unlist(lapply(thisData, getCol))
      plot(shape)
      # map('nzHires',add=TRUE,col="black",lwd=2)
      # map.axes()
      for(b in 1:dim(labeldf)[1]){
        polygon(sdata$shp$shp[[b]]$points,col=thisColors[b])
      }
      mtext(thisCode,side=3,adj=0)
    }
  }
}
dev.off()


pdf(paste(plotPath,"InitialDensMaps.pdf", sep=""), height=10)
par(lend=1, mfrow=c(8,2), mar=c(2,0.5,1.5,1.5))
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    if(sum(storeInitialSpatial[g,], na.rm=TRUE)>0){
      thisData<-storeInitialSpatial[g,]; thisMax<-max(thisData, na.rm=TRUE)
      thisColors<-unlist(lapply(thisData, getCol))
      plot(shape)
      # map('nzHires',add=TRUE,col="black",lwd=2)
      # map.axes()
      for(b in 1:dim(labeldf)[1]){
        polygon(sdata$shp$shp[[b]]$points,col=thisColors[b])
      }
      mtext(thisCode,side=3,adj=0)
    }
  }
}
dev.off()

#add up tonnes per m^2 from tracers at a timestep
this_nts<-min(nts-1,test_nts)
storeTracerSpatial<-array(NA, dim=c(ng,nboxes))
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]; thisName<-str_trim(groupsDF$Name[g], side="both")
  thisTracer<-paste(thisName,"_N", sep=""); thisData<-ncvar_get(ThisNC.nc, thisTracer)
  if(length(dim(thisData))==3){
    thisTonnes<-apply(thisData*thisVol, c(2,3), sum, na.rm=TRUE)[,this_nts]
  } else{
    thisTonnes<-apply(thisData*thisVol[6,,],c(1,2),sum,na.rm=TRUE)[,this_nts]
  }
  thisTonnesPm2ByBox<-thisTonnes/areaByBox
  storeTracerSpatial[g,]<-thisTonnesPm2ByBox
}



pdf(paste(plotPath,"DensMapsTS_",this_nts,".pdf", sep=""), height=10)
par(lend=1, mfrow=c(8,2), mar=c(2,0.5,1.5,1.5))
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
  # if(thisNumCohorts>1){
    if(sum(storeTracerSpatial[g,], na.rm=TRUE)>0){
      thisData<-storeTracerSpatial[g,]
      thisData[!dynBoxIndex]<-0; thisMax<-max(thisData, na.rm=TRUE)
      thisColors<-unlist(lapply(thisData, getCol))
      plot(shape)
      # map('nzHires',add=TRUE,col="black",lwd=2)
      # map.axes()
      for(b in 1:dim(labeldf)[1]){
        polygon(sdata$shp$shp[[b]]$points,col=thisColors[b])
      }
      mtext(thisCode,side=3,adj=0)
    }
  # }
}
dev.off()


pdf(paste(plotPath,"DensMapsTSnoAge_",this_nts,".pdf", sep=""), height=11)
par(lend=1, mfrow=c(9,2), mar=c(2,0.5,1.5,1.5))
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts==1){
  if(sum(storeTracerSpatial[g,], na.rm=TRUE)>0){
    thisData<-storeTracerSpatial[g,]
    thisData[!dynBoxIndex]<-0; thisMax<-max(thisData, na.rm=TRUE)
    thisColors<-unlist(lapply(thisData, getCol))
    plot(shape)
    # map('nzHires',add=TRUE,col="black",lwd=2)
    # map.axes()
    for(b in 1:dim(labeldf)[1]){
      polygon(sdata$shp$shp[[b]]$points,col=thisColors[b])
    }
    mtext(thisCode,side=3,adj=0)
  }
  }
}
dev.off()

bpByBox<-apply(storeTracerSpatial[groupsDF$NumCohorts==1, ],2, sum, na.rm=TRUE)

thisMax<-max(bpByBox, na.rm=TRUE)
plotColors<-unlist(lapply(bpByBox, getCol))
plot(shape)
# map('nzHires',add=TRUE,col="black",lwd=2)
# map.axes()
for(b in 1:dim(labeldf)[1]){
  polygon(sdata$shp$shp[[b]]$points,col=plotColors[b])
}


plot(seq(0,1), type="n",ylab="Observed", xlab="Predicted")
for(g in 1:ng){
  thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    par(new=TRUE)
    plot(x=storeInitialSpatial[g,], y=storeTracerSpatial[g,], pch=20, col=myBlue_trans,xlab="",ylab="",xaxt="n",yaxt="n")
  }
}

## spatial distribution by trophic level
#read in trophic levels
ppGroups<-as.character(read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ppGroups.csv",sep=""))[,1]);
npg<-length(ppGroups)

groupsTL<-read.csv(paste(basePath,"inputs\\biol_prm\\CRAM_trophicLevels_isotopes.csv",sep=""))
ppTLs<-rep(NA,npg)
for(g in 1:npg){
  thisVar<-ppGroups[g]
  thisCode<-gsub("juv|ad","",thisVar)
  thisAge<-gsub(thisCode,"",thisVar)
  thisTLrow<-groupsTL[groupsTL$Code==thisCode,]
  TL_A<-thisTLrow$Isotope
  if(thisAge=="juv"){
    thisTL<-max(thisTLrow$TrophicLevel1,TL_A-0.5,na.rm=TRUE)
  } else { #either adult or un-aged here
    thisTL<-TL_A
    if(is.na(TL_A)){thisTL<-thisTLrow$TrophicLevel2}
  }
  if(is.na(thisTL)){thisTL<-1}
  ppTLs[g]<-thisTL
}

getCol<-function(x){
  thisCol<-"white"
  y=round(((x)/(thisMax)),2)*100+1
  if(sum(x,na.rm=TRUE)>0){
    thisCol<-thisColRamp[y]
  }
  # return(y)
  return(thisCol)
}
gTLs<-rep(NA, ng)
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]
  thisTLrow<-groupsTL[groupsTL$Code==thisCode,]
  TL_A<-thisTLrow$Isotope
  if(thisAge=="juv"){
    thisTL<-max(thisTLrow$TrophicLevel1,TL_A-0.5,na.rm=TRUE)
  } else { #either adult or un-aged here
    thisTL<-TL_A
    if(is.na(TL_A)){thisTL<-thisTLrow$TrophicLevel2}
  }
  if(is.na(thisTL)){thisTL<-1}
  gTLs[g]<-thisTL
}

tlCols<-c(myYellow,myGreen,myBlue,myPurple,myRed,myGrey)

#set up biomass by trophic level
pdf(paste(plotPath,"TLspatial.pdf",sep=""))
par(mfrow=c(3,2))
for(t in 0:5){
  thisTLindex<-round(gTLs,0)==t
  thisBiomassTracers<-storeTracerSpatial[thisTLindex,]
  thisBiomassByBox<-apply(thisBiomassTracers,2,sum,na.rm=TRUE)
  thisCol<-tlCols[t+1]
  thisMax<-max(thisBiomassByBox,na.rm=TRUE); thisColRamp<-colorRampPalette(colors=c("white",thisCol,"black"))(101)

  thisColors<-unlist(lapply(thisBiomassByBox,getCol))
  plot(shape)
  # map('nzHires',add=TRUE,col="black",lwd=2)
  # map.axes()
  for(b in 1:dim(labeldf)[1]){
    polygon(sdata$shp$shp[[b]]$points,col=thisColors[b])
  }
  mtext(paste("Trophic level ", t),side=3,adj=0)
}
dev.off()



##################
## plot IVS with catch distribution
IVSdata<-read.csv(paste(DIR$'Data',"fisheries\\scampi_target2017.csv",sep=""))
plotData<-IVSdata[,c("start_latitude", "start_longitude")]; 
shortData<-distinct(plotData[plotData$start_latitude>-44.7 & plotData$start_latitude<(-42.6) & plotData$start_longitude>174, ])
thisCode<-"IVS"; codeIndex<-groupsDF$Code==thisCode

pdf(paste(plotPath,"SpatialCompare_",thisCode,".pdf", sep=""), height=4, width=5)
par(mar=c(3,3,1,1))
thisData<-storeTracerSpatial[codeIndex,]
thisData[!dynBoxIndex]<-0; thisMax<-max(thisData, na.rm=TRUE)
thisColRamp<-colorRampPalette(colors=c(myVeryLightGrey,myLightGrey,myGrey,"black"))(101)
thisColors<-unlist(lapply(thisData, getCol))
plot(shape)
map('nzHires',add=TRUE,col="black",lwd=2)
map.axes()
for(b in 1:dim(labeldf)[1]){
  polygon(sdata$shp$shp[[b]]$points,col=thisColors[b])
}
# mtext(thisCode,side=3,adj=0)
points(x=shortData$start_longitude, y=shortData$start_latitude, pch=20, col=myOrange_trans)
dev.off()


splitString<-function(x,split=" "){
  y<-as.double(unlist(str_split(x,split)))
  return(y)
}
trawlBiomass<-c(0,splitString("1000	2500	2060	1522	4610	3333	1493	193	1149	5642	3208	2983	2434	516	180	132	50	58	50	50	50	50	70	70
", split="\t"), rep(0,5))
IVStrawlSurvey<-trawlBiomass/areaByBox
thisMax<-max(IVStrawlSurvey); thisColRamp<-colorRampPalette(colors = c(myLightGreen,myGreen,myDarkGreen))(101)
tsColors<-unlist(lapply(IVStrawlSurvey, getCol))
pdf(paste(plotPath,"SpatialTrawlSurvey_",thisCode,".pdf", sep=""), height=4, width=5)
par(mar=c(3,3,1,1))
plot(shape)
map('nzHires',add=TRUE,col="black",lwd=2)
map.axes()
for(b in 1:dim(labeldf)[1]){
  polygon(sdata$shp$shp[[b]]$points,col=tsColors[b])
}
dev.off()