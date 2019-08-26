# read in catch nc file, and plot spatially

this_out<-"SNA3"; runFolder="TBGB_JP2"
thisDesc <- paste(runFolder, this_out,sep="")

this_path = paste(DIR$'Base', "TBGB\\",runFolder,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")

plotPath<-paste(DIR$'Base',"TBGB\\Figures\\Testing\\FIshing\\",sep="")

plotYearly<-TRUE

biolLines <- readLines(paste(this_path,"TBGB_biol.prm", sep=""))

#read in catches, as fomated and such in summaryCatchByFleet_TBGB.R
load(paste(this_path,"Catch_history\\storeCatchByFleet", sep=""))
nfg <- dim(fishedGroups)[1]

ThisNC.nc<-nc_open(paste(outPath,thisRun,".nc",sep=""))
thisVol <- ncvar_get(ThisNC.nc, "volume"); nts <- dim(thisVol)[3]
daysTimeStep<-73
numStepsPerYear<-365/daysTimeStep
year0<-1850
fishingStartYear<-1900
modelStartYear<-1950
newVol <- ncvar_get(newNC.nc, "volume")

xLabsTemp<-seq(burnin,(nts*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0
modelTime <- year0 + (1:nts)/numStepsPerYear

burnin <- modelStartYear-year0

lastYear <- max(modelTime)
plotYearlyIndex <- seq(burnin*numStepsPerYear, nts, by=numStepsPerYear)
# plotYearlyIndex <- seq(burnin*numStepsPerYear,nts)

getColor<-function(x,thisMax){
  thisCol<-"white"
  if(!is.na(x) & x>0){
    if(x<1){
      y<-1
    } else{
      y<-round(log10(x)/log10(thisMax),2)*100+1
    }
    thisCol<-thisColRamp[y]
  }
  return(thisCol)
}
thisColRamp <- colorRampPalette(colors=c(myLightBlue, myBlue, "midnightblue"))(101)

plotGrid<-function(x,y){
  thisX<-c(x-0.5,x-0.5,x+0.5,x+0.5); thisY<-c(y-0.5,y+0.5,y+0.5,y-0.5)
  thisCol<-plotColour[x,y]
  if(length(thisCol)>0){
    polygon(x=thisX,y=thisY,col=thisCol,border=NA)
  }
  return(NULL)
}


# check out spatial catches
catchNC.nc <- nc_open(paste(outPath, "outputCATCH.nc", sep=""))
catchVars <- sort(unique(names(catchNC.nc$var)))

thisCode <- "SCA"; thisName <- str_trim(groupsDF$Name[groupsDF$Code==thisCode])
thisVars <- catchVars[grep(paste(thisCode,"_Catch", sep=""), catchVars)]

thisCVar <- thisVars[1]
thisCdata <- ncvar_get(catchNC.nc, thisCVar)
for(i in 2:length(thisVars)){
  thisCVar <- thisVars[i]
  thisCdata <- thisCdata + ncvar_get(catchNC.nc, thisCVar)
}
ncatchYears <- dim(thisCdata)[2]; catchYears <- year0:(year0+ncatchYears-1)

biomassData <- ncvar_get(ThisNC.nc, paste(thisName,"_N", sep=""))
thisBiomassByts <- apply(biomassData * thisVol, 3, sum) * mg_2_tonne * X_CN
thisBiomassByYear<-thisBiomassByts[plotYearlyIndex]
# create as heatmap
plot(x=catchYears, y=catchByYear, type="h", lwd=5, col=myGrey, xlim=c(modelStartYear, max(catchYears)))
par(new=TRUE)
plot(x=modelTime[plotYearlyIndex], y=thisBiomassByts[plotYearlyIndex], type="l")  


# plotData<-t(tracersArray[t,,]); 
catchByYear <- apply(thisCdata, 2, sum); nonzeroIndex <- catchByYear>0
par(mfcol=c(2,2), mar=c(4,4,1,1))
xx<-thisBiomassByYear[nonzeroIndex]
xxx <- xx[seq(1,length(xx), by=numStepsPerYear)]
plotTime <- modelTime[seq(1,nts*numStepsPerYear, by=numStepsPerYear)][nonzeroIndex]
plotTime <- plotTime[!is.na(plotTime)]
prettyTime <- pretty(plotTime); prettyTimeAt <- prettyTime- prettyTime[1]+1

  plotData <- t(thisCdata[,nonzeroIndex]); thisMax <- max(plotData, na.rm=TRUE)
  plotColour<-apply(plotData,c(1,2),getColor,thisMax)
  tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotData)[1]),dim(plotData)[2]),"y"=sort(rep(seq(1,dim(plotData)[2]),dim(plotData)[1]))))
  par(mar=c(6,4,1,1))
  plot(x=seq(1,dim(plotData)[1]),y=rep(dim(plotData)[2],dim(plotData)[1]),type="n",xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,dim(plotData)[2]))
  mtext(allPlotTracers[t],side=3, adj=0)
  axis(at=prettyTimeAt,labels = prettyTime,side=1,las=2)
  axis(at=seq(1,dim(plotData)[2]),labels=colnames(plotData),side=2,las=1)
  temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
  box()
  mtext("Box",side=2,adj=0.5,line=3)
   plot(x=plotTime, y=xxx, type="l", ylab="Catch (tonnes)", xlab="Year", ylim=c(0, max(xxx)), col="midnightblue", lwd=2); 

  
  for(thisTime in timeSteps){
    plot(shape)
    map('nzHires',add=TRUE,col="black",lwd=2)
    map.axes()
    for(plotB in 1:dim(labeldf)[1]){
      Blab <- plotB +1
      if(Blab==(dim(labeldf)[1]+1)){Blab<-1}
      thisCol<-plotColour[thisTime, Blab]
      test<-sdata$shp$shp[[plotB]]$points
      thisX<-mean(test$'X'); thisY<-mean(test$'Y')
      text(Blab, x=thisX, y=thisY)
      polygon(sdata$shp$shp[[plotB]]$points,col=thisCol)
    }
    mtext(paste(allPlotTracers[t],", timestep ",thisTime, sep=""), side=3, adj=0)
  }
}




