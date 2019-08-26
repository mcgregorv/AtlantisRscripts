this_run<-"Chaos"
this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
thisDesc<-"DChaosUpDownBASE"
thisDesc <-"ChaosSampleA"

thisRunTemplate <- "ChaosFISH"
thisRunTemplate <- "Chaos|short"
thisRunTemplate <- "DChaos"; baseRunIndex<-5
thisRunTemplate <- "ChaosSampleA"; baseRunIndex<-1



nlayers<-6

groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
plotPath<-paste(this_path,"..\\Figures\\", thisDesc,sep="")

# brings in "relTracersArray", "tracers2plot", "tracersArray", "this_out", "runCols", "timeList", "nts_list"
load(paste(this_path,thisDesc,"modelTracers",sep=""));
max_nts <- max(nts_list)
nruns<-length(this_out)
runCols<-colorRampPalette(colors=c("midnightblue",myBlue,myAqua,myGold,  myOrange, "red"))(nruns)
base_nts <- nts_list[[baseRunIndex]]

thisMin<-min(relTracersArray, na.rm=TRUE); thisMax <- max(relTracersArray, na.rm=TRUE)

timeGroups <- seq(50,max_nts, by=1); ntgs<-length(timeGroups)
plot()
for(t in 2:ntgs){
  thisTimeMax<-timeGroups[t]; thisTimeMin <- max(c(1, timeGroups[t-1]), na.rm=TRUE)
  thisData <- relTracersArray[,,thisTimeMin:(thisTimeMax-1)]
  hist(thisData, main=paste("Years ", thisTimeMin+1865,":",thisTimeMax-1+1865, sep=""))
  # thisHist<-hist(thisData, plot=FALSE)
  # thisY <- lowess(thisHist$counts)$y
  # thisY <- thisHist$counts
  # points(x=thisHist$breaks[-1]-0.5, y=thisY, type="l", col=colByTime[t])
}


#read in biol.prm file so can subset by lifespan 
biolLines <- readLines(paste(this_path,"..\\CRAM_BH_hybrid_biol.prm", sep=""))
getLifeSpan<-function(tracer){
  thisName<-gsub("_N", "", tracer)
  thisCode <- groupsDF$Code[grep(thisName, groupsDF$Name)]
  if(length(thisCode)>1){thisCode <- groupsDF$Code[groupsDF$Name==thisName]}
  thisNumCohorts <- groupsDF$NumCohorts[groupsDF$Code==thisCode]
  thisLifespan <- NA
  if(length(thisCode)>0){
    if(thisNumCohorts>1){
      thisVar <- paste(thisCode, "_AgeClassSize", sep="")
      thisACS<-get_first_number(biolLines[grep(thisVar, biolLines)])
      thisLifespan <- thisNumCohorts * thisACS
    }
  }
  return(thisLifespan)
}  
tracerLifespans <- unlist(lapply(tracers2plot, getLifeSpan))
tracerLifespans[is.na(tracerLifespans)]<-0
sort(unique(tracerLifespans))

lifeSpans <- sort(unique(tracerLifespans)); nlifes <- length(lifeSpans)
## want them by lifespan - as timeseries..?
## give a heat map a whirl
# set up df - nlifespans by nts -just base it on means for now
## now to plot as heat map
getColor<-function(x,thisMax){
  thisCol<-"white"
  if(!is.na(x) & x>0){
    y<-round((x-thisMin)/(thisMax-thisMin),2)*100+1
    thisCol<-thisColRamp[y]
  }
  return(thisCol)
}
getColorLOG<-function(x, thisMax=100){
  thisCol<-"white"
  if(!is.na(x)){
    z<-(log(abs(x), base=10) -log(0.001, base=10))  / (log(100, base=10) - log(0.001, base=10))
    if(z<0){z<-0}
    y<-100 * z + 1
    thisCol<-thisColRamp[y]
      if(x==0){thisCol="white"}
    if(x>thisMax){thisCol=thisColRamp[101]}
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
thisColRamp<-colorRampPalette(colors=c(myLightAqua,myAqua,"midnightblue"))(101)

thisRunIndex<-c(3,7)

for( r in 1:nruns){
  if(r !=baseRunIndex){

    meanByLifeSpan <- data.frame(matrix(NA, ncol=(nlifes-1), nrow=base_nts));
    for(l in 2:nlifes){
      thisLifespan <- lifeSpans[l]
      thisLifeIndex <- tracerLifespans==thisLifespan
      thisData<-relTracersArray[r,thisLifeIndex,1:base_nts]
      if(length(dim(thisData))==2){
        thisMeanByTime <- apply(abs(thisData), 2, mean, na.rm=TRUE)
      } else{
        thisMeanByTime <- abs(thisData)
      }
      meanByLifeSpan[,l-1]<- thisMeanByTime
    }
    # 
    # t=grep("Pelagic_fish_sml_N", tracers2plot)
    # thisBaseData<-tracersArray[baseRunIndex, t,]; thisRunData<-tracersArray[4,t,]
    # plot(thisBaseData, type="l")
    # points(thisRunData, type="l", lty=2)
    # test<-100*abs(thisBaseData - thisRunData)/thisBaseData
    # plot(test, type="l")
    
    plotData<-meanByLifeSpan
    thisMax<-max(plotData, na.rm=TRUE); thisMin<-min(plotData, na.rm=TRUE)
    
    plotColour<-apply(plotData,c(1,2),getColor,thisMax)
    tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotData)[1]),dim(plotData)[2]),"y"=sort(rep(seq(1,dim(plotData)[2]),dim(plotData)[1]))))
    # pdf(paste(plotPath,"KemptonsQ_SSRsens.pdf",sep=""),height=4,width=10)
    par(mar=c(6,4,1.5,1))
    plot(x=seq(1,dim(plotData)[1]),y=rep(dim(plotData)[2],dim(plotData)[1]),type="n",xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,dim(plotData)[2]))
    axis(at=1:base_nts,labels = timeList[[baseRunIndex]],side=1,las=2)
    axis(at=1:(nlifes-1),labels=lifeSpans[-1],side=2,las=1)
    temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
    box()
    mtext(signif(thisMax,2), side=3, adj=1, col=runCols[r])
  }
}

## swap perspective; loop through lifespans, and plot for each run
pdf(paste(plotPath, "meanByRunLifespan.pdf", sep=""), height=7, width=7)
par(mfrow=c(4,2))
for(l in 2:nlifes){
  thisLifespan <- lifeSpans[l]
  thisLifeIndex <- tracerLifespans==thisLifespan
    
    meanByRun <- data.frame(matrix(NA, ncol=(nruns), nrow=base_nts));
    for(r in 1:nruns){
      if(r !=baseRunIndex){
        

      thisData<-relTracersArray[r,thisLifeIndex,1:base_nts]
      if(length(dim(thisData))==2){
        thisMeanByTime <- apply(abs(thisData), 2, mean, na.rm=TRUE)
      } else{
        thisMeanByTime <- abs(thisData)
      }
      meanByRun[,r]<- thisMeanByTime
      }
    }
    plotData<-meanByRun
    thisMax<-max(plotData, na.rm=TRUE); thisMin<-min(plotData, na.rm=TRUE)
    
    plotColour<-apply(plotData,c(1,2),getColor,thisMax)
    tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotData)[1]),dim(plotData)[2]),"y"=sort(rep(seq(1,dim(plotData)[2]),dim(plotData)[1]))))
    # pdf(paste(plotPath,"KemptonsQ_SSRsens.pdf",sep=""),height=4,width=10)
    par(mar=c(6,4,1.5,1))
    plot(x=seq(1,dim(plotData)[1]),y=rep(dim(plotData)[2],dim(plotData)[1]),type="n",xlab="",ylab="Run",xaxt="n",yaxt="n",ylim=c(0,dim(plotData)[2]))
    axis(at=(1:base_nts)[seq(1,base_nts, by=20)],labels = (timeList[[baseRunIndex]][seq(1,base_nts, by=20)]-1),side=1,las=2)
    axis(at=1:nruns,labels=gsub(thisRunTemplate, "",  this_out),side=2,las=1)
    temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
    box()
    mtext(paste("Lifespan:", thisLifespan,"Max =", signif(thisMax,2), sep=" "), side=3, adj=1)
}
dev.off()

legendText<-c(1,5,10,20,50,100)
legendColor<-unlist(lapply(legendText, getColorLOG, thisMax=thisMax))
legendTextConverted <- unlist(lapply(legendText, FUN=function(x){10^x}))
makeBlankPlot()
legend(legend=legendText, col=legendColor, x="center", pch=15, bty="n")

## read in trophic levels
groupsTL<-read.csv(paste(this_path,"..\\inputs\\biol_prm\\CRAM_trophicLevels_isotopes.csv",sep=""))
getTL<-function(tracer){
  thisTL<-NA
  thisName<-gsub("_N", "", tracer)
  thisCode <- groupsDF$Code[grep(thisName, groupsDF$Name)]
  if(length(thisCode)>0){
    if(length(thisCode)>1){thisCode <- groupsDF$Code[groupsDF$Name==thisName]}
    thisTL<-groupsTL$Isotope[groupsTL$Code==thisCode]
    if(is.na(thisTL)){thisTL<-groupsTL$TrophicLevel2[groupsTL$Code==thisCode]}
  }
  return(thisTL)
}
tracerTLs <- unlist(lapply(tracers2plot, getTL))
tracerRoundTLs<-round(tracerTLs)
rTLs<-sort(unique(tracerRoundTLs)); nTLs <- length(rTLs)

pdf(paste(plotPath, "meanByRunTL.pdf", sep=""), height=10, width=5)
par(mfrow=c(6,1), mar=c(3,8,1,1))
for(l in 1:length(rTLs)){
  thisTL <- rTLs[l]
  thisTLIndex <- tracerRoundTLs==thisTL
  
  meanByRun <- data.frame(matrix(NA, ncol=(nruns), nrow=base_nts));
  for(r in 1:nruns){
    if(r !=baseRunIndex){
      
      
      thisData<-relTracersArray[r,thisTLIndex,1:base_nts]
      if(length(dim(thisData))==2){
        thisMeanByTime <- apply(abs(thisData), 2, mean, na.rm=TRUE)
      } else{
        thisMeanByTime <- abs(thisData)
      }
      meanByRun[,r]<- thisMeanByTime
    }
  }
  plotData<-meanByRun
  thisMax<-max(plotData, na.rm=TRUE); thisMin<-min(plotData, na.rm=TRUE)
  
  plotColour<-apply(plotData,c(1,2),getColor,thisMax)
  tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotData)[1]),dim(plotData)[2]),"y"=sort(rep(seq(1,dim(plotData)[2]),dim(plotData)[1]))))
  # pdf(paste(plotPath,"KemptonsQ_SSRsens.pdf",sep=""),height=4,width=10)
  par(mar=c(6,4,1.5,1))
  plot(x=seq(1,dim(plotData)[1]),y=rep(dim(plotData)[2],dim(plotData)[1]),type="n",xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,dim(plotData)[2]))
  axis(at=(1:base_nts)[seq(1,base_nts, by=20)],labels = (timeList[[baseRunIndex]][seq(1,base_nts, by=20)]-1),side=1,las=2)
  axis(at=1:nruns,labels=this_out,side=2,las=1)
  temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
  box()
  mtext(paste("Trophic level:", thisTL,"Max =", signif(thisMax,2), sep=" "), side=3, adj=1)
}
dev.off()


## how about the half life? time taken for difference from base to be half what it was at the start
# thisTimeSeries<-tracersArray[1,1,]; baseTimeSeries<-tracersArray[baseRunIndex,1,]
calcHalfLife<-function(thisTimeSeries, baseTimeSeries){
  xx<-abs(baseTimeSeries-thisTimeSeries)/baseTimeSeries
  thisStartDiff <- xx[1]
  x <- min((1:length(thisTimeSeries))[xx<(thisStartDiff/2)], na.rm=TRUE)
  ## check if it leaves this bound
  test<-xx>(thisStartDiff/2); testMin<-min((seq(1,length(thisTimeSeries))[test]), na.rm=TRUE)
  if(testMin>x){x<-NA}
  return(x)
}

storeHalfLifes <- array(NA, dim=dim(tracersArray)[1:2])
for(r in 1:nruns){
  if(r != baseRunIndex){
    for(t in 1:dim(tracersArray)[2]){
      thisTimeSeries <- tracersArray[r,t,]; baseTimeSeries<-tracersArray[baseRunIndex, t, ]
      thisHalfLife <- calcHalfLife(thisTimeSeries, baseTimeSeries)
      if(thisHalfLife=="Inf"){thisHalfLife<-NA}
      storeHalfLifes[r,t]<- thisHalfLife
    }
  }
}
## take out baseRun (all NAs)
storeHalfLifes <- storeHalfLifes[-baseRunIndex,]
toPlot<-melt(storeHalfLifes); colnames(toPlot)<-c("run", "tracerNum", "halflife")
toPlot$tracer <- tracers2plot[toPlot$tracerNum]
toPlot$lifespan <- unlist(lapply(toPlot$tracer, getLifeSpan))
toPlot$TL <- unlist(lapply(toPlot$tracer, getTL))
toPlot$run <- as.factor(toPlot$run)
toPlot$product <- toPlot$TL * toPlot$lifespan




#############################
TLcolors<-colorRampPalette(colors=c(myGold,myGreen,myAqua,myBlue,"red"))(11)
getTLcolor<-function(x){
  thisCol<-"white"
  y<-round(2*x)
  if(y %in% 1:11){
    thisCol<-TLcolors[y]
  }
  return(thisCol)
}
toPlot$TLcolor<-unlist(lapply(toPlot$TL, getTLcolor))

plot(x=toPlot$lifespan, y=toPlot$halflife, pch=20, col=toPlot$TLcolor)
plot(x=toPlot$TL, y=toPlot$halflife, pch=20)
plot(x=(toPlot$product)^(1/2), y=(toPlot$halflife), pch=20, ylim=c(0,70), xlim=c(0,70))
points(x=c(0,80), y=c(0,80), type="l", col="red")

plot(x=(toPlot$product), y=(toPlot$halflife), pch=20, ylim=c(0,70))

toPlot$lifespan[is.na(toPlot$lifespan)]<-0
toPlot$lifespan <- as.factor(toPlot$lifespan)
toPlot$TL<-as.factor(myRounding(toPlot$TL, 0.5))
bp<-ggplot(data = toPlot, aes(x = run, fill = TL, y = halflife)) + 
  geom_bar(stat = 'identity')
bp

bp<-ggplot(data = toPlot, aes(x = run, fill = lifespan, y = halflife)) + 
  geom_bar(stat = 'identity')
bp

thisIndex <- toPlot$run==1 
hist(toPlot$halflife[thisIndex])

bp<-ggplot(data = toPlot, aes(x = lifespan, fill = run, y = halflife)) + 
  geom_bar(stat = 'identity')
bp
bp<-ggplot(data = toPlot, aes(x = TL, fill = run, y = halflife)) + 
  geom_bar(stat = 'identity')
bp
bp<-ggplot(data = toPlot, aes(x = product, fill = run, y = halflife)) + 
  geom_bar(stat = 'identity')
bp

for(r in 1:(nruns-1)){
  thisIndex <- toPlot$run==r
  bp<-ggplot(data = toPlot[thisIndex,], aes(x = lifespan, fill = run, y = halflife)) + 
    geom_bar(stat = 'identity')
  bp
}

# pdf(paste(plotPath,"DietSummary.pdf", sep=""),height=7,width=8)
# par(mar=c(4,3.5,1,1))
# bp + coord_flip()  + scale_fill_manual(values=preyGroupColours) + labs(y="Proportion of diet", x="") + theme_igray() + 
#   theme(axis.text=element_text(size=12, angle=0),axis.title=element_text(size=12)) + guides(fill=guide_legend(title="")) 
# dev.off()







