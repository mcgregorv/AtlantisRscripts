## have 3 runs per group - read these in, quantify the change in focus groups biomass, then quantifiy the changes in the other groups relative to that in focus group
# read them all in first and store, then do plots
this_run<-"base"
basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
outPath<-paste(basePath, "base\\IntSims\\", sep="")
plotPath<-paste(basePath,"Figures\\intEffects\\",sep="")

dataOutPath<-outPath

mg_2_tonne<-2e-8; X_CN<-7.5

nts<-50; burnin<-1 ## these are all 50 year runs without 35 year burnin
nruns<-4 # each focus group has 3 runs with additional mL and 4th one is the base


## read in groups file and biol.prm lines
groupsDF<-read.csv(paste(basePath, "CRAM_groups.csv", sep="")); 
ng<-dim(groupsDF)[1]

biomassByFocusGroup<-array(NA, dim=c(ng, nruns, ng, nts))

calcAE<-function(P,O){
  a<-0; n<-min(length(P), length(O))
  for(i in 1:n){a<-a + (P[i]-O[i])}
  b<-a/n
  return(b)
}

BaseModel<-nc_open(paste(basePath,"base\\outputBase\\output.nc", sep="")); baseVol<-ncvar_get(BaseModel, "volume")

for(f in 1:ng){
  thisFCode<-as.character(groupsDF$Code[f]); thisFNumCohorts<-groupsDF$NumCohorts[f]
  if(thisFNumCohorts>1){
    cat(thisFCode,"--")
    for(r in 1:(nruns-1)){
      this_out<-paste("outputXXX_mLSens", thisFCode, r,sep=""); thisOutFile<-paste(outPath, this_out, "\\output.nc", sep="")
      if(file.exists(thisOutFile)){
        ThisNC.nc<-nc_open(thisOutFile)
        thisVol<-ncvar_get(ThisNC.nc, "volume")
        #now get all the biomass tracers for this run
        for(g in 1:ng){
          as.character(thisCode<-groupsDF$Code[g]); thisNumCohorts<-groupsDF$NumCohorts[g]
          if(thisNumCohorts>1){
            thisName<-str_trim(groupsDF$Name[g], side="both"); thisTracer<-paste(thisName,"_N", sep="")
            thisData<-ncvar_get(ThisNC.nc, thisTracer)
            this_tsdata<-apply(thisData * thisVol, 3, sum)* mg_2_tonne *X_CN
            biomassByFocusGroup[f,r,g,]<-this_tsdata[1:nts]
          }
        }
      }
    }
  }
}
#populate base model as well
r=4
for(f in 1:ng){
  thisFCode<-as.character(groupsDF$Code[f]); thisFNumCohorts<-groupsDF$NumCohorts[f]
  if(thisFNumCohorts>1){
    cat(thisFCode,"--")
    #now get all the biomass tracers for this run
    for(g in 1:ng){
      as.character(thisCode<-groupsDF$Code[g]); thisNumCohorts<-groupsDF$NumCohorts[g]
      if(thisNumCohorts>1){
        thisName<-str_trim(groupsDF$Name[g], side="both"); thisTracer<-paste(thisName,"_N", sep="")
        thisData<-ncvar_get(BaseModel, thisTracer)
        this_tsdata<-apply(thisData * baseVol, 3, sum)* mg_2_tonne *X_CN
        biomassByFocusGroup[f,r,g,]<-this_tsdata[1:nts]
      }
    }
  }
}

## now fit the bumped runs wrt base runs for each group rel to focus group
storeModels<-NULL; storeMs<-array(NA, dim=c(ng, ng, 4)); storeBs<-storeMs; storeR2s<-storeMs; storeRCorrS<-storeMs; storeRCorrP<-storeMs #one for up, down, and combined
for(f in 1:ng){
  thisFCode<-as.character(groupsDF$Code[f]); thisFNumCohorts<-groupsDF$NumCohorts[f]
  if(thisFNumCohorts>1){
    cat(thisFCode,"--")
    thisBaseData<-biomassByFocusGroup[f,4,,]
    alt1Data<-biomassByFocusGroup[f,1,,]; alt2Data<-biomassByFocusGroup[f,2,,]; alt3Data<-biomassByFocusGroup[f,3,,]
    if(sum(thisBaseData, na.rm=TRUE)>0 & sum(alt1Data, na.rm=TRUE)>0 & sum(alt2Data, na.rm=TRUE)>0){
      for(g in 1:ng){
        thisCode<-as.character(groupsDF$Code[g]); thisNumCohorts<-groupsDF$NumCohorts[g]
        if(sum(thisBaseData[g,], na.rm=TRUE)>0 & sum(alt1Data[g,], na.rm=TRUE)>0 & sum(alt2Data[g,], na.rm=TRUE)>0){
          fitData1<-data.frame(cbind("g"=(alt1Data[g,]-thisBaseData[g,])/thisBaseData[g,], "f"=(alt1Data[f,]-thisBaseData[f,])/thisBaseData[f,]))
          fitData2<-data.frame(cbind("g"=(alt2Data[g,]-thisBaseData[g,])/thisBaseData[g,], "f"=(alt2Data[f,]-thisBaseData[f,])/thisBaseData[f,]))
          fitData3<-data.frame(cbind("g"=(alt3Data[g,]-thisBaseData[g,])/thisBaseData[g,], "f"=(alt3Data[f,]-thisBaseData[f,])/thisBaseData[f,]))
          combFitData<-rbind(fitData1, fitData2, fitData3)
          model1<-lm(g ~ f, data=fitData1); model2<-lm(g ~ f, data=fitData2); model3<-lm(g ~ f, data=fitData3)
          thisRs_1<-summary(model1)$r.squared; thisRs_2<-summary(model2)$r.squared; thisRs_3<-summary(model3)$r.squared
          storeR2s[f,g,1]<-thisRs_1; storeR2s[f,g,2]<-thisRs_2; storeR2s[f,g,3]<-thisRs_3
          storeMs[f,g,1]<-summary(model1)$coefficients[1]; storeBs[f,g,1]<-as.double(summary(model1)$coefficients[2])
          storeMs[f,g,2]<-summary(model2)$coefficients[1]; storeBs[f,g,2]<-summary(model2)$coefficients[2]
          storeMs[f,g,3]<-summary(model3)$coefficients[1]; storeBs[f,g,3]<-as.double(summary(model3)$coefficients[2])
          modelIndex1<-paste(f,g,1, sep="-"); modelIndex2<-paste(f,g,2,sep="-"); modelIndex3<-paste(f,g,3,sep="-")
          storeModels[[modelIndex1]]<-model1; storeModels[[modelIndex2]]<-model2; storeModels[[modelIndex3]]<-model3
          ## combined
          modelc<-lm(g ~ f, data=combFitData); 
          thisRs<-summary(modelc)$r.squared; storeR2s[f,g,4]<-thisRs
          storeMs[f,g,4]<-summary(modelc)$coefficients[1]; storeBs[f,g,4]<-as.double(summary(modelc)$coefficients[2])
          modelIndexc<-paste(f,g,4, sep="-");
          storeModels[[modelIndex1]]<-modelc; 
          
          # calc pearsons spearman on f and g
          tempCor<-rcorr(x=fitData1$f, y=fitData1$g, type="spearman"); storeRCorrS[f,g,1]<-tempCor$r[2,1]
          tempCor<-rcorr(x=fitData2$f, y=fitData2$g, type="spearman"); storeRCorrS[f,g,2]<-tempCor$r[2,1]
          tempCor<-rcorr(x=fitData3$f, y=fitData3$g, type="spearman"); storeRCorrS[f,g,3]<-tempCor$r[2,1]
          tempCor<-rcorr(x=combFitData$f, y=combFitData$g, type="spearman"); storeRCorrS[f,g,4]<-tempCor$r[2,1]
          
          # calc pearsons correlation on f and g
          tempCor<-rcorr(x=fitData1$f, y=fitData1$g, type="pearson"); storeRCorrP[f,g,1]<-tempCor$r[2,1]
          tempCor<-rcorr(x=fitData2$f, y=fitData2$g, type="pearson"); storeRCorrP[f,g,2]<-tempCor$r[2,1]
          tempCor<-rcorr(x=fitData3$f, y=fitData3$g, type="pearson"); storeRCorrP[f,g,3]<-tempCor$r[2,1]
          tempCor<-rcorr(x=combFitData$f, y=combFitData$g, type="pearson"); storeRCorrP[f,g,4]<-tempCor$r[2,1]
          
        }
      }
    }
   
  }
}

storeRCorr<-storeRCorrP
plotIndex<-groupsDF$NumCohorts>1

rcorrMin<-min(storeRCorr,na.rm=TRUE); rcorrMax<-max(storeRCorr,na.rm=TRUE)

negColRamp<-colorRampPalette(colors=c("lightyellow", myOrange, "red"))(101)
posColRamp<-colorRampPalette(colors=c(myLightAqua, myBlue, "midnightblue"))(101)

getRCorrpColor<-function(x){
  thisCol<-"white"
  if(!is.na(x)){
    if(x<0){
      y<-round(x/rcorrMin, 2) *100 +1
      thisCol<-posColRamp[y]
    } else{
      y<-round(x/rcorrMax,2)*100 + 1
      thisCol<-negColRamp[y]
    }
    
    if(x==1){thisCol<-negColRamp[101]}
    if(x==(-1)){thisCol<-posColRamp[101]}
    if(x==0){thisCol="white"}
    if(x>rcorrMax){thisCol=negColRamp[101]}
    if(x<rcorrMin){thisCol=posColRamp[101]}
  }
  return(thisCol)
}

# rcorrMin<-(-1); rcorrMax<-1

RCorrCol1<-apply(storeRCorr[plotIndex,plotIndex,1],c(1,2), getRCorrpColor)
RCorrCol2<-apply(storeRCorr[plotIndex,plotIndex,2],c(1,2), getRCorrpColor)
RCorrCol3<-apply(storeRCorr[plotIndex,plotIndex,3],c(1,2), getRCorrpColor)
RCorrColCom<-apply(storeRCorr[plotIndex, plotIndex, 4], c(1,2), getRCorrpColor)

## plot using Map
plotGrid<-function(x,y){
  thisX<-c(x-0.5,x-0.5,x+0.5,x+0.5); thisY<-c(y-0.5,y+0.5,y+0.5,y-0.5)
  thisCol<-plotColour[x,y]
  if(length(thisCol)>0){
    polygon(x=thisX,y=thisY,col=thisCol,border=NA)
  }
  return(NULL)
}
plotColour<-RCorrCol1
plotColour<-RCorrCol2
plotColour<-RCorrColCom

listColors<-list(RCorrCol1, RCorrCol2, RCorrCol3, RCorrColCom)
for(p in 1:4){
  plotColour<-listColors[[p]]
  # pdf(paste(plotPath,"spearmansCorrelationRun",p,".pdf", sep=""), height=7, width=7)
  pdf(paste(plotPath,"PearsonsCorrelationRun",p,".pdf", sep=""), height=7, width=7)
  
par(mar=c(4,4,1,1))
  tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotColour)[1]),dim(plotColour)[2]),"y"=sort(rep(seq(1,dim(plotColour)[2]),dim(plotColour)[1]))))
  
  par(mar=c(6,4,1.5,1))
  plot(x=seq(1,dim(plotColour)[1]),y=rep(dim(plotColour)[2],dim(plotColour)[1]),type="n",xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,dim(plotColour)[2]))
  mtext("Focus group", side=1, adj=0.5,line=4)
  axis(at=seq(1,dim(plotColour)[1]),labels = as.character(groupsDF$Code[plotIndex]),side=1,las=2)
  axis(at=seq(1,dim(plotColour)[2]),labels=as.character(groupsDF$Code[plotIndex]),side=2,las=1)
  temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
  mtext("Response group", side=2, adj=0.5, line=3)
  box()
  dev.off()
}

##plot legend
pdf(paste(plotPath,"SpearmasCorrelationRunLEGEND.pdf", sep=""), height=4, width=2.5)
par(mar=c(0,0,0,0))
legend<-seq(-1,1,by=0.2); legendColors<-unlist(lapply(legend, getRCorrpColor))
makeBlankPlot()
legend(legend=legend, col=legendColors, x="center", bty="n", pch=15, pt.cex=1.5, title="Spearman's\n correlation")
dev.off()

###############
##group by trophic level
temp<-read.csv(paste(basePath,"inputs\\biol_prm\\CRAM_trophicLevels_isotopes.csv",sep=""))
groupsTL<-temp[plotIndex,]
groupsTL$Isotope[is.na(groupsTL$Isotope)]<-groupsTL$TrophicLevel2[is.na(groupsTL$Isotope)]
TLindex<-order(groupsTL$Isotope)

labelNames<-gsub("_","", groupsDF$LongName[plotIndex][TLindex])
labelNames<-gsub("_"," ", groupsDF$Name[plotIndex][TLindex])

RCorrColTL1<-apply(storeRCorr[plotIndex, plotIndex, 1][TLindex,TLindex], c(1,2), getRCorrpColor)
RCorrColTL2<-apply(storeRCorr[plotIndex, plotIndex, 2][TLindex,TLindex], c(1,2), getRCorrpColor)
RCorrColTL3<-apply(storeRCorr[plotIndex, plotIndex, 3][TLindex,TLindex], c(1,2), getRCorrpColor)
RCorrColTL<-apply(storeRCorr[plotIndex, plotIndex, 4][TLindex,TLindex], c(1,2), getRCorrpColor)
listColors<-list(RCorrColTL1, RCorrColTL2, RCorrColTL3, RCorrColTL)
for(p in 1:4){
  plotColour<-listColors[[p]]
  # pdf(paste(plotPath,"spearmansCorrelationRunByTL",p,".pdf", sep=""), height=8, width=8)
  pdf(paste(plotPath,"PearsonsCorrelationRunByTL",p,".pdf", sep=""), height=8, width=8)
  
  par(mar=c(4,4,1,1))
  tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotColour)[1]),dim(plotColour)[2]),"y"=sort(rep(seq(1,dim(plotColour)[2]),dim(plotColour)[1]))))
  
  par(mar=c(10,10,1.5,1))
  plot(x=seq(1,dim(plotColour)[1]),y=rep(dim(plotColour)[2],dim(plotColour)[1]),type="n",xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,dim(plotColour)[2]))
  mtext("Focus group", side=1, adj=0.5,line=9)
  axis(at=seq(1,dim(plotColour)[1]),labels = labelNames,side=1,las=2)
  axis(at=seq(1,dim(plotColour)[2]),labels=labelNames,side=2,las=1)
  temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
  mtext("Response group", side=2, adj=0.5, line=9)
  box()
  dev.off()
}

############### get lowest point for each focus group under each scenario
lowestPoint<-array(NA, dim=c(ng, 3)); highestPoint<-0*lowestPoint; endPoint<-0*highestPoint
allLowestPoints<-array(NA, dim=c(ng, ng, 3)); allHighestPoints<-0*allLowestPoints; allEndPoints<-0*allHighestPoints
for(f in 1:ng){
  thisFCode<-as.character(groupsDF$Code[f]); thisFNumCohorts<-groupsDF$NumCohorts[f]
  if(thisFNumCohorts>1){
    cat(thisFCode,"--")
    thisBaseData<-biomassByFocusGroup[f,4,,]
    alt1Data<-biomassByFocusGroup[f,1,,]; alt2Data<-biomassByFocusGroup[f,2,,]; alt3Data<-biomassByFocusGroup[f,3,,]
    if(sum(thisBaseData, na.rm=TRUE)>0 & sum(alt1Data, na.rm=TRUE)>0 & sum(alt2Data, na.rm=TRUE)>0){
          fitData1<-(alt1Data[f,]-thisBaseData[f,])/thisBaseData[f,]
          fitData2<-(alt2Data[f,]-thisBaseData[f,])/thisBaseData[f,]
          fitData3<-(alt3Data[f,]-thisBaseData[f,])/thisBaseData[f,]
          lowestPoint[f,]<-signif(c(min(fitData1), min(fitData2), min(fitData3)),2)*100
          highestPoint[f,]<-signif(c(max(fitData1), max(fitData2), max(fitData3)),2)*100
          endPoint[f,]<-signif(c(fitData1[nts], fitData2[nts], fitData3[nts]),2)*100
          for(g in 1:ng){
            thisCode<-as.character(groupsDF$Code[g]); thisNumCohorts<-groupsDF$NumCohorts[g]
            if(sum(thisBaseData[g,], na.rm=TRUE)>0 & sum(alt1Data[g,], na.rm=TRUE)>0 & sum(alt2Data[g,], na.rm=TRUE)>0){
              fitData1<-(alt1Data[g,]-thisBaseData[g,])/thisBaseData[g,]
              fitData2<-(alt2Data[g,]-thisBaseData[g,])/thisBaseData[g,]
              fitData3<-(alt3Data[g,]-thisBaseData[g,])/thisBaseData[g,]
              allLowestPoints[f,g,]<-signif(c(min(fitData1, na.rm=TRUE), min(fitData2, na.rm=TRUE), min(fitData3, na.rm=TRUE)),2)*100
              allHighestPoints[f,g,]<-signif(c(max(fitData1, na.rm=TRUE), max(fitData2, na.rm=TRUE), max(fitData3, na.rm=TRUE)),2)*100
              allEndPoints[f,g,]<-signif(c(fitData1[nts],fitData2[nts], fitData3[nts]),2)*100
            }
          }
    }
  }
}


plot(1,ylim=c(-100,100), xlim=c(-100,0), type="n", ylab="Response", xlab="Focus")
for(f in 1:ng){
  for(r in 1:3){
    fmin<-endPoint[f,r]
    if(!is.na(fmin)){
      if(fmin<0){
         datamin<-(allEndPoints[f,,r])
         points(x=rep(fmin, length(datamin)), y=datamin, pch=20, col=myGrey_trans)
      }
    }
  }
}

## plot separately for each run
par(mfrow=c(3,1), mar=c(4,4,1,1))
for(r in 1:3){
  plot(1,ylim=c(-100,100), xlim=c(-100,0), type="n", ylab="Response", xlab="Focus")
  mtext(r,side=3,adj=0)
  for(f in 1:ng){
      fmin<-endPoint[f,r]
      if(!is.na(fmin)){
        if(fmin<0){
          datamin<-(allEndPoints[f,,r])
          points(x=rep(fmin, length(datamin)), y=datamin, pch=20, col=myGrey_trans)
        }
      }
    }
}
### color by TL - use colors from foodweb plot
colByTLtrans<-c(myBlue_trans, myPurple_trans, myRed_trans, myGrey_trans)
colByTL<-c(myBlue, myPurple, myRed, myGrey)
getTLcolor<-function(x){
  thisCol<-"white"
  y<-trunc(x)-1; 
  if(y>0){thisCol<-colByTL[y]}
  return(thisCol)
}
TLcolByGroupCode<-unlist(lapply(groupsTL$Isotope, getTLcolor))
## just want age-structured end points
asIndex<-groupsDF$NumCohorts>1; 
ageStructuredCodes<-groupsDF$Code[asIndex]; nag<-length(ageStructuredCodes)
r=3
plot(1,ylim=c(-100,100), xlim=c(-100,0), type="n", ylab="Response", xlab="Focus")
mtext(r,side=3,adj=0)
for(f in 1:nag){
  thisCode<-ageStructuredCodes[f]; codeIndex<-groupsDF$Code==thisCode
  fmin<-endPoint[codeIndex,r]
  if(!is.na(fmin)){
    if(fmin<0){
      datamin<-(allEndPoints[codeIndex,,r])
      points(x=rep(fmin, length(datamin)), y=datamin, pch=20, col=myGrey_trans, cex=thisCex)
    }
  }
}
ageStructuredNames<-gsub("_", " ", groupsDF$Name[asIndex])
ageByCode<-(groupsDF$NumAgeClassSize * groupsDF$NumCohorts) [asIndex]
ageLengthIndex<-order(ageByCode)
par(mar=c(10,4,1,1))
plot(endPoint[plotIndex, r][ageLengthIndex], type="h", lwd=5, xaxt="n", ylab="Change in focus group (%)", xlab="", col=TLcolByGroupCode[ageLengthIndex])
par(las=2)
axis(at=seq(1,nag), labels=ageStructuredNames[ageLengthIndex], side=1)

###########################################################
## based on end point, scale each response var by each focus var - just the most extreme one
#############################################################
relativeEndPoint<-0*allEndPoints[,,3]
for(f in 1:ng){
  thisEndPoint<-endPoint[f,3]
  thisResponses<-allEndPoints[f,,3]/thisEndPoint
  relativeEndPoint[f,]<-thisResponses
}
thisMin<-min(relativeEndPoint, na.rm=TRUE); thisMax<-max(relativeEndPoint, na.rm=TRUE)
getColorForMap<-function(x){
  thisCol<-"white"
  if(!is.na(x)){
    if(x<0){
      y<-round(x/thisMin, 2) *100 +1
      thisCol<-negColRamp[y]
    } else{
      y<-round(x/thisMax,2)*100 + 1
      thisCol<-posColRamp[y]
    }
    # 
    # if(x==1){thisCol<-posColRamp[101]}
    # if(x==(-1)){thisCol<-negColRamp[101]}
    if(x==0){thisCol="white"}
    if(x>thisMax){thisCol=posColRamp[101]}
    if(x<thisMin){thisCol=negColRamp[101]}
  }
  return(thisCol)
}

## do heat map
plotColour<-apply(relativeEndPoint[plotIndex, plotIndex][TLindex,TLindex], c(1,2), getColorForMap)
  pdf(paste(plotPath,"RelativeEndPointOrderedByTL.pdf", sep=""), height=8, width=8)
  
  par(mar=c(4,4,1,1))
  tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotColour)[1]),dim(plotColour)[2]),"y"=sort(rep(seq(1,dim(plotColour)[2]),dim(plotColour)[1]))))
  
  par(mar=c(10,10,1.5,1))
  plot(x=seq(1,dim(plotColour)[1]),y=rep(dim(plotColour)[2],dim(plotColour)[1]),type="n",xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,dim(plotColour)[2]))
  mtext("Focus group", side=1, adj=0.5,line=9)
  axis(at=seq(1,dim(plotColour)[1]),labels = labelNames,side=1,las=2)
  axis(at=seq(1,dim(plotColour)[2]),labels=labelNames,side=2,las=1)
  temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
  mtext("Response group", side=2, adj=0.5, line=9)
  box()
  dev.off()

  ##plot legend
  pdf(paste(plotPath,"RelativeEndPointLEGEND.pdf", sep=""), height=4, width=2.5)
  par(mar=c(0,0,0,0))
  legend<-seq(-6,2,by=1); legendColors<-unlist(lapply(legend, getColorForMap))
  makeBlankPlot()
  legend(legend=legend, col=legendColors, x="center", bty="n", pch=15, pt.cex=1.5, title="Relative change \nin biomass")
  dev.off()
  

  thisMin<- -100; thisMax<-100
  for(p in 1:3){
    plotColour<-apply(allEndPoints[plotIndex, plotIndex,p][TLindex,TLindex], c(1,2), getColorForMap)
    pdf(paste(plotPath,"EndPointOrderedByTL",p,".pdf", sep=""), height=8, width=8)
    
    par(mar=c(4,4,1,1))
    tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotColour)[1]),dim(plotColour)[2]),"y"=sort(rep(seq(1,dim(plotColour)[2]),dim(plotColour)[1]))))
    
    par(mar=c(10,10,1.5,1))
    plot(x=seq(1,dim(plotColour)[1]),y=rep(dim(plotColour)[2],dim(plotColour)[1]),type="n",xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,dim(plotColour)[2]))
    mtext("Focus group", side=1, adj=0.5,line=9)
    axis(at=seq(1,dim(plotColour)[1]),labels = labelNames,side=1,las=2)
    axis(at=seq(1,dim(plotColour)[2]),labels=labelNames,side=2,las=1)
    temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
    mtext("Response group", side=2, adj=0.5, line=9)
    box()
    dev.off()
  } 
  
  ##plot legend
  pdf(paste(plotPath,"EndPointLEGEND.pdf", sep=""), height=4, width=2.5)
  par(mar=c(0,0,0,0))
  legend<-seq(-100,100,by=20); legendColors<-unlist(lapply(legend, getColorForMap))
  makeBlankPlot()
  legend(legend=legend, col=legendColors, x="center", bty="n", pch=15, pt.cex=1.5, title="Proportional change \nin biomass (%)")
  dev.off()
  
### same thing, but log scale
  getColorForMapLOG<-function(x){
    thisCol<-"white"
    if(!is.na(x)){
      z<-(log(abs(x), base=10) -log(0.001, base=10))  / (log(100, base=10) - log(0.001, base=10))
      if(z<0){z<-0}
      y<-100 * z + 1
      if(x<0){
        thisCol<-negColRamp[y]
      } else{
         thisCol<-posColRamp[y]
      }
      # 
      # if(x==1){thisCol<-posColRamp[101]}
      # if(x==(-1)){thisCol<-negColRamp[101]}
      if(x==0){thisCol="white"}
      if(x>thisMax){thisCol=posColRamp[101]}
      if(x<thisMin){thisCol=negColRamp[101]}
    }
    return(thisCol)
  }
 
  ##plot legend
  pdf(paste(plotPath,"EndPointLogLEGEND.pdf", sep=""), height=4, width=2.5)
  par(mar=c(0,0,0,0))
  legend<-c(-100,-10,-1,-0.1,-0.01,0,0.01,0.1,1,10,100); legendColors<-unlist(lapply(legend, getColorForMapLOG))
  makeBlankPlot()
  legend(legend=legend, col=legendColors, x="center", bty="n", pch=15, pt.cex=1.5, title="Proportional change \nin biomass (%)")
  dev.off()
 
  for(p in 1:3){
    plotColour<-apply(allEndPoints[plotIndex, plotIndex,p][TLindex,TLindex], c(1,2), getColorForMapLOG)
    pdf(paste(plotPath,"EndPointLogOrderedByTL",p,".pdf", sep=""), height=8, width=8)
    
    par(mar=c(4,4,1,1))
    tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotColour)[1]),dim(plotColour)[2]),"y"=sort(rep(seq(1,dim(plotColour)[2]),dim(plotColour)[1]))))
    
    par(mar=c(10,10,1.5,1))
    plot(x=seq(1,dim(plotColour)[1]),y=rep(dim(plotColour)[2],dim(plotColour)[1]),type="n",xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,dim(plotColour)[2]))
    mtext("Focus group", side=1, adj=0.5,line=9)
    axis(at=seq(1,dim(plotColour)[1]),labels = labelNames,side=1,las=2)
    axis(at=seq(1,dim(plotColour)[2]),labels=labelNames,side=2,las=1)
    temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
    mtext("Response group", side=2, adj=0.5, line=9)
    box()
    dev.off()
  } 
  
## write some stuff out
  this_out<-c("allEndPoints", "endPoint", "groupsDF", "asIndex", "ageStructuredNames", "TLindex", "biomassByFocusGroup")
  ## allEndPoints and endPoint are relative to the end point in the base model (%), so -20 means it ended 20% lower than the base model
  ## allEndPoints has dimensions ngroups, ngroups, nruns where ngroups is the total number of species groups (55), nruns is the number of model runs (3)
  ## first dim is the focus group, second is the response group
  ## because i found it helpful, i also have endPoint which is just the end point for the focus groups for each of the 3 runs (dim = 55, 3)
  ##groupsDF is just the data defining groups, i think you have this already, but i put it here incase its helpful
  ## asIndex is the age structured index which is helpful as I only did this for the age structured groups
  ## TLindex puts age structured groups into trophic level order (roughly, as the trophic levels defined are a bit sketchy)
  ## to get endPoint ordered by TL, do endPoint[asIndex,][TLindex,] (you need the asIndex first)
  ##  biomassByFocusGroup has all the biomass tracers for all groups and all runs including the base run
  ## has dimensions ngroups, nruns, ngroup, ntimesteps (55, 4, 55, 50) 
  ## first dim is focus group, second is runs with 4th the base run and 1-3 increasing additional mortality, 3rd dimension is response group, fith is time (years)
  
  #write them out - commented out so only write out if intending to
  # save(list=this_out,file=paste(dataOutPath,"interactionData",sep=""))
  
########################################################################
#which groups have small changes and big effects..?
storef<-c(); storer<-c()
for(f in 1:ng){
  for(r in 1:3){
    fmin<-endPoint[f,r]
    if(!is.na(fmin)){
      if(fmin<0){
        datamin<-(allEndPoints[f,,r])
        if(fmin<20 & max(datamin, na.rm=TRUE)>50){storef<-c(storef,f); storer<-c(storer,r)}
      }
    }
  }
}

## get out total influence - then rank. 
totalInfl<-array()

plot(fitData1,type="l", ylim=c(min(fitData3), max(fitData1)))
points(fitData2, type="l", lty=2)          
points(fitData3, type="l", lty=3)    

par(lend=1, mar=c(9,4,1,1))
plot(abs(lowestPoint[plotIndex,][TLindex,][,1]), type="h", lwd=5, xaxt="n", ylab="Percentage decrease", xlab="", ylim=c(0,100))
par(las=2)
axis(at=seq(1,length(labelNames)), labels=labelNames, side=1)

plot(abs(lowestPoint[plotIndex,][TLindex,][,2]), type="h", lwd=5, xaxt="n", ylab="Percentage decrease", xlab="", ylim=c(0,100))
par(las=2)
axis(at=seq(1,length(labelNames)), labels=labelNames, side=1)

plot(abs(lowestPoint[plotIndex,][TLindex,][,3]), type="h", lwd=5, xaxt="n", ylab="Percentage decrease", xlab="", ylim=c(0,100))
par(las=2)
axis(at=seq(1,length(labelNames)), labels=labelNames, side=1)

## sim 1's only have one that reduces by at least 20% (IVS) and 2 others that reduced by 5% or more (PFS, EIS)


#############################################
#############################################

## plot DPI and IVH
par(mfrow=c(2,1), las=1, mar=c(4,4,1,1))
thisFCode<-"SND"; thisCode<-"HAK"
thisFCode<-"ASQ"; thisCode<-"IVS"
f<-seq(1,ng)[groupsDF$Code==thisFCode]; g<-seq(1,ng)[groupsDF$Code==thisCode]
thisFNumCohorts<-groupsDF$NumCohorts[f]
thisBaseData<-biomassByFocusGroup[f,2,f,]
alt1Data<-biomassByFocusGroup[f,1,f,]; alt2Data<-biomassByFocusGroup[f,2,f,]; alt3Data<-biomassByFocusGroup[f,3,f,]
thisMax<-max(max(thisBaseData, na.rm=TRUE), max(alt1Data, na.rm=TRUE), max(alt2Data, na.rm=TRUE))
plot(thisBaseData, type="l", xlab="Years simulated after burnin", lwd=2, ylim=c(0,thisMax), ylab="Biomass (tonnes)")
points(alt1Data, type="l", col=myGreen,lty=2, lwd=2)
points(alt2Data, type="l", col=myPurple, lty=4, lwd=2)
points(alt3Data, type="l", col=myBlue, lty=3, lwd=2)
# legend(legend=c("Up", "Middle", "Down"), col=c(myGreen, "black", myPurple), lty=c(2,1,4), lwd=2, seg.len=4,x="bottomright", bty="n")
mtext(thisFCode,side=3, adj=0, font=2)
##
thisBaseData<-biomassByFocusGroup[f,2,g,]
alt1Data<-biomassByFocusGroup[f,1,g,]; alt2Data<-biomassByFocusGroup[f,2,g,]; alt3Data<-biomassByFocusGroup[f,3,g,]
thisMax<-max(max(thisBaseData, na.rm=TRUE), max(alt1Data, na.rm=TRUE), max(alt2Data, na.rm=TRUE))
plot(thisBaseData, type="l", xlab="Years simulated after burnin", lwd=2, ylim=c(0,thisMax), ylab="Biomass (tonnes)")
points(alt1Data, type="l", col=myGreen,lty=2, lwd=2)
points(alt2Data, type="l", col=myPurple, lty=4, lwd=2)
points(alt3Data, type="l", col=myBlue, lty=3, lwd=2)
# legend(legend=c("Up", "Middle", "Down"), col=c(myGreen, "black", myPurple), lty=c(2,1,4), lwd=2, seg.len=4,x="bottomright", bty="n")
mtext(thisCode,side=3, adj=0, font=2)

###############################

# #plot the focus groups biomass tracers for each run
# pdf(paste(plotPath,"focusGroupBiomassUnderSimRunsALLGROUPS.pdf", sep=""))
# par(mfrow=c(5,2), mar=c(4,4,1.5,1), oma=c(0,0,0,0))
# for(f in 1:ng){
#   thisFCode<-as.character(groupsDF$Code[f]); thisFNumCohorts<-groupsDF$NumCohorts[f]
#   if(thisFNumCohorts>1){
#     cat(thisFCode,"--")
#     thisBaseData<-biomassByFocusGroup[f,2,f,]
#     alt1Data<-biomassByFocusGroup[f,1,f,]; alt2Data<-biomassByFocusGroup[f,3,f,]
#     if(sum(thisBaseData, na.rm=TRUE)>0 & sum(alt1Data, na.rm=TRUE)>0 & sum(alt2Data, na.rm=TRUE)>0){
#       thisMax<-max(max(thisBaseData, na.rm=TRUE), max(alt1Data, na.rm=TRUE), max(alt2Data, na.rm=TRUE))
#       plot(thisBaseData, type="l", xlab="Years simulated after burnin", lwd=2, ylim=c(0,thisMax), ylab="Biomass (tonnes)")
#       points(alt1Data, type="l", col=myGreen,lty=2, lwd=2)
#       points(alt2Data, type="l", col=myPurple, lty=4, lwd=2)
#       # legend(legend=c("Up", "Middle", "Down"), col=c(myGreen, "black", myPurple), lty=c(2,1,4), lwd=2, seg.len=4,x="bottomright", bty="n")
#       mtext(thisFCode,side=3, adj=0, font=2)
#     }
#   }
# }
# dev.off()



#plot the focus groups biomass tracers for each run
pdf(paste(plotPath,"focusGroupRelBiomassUnderSimRunsALLGROUPS.pdf", sep=""))
par(mfrow=c(5,2), mar=c(4,4,1.5,1), oma=c(0,0,0,0))
for(f in 1:ng){
  thisFCode<-as.character(groupsDF$Code[f]); thisFNumCohorts<-groupsDF$NumCohorts[f]
  if(thisFNumCohorts>1){
    cat(thisFCode,"--")
    thisBaseData<-biomassByFocusGroup[f,2,f,]
    alt1Data<-biomassByFocusGroup[f,1,f,]; alt2Data<-biomassByFocusGroup[f,3,f,]
    rel1<-alt1Data/thisBaseData; rel2<-alt2Data/thisBaseData
    if(sum(thisBaseData, na.rm=TRUE)>0 & sum(alt1Data, na.rm=TRUE)>0 & sum(alt2Data, na.rm=TRUE)>0){
      thisMax<-max(max(rel1, na.rm=TRUE), max(rel2, na.rm=TRUE))
      plot(rel1, type="l", xlab="Years simulated after burnin", col=myGreen, lty=2, lwd=2, ylim=c(0,thisMax), ylab="Simulation run biomass/base run biomass")
      points(rel2, type="l", col=myPurple, lty=4, lwd=2)
      # legend(legend=c("Up", "Middle", "Down"), col=c(myGreen, "black", myPurple), lty=c(2,1,4), lwd=2, seg.len=4,x="bottomright", bty="n")
      mtext(thisFCode,side=3, adj=0, font=2)
    }
  }
}
dev.off()


baseDiffs<-c(rep(1.4,50), rep(-1.4,50))
testDiffs<-c(rep(-1.4,50), rep(1.4,50))

rcorr(x=baseDiffs, y=testDiffs, type="pearson")
# x  y
# x  1 -1
# y -1  1

testDiffs<-c(rep(1.4,50), runif(50))
rcorr(x=baseDiffs, y=testDiffs, type="pearson")




