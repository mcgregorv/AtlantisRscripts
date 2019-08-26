this_out<-"TestsSCA5"; runFolder="TBGB_JP2"
thisDesc <- paste(runFolder, this_out,sep="")

this_path = paste(DIR$'Base', "TBGB\\",runFolder,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")

plotYearly<-FALSE

biolLines <- readLines(paste(this_path,"TBGB_biol.prm", sep=""))

#read in catches, as fomated and such in summaryCatchByFleet_TBGB.R
load(paste(this_path,"Catch_history\\storeCatchByFleet", sep=""))
nfg <- dim(fishedGroups)[1]

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol <- ncvar_get(ThisNC.nc, "volume"); nts <- dim(thisVol)[3]
daysTimeStep<-73
numStepsPerYear<-365/daysTimeStep
year0<-1865
fishingStartYear<-1900
modelStartYear<-1900
newVol <- ncvar_get(newNC.nc, "volume")

burnin <- (modelStartYear-year0)*numStepsPerYear

xLabsTemp<-seq(0,(nts*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0
modelTime <- year0 + (1:nts)/numStepsPerYear


lastYear <- max(modelTime)
if(plotYearly==TRUE){
plotYearlyIndex <- seq(1, nts, by=numStepsPerYear)
plotPath<-paste(DIR$'Base',"TBGB\\Figures\\Testing\\FIshing\\plotYearly_",sep="")
} else {
plotYearlyIndex <- seq(1,nts)
plotPath<-paste(DIR$'Base',"TBGB\\Figures\\Testing\\FIshing\\",sep="")
}

# read in trawl survey estimates
tsData <- read.csv("C:\\Projects\\2017\\ATLANTISChathamRiseRepository\\TBGB\\Data\\SurveysAndSuch\\AbsoluteBiomass_surveys.csv")
surveyYears<-as.double(gsub("X","",colnames(tsData))); surveyYears<-surveyYears[!is.na(surveyYears)]
## grab the SSB for each group - this is the adult biomass

get_age_mat<-function(x){
  #x is species group code
  thisVar<-paste(x,"_age_mat", sep="")
  y<-grep(thisVar,biolLines); z<-get_first_number(biolLines[y])
  return(z)
}
## store SSB tracers for each group -only fill with age-structured
storeSSB<-array(NA, dim=c(ng,nts)); storeBiomass <- 0*storeSSB
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
  thisName <- str_trim(groupsDF$Name[g], side="both")
  if(thisNumCohorts>1){
    cat(as.character(thisCode))
    thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode],side="both")
    ##get age mature so know which cohorts are part of SSB
    this_age_mat<-get_age_mat(thisCode)
    #if only 2 cohorts, just the adults are SSB
    if(thisNumCohorts==2){
      thisVar<-paste(thisName,"2_Nums", sep="");  fishNums<-ncvar_get(ThisNC.nc,thisVar);
      thisVar<-paste(thisName,"2_ResN", sep="");  fishResN<-ncvar_get(ThisNC.nc,thisVar);
      thisVar<-paste(thisName,"2_StructN", sep="");  fishStructN<-ncvar_get(ThisNC.nc,thisVar);
      tempFishSSB<-apply(fishNums*(fishResN + fishStructN),3,sum)*mg_2_tonne*X_CN;
      storeSSB[g,]<-tempFishSSB
    } else{
      thisSSB<-rep(0,nts); fishSSB<-thisSSB
      for(c in (this_age_mat+1):thisNumCohorts){
        thisVar<-paste(thisName,c,"_Nums", sep="");  fishNums<-ncvar_get(ThisNC.nc,thisVar);
        thisVar<-paste(thisName,c,"_ResN", sep=""); fishResN<-ncvar_get(ThisNC.nc,thisVar);
        thisVar<-paste(thisName,c,"_StructN", sep=""); fishStructN<-ncvar_get(ThisNC.nc,thisVar);
        
        tempFishSSB<-apply(fishNums*(fishResN + fishStructN),3,sum)*mg_2_tonne*X_CN;
        fishSSB<-fishSSB + tempFishSSB
      }
      storeSSB[g,]<-fishSSB
    }
    # total biomass
    thisTracer<-paste(thisName,"_N", sep=""); thisData <- ncvar_get(ThisNC.nc, thisTracer)
    thisBiomass <- apply(thisData * thisVol, 3, sum) * mg_2_tonne * X_CN
    storeBiomass[g,] <- thisBiomass
  }
}

## SCA only check
g=grep("SCA", fishedCodes)
thisCode <- as.character(fishedGroups$Code[g]); gIndex <- groupsDF$Code==thisCode
thisTracerName <- str_trim(fishedGroups$Name[g], side="both"); thisName <- gsub("_", " ", thisTracerName)
thisCatchHistory <- storeCatchByStock[g,]; thisF <- TBGB_F[g,]; 
thisFaxis <- pretty(seq(0,max(thisF), length.out=10)); thisCatchAxis <- pretty(0:max(thisCatchHistory))
thisSSB <- tsData[grep(thisCode,tsData$species),2:(length(surveyYears)+1)]
# get model SSB (roughly - inclu juv for now)
modelSSB <- storeSSB[gIndex,]
SCAPlotYearlyIndex <- seq(5, nts, by=numStepsPerYear)
modelSSB2plot <- modelSSB[SCAPlotYearlyIndex]
modelTime2plot <- modelTime[SCAPlotYearlyIndex]
index1 <- as.double(thisSSB[1,]); time1 <- as.double(gsub("X","",names(thisSSB)))[!is.na(index1)]
mean1 <- mean(index1, na.rm=TRUE); meanS1 <- mean(modelSSB2plot[round(modelTime2plot) %in% time1])
sIndex1 <- (index1/mean1)*meanS1
# second one
index2 <- as.double(thisSSB[2,]); time2 <- as.double(gsub("X","",names(thisSSB)))[!is.na(index2)]
mean2 <- mean(index2, na.rm=TRUE); meanS2 <- mean(modelSSB2plot[round(modelTime2plot) %in% time2])
sIndex2 <- (index2/mean2)*meanS2
thisBmax <- max(c(sIndex1, sIndex2, modelSSB2plot[modelTime2plot>=thisStartYear]), na.rm=TRUE)
pdf(paste(plotPath, "TBGB_SCAonlyCatchCheck_", thisDesc,".pdf", sep=""), height=3, width=5.5)
par( mar=c(3,4.5,1.5,5), lend=1,  par(las=0))
plot(x=catchYears, y=thisCatchHistory, type="h", col=myGrey_trans, lwd=3, yaxt="n", xaxt="n", ylab="", xlab="", xlim=c(1960, max(catchYears)))
axis(at=thisCatchAxis, labels=thisCatchAxis, side=4, line=1, col.axis=myGrey, col.ticks=myGrey, col=myGrey)
par(new=TRUE)
plot(x=modelTime2plot, y=modelSSB2plot, type="l", col=myBlue, xlim=c(1960, 2015), ylim=c(0, max(modelSSB2plot[modelTime2plot>1955])), xlab="", ylab="Biomass (tonnes)")
     points(x=time1, y=sIndex1[!is.na(index1)],pch=4, lwd=1.2, ylab="Biomass (tonnes)", xlab="", ylim=c(0, thisBmax*1.2), xlim=c(thisStartYear, max(catchYears)), cex.axis=thisCex, cex.lab=thisCex)
mtext("Catch (tonnes)", side=4, line=3, col=myGrey, font=2)
mtext(thisName, side=3, adj=0, font=2)
dev.off()

###########################
# plot them to check look right
pdf(paste(plotPath, "TBGB_CatchHistoriesAndF_fromModel",thisDesc,"_withTS.pdf", sep=""), height=10, width=10)
par(mfrow=c(5,2), mar=c(4,4.5,1.5,8), lend=1,  par(las=0))
for(g in 1:nfg){
  thisCode <- as.character(fishedGroups$Code[g]); gIndex <- groupsDF$Code==thisCode
  thisTracerName <- str_trim(fishedGroups$Name[g], side="both"); thisName <- gsub("_", " ", thisTracerName)
  thisCatchHistory <- storeCatchByStock[g,]; thisF <- TBGB_F[g,]; 
  thisFaxis <- pretty(seq(0,max(thisF), length.out=10)); thisCatchAxis <- pretty(0:max(thisCatchHistory))
  
  thisSSB <- tsData[grep(thisCode,tsData$species),2:(length(surveyYears)+1)]
  
  # get model SSB (roughly - inclu juv for now)
  modelSSB <- storeSSB[gIndex,]
  modelSSB2plot <- modelSSB[plotYearlyIndex]
  modelTime2plot <- modelTime[plotYearlyIndex]
  # if(dim(thisSSB)[1]==2){
  #   test<-apply(thisSSB[1:2,],2,sum, na.rm=TRUE)
  #   thisStartYear <- max(c(catchYears, min(gsub("X","",names(test)))))
  # else if(dim(thisSSB)[1]==1){
  if(dim(thisSSB)[1]>0){
    thisStartYear <- min(c(min(catchYears[(thisCatchHistory>0)]), min(as.double(gsub("X","",names(thisSSB))))))
  } else{
    thisStartYear <- min(catchYears[(thisCatchHistory>0)])
 }

  plot(x=catchYears, y=thisCatchHistory, type="h", col=myGrey_trans, lwd=3, yaxt="n", xaxt="n", ylab="", xlab="", xlim=c(thisStartYear, max(catchYears)))
  axis(at=thisCatchAxis, labels=thisCatchAxis, side=4, line=1, col.axis=myGrey, col.ticks=myGrey, col=myGrey)
  par(new=TRUE)
  if(dim(thisSSB)[1]>0){
    thisBiomassMax <- max(c(max(thisSSB, na.rm=TRUE),  max(modelSSB2plot[modelTime2plot>=modelStartYear], na.rm=TRUE), na.rm=TRUE))
    if(dim(thisSSB)[1]==2){
      index1 <- as.double(thisSSB[1,]); time1 <- as.double(gsub("X","",names(thisSSB)))[!is.na(index1)]
      mean1 <- mean(index1, na.rm=TRUE); meanS1 <- mean(modelSSB2plot[round(modelTime2plot) %in% time1])
      sIndex1 <- (index1/mean1)*meanS1
      # second one
      index2 <- as.double(thisSSB[2,]); time2 <- as.double(gsub("X","",names(thisSSB)))[!is.na(index2)]
      mean2 <- mean(index2, na.rm=TRUE); meanS2 <- mean(modelSSB2plot[round(modelTime2plot) %in% time2])
      sIndex2 <- (index2/mean2)*meanS2
      thisBmax <- max(c(sIndex1, sIndex2, modelSSB2plot[modelTime2plot>=thisStartYear]), na.rm=TRUE)
      plot(x=time1, y=sIndex1[!is.na(index1)],pch=4, lwd=1.2, ylab="Biomass (tonnes)", xlab="", ylim=c(0, thisBmax*1.2), xlim=c(thisStartYear, max(catchYears)), cex.axis=thisCex, cex.lab=thisCex)
      # points(x=time2, y=sIndex2[!is.na(index2)],pch=1, lwd=1.2, ylab="Biomass (tonnes)", xlab="", ylim=c(0, thisBiomassMax*1.1), xlim=c(modelStartYear, max(catchYears)), cex.axis=thisCex, cex.lab=thisCex)
    }else{
      index1 <- as.double(thisSSB); time1 <- as.double(gsub("X","",names(thisSSB)))[!is.na(index1)]
      mean1 <- mean(index1, na.rm=TRUE); meanS1 <- mean(modelSSB2plot[round(modelTime2plot) %in% time1])
      sIndex1 <- (index1/mean1)*meanS1
      thisBmax <- max(c(sIndex1,  modelSSB2plot[modelTime2plot>=thisStartYear]), na.rm=TRUE)
      plot(x=time1, y=sIndex1[!is.na(index1)],pch=4, lwd=1.2, ylab="Biomass (tonnes)", xlab="", ylim=c(0, thisBmax*1.2), xlim=c(thisStartYear, max(catchYears)), cex.axis=thisCex, cex.lab=thisCex)
     }
  }else{
    thisBiomassMax <- max(c(max(modelSSB2plot[modelTime2plot>=modelStartYear], na.rm=TRUE), na.rm=TRUE))
    plot(x=surveyYears, y=rep(0, length(surveyYears)),type="n", lty=2, lwd=2, ylab="Biomass (tonnes)", xlab="", ylim=c(0, thisBiomassMax*1.1), xlim=c(thisStartYear, max(catchYears)), cex.axis=thisCex, cex.lab=thisCex)
    
  }
  points(x=modelTime2plot, y=modelSSB2plot, type="l", col=myBlue, lwd=1, lty=1)
  par(las=0)
  mtext("Catch (tonnes)", side=4, line=3, col=myGrey, font=2)
  
  # 
  # par(new=TRUE)
  # plot(x=catchYears, y=thisF, type="l", col=myOrange, lwd=2, yaxt="n", xaxt="n", ylab="", xlab="")
  # axis(at=thisFaxis, labels=thisFaxis, side=4, line=5, col.axis=myOrange, col.ticks=myOrange, col=myOrange)
  # mtext("F", side=4, line=7, col=myOrange, font=2)
  
  mtext(thisName, side=3, adj=0, font=2)

}
dev.off()


plotStartYear<-1960

# plot them to check look right
pdf(paste(plotPath, "TBGBSSB_fromModel",thisDesc,"_withTS.pdf", sep=""), height=10, width=10)
par(mfrow=c(5,2), mar=c(4,4.5,1.5,8), lend=1,  par(las=0))
for(g in 1:nfg){
  thisCode <- as.character(fishedGroups$Code[g]); gIndex <- groupsDF$Code==thisCode
  thisTracerName <- str_trim(fishedGroups$Name[g], side="both"); thisName <- gsub("_", " ", thisTracerName)
  thisCatchHistory <- storeCatchByStock[g,]; thisF <- TBGB_F[g,]; 
  thisFaxis <- pretty(seq(0,max(thisF), length.out=10)); thisCatchAxis <- pretty(0:max(thisCatchHistory))
  
  thisSSB <- tsData[grep(thisCode,tsData$species),2:(length(surveyYears)+1)]
  
  # get model SSB (roughly - inclu juv for now)
  modelSSB <- storeSSB[gIndex,]
  
  modelSSB2plot <- modelSSB[plotYearlyIndex]
  modelTime2plot <- modelTime[plotYearlyIndex]
  
  plot(x=catchYears, y=thisCatchHistory, type="h", col=myGrey_trans, lwd=3, yaxt="n", xaxt="n", ylab="", xlab="", xlim=c(plotStartYear, max(catchYears)))
  axis(at=thisCatchAxis, labels=thisCatchAxis, side=4, line=1, col.axis=myGrey, col.ticks=myGrey, col=myGrey)
  par(new=TRUE)
  
   if(dim(thisSSB)[1]>0){
    thisBiomassMax <- max(c(max(thisSSB, na.rm=TRUE),  max(modelSSB2plot[modelTime2plot>=modelStartYear], na.rm=TRUE), na.rm=TRUE))
    if(dim(thisSSB)[1]==2){
      plot(x=surveyYears, y=thisSSB[1,],pch=4, lwd=1.2, ylab="Biomass (tonnes)", xlab="", xlim=c(plotStartYear, max(catchYears)), ylim=c(0, thisBiomassMax*1.2),  cex.axis=thisCex, cex.lab=thisCex)
      # points(x=surveyYears, y=thisSSB[2,],pch=1, lwd=1.2, ylab="Biomass (tonnes)", xlab="", ylim=c(0, thisBiomassMax*1.2), xlim=c(modelStartYear, max(catchYears)), cex.axis=thisCex, cex.lab=thisCex)
    }else{
      plot(x=surveyYears, y=thisSSB,pch=4, lwd=1.2, ylab="Biomass (tonnes)", xlab="", xlim=c(plotStartYear, max(catchYears)), ylim=c(0, thisBiomassMax*1.2), cex.axis=thisCex, cex.lab=thisCex)
    }
  }else{
    thisBiomassMax <- max(c(max(modelSSB2plot[modelTime2plot>=modelStartYear], na.rm=TRUE), na.rm=TRUE))
    plot(x=surveyYears, y=rep(0, length(surveyYears)),type="n", lty=2, lwd=2, ylab="Biomass (tonnes)", xlab="", ylim=c(0, thisBiomassMax*1.2), xlim=c(plotStartYear, max(catchYears)), cex.axis=thisCex, cex.lab=thisCex)
    
  }
  points(x=modelTime2plot, y=modelSSB2plot, type="l", col=myBlue, lwd=1, lty=1)

  mtext(thisName, side=3, adj=0, font=2)
}
dev.off()


plotStartYear<-1900

# plot them to check look right - alt version using total biomass
pdf(paste(plotPath, "TBGBBiomass_fromModel",thisDesc,"_withTS.pdf", sep=""), height=10, width=10)
par(mfrow=c(5,2), mar=c(4,4.5,1.5,8), lend=1,  par(las=0))
for(g in 1:nfg){
  thisCode <- as.character(fishedGroups$Code[g]); gIndex <- groupsDF$Code==thisCode
  thisTracerName <- str_trim(fishedGroups$Name[g], side="both"); thisName <- gsub("_", " ", thisTracerName)
  thisCatchHistory <- storeCatchByStock[g,]; thisF <- TBGB_F[g,]; 
  thisFaxis <- pretty(seq(0,max(thisF), length.out=10)); thisCatchAxis <- pretty(0:max(thisCatchHistory))
  
  thisSSB <-tsData[grep(thisCode,tsData$species),2:(length(surveyYears)+1)]
  
  # get model SSB (roughly - inclu juv for now)
  modelSSB <- storeBiomass[gIndex,]
  
  modelSSB2plot <- modelSSB[plotYearlyIndex]
  modelTime2plot <- modelTime[plotYearlyIndex]
  
  plot(x=catchYears, y=thisCatchHistory, type="h", col=myGrey_trans, lwd=3, yaxt="n", xaxt="n", ylab="", xlab="", xlim=c(plotStartYear, max(catchYears)))
  axis(at=thisCatchAxis, labels=thisCatchAxis, side=4, line=1, col.axis=myGrey, col.ticks=myGrey, col=myGrey)
  par(new=TRUE)
  
  if(dim(thisSSB)[1]>0){
    thisBiomassMax <- max(c(max(thisSSB, na.rm=TRUE),  max(modelSSB2plot[modelTime2plot>=modelStartYear], na.rm=TRUE), na.rm=TRUE))
    if(dim(thisSSB)[1]==2){
      plot(x=surveyYears, y=thisSSB[1,],pch=4, lwd=1.2, ylab="Biomass (tonnes)", xlab="", xlim=c(plotStartYear, max(catchYears)), ylim=c(0, thisBiomassMax*1.2),  cex.axis=thisCex, cex.lab=thisCex)
      # points(x=surveyYears, y=thisSSB[2,],pch=1, lwd=1.2, ylab="Biomass (tonnes)", xlab="", ylim=c(0, thisBiomassMax*1.2), xlim=c(modelStartYear, max(catchYears)), cex.axis=thisCex, cex.lab=thisCex)
    }else{
      plot(x=surveyYears, y=thisSSB,pch=4, lwd=1.2, ylab="Biomass (tonnes)", xlab="", xlim=c(plotStartYear, max(catchYears)), ylim=c(0, thisBiomassMax*1.2), cex.axis=thisCex, cex.lab=thisCex)
    }
  }else{
    thisBiomassMax <- max(c(max(modelSSB2plot[modelTime2plot>=modelStartYear], na.rm=TRUE), na.rm=TRUE))
    plot(x=surveyYears, y=rep(0, length(surveyYears)),type="n", lty=2, lwd=2, ylab="Biomass (tonnes)", xlab="", ylim=c(0, thisBiomassMax*1.2), xlim=c(plotStartYear, max(catchYears)), cex.axis=thisCex, cex.lab=thisCex)
    
  }
  points(x=modelTime2plot, y=modelSSB2plot, type="l", col=myBlue, lwd=1, lty=1)
  
  mtext(thisName, side=3, adj=0, font=2)
}
dev.off()

