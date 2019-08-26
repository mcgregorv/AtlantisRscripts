# plot tracers rel to base runs - tracers were read in and stored in readInChaosRunsAndSaveForAnalyses.R
# have 2 main sets of runs so far - SampleA and SampleKEY - both have fished and non-fished version
# Use DChaosBASE and DChaosFISH to compare
basePath <-paste(DIR$'Base',"AtlantisModels\\chaos\\",sep="")
thisChaosVersion <- "SampleA";
thisChaosVersion <- "FISHSampleA";
thisChaosVersion <- "SampleKEY";
thisChaosVersion <- "FISHSampleKEY";
nChaosRuns <- 10

mg_2_tonne <- 2e-8; X_CN<-5.7

thisChaosVersion <- "SampleKEYPlusBP"
nChaosRuns <- 35

thisCex=1.65

plotPath <- paste(DIR$'Reports',"Chaos\\Figures\\", sep="") ## paper version
plotPath <- paste(DIR$'Figures',"Chaos\\", sep="")

chaosVersions <- c("SampleA", "FISHSampleA", "SampleKEY", "FISHSampleKEY"); nchaosVersions <- length(chaosVersions)
chaosRunDescr <- c("A: Weighted randon - no fishing", "B: Weighted random - with fishing", 
                   "C: Weighted Keystone - no fishing", "D: Weighted keystone - with fishing")

chaosVersions <- c("SampleInformanceAll", "FISHSampleInformanceAll", "SampleKEYPlusBP", "FISHSampleKEYPlusBP"); nchaosVersions <- length(chaosVersions)
chaosRunDescr <- c("A: Weighted informance - no fishing", "B: Weighted informance - with fishing", 
                   "C: Weighted Keystone - no fishing", "D: Weighted keystone - with fishing")

# chaosVersions <- c("DChaosUpDownBASEmodelTracers", "DChaosNtracersFISH", paste("ChaosNtracers",c("SampleInformanceAll", "FISHSampleInformanceAll", "SampleKEYPlusBP", "FISHSampleKEYPlusBP"), sep="")); nchaosVersions <- length(chaosVersions)
# chaosRunDescr <- c("A: All up or down - no fishing", "B: All up or down - with fishing", "C: Weighted randon - no fishing", "D: Weighted random - with fishing", 
#                    "E: Weighted Keystone - no fishing", "F: Weighted keystone - with fishing")
chaosVersions <- c("DChaosNtracers", "DChaosNtracersFISH", paste("ChaosNtracers",c("SampleInformanceAll", "FISHSampleInformanceAll", "SampleKEYPlusBP", "FISHSampleKEYPlusBP"), sep="")); nchaosVersions <- length(chaosVersions)
chaosRunDescr <- c("A: All up or down - no fishing", "B: All up or down - with fishing", "C: Weighted randon - no fishing", "D: Weighted random - with fishing", 
                   "E: Weighted Keystone - no fishing", "F: Weighted keystone - with fishing")



groupsDF <- read.csv(paste(basePath, "..\\CRAM_groups.csv",sep="")); ng <- dim(groupsDF)[1]
nts <- 151; nlayers <-6
thisGroupsDF <- groupsDF[groupsDF$Code != "DC",]; ng <- dim(thisGroupsDF)[1]
keepGindex <- grep("DC", groupsDF$Code, invert = TRUE)

thisYmin <- -100; thisYmax <- 100
allBreaks <-seq(floor((thisYmin/10))*10, ceiling((thisYmax/10))*10, by=10)

# thisYmin <- -70; thisYmax <- 100
# allBreaks <-seq(floor((thisYmin/10))*10, ceiling((thisYmax/10))*10, by=5)
plotVersion<-"A"

v1s <- c(1,0.5,1,0.5); v2s <- v1s-0.5
showTimes <- c(1,50,100,150); nt<-length(showTimes)
showTimeLabels <- c("Initial", "1915", "1956","2015")
f1s <- seq(0,1, length.out=(nt+1))[1:nt]; f2s <- seq(0,1, length.out=(nt+1))[-1]

showTimes <- c(1,150); nt<-length(showTimes)
showTimeLabels <- c("Initial (1865)","2015")

colByTime <- colorRampPalette(colors=c(myBlue, myAqua))(nt); colByTime_trans <- paste(colByTime,"88", sep="")
nruns<-35 # this is the max number of runs in a set if they have different numbers

# load baseArray - has tracers from Chaos BASE and FISH runs (c("outputDChaosFISH", "outputDChaosBASE"))
load(paste(basePath,"ChaosNtracersBASEandBaseFISH_baseAarray",sep=""))

# store the hists
storeHists <- array(NA, dim=c(nchaosVersions, nts, length(allBreaks)-1))
storeAllRel <- array(NA, dim=c(nchaosVersions, nruns, ng, nts))
storeAllAbs <- storeAllRel
for(c in 1:nchaosVersions){
  rm(storeNarray, tracersArray, thisBaseArray)
  thisChaosVersion <- chaosVersions[c]
  this_v1 <- v1s[c]; this_v2 <- v2s[c]
  
  if(length(grep("with fishing",chaosRunDescr[c]))>0){
    thisBaseArray <- baseArray[1,,]
  } else{
    thisBaseArray <- baseArray[2,,]
  }
   # bring in storeNarray
  # load(paste(basePath,"ChaosNtracers",thisChaosVersion,sep=""))
  load(paste(basePath,thisChaosVersion,sep=""))
  # if(!exists("storeNarray")){
  #   storeNarray <- apply(tracersArray, c(1,2,3), FUN=function(x){x*mg_2_tonne*X_CN})
  #   thisBaseArray <- tracersArray[5,,]
  #   storeNarray <- tracersArray[-5,,]
  # }
  if(dim(storeNarray)[2]==(ng+1)){storeNarray<-storeNarray[,keepGindex,]}
  if(dim(thisBaseArray)[1]==(ng+1)){thisBaseArray<-thisBaseArray[keepGindex,]}
  # # TESTING
  # this_g <- 1
  # for(this_g in 1:ng){
  #   thisMax <- max(c(max(storeNarray[,this_g,], na.rm=TRUE), max(thisBaseArray[this_g, ], na.rm=TRUE)), na.rm=TRUE);
  #   thisMin <- min(c(min(storeNarray[,this_g,], na.rm=TRUE), min(thisBaseArray[this_g, ], na.rm=TRUE)), na.rm=TRUE)
  #  par(mfrow=c(2,1), mar=c(4,4,1,1))
  #    plot(storeNarray[1,this_g,], type="l", ylim=c(thisMin, thisMax), col=myGreen_trans, ylab="")
  #   for(r in 1:this_nruns){
  #      points(storeNarray[r,this_g,], type="l", col=myGreen)
  #   }
  #   points(thisBaseArray[this_g, ], type="l", col="red", ylab="")
  #  par(las=1); mtext(paste(this_g,": ",as.character(thisGroupsDF$Code[this_g]), sep=""), side=3, adj=0)
  #   tempRatio <- storeNarray[,this_g,]/thisBaseArray[this_g,]
  #   thisMax <- max(tempRatio, na.rm=TRUE); thisMin <- min(tempRatio, na.rm=TRUE)
  #   plot(tempRatio[1,], type="l", col=myGreen_trans, ylim=c(thisMin, thisMax))
  #   for(r in 1:this_nruns){
  #     points(tempRatio[r,], type="l", col=myGreen)
  #   }
  # 
  # }
  #
  #summarise the largest proportional differences - perhaps from 50 years, 100 years, and at the end (151 years)
  allRelTracers <- 0*storeNarray
  for(g in 1:ng){
    thisCode<-as.character(groupsDF$Code[g]); 
    thisBaseTracer <- thisBaseArray[g,]
    theseChaosTracers <- storeNarray[,g,]
    
    relTracers <- 0*theseChaosTracers; this_nruns <- dim(theseChaosTracers)[1]
    for(r in 1:this_nruns){
      allRelTracers[r,g,]<-(theseChaosTracers[r,] - thisBaseTracer)/thisBaseTracer
      storeAllRel[c,r,g,]<- (theseChaosTracers[r,] - thisBaseTracer)/thisBaseTracer
      storeAllAbs[c,r,g,]<- theseChaosTracers[r,]
    }
  }
  for(t in 1:nts){
    this_t <- t
    thisValues <- 100*as.double(allRelTracers[,,this_t])
    # plus bins
    thisValues[thisValues>thisYmax]<- thisYmax
    thisValues[thisValues<thisYmin]<- thisYmin
    h1 <- hist(thisValues,  breaks=allBreaks, plot=FALSE)
    storeHists[c,t,]<- h1$density
    if(c==1 & t==1){storeMids <- h1$mids}
  }
  rm(storeNarray, tracersArray)
}

overallHists <- array(NA, dim=c(nts, dim(storeHists)[3]))
for(t in 1:nts){
  this_t <- t
  thisValues <- 100*as.double(storeAllRel[,,,this_t])
  # plus bins
  thisValues[thisValues>thisYmax]<- thisYmax
  thisValues[thisValues<thisYmin]<- thisYmin
  h1 <- hist(thisValues,  breaks=allBreaks, plot=FALSE)
  overallHists[t,]<- h1$density
}
storeHists_as <- 0*storeHists
for(c in 1:nchaosVersions){
  for(t in 1:nts){
    this_t <- t
    thisValues <- 100*as.double(asRelTracers[c,,,this_t])
    # plus bins
    thisValues[thisValues>thisYmax]<- thisYmax
    thisValues[thisValues<thisYmin]<- thisYmin
    h1 <- hist(thisValues,  breaks=allBreaks, plot=FALSE)
    storeHists_as[c,t,]<- h1$density
    if(c==1 & t==1){storeMids <- h1$mids}
  }
}
timesteps <- seq(1865,2015)
# par(fig=c(0,1,0.5,1), mar=c(4,4,1,1))
# plot(x=storeMids, y=overallHists[1,], type="h", lwd=10, col=myGrey, ylim=c(0,max(overallHists, na.rm=TRUE)))
# points(x=storeMids, y=overallHists[151,], type="h", lwd=7, col=myAqua)


shiftDown<-0.7
fig_list<-NULL
fig_list[["topleft"]] <- c(0,0.5,0.66*shiftDown,1*shiftDown); fig_list[["topright"]] <- c(0.5,1,0.66*shiftDown,1*shiftDown)
fig_list[["midleft"]] <- c(0,0.5,0.33*shiftDown,0.66*shiftDown); fig_list[["midright"]] <- c(0.5,1,0.33*shiftDown,0.66*shiftDown)
fig_list[["botleft"]] <- c(0,0.5,0,0.33*shiftDown); fig_list[["botright"]] <- c(0.5,1,0,0.33*shiftDown)

pdf(paste(plotPath,"HistCompareChaosRuns.pdf", sep=""), height=7, width=8)
par(mar=c(3,1,0,1), oma=c(0,3,1,1), lend=1, las=1)
par(fig=c(0,1,shiftDown,1))
plot(x=storeMids, y=overallHists[1,], type="h", lwd=20, col=myGrey, ylim=c(0,max(overallHists, na.rm=TRUE)*1.5))
points(x=storeMids, y=overallHists[151,], type="h", lwd=14, col=myAqua)
legend(legend=c("Initial","Final"), col=c(myGrey, myAqua), lwd=c(15,12), x="topright", bty="n", seg.len=3)
mtext("A: Overall", side=3, adj=0.01, line=-1, cex=1.1)
par(fig=fig_list[["topleft"]])
par(new=TRUE)
for(c in 1:nchaosVersions){
  thisLetter <- c("B", "C", "D", "E", "F", "G")[c]
  thisDesc<-paste(thisLetter,": ",unlist(str_split(chaosRunDescr[c],":"))[2], sep="")
  thisYaxt<-"s"
  if(round(c/2)==(c/2)){thisYaxt="n"} #supress for odd
  plot(x=storeMids, y=storeHists[c,1,], type="h", lwd=10, col=myGrey, ylim=c(0,max(overallHists, na.rm=TRUE)*1.5), yaxt=thisYaxt)
  # points(x=storeMids, y=storeHists[c,50,], type="h", lwd=8, col=myBlue)
  points(x=storeMids, y=storeHists[c,151,], type="h", lwd=7, col=myAqua)
  mtext(thisDesc, side=3, adj=0.01, line=-1, cex=0.8)
  #set position for next plot
  if(c<nchaosVersions){
    par(fig=fig_list[[(c+1)]])
    par(new=TRUE)
  }
}
mtext("Relative change (%)", side=1, outer=TRUE, adj=0.5, line=-1.1)
par(las=0)
mtext("Density", side=2, outer=TRUE, adj=0.5, line=2)
dev.off()

pdf(paste(plotPath,"HistCompareChaosRuns_incl50yearsIn.pdf", sep=""), height=7, width=8)
par(mar=c(3,1,0,1), oma=c(0,3,1,1), lend=1, las=1)
par(fig=c(0,1,shiftDown,1))
plot(x=storeMids, y=overallHists[1,], type="h", lwd=20, col=myLightGrey, ylim=c(0,max(overallHists, na.rm=TRUE)*1.5))
points(x=storeMids, y=overallHists[50,], type="h", lwd=14, col=myBlue)
points(x=storeMids, y=overallHists[151,], type="h", lwd=10, col=myAqua)
legend(legend=c("Initial","50 years","Final"), col=c(myLightGrey, myBlue, myAqua), lwd=c(14,11,10), x="topright", bty="n", seg.len=3)
mtext("A: Overall", side=3, adj=0.01, line=0.05, cex=1.1)
par(fig=fig_list[["topleft"]])
par(new=TRUE)
for(c in 1:nchaosVersions){
  thisLetter <- c("B", "C", "D", "E", "F", "G")[c]
  thisDesc<-paste(thisLetter,": ",unlist(str_split(chaosRunDescr[c],":"))[2], sep="")
  thisYaxt<-"s"
  if(round(c/2)==(c/2)){thisYaxt="n"} #supress for odd
  plot(x=storeMids, y=storeHists[c,1,], type="h", lwd=10, col=myLightGrey, ylim=c(0,max(overallHists, na.rm=TRUE)*1.5), yaxt=thisYaxt)
  points(x=storeMids, y=storeHists[c,50,], type="h", lwd=7, col=myBlue)
  points(x=storeMids, y=storeHists[c,151,], type="h", lwd=5, col=myAqua)
  mtext(thisDesc, side=3, adj=0.01, line=0.05, cex=0.8)
  #set position for next plot
  if(c<nchaosVersions){
    par(fig=fig_list[[(c+1)]])
    par(new=TRUE)
  }
}
mtext("Relative change (%)", side=1, outer=TRUE, adj=0.5, line=-1.1)
par(las=0)
mtext("Density", side=2, outer=TRUE, adj=0.5, line=2)
dev.off()

########################################################
# how does it look if we limit to only age-structured groups?
asIndex <- groupsDF$NumCohorts[keepGindex]>1
asRelTracers <- storeAllRel[,,asIndex,]

overallHists_as <- array(NA, dim=c(nts, dim(storeHists)[3]))
for(t in 1:nts){
  this_t <- t
  thisValues <- 100*as.double(asRelTracers[,,,this_t])
  # plus bins
  thisValues[thisValues>thisYmax]<- thisYmax
  thisValues[thisValues<thisYmin]<- thisYmin
  h1 <- hist(thisValues,  breaks=allBreaks, plot=FALSE)
  overallHists_as[t,]<- h1$density
}

pdf(paste(plotPath,"HistCompareChaosRuns_ageStructured.pdf", sep=""), height=7, width=8)
par(mar=c(3,1,0,1), oma=c(0,3,1,1), lend=1, las=1)
par(fig=c(0,1,shiftDown,1))
plot(x=storeMids, y=overallHists_as[1,], type="h", lwd=20, col=myGrey, ylim=c(0,max(overallHists, na.rm=TRUE)*1.5))
points(x=storeMids, y=overallHists_as[151,], type="h", lwd=14, col=myAqua)
legend(legend=c("Initial","Final"), col=c(myGrey, myAqua), lwd=c(15,12), x="topright", bty="n", seg.len=3)
mtext("A: Overall", side=3, adj=0.01, line=-1, cex=1.1)
par(fig=fig_list[["topleft"]])
par(new=TRUE)
for(c in 1:nchaosVersions){
  thisLetter <- c("B", "C", "D", "E", "F", "G")[c]
  thisDesc<-paste(thisLetter,": ",unlist(str_split(chaosRunDescr[c],":"))[2], sep="")
  thisYaxt<-"s"
  if(round(c/2)==(c/2)){thisYaxt="n"} #supress for odd
  plot(x=storeMids, y=storeHists_as[c,1,], type="h", lwd=10, col=myGrey, ylim=c(0,max(overallHists, na.rm=TRUE)*1.5), yaxt=thisYaxt)
  points(x=storeMids, y=storeHists_as[c,151,], type="h", lwd=7, col=myAqua)
  mtext(thisDesc, side=3, adj=0.01, line=-1, cex=0.8)
  #set position for next plot
  if(c<nchaosVersions){
    par(fig=fig_list[[(c+1)]])
    par(new=TRUE)
  }
}
mtext("Relative change (%)", side=1, outer=TRUE, adj=0.5, line=-1.1)
par(las=0)
mtext("Density", side=2, outer=TRUE, adj=0.5, line=2)
dev.off()


## SOME OTHER PLOTS - HEATMAPS AND BOXPLOTS, THAT I DIDN'T END UP USING...
# 
# # turn to heatmap
# colByHeat <- colorRampPalette(colors=c("snow", myGrey, "black"))(101)
# getColor <- function(x, log=FALSE, base=10){
#   thisCol<-"white"
#   if(!is.na(x)){
#     thisIndex <- round(x/thisMax,2)*100 +1
#     if(log==TRUE){
#       thisIndex <- round((log(x, base=base)/log(thisMax, base=base)),2)*100 +1
#     }
#     if(length(thisIndex)>0){
#       if(thisIndex > 101){
#         thisIndex <- 101
#       }
#       if(sum(x, na.rm=TRUE)>0 & thisIndex>0){
#         thisCol<-colByHeat[thisIndex]
#       }
#     }
#   }
#   return(thisCol)
# }
# plotGrid<-function(x,y){
#   thisX<-c(x-0.5,x-0.5,x+0.5,x+0.5); thisY<-c(y-0.5,y+0.5,y+0.5,y-0.5)
#   thisCol<-plotColour[x,y]
#   if(length(thisCol)>0){
#     polygon(x=thisX,y=thisY,col=thisCol,border=NA)
#   }
#   return(NULL)
# }
# thisMax <- max(storeHists, na.rm=TRUE); 
# histByCol <- apply(storeHists, c(1,2,3), getColor, log=TRUE, base=10)
# 
# colByTime <- colorRampPalette(colors=c("Black", "red"))(nt); colByTime_trans<-paste(colByTime,"88", sep="")
# 
# # par(mfrow=c(2,2), mar=c(3,3,1,1))
# 
# fh1ByC <- c(0,0.65,0,0.65); fh2ByC <- c(0.35,1,0.35,1)
# fh3ByC <- c(0.35,0.35,0.35,0.35); fh4ByC <- c(0.65,0.65,0.65,0.65)
# fv1ByC <- c(0.5,0.5,0,0); fv2ByC <- c(1,1,0.5, 0.5);
# fv4ByC <- rev(c(0.25,0.5,0.75,1)); fv3ByC <- rev(c(0,0.25,0.5,0.75))
# fv4ByC <- rev(c(0.2,0.45,0.75,1)); fv3ByC <- rev(c(0,0.25,0.55,0.8))
# 
# endColor <- rep(c(myBlue, myRed),2); endColor_trans <- paste(endColor,"88", sep="")
# endColor_trans <- endColor
# startColor <- "midnightblue"; 
# pdf(paste(plotPath, "VarianceByTimeAndChaosRun",plotVersion,".pdf",sep=""), height=7, width=10)
# for(c in 1:nchaosVersions){
#   h1 <- fh1ByC[c]; h2 <- fh2ByC[c]; h3 <- fh3ByC[c]; h4 <- fh4ByC[c]
#   v1 <- fv1ByC[c]; v2 <- fv2ByC[c]; v3 <- fv3ByC[c]; v4 <- fv4ByC[c]
#   par(fig=c(h1,h2,v1,v2), mar=c(2,2,2,1), oma=c(0,0,0,0))
#   par(lend=1)
#   thisEndColor <- endColor[c]; thisEndColor_trans <- endColor_trans[c]
#   colByTime_trans<-c(startColor, thisEndColor_trans); colByTime <- c(startColor, thisEndColor)
#   if(c>1){par(new=TRUE)}
#   plotColour <- histByCol[c,,]
#   tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotColour)[1]),dim(plotColour)[2]),"y"=sort(rep(seq(1,dim(plotColour)[2]),dim(plotColour)[1]))))
#   # pdf(paste(plotPath,"KemptonsQ_SSRsens.pdf",sep=""),height=4,width=10)
#   plot(x=seq(1,dim(plotColour)[1]),y=rep(dim(plotColour)[2],dim(plotColour)[1]),type="n",xlab="",ylab="",yaxt="n", xaxt="n",ylim=c(0,dim(plotColour)[2]))
#   temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
#   box()
#   par(las=1)
#   zz <- pretty(allBreaks); axisBreaks <- zz[zz<=max(allBreaks) & zz >= min(allBreaks)]; xx <- match(axisBreaks, allBreaks); axisAt <- xx[!is.na(xx)]
#   axis(labels=axisBreaks, at=axisAt, side=2, cex.axis=0.5, cex.lab=0.5)
#   yy <- pretty(seq(1865,2015)); axisYears<-yy[yy %in% seq(1865,2015)]; xx <- match(axisYears, seq(1865,2015)); atYears <- xx[!is.na(xx)]
#   axis(at=atYears, labels = axisYears, side=1, cex.axis=0.5, cex.lab=0.5*thisCex, line=-1, lwd=0)
#   abline(v=atYears, col=myGrey_trans)
#   mtext(chaosRunDescr[c], side=3, adj=0, line=0, cex=thisCex*0.4, font=1)
#   abline(v=1,col=colByTime[1], lty=2, lwd=2.5); abline(v=151, col=colByTime[2], lty=2, lwd=2.5)
#   
#   par(fig=c(h3,h4,v3,v4), mar=c(2,2,0,0), oma=c(1,0,3,0))
#   par(new=TRUE)
#   par(lend=1)
#   plot(x=allBreaks, y=rep(0, length(allBreaks)), type="n", ylim=c(0,750), xlab="Percentage change", cex.axis=0.5, cex.lab=0.5)
#   lwds<-c(7,5)
#   for(t in 1:2){
#     this_t <- showTimes[t]
#     thisValues <- storeHists[c,this_t,]
#     points(x=storeMids, y=thisValues, type="h", lwd=lwds[t], col=colByTime_trans[t])
#   }
# }
# dev.off()
# 
# 
# fh1ByC <- c(0,0.65,0,0.65,0.65); fh2ByC <- c(0.35,1,0.35,1,0.35,1)
# fh3ByC <- c(0.35,0.35,0.35,0.35,0.35,0.35); fh4ByC <- c(0.65,0.65,0.65,0.65,0.65,0.65)
# fv1ByC <- c(0.66,0.66,0.33,0.33,0,0); fv2ByC <- c(1,1,0.66, 0.66,0.33,0.33);
# # fv4ByC <- rev(c(0.25,0.5,0.75,1)); fv3ByC <- rev(c(0,0.25,0.5,0.75))
# fv4ByC <- rev(c(0.3,0.6,0.75,1)); fv3ByC <- rev(c(0,0.17,0.34,0.55, 0.72))
# 
# pdf(paste(plotPath, "VarianceByTimeAndChaosRun6Versions",plotVersion,".pdf",sep=""), height=10, width=10)
# for(c in 1:nchaosVersions){
#   h1 <- fh1ByC[c]; h2 <- fh2ByC[c]; h3 <- fh3ByC[c]; h4 <- fh4ByC[c]
#   v1 <- fv1ByC[c]; v2 <- fv2ByC[c]; v3 <- fv3ByC[c]; v4 <- fv4ByC[c]
#   par(fig=c(h1,h2,v1,v2), mar=c(2,2,2,1), oma=c(0,0,0,0))
#   par(lend=1)
#   thisEndColor <- endColor[c]; thisEndColor_trans <- endColor_trans[c]
#   colByTime_trans<-c(startColor, thisEndColor_trans); colByTime <- c(startColor, thisEndColor)
#   if(c>1){par(new=TRUE)}
#   plotColour <- histByCol[c,,]
#   tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotColour)[1]),dim(plotColour)[2]),"y"=sort(rep(seq(1,dim(plotColour)[2]),dim(plotColour)[1]))))
#   # pdf(paste(plotPath,"KemptonsQ_SSRsens.pdf",sep=""),height=4,width=10)
#   plot(x=seq(1,dim(plotColour)[1]),y=rep(dim(plotColour)[2],dim(plotColour)[1]),type="n",xlab="",ylab="",yaxt="n", xaxt="n",ylim=c(0,dim(plotColour)[2]))
#   temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
#   box()
#   par(las=1)
#   zz <- pretty(allBreaks); axisBreaks <- zz[zz<=max(allBreaks) & zz >= min(allBreaks)]; xx <- match(axisBreaks, allBreaks); axisAt <- xx[!is.na(xx)]
#   axis(labels=axisBreaks, at=axisAt, side=2, cex.axis=0.5, cex.lab=0.5)
#   yy <- pretty(seq(1865,2015)); axisYears<-yy[yy %in% seq(1865,2015)]; xx <- match(axisYears, seq(1865,2015)); atYears <- xx[!is.na(xx)]
#   axis(at=atYears, labels = axisYears, side=1, cex.axis=0.5, cex.lab=0.5*thisCex, line=-1, lwd=0)
#   abline(v=atYears, col=myGrey_trans)
#   mtext(chaosRunDescr[c], side=3, adj=0, line=0, cex=thisCex*0.4, font=1)
#   abline(v=1,col=colByTime[1], lty=2, lwd=2.5); abline(v=151, col=colByTime[2], lty=2, lwd=2.5)
#   
#   par(fig=c(h3,h4,v3,v4), mar=c(2,2,0,0), oma=c(1,0,3,0))
#   par(new=TRUE)
#   par(lend=1)
#   plot(x=allBreaks, y=rep(0, length(allBreaks)), type="n", ylim=c(0,750), xlab="Percentage change", cex.axis=0.5, cex.lab=0.5)
#   lwds<-c(7,5)
#   for(t in 1:2){
#     this_t <- showTimes[t]
#     thisValues <- storeHists[c,this_t,]
#     points(x=storeMids, y=thisValues, type="h", lwd=lwds[t], col=colByTime_trans[t])
#   }
# }
# dev.off()
# 
# 
# ## check out D runs - perhaps can include them here too - this brings in tracersArray
# ## these are all up or all down; no fishing
# load(paste(basePath, "DChaosUpDownBASEmodelTracers", sep="")) # brings in tracersArray
# 
# nChaosRuns <- 9
# baseIndex <- 5
# 
# this_ng <- ng-1 # don't have DC included (carrion) as not used
# #summarise the largest proportional differences - perhaps from 50 years, 100 years, and at the end (151 years)
# upDownRelTracers <- 0*tracersArray
# for(g in 1:this_ng){
#   thisCode<-as.character(groupsDF$Code[g]); 
#   thisBaseTracer <- tracersArray[baseIndex,g,]
#   theseChaosTracers <- tracersArray[,g,]
#   
#   relTracers <- 0*theseChaosTracers; nruns <- dim(theseChaosTracers)[1]
#   for(r in 1:nruns){
#     upDownRelTracers[r,g,]<-(theseChaosTracers[r,] - thisBaseTracer)/thisBaseTracer
#   }
#   
# }
# # the fished version
# load(paste(basePath, "DChaosNtracersFISH", sep="")) #brings in storeNarray
# 
# nChaosRuns <- 9
# 
# baseIndex <- 1
# 
# this_ng <- ng-1 # don't have DC included (carrion) as not used
# #summarise the largest proportional differences - perhaps from 50 years, 100 years, and at the end (151 years)
# upDownFISHRelTracers <- 0*storeNarray
# for(g in 1:this_ng){
#   thisCode<-as.character(groupsDF$Code[g]); 
#   thisBaseTracer <- storeNarray[baseIndex,g,]
#   theseChaosTracers <- storeNarray[,g,]
#   
#   relTracers <- 0*theseChaosTracers; nruns <- dim(theseChaosTracers)[1]
#   for(r in 1:nruns){
#     upDownFISHRelTracers[r,g,]<-(theseChaosTracers[r,] - thisBaseTracer)/thisBaseTracer
#   }
#   
# }
# # pdf(paste(plotPath, "PercentageDifferenceByICRun.pdf", sep=""), height=6.5, width=10)
# pdf(paste(plotPath, "PercentageDifferenceByICRun_35Runs.pdf", sep=""), height=8, width=10)
# par(mfrow=c(3,2), mar=c(3,4.5,1.7,1))
# # the UpDown runs
# allRelTracers <- upDownRelTracers
# thisYmin<-100*min(allRelTracers, na.rm=TRUE); thisYmax <- 100*max(allRelTracers, na.rm=TRUE)
# cat("Min ", thisYmin, " Max ", thisYmax)
# thisYmin <- -30; thisYmax <- 30
# thisOutLine<-FALSE
# plot(0:4, type="n", ylim=c(thisYmin, thisYmax), xaxt="n", xlab="", ylab="Relative difference (%)", cex.lab=thisCex, cex.axis=thisCex)
# abline(h=0, col="red", lty=2, lwd=1.5)
# boxplot(100*as.double(allRelTracers[,,1]), at=1.5,add=TRUE, outline=thisOutLine, col=myBlue_trans, border=myBlue, pch=20, cex=0.5*thisCex, yaxt="n")
# boxplot(100*as.double(allRelTracers[,,50]), at=2.5,add=TRUE, outline=thisOutLine, col=myBlue_trans, border=myBlue, pch=20, cex=0.5*thisCex, yaxt="n")
# boxplot(100*as.double(allRelTracers[,,100]), at=3.5,add=TRUE, outline=thisOutLine, col=myBlue_trans, border=myBlue, pch=20, cex=0.5*thisCex, yaxt="n")
# boxplot(100*as.double(allRelTracers[,,150]), at=4.5,add=TRUE, outline=thisOutLine, col=myBlue_trans, border=myBlue, pch=20, cex=0.5*thisCex, yaxt="n")
# mtext("All up or down - no fishing", side=3, adj=0, line=0.05, cex=thisCex*0.85, font=1)
# axis(at=c(1.5,2.5,3.5,4.5), labels=c("Initial", "1915", "1956","2015"), side=1, cex.axis=thisCex)
# # with fishing
# allRelTracers <- upDownFISHRelTracers
# thisOutLine<-FALSE
# plot(0:4, type="n", ylim=c(thisYmin, thisYmax), xaxt="n", xlab="", ylab="Relative difference (%)", cex.lab=thisCex, cex.axis=thisCex)
# abline(h=0, col="red", lty=2, lwd=1.5)
# boxplot(100*as.double(allRelTracers[,,1]), at=1.5,add=TRUE, outline=thisOutLine, col=myBlue_trans, border=myBlue, pch=20, cex=0.5*thisCex, yaxt="n")
# boxplot(100*as.double(allRelTracers[,,50]), at=2.5,add=TRUE, outline=thisOutLine, col=myBlue_trans, border=myBlue, pch=20, cex=0.5*thisCex, yaxt="n")
# boxplot(100*as.double(allRelTracers[,,100]), at=3.5,add=TRUE, outline=thisOutLine, col=myBlue_trans, border=myBlue, pch=20, cex=0.5*thisCex, yaxt="n")
# boxplot(100*as.double(allRelTracers[,,150]), at=4.5,add=TRUE, outline=thisOutLine, col=myBlue_trans, border=myBlue, pch=20, cex=0.5*thisCex, yaxt="n")
# mtext("All up or down - with fishing", side=3, adj=0, line=0.05, cex=thisCex*0.85, font=1)
# axis(at=c(1.5,2.5,3.5,4.5), labels=c("Initial", "1915", "1956","2015"), side=1, cex.axis=thisCex)
# for(c in 1:nchaosVersions){
#   #summarise the largest proportional differences - perhaps from 50 years, 100 years, and at the end (151 years)
#   allRelTracers <- storeAllRel[c,,,]
#   thisOutLine<-FALSE
#   plot(0:4, type="n", ylim=c(thisYmin, thisYmax), xaxt="n", xlab="", ylab="Relative difference (%)", cex.lab=thisCex, cex.axis=thisCex)
#   abline(h=0, col="red", lty=2, lwd=1.5)
#   boxplot(100*as.double(allRelTracers[,,1]), at=1.5,add=TRUE, outline=thisOutLine, col=myBlue_trans, border=myBlue, pch=20, cex=0.5*thisCex, yaxt="n")
#   boxplot(100*as.double(allRelTracers[,,50]), at=2.5,add=TRUE, outline=thisOutLine, col=myBlue_trans, border=myBlue, pch=20, cex=0.5*thisCex, yaxt="n")
#   boxplot(100*as.double(allRelTracers[,,100]), at=3.5,add=TRUE, outline=thisOutLine, col=myBlue_trans, border=myBlue, pch=20, cex=0.5*thisCex, yaxt="n")
#   boxplot(100*as.double(allRelTracers[,,150]), at=4.5,add=TRUE, outline=thisOutLine, col=myBlue_trans, border=myBlue, pch=20, cex=0.5*thisCex, yaxt="n")
#   mtext(chaosRunDescr[c], side=3, adj=0, line=0.05, cex=thisCex*0.85, font=1)
#   axis(at=c(1.5,2.5,3.5,4.5), labels=c("Initial", "1915", "1956","2015"), side=1, cex.axis=thisCex)
# 
# }
# dev.off()
# 
# 
# # for(c in 1:nchaosVersions){
# #   plotColour <- histByCol[c,,]
# #   tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotColour)[1]),dim(plotColour)[2]),"y"=sort(rep(seq(1,dim(plotColour)[2]),dim(plotColour)[1]))))
# #   # pdf(paste(plotPath,"KemptonsQ_SSRsens.pdf",sep=""),height=4,width=10)
# #   par(mar=c(6,4,1.5,1))
# #   plot(x=seq(1,dim(plotColour)[1]),y=rep(dim(plotColour)[2],dim(plotColour)[1]),type="n",xlab="",ylab="",yaxt="n", xaxt="n",ylim=c(0,dim(plotColour)[2]))
# #   temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
# #   box()
# #   par(las=1)
# #   axisBreaks <- pretty(allBreaks); xx <- match(axisBreaks, allBreaks); axisAt <- xx[!is.na(xx)]
# #   axis(labels=axisBreaks, at=axisAt, side=2)
# #   yy <- pretty(seq(1865,2015)); axisYears<-yy[yy %in% seq(1865,2015)]; xx <- match(axisYears, seq(1865,2015)); atYears <- xx[!is.na(xx)]
# #   axis(at=atYears, labels = axisYears, side=1)
# #   mtext(chaosRunDescr[c], side=3, adj=0, line=0, cex=thisCex*0.4, font=1)
# #   abline(v=1,col=colByTime[1], lty=2); abline(v=151, col=colByTime[2], lty=2)
# #   
# #   plot(x=allBreaks, y=rep(0, length(allBreaks)), type="n", ylim=c(0,850), xlab="Percentage change")
# #   for(t in 1:nt){
# #     
# #     this_t <- showTimes[t]
# #     thisValues <- 100*as.double(allRelTracers[,,this_t])
# #     # plus bins
# #     thisValues[thisValues>thisYmax]<- thisYmax
# #     thisValues[thisValues<thisYmin]<- thisYmin
# #     h1 <- hist(thisValues,  breaks=allBreaks, plot=FALSE)
# #     
# #     hist(thisValues,  breaks=allBreaks, col=colByTime_trans[t], border=colByTime[t], main="",xlab="", yaxt="n", xaxt="n", ylab="", add=TRUE, width=0.1)
# #   }
# #   mtext(chaosRunDescr[c], side=3, adj=0, line=0, cex=thisCex*0.4, font=1)
# #   legend(legend=c("Initial", "Final"), fill=colByTime_trans, pch=15, x="topright", border=colByTime, bty="n", col=NA)
# # }
# # par(las=0)
# # mtext("Percentage change", side=2, outer=TRUE, adj=0.5, line=-1)
# 
# 
# par(mfrow=c(3,3), mar=c(3,4,1.2,0.2))
#   thisYmin<-100*min(allRelTracers, na.rm=TRUE); thisYmax <- 100*max(allRelTracers, na.rm=TRUE)
#   cat("Min ", thisYmin, " Max ", thisYmax)
#   thisYmin <- -50; thisYmax <- 50
#   thisOutLine<-FALSE
#   plot(0:4, type="n", ylim=c(thisYmin, thisYmax), xaxt="n", xlab="", ylab="Relative difference (%)", cex.lab=thisCex, cex.axis=thisCex)
#   abline(h=0, col="red", lty=2, lwd=1.5)
#   boxplot(100*as.double(allRelTracers[,,1]), at=1.5,add=TRUE, outline=thisOutLine, col=myBlue_trans, border=myBlue, pch=20, cex=0.5*thisCex, yaxt="n")
#   boxplot(100*as.double(allRelTracers[,,50]), at=2.5,add=TRUE, outline=thisOutLine, col=myBlue_trans, border=myBlue, pch=20, cex=0.5*thisCex, yaxt="n")
#   boxplot(100*as.double(allRelTracers[,,100]), at=3.5,add=TRUE, outline=thisOutLine, col=myBlue_trans, border=myBlue, pch=20, cex=0.5*thisCex, yaxt="n")
#   boxplot(100*as.double(allRelTracers[,,150]), at=4.5,add=TRUE, outline=thisOutLine, col=myBlue_trans, border=myBlue, pch=20, cex=0.5*thisCex, yaxt="n")
#   mtext(plotLetters[r], side=3, adj=0, line=0.05, cex=thisCex*0.85, font=1)
#   axis(at=c(1.5,2.5,3.5,4.5), labels=c("Initial", "1915", "1956","2015"), side=1, cex.axis=thisCex)
# 
# 
# 
# plotLetters <- c("A", "B", "C", "D", "E", "F", "G", "H","I")
# pdf(paste(plotPath, "AllUpOrDown.pdf", sep=""))
# par(mfrow=c(3,3), mar=c(3,4,1.2,0.2))
# for(r in c(1:4,6:9)){
#   thisYmin<-100*min(allRelTracers, na.rm=TRUE); thisYmax <- 100*max(allRelTracers, na.rm=TRUE)
#   cat("Min ", thisYmin, " Max ", thisYmax)
#   thisYmin <- -50; thisYmax <- 50
#   thisOutLine<-FALSE
#   plot(0:4, type="n", ylim=c(thisYmin, thisYmax), xaxt="n", xlab="", ylab="Relative difference (%)", cex.lab=thisCex, cex.axis=thisCex)
#   abline(h=0, col="red", lty=2, lwd=1.5)
#   boxplot(100*as.double(allRelTracers[r,,1]), at=1.5,add=TRUE, outline=thisOutLine, col=myBlue_trans, border=myBlue, pch=20, cex=0.5*thisCex, yaxt="n")
#   boxplot(100*as.double(allRelTracers[r,,50]), at=2.5,add=TRUE, outline=thisOutLine, col=myBlue_trans, border=myBlue, pch=20, cex=0.5*thisCex, yaxt="n")
#   boxplot(100*as.double(allRelTracers[r,,100]), at=3.5,add=TRUE, outline=thisOutLine, col=myBlue_trans, border=myBlue, pch=20, cex=0.5*thisCex, yaxt="n")
#   boxplot(100*as.double(allRelTracers[r,,150]), at=4.5,add=TRUE, outline=thisOutLine, col=myBlue_trans, border=myBlue, pch=20, cex=0.5*thisCex, yaxt="n")
#   mtext(plotLetters[r], side=3, adj=0, line=0.05, cex=thisCex*0.85, font=1)
#   axis(at=c(1.5,2.5,3.5,4.5), labels=c("Initial", "1915", "1956","2015"), side=1, cex.axis=thisCex)
# }
# dev.off()


### prob delete
# par(mfrow=c(2,2), mar=c(3,3,1,1))
# for(c in 1:4){
#   thisChaosVersion <- chaosVersions[c]
#   this_v1 <- v1s[c]; this_v2 <- v2s[c]
#   
#   groupsDF <- read.csv(paste(basePath, "..\\CRAM_groups.csv",sep="")); ng <- dim(groupsDF)[1]
#   nts <- 151; nlayers <-6
#   
#   # load baseArray - has tracers from Chaos BASE and FISH runs (c("outputDChaosFISH", "outputDChaosBASE"))
#   load(paste(basePath,"ChaosNtracersBASEandBaseFISH_baseAarray",sep=""))
#   
#   # bring in storeNarray
#   load(paste(basePath,"ChaosNtracers",thisChaosVersion,sep=""))
#   
#   #summarise the largest proportional differences - perhaps from 50 years, 100 years, and at the end (151 years)
#   allRelTracers <- 0*storeNarray
#   for(g in 1:ng){
#     thisCode<-as.character(groupsDF$Code[g]); 
#     thisBaseTracer <- baseArray[1,g,]
#     theseChaosTracers <- storeNarray[,g,]
#     
#     relTracers <- 0*theseChaosTracers; nruns <- dim(theseChaosTracers)[1]
#     for(r in 1:nruns){
#       allRelTracers[r,g,]<-(theseChaosTracers[r,] - thisBaseTracer)/thisBaseTracer
#     }
#   }
#   plot(1, type="n", ylim=c(0,850), xlim=c(thisYmin, thisYmax) ,xlab="Percentage change")
#   for(t in 1:nt){
#     
#     this_t <- showTimes[t]
#     thisValues <- 100*as.double(allRelTracers[,,this_t])
#     # plus bins
#     thisValues[thisValues>thisYmax]<- thisYmax
#     thisValues[thisValues<thisYmin]<- thisYmin
#     h1 <- hist(thisValues,  breaks=allBreaks, plot=FALSE)
#     
#     hist(thisValues,  breaks=allBreaks, col=colByTime_trans[t], border=colByTime[t], main="",xlab="", yaxt="n", xaxt="n", ylab="", add=TRUE, width=0.1)
#     par(xpd=TRUE)
#   }
#   mtext(chaosRunDescr[c], side=3, adj=0, line=0, cex=thisCex*0.4, font=1)
# }


# pdf(paste(plotPath, "PercentageDifferenceByICRun_35RunsHIST.pdf", sep=""), height=6.5, width=10)
# # par(mfrow=c(2,2), mar=c(3,4.2,1.3,1))
# par(mfrow=c(1,1), mar=c(3,4.2,1.3,1))
# for(c in 1:nchaosVersions){
#   thisChaosVersion <- chaosVersions[c]
#   this_v1 <- v1s[c]; this_v2 <- v2s[c]
#   
#   groupsDF <- read.csv(paste(basePath, "..\\CRAM_groups.csv",sep="")); ng <- dim(groupsDF)[1]
#   nts <- 151; nlayers <-6
#   
#   # load baseArray - has tracers from Chaos BASE and FISH runs (c("outputDChaosFISH", "outputDChaosBASE"))
#   load(paste(basePath,"ChaosNtracersBASEandBaseFISH_baseAarray",sep=""))
#   
#   # bring in storeNarray
#   load(paste(basePath,"ChaosNtracers",thisChaosVersion,sep=""))
#   
#   #summarise the largest proportional differences - perhaps from 50 years, 100 years, and at the end (151 years)
#   allRelTracers <- 0*storeNarray
#   for(g in 1:ng){
#     thisCode<-as.character(groupsDF$Code[g]); 
#     thisBaseTracer <- baseArray[1,g,]
#     theseChaosTracers <- storeNarray[,g,]
#     
#     relTracers <- 0*theseChaosTracers; nruns <- dim(theseChaosTracers)[1]
#     for(r in 1:nruns){
#       allRelTracers[r,g,]<-(theseChaosTracers[r,] - thisBaseTracer)/thisBaseTracer
#     }
#   }
#   # thisYmin<-100*min(allRelTracers, na.rm=TRUE); thisYmax <- 100*max(allRelTracers, na.rm=TRUE)
#   
#   cat("Min ", thisYmin, " Max ", thisYmax)
#   if(c %in% c(1,2)){
#     this_f1s <- f1s/2; this_f2s <- f2s/2
#   } else{
#     this_f1s <- 0.5+f1s/2; this_f2s <- 0.5+f2s/2
#   }
#   
#   for(t in 1:nt){
#     if(t==1){
#       par(fig=c(min(this_f1s),max(this_f2s),this_v2,this_v1), mar=c(0,0,0,0))
#       if(c > 1){
#         par(new=TRUE)
#       }
#       makeBlankPlot()
#       mtext(chaosRunDescr[c], side=3, adj=0, line=-1, cex=thisCex*0.7, font=1)
#     }
#     this_t <- showTimes[t]
#     thisValues <- 100*as.double(allRelTracers[,,this_t])
#     # plus bins
#     thisValues[thisValues>thisYmax]<- thisYmax
#     thisValues[thisValues<thisYmin]<- thisYmin
#     
#     t1<-hist(thisValues, plot=FALSE, breaks=allBreaks)
#     this_x <- t1$counts 
#     
#     HBreaks <- t1$breaks
#     HBreak1 <- t1$breaks[1]
#     hpos <- function(Pos) (Pos-HBreak1)*(length(HBreaks)-1)/ diff(range(HBreaks))
#     this_ylim<-hpos(c(thisYmin, thisYmax))
#     f1<-this_f1s[t]; f2 <- this_f2s[t]
#     par(fig=c(f1,f2,this_v2,this_v1), mar=c(4,1,1,1), oma=c(0,2,0,0))
#     
#     
#     # par(fig=c(f1,f2,0,1), mar=c(4,4,1,1))
#     # if(t>1 | c>1){ par(new=TRUE) }
#     par(new=TRUE)
#     barplot(this_x, space=0, horiz=TRUE, ylim=this_ylim, xaxt="n")
#     axis(at=min(this_x), labels = showTimeLabels[t], side=1)
#     if(c %in% c(1,2) & t==1){
#       par(las=1)
#       axis(at=hpos(t1$mids), labels=t1$mids, side=2)
#     }
#   }
#   # mtext(chaosRunDescr[c], side=3, adj=0, line=-1, cex=thisCex*0.7, font=1, outer=TRUE)
#   # axis(at=seq(1,nt), labels=c("Initial", "1915", "1956","2015"), side=1, cex.axis=thisCex)
#   
# }
# dev.off()
# 
# 

# pdf(paste(plotPath, "PercentageDifferenceByICRun_35RunsNoOutliers.pdf", sep=""), height=6.5, width=10)
# par(mfrow=c(2,2), mar=c(3,4.2,1.3,1))
# for(c in 1:nchaosVersions){
#   thisChaosVersion <- chaosVersions[c]
#   
#   
#   groupsDF <- read.csv(paste(basePath, "..\\CRAM_groups.csv",sep="")); ng <- dim(groupsDF)[1]
#   nts <- 151; nlayers <-6
#   
#   # load baseArray - has tracers from Chaos BASE and FISH runs (c("outputDChaosFISH", "outputDChaosBASE"))
#   load(paste(basePath,"ChaosNtracersBASEandBaseFISH_baseAarray",sep=""))
#   
#   # bring in storeNarray
#   load(paste(basePath,"ChaosNtracers",thisChaosVersion,sep=""))
#   
#   #summarise the largest proportional differences - perhaps from 50 years, 100 years, and at the end (151 years)
#   allRelTracers <- 0*storeNarray
#   for(g in 1:ng){
#     thisCode<-as.character(groupsDF$Code[g]); 
#     thisBaseTracer <- baseArray[1,g,]
#     theseChaosTracers <- storeNarray[,g,]
#     
#     relTracers <- 0*theseChaosTracers; nruns <- dim(theseChaosTracers)[1]
#     for(r in 1:nruns){
#       allRelTracers[r,g,]<-(theseChaosTracers[r,] - thisBaseTracer)/thisBaseTracer
#     }
#   }
#   thisYmin<--50; thisYmax <- 50
#   cat("Min ", thisYmin, " Max ", thisYmax)
#   thisOutLine<-FALSE
#   plot(0:4, type="n", ylim=c(thisYmin, thisYmax), xaxt="n", xlab="", ylab="Relative difference (%)", cex.lab=thisCex, cex.axis=thisCex)
#   abline(h=0, col="red", lty=2, lwd=1.5)
#   boxplot(100*as.double(allRelTracers[,,1]), at=1.5,add=TRUE, outline=thisOutLine, col=myBlue_trans, border=myBlue, pch=20, cex=0.5*thisCex, yaxt="n")
#   boxplot(100*as.double(allRelTracers[,,50]), at=2.5,add=TRUE, outline=thisOutLine, col=myBlue_trans, border=myBlue, pch=20, cex=0.5*thisCex, yaxt="n")
#   boxplot(100*as.double(allRelTracers[,,100]), at=3.5,add=TRUE, outline=thisOutLine, col=myBlue_trans, border=myBlue, pch=20, cex=0.5*thisCex, yaxt="n")
#   boxplot(100*as.double(allRelTracers[,,150]), at=4.5,add=TRUE, outline=thisOutLine, col=myBlue_trans, border=myBlue, pch=20, cex=0.5*thisCex, yaxt="n")
#   mtext(chaosRunDescr[c], side=3, adj=0, line=0.05, cex=thisCex*0.85, font=1)
#   axis(at=c(1.5,2.5,3.5,4.5), labels=c("Initial", "1915", "1956","2015"), side=1, cex.axis=thisCex)
#   
# }
# dev.off()
