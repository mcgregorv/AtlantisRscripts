# check out resp to IC perturbations for BP, ALL, and AS groups
# bring in c("storeRelResponse", "storeRatioResponse", "storeAllTracersFromALlRuns"), created in modelRelationshipsAt50years_groupStructure_vs_stability.R
basePath <-paste(DIR$'Base',"AtlantisModels\\chaos\\",sep="")
plotPath <- paste(DIR$'Reports',"Chaos\\Figures\\", sep="") ## paper version
plotPath <- paste(DIR$'Figures',"Chaos\\", sep="")

source(paste(DIR$'General functions',"chaosFunctions.R", sep=""))

load(paste(basePath, "AllTracersFromAllCHAOSRuns", sep=""))

groupsDF<- read.csv(paste(basePath,"..\\CRAM_groups.csv", sep="")); ng<-dim(groupsDF)[1]
thisGroupsDF <- groupsDF[groupsDF$Code != "DC",]; ng <- dim(thisGroupsDF)[1]
keepGindex <- grep("DC", groupsDF$Code, invert = TRUE)
thisCex=1.65

plotPath <- paste(DIR$'Reports',"Chaos\\Figures\\", sep="") ## paper version
plotPath <- paste(DIR$'Figures',"Chaos\\", sep="")

chaosVersions <- c("DChaosNtracers", "DChaosNtracersFISH", paste("ChaosNtracers",c("SampleInformanceAll", "FISHSampleInformanceAll", "SampleKEYPlusBP", "FISHSampleKEYPlusBP"), sep="")); nchaosVersions <- length(chaosVersions)
chaosRunDescr <- c("A: All up or down - no fishing", "B: All up or down - with fishing", "C: Weighted randon - no fishing", "D: Weighted random - with fishing", 
                   "E: Weighted Keystone - no fishing", "F: Weighted keystone - with fishing")

nruns <- 35
nts <- 151; nlayers <-6

# load baseArray - has tracers from Chaos BASE and FISH runs (c("outputDChaosFISH", "outputDChaosBASE"))
load(paste(basePath,"ChaosNtracersBASEandBaseFISH_baseAarray",sep=""))
keepBaseArray<- baseArray[,keepGindex,]

colByFished<-rep(c(myOrange_trans, myBlue_trans), 3)
allTime<-1865:2015


cvForNotFIshed <- apply(storeAllTracersFromALlRuns[seq(1, 5, by=2),,,], c(3,4), calc_cv)
cvForFIshed <- apply(storeAllTracersFromALlRuns[seq(2, 6, by=2),,,], c(3,4), calc_cv)


medCVfished <- apply(cvForFIshed, 2, median, na.rm=TRUE)
medCVNotFished <- apply(cvForNotFIshed, 2, median, na.rm=TRUE)
plot(medCVfished, type="l", ylim=c(0, 0.35))
points(medCVNotFished, type="l", lty=2)
abline(v=35, col="red")

BPindex <- thisGroupsDF$NumCohorts==1
BPgroups <- thisGroupsDF$Code[BPindex]; this_ng <- length(BPgroups)

medCVfished <- apply(cvForFIshed[BPindex,], 2, median, na.rm=TRUE)
medCVNotFished <- apply(cvForNotFIshed[BPindex,], 2, median, na.rm=TRUE)
plot(medCVfished, type="l", ylim=c(0, 0.35))
points(medCVNotFished, type="l", lty=2)
abline(v=35, col="red")

medCVfished <- apply(cvForFIshed[!BPindex,], 2, median, na.rm=TRUE)
medCVNotFished <- apply(cvForNotFIshed[!BPindex,], 2, median, na.rm=TRUE)
plot(medCVfished, type="l", ylim=c(0, 0.35))
points(medCVNotFished, type="l", lty=2)
abline(v=35, col="red")


axisCex<-0.8
options(scipen = -2)
pdf(paste(plotPath,"BiomassPoolTracersALL.pdf", sep=""), height=10, width=10)
par(mfrow=c(5,2), mar=c(3.5,4.2,1.5,4), oma=c(0,1,0,0))
for(g in 1:this_ng){
  thisCode <- BPgroups[g]
  gIndex <- grep(thisCode, thisGroupsDF$Code);
  if(length(gIndex>1)){gIndex <- gIndex[thisGroupsDF$Code[gIndex]==thisCode]}
  thisName<- gsub("_", " ", str_trim(thisGroupsDF$Name[gIndex]))
  thisMax <- max(c(max(storeAllTracersFromALlRuns[,,gIndex,], na.rm=TRUE), max(keepBaseArray[,gIndex,])))
  thisBaseTracer <- keepBaseArray[2, gIndex,1:dim(keepBaseArray)[3]]; thisBaseFishTracer <- keepBaseArray[1,gIndex,1:dim(keepBaseArray)[3]]
  par(las=1)
  plot(x=allTime[1:nts],y=thisBaseTracer, ylim=c(0,thisMax*1.1), type="n", ylab="", xlab="", cex=0.5, cex.axis=axisCex)
  par(las=0)
  mtext("Biomass (tonnes)", side=2, adj=0.5, line=3.5, cex=0.7)
  for(c in 1:dim(storeAllTracersFromALlRuns)[1]){
    for(r in 1:dim(storeAllTracersFromALlRuns)[2]){
      points(x=allTime[1:nts],y=storeAllTracersFromALlRuns[c,r,gIndex,1:dim(keepBaseArray)[3]], col=colByFished[c], type="l")
    }
  }
  par(las=0)
  mtext(thisName, side=3, adj=0, cex.lab=axisCex)
  thisFishCV <- cvForFIshed[gIndex,]; thisNotFishedCV <- cvForNotFIshed[gIndex,]; thisMax <- max(c(thisFishCV, thisNotFishedCV), na.rm=TRUE)
  thisYaxis <- pretty(seq(0, 1, length.out=5))
  par(new=TRUE)
  par(las=1)
  plot(x=allTime, y=thisFishCV, ylim=c(0, 1), cex=1, pch=4, col=myAqua,  ylab="", xlab="", yaxt="n", xaxt="n")
  points(x=allTime, y=thisNotFishedCV, cex=1, pch=3, col=myRed)
  axis(at=thisYaxis, labels = 100*thisYaxis, side=4, cex.lab=axisCex, cex.axis=axisCex); 
  par(las=0)
  mtext("CV (%)", side=4, adj=0.5, line=2, cex=0.5)
}
makeBlankPlot()
par(lend=1)
legend(legend=c("Fished models","Not fished models", "Fished CVs", "Not fished CVs"), x="center", col=c(myBlue, myOrange, myAqua, myRed), 
       pch=c(NA, NA, 4, 3), lty=c(1,1,NA, NA), lwd=c(2,2,1,1), bty="n", cex=thisCex)
dev.off()


ASgroups <- thisGroupsDF$Code[!BPindex]; this_ng <- length(ASgroups)
pdf(paste(plotPath,"AgeStructuredTracersALL_p1.pdf", sep=""), height=10, width=10)
par(mfrow=c(5,2), mar=c(3.5,4.2,1.5,4), oma=c(0,1,0,0))
for(g in 1:this_ng){
  thisCode <- ASgroups[g]
  gIndex <- grep(thisCode, thisGroupsDF$Code);
  if(length(gIndex>1)){gIndex <- gIndex[thisGroupsDF$Code[gIndex]==thisCode]}
  thisName<- gsub("_", " ", str_trim(thisGroupsDF$Name[gIndex]))
  thisMax <- max(c(max(storeAllTracersFromALlRuns[,,gIndex,], na.rm=TRUE), max(keepBaseArray[,gIndex,])))
  thisBaseTracer <- keepBaseArray[2, gIndex,1:dim(keepBaseArray)[3]]; thisBaseFishTracer <- keepBaseArray[1,gIndex,1:dim(keepBaseArray)[3]]
  par(las=1)
  plot(x=allTime[1:nts],y=thisBaseTracer, ylim=c(0,thisMax*1.1), type="n", ylab="", xlab="", cex=0.5, cex.axis=axisCex)
  par(las=0)
  mtext("Biomass (tonnes)", side=2, adj=0.5, line=3.5, cex=0.7)
  for(c in 1:dim(storeAllTracersFromALlRuns)[1]){
    for(r in 1:dim(storeAllTracersFromALlRuns)[2]){
      points(x=allTime[1:nts],y=storeAllTracersFromALlRuns[c,r,gIndex,1:dim(keepBaseArray)[3]], col=colByFished[c], type="l")
    }
  }
  # points(x=allTime[1:nts],y=thisBaseTracer, type="l", col="red")
  # points(x=allTime[1:nts],y=thisBaseFishTracer, type="l", lty=2, col="midnightblue")
  par(las=0)
  mtext(thisName, side=3, adj=0, cex.lab=axisCex)
  thisFishCV <- cvForFIshed[gIndex,]; thisNotFishedCV <- cvForNotFIshed[gIndex,]; thisMax <- max(c(thisFishCV, thisNotFishedCV), na.rm=TRUE)
  thisYaxis <- pretty(seq(0, 1, length.out=5))
  par(new=TRUE)
  par(las=1)
  plot(x=allTime, y=thisFishCV, ylim=c(0, 1), pch=4, col=myAqua, cex=1, ylab="", xlab="", yaxt="n", xaxt="n")
  points(x=allTime, y=thisNotFishedCV, pch=3, col=myRed,cex=1)
  axis(at=thisYaxis, labels = 100*thisYaxis, side=4, cex.lab=axisCex, cex.axis=axisCex); 
  par(las=0)
  mtext("CV (%)", side=4, adj=0.5, line=2, cex=0.5)
}
makeBlankPlot()
par(lend=1)
legend(legend=c("Fished models","Not fished models", "Fished CVs", "Not fished CVs"), x="center", col=c(myBlue, myOrange, myAqua, myRed), 
       pch=c(NA, NA, 4, 3), lty=c(1,1,NA, NA), lwd=c(2,2,1,1), bty="n")
dev.off()

axisCex<-0.8
ASgroups <- thisGroupsDF$Code[!BPindex]; this_ng <- length(ASgroups)
pdf(paste(plotPath,"AgeStructuredTracersALL_p2.pdf", sep=""), height=8, width=10)
par(mfrow=c(4,4), mar=c(3.5,4.2,1.5,4), oma=c(0,1,0,0))
for(g in 24:this_ng){
  thisCode <- ASgroups[g]
  gIndex <- grep(thisCode, thisGroupsDF$Code);
  if(length(gIndex>1)){gIndex <- gIndex[thisGroupsDF$Code[gIndex]==thisCode]}
  thisName<- gsub("_", " ", str_trim(thisGroupsDF$Name[gIndex]))
  thisMax <- max(c(max(storeAllTracersFromALlRuns[,,gIndex,], na.rm=TRUE), max(keepBaseArray[,gIndex,])))
  thisBaseTracer <- keepBaseArray[2, gIndex,1:dim(keepBaseArray)[3]]; thisBaseFishTracer <- keepBaseArray[1,gIndex,1:dim(keepBaseArray)[3]]
  par(las=1)
  plot(x=allTime[1:nts],y=thisBaseTracer, ylim=c(0,thisMax*1.1), type="n", ylab="", xlab="", cex=0.5, cex.axis=axisCex)
  par(las=0)
  mtext("Biomass (tonnes)", side=2, adj=0.5, line=3.5, cex=0.7)
  for(c in 1:dim(storeAllTracersFromALlRuns)[1]){
    for(r in 1:dim(storeAllTracersFromALlRuns)[2]){
      points(x=allTime[1:nts],y=storeAllTracersFromALlRuns[c,r,gIndex,1:dim(keepBaseArray)[3]], col=colByFished[c], type="l")
    }
  }
  # points(x=allTime[1:nts],y=thisBaseTracer, type="l", col="red")
  # points(x=allTime[1:nts],y=thisBaseFishTracer, type="l", lty=2, col="midnightblue")
  par(las=0)
  mtext(thisName, side=3, adj=0, cex.lab=axisCex)
  thisFishCV <- cvForFIshed[gIndex,]; thisNotFishedCV <- cvForNotFIshed[gIndex,]; thisMax <- max(c(thisFishCV, thisNotFishedCV), na.rm=TRUE)
  thisYaxis <- pretty(seq(0, 1, length.out=5))
  par(new=TRUE)
  par(las=1)
  plot(x=allTime, y=thisFishCV, ylim=c(0, 1), pch=4, col=myAqua, cex=0.2, ylab="", xlab="", yaxt="n", xaxt="n")
  points(x=allTime, y=thisNotFishedCV, pch=3, col=myRed,cex=0.2)
  axis(at=thisYaxis, labels = 100*thisYaxis, side=4, cex.lab=axisCex, cex.axis=axisCex); 
  par(las=0)
  mtext("CV (%)", side=4, adj=0.5, line=2, cex=0.5)
}
makeBlankPlot()
par(lend=1)
legend(legend=c("Fished models","Not fished models", "Fished CVs", "Not fished CVs"), x="center", col=c(myBlue, myOrange, myAqua, myRed), 
       pch=c(NA, NA, 4, 3), lty=c(1,1,NA, NA), lwd=c(2,2,1,1), bty="n")
dev.off()
# plot seperately examples to have in main body of text
axisCex=1
pdf(paste(plotPath, "CVsByTimeExamples.pdf", sep=""), height=10, width=8)
par(mfrow=c(2,1), mar=c(4,6,1.5,5), oma=c(5,0,0,0))
thisCode <- "DR"
gIndex <- grep(thisCode, thisGroupsDF$Code);
if(length(gIndex>1)){gIndex <- gIndex[thisGroupsDF$Code[gIndex]==thisCode]}
thisName<- gsub("_", " ", str_trim(thisGroupsDF$Name[gIndex]))
thisMax <- max(c(max(storeAllTracersFromALlRuns[,,gIndex,], na.rm=TRUE), max(keepBaseArray[,gIndex,])))
thisBaseTracer <- keepBaseArray[2, gIndex,1:dim(keepBaseArray)[3]]; thisBaseFishTracer <- keepBaseArray[1,gIndex,1:dim(keepBaseArray)[3]]
par(las=1)
plot(x=allTime[1:nts],y=thisBaseTracer, ylim=c(0,thisMax*1), type="n", ylab="", xlab="", cex=0.5, cex.axis=axisCex)
for(c in 1:dim(storeAllTracersFromALlRuns)[1]){
  for(r in 1:dim(storeAllTracersFromALlRuns)[2]){
    points(x=allTime[1:nts],y=storeAllTracersFromALlRuns[c,r,gIndex,1:dim(keepBaseArray)[3]], col=colByFished[c], type="l")
  }
}
par(las=0)
mtext("Biomass (tonnes)", side=2, adj=0.5, line=4)
par(las=0)
mtext(paste("A: ", thisName, sep=""), side=3, adj=0, cex.lab=axisCex)
thisFishCV <- cvForFIshed[gIndex,]; thisNotFishedCV <- cvForNotFIshed[gIndex,]; thisMax <- max(c(thisFishCV, thisNotFishedCV), na.rm=TRUE)
thisYaxis <- pretty(seq(0, 1, length.out=5))
par(new=TRUE)
par(las=1)
plot(x=allTime, y=thisFishCV, ylim=c(0, 1), pch=4, col=myAqua, cex=1, ylab="", xlab="", yaxt="n", xaxt="n")
points(x=allTime, y=thisNotFishedCV, pch=3, col=myRed,cex=1)
axis(at=thisYaxis, labels = 100*thisYaxis, side=4, cex.lab=axisCex, cex.axis=axisCex); 
par(las=0)
mtext("CV (%)", side=4, adj=0.5, line=2.5, cex=1)
## hoki
thisCode <- "HOK"
gIndex <- grep(thisCode, thisGroupsDF$Code);
if(length(gIndex>1)){gIndex <- gIndex[thisGroupsDF$Code[gIndex]==thisCode]}
thisName<- gsub("_", " ", str_trim(thisGroupsDF$Name[gIndex]))
thisMax <- max(c(max(storeAllTracersFromALlRuns[,,gIndex,], na.rm=TRUE), max(keepBaseArray[,gIndex,])))
thisBaseTracer <- keepBaseArray[2, gIndex,1:dim(keepBaseArray)[3]]; thisBaseFishTracer <- keepBaseArray[1,gIndex,1:dim(keepBaseArray)[3]]
par(las=1)
plot(x=allTime[1:nts],y=thisBaseTracer, ylim=c(0,thisMax*1), type="n", ylab="", xlab="", cex=0.5, cex.axis=axisCex)
for(c in 1:dim(storeAllTracersFromALlRuns)[1]){
  for(r in 1:dim(storeAllTracersFromALlRuns)[2]){
    points(x=allTime[1:nts],y=storeAllTracersFromALlRuns[c,r,gIndex,1:dim(keepBaseArray)[3]], col=colByFished[c], type="l")
  }
}
par(las=0)
mtext("Biomass (tonnes)", side=2, adj=0.5, line=4)
par(las=0)
mtext(paste("B: ", thisName, sep=""), side=3, adj=0, cex.lab=axisCex)
thisFishCV <- cvForFIshed[gIndex,]; thisNotFishedCV <- cvForNotFIshed[gIndex,]; thisMax <- max(c(thisFishCV, thisNotFishedCV), na.rm=TRUE)
thisYaxis <- pretty(seq(0, 1, length.out=5))
par(new=TRUE)
par(las=1)
plot(x=allTime, y=thisFishCV, ylim=c(0, 1), pch=4, col=myAqua, cex=1, ylab="", xlab="", yaxt="n", xaxt="n")
points(x=allTime, y=thisNotFishedCV, pch=3, col=myRed,cex=1)
axis(at=thisYaxis, labels = 100*thisYaxis, side=4, cex.lab=axisCex, cex.axis=axisCex); 
par(las=0)
mtext("CV (%)", side=4, adj=0.5, line=2.5, cex=1)
par(xpd="NA")
par(lend=1)
legend(legend=c("Fished models","Not fished models", "Fished CVs", "Not fished CVs"), x="bottom", inset=-0.4, col=c(myBlue, myOrange, myAqua, myRed), 
       pch=c(NA, NA, 4, 3), lty=c(1,1,NA, NA), lwd=c(2,2,1.5,1.5), bty="n", ncol=2)
dev.off()


## do CVs by time, for all species, for BP species, and for AS species
thisMaxCV <- 0.6
plot(x=allTime[1:nts],y=rep(0,nts), type="n", ylim=c(0, thisMaxCV), xlab="", ylab="CV (%)")
for(t in 1:nts){
  d1<-cvForNotFIshed[,t]; d2 <- cvForFIshed[,t]
  boxplot(d1, at=allTime[t]-0.1, add=TRUE, col=myOrange_trans, border=myOrange, outline=FALSE, pch=20, cex=0.2)
  boxplot(d2, at=allTime[t]+0.1, add=TRUE, col=myBlue_trans, border=myBlue, outline=FALSE, pch=20, cex=0.2)
}

plot(x=allTime[1:nts],y=rep(0,nts), type="n", ylim=c(0, thisMaxCV), xlab="", ylab="CV (%)")
for(t in 1:nts){
  d1<-cvForNotFIshed[BPindex,t]; d2 <- cvForFIshed[BPindex,t]
  boxplot(d1, at=t-0.1, add=TRUE, col=myOrange_trans, border=myOrange, outline=FALSE, pch=20, cex=0.2)
  boxplot(d2, at=t+0.1, add=TRUE, col=myBlue_trans, border=myBlue, outline=FALSE, pch=20, cex=0.2)
}
# pdf(paste(plotPath, "CVbyTime.pdf", sep=""), height=4, width=12)
# par(mfrow=c(1,2), mar=c(4,4,1.5,1), oma=c(0,0,0,8))
# plot(x=allTime[1:nts],y=rep(0,nts), type="n", ylim=c(0, thisMaxCV), xlab="", ylab="CV (%)")
# for(t in 1:nts){
#   d1<-cvForNotFIshed[!BPindex,t]; d2 <- cvForFIshed[!BPindex,t]
#   boxplot(d1, at=allTime[t]-0.1, add=TRUE, col=myOrange_trans, border=myOrange, outline=FALSE, pch=20, cex=0.2, xaxt="n", yaxt="n")
#   boxplot(d2, at=allTime[t]+0.1, add=TRUE, col=myBlue_trans, border=myBlue, outline=FALSE, pch=20, cex=0.2, xaxt="n", yaxt="n")
# }
# mtext("A: Age-structured species groups", side=3, adj=0)
# plot(x=allTime[1:nts],y=rep(0,nts), type="n", ylim=c(0, thisMaxCV), xlab="", ylab="CV (%)")
# for(t in 1:nts){
#   d1<-cvForNotFIshed[BPindex,t]; d2 <- cvForFIshed[BPindex,t]
#   boxplot(d1, at=allTime[t]-0.1, add=TRUE, col=myOrange_trans, border=myOrange, outline=FALSE, pch=20, cex=0.2, xaxt="n", yaxt="n")
#   boxplot(d2, at=allTime[t]+0.1, add=TRUE, col=myBlue_trans, border=myBlue, outline=FALSE, pch=20, cex=0.2, xaxt="n", yaxt="n")
# }
# mtext("B: Biomass-pool species groups", side=3, adj=0)
# par(xpd=NA, lend=1)
# legend(legend=c("Fished", "Not fished"), col=c(myBlue, myOrange), seg.len=2.5, lwd=5, x="right", inset=-0.4, bty="n")
# dev.off()

## alt: just plot med, CIs for each
calcCI <- function(x, k){
  y <- sort(x); n <- length(x); i <- max(1,round(n * k))
  thisOut <- y[i]
  return(thisOut)
}
thisCI <- 0.5 # 0.5 gives upper and lower quartiles (0.25 and 0.75) (approximately because of rounding the position here)
# age-structured unfished
ASmed_unfished <- apply(cvForNotFIshed[!BPindex,],2,median, na.rm=TRUE)
ASlq_unfished <- apply(cvForNotFIshed[!BPindex,],2,calcCI, k = (1-thisCI)*0.5); ASuq_unfished <- apply(cvForNotFIshed[!BPindex,],2,calcCI, k = (1-(1-thisCI)*0.5))
# age-structured fished
ASmed_fished <- apply(cvForFIshed[!BPindex,],2,median, na.rm=TRUE)
ASlq_fished <- apply(cvForFIshed[!BPindex,],2,calcCI, k = (1-thisCI)*0.5); ASuq_fished <- apply(cvForFIshed[!BPindex,],2,calcCI, k = (1-(1-thisCI)*0.5))
##
# age-structured unfished
BPmed_unfished <- apply(cvForNotFIshed[BPindex,],2,median, na.rm=TRUE)
BPlq_unfished <- apply(cvForNotFIshed[BPindex,],2,calcCI, k = (1-thisCI)*0.5); BPuq_unfished <- apply(cvForNotFIshed[BPindex,],2,calcCI, k = (1-(1-thisCI)*0.5))
# age-structured fished
BPmed_fished <- apply(cvForFIshed[BPindex,],2,median, na.rm=TRUE)
BPlq_fished <- apply(cvForFIshed[BPindex,],2,calcCI, k = (1-thisCI)*0.5); BPuq_fished <- apply(cvForFIshed[BPindex,],2,calcCI, k = (1-(1-thisCI)*0.5))

thisMaxCV <- max(c(BPuq_fished, BPuq_unfished, ASuq_fished, ASuq_unfished))*1.1

pdf(paste(plotPath, "CVbyTime.pdf", sep=""), height=4, width=12)
par(mfrow=c(1,2), mar=c(4,4,1.5,1), oma=c(0,0,0,8))
## AS
plot(x=allTime[1:nts],y=rep(0,nts), type="n", ylim=c(0, thisMaxCV), xlab="", ylab="CV (%)")
points(x=allTime, y=ASmed_unfished, type="l", lwd=2, col=myOrange)
points(x=allTime, y=ASlq_unfished, type="l", lwd=1.5, col=myOrange, lty=4)
points(x=allTime, y=ASuq_unfished, type="l", lwd=1.5, col=myOrange, lty=4)
# fished
points(x=allTime, y=ASmed_fished, type="l", lwd=2, col=myBlue)
points(x=allTime, y=ASlq_fished, type="l", lwd=1.5, col=myBlue, lty=4)
points(x=allTime, y=ASuq_fished, type="l", lwd=1.5, col=myBlue, lty=4)
mtext("A: Age-structured species groups", side=3, adj=0)
## BP
plot(x=allTime[1:nts],y=rep(0,nts), type="n", ylim=c(0, thisMaxCV), xlab="", ylab="CV (%)")
points(x=allTime, y=BPmed_unfished, type="l", lwd=2, col=myOrange)
points(x=allTime, y=BPlq_unfished, type="l", lwd=1.5, col=myOrange, lty=4)
points(x=allTime, y=BPuq_unfished, type="l", lwd=1.5, col=myOrange, lty=4)
# fished
points(x=allTime, y=BPmed_fished, type="l", lwd=2, col=myBlue)
points(x=allTime, y=BPlq_fished, type="l", lwd=1.5, col=myBlue, lty=4)
points(x=allTime, y=BPuq_fished, type="l", lwd=1.5, col=myBlue, lty=4)
mtext("B: Biomass-pool species groups", side=3, adj=0)
par(xpd=NA, lend=1)
legend(legend=c("Fished", "Not fished"), col=c(myBlue, myOrange), seg.len=2.5, lwd=2.5, x="right", inset=-0.4, bty="n")
dev.off()



plot(x=1:nts, y=rep(0,nts), type="n", ylim=c(0, thisMaxCV))
for(t in 1:nts){
  d1<-cvForNotFIshed[!BPindex,t]; d2 <- cvForNotFIshed[BPindex,t]
  boxplot(d1, at=t-0.1, add=TRUE, col=myGrey_trans, border=myAqua, outline=FALSE, pch=20, cex=0.2)
  boxplot(d2, at=t+0.1, add=TRUE, col=myPurple_trans, border=myPurple, outline=FALSE, pch=20, cex=0.2)
}

#####################################################################################
## fit models to CVs
######################################################################################
# first, set up the data frame that has explanatory variables as well as CVs
allTheRankings <- read.csv(paste(DIR$'Tables', "allTheRankings.csv", sep="")); nrankings <- dim(allTheRankings)[2]
rankDF<-read.csv(paste(DIR$'Tables',"interactionEffectsRANKINGs.csv",sep=""))
thisCex=1.65

## these are the factors that can be used for all groups - can redefine later for age-structured factors
thisFactors <- c("ChaosAlt", "TL",  "NumL1cons", "Lifespan",  "B0", "PropByTopPrey") # note, perturbMethod has chaosrun
# note Linf  cuts out whales and such, might want to take it out

allInts <- getIntVars(thisFactors); all_factors <- c(thisFactors, allInts)



chaosVersions <- c("DChaosNtracers", "DChaosNtracersFISH", paste("ChaosNtracers",c("SampleInformanceAll", "FISHSampleInformanceAll", "SampleKEYPlusBP", "FISHSampleKEYPlusBP"), sep="")); nchaosVersions <- length(chaosVersions)
chaosRunDescr <- c("A: All up or down - no fishing", "B: All up or down - with fishing", "C: Weighted randon - no fishing", "D: Weighted random - with fishing", 
                   "E: Weighted Keystone - no fishing", "F: Weighted keystone - with fishing")

test_timeSteps<-36:nts; ntts <- length(test_timeSteps)
allData2fit <- data.frame(matrix(NA, ncol=(nrankings+ntts+2), nrow=(nchaosVersions * ng)))
thisColnames <- c("ChaosAlt",  "Informance", colnames(allTheRankings),paste("ts_",test_timeSteps, sep=""))
colnames(allData2fit)<- thisColnames
chaosLetters <- c("A", "B", "C", "D", "E", "F")
for(c in 1:nchaosVersions){
    for(g in 1:ng){
      thisRow <- ng*(c-1)  + g
      # cat("c", c,", g, ",g, ", thisRow, ", thisRow,"\n") # always worth checking ;-)
      allData2fit$ChaosAlt[thisRow] <- chaosLetters[c];
      thisSimData <- storeAllTracersFromALlRuns[c,,g,]; thisCVs <- apply(thisSimData, 2, calc_cv)
      allData2fit[thisRow, grep("^ts", colnames(allData2fit))] <- thisCVs[test_timeSteps]

      thisCode <- as.character(thisGroupsDF$Code[g])
      allData2fit[thisRow,c(colnames(allTheRankings))] <- allTheRankings[allTheRankings$Code == thisCode,]
      allData2fit$Code[thisRow] <- thisCode
      thisNumCohorts <- thisGroupsDF$NumCohorts[g]
      if(thisNumCohorts==1){
        allData2fit$Keystone[thisRow] <- 0
      }
      thisInf <- rankDF$informPerformRank[rankDF$Code==thisCode]
      if(length(thisInf)==0){thisInf <- 16}
      allData2fit$Informance[thisRow] <- thisInf
    }
}


breaks <- seq(0,1,by=0.1)
colByTime <- paste(colorRampPalette(colors=c(myLightBlue, myBlue, "midnightblue"))(ntts), "88", sep="")
h1 <- hist(allData2fit[,c("ts_36")], plot=FALSE, breaks=breaks)
plot(x=h1$mids, y=h1$counts, col=colByTime[1], pch=8, lwd=2, ylim=c(0, 500))
for(t in 2:ntts){
  h1 <- hist(allData2fit[,c(paste("ts_",test_timeSteps[t], sep=""))], plot=FALSE, breaks=breaks)
  points(x=h1$mids, y=h1$counts, col=colByTime[t], pch=8, lwd=2)
  
}
storeModelDFs <- NULL; storeFinalR <- list()
for(t in 1:ntts){
  thisRespVar <- paste("ts_",test_timeSteps[t],sep="")
  this_df <- getBestLModel(thisVar = thisRespVar, thisFactors = all_factors, cutoff = 0.1, data2fit = allData2fit)
  storeModelDFs[[t]]<-this_df
  storeFinalR[[t]]<-max(this_df$Rsq)
  cat(t,", ")
}

#split by fished or not fished
storeUnfishedModelDFs <- NULL; storeUnfishedFinalR <- list(); storeUnfishedNfactors<-list()
allLogData <- allData2fit; allLogData[,grep("^ts", colnames(allLogData))]<- (allData2fit[,grep("^ts", colnames(allData2fit))])^(1/3)
# check chaosalt is still a factor, make it so in any case
allLogData$ChaosAlt <- factor(allLogData$ChaosAlt, levels=c("A", "B", "C", "D", "E", "F"))

unfishedIndex <- seq(1,6,by=2)
allUnfishedLogData <- allLogData[allLogData$ChaosAlt %in% chaosLetters[unfishedIndex],]
for(t in 1:ntts){
  thisRespVar <- paste("ts_",test_timeSteps[t],sep="")
  this_df <- getBestLModel(thisVar = thisRespVar, thisFactors = all_factors, cutoff = 0.1, 
                           data2fit =allUnfishedLogData)
  storeUnfishedModelDFs[[t]]<-this_df
  storeUnfishedFinalR[[t]]<-max(this_df$Rsq); storeUnfishedNfactors[[t]]<-dim(this_df)[1]
  cat(t,", ")
}
plot(unlist(storeUnfishedFinalR))
plot(unlist(storeUnfishedNfactors))

test<-unlist(storeUnfishedModelDFs)
table(names(test))
table(test[grep("Factor",names(test))])

# fished
# the logged version
allFishedLogData <- allLogData[!(allLogData$ChaosAlt %in% chaosLetters[unfishedIndex]),]
storeFishedModelDFs <- NULL; storeFishedFinalR <- list(); storeNfactos <- list()
fishedIndex <- seq(2,6,by=2)
for(t in 1:ntts){
  thisRespVar <- paste("ts_",test_timeSteps[t],sep="")
  this_df <- getBestLModel(thisVar = thisRespVar, thisFactors = all_factors, cutoff = 0.1, 
                           data2fit = allFishedLogData)
  storeFishedModelDFs[[t]]<-this_df
  storeFishedFinalR[[t]]<-max(this_df$Rsq)
  storeNfactos[[t]] <- dim(this_df)[1]
  cat(t,", ")
}
plot(unlist(storeFishedFinalR))
plot(unlist(storeNfactos))

test<-unlist(storeFishedModelDFs)
table(names(test))
table(test[grep("Factor",names(test))])

plot(x=allFishedLogData[,"ts_100"], y=allFishedLogData[,"ts_100"], type="n", ylim=c(0, 1), xlim=c(0, 1), xlab="Observed", ylab="Fitted")
for(t in 36:116){
  thisRespVar<-paste("ts_",t,sep="")
  this_formula <- as.formula(paste(thisRespVar, "~", paste(storeFishedModelDFs[[100]]$Factor, collapse=" + "), sep=""))
  this_model <- glm(this_formula, data=allFishedLogData)
  testFit <- predict.glm(this_model, newdata = allFishedLogData)
  par(mfrow=c(2,1))
  hist_fit <- hist(testFit, plot=FALSE); hist_data <- hist(allFishedLogData[,c(thisRespVar)], plot=FALSE)
  par(mfrow=c(1,1))
  points(x=allFishedLogData[,c(thisRespVar)], y=testFit, pch=20, col=colByTime[t], lwd=2, lty=2)
}
thisLine <- seq(0, 1, length.out = 100)
points(x=thisLine, y=thisLine, type="l", col="red", lwd=2, lty=2)
mtext("Models with fishing", side=3, adj=0)

#  all models whether fished or not
storeAllModelDFs <- NULL; storeFinalR <- list(); storeNfactos <- list()
fishedIndex <- seq(2,6,by=2)
for(t in 1:ntts){
  thisRespVar <- paste("ts_",test_timeSteps[t],sep="")
  this_df <- getBestLModel(thisVar = thisRespVar, thisFactors = all_factors, cutoff = 0.1, data2fit = allLogData)
  storeAllModelDFs[[t]]<-this_df
  storeFinalR[[t]]<-max(this_df$Rsq)
  storeNfactos[[t]] <- dim(this_df)[1]
  cat(t,", ")
}

ALLtest<-c(unlist(storeAllModelDFs), unlist(storeUnfishedModelDFs), unlist(storeFishedModelDFs))
ALLtestFactors <- unique(sort(ALLtest[grep("^Factor", names(ALLtest))]))



plot(unlist(storeFinalR))
plot(unlist(storeNfactos))
# which year does it shift from 4 factors to 3?
min((1:length(storeNfactos))[as.double(unlist(storeNfactos))<4])+1899

test<-unlist(storeAllModelDFs)
table(names(test))
table(test[grep("Factor",names(test))])

## plot explanatory vars by time and power
factorColors <- colorRampPalette(colors=c(myBlue,  myAqua, myGreen, myGold,myOrange,"red", myRed,  myPurple, myGrey,"black"))(length(all_factors))
factorNumColors <-rev(colorRampPalette(colors = c(myPurple,"midnightblue",myBlue, myAqua, myLightBlue))(6))

## Age structured explanatory variables
ASfactors <- c("ChaosAlt", "TL",  "NumL1cons", "Lifespan",  "B0", "Informance",     "Keystone" ,  "propAdM"  ,     "propJuvM"  , "PropByTopPrey" ,      "Linf" ) 
allInts <- getIntVars(ASfactors); all_ASfactors <- c(ASfactors, allInts)

# if already fitted the models and going back to sort out colours, can use the explanatory variables actually used, otherwise start with this defin. 
if(is.null(AStestFactors)==TRUE){
    FACTORSUSED <-sort(unique( c("NumL1cons", "ChaosAlt:TL", "ChaosAlt:NumL1cons", "ChaosAlt:Lifespan", "ChaosAlt:Informance", "TL:NumL1cons", "TL:Informance", 
  "TL:propAdM", "NumL1cons:Linf", "Lifespan:Linf", "propAdM:PropByTopPrey", "propAdM:Linf", "propJuvM:Linf", "PropByTopPrey:Linf", 
  "ChaosAlt:Linf", "B0:Linf", "ChaosAlt:propAdM", "Keystone:propAdM", "NumL1cons:PropByTopPrey", "TL:Linf", "Informance:propAdM", "NumL1cons:B0", "ChaosAlt:B0", 
  "TL:B0", "Keystone:propAdM", "Keystone:Linf", "TL:PropByTopPrey"))); 
} else{
  FACTORSUSED <- sort(unique(c(ALLtestFactors, BPtestFactors, AStestFactors)))
}

nFU <- length(FACTORSUSED)
FACTCRSUSEDColours <- colorRampPalette(colors=c("midnightblue",myBlue,  myAqua, myGreen, myYellow,  myOrange,"red", "magenta", myPurple))(nFU)
FACTCRSUSEDColours <- colorRampPalette(colors=c("midnightblue",myBlue,  myAqua, myGreen, myYellow,  "red", myRed))(nFU)
plot(1:nFU, col=FACTCRSUSEDColours, pch=20, cex=2)

# ALLFACTORS <- all_ASfactors; nALLFactors <- length(ALLFACTORS)
# ALLFactorColours <- colorRampPalette(colors=c("midnightblue",myBlue, myAqua,myLightBlue,  myDarkGreen,myGreen, "wheat",myYellow,  myGold,myOrange,"red", myRed, "pink","magenta", myPurple, myGrey,"black","brown","wheat"))(nALLFactors)
# plot(1:nALLFactors, col=ALLFactorColours, pch=20)
plotFactorsByTime <- function(thisData2plot, modelLabel="", all_factors=all_factors){
  temp<- unlist(thisData2plot)
  factorIndex <- grep("Factor",names(temp)); R2Index <- grep("Rsq", names(temp)); this_nrows <- length(factorIndex)
  toPlot <- data.frame(matrix(NA, ncol=4, nrow=this_nrows))
  colnames(toPlot)<-c("Factor", "FactorOrder", "R2", "Timestep")
  toPlot$Factor <- temp[factorIndex]
  toPlot$R2 <- as.double(temp[R2Index])
  toPlot$FactorOrder<-names(temp)[factorIndex]
  toPlot$Timestep[1]<-1; toPlot$increasedR2 <- NA; toPlot$increasedR2[1]<-toPlot$R2[1]
  for(t in 2:this_nrows){
    prevForder <- as.double(gsub("Factor","",toPlot$FactorOrder[(t-1)]))
    thisForder <- as.double(gsub("Factor", "", toPlot$FactorOrder[(t)]))
    if(is.na(thisForder)){thisForder<-1}
    if(is.na(prevForder)){prevForder<-1}
    if(prevForder<thisForder){
      toPlot$Timestep[t]<-toPlot$Timestep[(t-1)]
      toPlot$increasedR2[t] <- toPlot$R2[t]-toPlot$R2[(t-1)]
    } else{
      toPlot$Timestep[t]<-toPlot$Timestep[(t-1)]+1
      toPlot$increasedR2[t] <- toPlot$R2[t]
    }
  }
  toPlot$year <- toPlot$Timestep + 1899
 
  XX <- table(toPlot$Factor)
  thisFactorsUsed <- sort(names(XX[XX>0]))
  toPlot$Factor <- factor(toPlot$Factor, levels=thisFactorsUsed)
  thisFactorColours <- FACTCRSUSEDColours[match(thisFactorsUsed, FACTORSUSED)]
  
  bp<-ggplot(data = toPlot[order(toPlot$FactorOrder),], aes(x = year, fill = Factor, y = increasedR2)) + 
    geom_bar(stat = 'identity')+ scale_fill_manual(values=thisFactorColours) + labs(y=expression(R^2), x="") + 
    theme(axis.text=element_text(size=12, angle=0),axis.title=element_text(size=12), plot.margin = margin(0.2,0.2,0,0.2, "cm"), 
          legend.text = element_text(size=8))  + guides(fill=guide_legend(title="")) + ggtitle(modelLabel) +guides(size = guide_legend(order=0.5))
  return(bp) 
} 

plotFactorNumberByTime <- function(thisData2plot, modelLabel="", all_factors=all_factors){
  temp<- unlist(thisData2plot)
  factorIndex <- grep("Factor",names(temp)); R2Index <- grep("Rsq", names(temp)); this_nrows <- length(factorIndex)
  toPlot <- data.frame(matrix(NA, ncol=4, nrow=this_nrows))
  colnames(toPlot)<-c("Factor", "FactorOrder", "R2", "Timestep")
  toPlot$Factor <- temp[factorIndex]
  toPlot$R2 <- as.double(temp[R2Index])
  toPlot$FactorOrder<-names(temp)[factorIndex]
  toPlot$Timestep[1]<-1; toPlot$increasedR2 <- NA; toPlot$increasedR2[1]<-toPlot$R2[1]
  for(t in 2:this_nrows){
    prevForder <- as.double(gsub("Factor","",toPlot$FactorOrder[(t-1)]))
    thisForder <- as.double(gsub("Factor", "", toPlot$FactorOrder[(t)]))
    if(is.na(prevForder)){prevForder<-1}
    if(is.na(thisForder)){thisForder<-1}
    if(prevForder<thisForder){
      toPlot$Timestep[t]<-toPlot$Timestep[(t-1)]
      toPlot$increasedR2[t] <- toPlot$R2[t]-toPlot$R2[(t-1)]
    } else{
      toPlot$Timestep[t]<-toPlot$Timestep[(t-1)]+1
      toPlot$increasedR2[t] <- toPlot$R2[t]
    }
  }
  toPlot$year <- toPlot$Timestep + 1899
  toPlot$FactorOrder <- factor(toPlot$FactorOrder, levels=paste("Factor",c("",1:10),sep=""))
  
  bp<-ggplot(data = toPlot, aes(x = year, fill = FactorOrder, y = increasedR2)) + 
    geom_bar(stat = 'identity')+ scale_fill_manual(values=factorNumColors) + labs(y=expression(R^2), x="") + 
    theme(axis.text=element_text(size=12, angle=0),axis.title=element_text(size=12), plot.margin = margin(0.2,1.5,0,0.2, "cm"), 
          legend.text = element_text(size=8)) + 
    guides(fill=guide_legend(title="")) + ggtitle(modelLabel) +guides(size = guide_legend(order=0.5))
  return(bp) 
} 

lab_cex=2

thisData2plot <- storeAllModelDFs
bpA <- plotFactorsByTime(thisData2plot, modelLabel="A: All model runs")
bpFnumsA <- plotFactorNumberByTime(thisData2plot, modelLabel="D: All model runs")

thisData2plot <- storeFishedModelDFs
bpB <- plotFactorsByTime(thisData2plot, modelLabel="B: Fished model runs")
bpFnumsB <- plotFactorNumberByTime(thisData2plot, modelLabel="E: Fished model runs")

thisData2plot <- storeUnfishedModelDFs
bpC <- plotFactorsByTime(thisData2plot, modelLabel="C: No-fishing model runs")
bpFnumsC <- plotFactorNumberByTime(thisData2plot, modelLabel="F: No-fishing model runs")

pdf(paste(plotPath, "FactorsByTimeAndR2.pdf", sep=""), height=10, width=8)
par(las=1)
plot_grid(bpA, bpB, bpC, nrow=3)
dev.off()


pdf(paste(plotPath, "FactorNumberByTimeAndR2.pdf", sep=""), height=10, width=8)
plot_grid(bpFnumsA, bpFnumsB, bpFnumsC, nrow=3)
dev.off()


pdf(paste(plotPath, "FactorAndNumberByTimeAndR2.pdf", sep=""), height=10, width=8)
plot_grid(bpA, bpB, bpC, bpFnumsA, bpFnumsB, bpFnumsC, nrow=2)
dev.off()

############
## repeat this with only biomass pool species (and same expl. vars), then only age-structured species (and full expl. vars)
## BIOMASS POOL GROUPS ##
BPfactors <- c("ChaosAlt", "TL",  "NumL1cons", "Lifespan",  "B0") # note, perturbMethod has chaosrun
allInts <- getIntVars(BPfactors); all_BPfactors <- c(BPfactors, allInts)

#split by fished or not fished
BPlogData <- allLogData[allLogData$Code %in% BPgroups, ]
BPlogData$ChaosAlt <- factor(BPlogData$ChaosAlt, levels=c("A", "B", "C", "D", "E", "F"))
storeBPUnfishedModelDFs <- NULL; storeBPUnfishedFinalR <- list(); storeBPUnfishedNfactors<-list()
allBPUnfishedLogData <- BPlogData[BPlogData$ChaosAlt %in% chaosLetters[unfishedIndex],]
for(t in 1:ntts){
  thisRespVar <- paste("ts_",test_timeSteps[t],sep="")
  this_df <- getBestLModel(thisVar = thisRespVar, thisFactors = all_BPfactors, cutoff = 0.1, 
                           data2fit =allBPUnfishedLogData)
  storeBPUnfishedModelDFs[[t]]<-this_df
  storeBPUnfishedFinalR[[t]]<-max(this_df$Rsq); storeBPUnfishedNfactors[[t]]<-dim(this_df)[1]
  cat(t,", ")
}

# fished ## its not actually logged anymore - it was; now it's cubed root because residuals look better
allBPFishedLogData <- BPlogData[BPlogData$ChaosAlt %in% chaosLetters[fishedIndex],]
storeBPFishedModelDFs <- NULL; storeBPFishedFinalR <- list(); storeBPFishedNfactos <- list()
fishedIndex <- seq(2,6,by=2)
for(t in 1:ntts){
  thisRespVar <- paste("ts_",test_timeSteps[t],sep="")
  this_df <- getBestLModel(thisVar = thisRespVar, thisFactors = all_BPfactors, cutoff = 0.1, 
                           data2fit = allBPFishedLogData)
  storeBPFishedModelDFs[[t]]<-this_df
  storeBPFishedFinalR[[t]]<-max(this_df$Rsq)
  storeBPFishedNfactos[[t]] <- dim(this_df)[1]
  cat(t,", ")
}

#  all models whether fished or not
storeAllBPModelDFs <- NULL; storeAllBPFinalR <- list(); storeAllBPNfactos <- list()
for(t in 1:ntts){
  thisRespVar <- paste("ts_",test_timeSteps[t],sep="")
  this_df <- getBestLModel(thisVar = thisRespVar, thisFactors = all_BPfactors, cutoff = 0.1, data2fit = BPlogData)
  storeAllBPModelDFs[[t]]<-this_df
  storeAllBPFinalR[[t]]<-max(this_df$Rsq)
  storeAllBPNfactos[[t]] <- dim(this_df)[1]
  cat(t,", ")
}

BPtest<-c(unlist(storeAllBPModelDFs), unlist(storeBPUnfishedModelDFs), unlist(storeBPFishedModelDFs))
BPtestFactors <- unique(sort(BPtest[grep("^Factor", names(BPtest))]))

## plot explanatory vars by time and power
thisData2plot <- storeAllBPModelDFs
bpA <- plotFactorsByTime(thisData2plot, modelLabel="A: All model runs")
bpFnumsA <- plotFactorNumberByTime(thisData2plot, modelLabel="D: All model runs")

thisData2plot <- storeBPFishedModelDFs
bpB <- plotFactorsByTime(thisData2plot, modelLabel="B: Fished model runs")
bpFnumsB <- plotFactorNumberByTime(thisData2plot, modelLabel="E: Fished model runs")

thisData2plot <- storeBPUnfishedModelDFs
bpC <- plotFactorsByTime(thisData2plot, modelLabel="C: No-fishing model runs")
bpFnumsC <- plotFactorNumberByTime(thisData2plot, modelLabel="F: No-fishing model runs")

pdf(paste(plotPath, "BPFactorsByTimeAndR2.pdf", sep=""), height=10, width=8)

plot_grid(bpA, bpB, bpC, nrow=3)
dev.off()


pdf(paste(plotPath, "BPFactorNumberByTimeAndR2.pdf", sep=""), height=10, width=8)
plot_grid(bpFnumsA, bpFnumsB, bpFnumsC, nrow=3)
dev.off()


pdf(paste(plotPath, "BPFactorAndNumberByTimeAndR2.pdf", sep=""), height=10, width=8)
plot_grid(bpA, bpB, bpC, bpFnumsA, bpFnumsB, bpFnumsC, nrow=2)
dev.off()

############
## AGE STRUCTURED GROUPS ##

#split by fished or not fished
ASlogData <- allLogData[!(allLogData$Code %in% BPgroups), ]
ASlogData$ChaosAlt <- factor(ASlogData$ChaosAlt, levels=c("A", "B", "C", "D", "E", "F"))
storeASUnfishedModelDFs <- NULL; storeASUnfishedFinalR <- list(); storeASUnfishedNfactors<-list()
ASUnfishedLogData <- ASlogData[ASlogData$ChaosAlt %in% chaosLetters[unfishedIndex],]
for(t in 1:ntts){
  thisRespVar <- paste("ts_",test_timeSteps[t],sep="")
  this_df <- getBestLModel(thisVar = thisRespVar, thisFactors = all_ASfactors, cutoff = 0.1, 
                           data2fit =ASUnfishedLogData)
  storeASUnfishedModelDFs[[t]]<-this_df
  storeASUnfishedFinalR[[t]]<-max(this_df$Rsq); storeASUnfishedNfactors[[t]]<-dim(this_df)[1]
  cat(t,", ")
}
  
# fished # the logged version
ASFishedData <- ASlogData[ASlogData$ChaosAlt %in% chaosLetters[fishedIndex],]
storeASFishedModelDFs <- NULL; storeASFishedFinalR <- list(); storeASNfactos <- list()
fishedIndex <- seq(2,6,by=2)
for(t in 1:ntts){
  thisRespVar <- paste("ts_",test_timeSteps[t],sep="")
  this_df <- getBestLModel(thisVar = thisRespVar, thisFactors = all_ASfactors, cutoff = 0.1, 
                           data2fit = ASFishedData)
  storeASFishedModelDFs[[t]]<-this_df
  storeASFishedFinalR[[t]]<-max(this_df$Rsq)
  storeASNfactos[[t]] <- dim(this_df)[1]
  cat(t,", ")
}

#  all models whether fished or not
storeASAllModelDFs <- NULL; storeASFinalR <- list(); storeASNfactos <- list()
for(t in 1:ntts){
  thisRespVar <- paste("ts_",test_timeSteps[t],sep="")
  this_df <- getBestLModel(thisVar = thisRespVar, thisFactors = all_ASfactors, cutoff = 0.1, data2fit = ASlogData)
  storeASAllModelDFs[[t]]<-this_df
  storeASFinalR[[t]]<-max(this_df$Rsq)
  storeASNfactos[[t]] <- dim(this_df)[1]
  cat(t,", ")
}

AStest<-c(unlist(storeASAllModelDFs), unlist(storeASUnfishedModelDFs), unlist(storeASFishedModelDFs))
AStestFactors <- unique(sort(AStest[grep("^Factor", names(AStest))]))
# factorColors <- colorRampPalette(colors=c(myBlue,  myAqua, myGreen,  myGold,myOrange,"red", myRed,  myPurple, myLightBlue, myBlue, myDarkGreen, myGrey,"black"))(length(all_ASfactors))
# factorNumColors <-rev(colorRampPalette(colors = c(myPurple,"midnightblue",myBlue, myAqua, myLightBlue))(10))

## plot explanatory vars by time and power
thisData2plot <- storeASAllModelDFs
bpA <- plotFactorsByTime(thisData2plot, modelLabel="A: All model runs", all_factors=all_ASfactors)
bpFnumsA <- plotFactorNumberByTime(thisData2plot, modelLabel="D: All model runs", all_factors=all_ASfactors)

thisData2plot <- storeASFishedModelDFs
bpB <- plotFactorsByTime(thisData2plot, modelLabel="B: Fished model runs", all_factors=all_ASfactors)
bpFnumsB <- plotFactorNumberByTime(thisData2plot, modelLabel="E: Fished model runs", all_factors=all_ASfactors)

thisData2plot <- storeASUnfishedModelDFs
bpC <- plotFactorsByTime(thisData2plot, modelLabel="C: No-fishing model runs", all_factors=all_ASfactors)
bpFnumsC <- plotFactorNumberByTime(thisData2plot, modelLabel="F: No-fishing model runs", all_factors=all_ASfactors)


pdf(paste(plotPath, "AS_bpA.pdf", sep=""), height=4, width=15)
bpA
dev.off()

pdf(paste(plotPath, "AS_bpB.pdf", sep=""), height=4, width=15)
bpB
dev.off()

pdf(paste(plotPath, "AS_bpC.pdf", sep=""), height=4, width=15)
bpC
dev.off()



pdf(paste(plotPath, "ASFactorsByTimeAndR2.pdf", sep=""), height=10, width=8)
plot_grid(bpA, bpB, bpC, nrow=3)
dev.off()


pdf(paste(plotPath, "ASFactorNumberByTimeAndR2.pdf", sep=""), height=4, width=15)
plot_grid(bpFnumsA, bpFnumsB, bpFnumsC, nrow=1)
dev.off()


pdf(paste(plotPath, "ASFactorAndNumberByTimeAndR2.pdf", sep=""), height=6, width=15)
plot_grid(bpA, bpB, bpC, bpFnumsA, bpFnumsB, bpFnumsC, nrow=2)
dev.off()

### set upi varDescr for plot labels
varFactors <- c("ChaosAlt", "TL", "NumL1cons",  "Lifespan",  "B0",   "PropByTopPrey")
varDesc <- data.frame(matrix(NA, nrow=length(varFactors), ncol=2))
colnames(varDesc)<- c("Variable", "Decr")
varDesc$Variable <- varFactors
thisDescriptions <- c("ChaosAlt", "Trophic level",  "Primary connections", "Lifespan (years)", "", "Diet proportion by top prey")
varDesc$Decr <- thisDescriptions

lab_cex=1.2

#########################################################################################################
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
#####################################################################################################################
#
# all data together, taking from 1925 as all versions (ALL, BP, AS) most consistent after this timestep
#########################################################################################################
#
testMelt <- melt(allLogData, id.vars=colnames(allLogData)[1:13])
## fit to all CV over time 1900:2015, or 1950:2015 - this is timestep 61, as starts in 1865
testMelt$timestep <- as.double(unlist(lapply(testMelt$variable,FUN=function(x){get_first_number(x)})))
testMelt$ChaosAlt <- factor(testMelt$ChaosAlt, levels=c("A", "B", "C", "D", "E", "F"))
this_startTimestep <- 61
this_endTimestep <- nts

thisRespVar <- "value"

thisMin <- 0 
thisIndex <- testMelt$timestep >= this_startTimestep & testMelt$timestep <= this_endTimestep
ALLmeltData <- testMelt[thisIndex,]
ALLmelt_df <- getBestLModel(thisVar = thisRespVar, thisFactors = all_factors, cutoff = 0.1, data2fit = ALLmeltData) ## reduce r^2 cutoff to check if just missing any expl. vars

ALLmeltFORMULA <- as.formula(paste(thisRespVar, "~", paste(ALLmelt_df$Factor, collapse=" + "), sep=""))
ALLmeltModel <- glm(ALLmeltFORMULA, data=ALLmeltData)
fittedALL <- ((predict.glm(ALLmeltModel, newdata =ALLmeltData)))
thisIntVar <- "ChaosAlt:TL"; this_model=ALLmeltModel; thisData = ALLmeltData
plotIntEffect_noChaosRun(thisIntVar, this_model=ALLmeltModel, thisMin=0, inputMax=1, updateMax=TRUE, thisData=ALLmeltData)
 
ALL_pearsonsResids <- residuals.glm(ALLmeltModel, type="pearson")
ALL_resids <- residuals.glm(ALLmeltModel, type="deviance")

pdf(paste(plotPath,"residualsALL.pdf", sep=""), height=6, width=10)
par(mfrow=c(2,2), mar=c(4.5,5,1.2,1))
plot(x=fittedALL, y=ALL_pearsonsResids, pch=20, col=myBlue_trans, cex.lab=1.5, ylim=c(-0.5,0.5), ylab="Pearson's residuals", xlab="Fitted") 
abline(h=0, col="red", lwd=1.5, lty=2)
mtext("A", side=3, adj=0, font=2)
# by trophic level
plot(x=ALLmeltData$TL, y=ALL_pearsonsResids, pch=20, col=myBlue_trans, ylab="Pearson's residuals", xlab="Trophic level", cex.lab=1.5)
mtext("B", side=3, adj=0, font=2)
# by chaosAlt
plot(x=1:nchaosVersions, y=rep(0, nchaosVersions), xaxt="n", xlab="ChaosAlt", ylab="Pearson's residuals", cex.lab=1.5, ylim=c(-0.5, 0.5), type="n", xlim=c(0.5, (nchaosVersions+0.5)))
axis(at=1:nchaosVersions, labels=chaosLetters, side=1)
for(c in 1:nchaosVersions){
  thisCindex <- ALLmeltData$ChaosAlt==chaosLetters[c]
  boxplot(ALL_pearsonsResids[thisCindex], at=c, col=myBlue_trans, border=myBlue, add=TRUE, pch=20, yaxt="n")
}
mtext("C", side=3, adj=0, font=2)
dev.off()



#########################################################################################################
  
#BP only
thisBPindex <- testMelt$Code %in% BPgroups & testMelt$timestep >= this_startTimestep
BPmeltData <- testMelt[thisBPindex,]
this_BPdf <- getBestLModel(thisVar = thisRespVar, thisFactors =all_BPfactors, cutoff = 0.1, data2fit = BPmeltData)

BPmeltFORMULA <- as.formula(paste(thisRespVar, "~", paste(this_BPdf$Factor, collapse=" + "), sep=""))
BPmeltModel <- glm(BPmeltFORMULA, data=BPmeltData)
fittedBP <- ((predict.glm(BPmeltModel, newdata =BPmeltData)))
thisIntVar <- "ChaosAlt:TL"; this_model=BPmeltModel; thisData = BPmeltData
plotIntEffect_noChaosRun(thisIntVar, this_model=BPmeltModel, thisMin=0, inputMax=1, updateMax=TRUE, thisData=BPmeltData)

BP_pearsonsResids <- residuals.glm(BPmeltModel, type="pearson")


pdf(paste(plotPath,"residualsBPall.pdf", sep=""), height=9, width=10)
par(mfrow=c(3,2), mar=c(4.5,5,1.2,1))
plot(x=fittedBP, y=BP_pearsonsResids, pch=20, col=myBlue_trans, cex.lab=1.5, ylim=c(-0.5,0.5), ylab="Pearson's residuals", xlab="Fitted") 
abline(h=0, col="red", lwd=1.5, lty=2)
mtext("A", side=3, adj=0, font=2)
# by trophic level
plot(x=BPmeltData$TL, y=BP_pearsonsResids, pch=20, col=myBlue_trans, ylab="Pearson's residuals", xlab="Trophic level", cex.lab=1.5)
mtext("B", side=3, adj=0, font=2)
# by chaosAlt
plot(x=1:nchaosVersions, y=rep(0, nchaosVersions), xaxt="n", xlab="ChaosAlt", ylab="Pearson's residuals", cex.lab=1.5, ylim=c(-0.5, 0.5), type="n", xlim=c(0.5, (nchaosVersions+0.5)))
axis(at=1:nchaosVersions, labels=chaosLetters, side=1)
for(c in 1:nchaosVersions){
  thisCindex <- BPmeltData$ChaosAlt==chaosLetters[c]
  boxplot(BP_pearsonsResids[thisCindex], at=c, col=myBlue_trans, border=myBlue, add=TRUE, pch=20, yaxt="n")
}
mtext("C", side=3, adj=0, font=2)
# by NumL1Cons
plot(x=BPmeltData$NumL1cons, y=BP_pearsonsResids, pch=20, col=myBlue_trans, ylab="Pearson's residuals", xlab="Number of primary trophic connections", cex.lab=1.5)
mtext("D", side=3, adj=0, font=2)
# by B0
plot(x=BPmeltData$B0, y=BP_pearsonsResids, pch=20, col=myBlue_trans, ylab="Pearson's residuals", xlab=expression(B[0]), cex.lab=1.5)
mtext("E", side=3, adj=0, font=2)
dev.off()
##############################################################################
# AS only
thisASindex <- !(testMelt$Code %in% BPgroups) & testMelt$timestep >= this_startTimestep & !is.na(testMelt$Linf)
## Linf doesn't get selected, so leave NAs in there as otherwise it knocks out Informance=16 (the most uninformed ones)
## take it back - if have NA's it ends up being selected, and messes up results
## go back to taking NA Linfs out
# thisASindex <- !(testMelt$Code %in% BPgroups) & testMelt$timestep >= this_startTimestep 
ASmeltData <- testMelt[thisASindex,]
this_ASdf <- getBestLModel(thisVar = thisRespVar, thisFactors =all_ASfactors, cutoff = 0.1, data2fit = ASmeltData)

ASmeltFORMULA <- as.formula(paste(thisRespVar, "~", paste(this_ASdf$Factor, collapse=" + "), sep=""))
ASmeltModel <- glm(ASmeltFORMULA, data=ASmeltData)
fittedAS <- (predict.glm(ASmeltModel, newdata =ASmeltData))
thisIntVar <- "ChaosAlt:TL"; this_model=ASmeltModel; thisData = ASmeltData
plotIntEffect_noChaosRun(thisIntVar, this_model=ASmeltModel, thisMin=0, inputMax=1, updateMax=TRUE, thisData=ASmeltData)

AS_pearsonsResids <- residuals.glm(ASmeltModel, type="pearson")

pdf(paste(plotPath,"residualsASall.pdf", sep=""), height=6, width=10)
par(mfrow=c(2,2), mar=c(4.5,5,1.2,1))
plot(x=fittedAS, y=AS_pearsonsResids, pch=20, col=myBlue_trans, cex.lab=1.5, ylim=c(-0.5,0.5), ylab="Pearson's residuals", xlab="Fitted") 
abline(h=0, col="red", lwd=1.5, lty=2)
mtext("A", side=3, adj=0, font=2)
# by trophic level
plot(x=ASmeltData$TL, y=AS_pearsonsResids, pch=20, col=myBlue_trans, ylab="Pearson's residuals", xlab="Trophic level", cex.lab=1.5)
mtext("B", side=3, adj=0, font=2)
# by chaosAlt
plot(x=1:nchaosVersions, y=rep(0, nchaosVersions), xaxt="n", xlab="ChaosAlt", ylab="Pearson's residuals", cex.lab=1.5, ylim=c(-0.5, 0.5), type="n", xlim=c(0.5, (nchaosVersions+0.5)))
axis(at=1:nchaosVersions, labels=chaosLetters, side=1)
for(c in 1:nchaosVersions){
  thisCindex <- ASmeltData$ChaosAlt==chaosLetters[c]
  boxplot(AS_pearsonsResids[thisCindex], at=c, col=myBlue_trans, border=myBlue, add=TRUE, pch=20, yaxt="n")
}
mtext("C", side=3, adj=0, font=2)
# by NumL1cons
plot(x=ASmeltData$NumL1cons, y=AS_pearsonsResids, pch=20, col=myBlue_trans, ylab="Pearson's residuals", xlab="Number of primary trophic connections", cex.lab=1.5)
mtext("D", side=3, adj=0, font=2)
dev.off()

###################
## do fished and unfished version of AS
## FISHED AS MELTED ##
ASFishedMeltData <- ASmeltData[ASmeltData$ChaosAlt %in% chaosLetters[fishedIndex],]
this_ASFisheddf <- getBestLModel(thisVar = thisRespVar, thisFactors =all_ASfactors, cutoff = 0.1, data2fit = ASFishedMeltData)

ASFishedmeltFORMULA <- as.formula(paste(thisRespVar, "~", paste(this_ASFisheddf$Factor, collapse=" + "), sep=""))
ASFishedmeltModel <- glm(ASFishedmeltFORMULA, data=ASFishedMeltData)
fittedASFished <- (predict.glm(ASFishedmeltModel, newdata =ASFishedMeltData))
thisIntVar <- "ChaosAlt:TL"; this_model=ASFishedmeltModel; thisData = ASFishedMeltData
plotIntEffect_noChaosRun(thisIntVar, this_model=ASFishedmeltModel, thisMin=0, inputMax=1, updateMax=TRUE, thisData=ASFishedMeltData)

ASfished_pearsonsResids <- residuals.glm(ASFishedmeltModel, type="pearson")

pdf(paste(plotPath,"residualsASFished.pdf", sep=""), height=6, width=10)
par(mfrow=c(2,2), mar=c(4.5,5,1.5,1))
plot(x=fittedASFished, y=ASfished_pearsonsResids, pch=20, col=myBlue_trans, cex.lab=1.5, ylim=c(-0.5,0.5), ylab="Pearson's residuals", xlab="Fitted") 
abline(h=0, col="red", lwd=1.5, lty=2)
mtext("A", side=3, adj=0, font=2)
# by trophic level
plot(x=ASFishedMeltData$TL, y=ASfished_pearsonsResids, pch=20, col=myBlue_trans, ylab="Pearson's residuals", xlab="Trophic level", cex.lab=1.5)
mtext("B", side=3, adj=0, font=2)
# by chaosAlt
plot(x=(1:nchaosVersions)[fishedIndex], y=rep(0, nchaosVersions)[fishedIndex], xaxt="n", xlab="ChaosAlt", ylab="Pearson's residuals", cex.lab=1.5, ylim=c(-0.5, 0.5), type="n", xlim=c(0.5, (nchaosVersions+0.5)))
axis(at=(1:nchaosVersions)[fishedIndex], labels=chaosLetters[fishedIndex], side=1)
for(c in 1:nchaosVersions){
  thisC <- chaosLetters[fishedIndex][c]
  if(!is.na(thisC)){
    thisCindex <- ASFishedMeltData$ChaosAlt==thisC
    boxplot(ASfished_pearsonsResids[thisCindex], at=grep(thisC, chaosLetters), col=myBlue_trans, border=myBlue, add=TRUE, pch=20, yaxt="n")
  }
}
mtext("C", side=3, adj=0, font=2)
# by NumL1cons
plot(x=ASFishedMeltData$NumL1cons, y=ASfished_pearsonsResids, pch=20, col=myBlue_trans, ylab="Pearson's residuals", xlab="Number of level 1 trophic connections", cex.lab=1.5)
mtext("D", side=3, adj=0, font=2)
dev.off()


## UNFISHED AS MELTED ##
ASUnfishedMeltData <- ASmeltData[!(ASmeltData$ChaosAlt %in% chaosLetters[fishedIndex]),]
this_ASUnfisheddf <- getBestLModel(thisVar = thisRespVar, thisFactors =all_ASfactors, cutoff = 0.1, data2fit = ASUnfishedMeltData)

ASUnfishedmeltFORMULA <- as.formula(paste(thisRespVar, "~", paste(this_ASUnfisheddf$Factor, collapse=" + "), sep=""))
ASUnfishedmeltModel <- glm(ASUnfishedmeltFORMULA, data=ASUnfishedMeltData)
fittedASUnfished <- (predict.glm(ASUnfishedmeltModel, newdata =ASUnfishedMeltData))
thisIntVar <- "ChaosAlt:TL"; this_model=ASUnfishedmeltModel; thisData = ASUnfishedMeltData
plotIntEffect_noChaosRun(thisIntVar, this_model=ASFishedmeltModel, thisMin=0, inputMax=1, updateMax=TRUE, thisData=ASFishedMeltData)

ASUnfished_pearsonsResids <- residuals.glm(ASUnfishedmeltModel, type="pearson")

pdf(paste(plotPath,"residualsASUnfished.pdf", sep=""), height=9, width=10)
par(mfrow=c(3,2), mar=c(4.5,5,1.5,1))
plot(x=fittedASUnfished, y=ASUnfished_pearsonsResids, pch=20, col=myBlue_trans, cex.lab=1.5, ylim=c(-0.5,0.5), ylab="Pearson's residuals", xlab="Fitted") 
abline(h=0, col="red", lwd=1.5, lty=2)
mtext("A", side=3, adj=0, font=2)
# by trophic level
plot(x=ASUnfishedMeltData$TL, y=ASUnfished_pearsonsResids, pch=20, col=myBlue_trans, ylab="Pearson's residuals", xlab="Trophic level", cex.lab=1.5)
mtext("B", side=3, adj=0, font=2)
# by chaosAlt
plot(x=(1:nchaosVersions)[unfishedIndex], y=rep(0, nchaosVersions)[unfishedIndex], xaxt="n", xlab="ChaosAlt", ylab="Pearson's residuals", cex.lab=1.5, ylim=c(-0.5, 0.5), type="n", xlim=c(0.5, (nchaosVersions+0.5)))
axis(at=(1:nchaosVersions)[unfishedIndex], labels=chaosLetters[unfishedIndex], side=1)
for(c in 1:nchaosVersions){
  thisC <- chaosLetters[unfishedIndex][c]
  if(!is.na(thisC)){
    thisCindex <- ASUnfishedMeltData$ChaosAlt==thisC
    boxplot(ASUnfished_pearsonsResids[thisCindex], at=grep(thisC, chaosLetters), col=myBlue_trans, border=myBlue, add=TRUE, pch=20, yaxt="n")
  }
}
mtext("C", side=3, adj=0, font=2)
# by NumL1cons
plot(x=ASUnfishedMeltData$NumL1cons, y=ASUnfished_pearsonsResids, pch=20, col=myBlue_trans, ylab="Pearson's residuals", xlab="Number of primary trophic connections", cex.lab=1.5)
mtext("D", side=3, adj=0, font=2)
# by Informance
temp<- table(ASUnfishedMeltData$Informance)
plot(x=names(temp), y=temp, type="n", col=myBlue_trans, ylab="Pearson's residuals", xlab="Informance", cex.lab=1.5, ylim=c(-0.5,0.5),
     xlim=c(as.double(min(names(temp)))-0.5, as.double(max(names(temp)))+0.5), xaxt="n")
axis(at= names(temp) , labels=c(2,3,4), side=1)
for(i in 1:length(temp)){
  thisIndex <- ASUnfishedMeltData$Informance==names(temp)[i]
  boxplot(ASUnfished_pearsonsResids[thisIndex], at=as.double(names(temp)[i]), add=TRUE, col=myBlue_trans, border=myBlue, pch=20, yaxt="n")
}
mtext("E", side=3, adj=0, font=2)
dev.off()

########################################################################################################################
## plot effects
########################################################################################################################

## TL:ChaosALt
pdf(paste(plotPath,"GLMeffectsAllBpAS.pdf", sep=""), height=12, width=8)
par(mfrow=c(3,1), mar=c(4.5,6,1.5,1))
plotIntEffect_noChaosRun(thisIntVar, this_model=ALLmeltModel, thisMin=0, inputMax=1, updateMax=TRUE, thisData=ALLmeltData)
mtext("A: ALL species groups", side=3, adj=0, cex=1.2)
plotIntEffect_noChaosRun(thisIntVar, this_model=BPmeltModel, thisMin=0, inputMax=1, updateMax=TRUE, thisData=BPmeltData)
mtext("B: BP species groups", side=3, adj=0, cex=1.2)
plotIntEffect_noChaosRun(thisIntVar, this_model=ASmeltModel, thisMin=0, inputMax=1, updateMax=TRUE, thisData=ASmeltData)
mtext("C: AS species groups", side=3, adj=0, cex=1.2)
dev.off()

# also do NumL1cons:B0 for BP, and chaosAlt:NumL1cons for AS model
pdf(paste(plotPath,"GLMeffectsBpAS.pdf", sep=""), height=8, width=8)
par(mfrow=c(2,1), mar=c(4.5,5.5,1.5,1))
plotIntEffect_noChaosRun("NumL1cons:B0", this_model=BPmeltModel, thisMin=0, inputMax=1, updateMax=TRUE, thisData=BPmeltData)
mtext("A: BP species groups", side=3, adj=0, cex=1.2)
plotIntEffect_noChaosRun("ChaosAlt:NumL1cons", this_model=ASmeltModel, thisMin=0, inputMax=1, updateMax=TRUE, thisData=ASmeltData)
mtext("B: AS species groups", side=3, adj=0, cex=1.2)
dev.off()

########################################################################################################################

thisChaosRuns<-c("A", "B","C", "D", "E", "F")
# do the end-of-run version
x <- grep(":", this_ASdf$Factor); intVars <- this_ASdf$Factor[x]; nivars <- length(intVars)
plotLetters<-c("A", "B", "C")
# pdf(paste(plotPath,"InteractionEffects_ALLdata150years.pdf", sep=""), height=4, width=12)
par(mfrow=c(1,2), mar=c(4,5.5,1.5,1), oma=c(0,0,0,10))
plotSingleEffect(thisVar="TL", this_model=finalBPAllModel)
mtext(plotLetters[1], side=3, adj=0, line=0.1, font=2)
for(i in 1:nivars){
  thisIntVar <- intVars[i]; 
  plotIntEffect_noChaosRun(thisIntVar = intVars[i], this_model=ASModel, inputMax=0.3, updateMax = FALSE, thisData = testMelt[thisASindex,])
  mtext(plotLetters[(i+1)], side=3, adj=0, line=0.1, font=2)
}
predLegendValues <- pretty(seq(0,0.3, length.out=5))
legendValues <- predLegendValues # absoulte, and percentage
legendCols <- unlist(lapply(predLegendValues, getColor, log=FALSE))
# makeBlankPlot()
par(xpd=NA)
legend(legend=legendValues, col=legendCols, x="right", bty="n", pch=15, pt.cex=1.5, title="Absolute relative\ndifference (%)", inset=-0.5)

# dev.off()




plot(x=allLogData[,"ts_100"], y=allLogData[,"ts_100"], type="n", ylim=c(0, 1), xlim=c(0, 1), xlab="Observed", ylab="Fitted")
for(t in 36:116){
  thisRespVar<-paste("ts_",t,sep="")
  this_formula <- as.formula(paste(thisRespVar, "~", paste(storeModelDFs[[100]]$Factor, collapse=" + "), sep=""))
  this_model <- glm(this_formula, data=allLogData)
  testFit <- predict.glm(this_model, newdata = allLogData)
  par(mfrow=c(2,1))
  hist_fit <- hist(testFit, plot=FALSE); hist_data <- hist(allLogData[,c(thisRespVar)], plot=FALSE)
  par(mfrow=c(1,1))
  points(x=allLogData[,c(thisRespVar)], y=testFit, pch=20, col=colByTime[t], lwd=2, lty=2)
}
thisLine <- seq(0, 1, length.out = 100)
points(x=thisLine, y=thisLine, type="l", col="red", lwd=2, lty=2)
mtext("All model runs", side=3, adj=0)



thisRespVar<-"value"
this_formula <- as.formula(paste(thisRespVar, "~", paste(this_df$Factor, collapse=" + "), sep=""))
this_model <- glm(this_formula, data=testMelt[thisIndex,])
testFit <- predict.glm(this_model, newdata = testMelt[thisIndex,])
par(mfrow=c(2,1))
hist_fit <- hist(testFit, plot=FALSE); hist_data <- hist(testMelt[thisIndex,c(thisRespVar)], plot=FALSE)
thisMin <- min(c(testFit, testMelt[thisIndex,c(thisRespVar)]), na.rm=TRUE)
thisMax <- max(c(testFit, testMelt[thisIndex,c(thisRespVar)]), na.rm=TRUE)
thisLine <- seq(thisMin, thisMax, length.out = 100)
par(mfrow=c(1,1))
plot(x=testMelt[thisIndex,c(thisRespVar)], y=testFit, pch=20, ylim=c(thisMin, thisMax), xlim=c(thisMin, thisMax))
points(x=thisLine, y=thisLine, type="l", col="red", lwd=2, lty=2)


