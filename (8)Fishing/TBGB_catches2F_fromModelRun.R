this_out<-"PBnew3"; runFolder="TBGB_JP2"
thisDesc <- paste(runFolder, this_out,sep="")

this_path = paste(DIR$'Base', "TBGB\\",runFolder,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")

plotPath<-paste(DIR$'Base',"TBGB\\Figures\\Testing\\FIshing\\",sep="")


#read in catches, as fomated and such in summaryCatchByFleet_TBGB.R
load(paste(this_path,"Catch_history\\storeCatchByFleet", sep=""))
nfg <- dim(fishedGroups)[1]

ThisNC.nc<-nc_open(paste(outPath,thisRun,".nc",sep=""))
thisVol <- ncvar_get(ThisNC.nc, "volume"); nts <- dim(thisVol)[3]
daysTimeStep<-73
numStepsPerYear<-365/daysTimeStep
year0<-1880
fishingStartYear<-1900
modelStartYear<-1900

xLabsTemp<-seq(0,(nts*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0
modelTime <- year0 + (1:nts)/numStepsPerYear


# plot them to check look right
pdf(paste(plotPath, "TBGB_CatchHistoriesAndF_fromModel",thisDesc,".pdf", sep=""), height=10, width=10)
par(mfrow=c(5,2), mar=c(4,4.5,1.5,8), lend=1,  par(las=0))
for(g in 1:nfg){
  thisTracerName <- str_trim(fishedGroups$Name[g], side="both"); thisName <- gsub("_", " ", thisTracerName)
  thisCatchHistory <- storeCatchByStock[g,]; thisF <- TBGB_F[g,]; thisSSB <- SSBbyfishedCode[g,]
  thisFaxis <- pretty(seq(0,max(thisF), length.out=10)); thisCatchAxis <- pretty(0:max(thisCatchHistory))
  
  # get model SSB (roughly - inclu juv for now)
  thisTracer <- paste(thisTracerName, "_N", sep=""); tempData <- ncvar_get(ThisNC.nc, thisTracer)
  modelSSB <- apply(tempData * thisVol, 3, sum) * mg_2_tonne * X_CN
  
  plot(x=catchYears, y=thisCatchHistory, type="h", col=myGrey_trans, lwd=3, yaxt="n", xaxt="n", ylab="", xlab="", xlim=c(modelStartYear, max(catchYears)))
  axis(at=thisCatchAxis, labels=thisCatchAxis, side=4, line=1, col.axis=myGrey, col.ticks=myGrey, col=myGrey)
  par(new=TRUE)
  thisBiomassMax <- max(c(max(thisSSB, na.rm=TRUE), max(modelSSB[modelTime>=modelStartYear], na.rm=TRUE), na.rm=TRUE))
  plot(x=catchYears, y=thisSSB,type="l", lty=2, lwd=1.5, ylab="Biomass (tonnes)", xlab="", ylim=c(0, thisBiomassMax*1.2), xlim=c(modelStartYear, max(catchYears)), cex.axis=thisCex, cex.lab=thisCex)
  points(x=modelTime, y=modelSSB, type="l", col=myBlue, lwd=1, lty=1)
  par(las=0)
  mtext("Catch (tonnes)", side=4, line=3, col=myGrey, font=2)
  
  
  par(new=TRUE)
  plot(x=catchYears, y=thisF, type="l", col=myOrange, lwd=2, yaxt="n", xaxt="n", ylab="", xlab="")
  axis(at=thisFaxis, labels=thisFaxis, side=4, line=5, col.axis=myOrange, col.ticks=myOrange, col=myOrange)
  mtext("F", side=4, line=7, col=myOrange, font=2)
  
  mtext(thisName, side=3, adj=0, font=2)
}
dev.off()



