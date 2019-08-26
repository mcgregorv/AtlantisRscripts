this_run<-"TBGB_JP"
this_out<-"FishON"

this_path<-paste(DIR$'Base',"TBGB\\",this_run,"\\",sep="")

#read in catches, as fomated and such in summaryCatchByFleet_TBGB.R
load(paste(this_path,"Catch_history\\storeCatchByFleet", sep=""))
nfg <- dim(fishedGroups)[1]

# plot them to check look right
pdf(paste(plotPath, "TBGB_CatchHistoriesAndF.pdf", sep=""), height=10, width=10)
par(mfrow=c(5,2), mar=c(4,4.5,1.5,8), lend=1,  par(las=0))
for(g in 1:nfg){
  thisName <- gsub("_", " ",fishedGroups$Name[g])
  thisCatchHistory <- storeCatchByStock[g,]; thisF <- TBGB_F[g,]; thisSSB <- SSBbyfishedCode[g,]
  thisFaxis <- pretty(seq(0,max(thisF), length.out=10)); thisCatchAxis <- pretty(0:max(thisCatchHistory))
  
  plot(x=catchYears, y=thisCatchHistory, type="h", col=myGrey_trans, lwd=3, yaxt="n", xaxt="n", ylab="", xlab="")
  axis(at=thisCatchAxis, labels=thisCatchAxis, side=4, line=1, col.axis=myGrey, col.ticks=myGrey, col=myGrey)
  par(new=TRUE)
  plot(x=catchYears, y=thisSSB,type="l", lty=2, lwd=1.5, ylab="Biomass (tonnes)", xlab="", ylim=c(0, max(thisSSB, na.rm=TRUE)), cex.axis=thisCex, cex.lab=thisCex)
  par(las=0)
  mtext("Catch (tonnes)", side=4, line=3, col=myGrey, font=2)
  
  
  par(new=TRUE)
  plot(x=catchYears, y=thisF, type="l", col=myOrange, lwd=2, yaxt="n", xaxt="n", ylab="", xlab="")
  axis(at=thisFaxis, labels=thisFaxis, side=4, line=5, col.axis=myOrange, col.ticks=myOrange, col=myOrange)
  mtext("F", side=4, line=7, col=myOrange, font=2)
  
  mtext(thisName, side=3, adj=0, font=2)
}
dev.off()



