# the analyses were carried out in TCorr_mL.R - this just reads in the values required to create the plots, and creates the plots
## original plots were created in TCorr_mL.R, but formatting and layout have changed here
plotPath <- paste(DIR$'Base',"Reports\\Chaos\\Figures\\", sep="")
inputsPath <- paste(DIR$'Base', "ATLANTISmodels\\inputs\\", sep="")

load(paste(inputsPath, "..\\storeChaosTemp\\storeValuesToPlot",sep=""))

# some other things to read in
groupsDF <- read.csv(paste(inputsPath,"..\\CRAM_Groups.csv", sep="")); ng <- dim(groupsDF)[1]
plotCodes <- groupsDF$Code[groupsDF$NumCohorts>1]; nplot<-length(plotCodes)
groupBioPars<-read.csv(paste(inputsPath,"..\\CRAM_B0.csv", sep=""))
groupsDFPaper<-read.csv(paste(inputsPath,"..\\CRAM_groupsPaper.csv", sep=""))
groupsDFOrderPaper <- groupsDFPaper[order(groupsDFPaper$Name),]
paperIndex <- groupsDFOrderPaper$Code %in% plotCodes

makeBoxPlotFromStats <- function(stats,at,col,border, width=1){
  x0<-at-0.5*width;  x1<-at+0.5*width
  y0<-stats[2]; y1<-stats[4]
  polygon(x=c(x0,x0,x1,x1), y=c(y0,y1,y1,y0), col=col, border=border)
  segments(x0=at, x1=at, y0=stats[1], y1=stats[2], col=border)
  segments(x0=at, x1=at, y0=stats[4], y1=stats[5], col=border)
  segments(x0=x0, x1=x1, y0=stats[3], y1=stats[3], col=border, lwd=2.5)
  x2 <- at-0.25*width;  x3<-at+0.25*width
  segments(x0=x2, x1=x3, y0=stats[5], y1=stats[5], col=border)
  segments(x0=x2, x1=x3, y0=stats[1], y1=stats[1], col=border)
}

# set up plotting region


par(mar=c(10,4,1,1))
plot(x=1:nplot, y=rep(1, nplot), type="n", ylim=c(0,1.5), ylab="Applied M / Model M", xlab="", xaxt="n")
abline(v=seq(0,nplot), col=myGrey_trans)
abline(h=seq(0,2,by=0.125), col=myGrey_trans)
abline(h=1, col="red", lty=2, lwd=1.5)
par(las=2)
axis(labels=gsub("_", " ",groupsDFOrderPaper$Name[paperIndex]), at=1:nplot, side=1)
for(g in 1:nplot){
  thisCode<-as.character(groupsDFOrderPaper$Code[paperIndex][g])
  thisGindex <- grep(thisCode, groupsDF$Code)
  thisModelMjuv <-  storeJuvMs[grep(thisCode, groupsDF$Code)]; thisModelMad <-  storeAdMs[grep(thisCode, groupsDF$Code)]      
  thisAdstats <- storeMboxplots_ad[[thisCode]]/thisModelMad; thisJuvstats <- storeMboxplots_juv[[thisCode]]/thisModelMjuv
  if(sum(storeMboxplots_ad[[thisCode]])==0){thisAdstats<-rep(0,5)}
  if(sum(storeMboxplots_juv[[thisCode]])==0){thisJuvstats<-rep(0,5)}
  makeBoxPlotFromStats(stats=thisJuvstats, at=g, col=myOrange_trans, border=myOrange, width=0.5)
  makeBoxPlotFromStats(stats=thisAdstats, at=g, col=myBlue_trans, border=myBlue, width=0.5)
}  

# do hist of median proportions for adults, and for females
hist(medianProps_ad, col=myBlue_trans, border=myBlue, xlab="Applied M / Model M", main="")
hist(medianProps_juv, col=myOrange_trans,  border=myOrange, xlab="Applied M / Model M", main="")

############ on one plot
pdf(paste(plotPath,"AppliedMoverModelM.pdf", sep=""), height=5, width=9)
par(mar=c(9,4,1,1), fig=c(0, (2/3), 0, 1), las=1)
plot(x=1:nplot, y=rep(1, nplot), type="n", ylim=c(0,1.5), ylab="Applied M / Model M", xlab="", xaxt="n")
abline(v=seq(0,nplot), col=myGrey_trans)
abline(h=seq(0,2,by=0.125), col=myGrey_trans)
abline(h=1, col="red", lty=2, lwd=1.5)
par(las=2)
axis(labels=gsub("_", " ",str_trim(groupsDFOrderPaper$Name[paperIndex], side="both")), at=1:nplot, side=1)
for(g in 1:nplot){
  thisCode<-as.character(groupsDFOrderPaper$Code[paperIndex][g])
  thisGindex <- grep(thisCode, groupsDF$Code)
  thisModelMjuv <-  storeJuvMs[grep(thisCode, groupsDF$Code)]; thisModelMad <-  storeAdMs[grep(thisCode, groupsDF$Code)]      
  thisAdstats <- storeMboxplots_ad[[thisCode]]/thisModelMad; thisJuvstats <- storeMboxplots_juv[[thisCode]]/thisModelMjuv
  if(sum(storeMboxplots_ad[[thisCode]])==0){thisAdstats<-rep(0,5)}
  if(sum(storeMboxplots_juv[[thisCode]])==0){thisJuvstats<-rep(0,5)}
  makeBoxPlotFromStats(stats=thisJuvstats, at=g, col=myOrange_trans, border=myOrange, width=0.5)
  makeBoxPlotFromStats(stats=thisAdstats, at=g, col=myBlue_trans, border=myBlue, width=0.5)
}  
# do hist of median proportions for adults, and for females
par(new=TRUE, fig=c((2/3), 1, (0.5), (1)), mar=c(5,4,2,0.5))
hist(medianProps_ad, col=myBlue_trans, border=myBlue, xlab="Median applied M / Model M", main="Adult")
par(new=TRUE, fig=c((2/3), 1, 0, (0.5)), mar=c(5,4,2,0.5))
hist(medianProps_juv, col=myOrange_trans,  border=myOrange, xlab="Median applied M / Model M", main="Juvenile")
dev.off()
############

# output median proportions by adult and juvenile
mLProps_df <- data.frame(matrix(NA, ncol=3, nrow=ng)); colnames(mLProps_df)<- c("Code", "propAd", "propJuv")
mLProps_df$Code <- as.character(groupsDF$Code)
for(g in 1:ng){
  thisCode <- mLProps_df$Code[g]
  thisNumCohorts <- groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    thisModelMjuv <-  storeJuvMs[grep(thisCode, groupsDF$Code)]; thisModelMad <-  storeAdMs[grep(thisCode, groupsDF$Code)] 
    if(is.na(thisModelMad) | thisModelMad==0){
      mLProps_df$propAd[g] <-0
    } else{
      mLProps_df$propAd[g] <- storeMboxplots_ad[[thisCode]][3]/thisModelMad
    }
    if(is.na(thisModelMjuv) | thisModelMjuv==0){
      mLProps_df$propJuv[g] <-0
    } else{
      mLProps_df$propJuv[g] <- storeMboxplots_juv[[thisCode]][3]/thisModelMjuv
    }
    
  }
}

# copy these over the allTheRankings.csv







