## if 60% of mortality is from mL, what does this mean in terms of numbers..?
# give couple of examples - read in initial conditions, and calculate how many die of other causes
inputsPath <- paste(DIR$'Base', "ATLANTISmodels\\inputs\\", sep="")
plotPath <- paste(DIR$'Base',"Reports\\Chaos\\Figures\\", sep="")

load(paste(inputsPath,"..\\Chaos\\storeNumbersByGroup",sep=""))

groupsDF <- read.csv(paste(inputsPath,"..\\CRAM_Groups.csv", sep="")); ng <- dim(groupsDF)[1]
biolLines <- readLines( paste(inputsPath,"..\\CRAM_BH_hybrid_biol.prm", sep=""))

#suppose mL is 10% for juveniles and 90% for adults - how many are available as food..?
# what is a typical N - for juv and adults..?
get_age_mat<- function(thisCode){
  thisVar<-paste(thisCode,"_age_mat", sep="")
  y<-grep(thisVar,biolLines); z<-get_first_number(biolLines[y])
  return(z)
}

thisCode<-"MAC"
g <- grep(thisCode, groupsDF$Code)
thisName <- str_trim(groupsDF$Name[groupsDF$Code==thisCode], side="both")
thisNumCohorts <- groupsDF$NumCohorts[g]
thisNumsByCohort <- storeNumbersByGroup[1,g,]
this_age_mat <- get_age_mat(thisCode)
thisJuvNums <- thisNumsByCohort[1:this_age_mat]
thisAdNums <- thisNumsByCohort[(this_age_mat+1): thisNumCohorts]

#suppose M total is 0.3
thisM <- 0.3
prop2die <- 1-exp(-thisM)
num2die <- signif(signif(sum(thisNumsByCohort, na.rm=TRUE),2)*prop2die,2)
#if mL was 0.8*thisM, what would be the difference in number to die?
thisBackgroundM <- 0.8*thisM
backgroundProp2die <- 1-exp(-thisBackgroundM)
backgroupM2die <- signif(signif(sum(thisNumsByCohort, na.rm=TRUE),2)*backgroundProp2die,2)
#difference
num2die - backgroupM2die

(num2die - backgroupM2die)/num2die

## how many for each species group..?
numbersByGroup <- apply(storeNumbersByGroup[1,,],1, sum, na.rm=TRUE)
hist(numbersByGroup)

numJuvsByGroup <- rep(NA, ng); numAdsByGroup <- numJuvsByGroup
for(g in 1:ng){
  thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    thisCOde <- groupsDF$Code[g]
    this_age_mat <- get_age_mat(thisCode)
    numJuvsByGroup[g] <- sum(storeNumbersByGroup[1,g,1:this_age_mat], na.rm=TRUE)
    numAdsByGroup[g] <- sum(storeNumbersByGroup[1,g,(this_age_mat + 1): thisNumCohorts])
  }
}
plot(x=c(1,2), y=c(1,1), ylim=c(0, 8.5e+2), type="n", xaxt="n", bty="n", xlab="", xlim=c(0.5,2.5), ylab="Number of individuals (millions)")
boxplot(numJuvsByGroup/(1e+6), outline=FALSE, col=myOrange_trans, border=myOrange,at=1, add=TRUE, bty="n")
boxplot(numAdsByGroup/(1e+6), outline=FALSE, col=myBlue_trans, border = myBlue,at=2, add=TRUE)
axis(at=c(1,2), labels=c("Juveniles", "Adults"), side=1)

#######################################################
## plot this with mortality proportions
# the analyses were carried out in TCorr_mL.R - this just reads in the values required to create the plots, and creates the plots
## original plots were created in TCorr_mL.R, but formatting and layout have changed here
load(paste(inputsPath, "..\\storeChaosTemp\\storeValuesToPlot",sep=""))

# some other things to read in
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

############ on one plot
pdf(paste(plotPath,"AppliedMoverModelM_withNumbers.pdf", sep=""), height=5, width=9)
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
par(new=TRUE, fig=c((2/3), 1, 0.3, (1)), mar=c(5,4,2,0.5))
plot(x=c(1,2), y=c(1,1), ylim=c(0, 8.5e+2), type="n", xaxt="n", bty="n", xlab="", xlim=c(0.5,2.5), ylab="Number of individuals (millions)")
boxplot(numJuvsByGroup/(1e+6), outline=FALSE, col=myOrange_trans, border=myOrange,at=1, add=TRUE, bty="n")
boxplot(numAdsByGroup/(1e+6), outline=FALSE, col=myBlue_trans, border = myBlue,at=2, add=TRUE)
axis(at=c(1,2), labels=c("Juveniles", "Adults"), side=1)
dev.off()
############

hist(numJuvsByGroup/numAdsByGroup)






