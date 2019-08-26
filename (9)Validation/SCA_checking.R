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
year0<-1880
fishingStartYear<-1900
modelStartYear<-1976
newVol <- ncvar_get(newNC.nc, "volume")

xLabsTemp<-seq(burnin,(nts*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0
modelTime <- year0 + (1:nts)/numStepsPerYear

burnin <- modelStartYear-year0

lastYear <- max(modelTime)
plotYearlyIndex <- seq(burnin*numStepsPerYear, nts, by=numStepsPerYear)
# plotYearlyIndex <- seq(burnin*numStepsPerYear,nts)

# read in biomass index - think from stock assessment, need to check this
temp <- read.csv(paste(DIR$'Base',"TBGB\\Data\\Scallops\\SCABiomassEWEM1.csv", sep=""))
biomIndexData <- temp[!is.na(temp$obs) & temp$date>1976,]

# get biomass by age class - juveniles and adults
storeAgeBiomass <- array(NA, dim=c(2,nts))
thisCode<-"SCA"; thisName <- "Scallops"
for(c in 1:2){
  thisVar<-paste(thisName,c,"_Nums", sep="");  fishNums<-ncvar_get(ThisNC.nc,thisVar);
  thisVar<-paste(thisName,c,"_ResN", sep="");  fishResN<-ncvar_get(ThisNC.nc,thisVar);
  thisVar<-paste(thisName,c,"_StructN", sep="");  fishStructN<-ncvar_get(ThisNC.nc,thisVar);
  tempFishSSB<-apply(fishNums*(fishResN + fishStructN),3,sum)*mg_2_tonne*X_CN;
  storeAgeBiomass[c,]<-tempFishSSB
}


# check out spatial catches
catchNC.nc <- nc_open(paste(outPath, "outputCATCH.nc", sep=""))
catchVars <- sort(unique(names(catchNC.nc$var)))
scaVars <- catchVars[grep("SCA", catchVars)]

scaCatchData <- ncvar_get(catchNC.nc, "TAR_Catch_FC1")
sum(scaCatchData)

catchByYear <- apply(scaCatchData, 2, sum)

catchYears <- year0:(year0+dim(scaCatchData)[2])
catchYearIndex <- catchYears %in% round(modelTime)


plotObs <- (biomIndexData$obs/mean(biomIndexData$obs))*mean(storeAgeBiomass[2,modelTime>=min(biomIndexData$date)])
pdf(paste(DIR$'Figures',"Validation_TBGB",this_out,"SCA_biomassCatch.pdf", sep=""), height=4, width=5)
par(mar=c(4,4,1,4))
plot(x=modelTime[plotYearlyIndex], y=storeAgeBiomass[2,plotYearlyIndex], ylab="Biomass (tonnes)", type="l", col=myBlue, ylim=c(0, 1.3*max(storeAgeBiomass[,plotYearlyIndex])), xlab="Year")
# points(x=modelTime[plotYearlyIndex], y=storeAgeBiomass[1,plotYearlyIndex], type="l", col=myOrange)
points(x=biomIndexData$date, y=plotObs, pch=20); points(x=biomIndexData$date, y=plotObs, lty=2, type="l")
# 
# par(new=TRUE)
# plot(x=modelTime[plotYearlyIndex], y=storeAgeBiomass[1,plotYearlyIndex], type="l", col=myOrange)
# what about the EwE model? Something weird about it - wait for new version
# ewe2plot <- (biomIndexData$est/mean(biomIndexData$est)) * mean(storeAgeBiomass[2,modelTime>=min(biomIndexData$date)])
# points(x=biomIndexData$date, y=ewe2plot, type="l", col=myGrey)

par(new=TRUE, lend=1)
plot(x=catchYears[catchYearIndex], y=catchByYear[catchYearIndex], type="h",lwd=5, col=myGrey_trans, xaxt="n", yaxt="n", ylab="", xlab="")
catchAxis <- pretty(1:max(catchByYear))
axis(at=catchAxis, labels=catchAxis, col=myGrey, col.axis=myGrey, side=4)
mtext("Catch (tonnes)", side=4, col=myGrey, font=1, adj=0.5, line=2.5)
dev.off()


catchNonZeroIndex <- catchByYear>0
plot(x=catchYears[catchNonZeroIndex], y=catchByYear[catchNonZeroIndex], type="h", lwd=5)


mean(c(130000, 86700, 86200, 81000, 100000, 111000, 127000, 67000, 149000))
