source(paste(DIR$'General functions',"make_grid.r",sep=""))

#read in outputDietCheck.txt
options(stringsAsFactors = FALSE)

WC_COL<-myBlue
SD_COL<-myRed
EP_COL<-myGreen

habitat_layers<-c("WC","SED","EPIBENTHIC")
habitat_cols<-c(WC_COL,SD_COL,EP_COL)

this_run<-"base"
this_out<-paste("FISH",sep="")
SApath<-paste(DIR$'Tables', "StockAssessment\\", sep=""); CPUEpath<-paste(DIR$'Tables',"CPUE\\", sep="")

# 
basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
outPath<-paste(basePath,this_run,"\\","output",this_out,"\\",sep="")
burnin<-110 # starts plots at 1980 (115)
## read in SSB arra
load(paste(basePath,this_run,"\\SSBtracers\\basefishrun",sep=""))
groupsDF<-read.csv(paste(basePath,"\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
ng<-dim(groupsDF)[1]

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1865+burnin
modelStartYear<-year0

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
ntsteps<-dim(thisVol)[3]

nts<-dim(thisVol)[3]-burnin #number of timesteps
cat(paste(nts,"\n"))

xLabsTemp<-seq(0,(nts*daysTimeStep-1),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear +1
xLabs<-xLabsTemp+year0

# plotPath<-paste(DIR$'Figures',"timeVaryingDynamics\\", sep="") ## paper version
plotPath<-paste(basePath,"\\Figures\\",this_run,"\\", this_out,sep="")

thisPrey<-"LIN"
##read in YCS from stock assessment
thisYCS<-read.csv(paste(DIR$'Tables', "StockAssessment\\",thisPrey,"_YCS.csv", sep=""))
thisYCS<-thisYCS[thisYCS$year>=year0,]
YCSyears<-thisYCS$year


thisNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisPrey]; thisPreyName<-str_trim(groupsDF$Name[groupsDF$Code==thisPrey])
thisNumbersArray<-array(NA, dim=c(nts,thisNumCohorts)); thisWeightsArray<-0*thisNumbersArray
for(c in 1:thisNumCohorts){
  thisTracer<-paste(thisPreyName,c,"_Nums", sep=""); thisData<-ncvar_get(ThisNC.nc, thisTracer)
  thisNums<-apply(thisData, 3, sum); thisNumbersArray[, c]<-thisNums[burnin:(nts+burnin-1)]
  thisTracer<-paste(thisPreyName,c,"_ResN", sep=""); thisData<-ncvar_get(ThisNC.nc, thisTracer)
  thisResNs<-apply(thisData, 3, nonZeroMean); 
  thisTracer<-paste(thisPreyName,c,"_StructN", sep=""); thisData<-ncvar_get(ThisNC.nc, thisTracer)
  thisStructNs<-apply(thisData, 3, nonZeroMean); 
  thisWeightsArray[, c]<-(thisResNs + thisStructNs)[burnin:(nts+burnin-1)]
}

## outputSpecificMort
outSpecMort<-read.csv(paste(outPath, "outputSpecificMort.txt", sep=""),skip=0,sep=" ",header=TRUE)
thisSpecMort<-outSpecMort[outSpecMort[,1]>=((burnin-1)*365), grep(thisPrey, colnames(outSpecMort))]
allSpecMort<-outSpecMort[, grep(thisPrey, colnames(outSpecMort))]
predMort<-thisSpecMort[,grep("M2", colnames(thisSpecMort))]
# these are mg N consumed per year as a proportin of numbers. 
# In the Atlantis code, the numerator should have been turned to numbers and it wasnt (remaind as mg N)
## but then was divided by the numbers
## hence, we can fix it by dividing by the weights of individuals here :) 
predMortNums<-(predMort/thisWeightsArray)
predM<-apply(predMortNums,c(1,2), FUN=function(x){(-1)*log(1-x)})
colnames(predM)<-seq(1,thisNumCohorts)
toPlot<-melt(predM)
mortYears<-seq(1865+burnin,1865+burnin+nts-1)
toPlot$Year<-rep(mortYears,thisNumCohorts)
colnames(toPlot)<-c("timestep", "Ageclass", "M", "Year")
toPlot$Ageclass<-as.factor(toPlot$Ageclass)

bp<-ggplot(data = toPlot, aes(x = Year, fill = Ageclass, y = M)) + 
  geom_bar(stat = 'identity')
bp 

## need background mortality for ling
otherMort<-allSpecMort[,grep("M1", colnames(thisSpecMort))]
potherMortNums<-(otherMort)
otherM<-apply(potherMortNums,c(1,2), FUN=function(x){(-1)*log(1-x)})
meanOtherMortByYear<-apply(otherM, 1, mean)

meanPredMortByYear<-apply(predM, 1, mean)


## other mort is only linear mortality so is actually a constant rate
## it looks like it changes in the outputs as the rate is approximated based on numbers at the start of the timestep,
## so when they are under a lot of other mortality (ie fishing) it is not so accurate and changes
# take the mean M(other) from the first 15 timesteps here
meanOther<-mean(meanOtherMortByYear[1:15])

meanMortByYear<-meanPredMortByYear + meanOther


thisCode<-thisPrey
g<-grep(thisCode,groupsDF$Code)
thisName<-gsub("_", " ", groupsDF$Name[g]); thisDesc<-groupsDF$LongName[g]
thisCPUE<-read.csv(paste(CPUEpath, thisCode,"_cpue.csv", sep=""))
thisSA<-read.csv(paste(SApath, thisCode, "_SSB.csv", sep=""))

firstPlotYear<-min(c(thisCPUE$year, thisSA$year))
firstPlotYear<-min(c(thisCPUE$year, thisSA$year, 1975))
altFirstPlotYear<-1975;


fullEstYears<-seq(1900,2015)
thisEst<-storeSSBfish[g,fullEstYears>=firstPlotYear] ; meanEst<-mean(thisEst, na.rm=TRUE)
estYears<-fullEstYears[fullEstYears>=firstPlotYear]
meanEst4CPUE<-mean(storeSSBfish[g,estYears %in% thisCPUE$year], na.rm=TRUE)
scaledSA<- meanEst * (thisSA$SSB / mean( thisSA$SSB, na.rm=TRUE))

thisYmax<-max(c(thisEst,  scaledSA))

plotYCS<-thisYCS$YCS[match(estYears, thisYCS$year)]
plotYCS<-(plotYCS)/mean(plotYCS, na.rm=TRUE)

YCSAxis<-pretty(seq(min(plotYCS, na.rm=TRUE), max(plotYCS, na.rm=TRUE), length.out = 5))


pdf(paste(plotPath, thisPrey,"mortyVsRecruitment.pdf", sep=""), height=6, width=8)
par(mar=c(0,5,0,4), las=1, mfrow=c(2,1), oma=c(4,1,1,1))
plot(x=estYears, y=thisEst, xlim=c(firstPlotYear, max(estYears)), type="l", ylim=c(0,thisYmax*1.05), ylab="", xlab="", 
     cex.axis=thisCex, cex.lab=thisCex, xaxt="n")
points(x=thisSA$year[thisSA$year %in% estYears], y=scaledSA[thisSA$year %in% estYears], type="l", lwd=2, lty=4, col="red")
points(x=estYears[estYears>=firstPlotYear], y=thisEst,lwd=3, type="l")
mtext(thisDesc, side=3, adj=0, cex=thisCex)
par(las=0)
mtext("SSB (tonnes)", side=2,adj=0.5, line=4.5, cex=thisCex)
legend(legend=c("Atlantis", "Stock assessment"), col=c("black", "red"), lwd=2, lty=c(1,4), x="bottomleft", bty="n", seg.len=3, cex=thisCex)
## YCS and M
par(las=1)
plot(as.double(meanMortByYear), type="l", col="black",lwd=2,xaxt="n",  ylab="", xlab="", cex.axis=thisCex, cex.lab=thisCex)
axis(at=xLabsAt,labels=xLabs,side=1, cex.axis=thisCex, cex.lab=thisCex)
par(las=0)
mtext("M", side=2, adj=0.5, line=4, col="black", cex=thisCex)
par(new=TRUE, las=1)
plot(plotYCS, type="l", lwd=2, col="red", yaxt="n",  xlab="", ylab="", lty=4, xaxt="n", cex.axis=thisCex, cex.lab=thisCex)
axis(at=YCSAxis, labels=YCSAxis, side=4, col.axis="red",  cex.axis=thisCex,cex.lab=thisCex)
par(las=0)
mtext("Year class strength", side=4, adj=0.5, line=3, col="red", cex=thisCex)
# axis(at=YCSAxis, labels = YCSAxis,side=4, col="red", col.axis="red")
dev.off()

