#bring in the run, the observed data, the catch history and plot together
source(paste(DIR$'General functions',"getCIfromCV.R", sep=""))
source(paste(DIR$'General functions',"getCIfromCVs.R", sep=""))
nboxes<-30

## for groups with a stock assesment and/or CPUE read in here and plot with biomass from Atlantis
SApath<-paste(DIR$'Tables', "StockAssessment\\", sep=""); CPUEpath<-paste(DIR$'Tables',"CPUE\\", sep="")

mg_2_tonne<-2e-8; X_CN<-5.7
catchYears<-seq(1900,2014); ny<-length(catchYears)

#################################### TRACERS
#plot all tracers for a given box and layer
this_run<-"base"
base_out<-"BASE"
this_out<-"FISH"

burnin<-35 #number of years to skip in plot

this_path<-paste(DIR$'Base',"ATLANTISmodels\\ArchivedModels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")
baseOutPath<-paste(this_path,"output",base_out,"\\", sep="")
r<-""
biolLines<-readLines(paste(this_path,"..\\..\\CRAM_BH_hybrid_biol.prm", sep=""))

#read in groups csv file to get codes for fished groups
groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
catchCodes<-groupsDF$Code[groupsDF$IsFished==1]; ncg<-length(catchCodes)

catchPath<-paste(this_path,"inputs\\catch_history\\",sep="")
#read in catch ARRAY
catch_array<-read.csv(paste(DIR$'Tables',"ALLcatchHistories.csv",sep=""))

# plotPath<-paste(DIR$'Figures',"Validation\\",sep="") ## use this to overwrite the figure used in the report
plotPath<-paste(DIR$'Figures',"\\Validation\\testing\\",this_out,sep="")

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1900
fishingStartYear<-1900
modelStartYear<-1900

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
ntsteps<-dim(thisVol)[3]

BaseNC.nc<-nc_open(paste(baseOutPath,"output.nc",sep=""))
baseVol<-ncvar_get(BaseNC.nc, "volume")

nts<-dim(thisVol)[3]-burnin #number of timesteps
cat(paste(nts,"\n"))
xLabsTemp<-seq(0,(nts*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0

#get all tracer names
allTracers<-sort(names(ThisNC.nc$var))

B0data<-read.csv(paste(this_path,"..\\CRAM_B0.csv",sep=""))

topCatchYears<-seq(1950,2014); biomassYears<-topCatchYears
topCatchYears<-seq(1900,2014); biomassYears<-topCatchYears

catchBarColor<-myGrey_trans

get_age_mat<-function(x){
  #x is species group code
  thisVar<-paste(x,"_age_mat", sep="")
  y<-grep(thisVar,biolLines); z<-get_first_number(biolLines[y])
  return(z)
}
## store SSB tracers for each group -only fill with age-structured
storeSSBbase<-array(NA, dim=c(ng,nts)); storeSSBfish<-storeSSBbase
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    cat(as.character(thisCode))
    thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode],side="both")
    ##get age mature so know which cohorts are part of SSB
    this_age_mat<-get_age_mat(thisCode)
    #if only 2 cohorts, just the adults are SSB
    if(thisNumCohorts==2){
      thisVar<-paste(thisName,"2_Nums", sep=""); tempNums<-ncvar_get(BaseNC.nc,thisVar); fishNums<-ncvar_get(ThisNC.nc,thisVar);
      thisVar<-paste(thisName,"2_ResN", sep=""); tempResN<-ncvar_get(BaseNC.nc,thisVar); fishResN<-ncvar_get(ThisNC.nc,thisVar);
      thisVar<-paste(thisName,"2_StructN", sep=""); tempStructN<-ncvar_get(BaseNC.nc,thisVar); fishStructN<-ncvar_get(ThisNC.nc,thisVar);
      thisSSB<-apply(tempNums*(tempResN + tempStructN),3,sum)*mg_2_tonne*X_CN; 
      tempFishSSB<-apply(fishNums*(fishResN + fishStructN),3,sum)*mg_2_tonne*X_CN;
      storeSSBbase[g,]<-thisSSB[burnin:(nts+burnin-1)]; storeSSBfish[g,]<-tempFishSSB[burnin:(nts+burnin-1)]
    } else{
      thisSSB<-rep(0,nts); fishSSB<-thisSSB
      for(c in (this_age_mat+1):thisNumCohorts){
        thisVar<-paste(thisName,c,"_Nums", sep=""); tempNums<-ncvar_get(BaseNC.nc,thisVar); fishNums<-ncvar_get(ThisNC.nc,thisVar);
        thisVar<-paste(thisName,c,"_ResN", sep=""); tempResN<-ncvar_get(BaseNC.nc,thisVar); fishResN<-ncvar_get(ThisNC.nc,thisVar);
        thisVar<-paste(thisName,c,"_StructN", sep=""); tempStructN<-ncvar_get(BaseNC.nc,thisVar); fishStructN<-ncvar_get(ThisNC.nc,thisVar);
        tempSSB<-apply(tempNums*(tempResN + tempStructN),3,sum)*mg_2_tonne*X_CN;
        tempFishSSB<-apply(fishNums*(fishResN + fishStructN),3,sum)*mg_2_tonne*X_CN;
        thisSSB<-thisSSB +  tempSSB[burnin:(nts+burnin-1)]; fishSSB<-fishSSB + tempFishSSB[burnin:(nts+burnin-1)]
      }
      storeSSBbase[g,]<-thisSSB; storeSSBfish[g,]<-fishSSB
    }
  }
}

#write these out so can read them in elsewhere
# save(list=c("storeSSBbase"), file=paste(this_path,"SSBtracers\\baserun",sep=""))
# save(list=c("storeSSBfish"), file=paste(this_path,"SSBtracers\\basefishrun",sep=""))

thisCex<-1.5

estYears<-seq(1900,(1900+dim(storeSSBfish)[2]-1))

pdf(paste(plotPath,"SSB_CPUE_SA_LEGEND.pdf", sep=""), height=1.5, width=2.5)
par(mar=c(0,0,0,0), lend=1)
makeBlankPlot()
legend(legend=c("CRAM", "CPUE", "Stock assessment"), col=c("black", myBlue,  myRed), lwd=c(3,2,2), lty=c(1,2, 4), x="center", bty="n", seg.len=3)
dev.off()


png(paste(DIR$'Figures',"SSB_CPUE_SA_LEGEND.png", sep=""), height=2000, width=2000, bg="transparent", type="cairo")
par(mar=c(0,0,0,0), lend=1)
makeBlankPlot()
legend(legend=c("Atlantis", "CPUE", "Stock assessment"), text.col = c("white", myAqua, myGold), cex=6, col=c("white", myAqua,  myGold), lwd=5, lty=c(1,2, 4), x="center", bty="n", seg.len=5)
dev.off()

thisCode<-"ORH"
#read in 2 and join them together
temp1<-read.csv(paste(SApath,"ORH_EastSouth.csv", sep="")); temp2<-read.csv(paste(SApath, "ORH_Northwest.csv", sep=""))
temp1$year<-as.double(gsub("SSB|\\[|]","",temp1$Year)); temp2$year<-as.double(gsub("SSB|\\[|]", "", temp2$Year))
allYears<-min(c(temp1$year, temp2$year), na.rm=TRUE) : max(c(temp1$year, temp2$year), na.rm=TRUE) 
thisSA<-data.frame(matrix(NA, ncol=2, nrow=length(allYears))); colnames(thisSA)<-c("year", "SSB")
thisSA$year<-allYears
thisSA$SSB <- temp1$SSB[match(allYears, temp1$year)] + temp2$SSB[match(allYears, temp2$year)]
g<-grep(thisCode,groupsDF$Code)
thisName<-gsub("_", " ", groupsDF$Name[g]); thisDesc<-groupsDF$LongName[g]
estYears<-seq(1900,2015)
firstPlotYear<-min(allYears)
altFirstPlotYear<-1980;
altFirstPlotYear<-firstPlotYear
longEst<-storeSSBfish[g,estYears>=altFirstPlotYear]
thisEst<-storeSSBfish[g,estYears>=firstPlotYear] ; meanEst<-mean(thisEst, na.rm=TRUE)
scaledSA<- meanEst * (thisSA$SSB / mean( thisSA$SSB, na.rm=TRUE))
absSA<-thisSA$SSB[match(estYears[estYears>=firstPlotYear], thisSA$year)]; compareAbs <- thisEst/absSA
thisMin<-min(c(absSA, thisEst), na.rm=TRUE); thisMax<-max(c(absSA, thisEst), na.rm=TRUE)
plot(x=absSA, y=thisEst, pch=20, col=myBlue); points(x=c(thisMin, thisMax), y=c(thisMin, thisMax), type="l", lty=2, col="red") 
plot(absSA,type="l", ylim=c(0, thisMax)); points(thisEst,type="l",lty=2)

compareAbsList<-NULL;
compareAbsList[[thisCode]]<-compareAbs


thisYmax<-max(c(thisEst,  scaledSA))

pdf(paste(plotPath,thisCode,"SSB_CPUE_SA.pdf", sep=""), height=3.5, width=5.5)
par(mar=c(2.5,4.5,1.5,1))
plot(x=estYears[estYears>=firstPlotYear], y=thisEst, xlim=c(firstPlotYear, max(estYears)), type="l", ylim=c(0,thisYmax), ylab="SSB (tonnes)", xlab="", 
     cex.axis=thisCex, cex.lab=thisCex)
points(x=thisSA$year, y=scaledSA, type="l", lwd=2, lty=4, col=myRed)
points(x=estYears[estYears>=firstPlotYear], y=thisEst,lwd=3, type="l")
mtext(thisDesc, side=3, adj=0, cex=thisCex)
subplot(boxplot(compareAbsList[[thisCode]], bty="n", ylab="", ylim=c(0,2), lwd=0.5, outline=FALSE), x=1930,y=100000, size=c(1,0.8))

dev.off()

## compare absolute

# ##redo as perc B0
# CRAMB0<-storeSSBfish[g,1]; CRAMpercB0<-100*(thisEst/CRAMB0)
# SAB0<-thisSA$SSB[1]; SApercB0<-100*(thisSA$SSB/SAB0)
# thisYmax<-max(c(CRAMpercB0, SApercB0), na.rm=TRUE)
# plot(x=estYears[estYears>=firstPlotYear], y=CRAMpercB0, xlim=c(firstPlotYear, max(estYears)), type="l", ylim=c(0,thisYmax), ylab=expression("SSB(%"~B[0]~")"), xlab="", 
#      cex.axis=thisCex, cex.lab=thisCex)
# points(x=thisSA$year, y=SApercB0, type="l", lwd=2, lty=4, col=myRed)
# points(x=estYears[estYears>=firstPlotYear], y=CRAMpercB0,lwd=3, type="l")
# mtext(thisDesc, side=3, adj=0, cex=thisCex)


thisCode<-"IVS"
thisCode<-"HAK"
# thisCode<-"LIN"
# thisCode<-"HOK"
# 
g<-grep(thisCode,groupsDF$Code)
thisName<-gsub("_", " ", groupsDF$Name[g]); thisDesc<-groupsDF$LongName[g]
thisCPUE<-read.csv(paste(CPUEpath, thisCode,"_cpue.csv", sep=""))
thisSA<-read.csv(paste(SApath, thisCode, "_SSB.csv", sep=""))

firstPlotYear<-min(c(thisCPUE$year, thisSA$year))
firstPlotYear<-min(c(thisCPUE$year, thisSA$year, 1980))
altFirstPlotYear<-1980;
longEst<-storeSSBfish[g,estYears>=altFirstPlotYear]
thisEst<-storeSSBfish[g,estYears>=firstPlotYear] ; meanEst<-mean(thisEst, na.rm=TRUE)
meanEst4CPUE<-mean(storeSSBfish[g,estYears %in% thisCPUE$year], na.rm=TRUE)

scaledCPUE<-meanEst4CPUE * (thisCPUE$index/mean(thisCPUE$index, na.rm=TRUE))
scaledSA<- meanEst * (thisSA$SSB / mean( thisSA$SSB, na.rm=TRUE))

thisYmax<-max(c(thisEst, scaledCPUE, scaledSA))


## compare absolute
absSA<-thisSA$SSB[match(estYears[estYears>=firstPlotYear], thisSA$year)]; 
if(thisCode=="HOK"){
  absSA<-thisSA$SSB[match(estYears[estYears>=firstPlotYear], thisSA$year)]*1000; # convert  
}
compareAbs <- thisEst/absSA
thisMin<-min(c(absSA, thisEst), na.rm=TRUE); thisMax<-max(c(absSA, thisEst), na.rm=TRUE)
plot(x=absSA, y=thisEst, pch=20, col=myBlue); points(x=c(thisMin, thisMax), y=c(thisMin, thisMax), type="l", lty=2, col="red") 
plot(absSA,type="l", ylim=c(0, thisMax)); points(thisEst,type="l",lty=2)
summary(compareAbs)
compareAbsList[[thisCode]]<-compareAbs

## plot boxplots of compare absolute values
plot(1:5, type="n", xaxt="n", xlab="", ylim=c(0,2))
for(i in 1:5){
  boxplot(x=compareAbsList[[i]], at=i, add=TRUE, outline=FALSE)
}

pdf(paste(plotPath,thisCode,"SSB_CPUE_SA.pdf", sep=""), height=3.5, width=5.5)
par(mar=c(2.5,4.5,1.5,1))
plot(x=estYears[estYears>=firstPlotYear], y=thisEst, xlim=c(firstPlotYear, max(estYears)), type="l", ylim=c(0,thisYmax), ylab="SSB (tonnes)", xlab="", 
     cex.axis=thisCex, cex.lab=thisCex)
points(x=thisCPUE$year, y=scaledCPUE, type="l", lwd=2, lty=2, col=myBlue)
points(x=thisSA$year, y=scaledSA, type="l", lwd=2, lty=4, col=myRed)
points(x=estYears[estYears>=firstPlotYear], y=thisEst,lwd=3, type="l")
if(thisCode=="IVS"){thisDesc <-"Invertebrate scavangers (commercial)"}
mtext(thisDesc, side=3, adj=0, cex=thisCex)
if(thisCode=="IVS"){
  subplot(boxplot(compareAbsList[[thisCode]], bty="n", ylab="", ylim=c(0,2), lwd=0.5), x=1987,y=5000, size=c(1,0.8))
}
if(thisCode=="HAK"){
  subplot(boxplot(compareAbsList[[thisCode]], bty="n", ylab="", ylim=c(0,2), lwd=0.5), x=1983,y=12000, size=c(1,0.8))
  
}
if(thisCode=="LIN"){
  subplot(boxplot(compareAbsList[[thisCode]], bty="n", ylab="", ylim=c(0,2), lwd=0.5,outline=FALSE), x=2010,y=70000, size=c(1,0.8))
  
}
if(thisCode=="HOK"){
  subplot(boxplot(compareAbsList[[thisCode]], bty="n", ylab="", ylim=c(0,2), lwd=0.5,outline=FALSE), x=2010,y=650000, size=c(1,0.8))
  
}
dev.off()

# for the poster
plotLineCol <- "white"
# png(paste(DIR$'Figures', "HakeSSB4poster.png", sep=""), height=1000, width=2000, bg="transparent")
# par(mar=c(6,20,1.5,1), las=1)
# plot(x=estYears[estYears>=firstPlotYear], y=thisEst, xlim=c(firstPlotYear, max(estYears)), type="l", ylim=c(0,thisYmax), ylab="", xlab="", 
#      cex.axis=thisCex, cex.lab=thisCex, xaxt="n", yaxt="n", bty="n")
# axis(at=seq(1980,2010, by=10), labels=seq(1980, 2010, by=10), 
#      side=1, col.axis=plotLineCol, col.ticks=plotLineCol, col.lab=plotLineCol, cex=10, font=1, cex.axis=4, col=plotLineCol)
# axis(at=seq(0,60000, by=10000), labels=seq(0, 60000, by=10000), 
#      side=2, col.axis=plotLineCol, col.ticks=plotLineCol, col.lab=plotLineCol, cex=10, font=1, cex.axis=4, lwd=5, col=plotLineCol)
# points(x=thisCPUE$year, y=scaledCPUE, type="l", lwd=7, lty=2, col=myAqua)
# points(x=thisSA$year, y=scaledSA, type="l", lwd=7, lty=4, col=myGold)
# points(x=estYears[estYears>=firstPlotYear], y=thisEst,lwd=7, type="l", col="white")
# par(las=0)
# mtext("Biomass (tonnes)", side=2, adj=0.5, line=15, cex=5, col=plotLineCol)
# dev.off()
thisCex<-8
plotFile<-paste(DIR$'Figures',"HakeSSB4poster.png",sep="")
png(plotFile, width=3000, height=2000, bg="transparent", type="cairo")
par(mar=c(15,30,10,0.5), las=1)
plot(x=estYears[estYears>=firstPlotYear], y=thisEst, xlim=c(firstPlotYear, max(estYears)), type="l", ylim=c(0,thisYmax), ylab="", xlab="", 
     cex.axis=thisCex, cex.lab=thisCex, xaxt="n", yaxt="n", bty="n")
axis(at=seq(1980,2010, by=10), labels=seq(1980, 2010, by=10), 
     side=1, col.axis=plotLineCol, col.ticks=plotLineCol, col.lab=plotLineCol, cex=thisCex, font=1, cex.axis=8, col=plotLineCol)
axis(at=seq(0,60000, by=10000), labels=seq(0, 60000, by=10000), 
     side=2, col.axis=plotLineCol, col.ticks=plotLineCol, col.lab=plotLineCol, cex=thisCex, font=1, cex.axis=8, lwd=5, col=plotLineCol)
points(x=thisCPUE$year, y=scaledCPUE, type="l", lwd=7, lty=2, col=myAqua)
points(x=thisSA$year, y=scaledSA, type="l", lwd=7, lty=4, col=myGold)
points(x=estYears[estYears>=firstPlotYear], y=thisEst,lwd=7, type="l", col="white")
par(las=0)
mtext(paste("Biomass trajectories - hake", sep=""), font=3, col=myOrange,side=3, adj=0, cex=thisCex*1.2)
mtext("Biomass (tonnes)", side=2, line=25, adj=0.5, cex=thisCex, col="white", font=1)
dev.off()


##################
thisCode<-"CRA"

g<-grep(thisCode,groupsDF$Code)
thisName<-gsub("_", " ", groupsDF$Name[g]); thisDesc<-groupsDF$LongName[g]
thisCPUE<-read.csv(paste(CPUEpath, thisCode,"_cpue.csv", sep=""))
index1<-thisCPUE$index1; index2<-thisCPUE$index2

firstPlotYear<-min(c(thisCPUE$year, thisSA$year))
altFirstPlotYear<-1980;
longEst<-storeSSBfish[g,estYears>=altFirstPlotYear]
thisEst<-storeSSBfish[g,estYears>=firstPlotYear] ; meanEst<-mean(thisEst, na.rm=TRUE)
rm(scaledSA)

scaledCPUE1<-meanEst * (index1/mean(index1, na.rm=TRUE))
scaledCPUE2<-meanEst * (index2/mean(index2, na.rm=TRUE))
# scaledSA<- meanEst * (thisSA$SSB / mean( thisSA$SSB, na.rm=TRUE))

thisYmax<-max(c(thisEst, scaledCPUE1, scaledCPUE2))

pdf(paste(plotPath,thisCode,"SSB_CPUE_SA.pdf", sep=""), height=3.5, width=6)
par(mar=c(2.5,4.5,1.5,1))
plot(x=estYears[estYears>=firstPlotYear], y=thisEst, type="l", ylim=c(0,thisYmax), ylab="SSB (tonnes)", xlab="", 
     cex.axis=thisCex, cex.lab=thisCex)
points(x=thisCPUE$year, y=scaledCPUE1, type="l", lwd=2, lty=2, col=myBlue)
points(x=thisCPUE$year, y=scaledCPUE2, type="l", lwd=2, lty=4, col=myPurple)
points(x=estYears[estYears>=firstPlotYear], y=thisEst,lwd=3, type="l")
mtext(thisDesc, side=3, adj=0, cex=thisCex)
dev.off()


# 
# par(mar=c(2.5,4.5,1.5,1))
# plot(x=estYears[estYears>=altFirstPlotYear], y=longEst, xlim=c(altFirstPlotYear, max(estYears)), type="l", ylim=c(0,thisYmax), ylab="SSB (tonnes)", xlab="", 
#      cex.axis=thisCex, cex.lab=thisCex)
# points(x=thisCPUE$year, y=scaledCPUE, type="l", lwd=2, lty=2, col=myBlue)
# # points(x=thisSA$year, y=scaledSA, type="l", lwd=2, lty=3, col=myRed)
# points(x=estYears[estYears>=altFirstPlotYear], y=longEst,lwd=3, type="l")
# mtext(thisDesc, side=3, adj=0, cex=thisCex)

# rm(thisSA)










