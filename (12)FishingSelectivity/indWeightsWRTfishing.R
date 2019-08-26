#read in fishing and non-fishing model runs
#compare weights of individuals by age-class as fishing kicks in
#this is a start - will want mult runs with varying selectivities and compare these

mg_2_tonne<- 0.00000002; X_CN<-5.7  

#plot all tracers for a given box and layer
this_run<-"base"
this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
plotPath<-paste(this_path,"..\\Figures\\CompareRuns_FishSENS_burnin_",sep="")

# this_out<-c("Mdiets", "Mdietsfish", "SENsel1", "SENsel2"); 
this_out<-c("SENSselBASE", "SENSselBASEfish", "SENsel1", "SENsel2", "SENsel3");
# this_out<-c("Mdiets", "Mdietsfish", "SENsel1", "SENsel2"); 
this_out<-c("Mdietsfish", "SENsel1", "SENsel2");
fish_out<-this_out[2]

nruns<-length(this_out)
burnin<-rep(35,nruns) #number of years to skip in plot
burnin<-c(0,35,35)

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1900
fishingStartYear<-1900
modelStartYear<-1900

nc_list<-NULL; nts_list<-NULL; min_nts<-1e+12
for(r in 1:nruns){
  outPath<-paste(this_path,"output",this_out[r],"\\",sep="")
  
  nc_list[[r]]<-nc_open(paste(outPath,"output.nc",sep=""))
  thisVol<-ncvar_get(nc_list[[r]],"volume")
  thisDz<-ncvar_get(nc_list[[r]],"dz")
  nts_list[[r]]<-dim(thisVol)[3]-burnin[r] #number of timesteps
  if(nts_list[[r]]<min_nts){min_nts<-nts_list[[r]]}
}

xLabsTemp<-seq(0,((min_nts-burnin[1])*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0

#doing for RFI first
thisCode<-"RFI"; thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode], side="both"); thisNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisCode]

numbers_array<-array(NA,dim=c(nruns,thisNumCohorts,min_nts-burnin[1]+1))
weight_array<-numbers_array

for(r in 1:nruns){
  for(c in 1:thisNumCohorts){
    thisTracer<-paste(thisName,c,"_ResN",sep="");   thisTemp<-ncvar_get(nc_list[[r]],thisTracer);   thisRN<-apply(thisTemp,3,nonZeroMean)
    thisTracer<-paste(thisName,c,"_StructN",sep="");   thisTemp<-ncvar_get(nc_list[[r]],thisTracer); thisSN<-apply(thisTemp,3,nonZeroMean)
    temp<-((thisRN+thisSN)*mg_2_tonne*X_CN)[burnin[r]:min_nts]
    weight_array[r,c,1:length(temp)]<-temp
    #numbers
    thisTracer<-paste(thisName,c,"_Nums",sep="");   thisTemp<-ncvar_get(nc_list[[r]],thisTracer);
    temp<-apply(thisTemp,3,sum,na.rm=TRUE)[burnin[r]:min_nts]
    numbers_array[r,c,1:length(temp)]<-temp
  }
}

##read in catches to add
fishoutPath<-paste(this_path,"output",fish_out,"\\",sep="")
TotCatchNC.nc<-nc_open(paste(outPath,"outputTOTCATCH.nc",sep=""))
temp<-names(TotCatchNC.nc$var); catchTracers<-temp[grep("_Catch",temp)]

testRFI<-ncvar_get(TotCatchNC.nc,"Tot_RFI_Catch")
testRFICatch<-apply(testRFI,2,sum)[burnin:min_nts]
catchAxis<-pretty(seq(0,max(testRFICatch), length.out=5))

##get ymax and min
thisYmax<-0; thisYmin<-1e+100
for(c in 1:thisNumCohorts){
  for(r in 2:4){
    test<-weight_array[r,c,]/weight_array[1,c,]
    thisMax<-max(test,na.rm=TRUE); thisMin<-min(test, na.rm=TRUE)
    if(thisMax>thisYmax){thisYmax<-thisMax}
    if(thisMin<thisYmin){thisYmin<-thisMin}
  }
}
colByRun<-c(myGrey,colorRampPalette(colors=c(myOrange, myGreen, myBlue))(nruns-1))


pdf(paste(plotPath, thisCode,"diffInWeights_SENSsel.pdf", sep=""), height=6, width=5)
par(mfrow=c(thisNumCohorts,1), mar=c(0,4,0,4), oma=c(3,1,1,1), las=1)
for(c in 1:thisNumCohorts){
  plot(weight_array[2,c,]/weight_array[1,c,], lty=2, col=colByRun[2],type="l", lwd=2, xaxt="n",ylab="", xlab="", ylim=c(thisYmin*0.998,thisYmax*1.02))
  abline(h=1, col=myGrey_trans, lwd=2)
  axis(at=xLabsAt, labels=xLabs, side=1)
  points(weight_array[3,c,]/weight_array[1,c,],type="l", lwd=2, col=colByRun[3], lty=3)
  points(weight_array[4,c,]/weight_array[1,c,],type="l", lwd=2, col=colByRun[4], lty=4)
  par(new=TRUE)
  plot(testRFICatch, type="h", lwd=5, col=myGrey_trans, lend=1, xaxt="n", yaxt="n", ylab="", xlab="")
  axis(at=catchAxis, labels=catchAxis, side=4, col.axis=myGrey)
}
par(las=0)
mtext("Weights|fishing/Weights|no-fishing", side=2, adj=0.5, outer=TRUE, line=-0.5)
mtext("Catches (tonnes)", side=4, adj=0.5, outer=TRUE, line=-0.5, col=myGrey)
dev.off()

pdf(paste(plotPath, thisCode,"diffInNums_SENsel.pdf", sep=""), height=6, width=5)
par(mfrow=c(thisNumCohorts,1), mar=c(0,4,0,4), oma=c(3,1,1,1), las=1)
for(c in 1:thisNumCohorts){
  plot(numbers_array[2,c,]/numbers_array[1,c,], col=colByRun[2], lty=2,type="l", lwd=2, xaxt="n",ylab="", xlab="", ylim=c(0,2))
  abline(h=1, col=myGrey_trans, lwd=2)
  points(numbers_array[3,c,]/numbers_array[1,c,],type="l", lwd=2, col=colByRun[3], lty=3)
  points(numbers_array[4,c,]/numbers_array[1,c,],type="l", lwd=2, col=colByRun[4], lty=4)
  par(new=TRUE)
  plot(testRFICatch, type="h", lwd=5, col=myGrey_trans, lend=1, xaxt="n", yaxt="n", ylab="", xlab="")
  axis(at=catchAxis, labels=catchAxis, side=4, col.axis=myGrey)
}
axis(at=xLabsAt, labels=xLabs, side=1)
par(las=0)
mtext("Numbers|fishing/Numbers|no-fishing", side=2, adj=0.5, outer=TRUE, line=-0.5)
mtext("Catches (tonnes)", side=4, adj=0.5, outer=TRUE, line=-0.5, col=myGrey)
dev.off()

## compare total biomass and total numbers between the runs
totalBiomass<-apply(weight_array*numbers_array, c(1,3), sum)
thisYmax<-max(totalBiomass, na.rm=TRUE)

pdf(paste(plotPath, thisCode,"diffInTotalBiomass_SENsel.pdf", sep=""), height=4, width=5)
par(mar=c(4,5,1,1), las=1)
plot(totalBiomass[1,], col=myGrey, lwd=2, ylim=c(0,thisYmax), xlab="Year", ylab="", xaxt="n", type="l", xlim=c(1,length(totalBiomass[1,])))
for(r in 2:nruns){
  points(totalBiomass[r,], col=colByRun[r],type="l", lwd=2, lty=r )
}
axis(at=xLabsAt, labels=xLabs, side=1)
par(las=0)
mtext("Biomass (tonnes)", side=2, line=3.5, adj=0.5)
dev.off()

thisYmax<-max(totalBiomass[,35:80], na.rm=TRUE)
colByRun<-c(myOrange,myGreen,myBlue)

pdf(paste(plotPath, thisCode,"diffInTotalBiomass_SENselZOOMIN.pdf", sep=""), height=4, width=5)
par(mar=c(4,5,1,1), las=1)
plot(totalBiomass[1,35:80], col=myOrange, lwd=2, ylim=c(0,thisYmax), xlab="Year", ylab="", xaxt="n", type="l")
for(r in 2:nruns){
  points(totalBiomass[r,35:80], col=colByRun[r],type="l", lwd=2, lty=1 )
}
axis(at=seq(35:80), labels=xLabs[36:81], side=1)
par(las=0)
mtext("Biomass (tonnes)", side=2, line=3.5, adj=0.5)
dev.off()

#plot the legend
pdf(paste(plotPath, thisCode,"diffInX_SENselLEGEND.pdf", sep=""), height=2, width=3)
par(mar=c(0,0,0,0))
makeBlankPlot()
legend(legend=c("No fishing", "Base selectivity", "Uniform selectivity", "Adult selectivity"), col=colByRun, lty=seq(1, nruns), lwd=2, x="center", bty="n", seg.len = 3)
dev.off()


totalNumbers<-apply(numbers_array, c(1,3), sum)
pdf(paste(plotPath, thisCode,"diffInTotalNumbers_SENsel.pdf", sep=""), height=4, width=5)
par(mar=c(4,5,1,1), las=1)
plot(totalNumbers[1,], col=myGrey, lwd=2, ylim=c(0,max(totalNumbers)), xlab="Year", ylab="", xaxt="n", type="l", xlim=c(1,length(totalNumbers[1,])))
for(r in 2:nruns){
  points(totalNumbers[r,], col=colByRun[r],type="l", lwd=2, lty=r )
}
axis(at=xLabsAt, labels=xLabs, side=1)
par(las=0)
mtext("Numbers", side=2, line=3.5, adj=0.5)
dev.off()

testWeight_array<-weight_array
for(r in 2:nruns){
  testWeight_array[r,,]<-weight_array[2,,]
}
testTotalBiomass<-apply(testWeight_array*numbers_array, c(1,3), sum)

pdf(paste(plotPath, thisCode,"diffInTESTTotalBiomass_SENsel.pdf", sep=""), height=4, width=5)
par(mar=c(4,5,1,1), las=1)
plot(testTotalBiomass[1,], col=myGrey, lwd=2, ylim=c(0,max(testTotalBiomass)), xlab="Year", ylab="", xaxt="n", type="l", xlim=c(1,length(testTotalBiomass[1,])))
for(r in 2:nruns){
  points(testTotalBiomass[r,], col=colByRun[r],type="l", lwd=2, lty=r )
}
axis(at=xLabsAt, labels=xLabs, side=1)
par(las=0)
mtext("Biomass (tonnes)", side=2, line=3.5, adj=0.5)
dev.off()

# 
# plot(testTotalBiomass[3,]/totalBiomass[3,], type="l")

##read in harvest file and plot catch by age distr
harvestFiles<-c("CRAM_harvest_short.prm", "CRAM_harvest_shortTestSel1.prm", "CRAM_harvest_shortTestSel2.prm"); nhr<-length(harvestFiles)
pdf(paste(plotPath, thisCode,"propCatchAtAge.pdf", sep=""), height=4, width=5)
par(mar=c(4,4,1,1), las=1)
for(h in 1:nhr){
  harvestFile<-paste(basePath, harvestFiles[h], sep="")
  harvestLines<-readLines(harvestFile)
  
  thisVar<-paste("CatchTS_agedistrib", thisCode, sep=""); temp<-harvestLines[grep(thisVar, harvestLines)+1]; thisCatchPropByAge<-get_first_number(temp,n="all")
  if(h==1){
    plot(x=seq(1,thisNumCohorts),y=thisCatchPropByAge, type="h", lwd=5, lend=1, col=colByRun[h+1], ylab="Proportion catch at age", xlab="Age class", xaxt="n")
    axis(at=seq(1,thisNumCohorts), labels=seq(1,thisNumCohorts), side=1)
  } else{
    points(x=seq(1,thisNumCohorts)+(h-1)/10,y=thisCatchPropByAge, type="h", lwd=5, lend=1, col=colByRun[h+1])
  }
}
dev.off()

pdf(paste(plotPath, thisCode,"propCatchAtAgeLEGEND.pdf", sep=""), height=2, width=4)
par(mar=c(0,0,0,0), las=1, lend=1)
makeBlankPlot()
legend(legend=c("Base selectivity", "Uniform selectivity", "Adult only selectivity"), col=c(myOrange, myGreen, myBlue), 
       lwd=5, seg.len=3, x="center", bty="n")

dev.off()

