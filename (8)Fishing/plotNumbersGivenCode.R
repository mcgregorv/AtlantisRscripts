#plot all tracers for a given box and layer
this_run<-"base"

this_out<-"SENSselBASEfish"
# this_out<-"CETmum"

mg_2_tonne<-2e-8; X_CN<-5.7

burnin<-35 #number of years to skip in plot

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")
r<-""

#read in B0's 
thisB0df<-read.csv(paste(this_path,"..\\CRAM_B0.csv",sep=""))

groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

plotPath<-paste(this_path,"..\\Figures\\",this_run,"\\",this_out,"",sep="")

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1900
fishingStartYear<-1900
modelStartYear<-1900

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
ntsteps<-dim(thisVol)[3]


nts<-dim(thisVol)[3]-burnin #number of timesteps
cat(paste(nts,"\n"))
xLabsTemp<-seq(0,(nts*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0

#get all tracer names
allTracers<-sort(names(ThisNC.nc$var))

cohortCols<-colorRampPalette(colors=c(myYellow,myGreen,myAqua,myBlue,myPurple,myRed,"red"))(10)

thisCode<-"EID"; thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode], side="both"); thisNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisCode]
#biomass first
thisTracer<-paste(thisName,"_N", sep="")

TotCatchNC.nc<-nc_open(paste(outPath,"outputTOTCATCH.nc",sep=""))
temp<-names(TotCatchNC.nc$var); catchTracers<-temp[grep("_Catch",temp)]

testRFI<-ncvar_get(TotCatchNC.nc,"Tot_HOK_Catch")
xx<-apply(testRFI,2,sum)
testRFICatch<-xx[(burnin+1):length(xx)]
catchAxis<-pretty(seq(0,max(testRFICatch), length.out=5))


plotsFile<-paste(plotPath,"BiomassTracer_",thisCode,".pdf",sep="")
pdf(plotsFile, height=4, width=5)
par(mar=c(3,5,2,5), las=1)
temp<-ncvar_get(ThisNC.nc,thisTracer)
xx<-apply(temp*thisVol,3,sum)*(mg_2_tonne*X_CN); thisymax<-max(xx)
plot(x=seq(year0,(year0+length(testRFICatch)-1)), y=testRFICatch, type="h", col=myGrey_trans, lwd=5, lend=1, xlim=c(1970,2014), 
           yaxt="n", xaxt="n", xlab="", ylab="")
axis(catchAxis, side=4)
par(new=TRUE)
plot(x=seq(year0,(year0+length(xx)-1)),y=xx,type="l",col=myGreen,lwd=2.5,ylim=c(thisymin,thisymax), xlim=c(1970,2014), xlab="", ylab="")
par(las=0)
mtext("Catch (tonnes)", side=4, adj=0.5, line=3)
mtext("Biomass (tonnes)", side=2, adj=0.5, line=3.2)
mtext(thisName,side=3,adj=0,font=2)
dev.off()

#now do numbers by cohort
numTracers<-paste(thisName,seq(1,thisNumCohorts),"_Nums", sep="")
plotsFile<-paste(plotPath,"NumberTracers_",thisCode,"_relative.pdf",sep="")
pdf(plotsFile, height=4, width=5)
par(mar=c(3,5,2,5), las=1)
temp<-ncvar_get(ThisNC.nc,numTracers[1])
xx<-apply(temp,3,sum); thisY<-xx/xx[35]; thisymax<-max(xx)
thisX<-seq(year0,(year0+length(testRFICatch)-1))
plot(x=thisX, y=testRFICatch, type="h", col=myGrey_trans, lwd=5, lend=1, xlim=c(1900,2014), yaxt="n", xaxt="n", xlab="", ylab="")
axis(catchAxis, side=4)
par(new=TRUE)
thisY<-xx/xx[35]
plot(x=thisX,y=thisY[burnin:(length(xx)-1)],type="l",col=cohortCols[1],lwd=2.5,ylim=c(0,3), xlab="", ylab="", xlim=c(1900,2014))
for(c in 2:thisNumCohorts){
  temp<-ncvar_get(ThisNC.nc,numTracers[c])
  xx<-apply(temp,3,sum); thisymax<-max(xx)
  thisY<-xx/xx[35]
  points(x=thisX,y=thisY[burnin:(length(xx)-1)],type="l",col=cohortCols[c],lwd=2.5)
}
par(las=0)
mtext("Catch (tonnes)", side=4, adj=0.5, line=3)
mtext("Numbers (relative)", side=2, adj=0.5, line=4.2)
mtext(thisName,side=3,adj=0,font=2)
dev.off()

plotsFile<-paste(plotPath,"NumberTracers_",thisCode,"_relative.pdf",sep="")
pdf(plotsFile, height=4, width=5)
par(mar=c(3,5,2,5), las=1)
temp<-ncvar_get(ThisNC.nc,numTracers[1])
xx<-apply(temp,3,sum); thisY<-xx/xx[35]; thisymax<-max(xx)
thisX<-seq(year0,(year0+length(testRFICatch)-1))
plot(x=thisX, y=testRFICatch, type="h", col=myGrey_trans, lwd=5, lend=1, xlim=c(1970,2014), yaxt="n", xaxt="n", xlab="", ylab="")
axis(catchAxis, side=4)
par(new=TRUE)
thisY<-xx/xx[105]
plot(x=thisX,y=xx[burnin:(length(xx)-1)],type="l",col=cohortCols[1],lwd=2.5,ylim=c(0,max(xx,na.rm=TRUE)), xlab="", ylab="", xlim=c(1970,2014))
for(c in 2:thisNumCohorts){
  temp<-ncvar_get(ThisNC.nc,numTracers[c])
  xx<-apply(temp,3,sum); thisymax<-max(xx)
  thisY<-xx/xx[105]
  points(x=thisX,y=xx[burnin:(length(xx)-1)],type="l",col=cohortCols[c],lwd=2.5)
}
par(las=0)
mtext("Catch (tonnes)", side=4, adj=0.5, line=3)
mtext("Numbers", side=2, adj=0.5, line=4.2)
mtext(thisName,side=3,adj=0,font=2)
dev.off()
