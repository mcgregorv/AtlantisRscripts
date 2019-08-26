#plot all tracers for a given box and layer
this_run<-"base"

this_out<-"SENselALL1"
# this_out<-"SENSselBASE"
# this_out<-"BALBH"
this_out<-"XXX_FSMGA1"
this_out<-"FISHmum_XXXFC2"

mg_2_tonne<-2e-8; X_CN<-5.7

burnin<-35 #number of years to skip in plot

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")

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

cohortCols<-colorRampPalette(colors=c(myGold,myGreen,myAqua,myBlue,myPurple,myRed,"red"))(10)
allTracers<-sort(unique(names(ThisNC.nc$var)))

thisCode<-"ORH"
##pul out the tracers and store first
storeNums<-array(NA,dim=c(ng,thisNumCohorts, nts)); storeResN<-storeNums; storeStructN<-storeResN

for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisCode]
  
  if(thisNumCohorts>1){
    thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode], side="both")
    
    thisNumTracers<-paste(thisName,seq(1,thisNumCohorts),"_Nums",sep="")
    thisResNtracers<-paste(thisName,seq(1,thisNumCohorts),"_ResN", sep="")
    thisStructNtracers<-gsub("ResN","StructN", thisResNtracers)
    
    #do numbers tracer first
    for(c in 1:thisNumCohorts){
      thisTracer<-thisNumTracers[c]; temp<-ncvar_get(ThisNC.nc, thisTracer)
      storeNums[g,c,]<-apply(temp,3,sum)[(burnin+1):(burnin+nts)]
      thisTracer<-thisResNtracers[c]; temp<-ncvar_get(ThisNC.nc, thisTracer)
      storeResN[g,c,]<-apply(temp,3,nonZeroMean)[(burnin+1):(burnin+nts)]
      thisTracer<-thisStructNtracers[c]; temp<-ncvar_get(ThisNC.nc, thisTracer)
      storeStructN[g,c,]<-apply(temp,3,nonZeroMean)[(burnin+1):(burnin+nts)]
    }
  }
}
##########################################
##this section does 4 plots for groups that are fished - numbers tracers for each age class, SSB, Bcur, selectivity
for(g in 1:ng){
thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
thisIsFished<-groupsDF$IsFished[g]
if(thisIsFished==1){
  thisMax<-max(storeNums[g,,], na.rm=TRUE)
  thisPlotFile<-paste(plotPath, "SSBandSuch",thisCode,".jpg", sep="")
  jpeg(thisPlotFile)
  par(mfrow=c(2,2), mar=c(4,4,1,1), oma=c(1,1,1,1))
  plot(x=seq(1,nts),y=storeNums[g,1,],type="l",col=cohortCols[1],lwd=2.5,ylim=c(0,thisMax),ylab="",xlab="Year",xaxt="n")
  axis(at=xLabsAt,labels=xLabs,side=1)
  for(c in 2:thisNumCohorts){
    points(x=seq(1,nts),y=storeNums[g,c,],type="l",col=cohortCols[c],lwd=2.5)
  }
  mtext("Numbers", side=2, adj=0.5, line=3)
  
  ##get age mature
  thisVar<-paste(thisCode, "_age_mat", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageMat<-get_first_number(temp)
  if(thisNumCohorts==2){
    ssbIndex<-2
    ssb<-((storeResN[g,ssbIndex, ] + storeStructN[g,ssbIndex, ]) * storeNums[g,ssbIndex, ] ) *mg_2_tonnes * X_CN
  }else{
    ssbIndex<-seq(1,thisNumCohorts)[seq(1,thisNumCohorts)>ageMat]
    ssb<-apply((storeResN[g,ssbIndex, ] + storeStructN[g,ssbIndex, ]) * storeNums[g,ssbIndex, ], 2, sum) *mg_2_tonnes * X_CN
  }
  
  plot(ssb, type="l", ylim=c(0,max(ssb)), xlab="Year", ylab="SSB (tonnes)",xaxt="n")
  axis(at=xLabsAt,labels=xLabs,side=1)
  bcur<-(ssb/ssb[1])*100
  thisMax<-100
  if(!is.na(sum(bcur))){thisMax<-max(bcur)}
  plot(bcur, type="l",ylim=c(0,thisMax), xlab="Year", ylab="Bcur(%B0)", xaxt="n")
  axis(at=xLabsAt,labels=xLabs,side=1)
  
  ##get selectivity from harvest file
  harvestFile<-paste(basePath,"CRAM_harvest_short.prm", sep="")
    harvestLines<-readLines(harvestFile)
    
    thisVar<-paste("CatchTS_agedistrib", thisCode, sep=""); temp<-harvestLines[grep(thisVar, harvestLines)+1]; thisCatchPropByAge<-get_first_number(temp,n="all")
      plot(x=seq(1,thisNumCohorts),y=thisCatchPropByAge[1:thisNumCohorts], type="h", lwd=5, lend=1, col=myBlue, ylab="Selectivity at age", xlab="Age class", xaxt="n")
      axis(at=seq(1,thisNumCohorts), labels=seq(1,thisNumCohorts), side=1)
      
      mtext(as.character(thisCode), side=3,outer=TRUE)
      dev.off()
  }
}  

######################################
## to SSB with catches and trawl survey indices
# plotPath<-paste(DIR$'Figures',"Validation\\",sep="") ## use this to overwrite the figure used in the report
lastModelYear<-2015
catchPath<-paste(this_path,"inputs\\catch_history\\",sep="")
#read in catch ARRAY
catch_array<-read.csv(paste(DIR$'Tables',"ALLcatchHistories.csv",sep=""))
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
  thisIsFished<-groupsDF$IsFished[g]
  if(thisIsFished==1){
    thisMax<-max(storeNums[g,,], na.rm=TRUE)
    # thisPlotFile<-paste(plotPath, "SSBvsOBS",thisCode,".jpg", sep="")
    # jpeg(thisPlotFile)
    thisPlotFile<-paste(plotPath,"EstObsAndTotalCatchbyYear",thisCode,sep="")
    pdf(paste(thisPlotFile,".pdf",sep=""),height=3.5, width=4.5)
    par(mfrow=c(1,1), mar=c(4,4,1,1), oma=c(1,1,1,1))
  
    ##get age mature
    thisVar<-paste(thisCode, "_age_mat", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageMat<-get_first_number(temp)
    if(thisNumCohorts==2){
      ssbIndex<-2
      ssb<-((storeResN[g,ssbIndex, ] + storeStructN[g,ssbIndex, ]) * storeNums[g,ssbIndex, ] ) *mg_2_tonnes * X_CN
    }else{
      ssbIndex<-seq(1,thisNumCohorts)[seq(1,thisNumCohorts)>ageMat]
      ssb<-apply((storeResN[g,ssbIndex, ] + storeStructN[g,ssbIndex, ]) * storeNums[g,ssbIndex, ], 2, sum) *mg_2_tonnes * X_CN
    }
    
    #the catch history
    tempCatchData<-catch_array[,c(as.character(thisCode))]
    catchYearIndex<-!is.na(tempCatchData); firstCatchYear<-min(seq(1900,lastModelYear)[catchYearIndex]); thisCatchYears<-seq(firstCatchYear,lastModelYear)
    thisCatchData<-tempCatchData[seq(1900,lastModelYear)>=firstCatchYear]
    #trawl survey
    thisObs<-rep(NA,length(thisCatchYears)); obsYears<-seq(1,length(thisCatchYears))
    TSfile<-paste(DIR$'Tables',"TrawlSurvey\\",thisCode,".csv",sep="")
    if(file.exists(TSfile)){
      thisTS<-read.csv(TSfile)
      thisObs<-thisTS$Biomass[thisTS$Year %in% thisCatchYears]; obsYears<-thisTS$Year[thisTS$Year %in% thisCatchYears]
      thisCVs<-thisTS$cv[thisTS$Year %in% thisCatchYears]/100; 
    }
    thisE<-ssb[seq(1900,lastModelYear) %in% obsYears] #this one used to scale obs, as overlap temporaly
    thisEst<-ssb[seq(1900,lastModelYear) %in% thisCatchYears]
    biomassYears<-xLabs
    biomassIndex<-biomassYears %in% thisCatchYears
    biomassMax<-max(thisEst)*1.5
    biomassAxis<-pretty(seq(0,biomassMax,length.out=5))
    
    thisYmax<-max(thisCatchData/1000, na.rm=TRUE)*1.1
    if(file.exists(TSfile)){
      if(sum(thisObs,na.rm=TRUE)>0){
        # thisScale<-mean(thisEst[match(biomassYears, obsYears)], na.rm=TRUE)/mean(thisObs,na.rm=TRUE)
        newScaledReal<-(thisObs/mean(thisObs,na.rm=TRUE))*mean(thisE, na.rm=TRUE)
        # newScaledReal<-(thisObs * thisScale)
        #get confidence intervals
        CIs<-getCIfromCVs(newScaledReal,thisCVs)
        # limitMax<-max(CIs$UpperCI, na.rm=TRUE)
        limitMax<-max(newScaledReal, na.rm=TRUE)
        if(limitMax>biomassMax){biomassMax<-limitMax}
      }
    }
    
    par(mar=c(3,4.5,2,4.5))
    par(las=0)
    # thisYmax<-90000
    plot(x=thisCatchYears,y=thisCatchData/1000,type="h",lwd=5,lend=1,xlab="",ylab="",col=catchBarColor,cex.axis=thisCex,ylim=c(0,thisYmax), xlim=c(firstCatchYear, lastYear))
    mtext(thisCode,side=3,adj=0,line=0.2,cex=thisCex)
    mtext("Catch (tonnes)",side=2,cex=thisCex,line=3)
    plotBiomassYears<-thisCatchYears[match(seq(1900,lastModelYear),thisCatchYears)]
    plotBiomassYears<-plotBiomassYears[!is.na(plotBiomassYears)]
    par(new=TRUE)
    plot(x=thisCatchYears,y=thisEst,type="l",col="black",cex=thisCex,ylim=c(0,biomassMax),
         xaxt="n",xlab="",lwd=3,yaxt="n",ylab="", xlim=c(firstCatchYear, lastYear), yaxt="n")
    thisBiomassAxis<-pretty(seq(0,biomassMax), length.out = 5)
    axis(at=thisBiomassAxis,labels=thisBiomassAxis,side=4,cex.axis=thisCex)
    mtext("SSB (tonnes)",side=4,adj=0.5,line=2.5,cex=thisCex)
    if(file.exists(TSfile)){
      if(sum(thisObs,na.rm=TRUE)>0){
        points(x=obsYears,y=newScaledReal,pch=20,col="red",cex=1.2)
        segments(y0=CIs$LowerCI, x0=obsYears, y1=CIs$UpperCI, x1=obsYears, col="red")
      }
    }
    dev.off()
  }
}
    

    



