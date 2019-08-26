#bring in the run, the observed data, the catch history and plot together

nboxes<-30

mg_2_tonne<-2e-8; X_CN<-5.7

catchYears<-seq(1900,2014); ny<-length(catchYears)

#################################### TRACERS
#plot all tracers for a given box and layer
this_run<-"base"

this_out<-"SENselALL2"
# this_out<-"SENSselBASEfish"

burnin<-35 #number of years to skip in plot

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")
r<-""


#read in groups csv file to get codes for fished groups
groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
catchCodes<-groupsDF$Code[groupsDF$IsFished==1]; ncg<-length(catchCodes)

catchPath<-paste(this_path,"inputs\\catch_history\\",sep="")
#read in catch ARRAY
catch_array<-read.csv(paste(DIR$'Tables',"ALLcatchHistories.csv",sep=""))


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

B0data<-read.csv(paste(this_path,"..\\CRAM_B0.csv",sep=""))

###########################################
##plot them all seperately too
###################################
topCatchYears<-seq(1950,2014); biomassYears<-topCatchYears
topCatchYears<-seq(1900,2014); biomassYears<-topCatchYears

pdf(paste(plotPath,"EstObsAndTotalCatchbyYear_ALLcodes.pdf",sep=""))
par(mar=c(3,4.5,2,4.5),mfrow=c(3,1))
for(g in 1:ncg){
  thisCode<-catchCodes[g]
  cat(as.character(thisCode))
  thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode],side="both")
  thisVar<-paste(thisName,"_N",sep="")
  tempData<-ncvar_get(ThisNC.nc,thisVar)
  thisEst<-apply(tempData*thisVol,3,sum)*mg_2_tonne*X_CN; 
  thisEst<-thisEst[burnin:length(thisEst)]
  biomassAxis<-pretty(seq(0,max(thisEst),length.out=5))
  #the catch history
  thisData<-catch_array[catchYears %in% topCatchYears,c(as.character(thisCode))]
  thisData[is.na(thisData)]<-0
  #trawl survey
  thisObs<-rep(NA,length(topCatchYears)); obsYears<-seq(1,length(topCatchYears))
  TSfile<-paste(DIR$'Tables',"TrawlSurvey\\",thisCode,".csv",sep="")
  if(file.exists(TSfile)){
    thisTS<-read.csv(TSfile, sep="\t")
    thisObs<-thisTS$Biomass[thisTS$Year %in% topCatchYears]; obsYears<-thisTS$Year[thisTS$Year %in% topCatchYears]
  }
  thisE<-thisEst[biomassYears %in% obsYears]
  biomassYears<-xLabs
  biomassIndex<-biomassYears %in% topCatchYears
  
  thisB0<-B0data$B0[B0data$Code==thisCode]
  
  newBiomassAxis<-(biomassAxis/max(thisEst,na.rm=TRUE))*max(thisData/1000)
  
  thisYmax<-max(max(thisData/1000),newBiomassAxis,(thisB0/max(thisEst))*max(thisData/1000))
  
  # pdf(paste(plotPath,"EstObsAndTotalCatchbyYear_",thisCode,".pdf",sep=""),height=4)
  # par(mar=c(3,4.5,2,4.5))
  par(las=0)
  plot(x=topCatchYears,y=thisData/1000,type="h",lwd=5,lend=1,xlab="",ylab="",col=myBlue,cex.axis=thisCex,ylim=c(0,thisYmax))
  mtext(thisCode,side=3,adj=0,line=0.2,cex=thisCex)
  mtext("Catch (tonnes)",side=2,cex=thisCex,line=3)
  
  plotBiomassYears<-topCatchYears[match(biomassYears,topCatchYears)]
  plotBiomassYears<-plotBiomassYears[!is.na(plotBiomassYears)]
  
  points(x=plotBiomassYears,y=(thisEst[biomassIndex]/max(thisEst,na.rm=TRUE))*max(thisData/1000),type="l",col="black",cex=thisCex,ylim=c(0,max(thisEst)),xaxt="n",xlab="",lwd=3,yaxt="n",ylab="")
  
  
  # par(new=TRUE) #this means the next plot will overlay. But got to be carefull the x-axes line up
  # plot(x=plotBiomassYears,y=thisEst[biomassIndex],type="l",col="black",cex=thisCex,ylim=c(0,max(thisEst)),xaxt="n",xlab="",lwd=3,yaxt="n",ylab="")
  axis(at=newBiomassAxis,labels=biomassAxis,side=4,cex.axis=thisCex)
  mtext("Biomass (tonnes)",side=4,adj=0.5,line=2.5,cex=thisCex)
  abline(h=(thisB0/max(thisEst,na.rm=TRUE))*max(thisData/1000),col=myRed,lwd=2,lty=2)
  if(file.exists(TSfile)){
    # points(x=thisTS$Year,y=thisTS$Biomass,pch=20)
    
    if(sum(thisObs,na.rm=TRUE)>0){
      ##scale biomass
      # lmfit=nls(thisE~b+a*thisObs,algorithm="port",start=c(a=1,b=0),lower=c(a=0,b=0),upper=c(a=Inf,b=0))
      # ScaledReal = abs(coef(lmfit)[1]) * thisObs[!is.na(thisObs)] + coef(lmfit)[2]
      # 
      # newScaledReal<-(ScaledReal/max(thisEst,na.rm=TRUE))*max(thisData/1000)
      temp<-(thisObs/mean(thisObs,na.rm=TRUE))*mean(thisEst)
      
      newScaledReal<-(temp/max(temp,na.rm=TRUE))*max(thisData/1000)
      
      points(x=obsYears,y=newScaledReal,pch=20,col="red")
      points(x=obsYears,y=newScaledReal,type="l",col="red")
    }
  }
  par(las=1)
  axis(at=(thisB0/max(thisEst,na.rm=TRUE))*max(thisData/1000),labels = expression(B[0]),side=4,col=myRed,col.axis=myRed,cex.axis=thisCex,tick=FALSE,line=-1)
  
  # dev.off()
}
dev.off()

##do individual plot
g=17
g=5
for(g in 1:ncg){
  thisCode<-catchCodes[g]
  cat(as.character(thisCode))
  thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode],side="both")
  thisVar<-paste(thisName,"_N",sep="")
  tempData<-ncvar_get(ThisNC.nc,thisVar)
  thisEst<-apply(tempData*thisVol,3,sum)*mg_2_tonne*X_CN; 
  thisEst<-thisEst[burnin:length(thisEst)]
  biomassAxis<-pretty(seq(0,max(thisEst),length.out=5))
  #the catch history
  thisData<-catch_array[catchYears %in% topCatchYears,c(as.character(thisCode))]
  thisData[is.na(thisData)]<-0
  #trawl survey
  thisObs<-rep(NA,length(topCatchYears)); obsYears<-seq(1,length(topCatchYears))
  TSfile<-paste(DIR$'Tables',"TrawlSurvey\\",thisCode,".csv",sep="")
  if(file.exists(TSfile)){
    thisTS<-read.csv(TSfile)
    thisObs<-thisTS$Biomass[thisTS$Year %in% topCatchYears]; obsYears<-thisTS$Year[thisTS$Year %in% topCatchYears]
  }
  thisE<-thisEst[biomassYears %in% obsYears]
  biomassYears<-xLabs
  biomassIndex<-biomassYears %in% topCatchYears
  
  thisB0<-B0data$B0[B0data$Code==thisCode]
  
  newBiomassAxis<-(biomassAxis/max(thisEst,na.rm=TRUE))*max(thisData/1000)
  
  thisYmax<-max(max(thisData/1000),newBiomassAxis,(thisB0/max(thisEst))*max(thisData/1000))
  
  firstCatchYear<-min(topCatchYears[thisData>0]); lastYear<-max(topCatchYears)
  
  thisPlotFile<-paste(DIR$'Figures',"\\Validation\\EstObsAndTotalCatchbyYear_",thisCode,sep="")
  # pdf(paste(thisPlotFile,".pdf",sep=""),height=4)
  jpeg(paste(thisPlotFile, ".jpg", sep=""), quality=3000)
  par(mar=c(3,4.5,2,4.5))
  par(las=0)
  # thisYmax<-90000
  plot(x=topCatchYears,y=thisData/1000,type="h",lwd=5,lend=1,xlab="",ylab="",col=myBlue,cex.axis=thisCex,ylim=c(0,thisYmax), xlim=c(firstCatchYear, lastYear))
  mtext(thisCode,side=3,adj=0,line=0.2,cex=thisCex)
  mtext("Catch (tonnes)",side=2,cex=thisCex,line=3)
  plotBiomassYears<-topCatchYears[match(biomassYears,topCatchYears)]
  plotBiomassYears<-plotBiomassYears[!is.na(plotBiomassYears)]
  par(new=TRUE)
  plot(x=plotBiomassYears,y=(thisEst[biomassIndex]),type="l",col="black",cex=thisCex,ylim=c(0,max(thisEst)*1.5),
       xaxt="n",xlab="",lwd=3,yaxt="n",ylab="", xlim=c(firstCatchYear, lastYear), yaxt="n")
  thisBiomassAxis<-pretty(seq(0,max(thisEst)*1.5, length.out = 5))
  axis(at=thisBiomassAxis,labels=thisBiomassAxis,side=4,cex.axis=thisCex)
  mtext("Biomass (tonnes)",side=4,adj=0.5,line=2.5,cex=thisCex)
  if(file.exists(TSfile)){
    if(sum(thisObs,na.rm=TRUE)>0){
      newScaledReal<-(thisObs/max(thisObs,na.rm=TRUE))*max(thisEst)
      #get confidence intervals
      CIs<-getCIfromCV(newScaledReal,0.3)
      points(x=obsYears,y=newScaledReal,pch=20,col="red")
      points(x=obsYears,y=newScaledReal,type="l",col="red")
      polygon(x=c(obsYears, rev(obsYears)), y=c(CIs$LowerCI, rev(CIs$UpperCI)), col=myRed_trans, border=NA)
    }
  }
  legend(legend=c("Catches (forced)", "Model estimated biomass", "Trawl survey estimated biomass"), col=c(myBlue, "black", "red"), lwd=c(4,3,1.5), 
         seg.len=2.5,  lty=1, pch=c(NA, NA, 20), x="topleft", bty="n")
  
  dev.off()
  
}
################################################


