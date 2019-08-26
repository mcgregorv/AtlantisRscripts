#bring in the run, the observed data, the catch history and plot together
source(paste(DIR$'General functions',"getCIfromCV.R", sep=""))
source(paste(DIR$'General functions',"getCIfromCVs.R", sep=""))
nboxes<-30

mg_2_tonne<-2e-8; X_CN<-5.7

catchYears<-seq(1900,2014); ny<-length(catchYears)

#################################### TRACERS
#plot all tracers for a given box and layer
this_run<-"base"

# base_out<-"ORHmum4base"
# this_out<-"ORHmum4fish"

base_out<-"BASE4"
this_out<-"FISH4"
# 
# this_out<-"TEST150yrfish"
# base_out<-"TEST150yrbase"

burnin<-35 #number of years to skip in plot

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")
baseOutPath<-paste(this_path,"output",base_out,"\\", sep="")
r<-""


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

###########################################
##plot them all seperately too
###################################
topCatchYears<-seq(1950,2014); biomassYears<-topCatchYears
topCatchYears<-seq(1900,2014); biomassYears<-topCatchYears

catchBarColor<-myGrey_trans
##do individual plot
g=17
g=25
for(g in 1:ncg){
  thisCode<-catchCodes[g]
  cat(as.character(thisCode))
  thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode],side="both")
  thisVar<-paste(thisName,"_N",sep="")
  tempData<-ncvar_get(ThisNC.nc,thisVar)
  thisEst<-apply(tempData*thisVol,3,sum)*mg_2_tonne*X_CN; 
  thisEst<-thisEst[burnin:length(thisEst)]
  #the catch history
  thisData<-catch_array[catchYears %in% topCatchYears,c(as.character(thisCode))]
  thisData[is.na(thisData)]<-0
  #trawl survey
  thisObs<-rep(NA,length(topCatchYears)); obsYears<-seq(1,length(topCatchYears))
  TSfile<-paste(DIR$'Tables',"TrawlSurvey\\",thisCode,".csv",sep="")
  if(file.exists(TSfile)){
    thisTS<-read.csv(TSfile)
    thisObs<-thisTS$Biomass[thisTS$Year %in% topCatchYears]; obsYears<-thisTS$Year[thisTS$Year %in% topCatchYears]
    thisCVs<-thisTS$cv[thisTS$Year %in% topCatchYears]/100; 
  }
  thisE<-thisEst[biomassYears %in% obsYears]
  biomassYears<-xLabs
  biomassIndex<-biomassYears %in% topCatchYears
  biomassMax<-max(thisEst)*1.5
  biomassAxis<-pretty(seq(0,biomassMax,length.out=5))
  
  thisB0<-B0data$B0[B0data$Code==thisCode]
  
  # newBiomassAxis<-(biomassAxis/max(thisEst,na.rm=TRUE))*max(thisData/1000)
  
  thisYmax<-max(thisData/1000)*1.1
  if(file.exists(TSfile)){
    if(sum(thisObs,na.rm=TRUE)>0){
      # thisScale<-mean(thisEst[match(biomassYears, obsYears)], na.rm=TRUE)/mean(thisObs,na.rm=TRUE)
      newScaledReal<-(thisObs/mean(thisObs,na.rm=TRUE))*mean(thisEst[biomassYears %in% obsYears])
      # newScaledReal<-(thisObs * thisScale)
      #get confidence intervals
      CIs<-getCIfromCVs(newScaledReal,thisCVs)
      # limitMax<-max(CIs$UpperCI, na.rm=TRUE)
      limitMax<-max(newScaledReal, na.rm=TRUE)
      if(limitMax>biomassMax){biomassMax<-limitMax}
    }
  }
  
  firstCatchYear<-min(topCatchYears[thisData>0]); lastYear<-max(topCatchYears)
  
  thisPlotFile<-paste(plotPath,"EstObsAndTotalCatchbyYear",thisCode,sep="")
  pdf(paste(thisPlotFile,".pdf",sep=""),height=3.5, width=4.5)
  # jpeg(paste(thisPlotFile, ".jpg", sep=""), quality=3000)
  par(mar=c(3,4.5,2,4.5))
  par(las=0)
  # thisYmax<-90000
  plot(x=topCatchYears,y=thisData/1000,type="h",lwd=5,lend=1,xlab="",ylab="",col=catchBarColor,cex.axis=thisCex,ylim=c(0,thisYmax), xlim=c(firstCatchYear, lastYear))
  mtext(thisCode,side=3,adj=0,line=0.2,cex=thisCex)
  mtext("Catch (tonnes)",side=2,cex=thisCex,line=3)
  plotBiomassYears<-topCatchYears[match(biomassYears,topCatchYears)]
  plotBiomassYears<-plotBiomassYears[!is.na(plotBiomassYears)]
  par(new=TRUE)
  plot(x=plotBiomassYears,y=(thisEst[biomassIndex]),type="l",col="black",cex=thisCex,ylim=c(0,biomassMax),
       xaxt="n",xlab="",lwd=3,yaxt="n",ylab="", xlim=c(firstCatchYear, lastYear), yaxt="n")
  thisBiomassAxis<-pretty(seq(0,biomassMax), length.out = 5)
  axis(at=thisBiomassAxis,labels=thisBiomassAxis,side=4,cex.axis=thisCex)
  mtext("Biomass (tonnes)",side=4,adj=0.5,line=2.5,cex=thisCex)
  if(file.exists(TSfile)){
    if(sum(thisObs,na.rm=TRUE)>0){
      points(x=obsYears,y=newScaledReal,pch=20,col="red",cex=1.2)
      segments(y0=CIs$LowerCI, x0=obsYears, y1=CIs$UpperCI, x1=obsYears, col="red")
    }
  }
  # legend(legend=c("Catches (forced)", "Model estimated biomass", "Trawl survey estimated biomass"), col=c(catchBarColor, "black", "red"), lwd=c(4,3,1.5),
         # seg.len=2.5,  lty=1, pch=c(NA, NA, 20), x="topleft", bty="n")
  
  dev.off()
  
}
## do legend seperately
thisPlotFile<-paste(plotPath,"EstObsAndTotalCatchbyYearLEGEND",sep="")
pdf(paste(thisPlotFile,".pdf",sep=""),height=2, width=3)
# jpeg(paste(thisPlotFile, ".jpg", sep=""), quality=3000)
par(mar=c(0,0,0,0), lend=1)
makeBlankPlot()
legend(legend=c("Catches (forced)", "Model estimated biomass", "Trawl survey estimated biomass"), col=c(catchBarColor, "black", "red"), lwd=c(4,3,1.5), 
       seg.len=2.5,  lty=1, pch=c(NA, NA, 20), x="center", bty="n")
dev.off()

## create tex insert
orderByLName<-match(sort(groupsDF$LongName), groupsDF$LongName)
codesByLName<-groupsDF$Code[orderByLName]

count<-1; nperpage<-15
skipGroups<-c("CRA", "IVH", "IVS", "PFL","PFS") ## no trawl survey
# skipGroups<-c()
texFile<-paste(DIR$'Reports',"(01)BaseReport\\ObsVsEst_figures.tex", sep="")
cat("", file=texFile, append=FALSE)
thisCaption<-paste("Observed biomass estimated from trawl surveys (red), estimated biomass from CRAM (black) and forced catch history (grey) for all groups with trawl survey estimates.",sep="")
thisLab<-"ObsVsEst"
thisFigText<-paste("\\begin{figure}[H]
\\centering", sep="")
cat(thisFigText, file=texFile, append=TRUE)
for(g in 1:ng){
  thisCode<-codesByLName[g];  thisFig<-paste("EstObsAndTotalCatchbyYear", thisCode,".pdf", sep="")
  cat(as.character(thisCode), count, "\n")
  if(!(thisCode %in% skipGroups) & thisCode %in% catchCodes){
    if(count==(nperpage+1)){
      #start a new figure
      texFile<-paste(DIR$'Reports',"(01)BaseReport\\ObsVsEst_figuresPart2.tex", sep="")
      cat("", file=texFile, append=FALSE)
      thisCaption<-paste("Observed biomass estimated from trawl surveys (red), estimated biomass from CRAM (black) and forced catch history (grey) for all groups with trawl survey estimates.",sep="")
      thisLab<-"ObsVsEstP2"
      thisFigText<-paste("\\begin{figure}[H]
\\centering", sep="")
      cat(thisFigText, file=texFile, append=TRUE)
    }
    thisFigText<-paste("\\includegraphics[width=5cm]{", thisFig,"}\n", sep="")
    cat(thisFigText, file=texFile, append=TRUE)
    if(count==nperpage){
      #close off previous figure
      thisFigText<-paste("\\caption{",thisCaption,"}\\label{gif:", thisLab,"}
                         \\end{figure}\n",sep="")
      cat(thisFigText, file=texFile, append=TRUE)
    }
    count<-count+1
  }
}
thisFigText<-paste("\\caption{",thisCaption,"}\\label{gif:", thisLab,"}
\\end{figure}\n",sep="")
cat(thisFigText, file=texFile, append=TRUE)


## if done individually
# g=1
# thisCode<-catchCodes[g]; thisFig<-paste("EstObsAndTotalCatchbyYear", thisCode,".pdf", sep="")
# thisCaption<-paste("Observed biomass estimated from trawl surveys (red), estimated biomass from CRAM (black) and forced catch history (grey) for ", thisCode,".",sep="")
# thisLab<-paste("ObsVsEst_",thisCode,sep="")
# thisFigText<-paste("\\begin{figure}[H]
# \\centering
# \\includegraphics[width=11cm]{", thisFig,"}
# \\caption{",thisCaption,"}\\label{gif:", thisLab,"}
# \\end{figure}\n",sep="")
# cat(thisFigText, file=texFile, append=TRUE)
###############################################
## make each plot again, but from 1900 and with base run behind it

for(g in 1:ncg){
  thisCode<-catchCodes[g]
  cat(as.character(thisCode))
  thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode],side="both")
  thisVar<-paste(thisName,"_N",sep="")
  tempData<-ncvar_get(ThisNC.nc,thisVar)
  thisEst<-apply(tempData*thisVol,3,sum)*mg_2_tonne*X_CN; 
  thisEst<-thisEst[burnin:length(thisEst)]
  biomassAxis<-pretty(seq(0,max(thisEst),length.out=5))
  ## get the base biomass from model
  tempData<-ncvar_get(BaseNC.nc,thisVar)
  baseEst<-apply(tempData*baseVol,3,sum)*mg_2_tonne*X_CN; 
  baseEst<-baseEst[burnin:length(baseEst)]
 
  #the catch history
  thisData<-catch_array[catchYears %in% topCatchYears,c(as.character(thisCode))]
  thisData[is.na(thisData)]<-0
  #trawl survey
  thisObs<-rep(NA,length(topCatchYears)); obsYears<-seq(1,length(topCatchYears))
  TSfile<-paste(DIR$'Tables',"TrawlSurvey\\",thisCode,".csv",sep="")
  if(file.exists(TSfile)){
    thisTS<-read.csv(TSfile)
    thisObs<-thisTS$Biomass; obsYears<-thisTS$Year
    CIs<-getCIfromCVs(newScaledReal,thisCVs)
    # thisObs<-thisTS$Biomass[thisTS$Year %in% topCatchYears]; obsYears<-thisTS$Year[thisTS$Year %in% topCatchYears]
  }
  # thisE<-thisEst[biomassYears %in% obsYears]
  thisE<-thisEst
  biomassYears<-xLabs
  biomassIndex<-biomassYears %in% topCatchYears
  
  thisB0<-B0data$B0[B0data$Code==thisCode]
  
  newBiomassAxis<-(biomassAxis/max(thisEst,na.rm=TRUE))*max(thisData/1000)
  
  thisYmax<-max(max(thisData/1000),newBiomassAxis,(thisB0/max(thisEst))*max(thisData/1000))
  
  firstCatchYear<-min(topCatchYears[thisData>0]); lastYear<-max(topCatchYears)
  
  thisPlotFile<-paste(DIR$'Figures',"\\Validation\\testing\\",this_out,"EstObsAndTotalCatchbyYear_",thisCode,"_fullTimeSeries",sep="")
  # pdf(paste(thisPlotFile,".pdf",sep=""),height=4)
  jpeg(paste(thisPlotFile, ".jpg", sep=""), quality=3000)
  par(mar=c(3,4.5,2,4.5))
  par(las=0)
  # thisYmax<-90000
  plot(x=topCatchYears,y=thisData/1000,type="h",lwd=5,lend=1,xlab="",ylab="",col=catchBarColor,cex.axis=thisCex,ylim=c(0,thisYmax))
  mtext(thisCode,side=3,adj=0,line=0.2,cex=thisCex)
  mtext("Catch (tonnes)",side=2,cex=thisCex,line=3)
  plotBiomassYears<-topCatchYears[match(biomassYears,topCatchYears)]
  plotBiomassYears<-plotBiomassYears[!is.na(plotBiomassYears)]
  # plotBiomassYears<-biomassYears
  biomassMax<-max(max(thisEst), max(baseEst))*1.5
  if(file.exists(TSfile)){
    if(sum(thisObs,na.rm=TRUE)>0){
      newScaledReal<-(thisObs/mean(thisObs,na.rm=TRUE))*mean(thisEst[biomassYears %in% obsYears])
      obsMax<-max(newScaledReal, na.rm=TRUE)
      if(obsMax>biomassMax){biomassMax<-obsMax*1.2}
    }
  }
  
  par(new=TRUE, lend=1)
  plot(x=plotBiomassYears, y=baseEst[biomassIndex][1:length(thisEst[biomassIndex])], type="l", lwd=2, lty=2, col=myGrey, cex=thisCex,ylim=c(0,biomassMax),  xaxt="n",xlab="",yaxt="n",ylab="", yaxt="n")
  points(x=plotBiomassYears,y=thisEst[biomassIndex],type="l",col="black",lwd=2)
  thisBiomassAxis<-pretty(seq(0,max(thisEst)*1.5, length.out = 5))
  axis(at=thisBiomassAxis,labels=thisBiomassAxis,side=4,cex.axis=thisCex)
  mtext("Biomass (tonnes)",side=4,adj=0.5,line=2.5,cex=thisCex)
  if(file.exists(TSfile)){
    if(sum(thisObs,na.rm=TRUE)>0){
       #get confidence intervals
      CIs<-getCIfromCVs(newScaledReal,thisCVs)
      points(x=obsYears,y=newScaledReal,pch=20,col="red", cex=1.2)
      segments(y0=CIs$LowerCI, x0=obsYears, y1=CIs$UpperCI, x1=obsYears, col="red")
      
      # points(x=obsYears,y=newScaledReal,type="l",col="red")
      # polygon(x=c(obsYears, rev(obsYears)), y=c(CIs$LowerCI, rev(CIs$UpperCI)), col=myRed_trans, border=NA)
    }
  }
  # legend(legend=c("Catches (forced)", "Atlantis model", "Trawl survey"), col=c(myBlue, "black", "red"), lwd=c(4,3,1.5), 
         # seg.len=2.5,  lty=1, pch=c(NA, NA, 20), x="topleft", bty="n")
  
  legend(legend=c("Catches (forced)", "Atlantis model", "Atlantis no-fishing", "Trawl survey"), col=c(catchBarColor, "black", myGrey, "red"), lwd=c(4,2,2,NA), 
         seg.len=2.5,  pch=c(NA,NA, NA, 20), x="topleft", bty="n", lty=c(1,1,2,NA))
  
  
  dev.off()
  
}

##############################################
##do again for RFI with focus on trawl survey years
g<-29
thisCode<-catchCodes[g]
cat(as.character(thisCode))
thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode],side="both")
thisVar<-paste(thisName,"_N",sep="")
tempData<-ncvar_get(ThisNC.nc,thisVar)
thisEst<-apply(tempData*thisVol,3,sum)*mg_2_tonne*X_CN; 
thisEst<-thisEst[burnin:length(thisEst)]
#trawl survey
thisObs<-rep(NA,length(topCatchYears)); obsYears<-seq(1,length(topCatchYears))
TSfile<-paste(DIR$'Tables',"TrawlSurvey\\",thisCode,".csv",sep="")
if(file.exists(TSfile)){
  thisTS<-read.csv(TSfile)
  thisObs<-thisTS$Biomass[thisTS$Year %in% topCatchYears]; obsYears<-thisTS$Year[thisTS$Year %in% topCatchYears]
}
thisE<-thisEst[biomassYears %in% obsYears]
topCatchYears<-obsYears
biomassYears<-xLabs
biomassIndex<-biomassYears %in% topCatchYears
biomassMax<-max(thisEst)*5
biomassAxis<-pretty(seq(0,biomassMax,length.out=5))

thisB0<-B0data$B0[B0data$Code==thisCode]

#the catch history
thisData<-catch_array[catchYears %in% topCatchYears,c(as.character(thisCode))]
thisData[is.na(thisData)]<-0

thisYmax<-max(thisData/1000)*1.1

firstCatchYear<-min(obsYears); lastYear<-max(obsYears)

thisPlotFile<-paste(DIR$'Figures',"\\Validation\\testing\\",this_out,"EstObsAndTotalCatchbyYear_",thisCode, "zoomedIn",sep="")
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
plot(x=plotBiomassYears,y=(thisEst[biomassIndex]),type="l",col="black",cex=thisCex,ylim=c(0,biomassMax),
     xaxt="n",xlab="",lwd=3,yaxt="n",ylab="", xlim=c(firstCatchYear, lastYear), yaxt="n")
thisBiomassAxis<-pretty(seq(0,biomassMax), length.out = 5)
axis(at=thisBiomassAxis,labels=thisBiomassAxis,side=4,cex.axis=thisCex)
mtext("Biomass (tonnes)",side=4,adj=0.5,line=2.5,cex=thisCex)
if(file.exists(TSfile)){
  if(sum(thisObs,na.rm=TRUE)>0){
    thisScale<-mean(thisEst[match(biomassYears, obsYears)], na.rm=TRUE)/mean(thisObs,na.rm=TRUE)
    newScaledReal<-(thisObs * thisScale)
    #get confidence intervals
    CIs<-getCIfromCV(newScaledReal,0.3)
    points(x=obsYears,y=newScaledReal,pch=20,col="red", cex=1.5)
    points(x=obsYears,y=newScaledReal,type="l",col="red")
    # polygon(x=c(obsYears, rev(obsYears)), y=c(CIs$LowerCI, rev(CIs$UpperCI)), col=myRed_trans, border=NA)
  }
}
legend(legend=c("Catches (forced)", "Model estimated biomass", "Trawl survey estimated biomass"), col=c(myBlue, "black", "red"), lwd=c(4,3,1.5), 
       seg.len=2.5,  lty=1, pch=c(NA, NA, 20), x="topleft", bty="n")

dev.off()

#
##################################
#

# pdf(paste(plotPath,"EstObsAndTotalCatchbyYear_ALLcodes.pdf",sep=""))
# par(mar=c(3,4.5,2,4.5),mfrow=c(3,1))
# for(g in 1:ncg){
#   thisCode<-catchCodes[g]
#   cat(as.character(thisCode))
#   thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode],side="both")
#   thisVar<-paste(thisName,"_N",sep="")
#   tempData<-ncvar_get(ThisNC.nc,thisVar)
#   thisEst<-apply(tempData*thisVol,3,sum)*mg_2_tonne*X_CN; 
#   thisEst<-thisEst[burnin:length(thisEst)]
#   biomassAxis<-pretty(seq(0,max(thisEst),length.out=5))
#   #the catch history
#   thisData<-catch_array[catchYears %in% topCatchYears,c(as.character(thisCode))]
#   thisData[is.na(thisData)]<-0
#   #trawl survey
#   thisObs<-rep(NA,length(topCatchYears)); obsYears<-seq(1,length(topCatchYears))
#   TSfile<-paste(DIR$'Tables',"TrawlSurvey\\",thisCode,".csv",sep="")
#   if(file.exists(TSfile)){
#     thisTS<-read.csv(TSfile, sep="\t")
#     thisObs<-thisTS$Biomass[thisTS$Year %in% topCatchYears]; obsYears<-thisTS$Year[thisTS$Year %in% topCatchYears]
#   }
#   thisE<-thisEst[biomassYears %in% obsYears]
#   biomassYears<-xLabs
#   biomassIndex<-biomassYears %in% topCatchYears
#   
#   thisB0<-B0data$B0[B0data$Code==thisCode]
#   
#   newBiomassAxis<-(biomassAxis/max(thisEst,na.rm=TRUE))*max(thisData/1000)
#   
#   thisYmax<-max(max(thisData/1000),newBiomassAxis,(thisB0/max(thisEst))*max(thisData/1000))
#   
#   # pdf(paste(plotPath,"EstObsAndTotalCatchbyYear_",thisCode,".pdf",sep=""),height=4)
#   # par(mar=c(3,4.5,2,4.5))
#   par(las=0)
#   plot(x=topCatchYears,y=thisData/1000,type="h",lwd=5,lend=1,xlab="",ylab="",col=myBlue,cex.axis=thisCex,ylim=c(0,thisYmax))
#   mtext(thisCode,side=3,adj=0,line=0.2,cex=thisCex)
#   mtext("Catch (tonnes)",side=2,cex=thisCex,line=3)
#   
#   plotBiomassYears<-topCatchYears[match(biomassYears,topCatchYears)]
#   plotBiomassYears<-plotBiomassYears[!is.na(plotBiomassYears)]
#   
#   points(x=plotBiomassYears,y=(thisEst[biomassIndex]/max(thisEst,na.rm=TRUE))*max(thisData/1000),type="l",col="black",cex=thisCex,ylim=c(0,max(thisEst)),xaxt="n",xlab="",lwd=3,yaxt="n",ylab="")
#   
#   
#   # par(new=TRUE) #this means the next plot will overlay. But got to be carefull the x-axes line up
#   # plot(x=plotBiomassYears,y=thisEst[biomassIndex],type="l",col="black",cex=thisCex,ylim=c(0,max(thisEst)),xaxt="n",xlab="",lwd=3,yaxt="n",ylab="")
#   axis(at=newBiomassAxis,labels=biomassAxis,side=4,cex.axis=thisCex)
#   mtext("Biomass (tonnes)",side=4,adj=0.5,line=2.5,cex=thisCex)
#   abline(h=(thisB0/max(thisEst,na.rm=TRUE))*max(thisData/1000),col=myRed,lwd=2,lty=2)
#   if(file.exists(TSfile)){
#     # points(x=thisTS$Year,y=thisTS$Biomass,pch=20)
#     
#     if(sum(thisObs,na.rm=TRUE)>0){
#       ##scale biomass
#       # lmfit=nls(thisE~b+a*thisObs,algorithm="port",start=c(a=1,b=0),lower=c(a=0,b=0),upper=c(a=Inf,b=0))
#       # ScaledReal = abs(coef(lmfit)[1]) * thisObs[!is.na(thisObs)] + coef(lmfit)[2]
#       # 
#       # newScaledReal<-(ScaledReal/max(thisEst,na.rm=TRUE))*max(thisData/1000)
#       temp<-(thisObs/mean(thisObs,na.rm=TRUE))*mean(thisEst)
#       
#       newScaledReal<-(temp/max(temp,na.rm=TRUE))*max(thisData/1000)
#       
#       points(x=obsYears,y=newScaledReal,pch=20,col="red")
#       points(x=obsYears,y=newScaledReal,type="l",col="red")
#     }
#   }
#   par(las=1)
#   axis(at=(thisB0/max(thisEst,na.rm=TRUE))*max(thisData/1000),labels = expression(B[0]),side=4,col=myRed,col.axis=myRed,cex.axis=thisCex,tick=FALSE,line=-1)
#   
#   # dev.off()
# }
# dev.off()



