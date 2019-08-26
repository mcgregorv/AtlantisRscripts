source(paste(DIR$'General functions',"getCIfromCVs.R", sep=""))

#plot all tracers for a given box and layer
this_run<-"base"
# this_run<-"ClimateChange"
this_run<-"Chaos"
this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
biolLines<-readLines(paste(this_path,"..\\CRAM_BH_hybrid_biol.prm", sep=""))

this_out <- paste("A",c("ChaosFISHUp05", "ChaosFISHUp02", "ChaosFISHUp01", "ChaosFISHUp005","ChaosFISH", "ChaosFISHDown005", "ChaosFISHDown01", "ChaosFISHDown02", "ChaosFISHDown05"), sep="");
plotDescrip<-"ChaosUpDownFISH"

this_out <- paste("D",c("ChaosFishUp05", "ChaosFishUp02", "ChaosFishUp01", "ChaosFishUp005","ChaosFish", "ChaosFishDown005", "ChaosFishDown01", "ChaosFishDown02", "ChaosFishDown05"), sep="");
plotDescrip<-"DChaosUpDownFish"

this_out <- paste("ChaosFISHSampleA", 1:10, sep="");
plotDescrip<-"ChaosFISHSampleA"


mg_2_tonne<-2e-8; X_CN<-5.7
catchYears<-seq(1900,2014); ny<-length(catchYears)
surveyOK<-c("BIS", "CBO", "EIS", "ELI", "HAK", "HOK", "JAV", "LDO", "LIN",  "SPE")

nlayers<-6
groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
catchCodes<-groupsDF$Code[groupsDF$IsFished==1]; ncg<-length(catchCodes)
groupsDFPaper<-read.csv(paste(this_path,"..\\CRAM_groupsPaper.csv", sep=""))

plotPath<-paste(this_path,"..\\Figures\\chaos\\", plotDescrip,sep="")


catchPath<-paste(this_path,"..\\inputs\\catch_history\\",sep="")
catch_array<-read.csv(paste(DIR$'Tables',"ALLcatchHistories.csv",sep=""))

nruns<-length(this_out)
burnin<-rep(1,nruns) #number of years to skip in plot

runCols <- c(colorRampPalette(colors=c(myBlue,myAqua,myGreen))(4), "black", colorRampPalette(colors=c(myYellow, myOrange,"red"))(4))
# runCols <- c( "black",colorRampPalette(colors=c(myBlue,myAqua,myGreen))(4))

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1865
fishingStartYear<-1865
modelStartYear<-1865

nc_list<-NULL; nts_list<-NULL; min_nts<-1e+12
for(r in 1:nruns){
  outPath<-paste(this_path,"output",this_out[r],"\\",sep="")
  nc_list[[r]]<-nc_open(paste(outPath,"output.nc",sep=""))
  thisVol<-ncvar_get(nc_list[[r]],"volume")
  thisDz<-ncvar_get(nc_list[[r]],"dz")
  nts_list[[r]]<-dim(thisVol)[3]-burnin[r] #number of timesteps
  if(showall==TRUE){nts_list[[r]]<-dim(thisVol)[3]}
  if(nts_list[[r]]<min_nts){min_nts<-nts_list[[r]]}
}

nts_list
max_nts<-max(nts_list, na.rm=TRUE)

timeList<-NULL; timeMin <- 30000; timeMax <- 0
for(r in 1:nruns){
  this_nts<-nts_list[[r]]; this_burnin <- burnin[r]
  thisYear0<-1865 - this_burnin + 1
  this_time <- thisYear0 : (thisYear0 + this_nts -1)
  timeList[[r]]<-this_time
  if(max(this_time) > timeMax){timeMax<-max(this_time)}
  if(min(this_time) < timeMin){timeMin <- min(this_time)}
}
xLabsTemp<-seq(0,(max_nts*daysTimeStep-1),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0+burnin[1]

allTracers<-sort(names(ThisNC.nc$var))

B0data<-read.csv(paste(this_path,"..\\CRAM_B0.csv",sep=""))

topCatchYears<-seq(1900,2014); biomassYears<-topCatchYears

catchBarColor<-myGrey_trans

get_age_mat<-function(x){
  #x is species group code
  thisVar<-paste(x,"_age_mat", sep="")
  y<-grep(thisVar,biolLines); z<-get_first_number(biolLines[y])
  return(z)
}
## store SSB tracers for each group -only fill with age-structured
storeSSBfish<-array(NA, dim=c(nruns,ng,max_nts)); 
for(r in 1:nruns){
  for(g in 1:ng){
    thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
    if(thisNumCohorts>1){
      cat(as.character(thisCode))
      thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode],side="both")
      ##get age mature so know which cohorts are part of SSB
      this_age_mat<-get_age_mat(thisCode)
      #if only 2 cohorts, just the adults are SSB
      if(thisNumCohorts==2){
        thisVar<-paste(thisName,"2_Nums", sep=""); fishNums<-ncvar_get(nc_list[[r]],thisVar);
        thisVar<-paste(thisName,"2_ResN", sep=""); fishResN<-ncvar_get(nc_list[[r]],thisVar);
        thisVar<-paste(thisName,"2_StructN", sep=""); fishStructN<-ncvar_get(nc_list[[r]],thisVar);
        tempFishSSB<-apply(fishNums*(fishResN + fishStructN),3,sum)*mg_2_tonne*X_CN;
        storeSSBfish[r,g,burnin[r]:(nts_list[r]+burnin[r]-1)]<-tempFishSSB[burnin[r]:(nts_list[r]+burnin[r]-1)]
      } else{
        fishSSB<-rep(0,nts_list[r]); 
        for(c in (this_age_mat+1):thisNumCohorts){
          thisVar<-paste(thisName,c,"_Nums", sep="");  fishNums<-ncvar_get(nc_list[[r]],thisVar);
          thisVar<-paste(thisName,c,"_ResN", sep="");  fishResN<-ncvar_get(nc_list[[r]],thisVar);
          thisVar<-paste(thisName,c,"_StructN", sep="");  fishStructN<-ncvar_get(nc_list[[r]],thisVar);
          tempFishSSB<-apply(fishNums*(fishResN + fishStructN),3,sum)*mg_2_tonne*X_CN;
          fishSSB<-fishSSB + tempFishSSB[burnin[r]:(nts_list[r]+burnin[r]-1)]
        }
        storeSSBfish[r,g,burnin[r]:(nts_list[r]+burnin[r]-1)]<-fishSSB[burnin[r]:(nts_list[r]+burnin[r]-1)]
      }
    }
  }
}

#run legend
thisRunTemplate <-"AChaosFISH"
runText<-rep(NA, nruns)
for(r in 1:nruns){
  thisRunName<-this_out[r]; 
  xx<-gsub(thisRunTemplate,"", thisRunName); 
  temp <- gsub("Up|Down", "", xx); thisDir<-gsub(temp, "", xx)
  if(nchar(temp)==2){
    thisNum <- as.double(temp)*10
  }else{
    thisNum <- as.double(temp)
  }
  runText[r]<-paste(thisNum, "% ", thisDir, sep="")
  
}
runText[baseRunIndex]<-"Base"
thisPlotFile<-paste(plotPath,"EstSSBObsAndTotalCatchLEGEND",sep="")
# pdf(paste(thisPlotFile,".pdf",sep=""),height=3.5, width=5)
jpeg(paste(thisPlotFile, ".jpg", sep=""), quality=3000)
par(mar=c(0,0,0,0))
makeBlankPlot()
legend(legend=runText, col=runCols, lty=1, lwd=2, x="center", bty="n")
dev.off()

for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    cat(as.character(thisCode))
    if(thisCode %in% catchCodes){
      thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode],side="both")
      thisEst<-storeSSBfish[,g,];
      thisDescription<-gsub("_|_N", " ",groupsDFPaper$Name[groupsDF$Code==thisCode])
      thisDescription<-gsub("-","\n",thisDescription)
      thisDescription<-gsub("\\(","\n",thisDescription); thisDescription<-gsub(")", "", thisDescription)
      
      
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
      biomassYears<-xLabs
      biomassIndex<-biomassYears %in% topCatchYears
      biomassMax<-max(thisEst[,biomassIndex], na.rm=TRUE)*1.5
      biomassAxis<-pretty(seq(0,biomassMax,length.out=5))
      
      thisB0<-B0data$B0[B0data$Code==thisCode]
      
      # newBiomassAxis<-(biomassAxis/max(thisEst,na.rm=TRUE))*max(thisData/1000)
      
      thisYmax<-max(thisData/1000)*1.1
      if(file.exists(TSfile)){
        if(sum(thisObs,na.rm=TRUE)>0){
          # thisScale<-mean(thisEst[match(biomassYears, obsYears)], na.rm=TRUE)/mean(thisObs,na.rm=TRUE)
          newScaledReal<-(thisObs/mean(thisObs,na.rm=TRUE))*mean(thisEst[,biomassYears %in% obsYears], na.rm=TRUE)
          # newScaledReal<-(thisObs * thisScale)
          #get confidence intervals
          CIs<-getCIfromCVs(newScaledReal,thisCVs)
          # limitMax<-max(CIs$UpperCI, na.rm=TRUE)
          limitMax<-max(newScaledReal, na.rm=TRUE)
          if(limitMax>biomassMax){biomassMax<-limitMax}
        }
      }
      
      firstCatchYear<-min(topCatchYears[thisData>0]); lastYear<-max(topCatchYears)
      
      thisPlotFile<-paste(plotPath,"EstSSBObsAndTotalCatchbyYear",thisCode,sep="")
      # pdf(paste(thisPlotFile,".pdf",sep=""),height=3.5, width=5)
      jpeg(paste(thisPlotFile, ".jpg", sep=""), quality=3000)
      par(mar=c(3,4.2,3,4.2))
      par(las=0)
      # thisYmax<-90000
      plot(x=topCatchYears,y=thisData/1000,type="h",lwd=5,lend=1,xlab="",ylab="",col=catchBarColor,cex.axis=thisCex,ylim=c(0,thisYmax), xlim=c(firstCatchYear, lastYear))
      mtext(thisDescription,side=3,adj=0,line=0.2,cex=thisCex)
      mtext("Catch (tonnes)",side=2,cex=thisCex,line=3)
      plotBiomassYears<-topCatchYears[match(biomassYears,topCatchYears)]
      plotBiomassYears<-plotBiomassYears[!is.na(plotBiomassYears)]
      par(new=TRUE)
      plot(x=plotBiomassYears,y=(thisEst[baseRunIndex,biomassIndex]),type="l",col="black",cex=thisCex,ylim=c(0,biomassMax),
           xaxt="n",xlab="",lwd=1.5,yaxt="n",ylab="", xlim=c(firstCatchYear, lastYear), yaxt="n")
      for(r in 1:nruns){
        points(x=plotBiomassYears,y=(thisEst[r,biomassIndex]),type="l", col=runCols[r], lwd=1.5, lty=r)
      }
      # points(x=plotBiomassYears, y=thisEstBase[biomassIndex], type="l",lty=1,lwd=1.5)
      thisBiomassAxis<-pretty(seq(0,biomassMax), length.out = 5)
      axis(at=thisBiomassAxis,labels=thisBiomassAxis,side=4,cex.axis=thisCex)
      mtext("SSB (tonnes)",side=4,adj=0.5,line=2.5,cex=thisCex)
      if(file.exists(TSfile) & thisCode %in% surveyOK){
        if(sum(thisObs,na.rm=TRUE)>0){
          points(x=obsYears,y=newScaledReal,pch=20,col="red",cex=1.2)
          segments(y0=CIs$LowerCI, x0=obsYears, y1=CIs$UpperCI, x1=obsYears, col="red")
        }
      }
      dev.off()
    }
    # } else {
    #   ## just plot biomass, and survey if it is there
    #   thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode],side="both")
    #   thisEst<-storeSSBfish[,g,]; 
    #   #trawl survey
    #   thisObs<-rep(NA,length(biomassYears)); obsYears<-seq(1,length(biomassYears))
    #   TSfile<-paste(DIR$'Tables',"TrawlSurvey\\",thisCode,".csv",sep="")
    #   if(file.exists(TSfile)  & thisCode %in% surveyOK){
    #     thisTS<-read.csv(TSfile)
    #     thisObs<-thisTS$Biomass[thisTS$Year %in% biomassYears]; obsYears<-thisTS$Year[thisTS$Year %in% biomassYears]
    #     thisCVs<-thisTS$cv[thisTS$Year %in% biomassYears]/100; 
    #   }
    #   biomassYears<-xLabs
    #   biomassIndex<-biomassYears %in% biomassYears
    #   biomassMax<-max(thisEst[,biomassIndex], na.rm=TRUE)*1.5
    #   biomassAxis<-pretty(seq(0,biomassMax,length.out=5))
    #   
    #   thisB0<-B0data$B0[B0data$Code==thisCode]
    #   
    #   if(file.exists(TSfile)){
    #     if(sum(thisObs,na.rm=TRUE)>0){
    #       # thisScale<-mean(thisEst[match(biomassYears, obsYears)], na.rm=TRUE)/mean(thisObs,na.rm=TRUE)
    #       newScaledReal<-(thisObs/mean(thisObs,na.rm=TRUE))*mean(thisEst[,biomassYears %in% obsYears], na.rm=TRUE)
    #       # newScaledReal<-(thisObs * thisScale)
    #       #get confidence intervals
    #       CIs<-getCIfromCVs(newScaledReal,thisCVs)
    #       # limitMax<-max(CIs$UpperCI, na.rm=TRUE)
    #       limitMax<-max(newScaledReal, na.rm=TRUE)
    #       if(limitMax>biomassMax){biomassMax<-limitMax}
    #     }
    #   }
    #   
    #   thisPlotFile<-paste(plotPath,"EstSSBObsAndTotalCatchbyYear",thisCode,sep="")
    #   # pdf(paste(thisPlotFile,".pdf",sep=""),height=3.5, width=5)
    #   jpeg(paste(thisPlotFile, ".jpg", sep=""), quality=3000)
    #   par(mar=c(3,4,3,4.5))
    #   par(las=0)
    #   # thisYmax<-90000
    #   plot(x=topCatchYears,y=thisData/1000,type="h",lwd=5,lend=1,xlab="",ylab="",col=catchBarColor,cex.axis=thisCex,ylim=c(0,thisYmax), xlim=c(firstCatchYear, lastYear))
    #   mtext(thisDescription,side=3,adj=0,line=0.2,cex=thisCex)
    #   plotBiomassYears<-biomassYears
    #   par(new=TRUE)
    #   plot(x=plotBiomassYears,y=(thisEst[baseRunIndex,biomassIndex]),type="l",col="black",cex=thisCex,ylim=c(0,biomassMax),
    #        xaxt="n",xlab="",lwd=3,yaxt="n",ylab="", xlim=c(firstCatchYear, lastYear), yaxt="n")
    #   for(r in 1:nruns){
    #     points(x=plotBiomassYears,y=(thisEst[r,biomassIndex]),type="l", col=runCols[r], lwd=1.5)
    #   }      
    #   thisBiomassAxis<-pretty(seq(0,biomassMax), length.out = 5)
    #   axis(at=thisBiomassAxis,labels=thisBiomassAxis,side=2,cex.axis=thisCex)
    #   mtext("SSB (tonnes)",side=2,adj=0.5,line=2.5,cex=thisCex)
    #   if(file.exists(TSfile)){
    #     if(sum(thisObs,na.rm=TRUE)>0){
    #       points(x=obsYears,y=newScaledReal,pch=20,col="red",cex=1.2)
    #       segments(y0=CIs$LowerCI, x0=obsYears, y1=CIs$UpperCI, x1=obsYears, col="red")
    #     }
    #   }
    #   dev.off()
    # }
  }
}
