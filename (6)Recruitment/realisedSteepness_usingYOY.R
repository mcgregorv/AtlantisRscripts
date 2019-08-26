#read in the tracers (numbers), calculate spawning stock abundance, calculate rectruited abundance, then calculate steepness
this_run<-"base"
# this_out<-paste("Base2",sep="")
this_out<-paste("TEST",sep="")
this_out<-"SENSselBASE"
# 
htext<-"_h"; thisSteepness<-0.75

htext<-"";  thisSteepness<-"";

basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
# outPath<-paste(basePath,this_run,"\\","output",this_out,"_h",thisSteepness,"\\",sep="")
outPath<-paste(basePath,this_run,"\\","output",this_out,htext,thisSteepness,"\\",sep="")


plotPath<-paste(basePath,"\\Figures\\Recruitment\\",this_out,sep="")

baseBiolFile<-paste(basePath,"CRAM_BH_hybrid_biol.prm",sep="")
# thisBiolFile<-paste(basePath,"CRAM_base_biol_h",thisSteepness,".prm",sep="")
# if(!file.exists(thisBiolFile)){
#   biolLines<-readLines(baseBiolFile)
# }else{
  biolLines<-readLines(baseBiolFile)
# }

groupsDF<-read.csv(paste(basePath,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1900
fishingStartYear<-1900
modelStartYear<-1900

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
nts<-dim(thisVol)[3]; nboxes<-dim(thisVol)[2]; nlayers<-dim(thisVol)[1]

allTracers<-sort(unique(names(ThisNC.nc$var)))
numTracers<-allTracers[grep("_Nums",allTracers)]

YOY_df<-read.csv(paste(outPath,"outputYOY.txt",sep=""),sep=" ")

getBHab<-function(Code){
  out<-list()
  thisPar<-paste("BHalpha_",Code,sep="")
  x<-biolLines[grep(thisPar,biolLines)]; this_a<-get_first_number(x)
  out$'a'<-this_a
  thisPar<-paste("BHbeta_",Code,sep="")
  x<-biolLines[grep(thisPar,biolLines)]; this_b<-get_first_number(x)
  out$'b'<-this_b
  return(out)
}
get_recruitmentModel<-function(Code){
  thisPar<-paste("flagrecruit",Code,sep="")
  thisOut<-get_first_number(biolLines[grep(thisPar,biolLines)])
  return(thisOut)
}
get_KWRR<-function(Code){
  thisPar<-paste("KWRR_",Code,sep="")
  thisOut<-get_first_number(biolLines[grep(thisPar,biolLines)])
  return(thisOut)
}
get_KWSR<-function(Code){
  thisPar<-paste("KWSR_",Code,sep="")
  thisOut<-get_first_number(biolLines[grep(thisPar,biolLines)])
  return(thisOut)
}
calcR<-function(SSB,a,b){
  R<-(a*SSB)/(b+SSB)
  return(R)
}

#set up array to store SS and R for each group for each ts
storeSS<-array(0,dim=c(ng,nts)); 

for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[g]); thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    blPar<-paste("FSPB_",thisCode,sep="")
    x<-biolLines[grep(blPar,biolLines)+1]; thisPropSpawning<-get_first_number(x,n="all")
    blPar<-paste(thisCode,"_AgeClassSize",sep="")
    x<-biolLines[grep(blPar,biolLines)]; thisAgeClassSize<-get_first_number(x)
    weight_array<-array(NA,dim=c(thisNumCohorts,nts)); numbers_array<-weight_array
    for(c in 1:thisNumCohorts){
      thisTracer<-paste(thisName,c,"_ResN",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer);   thisRN<-apply(thisTemp,3,nonZeroMean)
      thisTracer<-paste(thisName,c,"_StructN",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer); thisSN<-apply(thisTemp,3,nonZeroMean)
      weight_array[c,]<-(thisRN+thisSN)*mg_2_tonne*X_CN
      #numbers
      thisTracer<-paste(thisName,c,"_Nums",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer);
      numbers_array[c,]<-apply(thisTemp,3,sum,na.rm=TRUE)
    }
    thisSSB<-apply(numbers_array*weight_array,2,sum)
    storeSS[g,]<-thisSSB
  }
}

storeSSsum<-storeSS

mg_2_tonne<- 0.00000002; X_CN<-5.7

thisCex<-1.5
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[g]); thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    thisRM<-get_recruitmentModel(thisCode)
    if(thisRM %in% c(3,10,19)){
      jpeg(paste(plotPath,"BHCompareRealised_h",thisSteepness,"_",thisCode,".jpg",sep=""),quality=300,height=400,width=500)
      # par(mfrow=c(2,1))
      #take out the first one because it is calculated slightly differently
      testHokSS<-storeSSsum[g,]; thisYOY<-YOY_df[,grep(thisCode,colnames(YOY_df))]; testB0<-testHokSS[1]/(mg_2_tonne * X_CN)
      thisKWRR<-get_KWRR(thisCode); thisKWSR<-get_KWSR(thisCode); thisRecruitWeight<-thisKWRR+thisKWSR
      testHokR<-thisYOY/(thisRecruitWeight*X_CN * mg_2_tonne)
      this_a<-getBHab(thisCode)$'a'; this_b<-getBHab(thisCode)$'b'; this_b_tonnes<-this_b*mg_2_tonne * X_CN
      
      compareR<-unlist(lapply(testHokSS/(mg_2_tonne*X_CN),calcR,a=this_a,b=this_b))
      thisR0<-calcR(testB0, a=this_a, b=this_b)
      thisR20<-calcR(testB0*0.2, a=this_a, b=this_b); this_h<-thisR20/thisR0
      
      fitBHcurveSS<-seq(0,max(testHokSS/(mg_2_tonne*X_CN))*1.2,length.out=1e+4); BHcurveR<-unlist(lapply(fitBHcurveSS,calcR,a=this_a,b=this_b))
      BHcurveSS<-fitBHcurveSS*(mg_2_tonne*X_CN)
      thisYMax<-max(c(this_a,max(testHokR))); thisXMax<-max(c(max(BHcurveSS),max(testHokSS)))
      if(max(testHokR)>10*this_a){
        thisYMax<-max(c(this_a,max(testHokR[testHokR<max(testHokR)])))*1.2;
      }
      
         ## testing; read in localrecruit and add to plot
      thisVar<-paste("",thisCode, "_Time_Spawn", sep=""); temp<-biolLines[grep(thisVar, biolLines)]; thisfl<-get_first_number(temp)
      
      plot(BHcurveSS,BHcurveR,type="l",lwd=2,col=myBlue,ylim=c(0,thisYMax),xlim=c(0,thisXMax),ylab="Recruitment (numbers)",xlab="SSB (tonnes)",cex.axis=thisCex,cex.lab=thisCex,cex=thisCex)
      points(testHokSS,testHokR,pch=8,col=myRed,cex=thisCex)
      mtext(thisCode,side=3,cex=thisCex,adj=0); mtext(thisfl,side=3,cex=thisCex,adj=1); 
      abline(h=this_a,col=myGrey,lwd=2,lty=2)
      # points(testHokSS,this_calcR,col="red",pch=20,cex=thisCex)
      abline(v=this_b_tonnes,col=myGrey,lwd=2,lty=2)
      abline(h=thisR0, col=myGrey_trans, lwd=3); abline(v=testB0*mg_2_tonne*X_CN, col=myGrey_trans, lwd=3)
      # 
      # plot(testHokR/compareR, type="l")
      # mtext(signif((testHokR/compareR)[length(testHokR/compareR)],2), side=3, adj=1)
      dev.off()
    }
  }
}


##read in Time_Spawn and Spawn_Period, add them together and print out those that stradle 365
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    thisVar<-paste(thisCode, "_Time_Spawn", sep=""); temp<-biolLines[grep(thisVar, biolLines)]; thisTimeSpawn<-get_first_number(temp)
    thisVar<-paste(thisCode, "_spawn_period", sep=""); temp<-biolLines[grep(thisVar, biolLines)]; thisSpawnPeriod<-get_first_number(temp)
    endSpawn<-thisTimeSpawn + thisSpawnPeriod
    if(thisTimeSpawn<365 & endSpawn>365){cat(as.character(thisCode), "--", thisSpawnPeriod,", ", endSpawn,"\n")}
    
  }
}

##linear recruitment
##Linear
get_KDENR<-function(Code){
  #this is the linear recruitment coefficient
  thisPar<-paste("KDENR_",Code,sep="")
  thisOut<-get_first_number(biolLines[grep(thisPar,biolLines)+1])
  return(thisOut)
}

#what if divide by recruit period?
get_RP<-function(Code){
  #this is the linear recruitment coefficient
  thisPar<-paste("Recruit_Period_",Code,sep="")
  thisOut<-get_first_number(biolLines[grep(thisPar,biolLines)+1])
  return(thisOut)
}

for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[g]); thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    thisRM<-get_recruitmentModel(thisCode)
    if(thisRM==12){
      this_KDENR<-get_KDENR(thisCode)
      jpeg(paste(plotPath,"LinearCompareRealised_h",thisSteepness,"_",thisCode,".jpg",sep=""),quality=300,height=300)
      #take out the first one because it is calculated slightly differently
      testHokSS<-storeSSsum[g,-1]; thisYOY<-YOY_df[-1,grep(thisCode,colnames(YOY_df))]
      thisRP<-get_RP(thisCode)
      thisKWRR<-get_KWRR(thisCode); thisKWSR<-get_KWSR(thisCode); thisRecruitWeight<-thisKWRR+thisKWSR
      testHokR<-(thisYOY/(thisRecruitWeight*X_CN * mg_2_tonne))
      linearSS<-seq(min(testHokSS),max(testHokSS),length.out=1000); linearR<-this_KDENR*linearSS
      thisYMax<-max(c(max(linearR),max(testHokR))); thisXMax<-max(c(max(BHcurveSS),max(testHokSS)))

      plot(linearSS,linearR,type="l",lwd=2,col=myBlue,ylim=c(0,thisYMax),ylab="Recruitment (numbers)",xlab="Spawning Stock Abundance (numbers)",cex.axis=thisCex,cex.lab=thisCex,cex=thisCex)
      points(testHokSS,testHokR,pch=20,col="red",cex=thisCex)
      mtext(thisCode,side=3,cex=thisCex,adj=0)
      
      # par(new=TRUE)
      # ##plot instead the difference between intended linear recruitment and actual recruitment
      # expectedRec<-this_KDENR*testHokSS
      # diffRec<-testHokR-expectedRec
      # plot(x=testHokSS,y=diffRec,pch=20,col=myGreen,yaxt="n",xaxt="n")
      # thisAxis<-signif(seq(1,max(diffRec),length.out = 5),2)
      # axis(at=thisAxis,labels = thisAxis,side=4,col=myGreen)
      
      dev.off()
    }
  }
}
##########

##########
#estimate scaling for KDENR to get the resulting recruitment what it should be
newKDENRvalues<-data.frame(matrix(NA,ncol=3,nrow=0)); colnames(newKDENRvalues)<-c("Code","KDENR","newKDENR")
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[g]); thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    thisRM<-get_recruitmentModel(thisCode)
    if(thisRM==12){
      this_KDENR<-get_KDENR(thisCode)
      #take out the first one because it is calculated slightly differently
      testHokSS<-storeSSsum[g,-1]; thisYOY<-YOY_df[-1,grep(thisCode,colnames(YOY_df))]
      thisRP<-get_RP(thisCode)
      thisKWRR<-get_KWRR(thisCode); thisKWSR<-get_KWSR(thisCode); thisRecruitWeight<-thisKWRR+thisKWSR
      testHokR<-(thisYOY/(thisRecruitWeight*X_CN * mg_2_tonne))[1:length(testHokSS)]
      linearSS<-seq(min(testHokSS),max(testHokSS),length.out=1000); linearR<-this_KDENR*linearSS
      
      thisScale<- max(testHokR/linearR)
      newKDENR<-signif(this_KDENR/thisScale,3)
      
      newKDENRvalues[g,]<-c(as.character(thisCode),this_KDENR,newKDENR)
    }
  }
}

newKDENRvalues<-newKDENRvalues[!is.na(newKDENRvalues$Code),]

write.csv(newKDENRvalues,file=paste(basePath,"CRAM_newKDENRpars.csv",sep=""))


##test steepness
thisCode<-"HOK"; this_a<-getBHab(thisCode)$'a'; this_b<-getBHab(thisCode)$'b'
hokR0<-745264763; hokSSB0<-734584744
hokSSB0.2<-0.2*hokSSB0
hokR0.2<-calcR(hokSSB0.2,this_a,this_b)
hokTesth<-hokR0.2/hokR0

######### http://www2.kenyon.edu/Depts/Math/hartlaub/Math224%20Fall2008/Markov-Sample2.pdf
T<-matrix(c(.7,.1,.3,.2,.8,.3,.1,.1,.4), 3, 3, byrow = TRUE)
ev<-eigen(T)
P<-c(0.48,0.51,0.01)
T %*% P

ev1<-ev$vectors