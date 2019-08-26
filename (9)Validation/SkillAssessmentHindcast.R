#bring in the run, the observed data, the catch history and plot together
source(paste(DIR$'General functions',"getCIfromCV.R", sep=""))
source(paste(DIR$'General functions',"getSNRN.R", sep=""))

nboxes<-30

mg_2_tonne<-2e-8; X_CN<-5.7

#################################### TRACERS
#plot all tracers for a given box and layer
this_run<-"base"; base_out<-this_run

this_out<-"FISH"

# define groups for which survey is appropriate for estimating abundance
surveyOK<-c("BIS", "CBO",  "ELI", "HAK", "HOK", "JAV", "LDO", "LIN",  "SPE")

burnin<-35 #number of years to skip in plot

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")
baseOutPath<-paste(this_path,"output",base_out,"\\", sep="")
biolLines<-readLines(paste(this_path,"..\\CRAM_BH_hybrid_biol.prm", sep=""))

#read in groups csv file to get codes for fished groups
groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

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

modelYears<-seq(year0,(year0+nts-1))

thisCex<-1.1

#get all tracer names
allTracers<-sort(names(ThisNC.nc$var))

##store biomass tracers
storeSSB<-array(NA,dim=c(ng,nts))
storeObs<-array(NA,dim=c(ng, nts))
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    cat(as.character(thisCode))
    thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode],side="both")
    if(thisNumCohorts>2){
      #cohort mature
      thisVar<-paste(thisCode,'_age_mat', sep=""); x<-grep(thisVar, biolLines)
      ageMat<-get_first_number(biolLines[x])
      SSB<-rep(0,nts)
      #get numbers, resn, structn, then multiply for biomass
      for(c in (ageMat+1):thisNumCohorts){
        tempvar<-paste(thisName,c,"_Nums",sep="")
        temp<-ncvar_get(ThisNC.nc,tempvar); tempNums<-apply(temp,3,sum)
        tempvar<-paste(thisName,c,"_ResN",sep="")
        temp<-ncvar_get(ThisNC.nc,tempvar); tempResN<-apply(temp,3,nonZeroMean)
        tempvar<-paste(thisName,c,"_StructN",sep="")
        temp<-ncvar_get(ThisNC.nc,tempvar); tempStructN<-apply(temp,3,nonZeroMean)
        
        tempBiomass<-tempNums*(tempResN+tempStructN)
        SSB=SSB+tempBiomass[(burnin+1):(burnin+nts)]
      }
    } else{
        c=2
        tempvar<-paste(thisName,c,"_Nums",sep="")
        temp<-ncvar_get(ThisNC.nc,tempvar); tempNums<-apply(temp,3,sum)
        tempvar<-paste(thisName,c,"_ResN",sep="")
        temp<-ncvar_get(ThisNC.nc,tempvar); tempResN<-apply(temp,3,nonZeroMean)
        tempvar<-paste(thisName,c,"_StructN",sep="")
        temp<-ncvar_get(ThisNC.nc,tempvar); tempStructN<-apply(temp,3,nonZeroMean)
        
        tempBiomass<-tempNums*(tempResN+tempStructN) *mg_2_tonne * X_CN
        SSB=tempBiomass[(burnin+1):(burnin+nts)]
        
    }
    storeSSB[g,]<-SSB

    #trawl survey
    TSfile<-paste(DIR$'Tables',"TrawlSurvey\\",thisCode,".csv",sep="")
    if(file.exists(TSfile)){
      thisTS<-read.csv(TSfile)
      obsYears<-thisTS$Year; thisObs<-thisTS$Biomass
      yearIndex<-modelYears %in% obsYears
      newScaledReal<-(thisObs/mean(thisObs,na.rm=TRUE))*mean(SSB[yearIndex])
      storeObs[g,yearIndex]<-newScaledReal
    }
  }
}

# calc model efficiency
calcModEff<-function(O,P){
  A<-0; B<-0; n<-length(O)
  mO<-mean(O, na.rm=TRUE);
  for(i in 1:n){
    A=A + (O[i] - mO)^2; B=B + (P[i] - O[i])^2
  }
  mE<-(A-B)/A
  return(mE)
}
calcRMSE<-function(O,P){
  A<-0; n<-length(O)
  for(i in 1:n){
    A = A + (P[i] - O[i])^2
  }
  rmse<-sqrt(A/n)
  return(rmse)
}

calcRI<-function(O, P){
  A=0; n<-length(O)
  for(i in 1:n){
    A = A + (log(O[i]/P[i]))^2
  }
  RI<-exp(sqrt(A/n))
  return(RI)
}
calcAAE<-function(O, P){
  A=0; n=length(O)
  for(i in 1:n){
    A <- A + abs(O[i] - P[i])
  }
  AAE<-A/n
  return(AAE)
}

storeModEffs<-rep(NA,ng); storePcorr<-rep(NA,ng); storeRMSE<-rep(NA,ng); storeRI<-rep(NA,ng); storeAAE<-rep(NA,ng)
for(g in 1:ng){
  thisObs<-storeObs[g,]; thisEst<-storeSSB[g,]
  thisCode<-groupsDF$Code[g]
  if(thisCode %in% surveyOK){
    if(sum(thisObs, na.rm=TRUE)>0 & sum(thisEst, na.rm=TRUE)>0){
      index<-!is.na(thisObs)
      thisME<-calcModEff(O=thisObs[index], P=thisEst[index])
      storeModEffs[g]<-thisME
      thisCor<-rcorr(x=thisObs[index], y=thisEst[index])
      storePcorr[g]<-thisCor[[1]][2,1]
      thisRMSE<-calcRMSE(O=thisObs[index], P=thisEst[index])
      storeRMSE[g]<-thisRMSE
      storeRI[g]<-calcRI(O=thisObs[index], P=thisEst[index])
      storeAAE[g]<-calcAAE(O=thisEst[index],P= thisObs[index])
    }
  }
}

# testing RI
testO<-rep(10,10); testP<-c(rep(20,5), rep(5,5));
testLog<-log(testO/testP)

testP<-c(rep(10*1.5,5), rep((10*(2/3)),5));
testLog<-log(testO/testP)
calcRI(O=testO, P=testP)

testP<-c(rep(10*1.5,8), rep((10*(2/3)),3));
calcRI(O=testO, P=testP)

testP<-c(rep(10*1.5,6), rep((10*(2/3)),3),20);
calcRI(O=testO, P=testP)

# pdf(paste(plotPath,"SkillAssessment.pdf", sep=""),height=6, width=9)
# par(mfrow=c(2,2), mar=c(4,4.5,1,0.5))
# index<-!is.na(storeAAE)
# par(lend=1)
# plot(storeAAE[index], type="h", lwd=5, xaxt="n",ylab="AAE", xlab="")
# par(las=2)
# axis(at=seq(1,length(storeAAE[index])), labels=as.character(groupsDF$Code[index]), side=1)
index<-!is.na(storeRI); 
xx<-as.character(groupsDF$LongName[index])
thisLabels<-gsub("_|_N", " ", xx)
# edit ELI manually
x <- grep("Elasmobranchs - invertivores", thisLabels); thisLabels[x]<-"Elasmobranchs (invertivores)"

pdf(paste(plotPath,"SkillAssessmentRI.pdf", sep=""),height=5, width=5)
par(mar=c(14,4,1,1))

par(lend=1)
plot(storeRI[index], type="h", lwd=5, xaxt="n",ylab="RI", xlab="", ylim=c(0, max(storeRI, na.rm=TRUE)))
par(las=2)
axis(at=seq(1,length(storeRI[index])), labels=thisLabels, side=1)
abline(h=1,col="red",lty=1)
dev.off()


pdf(paste(plotPath,"SkillAssessmentMFE.pdf", sep=""),height=5, width=5)
par(mar=c(14,4,1,1))
index<-!is.na(storeModEffs); 
par(lend=1)
plot(storeModEffs[index], type="h", lwd=5, xaxt="n",ylab="MEF", xlab="")
par(las=2)
axis(at=seq(1,length(storeModEffs[index])), labels=thisLabels, side=1)
dev.off()

pdf(paste(plotPath,"SkillAssessmentPearsonsCorr.pdf", sep=""),height=5, width=4)
par(mar=c(14,4,1,1))
index<-!is.na(storePcorr);
par(lend=1)
plot(storePcorr[index], type="h", lwd=5, xaxt="n",ylab="Pearson's correlation", xlab="", cex.axis=thisCex, cex.lab=thisCex)
abline(h=0); par(las=2)
axis(at=seq(1,length(storePcorr[index])), labels=thisLabels, side=1, cex.axis=thisCex)
dev.off()

################### consider error bounds on observations
calcRI_CI<-function(O, P, CVs){
  A=0; n<-length(O)
  meanObs<-mean(O, na.rm=TRUE); CIlower<-O-2*meanObs*(CVs/100); CIupper<-O+2*meanObs * (CVs/100)
  for(i in 1:n){
    if(P[i]<=CIupper[i] & P[i]>=CIlower[[i]]){P[i]<-O[i]}
    C<-(log(O[i]/P[i]))^2
    A = A + C
  }
  RI<-exp(sqrt(A/n))
  return(RI)
}
calcModEff_CI<-function(O,P, CVs){
  meanObs<-mean(O, na.rm=TRUE); SEs<-meanObs*(CVs/100); 
  A<-0; B<-0; n<-length(O)
  mO<-mean(O, na.rm=TRUE);
  for(i in 1:n){
    A=A + (O[i] - mO)^2/SEs[i]; B=B + (P[i] - O[i])^2/SEs[i]
  }
  mE<-(A-B)/A
  return(mE)
}
## treats predictions within bounds as equal to obs
calcModEff_CI<-function(O,P, CVs){
  meanObs<-mean(O, na.rm=TRUE); CIlower<-O-2*meanObs*(CVs/100); CIupper<-O+2*meanObs * (CVs/100)
  A<-0; B<-0; n<-length(O)
  mO<-mean(O, na.rm=TRUE);
  for(i in 1:n){
    thisP<-P[i]; if(thisP>=CIlower[i] & thisP<=CIupper[i]){thisP<-O[i]}
    A=A + (O[i] - mO)^2; B=B + (thisP - O[i])^2
  }
  mE<-(A-B)/A
  return(mE)
}

storeModEffs_CI<-rep(NA,ng); storeRI_CI<-rep(NA,ng); 
for(g in 1:ng){
  thisObs<-storeObs[g,]; thisEst<-storeSSB[g,]
  thisCode<-groupsDF$Code[g]
  if(thisCode %in% surveyOK){
    if(sum(thisObs, na.rm=TRUE)>0 & sum(thisEst, na.rm=TRUE)>0){
      ## read in CIs for obs
      TSfile<-paste(DIR$'Tables',"TrawlSurvey\\",thisCode,".csv",sep="")
        thisTS<-read.csv(TSfile)
        thisCVs<-thisTS$cv
   
      index<-!is.na(thisObs)
      
      thisME<-calcModEff_CI(O=thisObs[index], P=thisEst[index], CV=thisCVs)
      storeModEffs_CI[g]<-thisME
      
      thisRI<-calcRI_CI(O=thisObs[index], P=thisEst[index], CVs=thisCVs)
      storeRI_CI[g]<-thisRI
    }
  }
}

pdf(paste(plotPath,"SkillAssessmentRI.pdf", sep=""),height=5, width=4)
par(mar=c(14,4,1,1))
index<-!is.na(storeRI); 
par(lend=1)
plot(x=seq(1,length(storeRI_CI[index]))-0.1,y=storeRI[index], type="h", lwd=5, xaxt="n",ylab="RI", xlab="", ylim=c(0, max(storeRI, na.rm=TRUE)), cex.axis=thisCex, cex.lab=thisCex)
par(las=2)
axis(at=seq(1,length(storeRI[index])), labels=thisLabels, side=1, cex.axis=thisCex)
abline(h=1,col=myGrey_trans, lwd=3,lty=1)
points(x=seq(1,length(storeRI_CI[index]))+0.1,y=storeRI_CI[index], type="h", lwd=5, col=myOrange)
dev.off()

thisYmin<-min(c(storeModEffs[index], storeModEffs_CI[index])); thisYmax<-max(c(storeModEffs[index], storeModEffs_CI[index]))
pdf(paste(plotPath,"SkillAssessmentMFE.pdf", sep=""),height=5, width=4)
par(mar=c(14,4,1,1))
index<-!is.na(storeModEffs); 
par(lend=1)
plot(x=seq(1,length(storeModEffs_CI[index]))-0.1,y=storeModEffs[index], type="h", lwd=5, xaxt="n",ylab="MEF", xlab="", ylim=c(thisYmin, thisYmax), cex.axis=thisCex, cex.lab=thisCex)
abline(h=0); abline(h=1,col=myGrey_trans, lwd=3)
par(las=2)
axis(at=seq(1,length(storeModEffs[index])), labels=thisLabels, side=1, cex.axis=thisCex)
points(x=seq(1,length(storeModEffs_CI[index]))+0.1,y=storeModEffs_CI[index], type="h", lwd=5, col=myOrange)
dev.off()
