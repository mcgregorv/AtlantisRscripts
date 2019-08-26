#plot all tracers for a given box and layer
this_run<-"chaos"

this_out<-"Short2"
this_out<-"TestFish"
# this_out<-"ZeroSomePred"
# this_out<-"ShortLonger"
this_out<-"FishNUMS4"
# this_out<-"Chla"
# this_out<-"NUMS4"
this_out<-"DChaosBASE"

burnin<-0 #number of years to skip in plot

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")

plotPath<-paste(this_path,"..\\Figures\\",this_run,"\\",this_out,"",sep="")

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1865
fishingStartYear<-1865
modelStartYear<-1865

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
numTracers<-allTracers[grep("_Nums",allTracers)]; nnt<-length(numTracers)

cohortCols<-colorRampPalette(colors=c(myYellow,myGreen,myAqua,myBlue,myPurple,myRed,"red"))(10)

groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep=""))
ageGroupsDF<-groupsDF[groupsDF$NumCohorts>1,]; nag<-dim(ageGroupsDF)[1]

thisBiolFile<-paste(this_path,"..\\CRAM_base_biol.prm",sep="")
biolLines<-readLines(thisBiolFile)


thisNumsDF<-array(NA,dim=c(nag,nts,10))
thisMortDF<-array(NA,dim=c(nag,nts,10))


for(g in 1:nag){
  thisCode<-ageGroupsDF$Code[g]
  thisName<-str_trim(ageGroupsDF$Name[g],side="both")
  thisNumCohorts<-ageGroupsDF$NumCohorts[g]
  yy<-paste(thisCode,"_AgeClassSize",sep=""); xxx<-biolLines[grep(yy,biolLines)]
  thisAgeClassSize<-get_first_number(xxx)
  #just doing for cohorts 2-numcohorts for now - need recruits to do c 1
  for(c in 1:thisNumCohorts){
    thisTracer<-paste(thisName,c,"_Nums",sep="")
    temp<-ncvar_get(ThisNC.nc,thisTracer)
    thisNumsDF[g,,c]<-apply(temp,3,sum)
    if(c>1){
      for(t in 2:nts){
        xx<-thisNumsDF[g,t-1,c]*((thisAgeClassSize-1)/thisAgeClassSize)+thisNumsDF[g,t-1,c-1]/thisAgeClassSize
        xxx<-thisNumsDF[g,t,c]/xx
        thisM<-(-1)*log(xxx,base=exp(1))
        thisMortDF[g,t,c]<-thisM
      }
    }
  }
}

## do alt mortality - simpler, I think
altMortDF<-array(NA,dim=c(nag,nts,10))

for(g in 1:nag){
  thisCode<-ageGroupsDF$Code[g]
  thisName<-str_trim(ageGroupsDF$Name[g],side="both")
  thisNumCohorts<-ageGroupsDF$NumCohorts[g]
  yy<-paste(thisCode,"_AgeClassSize",sep=""); xxx<-biolLines[grep(yy,biolLines)]
  thisAgeClassSize<-get_first_number(xxx)
  #just doing for cohorts 2-numcohorts for now - need recruits to do c 1
  for(c in 1:thisNumCohorts){
    thisTracer<-paste(thisName,c,"_Nums",sep="")
    temp<-ncvar_get(ThisNC.nc,thisTracer)
    thisNumsDF[g,,c]<-apply(temp,3,sum)
    if(c>1){
      for(t in (thisAgeClassSize+1):nts){
        xx<-thisNumsDF[g,t-thisAgeClassSize,c-1]
        xxx<-thisNumsDF[g,t,c]/xx
        thisM<-(-1/thisAgeClassSize)*log(xxx,base=exp(1))
        altMortDF[g,t,c]<-thisM
      }
    }
  }
}

##read in B0 and R0 and calc base M
B0df<-read.csv(paste(this_path,"..\\CRAM_B0.csv",sep=""))

calcM<-function(m_daily_prop){
  thisM<-0
  if(!is.na(m_daily_prop)){
    if(m_daily_prop>0){
      thisM<-(-1)*log(1-m_daily_prop,base=exp(1))*365
    }
  }
  return(thisM)
}

##get mortality from mL and mQ for each group and cohort
MfromML<-array(NA,dim=c(nag,nts,10)); MfromMQ<-MfromML
for(g in 1:nag){
  thisCode<-ageGroupsDF$Code[g]
  thisNumCohorts<-ageGroupsDF$NumCohorts[g]
  xx<-paste(thisCode,"_mL",sep=""); xxx<-biolLines[grep(xx,biolLines)+1]
  tempMLbyCohort<-get_first_number(xxx,n="all")
  xx<-paste(thisCode,"_mQ",sep=""); xxx<-biolLines[grep(xx,biolLines)+1]
  tempMQbyCohort<-get_first_number(xxx,n="all")
  #get age mature
  thisMat<-get_first_number(biolLines[grep(paste(thisCode,"_age_mat",sep=""),biolLines)]) 
  mLbyCohort<-c(rep(tempMLbyCohort[1],thisMat),rep(tempMLbyCohort[2],(thisNumCohorts-thisMat)))
  mQbyCohort<-c(rep(tempMQbyCohort[1],thisMat),rep(tempMQbyCohort[2],(thisNumCohorts-thisMat)))
  thisMfromML<-unlist(lapply(mLbyCohort,calcM))
  #MQ is more complicated as depends on abundance
  deathsFromMQ<-(mQbyCohort*thisNumsDF[g,,])^2
  propMortFromQ<-deathsFromMQ/thisNumsDF[g,,]
  thisMfromMQ<-apply(propMortFromQ,c(1,2),calcM)
  thisBackgroundML<-0*thisMfromMQ; thisBackgroundMQ<-thisBackgroundML
  for(c in 1:thisNumCohorts){
    thisBackgroundML[,c]<-0*thisBackgroundML[,c]+thisMfromML[c]
    thisBackgroundMQ[,c]<-0*thisBackgroundMQ[,c]+thisMfromMQ[,c]
  }
  MfromML[g,,]<-thisBackgroundML
  MfromMQ[g,,]<-thisBackgroundMQ
}
  

pdf(paste(plotPath,"RealisedMbyGroup.pdf",sep=""))
for(g in 1:nag){
  thisCode<-ageGroupsDF$Code[g]
  thisNumCohorts<-ageGroupsDF$NumCohorts[g]
  thisMort<-thisMortDF[g,,]
  thisBaseM<-B0df$M[B0df$Code==thisCode]
  thisMax<-max(max(thisMort,na.rm=TRUE),thisBaseM,na.rm=TRUE)
  thisMin<-min(thisMort,na.rm=TRUE)
  if(!is.na(thisMax)){
    #get recruitment method
    yy<-paste("flagrecruit",thisCode,sep=""); yyy<-biolLines[grep(yy,biolLines)]
    thisFlagRecruit<-get_first_number(yyy)
    thisFlagDescr<-ifelse(thisFlagRecruit==12,"Constant","BH")
    plot(x=seq(1,nts),y=rep(thisMax,nts),type="n",xlab="",ylab="M",xaxt="n",ylim=c(min(0,thisMin),thisMax))
    for(c in 1:thisNumCohorts){
      points(thisMort[,c],type="l",lwd=2,col=cohortCols[c])
      points(MfromML[g,,c],type="l",lwd=2,col=cohortCols[c],lty=2)
      points(MfromMQ[g,,c],type="l",lwd=2,col=cohortCols[c],lty=3)
    }
    axis(at=xLabsAt,labels=xLabs,side=1)
    mtext(thisCode,side=3,adj=0)
    mtext(thisFlagDescr,side=3,adj=1)
    abline(h=thisBaseM,col=myGrey,lwd=2,lty=2)
  }
}
dev.off()

##plot alt mort

pdf(paste(plotPath,"RealisedMbyGroupAlt.pdf",sep=""))
for(g in 1:nag){
  thisCode<-ageGroupsDF$Code[g]
  thisNumCohorts<-ageGroupsDF$NumCohorts[g]
  thisMort<-altMortDF[g,,]
  thisBaseM<-B0df$M[B0df$Code==thisCode]
  thisMax<-max(max(thisMort,na.rm=TRUE),thisBaseM,na.rm=TRUE)
  thisMin<-min(thisMort,na.rm=TRUE)
  if(!is.na(thisMax)){
    #get recruitment method
    yy<-paste("flagrecruit",thisCode,sep=""); yyy<-biolLines[grep(yy,biolLines)]
    thisFlagRecruit<-get_first_number(yyy)
    thisFlagDescr<-ifelse(thisFlagRecruit==12,"Constant","BH")
    plot(x=seq(1,nts),y=rep(thisMax,nts),type="n",xlab="",ylab="M",xaxt="n",ylim=c(min(0,thisMin),thisMax))
    for(c in 1:thisNumCohorts){
      points(thisMort[,c],type="l",lwd=2,col=cohortCols[c])
     }
    axis(at=xLabsAt,labels=xLabs,side=1)
    mtext(thisCode,side=3,adj=0)
    mtext(thisFlagDescr,side=3,adj=1)
    abline(h=thisBaseM,col=myGrey,lwd=2,lty=2)
  }
}
dev.off()

g=2
plot(thisNumsDF[g,,1],type="l",lwd=2,col=cohortCols[1],ylim=c(0,max(thisNumsDF[g,,1])))
points(thisNumsDF[g,,1]/3,type="l",lwd=2,col=cohortCols[1],lty=2)
for(c in 2:10){
  points(thisNumsDF[g,,c],type="l",lwd=2,col=cohortCols[c])
}

##read in outputMort to compare
# outputMort<-readLines(paste(outPath,"outputSpecificPredMort.txt",sep=""))
outputMort<-read.csv(paste(outPath,"outputSpecificPredMort.txt",sep=""),sep=" ")

test<-colSums(outputMort[,5:dim(outputMort)[2]])
index<-test>50
index<-test<1
toPlot<-test[index]
plot(toPlot,type="h",xaxt="n",xlab="")
par(las=2)
axis(at=seq(1,length(toPlot)),labels=names(toPlot),side=1)


