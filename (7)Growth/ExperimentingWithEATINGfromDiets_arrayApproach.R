source(paste(DIR$'General functions',"get_first_number.R",sep=""))
source(paste(DIR$'General functions',"myRounding.R",sep=""))
source(paste(DIR$'General functions',"getSNRN.R",sep=""))
source(paste(DIR$'General functions',"nonZeroMean.R",sep=""))
source(paste(DIR$'General functions',"getICNumbers.R",sep=""))
source(paste(DIR$'General functions',"getRequiredFoodPerDay_bySpace.R",sep=""))

this_run<-"Base"

# this_path<-paste(DIR$'Base',"\\ATLANTISmodels\\",sep="")
this_path<-paste(DIR$'Base',"\\ATLANTISmodels\\",sep="")
plotPath<-paste(this_path,"Figures\\",this_run,"\\",sep="")

##we'll need the groups df 
groupsDF<-read.csv(paste(this_path,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
thisInitialConditionsFile<-paste(this_path,"CRAM_input_short.nc",sep="")
ThisIC.nc<-nc_open(thisInitialConditionsFile)
thisBiolFile<-paste(this_path,"CRAM_BH_hybrid_biol.prm",sep="")
biolLines<-readLines(thisBiolFile)
groupsTL<-read.csv(paste(this_path,"inputs\\biol_prm\\CRAM_trophicLevels_isotopes.csv",sep=""))

#get volume from ICs
volume<-ncvar_get(ThisIC.nc,"volume")[,,1]
nlayers<-dim(volume)[1]; nboxes<-dim(volume)[2]

ppGroups<-as.character(read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ppGroups.csv",sep=""))[,1]);
npg<-length(ppGroups)

load(paste(DIR$"Data","\\eating\\interaction_spatial_array",sep="")); ##to bring in interaction_spatial_array. dim = npg, npg, nlayers, nboxes
load(paste(DIR$"Data","\\eating\\interaction_gapeSize_array",sep="")); ## dim = ng, 10, ng, 10
load(paste(DIR$"Data","\\eating\\interaction_trophicLevel_array",sep="")); ##dim = npg, npg

# sampledAvails<-data.matrix(read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ExpectedValues.csv",sep="")))
sampledAvails<-data.matrix(read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ExpectedValuesScaledForMax.csv",sep="")))[,-1]
# sampledAvails<-data.matrix(read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ExpectedValues_scaledUp.csv",sep="")))

########################################
##turn trophic level and spatial interactions to by cohort - much easier to use for eating!
trophicLevelByCohort<-array(0,dim=c(ng,10,ng,10)); spatialByCohort<-array(NA,dim=c(ng,10,ng,10,nlayers,nboxes))
sampledAvailsByCohort<-array(0,dim=c(ng,10,ng,10));
ageCohortLinking<-data.frame(matrix(NA,nrow=ng,ncol=10))
rownames(ageCohortLinking)<-as.character(groupsDF$Code); colnames(ageCohortLinking)<-seq(1,10)
for(g in 1:npg){
  thisVar<-ppGroups[g]
  thisCode<-gsub('juv|ad',"",thisVar)
  thisAge<-gsub(thisCode,"",thisVar)
  thisCodeIndex<-groupsDF$Code==thisCode
  thisNumCohorts<-groupsDF$NumCohorts[thisCodeIndex]
  if(thisNumCohorts==1){
    ageCohortLinking[thisCodeIndex,1]<-thisVar
  } else{
    #get age mature
    thisMat<-get_first_number(biolLines[grep(paste(thisCode,"_age_mat",sep=""),biolLines)]) 
    if(thisAge=="juv"){
      ageCohortLinking[thisCodeIndex,1:thisMat]<-thisVar
    }else{
      ageCohortLinking[thisCodeIndex,(thisMat+1):thisNumCohorts]<-thisVar
    }
  }
}

write.csv(ageCohortLinking,paste(DIR$'Tables',"SampledAvailabilities\\ageCohortLinking.csv",sep=""),row.names = TRUE)

for(g in 1:ng){
  thisPredCode<-groupsDF$Code[g]; thisPredRow<-ageCohortLinking[c(thisPredCode),]
  for(h in 1:ng){
    thisPreyCode<-groupsDF$Code[h]; thisPreyRow<-ageCohortLinking[c(thisPreyCode),]
    trophicLevelByCohort[g,,h,]<-interaction_trophicLevel_array[match(thisPredRow,ppGroups),match(thisPreyRow,ppGroups)]
    spatialByCohort[g,,h,,,]<-interaction_spatial_array[match(thisPredRow,ppGroups),match(thisPreyRow,ppGroups),,]
    sampledAvailsByCohort[g,,h,]<-sampledAvails[match(thisPredRow,ppGroups),match(thisPreyRow,ppGroups)]
  }
}

##############################
## will need biomass and numbers (the latter only for age-structured) #########
biomassIC<-array(NA,dim=c(ng,10,nlayers,nboxes))
NumbersIC<-array(NA,dim=c(ng,10,nlayers,nboxes))
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]
  thisName<-str_trim(groupsDF$Name[g],side="both")
  thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts==1){
    thisVar<-paste(thisName,"_N",sep="")
    temp<-ncvar_get(ThisIC.nc,thisVar)
    if(length(dim(temp))==3){
      biomassIC[g,1,,]<-temp[,,1]
    }else{
      biomassIC[g,1,nlayers,]<-temp[,1]
    }
  } else{
    for(c in 1:thisNumCohorts){
      thisNumsVar<-paste(thisName,c,"_Nums",sep="")
      xx<-ncvar_get(ThisIC.nc,thisNumsVar)
      if(length(dim(xx))==3){xx<-xx[,,1]}
      NumbersIC[g,c,,]<-xx
      thisRN<-getSNRN(name=thisName,cohort=c,whichN = "RN",ThisIC.nc)
      thisSN<-getSNRN(name=thisName,cohort=c,whichN = "SN",ThisIC.nc)
      biomassIC[g,c,,]<-NumbersIC[g,c,,]*(thisRN+thisSN)
    }
  }
}

##write biomassIC out
save(list=c("biomassIC"),file=paste(DIR$"Data","\\eating\\biomassIC",sep=""))
# load(paste(DIR$"Data","\\eating\\biomassIC",sep="")); ##to bring it back


##############################
## try some eating #########
#total available food for a given predator is needed by cohort and space
#need to turn NAs in gape size interaction to zeros
NAtozero<-function(x){
  if(is.na(x)){x<-0}
  return(x)
}
interaction_gapeSize_noNA<-apply(interaction_gapeSize_array,c(1,2,3,4),FUN=NAtozero)
trophicLevelByCohort_noNA<-apply(trophicLevelByCohort,c(1,2,3,4),FUN=NAtozero)

biomassIC_noNA<-apply(biomassIC,c(1,2,3,4),FUN=NAtozero)
tensorBiomass<-to.tensor(biomassIC_noNA)
availFoodTotal<-array(NA,dim=c(ng,10,nlayers,nboxes))
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]
  thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts==1){
    tensorSAvail<-to.tensor(sampledAvailsByCohort[g,1,,])
    tensorSAvail<-to.tensor(apply(tensorSAvail,c(1,2),NAtozero))
    thisPrey<-tensorSAvail %e% tensorBiomass
    availFoodTotal[g,1,,]<-thisPrey
  }else{
    for(c in 1:thisNumCohorts){
      tensorGape<-to.tensor(sampledAvailsByCohort[g,c,,] * interaction_gapeSize_noNA[g,c,,])
      tensorGape<-to.tensor(apply(tensorGape,c(1,2),NAtozero))
      availFoodTotal[g,c,,]<-mul.tensor(X=tensorGape,i=c(1,2),Y=tensorBiomass,j=c(1,2))
    }
  }
}

#need growth for each pred to get prop of avail eaten
reqFoodByGroupCohort<-array(NA,dim=c(ng,10,dim(volume)))
for(g in 1:ng){
  groupCode<-groupsDF$Code[g]
  thisRF<-getRequiredFoodPerDay_bySpace(groupCode,groupsDF,biolLines,ThisIC.nc)
  reqFoodByGroupCohort[g,,,]<-thisRF
}

 zerotoNA<-function(x){
   if(!is.na(x)){
     if(x==0){
       x<-NA
     }
   }
   return(x)
 } 

availFoodTotal<-apply(availFoodTotal,c(1,2,3,4),zerotoNA)
propEat<-reqFoodByGroupCohort/availFoodTotal
# may be worth checking doesn't go over 1 here
maxTo1<-function(x){
  x<-min(x,1)
  return(x)
}
naTo1<-function(x){
  if(is.na(x)){
    x<-1
  }
  return(x)
}
propEat<-apply(propEat,c(1,2,3,4),maxTo1)
propEat<-apply(propEat,c(1,2,3,4),naTo1)

##############################
## do actual eating
actualGrowth<-array(NA,dim=c(ng,10,nlayers,nboxes))
actualMort<-array(0,dim=c(ng,10,nlayers,nboxes))
actualMortByPred<-array(0,dim=c(ng,10,ng,10,nlayers,nboxes))


for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
  cat(as.character(thisCode,"--"))
  thisE<-get_first_number(biolLines[grep(paste("E_",thisCode,sep=""),biolLines)])
  if(length(thisE)>0){
      for(c in 1:thisNumCohorts){
       #sort prop eat out
       thisPropEat<-propEat[g,c,,]
       tensorPropEat<-to.tensor(thisPropEat)
       tempFilled<-rep(rep(tensorPropEat,ng,name="group"),10,name="cohort")
       filledPropEat<-reorder(tempFilled,i=c("group","cohort","I1","I2"))
       #sort samples out
       thisSamples<-to.tensor(apply(sampledAvailsByCohort[g,1,,],c(1,2),NAtozero))
       tempFilled<-rep(rep(thisSamples,nlayers,name="layer"),nboxes,name="box")
       filledSAvails<-reorder(tempFilled,i=c("I1","I2","layer","box"))
       
       # need 
       FilledSavaila_array<-array(eval(filledSAvails),dim=c(ng,10,nlayers,nboxes))
       FilledPropEat_array<-array(eval(filledPropEat),dim=c(ng,10,nlayers,nboxes))
       thisEat<-FilledSavaila_array*FilledPropEat_array*biomassIC
       
       actualGrowth[g,c,,]<-apply(thisE*thisEat,c(3,4),sum,na.rm=TRUE)
       
       actualMort<-actualMort+thisEat
       
       actualMortByPred[g,c,,,,]<-thisEat
     }
  }
}

## mortality needs to be turned into numbers for age-structured
weightByInd<-array(NA,dim=c(ng,10))
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    thisName<-str_trim(groupsDF$Name[g],side="both")
    for(c in 1:thisNumCohorts){
      thisSN<-getSNRN(name=thisName,cohort=c,whichN="SN",ThisIC.nc)
      thisRn<-getSNRN(name=thisName,cohort=c,whichN="RN",ThisIC.nc)
      weightByInd[g,c]<-thisSN+thisRN
    }
  }
}
tempFilled<-rep(rep(to.tensor(weightByInd),nlayers,name="layer"),nboxes,name="box")
filledWeightByInd<-reorder(tempFilled,i=c("I1","I2","layer","box"))
weightByInd_array<-array(eval(filledWeightByInd),dim=c(ng,10,nlayers,nboxes))

actualMortNumbers<-actualMort/weightByInd_array

calcM<-function(mday){
  M<-(-1)*365*log(1-mday,base=exp(1))
  return(M)
}

NumbersIC<-apply(NumbersIC,c(1,2,3,4),zerotoNA)
dailyPropM<-actualMortNumbers/NumbersIC

MfromPredation<-apply(dailyPropM,c(1,2,3,4),calcM)

g<-28
thisMs<-MfromPredation[g,,,]

meanMprops<-apply(dailyPropM,c(1),mean,na.rm=TRUE)
meanMprops[meanMprops=="NaN"]<-NA
toPlot<-as.double(meanMprops)[!is.na(meanMprops)]

pdf(paste(DIR$'Figures',"Eat\\meanPropMortPerDay.pdf",sep=""),height=5)
par(las=2)
plot(toPlot,col=myOrange,lwd=10,type="h",lend=1,xaxt="n",ylim=c(0,5))
axis(labels=as.character(groupsDF$Code[!is.na(meanMprops)]),at=seq(1,length(toPlot)),side=1)
dev.off()

testGroups<-c("PFS","IVS","IVH")
meanMprops<-apply(dailyPropM,c(1,2),mean,na.rm=TRUE)
meanMprops<-meanMprops[groupsDF$Code %in% testGroups,]

par(mfrow=c(3,1),mar=c(3,4.5,1,1),las=1)
for(t in 1:3){
  plot(meanMprops[t,],col=myOrange,lwd=10,type="h",lend=1,ylab="")
  mtext(testGroups[t],side=3)
}

meanM<-apply(MfromPredation,c(1,2),mean,na.rm=TRUE)
# meanM<-meanM[groupsDF$Code %in% testGroups,]

actualMs<-c(0.23, 0.036, 0.045)
pdf(paste(DIR$'Figures',"Eat\\meanMbyGroup.pdf",sep=""))
par(mfrow=c(3,1),mar=c(3,4.5,1,1),las=1)
for(t in 1:dim(meanM)[1]){
  if(sum(meanM[t,],na.rm=TRUE)>0){
    thisCol=myBlue
    if(max(meanM[t,],na.rm=TRUE)>1){thisCol<-"red"}
    plot(meanM[t,],col=thisCol,lwd=10,type="h",lend=1,ylab="",ylim=c(0,max(meanM[t,],na.rm=TRUE)))
    mtext(groupsDF$Code[t],side=3)
  }
  # abline(h=actualMs[t],col="red",lty=2)
}
dev.off()

#prop of required growth that we get
reqFoodByGroupCohort<-apply(reqFoodByGroupCohort,c(1,2,3,4),zerotoNA)
propGrowthMet<-actualGrowth/reqFoodByGroupCohort
meanGprops<-apply(propGrowthMet,c(1),mean,na.rm=TRUE)
meanGprops[meanGprops=="NaN"]<-NA
toPlot<-as.double(meanGprops)[!is.na(meanGprops)]
groupAxis<-as.character(groupsDF$Code)

# PFS 4.82192E+11	99319560270	20457356179
# 0.230
# IVS 330649344	11869454.09	9717887.082
# -log(1-(11869454.09/330649344),base=exp(1))
# 0.036
# #IVH 70582275	3118134.61	2309970.934
# -log(1-(3118134.61/70582275),base=exp(1))
# 0.045

# pdf(paste(DIR$'Figures',"Eat\\testingGrowthPropMean.pdf",sep=""),height=4)
pdf(paste(DIR$'Figures',"Eat\\testingGrowthPropMean_scaledUp.pdf",sep=""),height=4)
par(mar=c(4,4,0,0),c(0,0,0,0))
par(las=2)
plot(toPlot,col=myOrange,lwd=10,type="h",lend=1,xaxt="n",ylab="Mean proportion growth met")
abline(h=1,col="red",lwd=2,lty=2)
axis(labels=groupAxis,at=seq(1,length(groupAxis)),side=1)
dev.off()

# pdf(paste(DIR$'Figures',"Eat\\testingGrowthPropMean_zoomin.pdf",sep=""),height=4)
pdf(paste(DIR$'Figures',"Eat\\testingGrowthPropMean_zoomin_scaledUp.pdf",sep=""),height=4)
par(mar=c(4,4,0,0),c(0,0,0,0))
par(las=2)
plot(toPlot,col=myOrange,lwd=10,type="h",lend=1,xaxt="n",ylab="Mean proportion growth met",ylim=c(0,1))
abline(h=1,col="red",lwd=2,lty=2)
axis(labels=groupAxis,at=seq(1,length(groupAxis)),side=1)
dev.off()

temp<-apply(propGrowthMet,c(1,2),mean,na.rm=TRUE)
testGprops<-temp[groupsDF$Code %in% testGroups,]

minGroups<-apply(propGrowthMet,c(1),min,na.rm=TRUE)
minGroups[minGroups=="NaN"]<-NA
index<-!is.na(minGroups) & minGroups != "Inf"
toPlot<-as.double(minGroups)[index]
groupAxis<-as.character(groupsDF$Code)[index]


# pdf(paste(DIR$'Figures',"Eat\\testingGrowthPropMin_zoomin.pdf",sep=""),height=4)
pdf(paste(DIR$'Figures',"Eat\\testingGrowthPropMin_zoomin_scaledUp.pdf",sep=""),height=4)
par(mar=c(4,4,0,0),c(0,0,0,0))
par(las=2)
plot(toPlot,col=myOrange,lwd=10,type="h",lend=1,xaxt="n",ylab="Min proportion growth met",ylim=c(0,1))
abline(h=1,col="red",lwd=2,lty=2)
axis(labels=groupAxis,at=seq(1,length(groupAxis)),side=1)
dev.off()

## max

maxGroups<-apply(propGrowthMet,c(1),max,na.rm=TRUE)
maxGroups[maxGroups=="NaN"]<-NA
index<-!is.na(maxGroups) & maxGroups != "Inf" & maxGroups != "-Inf"
toPlot<-as.double(maxGroups)[index]
groupAxis<-as.character(groupsDF$Code)[index]


# pdf(paste(DIR$'Figures',"Eat\\testingGrowthPropMax_zoomin.pdf",sep=""),height=4)
pdf(paste(DIR$'Figures',"Eat\\testingGrowthPropMax_zoomin_scaledUp.pdf",sep=""),height=4)
par(mar=c(4,4,0,0),c(0,0,0,0))
par(las=2)
plot(toPlot,col=myOrange,lwd=10,type="h",lend=1,xaxt="n",ylab="Max proportion growth met",ylim=c(0,10))
abline(h=1,col="red",lwd=2,lty=2)
axis(labels=groupAxis,at=seq(1,length(groupAxis)),side=1)
dev.off()

A=mean(reqFoodByGroupCohort,na.rm=TRUE)
B=mean(availFoodTotal,na.rm=TRUE)
A/B
# [1] 0.0008148434
mean(propEat,na.rm=TRUE)
# [1] 0.7473358

testProp<-array(propEat[1,1,,])
testGrowReq<-array(reqFoodByGroupCohort[1,1,,])
testAvail<-array(availFoodTotal[1,1,,])
modData<-data.frame(cbind(testGrowReq,testAvail))

mod<-lm(testGrowReq~testAvail,data=modData)
this_k<-summary(mod)$'coefficients'[2,1]
expGR<-this_k * testAvail
thisCode<-as.character(groupsDF$Code[g])

colorByCohort<-colorRampPalette(colors=c(myGold,myGreen,myAqua,myBlue,myRed,"red"))(10)

# pdf(paste(DIR$'Figures',"Eat\\testingGrowth_allGroups.pdf",sep=""),height=7)
# pdf(paste(DIR$'Figures',"Eat\\testingGrowth_allGroups_KLPsReturned.pdf",sep=""),height=7)
pdf(paste(DIR$'Figures',"Eat\\testingGrowth_allGroups_KLPsReturned_scaledUp.pdf",sep=""),height=7)
par(mar=c(4,6,2,15),las=1,mfrow=c(3,1))
for(g in 1:ng){
  par(xpd="TRUE",las=1)
  thisCode<-as.character(groupsDF$Code[g])
  thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    maxFoodReq<-max(reqFoodByGroupCohort[g,,,],na.rm=TRUE)
    maxFoodAvail<-max(availFoodTotal[g,,,],na.rm=TRUE)
    plot(1,type="n",ylim=c(0,maxFoodReq),xlim=c(0,maxFoodAvail),ylab="",xlab="")
    mtext("Available prey",side=1,line=3)
    par(las=0)
    mtext("Food required",side=2,line=4.5)
    mtext(thisCode,side=3,adj=0)
    store_k<-rep(0,thisNumCohorts)
    for(c in 1:thisNumCohorts){
      thisCol<-colorByCohort[c]
      testGrowReq<-array(reqFoodByGroupCohort[g,c,,])
      testAvail<-array(availFoodTotal[g,c,,])
      if(sum(testAvail,na.rm=TRUE)>0){
        if(sum(testGrowReq,na.rm=TRUE)>0){
          modData<-data.frame(cbind(testGrowReq,testAvail))
          
          mod<-lm(testGrowReq~testAvail,data=modData)
          this_k<-summary(mod)$'coefficients'[2,1]
          expGR<-this_k * testAvail
          par(xpd=FALSE)
          points(x=testAvail,y=testGrowReq,pch=20,xlab="",ylab="",col=thisCol)
          points(x=testAvail,y=expGR,col=thisCol,type="l",lty=2,lwd=2)
        } else{
          this_k<-0
        }
      }else{
        this_k<-NA
      }
       store_k[c]<-this_k
    }
    legendText<-paste("Cohort ",seq(1,thisNumCohorts),", k = ",signif(store_k,2))
    par(xpd=NA)
    legend(legend=legendText,pch=20,bty="n",x="right",col=colorByCohort[1:thisNumCohorts],inset=-0.4)
  }
}
dev.off()

notEnoughFood<-c()
for(g in 1:ng){
  thisCode<-as.character(groupsDF$Code[g])
  thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    store_k<-rep(0,thisNumCohorts)
    for(c in 1:thisNumCohorts){
      thisCol<-colorByCohort[c]
      testGrowReq<-array(reqFoodByGroupCohort[g,c,,])
      testAvail<-array(availFoodTotal[g,c,,])
      if(sum(testAvail,na.rm=TRUE)>0){
        if(sum(testGrowReq,na.rm=TRUE)>0){
          modData<-data.frame(cbind(testGrowReq,testAvail))
          
          mod<-lm(testGrowReq~testAvail,data=modData)
          this_k<-summary(mod)$'coefficients'[2,1]
          expGR<-this_k * testAvail
        } else{
          this_k<-0
        }
      }else{
        this_k<-NA
      }
      store_k[c]<-this_k
    }
    if(max(store_k,na.rm=TRUE)>(10/365)){notEnoughFood<-c(notEnoughFood,thisCode)}
  }
}
notEnoughFood
#scaled up  "BEE" "IVH" "ORH" "PFS"
#returned KLP and KUPs 
#"BEE" "CRA" "GSH" "IVH" "IVS" "ORH" "PFS"
#previous
# [1] "ASQ" "BEE" "BIS" "CRA" "DPI" "EIS" "GSH" "HOK" "IVH" "IVS" "LDO" "LIN" "MAC" "MJE" "ORH" "PFM" "PFS" "RFI" "SPD" #1/365
# [1] "BEE" "CRA" "GSH" "HOK" "IVH" "IVS" "MJE" "ORH" "PFS" "SPD" # 10/365

##plot M
dim(MfromPredation)
meanMbyPreyC<-apply(MfromPredation,c(1,2),mean,na.rm=TRUE)
thisMaxM<-max(meanMbyPreyC,na.rm=TRUE)

# pdf(paste(DIR$'Figures',"Eat\\testingPredM_allGroups.pdf",sep=""),height=5)
# pdf(paste(DIR$'Figures',"Eat\\testingPredM_allGroups_KLPsReturned.pdf",sep=""),height=5)
pdf(paste(DIR$'Figures',"Eat\\testingPredM_allGroups_KLPsReturned_scaledUp.pdf",sep=""),height=5)
plot(seq(1,length(groupAxis)),ylim=c(0,thisMaxM),type="n",xlab="",ylab="M",xaxt="n")
abline(v=seq(,length(groupAxis)),col=myGrey_trans,lwd=1)
for(c in 1:10){
  thisCol<-colorByCohort[c]
  points(meanMbyPreyC[,c],col=thisCol,pch=20,cex=1.5)
}
par(las=2)
axis(labels=groupAxis,at=seq(1,length(groupAxis)),side=1)
dev.off()

# 
# pdf(paste(DIR$'Figures',"Eat\\testingPredM_allGroups_zoomin.pdf",sep=""),height=5)

# pdf(paste(DIR$'Figures',"Eat\\testingPredM_allGroups_zoomin_KLPsReturned.pdf",sep=""),height=5)
pdf(paste(DIR$'Figures',"Eat\\testingPredM_allGroups_zoomin_KLPsReturned_scaledUp.pdf",sep=""),height=5)
plot(seq(1,length(groupAxis)),ylim=c(0,2),type="n",xlab="",ylab="M",xaxt="n")
abline(v=seq(,length(groupAxis)),col=myGrey_trans,lwd=1)
for(c in 1:10){
  thisCol<-colorByCohort[c]
  points(meanMbyPreyC[,c],col=thisCol,pch=20,cex=1.5)
}
par(las=2)
axis(labels=groupAxis,at=seq(1,length(groupAxis)),side=1)
dev.off()

#########################################
## UP TO HERE UP TO HERE #############
##############################################
##finish checking out mortality

tooMuchMSometime<-groupsDF$Code[apply(MfromPredation,1,max,na.rm=TRUE)>2]

tooMuchMOften<-groupsDF$Code[apply(MfromPredation,1,mean,na.rm=TRUE)>2]
tooMuchMOften<-tooMuchMOften[!is.na(tooMuchMOften)]
# ASQ ELI ELP HAK LIN PFL
tooLittleMsometimes<-groupsDF$Code[apply(MfromPredation,1,min,na.rm=TRUE)<0.1]
# BEE BID DPI EID EIS GSH IVH IVS JAV LDO LIN MAC PFM PFS RFI SB  SND SPE


#####################################################
rescaleAvails<-array(1,dim=c(ng,ng))

##focus plots ## ##focus plots ## ##focus plots ## ##focus plots ##
focusGroup<-"SPE"
fg<-seq(1,ng)[groupsDF$Code==focusGroup]
temp<-actualMortByPred[,,fg,,,] #get mort by pred
# sum up over cells
mortByPred<-apply(temp,c(1,2,3),sum,na.rm=TRUE)
mortByPreyCohort<-apply(temp,c(3,4,5),sum,na.rm=TRUE)/weightByInd_array[fg,,,]
#to match, divide by numbers by cohort and box
mdayByPreyCohort<-mortByPreyCohort/NumbersIC[fg,,,]
#check out mort by cohort
mbycohort<-apply(mdayByPreyCohort,1,sum,na.rm=TRUE)
plot(mbycohort[mbycohort>0],type="h",lwd=10,col=myAqua,lend=1)
mbybox<-apply(mdayByPreyCohort,3,sum,na.rm=TRUE)
plot(mbybox,type="h",lwd=10,col=myAqua,lend=1)
mByPred<-apply(temp,c(1),sum,na.rm=TRUE)
plot(mByPred[mByPred>0],type="h",lwd=10,lend=1,col=myRed,xaxt="n",xlab="",ylab="Biomass predated")
npreds<-length(mByPred[mByPred>0]); axisAt<-seq(1,npreds)
axisLab<-as.character(groupsDF$Code[mByPred>0])
axis(at=axisAt,labels = axisLab,side=1)
thisDF<-data.frame(cbind(axisLab,signif(mByPred[mByPred>0],2)))
#just prey cohort 2
minPred<-0
minPred<-1e+8
mByPred<-apply(temp[,,c(1:3),,],c(1),sum,na.rm=TRUE)
plot(mByPred[mByPred>minPred],type="h",lwd=5,lend=1,col=myRed,xaxt="n",xlab="",ylab="Biomass predated")
npreds<-length(mByPred[mByPred>minPred]); axisAt<-seq(1,npreds)
axisLab<-as.character(groupsDF$Code[mByPred>minPred])
axis(at=axisAt,labels = axisLab,side=1)

##PFS PREDATORS
#mostly PFM in cohorts 3 and 4, by a long way - reduce this
#also reduce HOK, GSH, MJE, SPD, MAC, CRA predation
preyCode<-"SPE"; fg<-groupsDF$Code==preyCode
reduceCodes<-c("SPE","HOK")
fpgs<-groupsDF$Code %in% reduceCodes
rescaleAvails[fpgs,fg]<-0.01

preyCode<-"CRA"; fg<-groupsDF$Code==preyCode
reduceCodes<-c("GSH","HOK","MJE","SPD")
fpgs<-groupsDF$Code %in% reduceCodes
rescaleAvails[fpgs,fg]<-0.01

preyCode<-"RFI"; fg<-groupsDF$Code==preyCode
reduceCodes<-c("PFM","LDO","RFI")
fpgs<-groupsDF$Code %in% reduceCodes
rescaleAvails[fpgs,fg]<-0.01

preyCode<-"LDO"; fg<-groupsDF$Code==preyCode
reduceCodes<-c("LDO","SND","EID")
fpgs<-groupsDF$Code %in% reduceCodes
rescaleAvails[fpgs,fg]<-0.01

preyCode<-"LIN"; fg<-groupsDF$Code==preyCode
reduceCodes<-c("LIN","SND")
fpgs<-groupsDF$Code %in% reduceCodes
rescaleAvails[fpgs,fg]<-0.01

preyCode<-"JAV"; fg<-groupsDF$Code==preyCode
reduceCodes<-c("HOK","SND")
fpgs<-groupsDF$Code %in% reduceCodes
rescaleAvails[fpgs,fg]<-0.01

preyCode<-"PFS"; fg<-groupsDF$Code==preyCode
reduceCodes<-c("PFM", "HOK", "GSH", "MJE", "SPD", "MAC", "CRA")
fpgs<-groupsDF$Code %in% reduceCodes
rescaleAvails[fpgs,fg]<-0.01

preyCode<-"PFM"; fg<-groupsDF$Code==preyCode
reduceCodes<-c("PFM", "SND")
fpgs<-groupsDF$Code %in% reduceCodes
rescaleAvails[fpgs,fg]<-0.01

preyCode<-"BEE"; fg<-groupsDF$Code==preyCode
reduceCodes<-c("SND")
fpgs<-groupsDF$Code %in% reduceCodes
rescaleAvails[fpgs,fg]<-0.01

preyCode<-"BIS"; fg<-groupsDF$Code==preyCode
reduceCodes<-c("HOK","GSH")
fpgs<-groupsDF$Code %in% reduceCodes
rescaleAvails[fpgs,fg]<-0.01

preyCode<-"CBO"; fg<-groupsDF$Code==preyCode
reduceCodes<-c("HOK")
fpgs<-groupsDF$Code %in% reduceCodes
rescaleAvails[fpgs,fg]<-0.01

preyCode<-"DPI"; fg<-groupsDF$Code==preyCode
reduceCodes<-c("SND","EID")
fpgs<-groupsDF$Code %in% reduceCodes
rescaleAvails[fpgs,fg]<-0.01

preyCode<-"EIS"; fg<-groupsDF$Code==preyCode
reduceCodes<-c("PFM","HOK")
fpgs<-groupsDF$Code %in% reduceCodes
rescaleAvails[fpgs,fg]<-0.01

preyCode<-"IVS"; fg<-groupsDF$Code==preyCode
reduceCodes<-c("HOK")
fpgs<-groupsDF$Code %in% reduceCodes
rescaleAvails[fpgs,fg]<-0.01

#make new sampledAvails
sampledAvailsByCohort[g,,h,]<-sampledAvails[match(thisPredRow,ppGroups),match(thisPreyRow,ppGroups)]

newSampledAvails<-sampledAvails
for(g in 1:ng){
  thisPredCode<-groupsDF$Code[g]; 
  thisPredRow<-grep(thisPredCode,ppGroups)
  for(h in 1:ng){
    thisPreyCode<-groupsDF$Code[h]; thisPreyRow<-grep(thisPreyCode,ppGroups)
    thisScale<-rescaleAvails[g,h]
    newSampledAvails[thisPredRow,thisPreyRow]<-thisScale*sampledAvails[thisPredRow,thisPreyRow]

    if(thisPredCode==thisPreyCode){
      test<-newSampledAvails[thisPredRow,thisPreyRow]
      if(length(dim(test))>0){
        newSampledAvails[thisPredRow,thisPreyRow]<-apply(test,c(1,2),FUN=function(x){min(0.05,x)})
      }
    }
  }
}
#write it out
write.csv(newSampledAvails,paste(DIR$'Tables',"SampledAvailabilities\\ExpectedValues.csv",sep=""),row.names = FALSE)
