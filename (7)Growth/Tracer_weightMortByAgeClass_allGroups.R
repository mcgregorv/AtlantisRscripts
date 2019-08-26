#read in outputs from a model, then plot together with re-estimated growth curves
source(paste(DIR$'General functions',"nonZeroMean.R", sep=""))

##first plot here is summary of realised and intended mortality (M) for all groups. (blue and black bars)

##second plot has numbers in and out by cohort, but not sure this is quite right or useful

this_run<-"base"

this_out<-"MmumAll"
this_out<-"BHmumTEST150yr3base"
# this_out<-"TestMUMLong"

##read in M estimates from GrowthFromConstantM.R
storeInitialMs<-read.csv(paste(DIR$'Tables', "EstimatedMsFromICs.csv", sep=""))

mg_2_tonne<-2e-8; X_CN<-5.7
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

burnin<-0 #number of years to skip in plot

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")
r<-""

plotPath<-paste(this_path,"..\\Figures\\Mortality\\",this_out,"",sep="")

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
nts<-dim(thisVol)[3]


colByTime<-colorRampPalette((colors=c(myYellow,myGreen,myAqua,myBlue,"midnightblue")))(nts)
formatCol<-function(x){
  y<-as.character(gsub("X\\.","#",x))
  return(y)
}

groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

biolLines<-readLines(paste(this_path,"..\\CRAM_BH_hybrid_biol.prm",sep=""))

# nts<-72

storeDiffM<-array(NA, dim=c(ng, nts))

for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[g],side="both")
  cat(as.character(thisCode),"--")
  thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    #get ageclass size
    thisVar<-paste(thisCode, "_AgeClassSize", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageClassSize<-get_first_number(temp)
    testTimeStep<-ageClassSize*thisNumCohorts; testTimeStep<-min(nts, testTimeStep)
    
    weight_array<-array(NA,dim=c(thisNumCohorts,nts)); numbers_array<-weight_array
    for(c in 1:thisNumCohorts){
      thisTracer<-paste(thisName,c,"_ResN",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer);   thisRN<-apply(thisTemp,3,nonZeroMean)
      thisTracer<-paste(thisName,c,"_StructN",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer); thisSN<-apply(thisTemp,3,nonZeroMean)
      weight_array[c,]<-(thisRN+thisSN)*mg_2_tonne*X_CN
      #numbers
      thisTracer<-paste(thisName,c,"_Nums",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer);
      numbers_array[c,]<-apply(thisTemp,3,sum,na.rm=TRUE)
    }
   
    compareNumbers<-numbers_array[,testTimeStep]
    lastNumbers<-numbers_array[,nts]
    initialNumbers<-numbers_array[,1]
    
    numbersByTimeStep<-apply(numbers_array, 2, sum)
    
    ## fit to initial conditions
    #### fit curve y=x*exp(-M*ageclasssize) to final timestep tracer numbers and est realised M
    numsFromM_fn<-function(M){
      thisNumFromM<-rep(NA,thisNumCohorts); thisNumFromM[1]<-numbers_array[1,1]
      for(i in 2:thisNumCohorts){
        thisNumFromM[i]<-(thisNumFromM[i-1])*exp(-M*ageClassSize)
      }
      return(thisNumFromM)
    }
    ll_IC_fn<-function(pars){
      optimM<-pars[1]
      estNums<-numsFromM_fn(optimM)
      ss<-0
      for(c in 1:thisNumCohorts){
        ss<-ss+abs(obsNumbers[c]-estNums[c])
      }
      return(ss)
    }
    storeEstM<-rep(NA, nts)
    for(t in 1:nts){
      obsNumbers<-numbers_array[,t]
      storeEstM[t]<-optim(par=0.3,fn=ll_IC_fn, method="Brent", upper=2, lower=0)$par
    }
    diffM<-storeInitialMs[g,] - storeEstM
    
    #store diff in Ms
    storeDiffM[g,]<-diffM
   }
}

#plot final M difference by group
plot_nts<-nts
plot_nts<-50
index<-!is.na(storeDiffM[,nts])
thisYmax<-max(max(storeInitialMs[index,]), max(storeDiffM[,plot_nts][index])); thisYmin<-min(min(storeInitialMs[index,]), min(storeDiffM[,plot_nts][index]))
pdf(paste(plotPath, "MrealisedDiffALLgroups",plot_nts,".pdf", sep=""), width=10, height=4)
plot(storeDiffM[,plot_nts][index], type="h", lend=1, lwd=5, xaxt="n", ylim=c(thisYmin, thisYmax))
points(storeInitialMs[index,],type="h", lend=1, lwd=2, col=myBlue)
par(las=2)
axis(at=seq(1,length(seq(1,ng)[index])), as.character(groupsDF$Code[index]), side=1)
dev.off()


cbind(storeDiffM[,plot_nts], as.character(groupsDF$Code), storeInitialMs)

diffAs_mL<-unlist(lapply(storeDiffM[,nts], FUN=function(x){1-exp(-x/365)}))

##check out recruits coming in compared to those aging out
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[g],side="both")
  cat(as.character(thisCode),"--")
  thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    #get ageclass size
    thisVar<-paste(thisCode, "_AgeClassSize", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageClassSize<-get_first_number(temp)
    testTimeStep<-ageClassSize*thisNumCohorts; testTimeStep<-min(nts, testTimeStep)
    
    weight_array<-array(NA,dim=c(thisNumCohorts,nts)); numbers_array<-weight_array
    for(c in 1:thisNumCohorts){
      thisTracer<-paste(thisName,c,"_ResN",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer);   thisRN<-apply(thisTemp,3,nonZeroMean)
      thisTracer<-paste(thisName,c,"_StructN",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer); thisSN<-apply(thisTemp,3,nonZeroMean)
      weight_array[c,]<-(thisRN+thisSN)*mg_2_tonne*X_CN
      #numbers
      thisTracer<-paste(thisName,c,"_Nums",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer);
      numbers_array[c,]<-apply(thisTemp,3,sum,na.rm=TRUE)
    }
    SSB<-apply(weight_array * numbers_array, 2, sum)/(mg_2_tonne * X_CN)
    #get alpha and beta
    thisVar<-paste("BHalpha_", thisCode, sep=""); temp<-biolLines[grep(thisVar,biolLines)]; BHalpha<-get_first_number(temp)
    thisVar<-paste("BHbeta_", thisCode, sep=""); temp<-biolLines[grep(thisVar,biolLines)]; BHbeta<-get_first_number(temp)
    thisRec<-calcR(SSB=SSB[1], a=BHalpha, b=BHbeta)
    
    thisPlotFile<-paste(plotPath,"NumsInNumsOutByCohort_",thisCode,sep="")
    jpeg(paste(thisPlotFile, ".jpg", sep=""), quality=3000)
    par(mar=c(3,4.5,2,4.5))
    #number leaving
    plot(x=seq(1,thisNumCohorts), y=rep(thisRec, thisNumCohorts), type="n", ylim=c(0,1.2*thisRec))
    points(x=1, y=thisRec, type="h", lend=1, lwd=5); points(x=1, y=numbers_array[1,1]*(1/ageClassSize), type="h", lend=1, lwd=3, col=myOrange)
    for(c in 2:thisNumCohorts){
      nage_up<-numbers_array[c,1]*(1/ageClassSize); nageInto<-numbers_array[c-1,1]*(1/ageClassSize);
      points(x=c, y=nageInto, type="h", lend=1, lwd=5); points(x=c, y=nage_up, type="h", lend=1, lwd=3, col=myOrange)
    }
    c=2
    nage_up<-numbers_array[c,1]*(1/ageClassSize); nageInto<-numbers_array[c-1,1]*(1/ageClassSize);
    Mprop_add<-(nageInto-nage_up)/numbers_array[c,1]
    Madd<-(-1)*log(1-Mprop_add); mLadd<-signif(1-exp(-Madd/365),2)
    dev.off()
  }
    
}
# # saveDiffMs<-storeDiffM
# plot(storeDiffM[,plot_nts][index], type="h", lend=1, lwd=5, xaxt="n")
# points(storeInitialMs[index,],type="h", lend=1, lwd=2, col=myBlue)
# par(las=2)
# axis(at=seq(1,length(seq(1,ng)[index])), as.character(groupsDF$Code[index]), side=1)
# 
# 
