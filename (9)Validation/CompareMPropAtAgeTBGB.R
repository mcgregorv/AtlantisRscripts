##Read in tracers and summarise how much mortality each group has (estimated from numbers at age) and compare this with M from mL
mg_2_tonne<-2e-8; X_CN<-5.7

#################################### TRACERS
#plot all tracers for a given box and layer


this_out<-"BaseBurnin35"; runFolder<-"TBGB_JP2";
# this_out<-"output"; runFolder<-"TBGBReportBase"; 
# this_out<-"outputSCAlong4"; runFolder<-"TBGB_SI";
# this_out<-"MyRun_Fish1899_better1_codeupdate_rewriteDmatrix1"; runFolder<-"TBGBFish";
thisDesc <- paste(runFolder, this_out,sep="")

basePath<-  paste(DIR$'Base', "TBGB\\",runFolder,"\\",sep="")
outPath<-paste(basePath,"\\",this_out,"\\",sep="")
plotPath<-paste(basePath,"..\\Figures\\Testing\\PropAtAge\\",thisDesc, sep="")

burnin<-1 #number of years to skip in plot

groupsDF<-read.csv(paste(basePath, "\\TBGB_Groups.csv",sep="")); ng<-dim(groupsDF)[1]
ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))

doingMedian<-TRUE
# doingMedian<-FALSE
#read in 'obs' M
groupBioPars<-read.csv(paste(basePath,"..\\TBGB_B0.csv",sep=""))
groupBioPars$M <- as.double(as.character(groupBioPars$M))


biolLines <- readLines(paste(basePath,  "TBGB_biol.prm",sep=""))

daysTimeStep<-73
numStepsPerYear<-365/daysTimeStep
year0<-1865
fishingStartYear<-1899
modelStartYear<-1899

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

#grab _N tracers and store as tonnes for each
storeNumbersByGroup<-array(NA, dim=c(nts, ng, 10))
for(g in 1:ng){
  thisName<-str_trim(groupsDF$Name[g], side="both"); thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    for(c in 1:thisNumCohorts){
      thisTracer<-paste(thisName,c,"_Nums", sep="")
      temp<-ncvar_get(ThisNC.nc, thisTracer)
      thisNumbers<-apply(temp,3, sum)
      storeNumbersByGroup[,g,c]<-thisNumbers[(burnin+1):length(thisNumbers)]
    }
  }
}
##get proportions at age for each timestep and group, so can get 95% CIs
storeAllPropsAtAge<-array(NA,dim=c(nts,ng,10))
for(g in 1:ng){
  thisNumCohorts<-groupsDF$NumCohorts[g]
  thisNumByTS<-apply(storeNumbersByGroup[,g,],1,sum,na.rm=TRUE)
  if(thisNumCohorts>2){
    for(t in 1:nts){
      storeAllPropsAtAge[t,g,1:thisNumCohorts]<-storeNumbersByGroup[t,g,1:thisNumCohorts]/thisNumByTS[t]
    }
  }
}

numsICFromM_fn<-function(M){
  thisNumFromM<-rep(NA,length(thisNumbers)); thisNumFromM[1]<-thisNumbers[1]
  for(i in 2:length(thisNumbers)){
    thisNumFromM[i]<-(thisNumFromM[i-1])*exp(-M*ageClassSize)
  }
  return(thisNumFromM)
}
ll_IC_fn<-function(pars){
  optimM<-pars[1]
  estNums<-numsICFromM_fn(optimM)
  ss<-0
  for(c in 1:length(thisNumbers)){
    ss<-ss+abs(thisNumbers[c]-estNums[c])
  }
  return(ss)
}
t=nts
calc_P1<-function(M, thisNumCohorts, ageClassSize){
  #what should the first proportion be for all to add to 1, for a given number of age classes and a given M?
  xx=0
  for(i in 1:thisNumCohorts){
    xx = xx + exp(-(i-1) * M * ageClassSize)
  }
  p1<-1/xx
  return(p1)
}
getProportionsAtAge<-function(M, thisNumCohorts, ageClassSize){
  thisP1<-calc_P1(M, thisNumCohorts, ageClassSize)
  propsAtAge<-rep(NA, thisNumCohorts); propsAtAge[1]<-thisP1
  for(i in 2:thisNumCohorts){
    propsAtAge[i]<-propsAtAge[i-1]* exp(-(M * ageClassSize))
  }
  return(propsAtAge)
}
#now loop through and estimate M and create numbers at age plots
storeJuvMs<-rep(NA,ng); storeAdMs<-rep(NA,ng)
storeInitialJuvMs<-rep(NA, ng); storeInitialAdMs <- rep(NA, ng)
toPlot<-c()
plotFile<-paste(plotPath,"FittedM.pdf",sep="")
pdf(plotFile, height=8)
par(mfrow=c(3,2), mar=c(4,4,4,2))
for(g in 1:ng){
  thisName<-str_trim(groupsDF$Name[g], side="both"); thisNumCohorts<-groupsDF$NumCohorts[g]; thisCode<-as.character(groupsDF$Code[g])
  if(thisNumCohorts>2){

      thisVar<-paste(thisCode, "_AgeClassSize", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageClassSize<-get_first_number(temp)
      thisVar<-paste(thisCode, "_age_mat", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageMature<-get_first_number(temp)
      ## define props at age using M ('obs')
      thisM<-groupBioPars$M[groupBioPars$Code==thisCode]
      if(!is.na(thisM)){     # 'observed' props at age using M from lit
        obsPropsAtAge<-getProportionsAtAge(thisM, thisNumCohorts, ageClassSize)
      }
      thisData<-storeNumbersByGroup[,g,]
      totalNumsTracer<-apply(thisData,1,sum, na.rm=TRUE)

      thisLongName<-gsub("_"," ", groupsDF$Name[g])
      thisLongName<-gsub("-","\n",thisLongName)
      thisLongName<-gsub("\\(","\n",thisLongName); thisLongName<-gsub(")", "", thisLongName)
   
      
      propsByCohortTime<- 0 * thisData
      for(t in 1:dim(thisData)[1]){propsByCohortTime[t,]<-thisData[t,]/totalNumsTracer[t]}
      
      thisProps<-thisData[t,]/totalNumsTracer[t]
      # if(doingMedian==TRUE){thisProps<-propMedian}
      #do juveniles first
      thisNumbers<-thisProps[1:(ageMature+1)]
      juvFittedM<-optim(par=0.1,fn=ll_IC_fn, method="Brent", upper=2, lower=0)$par
      juvFittedNums<-numsICFromM_fn(juvFittedM)
      ## do adults
      thisNumbers<-thisProps[(ageMature+1):thisNumCohorts]
      adFittedM<-optim(par=0.1,fn=ll_IC_fn, method="Brent", upper=2, lower=0)$par
      adFittedNums<-numsICFromM_fn(adFittedM)
      
      thisYmax<- max(thisProps, na.rm=TRUE)*1.2
      
      if(thisYmax>0){
        
        toPlot<-c(toPlot, thisCode)
        # par(mar=c(4.5,4.5,3.5,0.5))
        bpdata<-melt(propsByCohortTime)
        boxplot(value ~ X2, data=bpdata, xlab="Age (years)", ylab="Proportion at age", xaxt="n",cex=thisCex, cex.axis=thisCex, cex.lab=thisCex,ylim=c(0,thisYmax), outline=FALSE)
         axis(at=seq(1,thisNumCohorts), labels=seq(1,thisNumCohorts)*ageClassSize, side=1, cex.axis=thisCex, cex.lab=thisCex)
        mtext(thisLongName,side=3, adj=0, cex=thisCex)
        abline(v=(ageMature+1), col=myGrey_trans, lwd=3)
        mtext(paste("adult M ", signif(adFittedM,2),"\nJuv M ",signif(juvFittedM,2),sep=""),side=3, adj=1, line=-2.5)
        if(!is.na(thisM)){     # 'observed' props at age using M from lit
          points(obsPropsAtAge, type="l", col=myBlue,lwd=2)
          mtext(paste("Obs M ", thisM, sep=""), side=3, col=myBlue,adj=1)
        }
        storeAdMs[g]<-adFittedM; storeJuvMs[g]<-juvFittedM
      } else{
        makeBlankPlot(); mtext(thisCode,side=3, adj=0)
      }
    # }
  }
  
}
dev.off()
