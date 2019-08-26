#read in outputs from a model, then plot together with re-estimated growth curves
this_run<-"TBGBMM"

this_out<-""

basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
# outPath<-paste(basePath,this_run,"\\","output",this_out,"_h",thisSteepness,"\\",sep="")
outPath<-paste(basePath,this_run,"\\","output",this_out,"\\",sep="")


plotPath<-paste(basePath,"\\Figures\\Recruitment\\TBGB",sep="")


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

burnin<-0 #number of years to skip in plot

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
nts<-dim(thisVol)[3]


colByTime<-colorRampPalette((colors=c(myYellow,myGreen,myAqua,myBlue,"midnightblue")))(nts)
formatCol<-function(x){
  y<-as.character(gsub("X\\.","#",x))
  return(y)
}


groupsDF<-read.csv(paste(outPath,"..\\TBGB_groups.csv",sep="")); ng<-dim(groupsDF)[1]

biolLines<-readLines(paste(outPath,"..\\TBGB_biol.prm",sep=""))

storeAllGs<-array(NA,dim=c(ng,10)); storeAllMs<-rep(NA,ng)

for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[g],side="both")
  cat(as.character(thisCode),"--")
  thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    
    weight_array<-array(NA,dim=c(thisNumCohorts,nts)); numbers_array<-weight_array
    
    for(c in 1:thisNumCohorts){
      thisTracer<-paste(thisName,c,"_ResN",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer);   thisRN<-apply(thisTemp,3,nonZeroMean)
      thisTracer<-paste(thisName,c,"_StructN",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer); thisSN<-apply(thisTemp,3,nonZeroMean)
      weight_array[c,]<-(thisRN+thisSN)*mg_2_tonne*X_CN
    }
    colnames(weight_array)<-colByTime
    tempArray<-data.frame(weight_array); tempArray$Cohort<-seq(1,thisNumCohorts)
    plotWeights<-melt(tempArray, id.var="Cohort")
    plotWeights$Color<-unlist(lapply(plotWeights$variable, FUN=formatCol))
    
    
    ## do array version
    storeGs<-array(NA, dim=c(thisNumCohorts,nyears))
    storeNumbers<-array(NA,dim=c(thisNumCohorts, nyears)); storeWeights<-storeNumbers 
    calcM<-rep(NA,thisNumCohorts)
    #fill weights and numbers from initial conditions (or first timestep from a run)
    for(c in 1:thisNumCohorts){
      thisTracer<-paste(thisName,c,"_ResN",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer);   thisRN<-apply(thisTemp,3,nonZeroMean)[1]
      thisTracer<-paste(thisName,c,"_StructN",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer); thisSN<-apply(thisTemp,3,nonZeroMean)[1]
      storeWeights[c,1]<-(thisRN+thisSN)
      ## numbers
      thisTracer<-paste(thisName,c,"_Nums",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer);
      storeNumbers[c,1]<-sum(thisTemp[,,1],na.rm=TRUE)
      numbers_array[c,]<-apply(thisTemp,3,sum,na.rm=TRUE)
      if(c>1){
        calcM[c]<-(1/ageClassSize)*(-1)*log(storeNumbers[c,1]/storeNumbers[c-1,1])
      }
    }
    
    colnames(numbers_array)<-colByTime
    tempArray<-data.frame(numbers_array); tempArray$Cohort<-seq(1,thisNumCohorts)
    plotTracerNumbers<-melt(tempArray, id.var="Cohort")
    plotTracerNumbers$Color<-unlist(lapply(plotTracerNumbers$variable, FUN=formatCol))
    
    
    
    ## fit to initial conditions
    #### fit curve y=x*exp(-M*ageclasssize) to final timestep tracer numbers and est realised M
    numsICFromM_fn<-function(M){
      thisNumFromM<-rep(NA,thisNumCohorts); thisNumFromM[1]<-numbers_array[1,1]
      for(i in 2:thisNumCohorts){
        thisNumFromM[i]<-(thisNumFromM[i-1])*exp(-M*ageClassSize)
      }
      return(thisNumFromM)
    }
    ll_IC_fn<-function(pars){
      optimM<-pars[1]
      estNums<-numsICFromM_fn(optimM)
      ss<-0
      for(c in 1:thisNumCohorts){
        ss<-ss+abs(numbers_array[c,1]-estNums[c])
      }
      return(ss)
    }
    
    
    testM<-mean(calcM,na.rm=TRUE); storeAllMs[g]<-testM
    fittedInitialM<-optim(par=testM,fn=ll_IC_fn, method="Brent", upper=2, lower=0)$par
    fittedInitialNums<-numsICFromM_fn(fittedInitialM)
    
    
    # NumsByM<-rep(NA,thisNumCohorts); NumsByM[1]<-storeNumbers[1]
    # for(c in 2:thisNumCohorts){
    #   NumsByM[c]<-NumsByM[c-1]*exp(-(testM*ageClassSize))
    # }
    
    
    #### fit curve y=x*exp(-M*ageclasssize) to final timestep tracer numbers and est realised M
    numsFromM_fn<-function(M){
      thisNumFromM<-rep(NA,thisNumCohorts); thisNumFromM[1]<-numbers_array[1,nts]
      for(i in 2:thisNumCohorts){
        thisNumFromM[i]<-(thisNumFromM[i-1])*exp(-M*ageClassSize)
      }
      return(thisNumFromM)
    }
    ll_fn<-function(pars){
      optimM<-pars[1]
      estNums<-numsFromM_fn(optimM)
      ss<-0
      for(c in 1:thisNumCohorts){
        ss<-ss+abs(numbers_array[c,nts]-estNums[c])
      }
      return(ss)
    }
    ##fit to final timestep
    fittedM<-optim(par=0.3,fn=ll_fn, method="Brent", upper=2, lower=0)$par
    fittedNums<-numsFromM_fn(fittedM)
    
    # testM<-0.12
    ##need to calculate growth for first cohort. 
    kwrr<-get_KWRR(thisCode); kwsr<-get_KWSR(thisCode)
    recWeight<-(kwsr+kwrr)
    ## starting G for cohort 1 doesn't take into account M, but OK place to start
    G1<-(log(storeWeights[2,1]) - log(storeWeights[1,1]))/(ageClassSize*365/2)
    G1=max(G1,0)
    G1<-1e-6
    
    for(t in 2:nyears){
      for(c in 1:(thisNumCohorts-1)){
        if(c>1){
          #second cohort and higher
          thisG<-calcG(Wi=storeWeights[c,1]  ,Wex=storeWeights[c,t-1], Nex=storeNumbers[c,t-1], Wrec=storeWeights[c-1,t], Nrec=storeNumbers[c-1,t]/ageClassSize, A=ageClassSize, M=testM)
          storeGs[c,t]<-thisG
          recWeight<-storeWeights[c-1,t]; recNums<-storeNumbers[c-1,t]
        } else{
          thisG<-G1
          recWeight<-storeWeights[1,1]; recNums<-storeNumbers[1,1]
        }
        newWeight<-storeWeights[c,t-1]*exp(thisG*365)
        newNums<-(storeNumbers[c,t-1]*(ageClassSize-1)/ageClassSize)*exp(-testM)
        totalNewWeight<-(recWeight*recNums) + (newWeight*newNums)
        totalWeightPerInd<-totalNewWeight/(newNums+recNums)
        storeWeights[c,t]<-totalWeightPerInd; storeNumbers[c,t]<-(newNums+recNums)
        
      }
    }
    storeAllGs[g,1:thisNumCohorts]<-storeGs[,nyears]
    
    thisCols<-colorRampPalette(colors=c(myGrey_trans,"black"))(nyears)
    
    jpeg(paste(plotPath,"GrowthCurve_",thisCode,"_re_estimated.jpg",sep=""),quality=300, width=800, height = 400)
    par(mfrow=c(1,2), mar=c(4,4,1,1), oma=c(1,1,1,1))
    
    plot(x=plotWeights$Cohort,y=plotWeights$value,col=plotWeights$Color,pch=20,cex=1.5,ylim=c(0,max(plotWeights$value)))
    points(x=seq(1,thisNumCohorts),y=weight_array[,1],col=myOrange,type="l",lwd=2)
    points(x=seq(1,thisNumCohorts),y=weight_array[,nts],col="midnightblue",type="l",lwd=2)
    for(c in 1:thisNumCohorts){
      points(x=rep(c,nyears), y=storeWeights[c,]*(mg_2_tonne*X_CN), pch=4, col=thisCols)
      points(x=c,y=storeWeights[c,nyears]*(mg_2_tonne*X_CN),pch=8,col="red")
    }
    plot(x=plotTracerNumbers$Cohort,y=plotTracerNumbers$value,col=plotTracerNumbers$Color,pch=20,cex=1.5,ylim=c(0,max(plotTracerNumbers$value)))
    #add in what the numbers should be based on M
    points(x=seq(1,thisNumCohorts),y=fittedInitialNums,type="l",col="red",lty=1)
    points(x=seq(1,thisNumCohorts),y=fittedNums,type="l",col=myGrey)
    mtext(paste("M initial ",signif(testM,2)),side=3,adj=1,line=1,col="red")
    mtext(paste("M fitted ",signif(fittedM,2)),side=3,adj=1,line=0,col=myGrey)
    dev.off()  
  }
  
}
