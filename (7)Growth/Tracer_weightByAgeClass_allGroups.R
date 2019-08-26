#read in outputs from a model, then plot together with re-estimated growth curves
this_run<-"base"

this_out<-"XXX_mLB5"
# this_out<-"SENSselBASE"
# this_out<-"TestMUMLong"

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

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")
r<-""

plotPath<-paste(this_path,"..\\Figures\\growth\\",this_out,"",sep="")

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
nts<-dim(thisVol)[3]
nyears<-nts


colByTime<-colorRampPalette((colors=c(myYellow,myGreen,myAqua,myBlue,"midnightblue")))(nts)
formatCol<-function(x){
  y<-as.character(gsub("X\\.","#",x))
  return(y)
}


groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

biolLines<-readLines(paste(this_path,"..\\CRAM_BH_hybrid_biol.prm",sep=""))

storeAllGs<-array(NA,dim=c(ng,10)); storeAllMs<-rep(NA,ng)

calcG<-function(Wi, Wex, Nex, Wrec, Nrec, A, M){
  Nnew<-Nex*((A-1)/A)*exp(-M)
  G<-(1/365)*log( (Wi*(Nnew + Nrec) - Wrec * Nrec) / (Wex * Nnew))
  return(G)
}

for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[g],side="both")
  cat(as.character(thisCode),"--")
  thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    
    weight_array<-array(NA,dim=c(thisNumCohorts,nts)); numbers_array<-weight_array
    #get ageclass size
    thisVar<-paste(thisCode, "_AgeClassSize", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageClassSize<-get_first_number(temp)
    
    
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
    calcM<-rep(NA,thisNumCohorts); calcMfinal<-calcM
    #fill weights and numbers from initial conditions (or first timestep from a run)
    for(c in 1:thisNumCohorts){
      thisTracer<-paste(thisName,c,"_ResN",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer);   thisRN<-apply(thisTemp,3,nonZeroMean)[1]
      thisTracer<-paste(thisName,c,"_StructN",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer); thisSN<-apply(thisTemp,3,nonZeroMean)[1]
      storeWeights[c,]<-(thisRN+thisSN)
      ## numbers
      thisTracer<-paste(thisName,c,"_Nums",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer);
      storeNumbers[c,]<-sum(thisTemp[,,1])
      numbers_array[c,]<-apply(thisTemp,3,sum,na.rm=TRUE)
      if(c>1){
        calcM[c]<-(1/ageClassSize)*(-1)*log(storeNumbers[c,1]/storeNumbers[c-1,1])
        calcMfinal[c]<-(1/ageClassSize)*(-1)*log(storeNumbers[c,nts]/storeNumbers[c-1,nts])
      }
    }
    
    colnames(numbers_array)<-colByTime
    tempArray<-data.frame(numbers_array); tempArray$Cohort<-seq(1,thisNumCohorts)
    plotTracerNumbers<-melt(tempArray, id.var="Cohort")
    plotTracerNumbers$Color<-unlist(lapply(plotTracerNumbers$variable, FUN=formatCol))
    
    
    numA<-numbers_array[1,1]; numB<-numbers_array[1,2]
    M1<-(-1)*log(numA/numB); M1day<-exp(-(M1/365))
    1-M1day
    
    numA<-numbers_array[2,1]; numB<-numbers_array[2,2]
    M1<-(-1)*log(numA/numB); M1day<-exp(-(M1/365))
    1-M1day
    ##BAL ran with mLj = 0.0008, then increased to 0.00257
    ##think latter will be too far
    #something not quite right about numbers by ageclass perhaps - for enough to age up (1/12) each year to keep the numbers in the next ageclass constant
    
    
    
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
    

    ##try taking into acount proportion aging up - every year 1/12 th of them age up (BAL)
    #have numbers_array from tracers
    ageingM<-rep(NA,thisNumCohorts)
    for(c in 2:thisNumCohorts){
      x<-(12*numbers_array[c,2])/(numbers_array[c-1,1] + 11*numbers_array[c,1])
      ageingM[c]<-(-1)*log(x)
    }
      
    
    agingDailyMl<-unlist(lapply(ageingM, FUN=function(x){signif(1-exp(-x/365),2)}))
    mean(agingDailyMl[3:10])
    
    totalNumbers<-apply(numbers_array,2,sum)
    
    plot(totalNumbers/totalNumbers[1],type="l")
    
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
    
    # #do the version of growth using Shane's fn 
    # calcGrowth<-function(x,out_a1=FALSE){
    #   #shanes version
    #   a1<-1-((log(lambda0,base=exp(1))-log(lambda1,base=exp(1)))/(log(maxWeight,base=exp(1))-log(minWeight,base=exp(1))))
    #   a0<-thisE*lambda0*x^(1-a1)
    #   mum<-a0*x^a1
    #   if(out_a1==TRUE){
    #     return(a1)
    #   }else{
    #     return(mum)
    #   }
    # }
    # 
    # lambda0<-0.3; lambda1<-0.1; thisE<-0.5
    # 
    # growthByCohort<-calcGrowth(storeWeights[,1])
    # growthByCohortProp<-growthByCohort/storeWeights[,1]+1
    # intrinsGrowthYear<-unlist(lapply(growthByCohortProp,FUN=function(x){365*ageClassSize*log(x)}))
    
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
    mtext(paste("Age class size: ", ageClassSize,sep=""),side=3,adj=0)
    plot(x=plotTracerNumbers$Cohort,y=plotTracerNumbers$value,col=plotTracerNumbers$Color,pch=20,cex=1.5,ylim=c(0,max(plotTracerNumbers$value)))
    #add in what the numbers should be based on M
    points(x=seq(1,thisNumCohorts),y=fittedInitialNums,type="l",col="red",lty=1)
    points(x=seq(1,thisNumCohorts),y=fittedNums,type="l",col=myGrey)
    mtext(paste("M initial ",signif(testM,2)),side=3,adj=1,line=1,col="red")
    mtext(paste("M fitted ",signif(fittedM,2)),side=3,adj=1,line=0,col=myGrey)
    dev.off()  
  }
       
}

# storeAllGs
# storeAllMs

##overwrite mum values in biolLines
#for each, multiple growth rate by weight (in mg N) per individual in each ageclass. do this bit first
growth_mgN<-0*storeAllGs
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]; thisName<-str_trim(groupsDF$Name[g])
  if(thisNumCohorts>2){
    temp<-exp(storeAllGs[g,]); temp[is.na(temp)]<-1; temp<-temp[1:thisNumCohorts]
    #get weight of individuals at initial conditions
    thisICweight<-rep(0,thisNumCohorts)
    for(c in 1:thisNumCohorts){
      thisTracer<-paste(thisName,c,"_ResN",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer);   thisRN<-apply(thisTemp,3,nonZeroMean)[1]
      thisTracer<-paste(thisName,c,"_StructN",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer); thisSN<-apply(thisTemp,3,nonZeroMean)[1]
      thisICweight[c]<-thisRN+thisSN
    }
    ##also need E to multiply by. just set to 0.5 for now as seems to be what it is set for for all
    # thisE<-0.5
    growthBy_mgN<-thisICweight * temp 
    growth_mgN[g,1:thisNumCohorts]<-growthBy_mgN
  }
}
# 
# # 
# for(g in 1:ng){
#   thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
#   if(thisNumCohorts>2){
#     thisVar<-paste("mum_", thisCode,sep=""); x<-grep(thisVar,biolLines)
#     curLine<-get_first_number(biolLines[x+1],n="all")
#     temp<-round(growth_mgN[g,],2); temp[is.na(temp)]<-1; temp<-temp[1:thisNumCohorts]
#     newLine<-paste(temp,collapse="\t")
#     biolLines[x+1]<-newLine
#   }
# }
# writeLines(biolLines,paste(this_path,"..\\CRAM_BH_hybrid_biol.prm",sep=""))
# 
# # 
# # 




