#read in outputs from a model, then plot together with re-estimated growth curves. use this for other linear rec groups too
this_run<-"base"

this_out<-"PreSENS2"
this_out<-"CETmum"
# this_out<-"TestMUMLong"

source(paste(DIR$'General functions',"numsFromBothMs_fn.R",sep=""))
source(paste(DIR$'General functions',"numsFromBothMsR0_fn.R",sep=""))

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


colByTime<-colorRampPalette((colors=c(myYellow,myGreen,myAqua,myBlue,"midnightblue")))(nts)
formatCol<-function(x){
  y<-as.character(gsub("X\\.","#",x))
  return(y)
}


groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

biolLines<-readLines(paste(this_path,"..\\CRAM_BH_hybrid_biol.prm",sep=""))

storeAllGs<-array(NA,dim=c(ng,10)); storeAllMs<-rep(NA,ng)

this_h<-0.85


rbind(seq(1,ng), as.character(groupsDF$Code))
# for(g in 1:ng){
### UPT TO HERE ELP ELP ELP
g=9
  thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[g],side="both")
  cat(as.character(thisCode),"--")
  thisNumCohorts<-groupsDF$NumCohorts[g]
  # if(thisNumCohorts>1){
  #get ageclass size
  thisVar<-paste(thisCode, "_AgeClassSize", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageClassSize<-get_first_number(temp)
    
    weight_array<-array(NA,dim=c(thisNumCohorts,nts)); numbers_array<-weight_array
    for(c in 1:thisNumCohorts){
      thisTracer<-paste(thisName,c,"_ResN",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer);   thisRN<-apply(thisTemp,3,nonZeroMean)
      thisTracer<-paste(thisName,c,"_StructN",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer); thisSN<-apply(thisTemp,3,nonZeroMean)
      weight_array[c,]<-(thisRN+thisSN)*mg_2_tonne*X_CN
      ##numbers
      thisVar<-paste(thisName,c,"_Nums",sep=""); thisTemp<-ncvar_get(ThisNC.nc, thisVar); thisNums<-apply(thisTemp,3,sum)
      numbers_array[c,]<-thisNums
    }
    colByCohort<-colorRampPalette(colors=c(myBlue, myAqua,myGreen,myGold,"red"))(10)
    plot(numbers_array[1,],type="l",ylim=c(0,max(numbers_array)))
    for(c in 2:thisNumCohorts){
      points(numbers_array[c,], col=colByCohort[c],type="l")
    }
    
    plot(weight_array[1,],type="l",ylim=c(0,max(weight_array)))
    for(c in 2:thisNumCohorts){
      points(weight_array[c,], col=colByCohort[c],type="l")
    }
  
    ##number of adults spawning
    ##actually, atlantis uses the whole population for BH10
    FSPB<-rep(1,thisNumCohorts)
    ##try with BH pars for recruitment
    calcR<-function(S,a,b){
      R<-(a*S)/(b+S)
      return(R)
    }
    calc_a<-function(h,R0){
      a<-(4*h*R0)/(5*h-1)
      return(a)
    }
    calc_b<-function(h,B0){
      b<-(B0*(1-h))/(5*h-1)
      return(b)
    }
    thisR0<-sum(numbers_array[,1])/ageClassSize # doesnt' take into account mortality within the ageclass, but OK as starting point
    thisB0<-sum(simNums_array[,1] * FSPB * weight_array[,1]) /(mg_2_tonne * X_CN)
    BHpars<-c(calc_a(this_h, thisR0), calc_b(this_h, thisB0))
    
    ########################################################################
    #do the aging and stuff in the function to optimise
    ##try optimising the R0 as well
     llinclR0_fn<-function(pars){
      optimMj<-pars[1]; optimMa<-pars[2]; optimR0<-pars[3]
      estNums<-numsFromBothMsR0_fn(optimMj, optimMa, optimR0)
      ss<-0
      for(c in 1:(thisNumCohorts)){
        ss<-ss+abs(numbers_array[c,1]-estNums[c])
      }
      return(ss)
     }
     thisOptim<-optim(par=c(0.1,0.1,thisR0),fn=llinclR0_fn)
     
     
     ll_fn<-function(pars){
       optimMj<-pars[1]; optimMa<-pars[2];
       estNums<-numsFromBothMs_fn(optimMj, optimMa)
       ss<-0
       for(c in 1:(thisNumCohorts)){
         ss<-ss+abs(numbers_array[c,1]-estNums[c])
       }
       return(ss)
     }
    
    thisOptim<-optim(par=c(0.01,0.1),fn=ll_fn)
    fittedBothM<-thisOptim$par
    fittedNums<-numsFromBothMs_fn(fittedBothM[1], fittedBothM[2])
    plot(numbers_array[,1],pch=20)
    points(fittedNums,type="l")
 
    BHpars<-c(calc_a(this_h,thisR0), calc_b(this_h,thisB0))
    
    signif(1-unlist(lapply(fittedBothM[1:2], FUN=function(x){exp(-x/365)})),2)
    
    #############################################################
    # now that have mortality and recruitment estimates, do growth estimates
    ########################################################
    
    
    ########################################################
 
    juvM<-fittedBothM[1]; adM<-fittedBothM[2]
    ##need to calculate growth for first cohort. 
    kwrr<-get_KWRR(thisCode); kwsr<-get_KWSR(thisCode)
    recWeight<-(kwsr+kwrr)
    ## starting G for cohort 1 doesn't take into account M, but OK place to start
    G1<-(log(weight_array[2,1]) - log(weight_array[1,1]))/(ageClassSize*365)
    # G1=max(G1,0)
    # G1<-1e-6
   
    
    storeWeights<-array(NA,dim=c(thisNumCohorts,nyears)); storeNumbers<-storeWeights

   mum_array<-rep(1,thisNumCohorts)
    
    weightsAndNumbers_fn<-function(mum_array, nyears=250){
      simNums_array<-array(NA,dim=c(thisNumCohorts, nyears)); simNums_array[,1]<-numbers_array[,1]
      simWeights_array<-simNums_array; simWeights_array[,1]<-weight_array[,1]
      for(y in 2:nyears) {
        for(c in 1:thisNumCohorts){
          thisM<-Ma; if(c<3){thisM<-juvM}
          if(c==1){
            #calc recruits
            thisSSB<-sum(simNums_array[, y-1] * FSPB * weight_array[,1]) / (mg_2_tonne * X_CN)
            thisRecruits<-calcR(S=thisSSB, a=BHpars[1], b=BHpars[2])
            thisRecWeightTotal<-(kwsr+kwrr)*thisRecruits * mg_2_tonne * X_CN
          } else{
            thisRecruits<-simNums_array[c-1, y-1]*((ageClassSize-1)/ageClassSize) #these are those that aged up from previous cohort
            thisRecWeightTotal<-(simWeights_array[c-1,y-1])*thisRecruits
          }
          # add the recruits, take out the aged, and take out the mortality
          thisAgeClassNums<-simNums_array[c,y-1]*((ageClassSize-1)/ageClassSize)
          simNums_array[c,y] <- (thisAgeClassNums + thisRecruits)*exp((-1)*(thisM*ageClassSize))
          thisAgeClassWeightTotal<-thisAgeClassNums * simWeights_array[c,y-1] * exp(mum_array[c])
          newAvWeight<-(thisAgeClassWeightTotal + thisRecWeightTotal)/(thisAgeClassNums + thisRecruits)
          simWeights_array[c,y]<-newAvWeight
          #if it's the last age class, kill any extra
          if(c==thisNumCohorts){
            simNums_array[c,y]<-simNums_array[c,1]
            simWeights_array[c,y]<-simWeights_array[c,1]
          }
          
        }
      }
      ##output the last timestep of numbers to compare with the initial conditions
      thisNumFromM<-simNums_array[,nyears]
      thisWeightsFromMum<-simWeights_array[,nyears]
      return(thisWeightsFromMum)
    }
    
    
    ll_fn<-function(pars){
      optimMum<-c(pars,0)
      estWeights<-weightsAndNumbers_fn(optimMum, nyears=150)
      ss<-0
      if(min(optimMum)<0){ss<-1e+6}
      for(c in 1:(thisNumCohorts)){
        ss<-ss+(weight_array[c,1]-estWeights[c])^2
      }
      return(ss)
    }
    
    initPars<-c(G1, 0.85, 0.2, rep(0.05,7))
    thisOptim<-optim(par=initPars,fn=ll_fn)
     optimMum<-thisOptim$par
    fittedWeights<-weightsAndNumbers_fn(optimMum)
    plot(weight_array[,1],fittedWeights)
    
    plot(fittedWeights,type="l", ylim=c(0,max(max(fittedWeights), max(weight_array[,1]))))
    points(weight_array[,1])
    
    #test with longer as well
    fittedWeights2<-weightsAndNumbers_fn(optimMum, nyears=250)
    points(fittedWeights2,type="l",col="red")
    fittedWeights3<-weightsAndNumbers_fn(optimMum, nyears=20)
    points(fittedWeights3,type="l",col="blue", lty=2)
    
    
    new_mum<-c(optimMum)
    
    signif(new_mum,4)
    
    
    #################################################
 
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
# 




