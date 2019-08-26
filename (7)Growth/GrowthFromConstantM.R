#read in tracer, fit and estimate M in initial conditions,  then fit growth curve using constant G over some power of ageclass
#read in outputs from a model, then plot together with re-estimated growth curves
source(paste(DIR$'General functions', "numsFromBothMs_fn.R", sep=""))

this_run<-"base"

this_out<-"PreSENS2"
this_out<-"SENSselBASE"
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

plotPath<-paste(this_path,"..\\Figures\\growth\\IC_",sep="")

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

storeAllGs<-array(NA,dim=c(ng,10)); storeAllMs<-rep(NA,ng); storeAllMums<-array(NA,dim=c(ng, 10))

for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[g],side="both")
  cat(as.character(thisCode),"--")
  thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    thisVar<-paste(thisCode, "_AgeClassSize", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageClassSize<-get_first_number(temp)
    
    
    weight_array<-array(NA,dim=c(thisNumCohorts,nts)); numbers_array<-weight_array
    
    for(c in 1:thisNumCohorts){
      thisTracer<-paste(thisName,c,"_ResN",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer);   thisRN<-apply(thisTemp,3,nonZeroMean)
      thisTracer<-paste(thisName,c,"_StructN",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer); thisSN<-apply(thisTemp,3,nonZeroMean)
      weight_array[c,]<-(thisRN+thisSN)*mg_2_tonne*X_CN
      thisTracer<-paste(thisName,c,"_Nums",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer);
      numbers_array[c,]<-apply(thisTemp, 3, sum)
    }
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
    
    fittedInitialM<-optim(par=0.1,fn=ll_IC_fn, method="Brent", upper=2, lower=0)$par
    fittedInitialNums<-numsICFromM_fn(fittedInitialM)
    jpeg(paste(plotPath,"fittedMortAndGrowth_", thisCode, ".jpg", sep=""))
    par(mfrow=c(1,2))
    plot(numbers_array[,1])
    points(fittedInitialNums,type="l")
    mtext(as.character(thisCode),side=3,adj=0)
    
    storeAllMs[g]<-fittedInitialM

    #now assume the numbers by age class are constant, and estimate growth by age class
    weightsICFromG_fn<-function(G, n){
      thisWeightFromG<-rep(NA,thisNumCohorts); thisWeightFromG[1]<-weight_array[1,1]
      for(i in 2:thisNumCohorts){
        thisWeightFromG[i]<-(thisWeightFromG[i-1])*exp((G/i^n)*ageClassSize)
      }
      return(thisWeightFromG)
    }
    ll_wIC_fn<-function(pars){
      optimG<-pars[1]; optim_n<-pars[2]
      estWeights<-weightsICFromG_fn(optimG, optim_n)
      ss<-0
      for(c in 1:thisNumCohorts){
        ss<-ss+abs(weight_array[c,1]-estWeights[c])
      }
      return(ss)
    }
    fittedInitialG<-optim(par=c(0.1, 1),fn=ll_wIC_fn)$par; fittedG<-fittedInitialG[1]; fitted_n<-fittedInitialG[2]
    fittedInitialWeights<-weightsICFromG_fn(fittedG, fitted_n)
    plot(weight_array[,1])
    points(fittedInitialWeights,type="l")
    dev.off()
    
    GbyAgeClass<-unlist(lapply(seq(1,thisNumCohorts), FUN=function(x){fittedG/(x*fitted_n)}))
    
    storeAllMums[g,1:thisNumCohorts]<-exp(GbyAgeClass/365)
  }
    
}

##update biol.prm
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[g],side="both")
  cat(as.character(thisCode),"--")
  thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    #do mum first
    thisVar<-paste("mum_", thisCode,sep=""); x<-grep(thisVar,biolLines)
    curLine<-get_first_number(biolLines[x+1],n="all")
    temp<-signif(storeAllMums[g,],5); temp[is.na(temp)]<-1; temp<-temp[1:thisNumCohorts]
    newLine<-paste(temp,collapse="\t")
    biolLines[x+1]<-newLine
    ## do mL - don't overwrite these, just enter as an extra line with max value
    thisVar<-paste(thisCode, "_mL",sep=""); x<-grep(thisVar,biolLines)
    curLine<-get_first_number(biolLines[x+1],n="all")
    temp<-storeAllMs[g]; maxMl<-signif(1-exp(-temp/365),4)
    newLine<-paste(paste(curLine, collapse=" "), "\n #Max value ",paste(rep(maxMl, 2), collapse=" "))
    biolLines[x+1]<-newLine
  }
}

writeLines(biolLines,paste(this_path,"..\\CRAM_BH_hybrid_biol.prm",sep=""))

#write Mum and Ms out
write.csv(storeAllMums, paste(DIR$'Tables', "EstimatedMumsFromICs.csv", sep=""), row.names = FALSE)
write.csv(storeAllMs, paste(DIR$'Tables', "EstimatedMsFromICs.csv", sep=""), row.names = FALSE)
