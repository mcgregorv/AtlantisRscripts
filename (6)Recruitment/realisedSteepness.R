#############################################
## NOT RIGHT - USE YOY VERSION INSTEAD
#############################################

#read in the tracers (numbers), calculate spawning stock abundance, calculate rectruited abundance, then calculate steepness
this_run<-"base"
# this_out<-paste("Base2",sep="")
this_out<-paste("Short2",sep="")
# 
basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
outPath<-paste(basePath,this_run,"\\","output",this_out,"\\",sep="")

plotPath<-paste(basePath,"\\Figures\\Recruitment\\",sep="")

biolLines<-readLines(paste(basePath,"CRAM_base_biol.prm",sep=""))
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
calcR<-function(SSB,a,b){
  R<-(a*SSB)/(b+SSB)
  return(R)
}

#set up array to store SS and R for each group for each ts
storeSS<-array(0,dim=c(ng,nlayers,nboxes,nts)); storeR<-storeSS

for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[g]); thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    blPar<-paste("FSPB_",thisCode,sep="")
    x<-biolLines[grep(blPar,biolLines)+1]; thisPropSpawning<-get_first_number(x,n="all")
    blPar<-paste(thisCode,"_AgeClassSize",sep="")
    x<-biolLines[grep(blPar,biolLines)]; thisAgeClassSize<-get_first_number(x)
    for(c in 1:thisNumCohorts){
      thisTracer<-paste(thisName,c,"_Nums",sep="")
      thisData<-ncvar_get(ThisNC.nc,thisTracer)
      thisProp<-thisPropSpawning[c]
      storeSS[g,,,]<-storeSS[g,,,]+thisData*thisProp
      if(c==1){
        storeR[g,,,]<-thisData/thisAgeClassSize
      }
    }
  }
}

##alt for reruitment
altStoreR<-0*storeR
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[g]); thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    blPar<-paste("FSPB_",thisCode,sep="")
    x<-biolLines[grep(blPar,biolLines)+1]; thisPropSpawning<-get_first_number(x,n="all")
    blPar<-paste(thisCode,"_AgeClassSize",sep="")
    x<-biolLines[grep(blPar,biolLines)]; thisAgeClassSize<-get_first_number(x)
    c=1
    thisTracer<-paste(thisName,c,"_Nums",sep="")
    thisData<-ncvar_get(ThisNC.nc,thisTracer)
    if(thisAgeClassSize==1){
      thisRecruits<-thisData
    }else{
      thisRecruits<-0*thisData
      # thisRecruits[,,dim(thisData)[3]]<-NA
      for(t in 2:dim(thisData)[3]){
        thisRecruits[,,(t-1)]<-thisData[,,t]*(thisAgeClassSize/(thisAgeClassSize-1))-thisData[,,(t-1)]
      }
    }
    altStoreR[g,,,]<-thisRecruits
  }
}

storeSSsum<-apply(storeSS,c(1,4),sum); storeRsum<-apply(storeR,c(1,4),sum);
storeRsum<-apply(altStoreR,c(1,4),sum); storeRsumOriginal<-apply(storeR,c(1,4),sum)

thisCex<-1.5
for(g in 2:ng){
  thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[g]); thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    thisRM<-get_recruitmentModel(thisCode)
    if(thisRM==19){
      jpeg(paste(plotPath,"BHCompareRealised_",thisCode,".jpg",sep=""),quality=300,height=300,width=340)
      testHokSS<-storeSSsum[g,]; testHokR<-storeRsum[g,]
      this_a<-getBHab(thisCode)$'a'; this_b<-getBHab(thisCode)$'b'
      BHcurveSS<-seq(0,max(testHokSS)*1.2,length.out=1e+4); BHcurveR<-unlist(lapply(BHcurveSS,calcR,a=this_a,b=this_b))
      thisYMax<-max(c(this_a,max(testHokR))); thisXMax<-max(c(max(BHcurveSS),max(testHokSS)))
      plot(BHcurveSS,BHcurveR,type="l",lwd=2,col=myBlue,ylim=c(0,thisYMax),xlim=c(0,thisXMax),ylab="Recruitment (numbers)",xlab="Spawning Stock Abundance (numbers)",cex.axis=thisCex,cex.lab=thisCex,cex=thisCex)
      points(testHokSS,testHokR,pch=8,col=myOrange,cex=thisCex)
      mtext(thisCode,side=3,cex=thisCex,adj=0)
      abline(h=this_a,col=myGrey,lwd=2,lty=2)
      this_calcR<-unlist(lapply(testHokSS,calcR,a=this_a,b=this_b))
      # points(testHokSS,this_calcR,col="red",pch=20,cex=thisCex)
      abline(v=this_b,col=myGrey,lwd=2,lty=2)
      dev.off()
    }
  }
}

###compare alt R with orignial for Hoki
g=28; thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[g]); thisNumCohorts<-groupsDF$NumCohorts[g]
blPar<-paste(thisCode,"_AgeClassSize",sep="")
x<-biolLines[grep(blPar,biolLines)]; thisAgeClassSize<-get_first_number(x)
testAlt<-storeRsum[g,]; testOrig<-storeRsumOriginal[g,]
###

##Linear
get_KDENR<-function(Code){
  #this is the linear recruitment coefficient
  thisPar<-paste("KDENR_",Code,sep="")
  thisOut<-get_first_number(biolLines[grep(thisPar,biolLines)+1])
  return(thisOut)
}

for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[g]); thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    thisRM<-get_recruitmentModel(thisCode)
    if(thisRM==12){
      this_KDENR<-get_KDENR(thisCode)
      jpeg(paste(plotPath,"LinearCompareRealised_",thisCode,".jpg",sep=""),quality=300,height=300)
      testHokSS<-storeSSsum[g,]; testHokR<-storeRsum[g,]
      linearSS<-seq(min(testHokSS),max(testHokSS),length.out=1000); linearR<-this_KDENR*linearSS
      thisYMax<-max(c(max(linearR),max(testHokR))); thisXMax<-max(c(max(BHcurveSS),max(testHokSS)))
      plot(linearSS,linearR,type="l",lwd=2,col=myBlue,ylim=c(0,thisYMax),ylab="Recruitment (numbers)",xlab="Spawning Stock Abundance (numbers)",cex.axis=thisCex,cex.lab=thisCex,cex=thisCex)
      points(testHokSS,testHokR,pch=20,col="red",cex=thisCex)
      mtext(thisCode,side=3,cex=thisCex,adj=0)
      dev.off()
    }
  }
}
#####################################################
##read in second model and add SSR values to plots
outPath2<-paste(basePath,this_run,"\\","outputBase2\\",sep="")

ThisNC2.nc<-nc_open(paste(outPath2,"output.nc",sep=""))
thisVol2<-ncvar_get(ThisNC2.nc,"volume")
nts2<-dim(thisVol2)[3]; nboxes<-dim(thisVol)[2]; nlayers<-dim(thisVol)[1]

#set up array to store SS and R for each group for each ts
storeSS2<-array(0,dim=c(ng,nlayers,nboxes,nts2)); storeR2<-storeSS2

for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[g]); thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    blPar<-paste("FSPB_",thisCode,sep="")
    x<-biolLines[grep(blPar,biolLines)+1]; thisPropSpawning<-get_first_number(x,n="all")
    blPar<-paste(thisCode,"_AgeClassSize",sep="")
    x<-biolLines[grep(blPar,biolLines)]; thisAgeClassSize<-get_first_number(x)
    for(c in 1:thisNumCohorts){
      thisTracer<-paste(thisName,c,"_Nums",sep="")
      thisData<-ncvar_get(ThisNC2.nc,thisTracer)
      thisProp<-thisPropSpawning[c]
      storeSS2[g,,,]<-storeSS2[g,,,]+thisData*thisProp
      if(c==1){
        storeR2[g,,,]<-thisData/thisAgeClassSize
      }
    }
  }
}

storeSSsum2<-apply(storeSS2,c(1,4),sum); storeRsum2<-apply(storeR2,c(1,4),sum)

combinedSSsum<-cbind(storeSSsum,storeSSsum2); combinedRsum<-cbind(storeRsum,storeRsum2)

thisCex<-1.5
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[g]); thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    thisRM<-get_recruitmentModel(thisCode)
    if(thisRM==19){
      jpeg(paste(plotPath,"BHCompareRealised_",thisCode,".jpg",sep=""),quality=300,height=300,width=340)
      testHokSS<-storeSSsum[g,]; testHokR<-storeRsum[g,]
      testHokSS2<-storeSSsum2[g,]; testHokR2<-storeRsum2[g,]
      this_a<-getBHab(thisCode)$'a'; this_b<-getBHab(thisCode)$'b'
      BHcurveSS<-seq(0,max(c(max(testHokSS),max(testHokSS)))*1.2,length.out=1e+4); BHcurveR<-unlist(lapply(BHcurveSS,calcR,a=this_a,b=this_b))
      thisYMax<-max(c(this_a,max(testHokR))); thisXMax<-max(c(max(BHcurveSS),max(testHokSS)))
      plot(BHcurveSS,BHcurveR,type="l",lwd=2,col=myBlue,ylim=c(0,thisYMax),xlim=c(0,thisXMax),ylab="Recruitment (numbers)",xlab="Spawning Stock Abundance (numbers)",cex.axis=thisCex,cex.lab=thisCex,cex=thisCex)
      points(testHokSS,testHokR,pch=8,col=myOrange,cex=thisCex)
      # points(testHokSS2,testHokR2,pch=8,col=myPurple,cex=thisCex)
      mtext(thisCode,side=3,cex=thisCex,adj=0)
      abline(h=this_a,col=myGrey,lwd=2,lty=2)
      this_calcR<-unlist(lapply(testHokSS,calcR,a=this_a,b=this_b))
      # points(testHokSS,this_calcR,col="red",pch=20,cex=thisCex)
      abline(v=this_b,col=myGrey,lwd=2,lty=2)
      dev.off()
    }
  }
}
#####################################################

###########################################
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[g]); thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    thisRM<-get_recruitmentModel(thisCode)
    if(thisRM==12){
      this_KDENR<-get_KDENR(thisCode)
      jpeg(paste(plotPath,"BHTestingLinear_",thisCode,".jpg",sep=""),quality=300,height=300)
      testHokSS<-storeSSsum[g,]; testHokR<-storeRsum[g,]
      linearSS<-seq(min(testHokSS),max(testHokSS),length.out=1000); linearR<-this_KDENR*linearSS
      this_calcR<-this_KDENR*testHokSS
      ratio<-testHokR/this_calcR
      plot(ratio,type="l",lwd=2,col=myBlue,ylim=c(0,max(ratio)),ylab="Ratio actual/calculated recruitment",xlab="Timestep (years)",cex.axis=thisCex,cex.lab=thisCex,cex=thisCex)
      mtext(thisCode,side=3,cex=thisCex,adj=0)
      mtext(this_KDENR,side=3,adj=1,cex=thisCex)
      dev.off()
    }
  }
}

###############################################################################################
# for(g in 1:ng){
#   thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[g]); thisNumCohorts<-groupsDF$NumCohorts[g]
#   if(thisNumCohorts>1){
#     thisRM<-get_recruitmentModel(thisCode)
#     if(thisRM==19){
#       jpeg(paste(plotPath,"BHCompareRealised_",thisCode,".jpg",sep=""),quality=300,height=300)
#       testHokSS<-storeSSsum[g,]; testHokR<-storeRsum[g,]
#       this_a<-getBHab(thisCode)$'a'; this_b<-getBHab(thisCode)$'b'
#       thisYMax<-max(c(this_a,max(testHokR)))
#       plot(testHokSS,testHokR,pch=20,col=myBlue,ylim=c(0,thisYMax),ylab="Recruitment (numbers)",xlab="Spawning Stock Abundance (numbers)",cex.axis=thisCex,cex.lab=thisCex,cex=thisCex)
#       mtext(thisCode,side=3,cex=thisCex,adj=0)
#       abline(h=this_a,col=myGrey,lwd=2,lty=2)
#       this_calcR<-unlist(lapply(testHokSS,calcR,a=this_a,b=this_b))
#       points(testHokSS,this_calcR,col="red",pch=20,cex=thisCex)
#       abline(v=this_b,col=myGrey,lwd=2,lty=2)
#       dev.off()
#     }
#   }
# }

######################################################################################################
#not sumarised over cells
# for(g in 1:ng){
#   thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[g]); thisNumCohorts<-groupsDF$NumCohorts[g]
#   if(thisNumCohorts>1){
#     thisRM<-get_recruitmentModel(thisCode)
#     if(thisRM==19){
#       jpeg(paste(plotPath,"BHCompareRealisedByCell_",thisCode,".jpg",sep=""),quality=300,height=300)
#       testHokSS<-as.vector(storeSS[g,,,]); testHokR<-as.vector(storeR[g,,,])
#       this_a<-getBHab(thisCode)$'a'; this_b<-getBHab(thisCode)$'b'
#       this_calcR<-unlist(lapply(testHokSS,calcR,a=this_a,b=this_b))
#       thisYMax<-max(c(max(this_calcR),max(testHokSS)))
#       plot(testHokSS,testHokR,pch=20,col=myBlue,ylim=c(0,thisYMax),ylab="Recruitment (numbers)",xlab="Spawning Stock Abundance (numbers)",cex.axis=thisCex,cex.lab=thisCex,cex=thisCex)
#       mtext(thisCode,side=3,cex=thisCex,adj=0)
#       abline(h=this_a,col=myGrey,lwd=2,lty=2)
#       points(testHokSS,this_calcR,col="red",pch=20,cex=thisCex)
#       abline(v=this_b,col=myGrey,lwd=2,lty=2)
#       dev.off()
#     }
#   }
# }
# 
# 
# testHokSS<-storeSSsum[28,]; testHokR<-storeRsum[28,]
# this_a<-getBHab("HOK")$'a'; this_b<-getBHab("HOK")$'b'
# thisYMax<-max(c(this_a,max(testHokSS)))
# plot(testHokSS,testHokR,pch=20,col=myBlue_trans,ylim=c(0,thisYMax))
# abline(h=this_a,col=myGrey,lwd=2,lty=2)
# this_calcR<-unlist(lapply(testHokSS,calcR,a=this_a,b=this_b))
# points(testHokSS,this_calcR,col="red",pch=20)
# abline(v=this_b,col=myGrey,lwd=2,lty=2)
# 
# 
# 
# 
# testHokSS<-as.vector(storeSS[28,,,]); testHokR<-as.vector(storeR[28,,,])
# this_a<-getBHab("HOK")$'a'; this_b<-getBHab("HOK")$'b'
# this_calcR<-unlist(lapply(testHokSS,calcR,a=this_a,b=this_b))
# thisYMax<-max(c(max(testHokSS),max(this_calcR)))
# plot(testHokSS,testHokR,pch=20,col=myBlue_trans,ylim=c(0,thisYMax))
# abline(h=this_a,col=myGrey,lwd=2,lty=2)
# points(testHokSS,this_calcR,col="red",pch=20)
# abline(v=this_b,col=myGrey,lwd=2,lty=2)
# 
# IC.nc<-nc_open(paste(basePath,"CRAM_input_short.nc",sep=""))
# ICvars<-unique(sort(names(IC.nc$var)))
# thisIC<-array(0,dim=dim(thisVol)[1:2])
# thisMature<-3
# for(c in 1:thisNumCohorts){
#   if(c>=thisMature){
#     thisVar<-paste(thisName,c,"_Nums",sep="")
#     temp<-ncvar_get(IC.nc,thisVar)[,,1]
#     thisIC<-thisIC+temp
#   }
# }
# abline(v=sum(thisIC),col=myOrange)
