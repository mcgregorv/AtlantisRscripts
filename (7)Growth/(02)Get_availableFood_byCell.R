#read in initial conditions and availability matrix then calculate how much food is available (in mg N) for each predator by cell (box/layer)

basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

plotPath<-paste(basePath,"Figures\\growth\\",sep="") #this is created in DIET_read_andCreateCsv.R

#created in DIET_read_andCreateCsv.R in setting_up/Testing/ RERUN this if any changes to pPREY 
pPREY_df<-read.csv(paste(basePath,"pPREY_OUT.csv",sep="")) 

biolFile<-paste(basePath,"CRAM_base_biol.prm",sep=""); biolLines<-readLines(biolFile)

initFile<-paste(basePath,"CRAM_input_short.nc",sep=""); initiNC<-nc_open(initFile)

groupsFile<-paste(basePath,"CRAM_groups.csv",sep=""); groupsDF<-read.csv(groupsFile); ng<-dim(groupsDF)[1]

max_nCohorts<-10; nboxes<-30; nlayers<-6

volume<-ncvar_get(initiNC,"volume")[,,1]

mumByGroupCohort<-read.csv(paste(basePath,"mumByCohort.csv",sep=""),header=TRUE)


allStructN<-NULL; allResN<-NULL; allDens<-NULL; blankArray<-array(NA,dim=c(nlayers,nboxes,max_nCohorts))

groupsDF$CohortMature<-NA

for(g in 1:ng){
  thisCode<-as.character(groupsDF$Code[g]); thisName<-str_trim(groupsDF$Name[g],side="both"); thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    #while looping here, also add cohort mature in here
    temp<-biolLines[grep(paste(thisCode,"_age_mat",sep=""),biolLines)]
    thisMature<-get_first_number(temp)
    groupsDF$CohortMature[g]<-thisMature
    SNarray<-blankArray; RNarray<-blankArray; DENarray<-blankArray
    for(c in 1:thisNumCohorts){
      thisVar<-paste(thisName,c,"_StructN",sep="")
      thisICdata<-ncvar_get(initiNC,thisVar)
      SNarray[,,c]<-thisICdata[,,1]
      #ResN
      thisVar<-paste(thisName,c,"_ResN",sep="")
      thisICdata<-ncvar_get(initiNC,thisVar)
      RNarray[,,c]<-thisICdata[,,1]
      #Density
      thisVar<-paste(thisName,c,"_Nums",sep="")
      thisICdata<-ncvar_get(initiNC,thisVar)
      DENarray[,,c]<-thisICdata[,,1]
    }
    allStructN[[as.character(thisCode)]]<-SNarray; allResN[[as.character(thisCode)]]<-RNarray; allDens[[as.character(thisCode)]]<-DENarray
  }
}

#get total mg N of each species/cohort group per cell
allBiomass<-NULL
for(g in 1:ng){
  thisCode<-as.character(groupsDF$Code[g]); thisName<-str_trim(groupsDF$Name[g],side="both"); thisNumCohorts<-groupsDF$NumCohorts[g]
  # cat(thisCode,", ")
  if(thisNumCohorts>1){
    #density times weight
    thisBiomass<-allDens[[thisCode]]*(allStructN[[thisCode]]+allResN[[thisCode]])
  }else{
    #read in _N and multiply by volume
    tempData<-ncvar_get(initiNC,paste(thisName,"_N",sep=""))
    if(length(dim(tempData))==2){
      tempData<-tempData[,1]
      thisBiomass<-tempData*volume[6,]
    }else{
      tempData<-tempData[,,1]
      thisBiomass<-tempData*volume
    }
  }
  allBiomass[[thisCode]]<-thisBiomass
}

zeroDF<-data.frame(matrix(0,nrow=nlayers,ncol=nboxes))
#calculate how much food there is in each box for each predator to eat (not taking into account clearance, handling time or max growth)
howMuchFoodAdults<-NULL; howMuchFoodJuv<-NULL
for(g in 1:ng){
  thisCode<-as.character(groupsDF$Code[g]); thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    #need to do seperate diet for adults and juveniles
    thispPrey<-pPREY_df[pPREY_df$PredatorCode==thisCode,]
    adultPrey<-zeroDF; juvPrey<-zeroDF
    for(p in 1:ng){
      thisPreyCode<-as.character(groupsDF$Code[p]); thisPreyNumCohorts<-groupsDF$NumCohorts[p]; 
      thisBiomass<-allBiomass[[thisPreyCode]]
      if(thisPreyNumCohorts>1){
        thisPreyAvail22<-thispPrey[thispPrey$PreyAge==2 & thispPrey$PredatorAge==2,c(thisPreyCode)]
        thisPreyAvail21<-thispPrey[thispPrey$PreyAge==2 & thispPrey$PredatorAge==1,c(thisPreyCode)]
        thisPreyAvail12<-thispPrey[(thispPrey$PreyAge==1 | is.na(thispPrey$PreyAge)) & thispPrey$PredatorAge==2,c(thisPreyCode)]
        thisPreyAvail11<-thispPrey[(thispPrey$PreyAge==1 | is.na(thispPrey$PreyAge)) & thispPrey$PredatorAge==1,c(thisPreyCode)]
        thisCM<-groupsDF$CohortMature[p]
        temp<-thisBiomass[,,(thisCM+1):thisPreyNumCohorts]
        if(length(dim(temp))==3){
          thisAdBiomass<-apply(temp,c(1,2),sum)
        }else{
          thisAdBiomass<-temp
        }
        temp<-thisBiomass[,,1:thisCM]
        if(length(dim(temp))==3){
          thisJuvBiomass<-apply(temp,c(1,2),sum)
        }else{
          thisJuvBiomass<-temp
        }
        adultPrey<-adultPrey+thisAdBiomass*thisPreyAvail22+thisJuvBiomass*thisPreyAvail12
        juvPrey<-juvPrey+thisAdBiomass*thisPreyAvail21+thisJuvBiomass*thisPreyAvail11
      } else {
        thisPreyAvail12<-thispPrey[(thispPrey$PreyAge==1 | is.na(thispPrey$PreyAge)) & thispPrey$PredatorAge==2,c(thisPreyCode)]
        thisPreyAvail11<-thispPrey[(thispPrey$PreyAge==1 | is.na(thispPrey$PreyAge)) & thispPrey$PredatorAge==1,c(thisPreyCode)]
        if(length(dim(thisBiomass))==2){
          adultPrey<-adultPrey+thisBiomass*thisPreyAvail12
          juvPrey<-juvPrey+thisBiomass*thisPreyAvail11
        }
      }
    }
    howMuchFoodAdults[[as.character(thisCode)]]<-adultPrey; howMuchFoodJuv[[as.character(thisCode)]]<-juvPrey
  }
}

#some plots of prey avaibility and prey required ##ADD in this loop to store min proportion of required prey avail
for(g in 14:ng){
  thisCode<-as.character(groupsDF$Code[g]); thisName<-str_trim(groupsDF$Name[g],side="both"); thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    thisNumbers<-allDens[[thisCode]]
    thisMum<-mumByGroupCohort[mumByGroupCohort$Code==thisCode,1:thisNumCohorts]
    thisDailyRequiredFood<-apply(thisNumbers[,,1:thisNumCohorts],c(1,2),FUN=function(x){as.double(x*thisMum)})
    thisCM<-groupsDF$CohortMature[g]
    tempAd<-thisDailyRequiredFood[(thisCM+1):thisNumCohorts,,]; 
    tempJuv<-thisDailyRequiredFood[1:thisCM,,]
    if(length(dim(tempAd))==3){
      adDailyFoodRequired<-apply(tempAd,c(2,3),sum)
    } else{
      adDailyFoodRequired<-tempAd
    }
    if(length(dim(tempJuv))==3){
      juvDailyFoodRequired<-apply(tempJuv,c(2,3),sum)
    } else{
      juvDailyFoodRequired<-tempJuv
    }
    
    adAxis<-pretty(seq(0,max(adDailyFoodRequired),length.out=5))
    juvAxis<-pretty(seq(0,max(juvDailyFoodRequired),length.out=5))
    
    foodRequiredColor<-myBlue; foodAvailableColor<-myGreen
    
    pdf(paste(plotPath,"\\",thisCode,"AdultFoodCompare.pdf",sep=""))
    par(lend=1,mar=c(0,4,1,4),mfrow=c(nlayers,1),oma=c(5,3,1,1),las=1)
    for(l in 1:nlayers){
      plot(adDailyFoodRequired[l,],type="h",lwd=5,col=foodRequiredColor,yaxt="n",xaxt="n",ylab="",xlab="")
      axis(at=adAxis,labels=adAxis,side=4)
      par(new=TRUE)
      plot(as.double(howMuchFoodAdults[[thisCode]][l,]),type="l",col=foodAvailableColor,lwd=2,xlab="",ylab="",xaxt="n")
      
    }
    axis(at=seq(1,nboxes),labels=seq(1,nboxes),side=1,outer=TRUE)
    mtext("Daily food required",side=4,adj=0.5,col=foodRequiredColor,line=-1.3,outer=TRUE,las=0)
    mtext("Food available (mg N)",side=2,adj=0.5,line=1,outer=TRUE,col=foodAvailableColor,las=0)
    mtext("Polygon",side=1,adj=0.5,line=2.5,outer=TRUE)
    dev.off()
    
    # if(sum(juvDailyFoodRequired,na.rm=TRUE)>0)
    pdf(paste(plotPath,"\\",thisCode,"JuvFoodCompare.pdf",sep=""))
    par(lend=1,mar=c(0,4,1,4),mfrow=c(nlayers,1),oma=c(5,3,1,1),las=1)
    for(l in 1:nlayers){
      plot(juvDailyFoodRequired[l,],type="h",lwd=5,col=foodRequiredColor,yaxt="n",xaxt="n",ylab="",xlab="")
      axis(at=juvAxis,labels=juvAxis,side=4)
      par(new=TRUE)
      plot(as.double(howMuchFoodJuv[[thisCode]][l,]),type="l",col=foodAvailableColor,lwd=2,xlab="",ylab="",xaxt="n")
      
      
    }
    axis(at=seq(1,nboxes),labels=seq(1,nboxes),side=1,outer=TRUE)
    mtext("Daily food required",side=4,adj=0.5,col=foodRequiredColor,line=-1.3,outer=TRUE,las=0)
    mtext("Food available (mg N)",side=2,adj=0.5,line=1,outer=TRUE,col=foodAvailableColor,las=0)
    mtext("Polygon",side=1,adj=0.5,line=2.5,outer=TRUE)
    dev.off()
  }
}

#calculate how many of each prey are eaten if each predator eats what is available to it, capped at its (max growth (mum) divided by E) multiplied by predator numbers. Do age-structured first as they are in numbers
howMuchDeathAdults<-NULL; howMuchDeathJuv<-NULL

#for a predator, get it's numbers
#get how much food is available to it for each cell. If it is less than required, it will eat it all and flag it
#if it is more than required, it will eat a proportion of each available prey prop=preyEatn/preyAvail
#store the realised diet for the predator, and remove the fatalities from the prey
#the numbers of prey are calculated from weight eaten divided by average weight of individuals
getFoodProportion<-function(foodReq,foodAvail){
  y<-NA
  test1<-foodReq>0; test2<-foodAvail>0
  if(test1==TRUE & (test2==FALSE | is.na(test2))){y<-"missing"}
  if(test1==TRUE & test2==TRUE){
    y<-foodReq/foodAvail
    
  }
  return(y)
}

numbersEaten11<-array(dim=c(ng,ng,nlayers,nboxes))
numbersEaten12<-numbersEaten11; numbersEaten21<-numbersEaten11; numbersEaten22<-numbersEaten11

for(g in 1:ng){
  thisCode<-as.character(groupsDF$Code[g]); thisName<-str_trim(groupsDF$Name[g],side="both"); thisNumCohorts<-groupsDF$NumCohorts[g]
  cat(thisCode,"--")
  if(thisNumCohorts>1){
    thisNumbers<-allDens[[thisCode]]
    #add together adult and juvenile numbers
    thisMC<-groupsDF$CohortMature[g]
    juvNumbersByCohort<-thisNumbers[,,1:thisMC]; adNumbersByCohort<-thisNumbers[,,(thisMC+1):thisNumCohorts]
    if(length(dim(juvNumbersByCohort))==3){
      juvNumbers<-apply(juvNumbersByCohort,c(1,2),sum)
    } else{
      juvNumbers<-juvNumbersByCohort
    }
    if(length(dim(adNumbersByCohort))==3){
      adNumbers<-apply(adNumbersByCohort,c(1,2),sum)
    } else{
      adNumbers<-adNumbersByCohort
    }
    #use average mum for adults and juveniles
    thisMum<-as.double(mumByGroupCohort[mumByGroupCohort$Code==as.character(thisCode),1:(dim(mumByGroupCohort)[2]-1)])
    thisAdMum<-thisMum[(thisCM):thisNumCohorts]; thisJuvMum<-thisMum[1:(thisCM-1)]
    if(length(thisAdMum)>1){thisAdMum<-mean(thisAdMum,na.rm=TRUE)}
    if(length(thisJuvMum)>1){thisJuvMum<-mean(thisJuvMum,na.rm=TRUE)}
    adTotalPrey<-howMuchFoodAdults[[thisCode]]
    juvTotalPrey<-howMuchFoodJuv[[thisCode]]
    adFoodRequired<-adNumbers*(thisAdMum/thisE); juvFoodRequired<-juvNumbers*(thisJuvMum/thisE)

    adPropToEat<-0*adFoodRequired;   juvPropToEat<-0*juvFoodRequired
    for(i in 1:nrow(adPropToEat)){
      for(j in 1:ncol(adPropToEat)){
        adPropToEat[i,j]<-getFoodProportion(foodReq=adFoodRequired[i,j],foodAvail=adTotalPrey[i,j])
        juvPropToEat[i,j]<-getFoodProportion(foodReq=juvFoodRequired[i,j],foodAvail=juvTotalPrey[i,j])
      }
    }
    if("missing" %in% adPropToEat){cat("Food missing for ",thisCode," adults\n")}
    if("missing" %in% juvPropToEat){cat("Food missing for ",thisCode," juveniles\n")}
    thispPrey<-pPREY_df[pPREY_df$PredatorCode==thisCode,]
    
    #now apply this proportion to each prey and turn it into numbers and store that
    for(p in 1:ng){
      thisPreyCode<-as.character(groupsDF$Code[p]); thisPreyNumCohorts<-groupsDF$NumCohorts[p]; 
      thisBiomass<-allBiomass[[thisPreyCode]]
      if(thisPreyNumCohorts>1){
        thisPreyAvail22<-thispPrey[thispPrey$PreyAge==2 & thispPrey$PredatorAge==2,c(thisPreyCode)]
        thisPreyAvail21<-thispPrey[thispPrey$PreyAge==2 & thispPrey$PredatorAge==1,c(thisPreyCode)]
        thisPreyAvail12<-thispPrey[(thispPrey$PreyAge==1 | is.na(thispPrey$PreyAge)) & thispPrey$PredatorAge==2,c(thisPreyCode)]
        thisPreyAvail11<-thispPrey[(thispPrey$PreyAge==1 | is.na(thispPrey$PreyAge)) & thispPrey$PredatorAge==1,c(thisPreyCode)]
        if(sum(thisPreyAvail11,thisPreyAvail12,thisPreyAvail21,thisPreyAvail22)>0){
          #get numbers by juv/ad and mean weight per idividual
          thisPreyMC<-groupsDF$CohortMature[p]
          thisPreyNumbers<-allDens[[thisPreyCode]]
          #add together adult and juvenile numbers
          juvPreyNumbersByCohort<-thisPreyNumbers[,,1:thisPreyMC]; adPreyNumbersByCohort<-thisPreyNumbers[,,(thisPreyMC+1):thisPreyNumCohorts]
          if(length(dim(juvPreyNumbersByCohort))==3){
            juvPreyNumbers<-apply(juvPreyNumbersByCohort,c(1,2),sum)
          } else{
            juvPreyNumbers<-juvPreyNumbersByCohort
          }
          if(length(dim(adPreyNumbersByCohort))==3){
            adPreyNumbers<-apply(adPreyNumbersByCohort,c(1,2),sum)
          } else{
            adPreyNumbers<-adPreyNumbersByCohort
          }
          ## ResN
          thisPreyNumbers<-allResN[[thisPreyCode]]
          #add together adult and juvenile numbers
          juvPreyNumbersByCohort<-thisPreyNumbers[,,1:thisPreyMC]; adPreyNumbersByCohort<-thisPreyNumbers[,,(thisPreyMC+1):thisPreyNumCohorts]
          if(length(dim(juvPreyNumbersByCohort))==3){
            juvPreyResN<-apply(juvPreyNumbersByCohort,c(1,2),sum)
          } else{
            juvPreyResN<-juvPreyNumbersByCohort
          }
          if(length(dim(adPreyNumbersByCohort))==3){
            adPreyResN<-apply(adPreyNumbersByCohort,c(1,2),sum)
          } else{
            adPreyResN<-adPreyNumbersByCohort
          }
          ## StructN
          thisPreyNumbers<-allStructN[[thisPreyCode]]
          #add together adult and juvenile numbers
          juvPreyNumbersByCohort<-thisPreyNumbers[,,1:thisPreyMC]; adPreyNumbersByCohort<-thisPreyNumbers[,,(thisPreyMC+1):thisPreyNumCohorts]
          if(length(dim(juvPreyNumbersByCohort))==3){
            juvPreyStructN<-apply(juvPreyNumbersByCohort,c(1,2),sum)
          } else{
            juvPreyStructN<-juvPreyNumbersByCohort
          }
          if(length(dim(adPreyNumbersByCohort))==3){
            adPreyStructN<-apply(adPreyNumbersByCohort,c(1,2),sum)
          } else{
            adPreyStructN<-adPreyNumbersByCohort
          }
          ##biomass by adult and juvenile
          juvBiomass<-apply(thisBiomass[,,1:thisPreyMC],c(1,2),sum,na.rm=TRUE); adBiomass<-apply(thisBiomass[,,(thisPreyMC+1):thisPreyNumCohorts],c(1,2),sum,na.rm=TRUE)
          #total weight
          adPreyWeightPerInd<-adPreyStructN+adPreyResN; juvPreyWeightPerInd<-juvPreyStructN+juvPreyResN
          ##
          #amount (mg N) eaten 
          amountEaten<-thisPreyAvail22*adBiomass*as.double(adPropToEat)
          #number of individuals required (round up)
          numberEaten<-amountEaten / adPreyWeightPerInd
          numbersEaten22[g,p,,]<-numberEaten
          ##
          ## adults eating juveniles (prey then predator)
          amountEaten<-thisPreyAvail12*juvBiomass*as.double(adPropToEat)
          #number of individuals required (round up)
          numberEaten<-amountEaten / juvPreyWeightPerInd
          numbersEaten12[g,p,,]<-numberEaten
          ##
          ## juveniles eating juveniles (prey then predator)
          amountEaten<-thisPreyAvail11*juvBiomass*as.double(juvPropToEat)
          #number of individuals required (round up)
          numberEaten<-amountEaten / juvPreyWeightPerInd
          numbersEaten11[g,p,,]<-numberEaten
          ##
          ## juveniles eating adults (prey then predator)
          amountEaten<-thisPreyAvail21*adBiomass*as.double(juvPropToEat)
          #number of individuals required (round up)
          numberEaten<-amountEaten / adPreyWeightPerInd
          numbersEaten21[g,p,,]<-numberEaten
        }
      }
    }
  }
}
 
#now plot the number eaten compared with the number that exist
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    juvEatenByPred<-numbersEaten11[,g,,]+numbersEaten12[,g,,]
    adEatenByPred<-numbersEaten21[,g,,]+numbersEaten22[,g,,]
    juvEaten<-apply(juvEatenByPred,c(2,3),sum,na.rm=TRUE)
    adEaten<-apply(adEatenByPred,c(2,3),sum,na.rm=TRUE)
    #get adult and juv numbers too
    thisNumbers<-allDens[[as.character(thisCode)]]
    #add together adult and juvenile numbers
    thisMC<-groupsDF$CohortMature[g]
    juvNumbersByCohort<-thisNumbers[,,1:thisMC]; adNumbersByCohort<-thisNumbers[,,(thisMC+1):thisNumCohorts]
    if(length(dim(juvNumbersByCohort))==3){
      juvNumbers<-apply(juvNumbersByCohort,c(1,2),sum)
    } else{
      juvNumbers<-juvNumbersByCohort
    }
    if(length(dim(adNumbersByCohort))==3){
      adNumbers<-apply(adNumbersByCohort,c(1,2),sum)
    } else{
      adNumbers<-adNumbersByCohort
    }

    pdf(paste(plotPath,"\\",thisCode,"JuvDeathByPredatorCompare.pdf",sep=""))
    par(lend=1,mar=c(0,4,1,4),mfrow=c(nlayers,1),oma=c(5,3,1,1),las=1)
    for(l in 1:nlayers){
      plot(juvNumbers[l,],type="h",lwd=5,col=foodRequiredColor,xaxt="n",ylab="",xlab="")
      points(juvEaten[l,]*365,type="l",col=foodAvailableColor,lwd=2,xlab="",ylab="",xaxt="n")
      
      
    }
    axis(at=seq(1,nboxes),labels=seq(1,nboxes),side=1,outer=TRUE)
    mtext("Numbers",side=2,adj=0.5,line=1,outer=TRUE,las=0)
    mtext("Polygon",side=1,adj=0.5,line=2.5,outer=TRUE)
    dev.off()
    
    pdf(paste(plotPath,"\\",thisCode,"adDeathByPredatorCompare.pdf",sep=""))
    par(lend=1,mar=c(0,4,1,4),mfrow=c(nlayers,1),oma=c(5,3,1,1),las=1)
    for(l in 1:nlayers){
      plot(adNumbers[l,],type="h",lwd=5,col=foodRequiredColor,xaxt="n",ylab="",xlab="")
      points(adEaten[l,]*365,type="l",col=foodAvailableColor,lwd=2,xlab="",ylab="",xaxt="n")
    }
    axis(at=seq(1,nboxes),labels=seq(1,nboxes),side=1,outer=TRUE)
    mtext("Numbers",side=2,adj=0.5,line=1,outer=TRUE,las=0)
    mtext("Polygon",side=1,adj=0.5,line=2.5,outer=TRUE)
    dev.off()
    
    cat("Max eaten of ",as.character(thisCode)," per year juv is ",max(juvEaten,na.rm=TRUE)*365, " Max numbers ",max(juvNumbers,na.rm=TRUE),"-- ")
    cat(", Max eaten of ",as.character(thisCode)," per year ad is ",max(adEaten,na.rm=TRUE)*365, " Max numbers ",max(adNumbers,na.rm=TRUE),"\n")
     
  }
  
}
  
#group them - those with >=100%, those with >=50% predation mortality, >= 10%, >=1% (based on max), <1%
codes1<-c(); codes2<-c(); codes3<-c(); codes4<-c(); codes5<-c()
codeNumbers1<-c(); codeNumbers2<-c(); codeNumbers3<-c(); codeNumbers4<-c(); codeNumbers5<-c()
for(g in 1:ng){
  #get this death by predators
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    juvEatenByPred<-numbersEaten11[,g,,]+numbersEaten12[,g,,]
    adEatenByPred<-numbersEaten21[,g,,]+numbersEaten22[,g,,]
    juvEaten<-apply(juvEatenByPred,c(2,3),sum,na.rm=TRUE)
    adEaten<-apply(adEatenByPred,c(2,3),sum,na.rm=TRUE)
    #get adult and juv numbers too
    thisNumbers<-allDens[[as.character(thisCode)]]
    #add together adult and juvenile numbers
    thisMC<-groupsDF$CohortMature[g]
    juvNumbersByCohort<-thisNumbers[,,1:thisMC]; adNumbersByCohort<-thisNumbers[,,(thisMC+1):thisNumCohorts]
    if(length(dim(juvNumbersByCohort))==3){
      juvNumbers<-apply(juvNumbersByCohort,c(1,2),sum)
    } else{
      juvNumbers<-juvNumbersByCohort
    }
    if(length(dim(adNumbersByCohort))==3){
      adNumbers<-apply(adNumbersByCohort,c(1,2),sum)
    } else{
      adNumbers<-adNumbersByCohort
    }
    juvProp<-(juvEaten*365)/juvNumbers
    juvPropMax<-max(juvProp,na.rm=TRUE)
    if(juvPropMax>=1){
      codes1<-c(codes1,paste(thisCode,1,sep=""))
      codeNumbers1<-c(codeNumbers1,signif(juvPropMax,2))
    } else if(juvPropMax>=0.5){
      codes2<-c(codes2,paste(thisCode,1,sep=""))
      codeNumbers2<-c(codeNumbers2,signif(juvPropMax,2))
    } else if(juvPropMax>=0.1){
      codes3<-c(codes3,paste(thisCode,1,sep=""))
      codeNumbers3<-c(codeNumbers3,signif(juvPropMax,2))
    } else if(juvPropMax>=0.01){
      codes4<-c(codes4,paste(thisCode,1,sep=""))
      codeNumbers4<-c(codeNumbers4,signif(juvPropMax,2))
    } else {
      codes5<-c(codes5,paste(thisCode,1,sep=""))
      codeNumbers5<-c(codeNumbers5,signif(juvPropMax,2))
    }
    adProp<-(adEaten*365)/adNumbers
    adPropMax<-max(adProp,na.rm=TRUE)
    if(adPropMax>=1){
      codes1<-c(codes1,paste(thisCode,2,sep=""))
      codeNumbers1<-c(codeNumbers1,signif(adPropMax,2))
    } else if(adPropMax>=0.5){
      codes2<-c(codes2,paste(thisCode,2,sep=""))
      codeNumbers2<-c(codeNumbers2,signif(adPropMax,2))
    } else if(adPropMax>=0.1){
      codes3<-c(codes3,paste(thisCode,2,sep=""))
      codeNumbers3<-c(codeNumbers3,signif(adPropMax,2))
    } else if(adPropMax>=0.01){
      codes4<-c(codes4,paste(thisCode,2,sep=""))
      codeNumbers4<-c(codeNumbers4,signif(adPropMax,2))
    } else {
      codes5<-c(codes5,paste(thisCode,2,sep=""))
      codeNumbers5<-c(codeNumbers5,signif(adPropMax,2))
    }
  }

}

# group then in terms of how much of their required food they have available
# at least 1000 times, 100 times, 10 times, 1 times, less than
gcodes1<-c(); gcodes2<-c(); gcodes3<-c(); gcodes4<-c(); gcodes5<-c()
gcodeNumbers1<-c(); gcodeNumbers2<-c(); gcodeNumbers3<-c(); gcodeNumbers4<-c(); gcodeNumbers5<-c()
for(g in 1:ng){
  #get this death by predators
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    thisAdultFoodAvail<-howMuchFoodAdults[[as.character(thisCode)]]
    thisJuvFoodAvail<-howMuchFoodJuv[[as.character(thisCode)]]
    ##daily food reaquired for max growth
    thisNumbers<-allDens[[as.character(thisCode)]]
    thisMum<-mumByGroupCohort[mumByGroupCohort$Code==thisCode,1:thisNumCohorts]
    thisDailyRequiredFood<-apply(thisNumbers[,,1:thisNumCohorts],c(1,2),FUN=function(x){as.double(x*thisMum)})
    thisCM<-groupsDF$CohortMature[g]
    tempAd<-thisDailyRequiredFood[(thisCM+1):thisNumCohorts,,]; 
    tempJuv<-thisDailyRequiredFood[1:thisCM,,]
    if(length(dim(tempAd))==3){
      adDailyFoodRequired<-apply(tempAd,c(2,3),sum)
    } else{
      adDailyFoodRequired<-tempAd
    }
    if(length(dim(tempJuv))==3){
      juvDailyFoodRequired<-apply(tempJuv,c(2,3),sum)
    } else{
      juvDailyFoodRequired<-tempJuv
    }
    #proportion of food avail over food required (per year)
    propJuv<-thisJuvFoodAvail/(juvDailyFoodRequired*365)
    thisJuvMinProp<-min(propJuv,na.rm=TRUE)
    propAd<-thisAdultFoodAvail/(adDailyFoodRequired*365)
    thisAdMinProp<-min(propAd,na.rm=TRUE)  
    if(thisJuvMinProp>=1000){
      gcodes1<-c(gcodes1,paste(thisCode,1,sep=""))
      gcodeNumbers1<-c(gcodeNumbers1,signif(thisJuvMinProp,2))
    } else if(thisJuvMinProp>=100){
      gcodes2<-c(gcodes2,paste(thisCode,1,sep=""))
      gcodeNumbers2<-c(gcodeNumbers2,signif(thisJuvMinProp,2))
    } else if(thisJuvMinProp>=10){
      gcodes3<-c(gcodes3,paste(thisCode,1,sep=""))
      gcodeNumbers3<-c(gcodeNumbers3,signif(thisJuvMinProp,2))
    } else if(thisJuvMinProp>=1){
      gcodes4<-c(gcodes4, paste(thisCode,1,sep=""))
      gcodeNumbers4<-c(gcodeNumbers4,signif(thisJuvMinProp,2))
    } else{
      gcodes5<-c(gcodes5,paste(thisCode,1,sep=""))
      gcodeNumbers5<-c(gcodeNumbers5,signif(thisJuvMinProp,2))
    }
    ##adults
    if(thisAdMinProp>=1000){
      gcodes1<-c(gcodes1,paste(thisCode,2,sep=""))
      gcodeNumbers1<-c(gcodeNumbers1,signif(thisAdMinProp,2))
    } else if(thisAdMinProp>=100){
      gcodes2<-c(gcodes2,paste(thisCode,2,sep=""))
      gcodeNumbers2<-c(gcodeNumbers2,signif(thisAdMinProp,2))
    } else if(thisAdMinProp>=10){
      gcodes3<-c(gcodes3,paste(thisCode,2,sep=""))
      gcodeNumbers3<-c(gcodeNumbers3,signif(thisAdMinProp,2))
    } else if(thisAdMinProp>=1){
      gcodes4<-c(gcodes4, paste(thisCode,2,sep=""))
      gcodeNumbers4<-c(gcodeNumbers4,signif(thisAdMinProp,2))
    } else{
      gcodes5<-c(gcodes5,paste(thisCode,2,sep=""))
      gcodeNumbers5<-c(gcodeNumbers5,signif(thisAdMinProp,2))
    }
    
  }
}
#now compare with recruitment too..?


##do a focus plot for a given predator to check out what is available to them to eat by box and layer
# and also the biomass of their prey by box and layer
thisCode<-"SSO"
thispPrey<-pPREY_df[pPREY_df$PredatorCode==thisCode,]
#gather the data
adultPrey<-NULL; juvPrey<-NULL
for(p in 1:ng){
  thisPreyCode<-as.character(groupsDF$Code[p]); thisPreyNumCohorts<-groupsDF$NumCohorts[p]; 
  thisBiomass<-allBiomass[[thisPreyCode]]
  if(thisPreyNumCohorts>1){
    thisPreyAvail22<-thispPrey[thispPrey$PreyAge==2 & thispPrey$PredatorAge==2,c(thisPreyCode)]
    thisPreyAvail21<-thispPrey[thispPrey$PreyAge==2 & thispPrey$PredatorAge==1,c(thisPreyCode)]
    thisPreyAvail12<-thispPrey[(thispPrey$PreyAge==1 | is.na(thispPrey$PreyAge)) & thispPrey$PredatorAge==2,c(thisPreyCode)]
    thisPreyAvail11<-thispPrey[(thispPrey$PreyAge==1 | is.na(thispPrey$PreyAge)) & thispPrey$PredatorAge==1,c(thisPreyCode)]
    thisCM<-groupsDF$CohortMature[p]
    temp<-thisBiomass[,,(thisCM+1):thisPreyNumCohorts]
    if(length(dim(temp))==3){
      thisAdBiomass<-apply(temp,c(1,2),sum)
    }else{
      thisAdBiomass<-temp
    }
    temp<-thisBiomass[,,1:thisCM]
    if(length(dim(temp))==3){
      thisJuvBiomass<-apply(temp,c(1,2),sum)
    }else{
      thisJuvBiomass<-temp
    }
    adultPrey[[paste(thisPreyCode,1,sep="")]]<-thisJuvBiomass*thisPreyAvail12
    adultPrey[[paste(thisPreyCode,2,sep="")]]<-thisAdBiomass*thisPreyAvail22
    
    juvPrey[[paste(thisPreyCode,1,sep="")]]<-thisJuvBiomass*thisPreyAvail11
    juvPrey[[paste(thisPreyCode,2,sep="")]]<-thisAdBiomass*thisPreyAvail21
    
  } else {
    thisPreyAvail12<-thispPrey[(thispPrey$PreyAge==1 | is.na(thispPrey$PreyAge)) & thispPrey$PredatorAge==2,c(thisPreyCode)]
    thisPreyAvail11<-thispPrey[(thispPrey$PreyAge==1 | is.na(thispPrey$PreyAge)) & thispPrey$PredatorAge==1,c(thisPreyCode)]
    if(length(dim(thisBiomass))==2){
      adultPrey[[thisPreyCode]]<-thisBiomass*thisPreyAvail12
      juvPrey[[thisPreyCode]]<-thisBiomass*thisPreyAvail11
    }
  }
}

#plot it. Get plot limits first
getCode<-function(x){
  yy<-gsub("\\d","",x,perl=TRUE)
  return(yy)
}
getAge<-function(x){
  yy<-gsub("\\d","",x,perl=TRUE)
  yyy<-gsub(yy,"",x)
  return(yyy)
}
thisBiomassLimScalar<-0.01 #use this if want to zoom in and see smaller prey groups - ie cut out zooplankton
#set to 1 if want to keep all in
np<-length(juvPrey);
colorByGroup<-colorRampPalette(colors=c(myGold,myOrange,"red",myRed,myPurple,myBlue,myAqua,myGreen,myDarkGreen,"black",myGrey))(np)
#can use total food for ylim
thisYlim<-max(howMuchFoodJuv[[as.character(thisCode)]],na.rm=TRUE)
pdf(paste(plotPath,thisCode,"PreyAvailableAndBiomass_juv.pdf",sep=""))
par(mfrow=c(2,1),mar=c(3,4.5,1,1),oma=c(5,1,1,1))
for(l in 1:nlayers){
  doPlot<-FALSE
  par(xpd=FALSE)
  plot(x=seq(1,nboxes),y=rep(0,nboxes),ylim=c(0,thisYlim),type="n",ylab="Prey available (mg N)",xlab="")
  storeY<-rep(0,nboxes)
  legendPreys<-c(); legendColors<-c()
  for(p in 1:length(juvPrey)){
    thisValues<-juvPrey[[p]][l,]
    
    for(b in 1:nboxes){
      polygon(x=c(b-0.5,b-0.5,b+0.5,b+0.5),y=c(storeY[b],storeY[b]+thisValues[b],storeY[b]+thisValues[b],storeY[b]),col=colorByGroup[p])
    }
    storeY<-storeY+thisValues
  }
  #for those with positive prey values, plot biomass for each box layer
  #get max first
  thisBiomassLim<-0
  for(p in 1:length(juvPrey)){
    thisValues<-juvPrey[[p]][l,]
    if(sum(thisValues,na.rm=TRUE)>0){
      legendPreys<-c(legendPreys,names(juvPrey)[p]); legendColors<-c(legendColors,colorByGroup[p])
      thisPreyCode<-getCode(names(juvPrey)[p])
      thisPreyAge<-getAge(names(juvPrey)[p])
      thisPreyCM<-groupsDF$CohortMature[groupsDF$Code==thisPreyCode]
      thisBiomass<-allBiomass[[thisPreyCode]]
      if(thisPreyNumCohorts>1){
        if(thisPreyAge=="1"){
          temp<-thisBiomass[,,1:thisPreyCM]
          if(length(dim(temp))==3){
            thisPreyBiomass<-apply(temp,c(1,2),sum)
          }else{
            thisPreyBiomass<-temp
          }
        } else{
          temp<-thisBiomass[,,(thisPreyCM+1):thisPreyNumCohorts]
          if(length(dim(temp))==3){
            thisPreyBiomass<-apply(temp,c(1,2),sum)
          }else{
            thisPreyBiomass<-temp
          }
        }
          
      } else{
        thisPreyBiomass<-thisBiomass
      }
      thisBiomassLim<-max(thisBiomassLim,max(thisPreyBiomass,na.rm=TRUE))
    }
  }
  par(xpd=FALSE)
  plot(x=seq(1,nboxes),y=rep(0,nboxes),ylim=c(0,thisBiomassLim*thisBiomassLimScalar),type="n",ylab="Prey biomass (mg N)",xlab="")
  for(p in 1:length(juvPrey)){
    thisValues<-juvPrey[[p]][l,]
    if(sum(thisValues,na.rm=TRUE)>0){
      doPlot<-TRUE
      thisPreyCode<-getCode(names(juvPrey)[p])
      thisPreyAge<-getAge(names(juvPrey)[p])
      thisPreyCM<-groupsDF$CohortMature[groupsDF$Code==thisPreyCode]; thisPreyNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisPreyCode]
      thisBiomass<-allBiomass[[thisPreyCode]]
      if(thisPreyNumCohorts>1){
        if(thisPreyAge=="1"){
          temp<-thisBiomass[,,1:thisPreyCM]
          if(length(dim(temp))==3){
            thisPreyBiomass<-apply(temp,c(1,2),sum)
          }else{
            thisPreyBiomass<-temp
          }
        } else{
          temp<-thisBiomass[,,(thisPreyCM+1):thisPreyNumCohorts]
          if(length(dim(temp))==3){
            thisPreyBiomass<-apply(temp,c(1,2),sum)
          }else{
            thisPreyBiomass<-temp
          }
        }
        
      } else{
        thisPreyBiomass<-thisBiomass
      }
      points(x=seq(1,nboxes),y=thisPreyBiomass[l,],type="l",lwd=2,col=colorByGroup[p])
    }
  }
  if(doPlot==TRUE){
    par(xpd=NA)
    legend(legend=legendPreys,col=legendColors,x="bottom",inset=-0.7,lwd=2,ncol=5)
  }
}
dev.off()
### ADULTS
np<-length(adultPrey);
colorByGroup<-colorRampPalette(colors=c(myGold,myOrange,"red",myRed,myPurple,myBlue,myAqua,myGreen,myDarkGreen,"black",myGrey))(np)
#can use total food for ylim
thisYlim<-max(howMuchFoodAdults[[as.character(thisCode)]],na.rm=TRUE)
pdf(paste(plotPath,thisCode,"PreyAvailableAndBiomass_adults.pdf",sep=""))
par(mfrow=c(2,1),mar=c(3,4.5,1,1),oma=c(5,1,1,1))
for(l in 1:nlayers){
  doPlot<-FALSE
  par(xpd=FALSE)
  plot(x=seq(1,nboxes),y=rep(0,nboxes),ylim=c(0,thisYlim),type="n",ylab="Prey available (mg N)",xlab="")
  storeY<-rep(0,nboxes)
  legendPreys<-c(); legendColors<-c()
  for(p in 1:length(adultPrey)){
    thisValues<-adultPrey[[p]][l,]
    
    for(b in 1:nboxes){
      polygon(x=c(b-0.5,b-0.5,b+0.5,b+0.5),y=c(storeY[b],storeY[b]+thisValues[b],storeY[b]+thisValues[b],storeY[b]),col=colorByGroup[p])
    }
    storeY<-storeY+thisValues
  }
  #for those with positive prey values, plot biomass for each box layer
  #get max first
  thisBiomassLim<-0
  for(p in 1:length(adultPrey)){
    thisValues<-adultPrey[[p]][l,]
    if(sum(thisValues,na.rm=TRUE)>0){
      legendPreys<-c(legendPreys,names(adultPrey)[p]); legendColors<-c(legendColors,colorByGroup[p])
      thisPreyCode<-getCode(names(adultPrey)[p])
      thisPreyAge<-getAge(names(adultPrey)[p])
      thisPreyCM<-groupsDF$CohortMature[groupsDF$Code==thisPreyCode]
      thisBiomass<-allBiomass[[thisPreyCode]]
      if(thisPreyNumCohorts>1){
        if(thisPreyAge=="1"){
          temp<-thisBiomass[,,1:thisPreyCM]
          if(length(dim(temp))==3){
            thisPreyBiomass<-apply(temp,c(1,2),sum)
          }else{
            thisPreyBiomass<-temp
          }
        } else{
          temp<-thisBiomass[,,(thisPreyCM+1):thisPreyNumCohorts]
          if(length(dim(temp))==3){
            thisPreyBiomass<-apply(temp,c(1,2),sum)
          }else{
            thisPreyBiomass<-temp
          }
        }
        
      } else{
        thisPreyBiomass<-thisBiomass
      }
      thisBiomassLim<-max(thisBiomassLim,max(thisPreyBiomass,na.rm=TRUE))
    }
  }
  par(xpd=FALSE)
  plot(x=seq(1,nboxes),y=rep(0,nboxes),ylim=c(0,thisBiomassLim*thisBiomassLimScalar),type="n",ylab="Prey biomass (mg N)",xlab="")
  for(p in 1:length(juvPrey)){
    thisValues<-adultPrey[[p]][l,]
    if(sum(thisValues,na.rm=TRUE)>0){
      doPlot<-TRUE
      thisPreyCode<-getCode(names(adultPrey)[p])
      thisPreyAge<-getAge(names(adultPrey)[p])
      thisPreyCM<-groupsDF$CohortMature[groupsDF$Code==thisPreyCode]; thisPreyNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisPreyCode]
      thisBiomass<-allBiomass[[thisPreyCode]]
      if(thisPreyNumCohorts>1){
        if(thisPreyAge=="1"){
          temp<-thisBiomass[,,1:thisPreyCM]
          if(length(dim(temp))==3){
            thisPreyBiomass<-apply(temp,c(1,2),sum)
          }else{
            thisPreyBiomass<-temp
          }
        } else{
          temp<-thisBiomass[,,(thisPreyCM+1):thisPreyNumCohorts]
          if(length(dim(temp))==3){
            thisPreyBiomass<-apply(temp,c(1,2),sum)
          }else{
            thisPreyBiomass<-temp
          }
        }
        
      } else{
        thisPreyBiomass<-thisBiomass
      }
      points(x=seq(1,nboxes),y=thisPreyBiomass[l,],type="l",lwd=2,col=colorByGroup[p])
    }
  }
  if(doPlot==TRUE){
    par(xpd=NA)
    legend(legend=legendPreys,col=legendColors,x="bottom",inset=-0.7,lwd=2,ncol=5)
  }
}
dev.off()

## and focus plots on predators of a given prey and compare with biomass (or numbers..?)
# thisCode<-"ZL"

gcodes5
gcodeNumbers5

rbind(gcodes4,gcodeNumbers4)

