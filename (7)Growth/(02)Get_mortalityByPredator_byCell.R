#read in initial conditions and availability matrix then calculate how much of each prey is eaten (in mg N then numbers) by predator by cell (box/layer)
##this has invert predators as well ((02)Get_availableFood_byCell.R doesn't)

basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

plotPath<-paste(basePath,"Figures\\growth\\",sep="") #this is created in DIET_read_andCreateCsv.R

#created in DIET_read_andCreateCsv.R in setting_up/Testing/ RERUN this if any changes to pPREY 
pPREY_df<-read.csv(paste(basePath,"pPREY_OUT.csv",sep="")) 

biolFile<-paste(basePath,"CRAM_base_biol.prm",sep=""); biolLines<-readLines(biolFile)

initFile<-paste(basePath,"CRAM_input_short.nc",sep=""); initiNC<-nc_open(initFile); 
initFileText<-paste(basePath,"CRAM_input_short.txt",sep=""); initLines<-readLines(initFileText)

groupsFile<-paste(basePath,"CRAM_groups.csv",sep=""); groupsDF<-read.csv(groupsFile); ng<-dim(groupsDF)[1]

max_nCohorts<-10; nboxes<-30; nlayers<-6

volume<-ncvar_get(initiNC,"volume")[,,1]

mumByGroupCohort<-read.csv(paste(basePath,"mumByCohort.csv",sep=""),header=TRUE)

#turn NAs to zeros
NAs2zeros<-function(x){
  y<-x
  if(is.na(x)){y<-0}
  return(y)
}

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
howMuchFoodAdults<-NULL; howMuchFoodJuv<-NULL; howMuchFoodInverts<-NULL
for(g in 1:ng){
  thisCode<-as.character(groupsDF$Code[g]); thisNumCohorts<-groupsDF$NumCohorts[g]
  thispPrey<-pPREY_df[pPREY_df$PredatorCode==thisCode,]
  if(thisNumCohorts>1){
    #need to do seperate diet for adults and juveniles
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
  } else{
    invertPrey<-zeroDF; 
    for(p in 1:ng){
      thisPreyCode<-as.character(groupsDF$Code[p]); thisPreyNumCohorts<-groupsDF$NumCohorts[p]; 
      thisBiomass<-allBiomass[[thisPreyCode]]
      if(thisPreyNumCohorts>1){
        temp<-thispPrey[,c(thisPreyCode)]
        if(length(temp)>0){
          if(temp>0){
            cat("Case of invert eating age-structured. Predator ",thisCode,", prey ",thisPreyCode,"\n")
          }
        }
      } else {
        thisPreyAvail<-thispPrey[,c(thisPreyCode)]
        if(length(dim(thisBiomass))==2){
          invertPrey<-invertPrey+thisBiomass*thisPreyAvail
         }
      }
    }
    howMuchFoodInverts[[as.character(thisCode)]]<-invertPrey; 
  }
}

##get biomass eaten for all groups, and numbers eaten for age structured groups

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

biomassEaten11<-numbersEaten11; biomassEaten12<-numbersEaten12; biomassEaten22<-numbersEaten22; biomassEaten21<-numbersEaten21

# trackBiomass<-allBiomass; trackDensity<-allDens #when prey are eaten, we'll remove from here so we know what's left

for(g in 1:ng){
  thisCode<-as.character(groupsDF$Code[g]); thisName<-str_trim(groupsDF$Name[g],side="both"); thisNumCohorts<-groupsDF$NumCohorts[g]
  cat(as.character(thisCode),"--")
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
          biomassEaten22[g,p,,]<-amountEaten
          ##
          ## adults eating juveniles (prey then predator)
          amountEaten<-thisPreyAvail12*juvBiomass*as.double(adPropToEat)
          #number of individuals required (round up)
          numberEaten<-amountEaten / juvPreyWeightPerInd
          numbersEaten12[g,p,,]<-numberEaten
          biomassEaten12[g,p,,]<-amountEaten
          ##
          ## juveniles eating juveniles (prey then predator)
          amountEaten<-thisPreyAvail11*juvBiomass*as.double(juvPropToEat)
          #number of individuals required (round up)
          numberEaten<-amountEaten / juvPreyWeightPerInd
          numbersEaten11[g,p,,]<-numberEaten
          biomassEaten11[g,p,,]<-amountEaten
          ##
          ## juveniles eating adults (prey then predator)
          amountEaten<-thisPreyAvail21*adBiomass*as.double(juvPropToEat)
          #number of individuals required (round up)
          numberEaten<-amountEaten / adPreyWeightPerInd
          numbersEaten21[g,p,,]<-numberEaten
          biomassEaten21[g,p,,]<-amountEaten
        }
      } else{
        thisPreyAvail22<-0
        thisPreyAvail12<-max(thispPrey[thispPrey$PredatorAge==2,c(thisPreyCode)])
        thisPreyAvail21<-0
        thisPreyAvail11<-max(thispPrey[thispPrey$PredatorAge==1,c(thisPreyCode)])
        if(sum(thisPreyAvail11,thisPreyAvail12,thisPreyAvail21,thisPreyAvail22)>0){
          #amount (mg N) eaten 
          ##adult predators
          if(length(dim(thisBiomass))==2){
            amountEaten<-thisPreyAvail12*thisBiomass*as.double(adPropToEat)
          } else{
            tempBiomass<-array(0,dim=c(nlayers,nboxes))
            tempBiomass[nlayers,]<-thisBiomass
            amountEaten<-thisPreyAvail12*tempBiomass*as.double(adPropToEat)
          }
          biomassEaten12[g,p,,]<-amountEaten
          biomassEaten22[g,p,,]<-0*amountEaten
          ##juvenile predators
          if(length(dim(thisBiomass))==2){
            amountEaten<-thisPreyAvail11*thisBiomass*as.double(juvPropToEat)
          } else{
            tempBiomass<-array(0,dim=c(nlayers,nboxes))
            tempBiomass[nlayers,]<-thisBiomass
            amountEaten<-thisPreyAvail11*tempBiomass*as.double(juvPropToEat)
          }
          biomassEaten11[g,p,,]<-amountEaten
          biomassEaten21[g,p,,]<-0*amountEaten
        }
      }
    }
  }else{ #this is non-agestructured predators (invert)
    #to get mum for inverts, read it in from biol.prm file
    thisVar<-paste("mum_",thisCode,sep=""); x<-grep(thisVar,biolLines); thisLine<-biolLines[x]
    #if 'T15" exists in the line, take it out. Then get first number
    thisLine<-gsub("T15","",thisLine); thisMum<-get_first_number(thisLine)
    #need to read in thisInvertE too 
    thisVar<-paste("E_",thisCode,sep=""); x<-grep(thisVar,biolLines); thisLine<-biolLines[x]
    thisInvertE<-get_first_number(thisLine)
    #if there is no E, they are not eating (like bacteria, so skipt)
    if(length(thisInvertE)>0){
      thisBiomass<-allBiomass[[as.character(thisCode)]]
      thisTotalPrey<-howMuchFoodInverts[[thisCode]]
      #if there are no prey, skip
      if(nrow(thisTotalPrey)>0){
        thisFoodRequired<-thisBiomass*(thisMum/thisInvertE)
        #if 1D, need to fill out other cells
        if(length(dim(thisFoodRequired))==0){
          temp<-thisFoodRequired; thisFoodRequired<-0*thisTotalPrey; thisFoodRequired[nlayers,]<-temp
        }
        invertPropToEat<-0*thisFoodRequired;  
        for(i in 1:nrow(invertPropToEat)){
          for(j in 1:ncol(invertPropToEat)){
            invertPropToEat[i,j]<-getFoodProportion(foodReq=thisFoodRequired[i,j],foodAvail=thisTotalPrey[i,j])
          }
        }
        if("missing" %in% invertPropToEat){cat("Food missing for ",thisCode,"\n")}
        thispPrey<-pPREY_df[pPREY_df$PredatorCode==thisCode,]
        #now apply this proportion to each prey and turn it into numbers and store that
        for(p in 1:ng){
          thisPreyCode<-as.character(groupsDF$Code[p]); thisPreyNumCohorts<-groupsDF$NumCohorts[p]; 
          thisBiomass<-allBiomass[[thisPreyCode]]
          if(thisPreyNumCohorts==1){
            thisPreyAvail<-max(thispPrey[,c(thisPreyCode)]) #there should only be one value, but just in case there is not, take max
             if(thisPreyAvail>0){
              #amount (mg N) eaten 
              ##adult predators
              if(length(dim(thisBiomass))==2){
                amountEaten<-thisPreyAvail*thisBiomass*as.double(adPropToEat)
              } else{
                tempBiomass<-array(0,dim=c(nlayers,nboxes))
                tempBiomass[nlayers,]<-thisBiomass
                amountEaten<-thisPreyAvail*tempBiomass*as.double(adPropToEat)
              }
              biomassEaten11[g,p,,]<-amountEaten
            }
            
          }
        }
      }
    }
  }
}

#for each group, plot biomass eatn by predator by cell. Include a line with initial biomass of the prey group
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
  thisBiomass<-allBiomass[[as.character(thisCode)]]
  if(thisNumCohorts>1){
    #need to split mature and juv
    thisCM<-groupsDF$CohortMature[g]
    #get adult and juv biomass
    juvBiomass<-apply(thisBiomass[,,1:thisCM],c(1,2),sum,na.rm=TRUE); adBiomass<-apply(thisBiomass[,,(thisCM+1):thisNumCohorts],c(1,2),sum,na.rm=TRUE)
    #get adult and juv biomass eaten
    temp11<-apply(biomassEaten11[,g,,],seq(1,length(dim(biomassEaten11[,g,,]))),NAs2zeros)
    temp12<-apply(biomassEaten12[,g,,],seq(1,length(dim(biomassEaten12[,g,,]))),NAs2zeros)
    temp21<-apply(biomassEaten21[,g,,],seq(1,length(dim(biomassEaten21[,g,,]))),NAs2zeros)
    temp22<-apply(biomassEaten22[,g,,],seq(1,length(dim(biomassEaten22[,g,,]))),NAs2zeros)
    
    adBiomassEaten<-temp21+temp22
    juvBiomassEaten<-temp11+temp12
    if(sum(juvBiomassEaten,na.rm=TRUE)>0){
      pdf(paste(plotPath,thisCode,"juvPredationByCellandPredator.pdf",sep=""),height=9)
      par(mar=c(0,5,1,1),mfrow=c(nlayers,1),oma=c(9,3,1,1),las=1)
      legendPreds<-c(); 
      for(l in 1:nlayers){
        thisJuvBiomassEaten<-juvBiomassEaten[,l,]; thisYlim<-max(colSums(thisJuvBiomassEaten,na.rm=TRUE)*365,max(juvBiomass[l,],na.rm=TRUE))
        storeY<-rep(0,nboxes)
        par(xpd=FALSE)
        plot(x=seq(1,nboxes),y=rep(0,nboxes),ylim=c(0,thisYlim),type="n",ylab="",xlab="",xaxt="n")
        for(p in 1:ng){
          thisValues<-thisJuvBiomassEaten[p,]*365; thisValues[is.na(thisValues)]<-0
          if(sum(thisValues)>0){
            legendPreds<-c(legendPreds,as.character(groupsDF$Code[p])); 
          }
          for(b in 1:nboxes){
            polygon(x=c(b-0.5,b-0.5,b+0.5,b+0.5),y=c(storeY[b],storeY[b]+thisValues[b],storeY[b]+thisValues[b],storeY[b]),col=colorByGroup[p])
          }
          storeY<-storeY+thisValues
          # cat(p,", ", max(storeY),"--")
        }
        points(x=seq(1,nboxes),y=juvBiomass[l,],type="l",col="midnightblue",lwd=2)
        
        maxPercEaten<-max((storeY/juvBiomass[l,])[juvBiomass[l,]>0])
        mtext(paste("Max. proportion eaten: ", signif(maxPercEaten,2),sep=""),side=3,adj=0,line=-1.2)
      }
      axis(at=seq(2,nboxes,by=2),labels = seq(2,nboxes,by=2),side=1,outer=TRUE)
      mtext("Box number",side=1,adj=0.5,line=2.5,outer=TRUE)
      mtext("Biomass eaten per year (mg N)",side=2,adj=0.5,outer=TRUE,las=0)
      par(xpd=NA)
      legendPreds<-unique(sort(legendPreds)); legendColors<-colorByGroup[match(legendPreds,groupsDF$Code)]
      par(lend=1)
      legend(legend=legendPreds,col=legendColors,x="bottom",inset=-0.9,lwd=5,ncol=8)
      dev.off()
    }
    ##adults
    if(sum(adBiomassEaten,na.rm=TRUE)>0){
      pdf(paste(plotPath,thisCode,"adPredationByCellandPredator.pdf",sep=""),height=9)
      par(mar=c(0,5,1,1),mfrow=c(nlayers,1),oma=c(9,3,1,1),las=1)
      legendPreds<-c(); 
      for(l in 1:nlayers){
        thisAdBiomassEaten<-adBiomassEaten[,l,]; thisYlim<-max(colSums(thisAdBiomassEaten,na.rm=TRUE)*365,max(adBiomass[l,],na.rm=TRUE))
        storeY<-rep(0,nboxes)
        par(xpd=FALSE)
        plot(x=seq(1,nboxes),y=rep(0,nboxes),ylim=c(0,thisYlim),type="n",ylab="",xlab="",xaxt="n")
        for(p in 1:ng){
          thisValues<-thisAdBiomassEaten[p,]*365; thisValues[is.na(thisValues)]<-0
          if(sum(thisValues)>0){
            legendPreds<-c(legendPreds,as.character(groupsDF$Code[p])); 
          }
          for(b in 1:nboxes){
            polygon(x=c(b-0.5,b-0.5,b+0.5,b+0.5),y=c(storeY[b],storeY[b]+thisValues[b],storeY[b]+thisValues[b],storeY[b]),col=colorByGroup[p])
          }
          storeY<-storeY+thisValues
          # cat(p,", ", max(storeY),"--")
        }
        points(x=seq(1,nboxes),y=adBiomass[l,],type="l",col="midnightblue",lwd=2)
        
        maxPercEaten<-max((storeY/adBiomass[l,])[adBiomass[l,]>0])
        mtext(paste("Max. proportion eaten: ", signif(maxPercEaten,2),sep=""),side=3,adj=0,line=-1.2)
      }
      axis(at=seq(2,nboxes,by=2),labels = seq(2,nboxes,by=2),side=1,outer=TRUE)
      mtext("Box number",side=1,adj=0.5,line=2.5,outer=TRUE)
      mtext("Biomass eaten per year (mg N)",side=2,adj=0.5,outer=TRUE,las=0)
      par(xpd=NA)
      legendPreds<-unique(sort(legendPreds)); legendColors<-colorByGroup[match(legendPreds,groupsDF$Code)]
      par(lend=1)
      legend(legend=legendPreds,col=legendColors,x="bottom",inset=-0.9,lwd=5,ncol=8)
      dev.off()
    }
  } else{

    temp11<-apply(biomassEaten11[,g,,],seq(1,length(dim(biomassEaten11[,g,,]))),NAs2zeros)
    temp12<-apply(biomassEaten12[,g,,],seq(1,length(dim(biomassEaten12[,g,,]))),NAs2zeros)
    temp21<-apply(biomassEaten21[,g,,],seq(1,length(dim(biomassEaten21[,g,,]))),NAs2zeros)
    temp22<-apply(biomassEaten22[,g,,],seq(1,length(dim(biomassEaten22[,g,,]))),NAs2zeros)
    thisBiomassEaten<-temp11+temp12+temp21+temp22

    if(sum(thisBiomassEaten,na.rm=TRUE)>0){
      pdf(paste(plotPath,thisCode,"PredationByCellandPredator.pdf",sep=""),height=9)
      par(mar=c(0,5,1,1),mfrow=c(nlayers,1),oma=c(9,3,1,1),las=1)
      legendPreds<-c(); 
      for(l in 1:nlayers){
        thisJuvBiomassEaten<-thisBiomassEaten[,l,]; thisYlim<-max(colSums(thisJuvBiomassEaten,na.rm=TRUE)*365,max(thisBiomass[l,],na.rm=TRUE))
        storeY<-rep(0,nboxes)
        par(xpd=FALSE)
        plot(x=seq(1,nboxes),y=rep(0,nboxes),ylim=c(0,thisYlim),type="n",ylab="",xlab="",xaxt="n")
        for(p in 1:ng){
          thisValues<-thisJuvBiomassEaten[p,]*365; thisValues[is.na(thisValues)]<-0
          if(sum(thisValues)>0){
            legendPreds<-c(legendPreds,as.character(groupsDF$Code[p])); 
          }
          for(b in 1:nboxes){
            polygon(x=c(b-0.5,b-0.5,b+0.5,b+0.5),y=c(storeY[b],storeY[b]+thisValues[b],storeY[b]+thisValues[b],storeY[b]),col=colorByGroup[p])
          }
          storeY<-storeY+thisValues
          # cat(p,", ", max(storeY),"--")
        }
        points(x=seq(1,nboxes),y=thisBiomass[l,],type="l",col="midnightblue",lwd=2)
        
        maxPercEaten<-max((storeY/thisBiomass[l,])[thisBiomass[l,]>0])
        mtext(paste("Max. proportion eaten: ", signif(maxPercEaten,2),sep=""),side=3,adj=0,line=-1.2)
      }
      axis(at=seq(2,nboxes,by=2),labels = seq(2,nboxes,by=2),side=1,outer=TRUE)
      mtext("Box number",side=1,adj=0.5,line=2.5,outer=TRUE)
      mtext("Biomass eaten per year (mg N)",side=2,adj=0.5,outer=TRUE,las=0)
      par(xpd=NA)
      legendPreds<-unique(sort(legendPreds)); legendColors<-colorByGroup[match(legendPreds,groupsDF$Code)]
      par(lend=1)
      legend(legend=legendPreds,col=legendColors,x="bottom",inset=-0.9,lwd=5,ncol=8)
      dev.off()
    }
  }
  ##NEED TO DO A VERSION THAT ADDS IN RECRUITMENT (OR GROWTH FOR THE UN AGE STRUCTURED)
  
}

##update spatial distribution of zoo plankton groups based on predation
updateZooIC<-FALSE #set this so don't update by mistake
small_number<-1e-8
ni<-dim(volume)[1]; nj<-dim(volume)[2]
if(updateZooIC==TRUE){
  Zcodes<-c("ZL","ZM","ZG","ZS"); nz<-length(Zcodes)
  for(z in 1:nz){
    thisCode<-zcodes[z]
    groupIndex<-groupsDF$Code==thisCode
    thisName<-str_trim(groupsDF$Name[groupIndex]); thisVar<-paste(thisName,"_N",sep="")
    thisIC<-ncvar_get(initiNC,thisVar)[,,1]
    thisIC_biomass<-thisIC*volume
    #get predation by cell
    temp11<-apply(biomassEaten11[,groupIndex,,],seq(1,length(dim(biomassEaten11[,groupIndex,,]))),NAs2zeros)
    temp12<-apply(biomassEaten12[,groupIndex,,],seq(1,length(dim(biomassEaten12[,groupIndex,,]))),NAs2zeros)
    temp21<-apply(biomassEaten21[,groupIndex,,],seq(1,length(dim(biomassEaten21[,groupIndex,,]))),NAs2zeros)
    temp22<-apply(biomassEaten22[,groupIndex,,],seq(1,length(dim(biomassEaten22[,groupIndex,,]))),NAs2zeros)
    thisBiomassEaten<-apply(temp11+temp12+temp21+temp22,c(2,3),sum)
    #max biomassEaten
    thisMaxEatenPerBiomass<-0; temp<-thisBiomassEaten/thisIC_biomass; thisMinEatenPerBiomass<-1e+10
    for(i in 1:ni){
      for(j in 1:nj){
        if(temp[i,j]!="NaN"){
          if(temp[i,j]>0){
            if(temp[i,j]>thisMaxEatenPerBiomass){thisMaxEatenPerBiomass<-temp[i,j]}
            if(temp[i,j]<thisMinEatenPerBiomass){thisMinEatenPerBiomass<-temp[i,j]}
          }
        }
      }
    }
    biomassPerm3_max<-max(thisIC); biomassPerm3_min<-biomassPerm3_max*(thisMinEatenPerBiomass/thisMaxEatenPerBiomass)
    newIC<-thisIC
    for(i in 1:ni){
      for(j in 1:nj){
        newIC[i,j]<-biomassPerm3_max*(thisBiomassEaten[i,j]/thisMaxEatenPerBiomass)
      }
    }
    thisIC_biomass<-newIC*volume

    # par(mfrow=c(nlayers,1))
    # for(l in 1:nlayers){
    #   plot(thisIC_biomass[l,],type="h",lwd=5,col=myGreen,lend=1)
    #   par(new=TRUE)
    #   plot(thisBiomassEaten[l,],col=myOrange,pch=8,xlab="",ylab="",xaxt="n",yaxt="n")
    # }
    
    #write over IC lines
    
    
  }
  
  
}



