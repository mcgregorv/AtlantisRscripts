#read in initial diets and plot intended predators for each prey - just uses proportion of diet with no scaling for predator abundance or prey abundance

basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
inputsPath<-paste(basePath, "\\inputs\\",sep="")
plotPath<-paste(DIR$'Figures',"Mortality\\", sep="")

diets<-read.csv(paste(inputsPath,"biol_prm\\pPREY\\CRAM_biol_pPREY_original.csv", sep=""))
preys<-diets[,1]; npreys<-length(preys)
getAgeColor<-function(x){
  thisCol<-myGrey
  if(is.na(x)){
    thisCol<-myGrey
  } else{
    if(x==2){thisCol="black"}
  }
  return(thisCol)
}
getAgeText<-function(x){
  thisText<-"Juv"
  if(is.na(x)){
    thisText<-""
  }else{
    if(x==2){thisText="Adult"}
  }
  return(thisText)
}

pdf(paste(plotPath,"AllIntendedPredators.pdf", sep="")); par(mfrow=c(3,2), mar=c(6,4,1,1), oma=c(0,0,0,0))
for(p in 1:npreys){
  temp<-diets[p,-1]
  preyCode<-preys[p]
  index<-temp>0
  thisPreyAvails<-temp[index]; thisPreds<-colnames(diets)[-1][index]
  if(length(thisPreyAvails)>0){
    ## do a juvenile and an adult plot. include not aged in juveniles
    juvIndex<-grep("pPREY2", thisPreds, invert=TRUE); adultIndex<-grep("pPREY2", thisPreds)
    juvPreds<-thisPreds[juvIndex]; juvAvails<-thisPreyAvails[juvIndex]
    predAges<-unlist(lapply(juvPreds, get_first_number, n=2));
    if(sum(predAges, na.rm=TRUE)==0){
      predAgeCols<-rep(myGrey, length(predAges))
      predAgeText<-rep("", length(predAges))
    }else{
      predAgeCols<-unlist(lapply(predAges, getAgeColor))
      predAgeText<-unlist(lapply(predAges, getAgeText))
    }
    if(length(juvPreds)>0){
      predCodes<-gsub("pPREY|1|2","", juvPreds); predLabels<-paste(predCodes, predAgeText, sep=" ")
      par(las=1, lend=1)
      plot(x=seq(1,length(juvPreds)), y=juvAvails, col=predAgeCols, type="h", lwd=5, xaxt="n", ylab="Proportion of diet", xlab="")
      mtext("Juvenile", side=3, adj=0.01,line=-1); mtext(preyCode,side=3,adj=0,font=2)
      par(las=2)
      axis(at=seq(1,length(juvPreds)), labels=predLabels, side=1)
    }
    if(length(adultIndex)>0){
      adPreds<-thisPreds[adultIndex]; adAvails<-thisPreyAvails[adultIndex]
      
      
      predAges<-unlist(lapply(adPreds, get_first_number, n=2)); predAgeCols<-unlist(lapply(predAges, getAgeColor))
      predAgeText<-unlist(lapply(predAges, getAgeText))
      predCodes<-gsub("pPREY|1|2","", adPreds); predLabels<-paste(predCodes, predAgeText, sep=" ")
      par(las=1, lend=1)
      plot(x=seq(1,length(adPreds)), y=adAvails, col=predAgeCols, type="h", lwd=5, xaxt="n", ylab="Proportion of diet", xlab="")
      mtext("Adult", side=3, adj=0.01,line=-1); mtext(preyCode,side=3,adj=0,font=2)
      par(las=2)
      axis(at=seq(1,length(adPreds)), labels=predLabels, side=1)
    }
  }
}
dev.off()


###########################################################################
## realised predators from a model run
this_run<-"XXX_FSMGDistribB3"; out_path<-paste(basePath,"base\\output",this_run,"\\", sep="")
dietCheck<-read.csv(paste(out_path,"outputDietCheck.txt", sep=""), sep=" ")
ThisNC.nc<-nc_open(paste(out_path,"output.nc", sep="")); thisVol<-ncvar_get(ThisNC.nc, "volume")

preys<-colnames(dietCheck)[6:dim(dietCheck)[2]]; npreys<-length(preys)

##need biol.prm for age mature, and groupsDF for number of cohorts
groupsDF<-read.csv(paste(basePath,"CRAM_groups.csv", sep="")); 
biolLines<-readLines(paste(basePath,"CRAM_BH_hybrid_biol.prm", sep=""))

pdf(paste(plotPath,"AllRealisedPredationMortality.pdf", sep="")); par(mfrow=c(3,2), mar=c(6,4,1,1), oma=c(0,0,0,0))
for(p in 1:npreys){
  thisPrey<-preys[p]; thisData<-dietCheck[,c(1:3,grep(thisPrey, colnames(dietCheck)))]
  tempData<-tapply(thisData[,c(thisPrey)], list(thisData$Time, thisData$Predator), sum, na.rm=TRUE)
  index<-colSums(tempData, na.rm=TRUE)>0
  thisPredData<-tempData[,index]
  if(length(dim(thisPredData))==0){
    sumPredData<-sum(thisPredData); predLabels<-colnames(tempData)[index]
  }else{
    sumPredData<-apply(thisPredData,2,sum,na.rm=TRUE); predLabels<-names(sumPredData)
  }
  if(sum(sumPredData)>0){
    par(las=1, lend=1)
    plot(sumPredData, type="h", lwd=5, xaxt="n", xlab="", ylab="Proportion predation mortality")
    mtext(thisPrey,side=3, adj=0, font=2)
    par(las=2)
    axis(at=seq(1,length(sumPredData)), labels =predLabels, side=1)
  }
}
dev.off()

## get realised biomass of the predators and plot proprtion mortalities divided by predator biomass so comparible with proportion of diets
## similar to proportion of diet multiplied by predator abundance, but quicker to do here
getBiomass<-function(code){
  thisPredName<-str_trim(groupsDF$Name[groupsDF$Code==code], side="both")
  thisTracer<-paste(thisPredName, "_N", sep=""); thisTracerData<-ncvar_get(ThisNC.nc, thisTracer); 
  if(length(dim(thisTracerData))==3){  
    thisBiomass<-apply(thisTracerData*thisVol, 3, sum)*X_CN * mg_2_tonne
  } else{
    thisBiomass<-apply(thisTracerData*thisVol[6,,], 2, sum)*X_CN * mg_2_tonne
  }
  return(thisBiomass)
}

pdf(paste(plotPath,"AllRealisedPredationMortalityDivByPredBiomass.pdf", sep="")); par(mfrow=c(3,2), mar=c(6,4,1,1), oma=c(0,0,0,0))
for(p in 1:npreys){
  thisPrey<-preys[p]; thisData<-dietCheck[,c(1:3,grep(thisPrey, colnames(dietCheck)))]
  tempData<-tapply(thisData[,c(thisPrey)], list(thisData$Time, thisData$Predator), sum, na.rm=TRUE)
  index<-colSums(tempData, na.rm=TRUE)>0
  thisPredData<-tempData[,index]
  if(sum(thisPredData, na.rm=TRUE)>0){
    if(length(dim(thisPredData))==0){
      predLabels<-colnames(tempData)[index]
      thisBiomass<-getBiomass(predLabels)
      sumPredData<-sum(thisPredData / thisBiomass[1:length(thisPredData)]); 
    }else{
      predLabels<-colnames(tempData)[index]
      test<-lapply(predLabels, getBiomass)
      sumPredData<-rep(0, length(predLabels))
      for(i in 1:length(predLabels)){
        xx<-thisPredData[,i] / test[[i]][1:dim(thisPredData)[1]]
        sumPredData[i]<-sum(xx, na.rm=TRUE)
      }
    }
    if(sum(sumPredData)>0){
      par(las=1, lend=1)
      plot(sumPredData, type="h", lwd=5, xaxt="n", xlab="", ylab="Est. realised proportion predators diet")
      mtext(thisPrey,side=3, adj=0, font=2)
      par(las=2)
      axis(at=seq(1,length(sumPredData)), labels =predLabels, side=1)
    }
  }
}
dev.off()

pdf(paste(plotPath,"AllRealisedPredationPredBiomass.pdf", sep="")); par(mfrow=c(3,2), mar=c(6,4,1,1), oma=c(0,0,0,0))
for(p in 1:npreys){
  thisPrey<-preys[p]; thisData<-dietCheck[,c(1:3,grep(thisPrey, colnames(dietCheck)))]
  tempData<-tapply(thisData[,c(thisPrey)], list(thisData$Time, thisData$Predator), sum, na.rm=TRUE)
  index<-colSums(tempData, na.rm=TRUE)>0
  thisPredData<-tempData[,index]
  if(sum(thisPredData, na.rm=TRUE)>0){
    if(length(dim(thisPredData))==0){
      predLabels<-colnames(tempData)[index]
      thisBiomass<-getBiomass(predLabels)
      sumPredData<-sum(thisBiomass[1:length(thisPredData)]); 
    }else{
      predLabels<-colnames(tempData)[index]
      test<-lapply(predLabels, getBiomass)
      sumPredData<-rep(0, length(predLabels))
      for(i in 1:length(predLabels)){
        xx<-test[[i]][1:dim(thisPredData)[1]]
        sumPredData[i]<-sum(xx, na.rm=TRUE)
      }
    }
    if(sum(sumPredData)>0){
      par(las=1, lend=1)
      plot(sumPredData, type="h", lwd=5, xaxt="n", xlab="", ylab="Est. realised proportion predators diet")
      mtext(thisPrey,side=3, adj=0, font=2)
      par(las=2)
      axis(at=seq(1,length(sumPredData)), labels =predLabels, side=1)
    }
  }
}
dev.off()



