##plots diets of each group by summary prey groupings (in preyGroups.csv)

source(paste(DIR$'General functions',"get_interaction_trophicLevel.R",sep=""))
source(paste(DIR$'General functions',"get_first_number.R",sep=""))
source(paste(DIR$'General functions',"myRounding.R",sep=""))
source(paste(DIR$'General functions',"getSNRN.R",sep=""))
source(paste(DIR$'General functions',"nonZeroMean.R",sep=""))
source(paste(DIR$'General functions',"getICNumbers.R",sep=""))

source(paste(DIR$'General functions',"get_interaction_spatial_bySpace.R",sep=""))
source(paste(DIR$'General functions',"get_interaction_gapeSize_byCohort.R",sep=""))

this_out<-"TestsSCA1"; runFolder<-"TBGB_JP2";
# this_out<-"output"; runFolder<-"TBGBReportBase";
# this_out<-"outputSG"; runFolder<-"TBGBReportBase"; 
# this_out<-"MyRun_Fish1899_better1_codeupdate_rewriteDmatrix1"; runFolder<-"TBGBFish"; 
thisDesc <- paste(runFolder, this_out,sep="")

basePath<-  paste(DIR$'Base', "TBGB\\",runFolder,"\\",sep="")
outPath<-paste(basePath,"\\output",this_out,"\\",sep="")
plotPath<-paste(basePath,"..\\Figures\\Testing\\DIET\\",thisDesc, sep="")

groupsDF<-read.csv(paste(basePath,"\\TBGB_Groups.csv",sep="")); ng<-dim(groupsDF)[1]

  year0<-1899; #this is when the model starts - includes 35 year burn-in that takes it up to the real start of 1900
  
  dietYears<-seq(1880,1960) #just replace with 1970 if only want one year, or give a range of years - will give average over them
  dietTimeSteps<-dietYears-year0 #to index dietYears
  
  preyGroupsDF<-read.csv(paste(basePath,"..\\TBGB_PreyGroups.csv",sep=""))
  preyGroups<-sort(unique(preyGroupsDF$PreyGroup)); nPreyGroups<-length(preyGroups)
  
   ## first plot intended diets
  predators<-as.character(groupsDF$Code[groupsDF$IsPredator==1]); npreds<-length(predators)
  
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
  
  preyGroupColours<-colorRampPalette(colors=c("black", myGrey,"midnightblue",myBlue,myLightBlue,myLightAqua, myAqua,myDarkGreen,myGreen,myGold,myOrange,"red","brown",myRed))(nPreyGroups)
    #####################################
  realisedDiets<-read.csv(paste(outPath,"\\outputDietCheck.txt", sep=""), sep=" ")
  realisedDiets$year_ts<-realisedDiets$Time/365
  yearIndex<-realisedDiets$year_ts %in% dietTimeSteps
  yearRelDiets<-realisedDiets[yearIndex, ]
  
  ## check ZL preds
  thisPrey <- "SCA"
  index<-realisedDiets[,c(thisPrey)]>0;
  test<-realisedDiets[index,c(1:4,grep(thisPrey,colnames(realisedDiets)))]
  timeIndex<-test$Time<500
  yy <- tapply(test[,c(thisPrey)], test$Predator, sum, na.rm=TRUE)
  yy <- yy[!is.na(yy)]
  par(mfrow=c(1,1), lend=1, las=2)
  plot(as.double(yy[yy>0]), type="h", lwd=5, xaxt="n"); axis(at=seq(1,length(yy[yy>0])), labels=names(yy[yy>0]), side=1)

  byTime <- tapply(test[,c(thisPrey)], list(test$Predator, test$Time), sum, na.rm=TRUE)
  par(mfrow=c(3,2), las=1)
  for(i in 1:length(yy)){
    plot(byTime[c(names(yy)[i]),], type="l"); mtext(names(yy)[i], side=3)
  }

   # test<-realisedDiets[index,]
   
   predByTime <- tapply(test[,c(thisPrey)], list(test[,c("Time")], test[,c("Predator")]), nonZeroMean)
  
  summaryRealisedDiets<-array(NA,dim=c(npreds, nPreyGroups))
  
  par(las=2, mfrow=c(1,1))
  thisPred <- "BC"
  test<-realisedDiets[grep(thisPred,realisedDiets$Predator ),]; 
  xx <- apply(test[,6:(dim(test)[2]-1)],2,mean, na.rm=TRUE)
  # par(las=2, mfrow=c(1,1))
  plot(as.double(xx[xx>0]), type="h", lwd=5, xaxt="n"); axis(at=seq(1,length(xx[xx>0])), labels=names(xx[xx>0]), side=1)
  par(las=1)
  mtext(thisPred, side=3)
  
  par(mfrow=c(2,1), las=2, lend=1)
  testAd <- test[test$Cohort>0,]
  xxAd <- apply(testAd[,6:(dim(testAd)[2]-1)],2,mean, na.rm=TRUE)
  # par(las=2, mfrow=c(1,1))
  plot(as.double(xxAd[xxAd>0]), type="h", lwd=5, xaxt="n", xlab="Prey", ylab="Proportion of diet"); axis(at=seq(1,length(xxAd[xxAd>0])), labels=names(xxAd[xxAd>0]), side=1)
  par(las=1)
  mtext(paste(thisPred," adults",sep=""), side=3)
  
  testJuv <- test[test$Cohort==0,]
  xxJuv <- apply(testJuv[,6:(dim(testJuv)[2]-1)],2,mean, na.rm=TRUE)
  # par(las=2, mfrow=c(1,1))
  par(las=2)
  plot(as.double(xxJuv[xxJuv>0]), type="h", lwd=5, xaxt="n", xlab="Prey", ylab="Proportion of diet"); axis(at=seq(1,length(xxJuv[xxJuv>0])), labels=names(xxJuv[xxJuv>0]), side=1)
  par(las=1)
  mtext(paste(thisPred," juv",sep=""), side=3)
  
  
  ### SEPERATE THESE INTO ADULT AND JUVENILE, IT MAKES A DIFFERENCE FOR SOME (SUCH AS ELP)
  pred2plot<-sort(names(yy[yy>0]), decreasing = TRUE)
  par(mfrow=c(3,2))
  for(p in 1:length(pred2plot)){
    thisPred <- pred2plot[p]
    test<-realisedDiets[grep(thisPred,realisedDiets$Predator ),];
    xx <- apply(test[,6:(dim(test)[2]-1)],2,max, na.rm=TRUE)
    par(las=2)
    plot(as.double(xx[xx>0]), type="h", lwd=5, xaxt="n"); axis(at=seq(1,length(xx[xx>0])), labels=names(xx[xx>0]), side=1)
    par(las=1)
    mtext(thisPred, side=3)
  }
  
  
  pdf(paste(plotPath,"TestingSummary.pdf", sep=""),height=10,width=7)
  par(mar=c(7,4,1,1), las=1, mfrow=c(5,2))
  for(p in 1:npreds){
    thisPred<-predators[p]; cat(thisPred)
    temp<-realisedDiets[realisedDiets$Predator==thisPred,]
    test<-apply(temp[,c(6:(dim(temp)[2]-1))],2,mean)
    thisPreydf<-data.frame(cbind("preyCode"=names(test), "preyProportion"=as.double(test))); 
    thisPreydf$preyGroup<-preyGroupsDF$PreyGroup[match(thisPreydf$preyCode,preyGroupsDF$Code)]
    # thisPreydf$preyProportion[thisPreydf$preyCode %in% c("BO","BD")]<-0
    thisPreySummary<-tapply(as.double(as.character(thisPreydf$preyProportion)),thisPreydf$preyGroup,sum)
    par(las=1)
    plot(thisPreySummary,type="h",xlab="", ylim=c(0,1), xaxt="n")
    mtext(thisPred,side=3)
    par(las=2)
    axis(at=seq(1,nPreyGroups), labels=names(thisPreySummary), side=1)
    
  } 
  dev.off()
  
  # preyGroups
  # Algae            Bacteria         Bird             Cetacea          Coelenterate     Crustacean       Detritus         Echinoderm       Elasmobranch    
  # [10] Mammal           Microzooplankton Mollusc          Phytoplankton    Polychaete       Teleost          Tunicate        
  
  
  thisPreydf[thisPreydf$preyGroup %in% c("Echinoderm"),]
  thisPreydf[thisPreydf$preyGroup %in% c("Polychaete"),]
  thisPreydf[thisPreydf$preyGroup %in% c("Coelenterate"),]
  
  
  # preyGroups
  # [1] Algae            Bacteria         Bird             Cetacea          Coelenterate     Crustacean       Detritus         Echinoderm       Elasmobranch    
  # [10] Mammal           Microzooplankton Mollusc          Phytoplankton    Polychaete       Teleost          Tunicate        
  allPrey<-colnames(realisedDiets)[c(6:(dim(realisedDiets)[2]-1))]; nAllPrey<-length(allPrey)
  
  realisedDietsByPreyGroup<-array(NA, dim=c(npreds, nPreyGroups)); colnames(realisedDietsByPreyGroup)<-preyGroups; rownames(realisedDietsByPreyGroup)<-predators
  realisedDietsByPrey<-array(NA, dim=c(npreds, nAllPrey)); colnames(realisedDietsByPrey)<-allPrey; rownames(realisedDietsByPrey)<-predators
  for(p in 1:npreds){
    thisPred<-predators[p]
    temp<-realisedDiets[realisedDiets$Predator==thisPred,]
    test<-temp[,6:(dim(temp)[2]-1)]
    posPreyIndex<-colSums(test, na.rm=TRUE)>0
    xx <- test[,posPreyIndex]
    thisDiet<-cbind("year"=temp$year_ts, xx)
    thisPreys<-colnames(test)[posPreyIndex]; nprey<-length(thisPreys)
    #get proportions of each prey
    if(dim(thisDiet)[2]>1){
      totalPrey<-sum(thisDiet[,-1], na.rm=TRUE); 
      if(length(thisPreys)==1){
        preySums<-sum(thisDiet[,2], na.rm=TRUE)
      }else{
        preySums<-apply(thisDiet[,-1],2,sum,na.rm=TRUE); 
      }
      
    } else{
      preySums<-apply(thisDiet[,-1],2,sum,na.rm=TRUE); 
    }
    preyProps<-preySums/totalPrey
    thisPreyGroupIndex<-match(thisPreys,preyGroupsDF$Code); thisPreyGroups<-preyGroupsDF$PreyGroup[thisPreyGroupIndex]
    tempDF<-data.frame(cbind(preyProps,"preyGroups"=as.character(thisPreyGroups))); preyByPreyGroup<-tapply(as.double(as.character(tempDF$preyProps)),tempDF$preyGroups,sum, na.rm=TRUE)
    realisedDietsByPreyGroup[p, match(names(preyByPreyGroup),preyGroups)]<-preyByPreyGroup
    ## allprey
    realisedDietsByPrey[p,match(names(preyProps), allPrey)]<-preyProps
  }
  
  
  rDietsByPreyGroup<-data.frame(realisedDietsByPreyGroup); 
  rDietsByPreyGroup$Predator<-groupsDF$Name[match(rownames(rDietsByPreyGroup), groupsDF$Code)]
  #skip cetacea to keep the colours lined up - not noticable amount eaten of them
  rDietsByPreyGroup<-rDietsByPreyGroup[,colnames(rDietsByPreyGroup)!="Cetacea"]
  
  tempplotDF<-melt(rDietsByPreyGroup, id.var="Predator")
  index<-tempplotDF$variable=="Mammal"; plotDF<-tempplotDF[!index,]
  bp<-ggplot(data = plotDF, aes(x = Predator, fill = variable, y = value)) + 
    geom_bar(stat = 'identity')
  
  
  pdf(paste(plotPath,"DietSummary.pdf", sep=""),height=7,width=8)
  par(mar=c(4,3.5,1,1))
  bp + coord_flip()  + scale_fill_manual(values=preyGroupColours) + labs(y="Proportion of diet", x="") + theme_igray() + 
    theme(axis.text=element_text(size=12, angle=0),axis.title=element_text(size=12)) + guides(fill=guide_legend(title="")) 
  dev.off()
  
  ## do intended diets
  ## colnames will have predators, with pPREY, prey age, PredCode, pred age
  predCode<-groupsDF$Code[groupsDF$IsPredator==1]; npreds <- length(predCode)
  ageIndex <- groupsDF$NumCohorts[groupsDF$IsPredator==1]>1
  agedVrs <- paste("pPREY",c(1,1,2,2), sort(rep(predCode[ageIndex],4)),c(1,2,1,2), sep="")
  noageVrs <- paste("pPREY",predCode[!ageIndex], sep="")
  
  predVars<-sort(c(agedVrs, noageVrs)); npredVrs<-length(predVars)
  nPrey <- dim(groupsDF)[1]
  
  dietDF <- data.frame(matrix(NA, ncol=npredVrs, nrow=nPrey))
  colnames(dietDF)<- predVars; rownames(dietDF)<- as.character(groupsDF$Code)
  dietsOriginal <- read.csv(paste(basePath, "..\\TBBOriginalDiet.csv",sep=""))
  
  preyCodesOrig<-colnames(dietsOriginal)[7:(7+ng-1)]; preyCodeIndex<-match(groupsDF$Code, preyCodesOrig)
  
  
  for(j in 1:npredVrs){
    thisPredVar <- predVars[j]
    preyAge<-get_first_number(thisPredVar, n=1); predAge<-get_first_number(thisPredVar, n=2); thisCode<-gsub(paste("pPREY", predAge,preyAge,sep="|"), "", thisPredVar)
    xx<-grep(thisCode, dietsOriginal$Code); 
    if(!is.na(preyAge) & !is.na(predAge)){
      xxx<-dietsOriginal[xx:(xx+3),]
      predAgeText<-c("Juvenile", "Adult")[predAge]; preyAgeText<-c("Juvenile", "Adult")[preyAge]
    } else{
      xxx<-dietsOriginal[xx,]
      predAgeText<-"All"; preyAgeText<-"All"
    }  
    index<-xxx$Predator.age==predAgeText & xxx$Prey.age==preyAgeText
    thisDiets<-xxx[index,]
    newColumn<-thisDiets[7:dim(dietsOriginal)[2]]
    dietDF[1:ng,j]<-as.double(newColumn[preyCodeIndex])
  
  }
  preys<-dietDF[,1]; npreys<-length(preys)
  diets<-dietDF
  dietsByPreyGroup<-array(NA, dim=c(npreds, nPreyGroups)); colnames(dietsByPreyGroup)<-preyGroups; rownames(dietsByPreyGroup)<-predators
  for(p in 1:npreds){
    thisPred<-predators[p]
    if(thisPred=="BO"){
      temp<-diets[,colnames(diets)==paste("pPREY",thisPred,sep="")]; thisDiet<-temp
    }else{
      temp<-diets[,grep(thisPred,colnames(diets))]; thisDiet<-temp
    }
    if(length(dim(temp))==2){
      if(dim(temp)[2]>1){
        thisDiet<-apply(temp,1,sum,na.rm=TRUE)
      }
    }
    index<-thisDiet>0
    nonZeroDiet<-thisDiet[index]
    nonZeroPreys<-rownames(diets)[index]
    thisPreyGroups<-preyGroupsDF$PreyGroup[match(nonZeroPreys,preyGroupsDF$Code)]
    temp<-data.frame(cbind(nonZeroDiet,"preyGroup"=as.character(thisPreyGroups)))
    dietByPreyGroup<-tapply(as.double(as.character(temp$nonZeroDiet)),temp$preyGroup, sum, na.rm=TRUE)
    propDiet<-dietByPreyGroup/sum(dietByPreyGroup, na.rm=TRUE); preyGroupIndex<-preyGroups %in% names(propDiet)
    dietsByPreyGroup[p,preyGroupIndex]<-propDiet
  }
  predatorNames<-gsub("_", " ", groupsDF$Name[match(predators, groupsDF$Code)])
  dietsByPreyGroup<-data.frame(dietsByPreyGroup); dietsByPreyGroup$Predator<-predatorNames
  
  plotDF<-melt(dietsByPreyGroup, id.var="Predator")
  bp<-ggplot(data = plotDF, aes(x = Predator, fill = variable, y = value)) + 
    geom_bar(stat = 'identity')
  preyGroupColours<-colorRampPalette(colors=c("black",myGrey,myPurple, myLightBlue,myBlue,myLightAqua, myAqua,myDarkGreen,myGreen,myGold,myOrange,"red",myRed,myPurple))(nPreyGroups)
  
  preyGroupColours<-colorRampPalette(colors=c("black", myGrey,"midnightblue",myBlue,myLightBlue,myLightAqua, myAqua,myDarkGreen,myGreen,myGold,myOrange,"red","brown",myRed))(nPreyGroups)
  
  # 
  # pdf(paste(plotPath,"DietSummary_intended.pdf", sep=""),height=7,width=9)
  # par(mar=c(4,3.5,1,1))
  bp + coord_flip()  + scale_fill_manual(values=preyGroupColours) + labs(y="Proportion of diet", x="") + theme_igray() +
    theme(axis.text=element_text(size=12, angle=0),axis.title=element_text(size=12)) + guides(fill=guide_legend(title=""))
  # dev.off()
  
  
  
  pdf(paste(plotPath,this_out,"TestingSummary_dietRealisedAndIntended.pdf", sep=""),height=10,width=7)
  par(mar=c(7,4,2,1), las=1, mfrow=c(4,2))
  for(p in 1:npreds){
    thisPred<-predators[p]; cat(thisPred)
    temp<-realisedDiets[realisedDiets$Predator==thisPred,]
    test<-apply(temp[,c(6:(dim(temp)[2]-1))],2,mean)
    thisPreydf<-data.frame(cbind("preyCode"=names(test), "preyProportion"=as.double(test))); 
    thisPreydf$preyGroup<-preyGroupsDF$PreyGroup[match(thisPreydf$preyCode,preyGroupsDF$Code)]
    # thisPreydf$preyProportion[thisPreydf$preyCode %in% c("BO","BD")]<-0
    thisPreySummary<-tapply(as.double(as.character(thisPreydf$preyProportion)),thisPreydf$preyGroup,sum)
    par(las=1)
    plot(thisPreySummary,type="n",xlab="", ylim=c(0,1), xaxt="n")
    points(x=seq(1,length(thisPreySummary))[thisPreySummary>0], y=thisPreySummary[thisPreySummary>0],pch=20,xlab="", ylim=c(0,1), xaxt="n", col="cornflowerblue", cex=1.5)
    mtext(thisPred,side=3,adj=0)
    par(las=2)
    axis(at=seq(1,nPreyGroups), labels=names(thisPreySummary), side=1)
    
    thisIntended <- dietsByPreyGroup[c(thisPred),colnames(dietsByPreyGroup)!="Predator"]
    points(as.double(thisIntended), pch=8, col="red")
    par(xpd=NA)
     legend(legend=c("Intended","realised"),pch=c(8,20), lwd=c(NA, NA), col=c("red", "cornflowerblue"), x="top", bty="n", ncol=2, inset=-0.2)
  } 
  dev.off()
 
  
  
  ## write the summary table as .csv
  yearsText<-paste(unique(c(min(dietYears),max(dietYears))),collapse="_") #will give year range if more than one year
  write.csv(realisedDietsByPreyGroup, paste(outPath, "\\realiseDietSummarySnapshot",yearsText,".csv", sep=""), row.names = predators)
  ## write the full diet summary
  write.csv(realisedDietsByPrey, paste(outPath, "\\realiseDietSnapshot", yearsText,".csv", sep=""), row.names = predators)
  
  
# }
# 
# ## who eats scallops
# test<-realisedDiets[,c("Time","Predator","Cohort","SCA")]; 
# index<-test$SCA>0
# scaPred <- test[index,]
# 
# meanByPred<-tapply(scaPred$SCA, scaPred$Predator, mean, na.rm=TRUE)
# 
# meanByPred[!is.na(meanByPred)]
# 
# test<-realisedDiets[,c("Time","Predator","Cohort","OYS")]; 
# index<-test$OYS>0
# scaPred <- test[index,]
# 
# meanByPred<-tapply(scaPred$OYS, scaPred$Predator, mean, na.rm=TRUE)
# 
# meanByPred[!is.na(meanByPred)]
# # meanByPred[!is.na(meanByPred)]
# # ELI          FLA          GUR          IVH           SB          TAR 
# # 1.933896e-01 3.275173e-01 1.182062e-02 5.741849e-05 6.213887e-02 1.055614e-01 
# # round(meanByPred[!is.na(meanByPred)],2)
# # ELI  FLA  GUR  IVH   SB  TAR 
# # 0.19 0.33 0.01 0.00 0.06 0.11 
# 
# maxByPred<-tapply(scaPred$SCA, scaPred$Predator, max, na.rm=TRUE)
# round(maxByPred[!is.na(maxByPred)],2)
# 
# ## fla cons as ts
# flaPred<-scaPred[scaPred$Predator=="FLA" & scaPred$Cohort==4,]
# plot(x=flaPred$Time/365, y=flaPred$SCA, type="l", ylab="Proportion of FLA age-class 4 diet", xlab="Year")
# 
# 
# ROMStemp<-nc_open(paste(basePath,"ROMS_inputs\\GoldenBay26_temp.nc", sep=""))
# ROMStempdata<-ncvar_get(ROMStemp,"temperature")

  # 
  # shQuote("C:/NIST08/AMDIS32/AMDIS_32.exe /S C:/Users/Ento/Documents/GCMS/test_cataglyphis_iberica/queens/CI23_Q_120828_01.CDF" , type = "cmd" )
  # 
  # 
  # atlantismain -i create_biol_input\TBGB_input.nc 0 -o output.nc -r TBGB_run_Fish.prm -f TBGB_force.prm -p TBGB_physics.prm -b TBGB_biol.prm -h TBGB_harvest_TS_setup.prm -s TBGB_Groups.csv -q TBGB_Fisheries.csv -d outputLongUnorderedGroups
  # 
  # 
  # shQuote("C:/NIST08/AMDIS32/AMDIS_32.exe /S C:/Users/Ento/Documents/GCMS/test_cataglyphis_iberica/queens/CI23_Q_120828_01.CDF" , type = "cmd" )