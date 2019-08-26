##plots diets of each group by summary prey groupings (in preyGroups.csv)

source(paste(DIR$'General functions',"get_interaction_trophicLevel.R",sep=""))
source(paste(DIR$'General functions',"get_first_number.R",sep=""))
source(paste(DIR$'General functions',"myRounding.R",sep=""))
source(paste(DIR$'General functions',"getSNRN.R",sep=""))
source(paste(DIR$'General functions',"nonZeroMean.R",sep=""))
source(paste(DIR$'General functions',"getICNumbers.R",sep=""))

source(paste(DIR$'General functions',"get_interaction_spatial_bySpace.R",sep=""))
source(paste(DIR$'General functions',"get_interaction_gapeSize_byCohort.R",sep=""))

this_run<-"Base"
this_out<-"FISHXXX_mLFB6"
dietTimeSteps<-seq(1,8) #set this to summarise diet over these timesteps - may want to use all, may want to use last 10 or 50 years for example

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
inputsPath<-paste(this_path,"inputs\\",sep="")
# plotPath<-paste(this_path,"Figures\\",this_run,"\\",sep="")
plotPath<-paste(DIR$'Figures',"Mortality\\test",this_out, sep="")
outPath<-paste(this_path,"base\\output",this_out,"\\",sep="")

preyGroupsDF<-read.csv(paste(this_path,"inputs\\supporting\\preyGroups.csv",sep=""))
preyGroups<-sort(unique(preyGroupsDF$PreyGroup)); nPreyGroups<-length(preyGroups)

## first plot intended diets
groupsDF<-read.csv(paste(this_path,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

predators<-as.character(groupsDF$Code[groupsDF$IsPredator==1]); npreds<-length(predators)

dietsOriginal<-read.csv(paste(inputsPath,"supporting\\CRAM_Diets.csv",sep=""))
diets<-read.csv(paste(inputsPath,"biol_prm\\pPREY\\CRAM_biol_pPREY_original.csv", sep=""))
preyCodesOrig<-colnames(dietsOriginal)[7:(7+ng-1)]; preyCodeIndex<-match(groupsDF$Code, preyCodesOrig)
##overwrite diets with values from diestOriginal. Diets has values ready for pPREY (proportions of prey available to pred)
# but want proportions of diets here
for(j in 2:dim(diets)[2]){
  thisVar<-colnames(diets)[j]; preyAge<-get_first_number(thisVar, n=1); predAge<-get_first_number(thisVar, n=2); thisCode<-gsub(paste("pPREY", predAge,preyAge,sep="|"), "", thisVar)
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
  diets[1:ng,j]<-as.double(newColumn[preyCodeIndex])
}

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
  nonZeroPreys<-diets[index,1]
  thisPreyGroups<-preyGroupsDF$PreyGroup[match(nonZeroPreys,preyGroupsDF$Code)]
  temp<-data.frame(cbind(nonZeroDiet,"preyGroup"=as.character(thisPreyGroups)))
  dietByPreyGroup<-tapply(as.double(as.character(temp$nonZeroDiet)),temp$preyGroup, sum, na.rm=TRUE)
  propDiet<-dietByPreyGroup/sum(dietByPreyGroup, na.rm=TRUE); preyGroupIndex<-preyGroups %in% names(propDiet)
  dietsByPreyGroup[p,preyGroupIndex]<-propDiet
}
dietsByPreyGroup<-data.frame(dietsByPreyGroup); dietsByPreyGroup$Predator<-predators

preyGroupColours<-colorRampPalette(colors=c("black", myGrey,"midnightblue",myBlue,myLightBlue,myLightAqua, myAqua,myDarkGreen,myGreen,myGold,myOrange,"red","brown",myRed))(nPreyGroups)


#####################################
realisedDiets<-read.csv(paste(outPath,"outputDietCheck.txt", sep=""), sep=" ")
realisedDiets$year_ts<-realisedDiets$Time/365
yearIndex<-realisedDiets$year_ts %in% dietTimeSteps
yearRelDiets<-realisedDiets[yearIndex, ]

summaryRealisedDiets<-array(NA,dim=c(npreds, nPreyGroups))

### SEPERATE THESE INTO ADULT AND JUVENILE, IT MAKES A DIFFERENCE FOR SOME (SUCH AS ELP)
p<-grep("EIS",predators)

# 
# pdf(paste(plotPath,"TestingSummary_realised.pdf", sep=""),height=10,width=7)
# par(mar=c(7,4,1,1), las=1, mfrow=c(5,2))
# for(p in 1:npreds){
  thisPred<-predators[p]; cat(thisPred)
  temp<-realisedDiets[realisedDiets$Predator==thisPred,]
  test<-apply(temp[,c(6:(dim(temp)[2]-1))],2,mean)
  thisPreydf<-data.frame(cbind("preyCode"=names(test), "preyProportion"=as.double(test)));
  thisPreydf$preyGroup<-preyGroupsDF$PreyGroup[match(thisPreydf$preyCode,preyGroupsDF$Code)]
  # thisPreydf$preyProportion[thisPreydf$preyCode %in% c("BO","BD")]<-0
  thisPreySummary<-tapply(as.double(as.character(thisPreydf$preyProportion)),thisPreydf$preyGroup,sum)
  par(las=1)
  plot(thisPreySummary,type="h",xlab="", ylim=c(0,1))
  mtext(thisPred,side=3)
  points(as.double(dietsByPreyGroup[p,]),pch=8,col="red")
  par(las=2)
  axis(at=seq(1,nPreyGroups), labels=names(thisPreySummary), side=1)
# 
# }
# dev.off()

# preyGroups
# Algae            Bacteria         Bird             Cetacea          Coelenterate     Crustacean       Detritus         Echinoderm       Elasmobranch    
# [10] Mammal           Microzooplankton Mollusc          Phytoplankton    Polychaete       Teleost          Tunicate        


thisPreydf[thisPreydf$preyGroup %in% c("Echinoderm"),]
thisPreydf[thisPreydf$preyGroup %in% c("Polychaete"),]

# which groups eat CEP and ASQ?
thisCode<-"HOK"
index<-realisedDiets[,grep(thisCode,colnames(realisedDiets))]>0
test<-realisedDiets[index,c(1:3,grep(thisCode,colnames(realisedDiets)))]
eatenByTime<-tapply(test[,c(thisCode)], test$Time, sum, na.rm=TRUE)
# plot(eatenByTime)
eatenByPredator<-tapply(test[,c(thisCode)], test$Predator, sum, na.rm=TRUE)
eatenByPredator<-eatenByPredator[!is.na(eatenByPredator)]
par(las=2)
plot(eatenByPredator,type="h",xaxt="n"); axis(at=seq(1,length(eatenByPredator)), labels = names(eatenByPredator), side=1)  
par(las=1)
mtext(thisCode,side=3)

eatenByCohort<-tapply(test[,c(thisCode)], test$Cohort, sum, na.rm=TRUE)

# preyGroups
# [1] Algae            Bacteria         Bird             Cetacea          Coelenterate     Crustacean       Detritus         Echinoderm       Elasmobranch    
# [10] Mammal           Microzooplankton Mollusc          Phytoplankton    Polychaete       Teleost          Tunicate        

realisedDietsByPreyGroup<-array(NA, dim=c(npreds, nPreyGroups)); colnames(realisedDietsByPreyGroup)<-preyGroups; rownames(dietsByPreyGroup)<-predators
for(p in 1:npreds){
  thisPred<-predators[p]
  temp<-realisedDiets[realisedDiets$Predator==thisPred,]
  test<-temp[,6:(dim(temp)[2]-1)]
  posPreyIndex<-colSums(test)>0
  thisDiet<-cbind("year"=temp$year_ts, test[,posPreyIndex])
  thisPreys<-colnames(test)[posPreyIndex]; nprey<-length(thisPreys)
  #get proportions of each prey
  # totalPrey<-sum(thisDiet[,-1]); preySums<-apply(thisDiet[,-1],2,sum,na.rm=TRUE); preyProps<-preySums/totalPrey
  # 
  # timeSums<-array(NA, dim=c(length(unique(thisDiet$year)),(dim(thisDiet)[2]-1)))
  # for(i in 1:length(thisPreys)){
  #   xx<-thisDiet[,c(1,(i+1))]; xxx<-tapply(xx[,2], xx[,1], sum)
  #   timeSums[,i]<-xxx
  # }
  thisPreyGroupIndex<-match(thisPreys,preyGroupsDF$Code); thisPreyGroups<-preyGroupsDF$PreyGroup[thisPreyGroupIndex]
  tempDF<-data.frame(cbind(preyProps,"preyGroups"=as.character(thisPreyGroups))); preyByPreyGroup<-tapply(as.double(as.character(tempDF$preyProps)),tempDF$preyGroups,sum, na.rm=TRUE)
  realisedDietsByPreyGroup[p, match(names(preyByPreyGroup),preyGroups)]<-preyByPreyGroup
}

rDietsByPreyGroup<-data.frame(realisedDietsByPreyGroup); rDietsByPreyGroup$Predator<-predators
#skip cetacea to keep the colours lined up - not noticable amount eaten of them
rDietsByPreyGroup<-rDietsByPreyGroup[,colnames(rDietsByPreyGroup)!="Cetacea"]

tempplotDF<-melt(rDietsByPreyGroup, id.var="Predator")
index<-tempplotDF$variable=="Mammal"; plotDF<-tempplotDF[!index,]
bp<-ggplot(data = plotDF, aes(x = Predator, fill = variable, y = value)) + 
  geom_bar(stat = 'identity')


pdf(paste(plotPath,"DietSummary_realised.pdf", sep=""),height=7,width=7)
par(mar=c(4,4,1,1))
bp + coord_flip()  + scale_fill_manual(values=preyGroupColours) + labs(y="Proportion of diet") + theme_igray() + 
  theme(axis.text=element_text(size=12, angle=0),axis.title=element_text(size=12)) + guides(fill=guide_legend(title="")) 
dev.off()

## create a table for the groupings

