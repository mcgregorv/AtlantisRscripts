#under initial conditions, using prey available, how much food does each predator have in the same place and how much does it vary?

load(paste(DIR$"Data","\\eating\\interaction_spatial_array",sep="")); ##to bring in interaction_spatial_array. 
# created in ExperimentingWithInformingDiets_arrayApproach.R
ppGroups<-as.character(read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ppGroups.csv",sep=""))[,1]); npg<-length(ppGroups)

basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
plotPath<-paste(basePath, "Figures\\Base\\DIET_check\\",sep="")

#read in original diet csv and plot diets
origDiets<-read.csv(paste(basePath,"inputs\\supporting\\CRAM_DietsForR.csv", sep=""))
predCodes<-unique(origDiets$Code); np<-length(predCodes)

plotPath<-paste(basePath, "Figures\\Base\\DIET_check\\",this_out,sep="")

diet_check<-read.csv(paste(outPath,"outputDietCheck.txt",sep=""),skip=0,sep=" ",header=TRUE)
predators<-unique(diet_check$Predator); np<-length(predators)

spatialInteractionByGroup<-apply(interaction_spatial_array,c(1,2), sum)
foodAvailByPredator<-rep(NA,npg)

for(p in 1:npg){
  thisPredVar<-ppGroups[p]
  thisPred<-gsub("juv|ad","",thisPredVar); thisPredAge<-gsub(thisPred,"",thisPredVar)
  if(thisPred %in% predators){
    thisPredLongAge<-"All"
    if(thisPredAge=="ad"){
      thisPredLongAge<-"Adult"
    }else if(thisPredAge=="juv"){
      thisPredLongAge<-"Juvenile"
    }
    temp<-origDiets[origDiets$Code==thisPred & origDiets$Predator.age==thisPredLongAge,]
    thisAvail<-temp[,colnames(temp) %in% groupsDF$Code]
    thisFoodAvail<-0
    if(thisPredAge==""){
      preyAge<-c("juv", "ad")[preyAgeIndex]
      thisPreys<-names(thisAvail)[thisAvail>0]; thisPreyAvails<-thisAvail[thisAvail>0]
      npreys<-length(thisPreys)
      for(f in 1:npreys){
        thisPrey<-thisPreys[f]
        thisAgePrey<-paste(thisPrey,preyAge, sep="")
        if(!(thisAgePrey %in% ppGroups)){thisAgePrey<-thisPrey}
        rowIndex<-seq(1,npg)[ppGroups==thisPredVar]; colIndex<-seq(1,npg)[ppGroups==thisAgePrey]
        thisSpaceAvail<-spatialInteractionByGroup[rowIndex, colIndex]*thisPreyAvails[f]
        thisFoodAvail<-thisFoodAvail+thisSpaceAvail
      }  
    } else{
      for(preyAgeIndex in 1:2){
        preyAge<-c("juv", "ad")[preyAgeIndex]
        thisPreys<-colnames(thisAvail)[thisAvail[preyAgeIndex,]>0]; thisPreyAvails<-thisAvail[preyAgeIndex,thisAvail[preyAgeIndex,]>0]
        npreys<-length(thisPreys)
        if(npreys>0){
          for(f in 1:npreys){
            thisPrey<-thisPreys[f]
            thisAgePrey<-paste(thisPrey,preyAge, sep="")
            if(!(thisAgePrey %in% ppGroups)){thisAgePrey<-thisPrey}
            rowIndex<-seq(1,npg)[ppGroups==thisPredVar]; colIndex<-seq(1,npg)[ppGroups==thisAgePrey]
            thisSpaceAvail<-spatialInteractionByGroup[rowIndex, colIndex]*thisPreyAvails[f]
            thisFoodAvail<-thisFoodAvail+thisSpaceAvail
          }  
        }
      }
    }
    foodAvailByPredator[p]<-as.double(thisFoodAvail)
  }
}

plotIndex<-!is.na(foodAvailByPredator)

plotIndex<-foodAvailByPredator>10 & !is.na(foodAvailByPredator)

plotIndex<-foodAvailByPredator>7 & foodAvailByPredator<=10 & !is.na(foodAvailByPredator)

plotIndex<-foodAvailByPredator>6 & foodAvailByPredator<=7 & !is.na(foodAvailByPredator)

plotIndex<-foodAvailByPredator>5 & foodAvailByPredator<=6 & !is.na(foodAvailByPredator)

par(lend=1)
plot(foodAvailByPredator[plotIndex], type="h", lwd=5, xaxt="n", xlab="", ylab="Saptial availbility")
par(las=2)
axis(at=seq(1,length(foodAvailByPredator[plotIndex])), labels=ppGroups[plotIndex], side=1)

##want spatial overlap from tracers
