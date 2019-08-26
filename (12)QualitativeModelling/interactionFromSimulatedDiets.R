##creates interaction matrix based on simulated diets from 1865-2015 model run (no fishing)
##based on some min cut-off and selected years

minCutoff<-0.01 #units: mg N eaten in one year
yearsOfInterest<-seq(1900,2015)

year0<-1865; #this is when the model starts - includes 35 year burn-in that takes it up to the real start of 1900

yoiTimeSteps<-yearsOfInterest-year0+1 #to index years of interest timesteps

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
inputsPath<-paste(this_path,"base\\QualMod\\",sep="")

## first plot intended diets
groupsDF<-read.csv(paste(inputsPath,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
predators<-as.character(groupsDF$Code[groupsDF$IsPredator==1]); npreds<-length(predators)


# getAgeText<-function(x){
#   thisText<-"Juv"
#   if(is.na(x)){
#     thisText<-""
#   }else{
#     if(x==2){thisText="Adult"}
#   }
#   return(thisText)
# }

#####################################
realisedDiets<-read.csv(paste(inputsPath,"outputDietCheck.txt", sep=""), sep=" ")
realisedDiets$year_ts<-realisedDiets$Time/365
yearIndex<-realisedDiets$year_ts %in% yoiTimeSteps
yearRelDiets<-realisedDiets[yearIndex, ]
allPrey<-colnames(realisedDiets)[c(6:(dim(realisedDiets)[2]-1))]; nAllPrey<-length(allPrey)

realisedDietsByPrey<-array(0, dim=c(npreds, nAllPrey)); colnames(realisedDietsByPrey)<-allPrey; rownames(realisedDietsByPrey)<-predators
for(p in 1:npreds){
  thisPred<-predators[p]
  temp<-yearRelDiets[yearRelDiets$Predator==thisPred,]
  test<-temp[,6:(dim(temp)[2]-1)]
  posPreyIndex<-colSums(test)>minCutoff
  thisPrey<-allPrey[posPreyIndex]
  ## allprey
  realisedDietsByPrey[p,match(thisPrey, allPrey)]<-1
}

## write the summary table as .csv
yearsText<-paste(unique(c(min(yearsOfInterest),max(yearsOfInterest))),collapse="_") #will give year range if more than one year
cutoffText<-gsub("\\.", "", signif(minCutoff,2))
## write the full diet summary
write.csv(realisedDietsByPrey, paste(inputsPath, "interactionFromSimDiets", yearsText,"_", cutoffText,".csv", sep=""), row.names = predators)

#check how many prey for each predator for given cut-off
test<-apply(realisedDietsByPrey, 1, sum)
test
# and how many predator groups eating each prey group
test<-apply(realisedDietsByPrey, 2, sum)
test

