#read in the catch.ts files for each box, join them together into an array, then get prop by box within each group for the last 5 years

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

catchPath<-paste(this_path,"inputs\\catch_history\\catchts\\",sep="")

dynBoxes<-seq(1,24) #in .ts files, there is boundary box (called boundary) then the dynamic boxes start at 1,
#the boxes 25-29 are boundaries too
nboxes<-length(dynBoxes)
catchYears<-seq(1900,2014); ny<-length(catchYears)
nts<-ny*12

#read in groups DF
groupsDF<-read.csv(paste(this_path,"CRAM_groups.csv",sep=""))
catchGroupsDF<-groupsDF[groupsDF$IsFished==1,]; ncg<-dim(catchGroupsDF)[1]

#set up arry
allCatches_array<-array(NA,dim=c(nboxes,ncg,nts))

timeSteps<-rep(0,nts)

for(b in 1:nboxes){
  thisTSfile<-paste(catchPath,"catch",b,".ts",sep="")
  thisTempLines<-readLines(thisTSfile)
  thisTSlines<-thisTempLines[grep("#", thisTempLines, invert = TRUE)]
  for(t in 1:nts){
    allCatches_array[b,,t]<-get_first_number(thisTSlines[t],n="all")[-1]
    timeSteps[t]<-get_first_number(thisTSlines[t],n=1)
  }
}
#get proportions over just the last 5 years - this will be the last 5*12 timesteps
startYear<-2010; start_ts<-(2010-1900)*12; end_ts<-nts
last5yearsCatchArray<-allCatches_array[,,start_ts:nts]

#for each group, get average prop in each box
aveBoxProp<-array(NA,dim=c(nboxes,ncg))

for(g in 1:ncg){
  thisData<-last5yearsCatchArray[,g,]
  temp<-apply(thisData,1,mean,na.rm=TRUE)
  aveBoxProp[,g]<-temp/sum(temp,na.rm=TRUE)
}

#write this out so can use it to create scenario forcing files
write.csv(aveBoxProp,paste(DIR$'Tables',"CatchHist_last5YearPropByBox.csv",sep=""),row.names = FALSE)

## write out hist catches, just by species and timestep
yearIndex<-ceiling(seq(1,nts)/12)
temp<-apply(allCatches_array, c(2,3), sum, na.rm=TRUE)
temp2<-melt(temp); 

temp2[,2]<-ceiling(temp2[,2]/12)
catchBySpeciesYear<-tapply(temp2$value,temp2[,2], sum, na.rm=TRUE)
