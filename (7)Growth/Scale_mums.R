#weights seem to be settling at around 0.5 for all. Will try multiplying mum's by 2, then if it works figger out why!
##uses csv created in Estimate_mum_byCohort.R
mumByGroupCohort<-read.csv(paste(basePath,"mumByCohort.csv",sep=""))
mumByGroupCohort[,1:10]<-2*mumByGroupCohort[,1:10]
#groups csv file
groupsFile<-paste(basePath,"CRAM_groups.csv",sep=""); groupsDF<-read.csv(groupsFile); ng<-dim(groupsDF)[1]
#biolprm file
biolFile<-paste(basePath,"CRAM_base_biol.prm",sep=""); biolCopyFile<-paste(basePath,"CRAM_base_biol_copy.prm",sep="")
#make a copy in case botch it up
file.copy(biolFile,biolCopyFile)
##copy it back if need to :-)
# file.copy(biolCopyFile,biolFile,overwrite=TRUE)

biolLines<-readLines(biolFile)

for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
  thisVar<-paste("mum_",thisCode,sep="")
  thisMums<-as.double(mumByGroupCohort[mumByGroupCohort$Code==as.character(thisCode),1:thisNumCohorts])
  biolIndex<-grep(thisVar,biolLines)
  newLine<-paste(thisMums,collapse="\t")
  biolLines[biolIndex+1]<-newLine
}

#write it out
writeLines(biolLines,biolFile)
