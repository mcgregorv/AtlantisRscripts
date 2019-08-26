# read in base model files, drop a group and create a new model without it
## testing if low keystoneness groups can be dropped without consequence
## would run faster :) 
basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

## run command
runFile<-paste(basePath,"RUN_base_short.bat", sep="")
temp<-readLines(runFile); runCommand<- temp[grep("atlantis", temp)]

getFileName<-function(inputCode){
  x<-unlist(str_split(runCommand, " ")); xx<-x[(grep(inputCode, x) +1)]
  return(xx)
}

baseICfile<-getFileName(inputCode="-i")
baseRunFile<-getFileName(inputCode="-r")
baseForceFile<-getFileName(inputCode="-f")
basePhysicsFile<-getFileName(inputCode="-p")
baseBiolFile<-getFileName(inputCode="-b")
baseGroupsFile<-getFileName(inputCode="-s")
baseFisheriesFile<-getFileName(inputCode="-q")
baseHarvestFile<-getFileName(inputCode="-h")
baseOutputFile<-getFileName(inputCode="-d")

groupsDF<-read.csv(paste(basePath, baseGroupsFile, sep="")); ng<-dim(groupsDF)[1]

DROPcode<-"CRA"

# turn this group off in the groups csv file
newGroupsDF<-groupsDF; newGroupsDF$IsTurnedOn[newGroupsDF$Code==DROPcode]<-0
newGroupsFile<-gsub("\\.csv",paste("DROP",DROPcode,".csv", sep=""),baseGroupsFile)
write.csv(newGroupsDF, file=paste(basePath, newGroupsFile, sep=""), row.names = FALSE)

## edit number of species recorded in run.prm file
newRunFile<-gsub("\\.prm", paste("DROP", DROPcode, ".prm", sep=""), baseRunFile)
runLines<-readLines(paste(basePath, baseRunFile, sep=""))
x<-grep("K_num_tot_sp", runLines); base_kNum<-get_first_number(runLines[x], n=1)
newkNum<-base_kNum-1
newkNum<-base_kNum
newLine<-paste("K_num_tot_sp", newkNum, collapse="\t")
runLines[x]<-newLine
writeLines(runLines, paste(basePath,newRunFile, sep=""))

newRunCommand<-gsub(baseGroupsFile, newGroupsFile, runCommand)
newRunCommand <- gsub(baseRunFile, newRunFile, newRunCommand)

newOutputFile<-paste("outputDROP",DROPcode,sep="")
newRunCommand <- gsub(baseOutputFile, newOutputFile, newRunCommand, fixed=TRUE)

newBatFile<-paste(basePath,"RUN_DROP",DROPcode, ".bat", sep="")
writeLines(newRunCommand, newBatFile)
