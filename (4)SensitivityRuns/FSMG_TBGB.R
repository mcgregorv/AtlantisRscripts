source(paste(DIR$'General functions',"createSensitivityRunTBGB.R",sep=""))
source(paste(DIR$'General functions',"get_biol_vector.R",sep=""))

this_run<-"base"

thisPath<-paste(DIR$'Base',"TBGB\\TBGB_SI\\",sep="")

#variables to edit
editVar<-"^XXX_FSMG"
version<-"A"

changes<-c(1,2,5,10,50)

proportionalChange<-TRUE #set to true if want the existing values to be multiplied by the change values, FALSE if the change values replace the existing ones
nextLine<-TRUE #some paramters are on the same line as the parameter name (in which case set this to FALSE), some are the next line

isVector<-TRUE

groupsFileName<-"TBGB_Groups"
biolFileName<-"TBGB_biol"
forceFile<-"TBGB_force"

varMax<-1e+20

groupsDF<-read.csv(paste(thisPath,groupsFileName,".csv",sep=""))
editGroups<-c("BAR") #V A-C

changeVectorIndex<-seq(1,2) #distribA


createSensitivityRunTBGB(this_run,runLength,thisPath,editVar,version,editGroups,runPrmFile="TBGB_run",changes,proportionalChange,nextLine,isVector,changeVectorIndex=changeVectorIndex,groupsFileName,biolFileName,varMax,ICfile="create_biol_input/TBGB_input.nc", forceFile=forceFile)
