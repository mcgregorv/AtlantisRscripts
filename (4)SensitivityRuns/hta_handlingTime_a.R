source(paste(DIR$'General functions',"createSensitivityRun.R",sep=""))
source(paste(DIR$'General functions',"get_biol_vector.R",sep=""))

this_run<-"base"

thisPath<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")

#variables to edit
editVar<-"hta_XXX"
version<-"A"  
version<-"B"  

changes<-c(1e-10,1e-5,1e-3,1e-2,1e-1) #version A
changes<-c(1e-100,1e-50,1e-30,1e-20,1e-15) #version B

proportionalChange<-TRUE #set to true if want the existing values to be multiplied by the change values, FALSE if the change values replace the existing ones
nextLine<-FALSE #some paramters are on the same line as the parameter name (in which case set this to FALSE), some are the next line

isVector<-FALSE

groupsFileName<-"..\\CRAM_Groups"
biolFileName<-"..\\CRAM_base_biol"

varMax<-1e+20

groupsDF<-read.csv(paste(thisPath,groupsFileName,".csv",sep=""))
editGroups<-groupsDF$Code[groupsDF$NumCohorts>1]

createSensitivityRun(thisRun,runLength,thisPath,editVar,version,editGroups,changes,proportionalChange,nextLine,isVector,changeVectorIndex=NULL,groupsFileName,biolFileName,varMax)

