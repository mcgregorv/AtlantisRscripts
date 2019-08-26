source(paste(DIR$'General functions',"createSensitivityRun.R",sep=""))
source(paste(DIR$'General functions',"get_biol_vector.R",sep=""))

this_run<-"base"

thisPath<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")

#variables to edit
# editVar<-"KI_XXX_T15"
editVar<-"KN_XXX"
version<-"A"
version<-"B"
version<-"C" #C is same as B but with KI_PL/PS=5, not 0.01 and higher predation on DR



nextLine<-FALSE #some paramters are on the same line as the parameter name (in which case set this to FALSE), some are the next line

isVector<-FALSE

groupsFileName<-"..\\CRAM_Groups"
biolFileName<-"..\\CRAM_base_biol"

varMax<-1000

changes<-c(0,0.01,0.1,0.5,1.0,5.0); proportionalChange<-FALSE #set to true if want the existing values to be multiplied by the change values, FALSE if the change values replace the existing ones #Version A, KI
changes<-c(0,0.01,0.1,0.5,1.0,5.0,10); proportionalChange<-TRUE #set to true if want the existing values to be multiplied by the change values, FALSE if the change values replace the existing ones  #Version A, KN
editGroups<-c("PL","PS","MB","DF","MA") #version A, KN 

#Version B, KN
changes<-c(0,0.01,0.1,0.5,1.0,5.0,10); proportionalChange<-TRUE #set to true if want the existing values to be multiplied by the change values, FALSE if the change values replace the existing ones  
editGroups<-c("PL","PS") #version B, KN 


createSensitivityRun(thisRun,runLength,thisPath,editVar,version,editGroups,changes,proportionalChange,nextLine,isVector,changeVectorIndex=NULL,groupsFileName,biolFileName,varMax,ICfile = "CRAM_input_short.nc")

