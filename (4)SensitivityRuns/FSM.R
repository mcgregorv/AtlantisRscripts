source(paste(DIR$'General functions',"createSensitivityRun.R",sep=""))
source(paste(DIR$'General functions',"get_biol_vector.R",sep=""))

this_run<-"base"; runLength<-"short"

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
biolLines<-readLines(paste(this_path,"..\\CRAM_BH_hybrid_biol.prm",sep=""))
groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

#FIRST check which groups migrate out of the model domain
doesMigrate<-rep(NA,ng)
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]
  thisVar<-paste("flag", thisCode, "Migrate", sep=""); xx<-biolLines[grep(thisVar, biolLines)]
  thisFlag<-get_first_number(xx)
  if(length(thisFlag)>0){ doesMigrate[g]<-thisFlag}
}

temp<-groupsDF$Code[doesMigrate==1]; groupsThatMigrate<-temp[!is.na(temp)]

#increase FSMG on these groups
groups1<-c("BAL", "ELP", "CET"); 
# groups1<-c("PIN", "HOK", "PFL", "SB", "CET", "MAC", "BAL")

ngroups<-length(groups1)

#variables to edit
editVar<-"XXX_FSMG"
version<-"A"

proportionalChange<-FALSE #set to true if want the existing values to be multiplied by the change values, FALSE if the change values replace the existing ones
nextLine<-TRUE #some paramters are on the same line as the parameter name (in which case set this to FALSE), some are the next line

isVector<-FALSE

groupsFileName<-"..\\CRAM_Groups"
biolFileName<-"..\\CRAM_BH_hybrid_biol"
forceFile<-"CRAM_forceBURNIN1865"


varMax<-10
changes<-c(1, 1.2, 1.5, 2, 5, 10) #Version A #changes them all at the same time, to this value
changes<-c(1,1.2,1.5,1.8,2,2.5) #distribA
changes<-c(0.15, seq(0.2,1,by=0.1)) #distribB

editGroups<-c("BAL") #version A
createSensitivityRun(this_run,runLength,thisPath,editVar,version,editGroups,runPrmFile="CRAM_base_run50yr",changes,proportionalChange,nextLine,isVector,changeVectorIndex=changeVectorIndex,groupsFileName,biolFileName,varMax,ICfile="CRAM_input_short.nc", forceFile=forceFile)

