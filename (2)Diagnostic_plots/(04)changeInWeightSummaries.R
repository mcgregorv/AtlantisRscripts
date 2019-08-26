##used code from pred prey interaction plots
this_run<-"base"
# this_out<-""

# this_out<-"Sensitivity\\KDENR_XXXE1"
this_out<-""
# # 

# outPath<-paste(DIR$'Base',"\\ATLANTIS\\Model_runs\\TBGB\\",this_run,"\\output",this_out,"\\",sep="")

# outPath<-paste(DIR$'Base',"\\ATLANTIS\\Model_runs\\TBGB\\",this_run,"\\Sensitivity\\KDENR\\outputKDENR_XXXD1\\",sep="")

thisPath<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(thisPath,"output",this_out,"\\",sep="")

plotPath<-paste(thisPath,"..\\Figures\\",this_run,"\\",this_out,"",sep="")


# daysTimeStep<-73
daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1880
fishingStartYear<-1900
modelStartYear<-1880

source(paste(DIR$'General functions',"create_change_in_weight.R",sep=""))
outdir=outPath
startyear=modelStartYear
toutinc=daysTimeStep
diet = TRUE
ncout<-"output"
fgfile<-paste(thisPath,"..\\CRAM_Groups.csv",sep="")
biolprm<-paste(thisPath,"..\\CRAM_base_biol",this_out,".prm",sep="")

changeInWeight<-create_change_in_weight(outdir, fgfile, biolprm, ncout, startyear, toutinc)

changeInWeightData<-changeInWeight$relWeight
# changeInWeightData<-output$relWeight



source(paste(DIR$'General functions',"plotChangeInWeight.R",sep=""))
source(paste(DIR$'General functions',"plotChangeInWeightSummary.R",sep=""))

# plotPath<-paste(thisPath,"Figures\\",sep="")

plotChangeInWeight(changeInWeightData,plotPath)
plotChangeInWeightSummary(changeInWeightData,plotPath)
