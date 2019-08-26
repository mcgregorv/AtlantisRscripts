#read in tracers from the base model and store ready to check out gape sizes
this_run<-"base"
this_out<-"Base"
source(paste(DIR$'General functions',"nonZeroMean.R", sep=""))

mg_2_tonne<-2e-8; X_CN<-5.7

burnin<-35 #number of years to skip for burnin

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
# outPath<-paste(this_path,"ZmQpermute\\",this_out,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")

groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1900
fishingStartYear<-1900
modelStartYear<-1900

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
nts<-dim(thisVol)[3]-burnin

#storetracers in an array - one for sn, one for Rn, one for numbrs, and one for tonnes
storeRNtracers<-array(NA, dim=c(nts, ng, 10)); storeSNtracers<-0*storeRNtracers
storeNumtracers<-0*storeRNtracers
storeTonnesTracers<-array(NA, dim=c(nts, ng))

#loop through the species groups and populate
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]; thisName<-str_trim(groupsDF$Name[g])
  # for all groups, get biomass in tonnes
  thisTracer<-paste(thisName,"_N", sep=""); thisData<-ncvar_get(ThisNC.nc, thisTracer)
  #if has 3 dimensions, then multiply by volume; else just by the area
  if(length(dim(thisData))==3){
    xx<-apply(thisData * thisVol,3, sum, na.rm=TRUE) * mg_2_tonne * X_CN
  } else{
    xx<-apply(thisData * thisVol[nlayers,,],2, sum, na.rm=TRUE) * mg_2_tonne * X_CN
  }
  #put in storetracers
  storeTonnesTracers[,g]<-xx[(burnin+1):(nts+burnin)]
  #if age structured, get numbers and weights as well
  if(thisNumCohorts > 1){
    for(c in 1:thisNumCohorts){
      thisTracer<- paste(thisName,c,"_Nums",sep=""); thisData<-ncvar_get(ThisNC.nc, thisTracer)
      xx<-apply(thisData,3,sum, na.rm=TRUE)
      storeNumtracers[,g,c]<-xx[(burnin+1):(burnin+nts)]
      #Res N
      thisTracer<-paste(thisName, c, "_ResN", sep=""); thisData<-ncvar_get(ThisNC.nc, thisTracer)
      #these are the same across all cells where they exist, so use nonzeromean to get the mean of the cells for which they are in
      xx<-apply(thisData, 3, nonZeroMean)
      storeRNtracers[,g,c]<-xx[(burnin+1):(burnin+nts)]
      # StructN
      thisTracer<-paste(thisName, c, "_StructN", sep=""); thisData<-ncvar_get(ThisNC.nc, thisTracer)
      #these are the same across all cells where they exist, so use nonzeromean to get the mean of the cells for which they are in
      xx<-apply(thisData, 3, nonZeroMean)
      storeSNtracers[,g,c]<-xx[(burnin+1):(burnin+nts)]
    }
  }
}

## also need diet summary by juv-ad
dietCheck<-read.csv(paste(outPath, "outputDietCheck.txt", sep=""), sep=" ")
#time is in days - turn to years, then cut off the burnin
dietYears<-dietCheck$Time/365
dietCheck<-dietCheck[dietYears>burnin,]
#summaries diets by cohort and group - they are proportions, so mean should be appropriate
tempData<-dietCheck[,c("Predator", "Cohort", as.character(groupsDF$Code))]
meltData<-data.frame(melt(tempData, id.vars=c("Predator", "Cohort")))
meltData$Cohort<-factor(meltData$Cohort)
meltData$value<-as.double(as.character(meltData$value))
dietSummary<-tapply(meltData$value,meltData[,c("Predator", "Cohort", "variable")], mean, na.rm=TRUE)

## read in KLP and KUP values from biol.prm file and store these too
biolFile<-paste(this_path,"..\\CRAM_BH_hybrid_biol.prm", sep="")
biolLines<-readLines(biolFile)
storeGapeScalars<-data.frame(matrix(NA, ncol=3, nrow=ng)); colnames(storeGapeScalars)<-c("Code", "KLP", "KUP")
storeGapeScalars$Code<-as.character(groupsDF$Code)
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    thisVar<-paste("KLP_", thisCode,sep=""); thisLine<-biolLines[grep(thisVar, biolLines)]; thisValue<-get_first_number(thisLine)
    storeGapeScalars$KLP[g]<-thisValue
    thisVar<-paste("KUP_", thisCode,sep=""); thisLine<-biolLines[grep(thisVar, biolLines)]; thisValue<-get_first_number(thisLine)
    storeGapeScalars$KUP[g]<-thisValue
    
  }
}

## write them out so can just read in.
storeTracersFile<-paste(this_path,"GapeSizes\\baseModelTracers", sep="")
save(list=c("storeTonnesTracers", "storeRNtracers", "storeSNtracers", "storeNumtracers", "dietSummary", "storeGapeScalars"),file=storeTracersFile)








