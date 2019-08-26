basePath <- paste(DIR$'Base',"AtlantisModels\\chaos\\", sep="")

allTheRankings <- read.csv(paste(DIR$'Tables', "allTheRankings.csv", sep="")); nrankings <- dim(allTheRankings)[2]
groupsDF<- read.csv(paste(basePath,"..\\CRAM_groups.csv", sep="")); ng<-dim(groupsDF)[1]

# read biol pars, and biomass pool indiv, sizes
bp_sizes <- read.csv(paste(DIR$'Tables', "BiomassPoolIndividualSizes.csv", sep=""))

## a and b are the length to weight conversion pars, the other 3 are the 3-parameter VB growth pars
lw_pars<-read.csv(paste(basePath,"..\\inputs\\supporting\\length2weights.csv", sep="")); 
colnames(lw_pars)<-c("Code", "a", "b","Linf", "k", "t0")

allTheRankings$Linf <- NA
for(g in 1:ng){
  thisCode <- as.character(allTheRankings$Code[g])
  thisNumCohorts <- groupsDF$NumCohorts[groupsDF$Code==thisCode]
  if(thisNumCohorts==1){
    thisLinf <- bp_sizes$Size[bp_sizes$Code==thisCode]
  } else{
    thisLinf <- lw_pars$Linf[lw_pars$Code == thisCode]
  }
  if(length(thisLinf)==0){thisLinf <- NA}
  allTheRankings$Linf[g] <- thisLinf
}

write.csv(allTheRankings, file=paste(DIR$'Tables', "allTheRankings.csv", sep=""), row.names=FALSE)
