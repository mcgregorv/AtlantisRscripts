basePath <-paste(DIR$'Base',"AtlantisModels\\chaos\\",sep="")

allData2fit <- read.csv(paste(DIR$'Tables',"Chaos_allData2fit.csv", sep=""))

groupsDF<- read.csv(paste(basePath,"..\\CRAM_groups.csv", sep="")); ng<-dim(groupsDF)[1]
keepGindex <- grep("DC", groupsDF$Code, invert = TRUE)
thisGroupsDF <- groupsDF[keepGindex,]; ng <- dim(thisGroupsDF)[1]
thisCex=1.65

plotPath <- paste(DIR$'Reports',"Chaos\\Figures\\", sep="") ## paper version
plotPath <- paste(DIR$'Figures',"Chaos\\", sep="")

chaosVersions <- c("DChaosNtracers", "DChaosNtracersFISH", paste("ChaosNtracers",c("SampleInformanceAll", "FISHSampleInformanceAll", "SampleKEYPlusBP", "FISHSampleKEYPlusBP"), sep="")); nchaosVersions <- length(chaosVersions)
chaosRunDescr <- c("A: All up or down - no fishing", "B: All up or down - with fishing", "C: Weighted randon - no fishing", "D: Weighted random - with fishing", 
                   "E: Weighted Keystone - no fishing", "F: Weighted keystone - with fishing")

nruns <- 35
nts <- 151; nlayers <-6

thisFactors <- colnames(allData2fit)[c(10, 12:21)]

# biomass pool factors
bp_factors <- c("TL",     "NumL1cons",   "B0",    "PropByTopPrey")

## test for correlation between explan, vars
nfactors <- length(bp_factors)
allCors <- array(NA, dim=c(nfactors, nfactors))
for(i in 1:nfactors){
  for(j in 1:nfactors){
    # if (i != j){
    y1 <- allData2fit[,c(bp_factors[i])]; y2 <- allData2fit[,c(bp_factors[j])]
    thisIndex <- !is.na(y1) & !is.na(y2)
    thisCor <- cor(x=y1[thisIndex], y=y2[thisIndex], method="spearman")
    allCors[i,j]<-thisCor
  }
  # }
}

rownames(allCors)<- bp_factors; colnames(allCors)<-bp_factors
write.csv(signif(allCors,2), paste(DIR$'Tables',"Chaos_BPCorrelations.csv", sep=""))

# just age-structured as some vars do not have values for biomass pool groups
asIndex <- thisGroupsDF$NumCohorts >1
asCodes <- thisGroupsDF$Code[asIndex]
asData2fit <- allData2fit[allData2fit$Code %in% asCodes,]

## test for correlation between explan, vars
nfactors <- length(thisFactors)
asCors <- array(NA, dim=c(nfactors, nfactors))
for(i in 1:nfactors){
  for(j in 1:nfactors){
    # if (i != j){
    y1 <- asData2fit[,c(thisFactors[i])]; y2 <- asData2fit[,c(thisFactors[j])]
    thisIndex <- !is.na(y1) & !is.na(y2)
    thisCor <- cor(x=y1[thisIndex], y=y2[thisIndex], method="spearman")
    asCors[i,j]<-thisCor
  }
  # }
}


rownames(asCors)<- thisFactors; colnames(asCors)<-thisFactors
write.csv(signif(asCors,2), paste(DIR$'Tables',"Chaos_ASCorrelations.csv", sep=""))

