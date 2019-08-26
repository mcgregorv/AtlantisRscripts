# read in tracers from runs with varie initial conditions, and caluclate ecological indicators, then asess how much these change wrt time and between runs
# use connectedmess, diversity (modified Kempton's Q), and mean trophic level

source(paste(DIR$'General functions',"calcQ.r", sep=""))
source(paste(DIR$'General functions',"calcTLindex.r", sep=""))
source(paste(DIR$'General functions',"chaosFunctions.R", sep=""))

# read in tracers from Chaos runs
basePath <-paste(DIR$'Base',"AtlantisModels\\chaos\\",sep="")
plotPath <- paste(DIR$'Reports',"Chaos\\Figures\\", sep="") ## paper version
plotPath <- paste(DIR$'Figures',"Chaos\\", sep="")

# bring in c("storeRelResponse", "storeRatioResponse", "storeAllTracersFromALlRuns"), created in modelRelationshipsAt50years_groupStructure_vs_stability.R
# using storeAlltracersFromALlRuns, which has dimensions ChaosAlt, run, group, timestep
load(paste(basePath, "AllTracersFromAllCHAOSRuns", sep=""))

groupsDF<- read.csv(paste(basePath,"..\\CRAM_groups.csv", sep="")); ng<-dim(groupsDF)[1]
thisGroupsDF <- groupsDF[groupsDF$Code!="DC",] # carrion not used

fishedIndex <- seq(2, 6, by=2); unfishedIndex <- seq(1,5,by=2)
## KEPMTON'S q #################################################################################
######### KEPMTON'S q ####################################################
##################### KEPMTON'S q ########################################################
########################### KEPMTON'S q ####################################################
Qtimes<-array(NA, dim=dim(storeAllTracersFromALlRuns)[c(1,2,4)])
nchaosAlts <- dim(storeAllTracersFromALlRuns)[1]; nruns <- dim(storeAllTracersFromALlRuns)[1]
for(c in 1:nchaosAlts){
  for(r in 1:nruns){
    xx<-storeAllTracersFromALlRuns[c,r,,]
    thisTest<-sum(xx, na.rm=TRUE)
    if(thisTest>0){
      NAtest<-apply(xx,2,sum, na.rm=TRUE); NAindex<-NAtest>0
      yy<-apply(xx[,NAindex],2, calcQ)
      # 
      # for(i in 1:dim(xx)[2]){
      #   thisQ<-calcQ(xx[,i])
      # }
      Qtimes[c,r,1:length(yy)]<-yy
    }
  }
}
colByC <- colorRampPalette(colors=c("red",myGold, myGreen,myAqua, myBlue,myPurple))(nchaosAlts)
thisMax <- max(Qtimes, na.rm=TRUE)
thisMin <- min(Qtimes, na.rm=TRUE)
plot(Qtimes[1,1,], type="l", col=myGrey_trans, ylim=c(thisMin, thisMax))
for(c in c(5,6)){
  for(r in 1:nruns){
    points(Qtimes[c,r,], type="l", col=colByC[c])
  }
}
CVbyChaosAlt <- apply(Qtimes, c(1,3), calc_cv)
CVbyTimestep <- apply(Qtimes, 3, calc_cv)



## MEAN TROPHIC INDEX ###############################################################################
####### MEAN TROPHIC INDEX ######################################################################
################### MEAN TROPHIC INDEX ######################################################
################################### MEAN TROPHIC INDEX ###########################################
#read in trophic levels
groupsTL<-read.csv(paste(basePath,"..\\inputs\\biol_prm\\CRAM_trophicLevels_isotopes.csv",sep=""))
groupsTL$TL<-groupsTL$Isotope; groupsTL$TL[is.na(groupsTL$Isotope)]<-groupsTL$TrophicLevel2[is.na(groupsTL$Isotope)]

asGroupIndex <- thisGroupsDF$NumCohorts>1
TLbyGroup <- groupsTL$TL[match(thisGroupsDF$Code, groupsTL$Code)]
TLbyGroup <- TLbyGroup[asGroupIndex]

TLindexArray<- array(NA, dim=c(dim(storeAllTracersFromALlRuns)[c(1,2,4)]))
for(c in 1:nchaosAlts){
  for(r in 1:nruns){
    temp<-storeAllTracersFromALlRuns[c,r,asGroupIndex,]
    TLindexArray[c,r,] <- calcTLindex(temp,TLbyGroup)
  }
}

TL_CVbyChaosAlt <- apply(TLindexArray, c(1,3), calc_cv)
TL_CVbyTimestep <- apply(TLindexArray, 3, calc_cv)



#########################
### PELAGIC BIOMASS / TOTAL BIOMASS
######################################
ageIndex <- thisGroupsDF$NumCohorts>1
biomByRunYear<-apply(storeAllTracersFromALlRuns[,,ageIndex,],c(1,2,4), sum, na.rm=TRUE)

test <- apply(storeAllTracersFromALlRuns[,,])

pelIndex<-thisGroupsDF$Code %in% c("PFS","PFL", "PFM", "MAC")
pelBiomass<-apply(storeAllTracersFromALlRuns[,,pelIndex,], c(1,2,4), sum)
# 
pelBiomOverCatch<-0* pelBiomass
for(c in 1:nchaosAlts){
  for(r in 1:nruns){
    thisBiom<-biomByRunYear[c,r,]
    if(sum(thisBiom, na.rm=TRUE)>0){
      thisB<-pelBiomass[c,r,]; 
      thisR<-thisB/thisBiom
      pelBiomOverCatch[c,r,]<-thisR
    }
  }
}
colByFished <- rep(c( myOrange, myBlue),3)
axisCex=1.5
pdf(paste(plotPath,"EcoIndicators.pdf",sep=""), height=10, width=8)
par(mfrow=c(3,1), mar=c(4,5,1.5,1), oma=c(5,0,0,0))
thisMax<-max(TLindexArray, na.rm=TRUE); thisMin <- min(TLindexArray,na.rm=TRUE)
plot(x=1865:2015,y=TLindexArray[1,1,], type="n", ylim=c(3.62, 3.72), ylab="Mean trophic level",xlab="", cex.axis=axisCex, cex.lab=2)
for(c in 1:nchaosAlts){
  thisCol<-paste(colByFished[c],"88", sep="")
  for(r in 1:nruns){
    points(x=1865:2015,y=TLindexArray[c,r,], type="l", col=thisCol)
  }
}
mtext("A:", side=3, adj=0)
thisMax<-max(Qtimes, na.rm=TRUE); thisMin <- min(Qtimes,na.rm=TRUE)
plot(x=1865:2015,y=Qtimes[1,1,], type="n", ylim=c(thisMin, thisMax), ylab="Kempton's Q",xlab="", cex.axis=axisCex, cex.lab=2)
for(c in 1:nchaosAlts){
  thisCol<-paste(colByFished[c],"88", sep="")
  for(r in 1:nruns){
    points(x=1865:2015,y=Qtimes[c,r,], type="l", col=thisCol)
  }
}
mtext("B:", side=3, adj=0)
thisMax<-max(pelBiomOverCatch, na.rm=TRUE); thisMin <- min(pelBiomOverCatch,na.rm=TRUE)
plot(x=1865:2015,y=pelBiomOverCatch[1,1,], type="n", ylim=c(0, 0.5), ylab="Biomass ratio:pelagic/all",xlab="", cex.axis=axisCex, cex.lab=2)
for(c in 1:nchaosAlts){
  thisCol<-paste(colByFished[c],"88", sep="")
  for(r in 1:nruns){
    points(x=1865:2015,y=pelBiomOverCatch[c,r,], type="l", col=thisCol)
  }
}
mtext("C:", side=3, adj=0)
par(xpd=NA, lend=1)
legend(legend=c("Unfished", "Fished"), ncol=2, col=colByFished[1:2], lwd=3, seg.len=3, lty=1, x="bottom", bty="n", inset=-0.5, cex=2)
dev.off()

