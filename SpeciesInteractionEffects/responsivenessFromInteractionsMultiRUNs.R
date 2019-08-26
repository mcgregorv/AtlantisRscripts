# like keystoneness but from the responding groups perspective
this_run<-"base"
basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
outPath<-paste(basePath, "base\\IntSims\\", sep="")
plotPath<-paste(basePath,"Figures\\intEffects\\",sep="")
# plotPath<-paste(DIR$'Figures',"\\Sensitivities\\", sep="") ## overwrites figure in paper

dataOutPath<-outPath

load(file=paste(dataOutPath,"interactionData",sep=""))
groupsDFPaper<-read.csv(paste(this_path,"..\\CRAM_groupsPaper.csv", sep=""))
ageStructuredNames<-gsub("_", " ", groupsDFPaper$Name[groupsDF$NumCohorts>1])


nruns=3

ng<-dim(endPoint)[1]

calcResponse<-function(x,f){
  ## f is the focus group, x is the vector of impacts
  y<-0
  for(g in 1:length(x)){
    if(g != f){
      if(!is.na(x[g])){
        y<-y + ( x[g]^2 * biomassPbyfocusGroup[g] )
      }
    }
  }
  z<-sqrt(y)
  return(z)
}
calcBiomassProportion<-function(x,f){
  thisBiom<-x[f]; totalBiom<-sum(x, na.rm=TRUE)
  thisP<-thisBiom/totalBiom
  return(thisP)
}

responsArray<-array(NA, dim=c(ng, nruns))
for(r in 2:3){
  ## use end point in biomass for impact
  allEndPointsBiomass<-0*allEndPoints[,,r]
  for(f in 1:ng){
    thisEPs<-allEndPoints[f,,r]
    thisBs<-biomassByFocusGroup[f,4,,50]
    thisEPBs<-thisEPs * thisBs
    allEndPointsBiomass[f,]<-thisEPBs
  }
  biomassPbyfocusGroup<-rep(NA, ng)
  for(f in 1:ng){
    biomassPbyfocusGroup[f]<-calcBiomassProportion(x=biomassByFocusGroup[f,4,,50], f) * (abs(endPoint[f,r])/100)
  }
  
  impactByResponseGroup<-rep(NA, ng); 
  for(f in 1:ng){
    impactByResponseGroup[f]<-calcResponse(x=allEndPoints[,f,r], f)
  }
  responsArray[,r]<-impactByResponseGroup
}
impactByResponseGroup<-apply(responsArray,1,mean, na.rm=TRUE)
keyStoneAS<-impactByResponseGroup[asIndex]

keystoneOrderNumber<-rev(order(keyStoneAS)); 
groupsDF$Code[asIndex][keystoneOrderNumber]

nag<-length(ageStructuredNames)
legendText<-c()
for(f in 1:nag){
  thisNum<-f; thisName<-ageStructuredNames[keystoneOrderNumber][f]
  thisLegend<-paste(thisNum, thisName,sep=" ")
  legendText<-c(legendText,thisLegend)
}


pdf(paste(plotPath,"responsivenessRUNS2and3.pdf",sep=""), width=7,height=5)
par(lend=1, mar=c(10.5,4,1,1))
plot(keyStoneAS[keystoneOrderNumber], type="h", lwd=5, xaxt="n", xlab="", ylab="Responsiveness", col=myGrey)
par(las=2)
axis(at=seq(1,length(keyStoneAS)), labels=legendText, side=1)
dev.off()


responsiveness<-keyStoneAS; 
save(list=c("responsiveness"),file=paste(dataOutPath,"responsiveness",sep=""))
