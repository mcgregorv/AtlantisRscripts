# we've done this for one run (run 3 with the most additional mortlaity)
## is there anything that can be done with the other 2 runs..? worth it..?
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

calcImpact<-function(x,f){
  ## f is the focus group, x is the vector of impacts
  y<-0
  for(g in 1:length(x)){
    if(g != f){
      if(!is.na(x[g])){
        y<-y + x[g]^2
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

getColorForMap<-function(x, doLog=FALSE){
  thisCol<-"white"
  if(!is.na(x)){
    y<-round((x-thisMin) / (thisMax - thisMin), 2) *100 +1
    if(doLog==TRUE){
      x1<-(x-thisMin)/ (thisMax-thisMin)
      z<-(log(abs(x1), base=10) -log(0.001, base=10))  / (log(1, base=10) - log(0.001, base=10))
      if(z<0){z<-0}
      y<-100 * z + 1
      
    }
    thisCol<-colRamp[y]
    # 
    # if(x==1){thisCol<-posColRamp[101]}
    # if(x==(-1)){thisCol<-negColRamp[101]}
    if(x==0){thisCol="white"}
    if(x>thisMax){thisCol=colRamp[101]}
    if(x<thisMin){thisCol=colRamp[1]}
  }
  return(thisCol)
}

## set up arrays to store effect and keystoneness and biomass proportions
effectArray<-array(NA, dim=c(ng,nruns)); keyArray<-effectArray; propBarray<-keyArray

for(r in 2:nruns){
  ## use end point in biomass for impact
  allEndPointsBiomass<-0*allEndPoints[,,r]
  for(f in 1:ng){
    thisEPs<-allEndPoints[f,,r]
    thisBs<-biomassByFocusGroup[f,4,,50]
    thisEPBs<-thisEPs * thisBs
    allEndPointsBiomass[f,]<-thisEPBs
  }
  
  impactByFocusGroup<-rep(NA, ng); impactByFocusGroup2<-rep(NA, ng)
  for(f in 1:ng){
    impactByFocusGroup[f]<-calcImpact(x=allEndPoints[f,,r], f)
  }
  # impactByFocusGroup<-impactByFocusGroup2
  
  biomassPbyfocusGroup<-rep(NA, ng)
  for(f in 1:ng){
    biomassPbyfocusGroup[f]<-calcBiomassProportion(x=biomassByFocusGroup[f,4,,50], f) * (abs(endPoint[f,r])/100)
  }
  
  keystoneByFocusGroup<-rep(NA, ng); 
  for(f in 1:ng){
    # keystoneByFocusGroup[f]<-log((1/100)*impactByFocusGroup[f] * (1/biomassPbyfocusGroup[f])^(1/1000)) #set exp to 1000 for ALL Runs, and for runs 2-3
    keystoneByFocusGroup[f]<-log((1/100)*impactByFocusGroup[f] * (1 - biomassPbyfocusGroup[f])^1)
  }
  effectArray[,r]<-impactByFocusGroup; keyArray[,r]<-keystoneByFocusGroup; propBarray[,r]<-biomassPbyfocusGroup
  
}

## add them up for each group
keystoneByFocusGroup<-apply(keyArray,1,mean, na.rm=TRUE)
impactByFocusGroup<-apply(effectArray,1,mean, na.rm=TRUE)
biomassPbyfocusGroup<-apply(propBarray,1,mean, na.rm=TRUE)

# plot(x=keystoneByFocusGroup, y=keystoneByFocusGroup2)
## plot impact on x and keystone on right
relImpact<-(impactByFocusGroup/max(impactByFocusGroup, na.rm=TRUE))[asIndex]
keyStoneAS<-keystoneByFocusGroup[asIndex]
biomProp<-biomassPbyfocusGroup[asIndex]


# save(list=c("keyStoneAS"),file=paste(dataOutPath,"keyStoneAS",sep=""))


keystoneOrderNumber<-rev(order(keyStoneAS)); 
groupsDF$Code[asIndex][keystoneOrderNumber]

nag<-length(ageStructuredNames)
legendText<-c()
for(f in 1:nag){
  thisNum<-f; thisName<-ageStructuredNames[keystoneOrderNumber][f]
  thisLegend<-paste(thisNum, thisName,sep=" ")
  legendText<-c(legendText,thisLegend)
}


## do a heat map(isH) that shows biomass and effect influences on keystoneness
colRamp<-colorRampPalette(colors=c(myGold,myOrange,myPurple,myBlue))(101)
temp<-data.frame(cbind("Biomass"=biomProp, "Effect"=relImpact, "Keystone"=keyStoneAS))
thisMin<-min(temp$Keystone, na.rm=TRUE); thisMax<-max(temp$Keystone, na.rm=TRUE)
ksCol<-unlist(lapply(temp$Keystone, getColorForMap, doLog=FALSE))

plot(x=temp$Biomass, y=temp$Effect,col=ksCol,pch=20,cex=1.2)

thisMax<-max(biomProp, na.rm=TRUE); thisMin<-min(biomProp, na.rm=TRUE)
colByBiom<-unlist(lapply(biomProp, getColorForMap,doLog=TRUE))
pdf(paste(plotPath,"keystoneByImpactRUNS2and3.pdf",sep=""), width=7,height=5)
par(mar=c(4.5,4.5,0.5,0.5))
plot(x=relImpact, y=keyStoneAS, type="n",xlab="Relative overall effect", ylab="Keystoneness", cex.lab=thisCex, cex.axis=thisCex)
# points(x=seq(0.01,1,by=0.01), y=log(seq(0.01,1,by=0.01)), type="l", lty=2, col=myGrey)
text(seq(1,length(keyStoneAS)),x=relImpact[keystoneOrderNumber],y=keyStoneAS[keystoneOrderNumber],cex=1.2,col=colByBiom[keystoneOrderNumber])
legendText

dev.off()

par(mar=c(4.5,4.5,0.5,0.5))
plot(x=biomProp, y=keyStoneAS, type="n",xlab="Biomass proportion", ylab="Keystoneness", cex.lab=thisCex, cex.axis=thisCex)
text(seq(1,length(keyStoneAS)),x=biomProp[keystoneOrderNumber],y=keyStoneAS[keystoneOrderNumber],cex=1.4)


legendCol<-c(rep(myBlue,4), rep("black", 33))
# pdf(paste(plotPath,"keystoneByImpactALLRUNSLEGEND.pdf",sep=""), width=10,height=6)
pdf(paste(plotPath,"keystoneByImpactRUNS2and3LEGEND.pdf",sep=""), width=10,height=6)
par(mar=c(0,0,0,0))
makeBlankPlot()
legend(legendText,pch=NA,lwd=NA,lty=NA,col=legendCol,x="center",ncol=3,cex=1.4, bty="n")
dev.off()


