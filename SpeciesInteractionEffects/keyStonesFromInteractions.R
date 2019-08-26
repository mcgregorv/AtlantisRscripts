## data created in plotIntEffects.R

this_run<-"base"
basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
outPath<-paste(basePath, "base\\IntSims\\", sep="")
plotPath<-paste(basePath,"Figures\\intEffects\\",sep="")
# plotPath<-paste(DIR$'Figures',"\\Sensitivities\\", sep="") ## overwrites figure in paper

dataOutPath<-outPath

load(file=paste(dataOutPath,"interactionData",sep=""))
# this_out<-c("allEndPoints", "endPoint", "groupsDF", "asIndex", "ageStructuredNames", "TLindex", "biomassByFocusGroup")
## allEndPoints and endPoint are relative to the end point in the base model (%), so -20 means it ended 20% lower than the base model
## allEndPoints has dimensions ngroups, ngroups, nruns where ngroups is the total number of species groups (55), nruns is the number of model runs (3)
## first dim is the focus group, second is the response group
## because i found it helpful, i also have endPoint which is just the end point for the focus groups for each of the 3 runs (dim = 55, 3)
##groupsDF is just the data defining groups, i think you have this already, but i put it here incase its helpful
## asIndex is the age structured index which is helpful as I only did this for the age structured groups
## TLindex puts age structured groups into trophic level order (roughly, as the trophic levels defined are a bit sketchy)
## to get endPoint ordered by TL, do endPoint[asIndex,][TLindex,] (you need the asIndex first)
##  biomassByFocusGroup has all the biomass tracers for all groups and all runs including the base run
## has dimensions ngroups, nruns, ngroup, ntimesteps (55, 4, 55, 50) 
## first dim is focus group, second is runs with 4th the base run and 1-3 increasing additional mortality, 3rd dimension is response group, fith is time (years)
## groups are in order of groups csv file
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
calcImpact2<-function(x,f){
  ## f is the focus group, x is the vector of impacts
  y<-0
  for(g in 1:length(x)){
    if(g != f){
      if(!is.na(x[g])){
        y<-y + abs(x[g]/100)
      }
    }
  }
  return(y)
}
calcBiomassProportion<-function(x,f){
  thisBiom<-x[f]; totalBiom<-sum(x, na.rm=TRUE)
  thisP<-thisBiom/totalBiom
  return(thisP)
}
## use end point in biomass for impact
allEndPointsBiomass<-0*allEndPoints[,,3]
for(f in 1:ng){
    thisEPs<-allEndPoints[f,,3]
    thisBs<-biomassByFocusGroup[f,4,,50]
    thisEPBs<-thisEPs * thisBs
    allEndPointsBiomass[f,]<-thisEPBs
}

impactByFocusGroup<-rep(NA, ng); impactByFocusGroup2<-rep(NA, ng)
for(f in 1:ng){
  impactByFocusGroup[f]<-calcImpact(x=allEndPoints[f,,3], f)
  impactByFocusGroup2[f]<-calcImpact2(x=allEndPointsBiomass[f,], f)
}
# impactByFocusGroup<-impactByFocusGroup2

biomassPbyfocusGroup<-rep(NA, ng)
for(f in 1:ng){
  biomassPbyfocusGroup[f]<-calcBiomassProportion(x=biomassByFocusGroup[f,4,,50], f) * (abs(endPoint[f,3])/100)
}

keystoneByFocusGroup<-rep(NA, ng); keystoneByFocusGroup2<-rep(NA, ng); keystoneByFocusGroup3<-rep(NA, ng)
for(f in 1:ng){
  keystoneByFocusGroup[f]<-log((1/100)*impactByFocusGroup[f] * (1 - biomassPbyfocusGroup[f])) #original 1/100 makes fraction rather than percentage
  # keystoneByFocusGroup2[f]<-log((1/100)*impactByFocusGroup[f] * (1 - biomassPbyfocusGroup[f])^10) # difference
  keystoneByFocusGroup2[f]<-log((1/100)*impactByFocusGroup[f] * (1/biomassPbyfocusGroup[f])^(1/10)) #quotient
  keystoneByFocusGroup2[f]<-log((1/100)*impactByFocusGroup[f] * ((1- biomassPbyfocusGroup[f])/biomassPbyfocusGroup[f])^(1/4)) #quotient
  
}
# plot(x=keystoneByFocusGroup, y=keystoneByFocusGroup2)

# keystoneByFocusGroup<-keystoneByFocusGroup2
## plot impact on x and keystone on right
relImpact<-(impactByFocusGroup/max(impactByFocusGroup, na.rm=TRUE))[asIndex]
keyStoneAS<-keystoneByFocusGroup[asIndex]
biomProp<-biomassPbyfocusGroup[asIndex]

relImpact2<-(impactByFocusGroup2/max(impactByFocusGroup2, na.rm=TRUE))[asIndex]
keyStoneAS2<-keystoneByFocusGroup2[asIndex]

keystoneOrderNumber<-rev(order(keyStoneAS)); keystoneOrderNumber2<-rev(order(keyStoneAS2))
groupsDF$Code[asIndex][keystoneOrderNumber]

nag<-length(ageStructuredNames)
legendText<-c()
for(f in 1:nag){
  thisNum<-f; thisName<-ageStructuredNames[keystoneOrderNumber][f]
  thisLegend<-paste(thisNum, thisName,sep=" ")
  legendText<-c(legendText,thisLegend)
}

## do a heat map(isH) that shows biomass and effect influences on keystoneness
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
colRamp<-colorRampPalette(colors=c(myGold,myOrange,myPurple,myBlue))(101)
temp<-data.frame(cbind("Biomass"=biomProp, "Effect"=relImpact, "Keystone"=keyStoneAS))
thisMin<-min(temp$Keystone, na.rm=TRUE); thisMax<-max(temp$Keystone, na.rm=TRUE)
ksCol<-unlist(lapply(temp$Keystone, getColorForMap, doLog=FALSE))

plot(x=temp$Biomass, y=temp$Effect,col=ksCol,pch=20,cex=1.2)

thisMax<-max(biomProp, na.rm=TRUE); thisMin<-min(biomProp, na.rm=TRUE)
colByBiom<-unlist(lapply(biomProp, getColorForMap,doLog=TRUE))
# pdf(paste(plotPath,"keystoneByImpact_quotient.pdf",sep=""), width=7,height=5)
# pdf(paste(plotPath,"keystoneByImpact_difference.pdf",sep=""), width=7,height=5)
# pdf(paste(plotPath,"keystoneByImpact_libralatoMethod.pdf",sep=""), width=7,height=5)
par(mar=c(4.5,4.5,0.5,0.5))
plot(x=relImpact, y=keyStoneAS, type="n",xlab="Relative overall effect", ylab="Keystoneness", cex.lab=thisCex, cex.axis=thisCex)
points(x=seq(0.01,1,by=0.01), y=log(seq(0.01,1,by=0.01)), type="l", lty=2, col=myGrey)
text(seq(1,length(keyStoneAS)),x=relImpact[keystoneOrderNumber],y=keyStoneAS[keystoneOrderNumber],cex=1.2,col=colByBiom[keystoneOrderNumber])
legendText

dev.off()

par(mar=c(4.5,4.5,0.5,0.5))
plot(x=biomProp, y=keyStoneAS, type="n",xlab="Biomass proportion", ylab="Keystoneness", cex.lab=thisCex, cex.axis=thisCex)
text(seq(1,length(keyStoneAS)),x=biomProp[keystoneOrderNumber],y=keyStoneAS[keystoneOrderNumber],cex=1.4)


pdf(paste(plotPath,"ImpactByBP_quotient.pdf",sep=""), width=7,height=5)
# pdf(paste(plotPath,"ImpactByBP_difference.pdf",sep=""), width=7,height=5)
# pdf(paste(plotPath,"ImpactByBP_libralatoMethod.pdf",sep=""), width=7,height=5)
par(mar=c(4.5,4.5,0.5,0.5))
plot(x=temp$Biomass, y=temp$Effect, type="n",xlab="Biomass proportion", ylab="Relative overall effect", cex.lab=thisCex, cex.axis=thisCex)
text(seq(1,length(keyStoneAS)),x=temp$Biomass[keystoneOrderNumber],y=temp$Effect[keystoneOrderNumber],cex=1.4,col=colByBiom[keystoneOrderNumber])
dev.off()


legendCol<-c(rep(myBlue,4), rep("black", 33))
pdf(paste(plotPath,"keystoneByImpact_quotientLEGEND.pdf",sep=""), width=10,height=6)
# pdf(paste(plotPath,"keystoneByImpact_differenceLEGEND.pdf",sep=""), width=10,height=6)
# pdf(paste(plotPath,"keystoneByImpact_libralatoMethodLEGEND.pdf",sep=""), width=10,height=6)
par(mar=c(0,0,0,0))
makeBlankPlot()
legend(legendText,pch=NA,lwd=NA,lty=NA,col=legendCol,x="center",ncol=3,cex=1.4, bty="n")
dev.off()


BPlegendText<-c(1e-4, 1e-3, 1e-2, 1e-1); legendCol<-unlist(lapply(BPlegendText, getColorForMap, doLog=TRUE))
pdf(paste(plotPath,"biomPropLEGEND.pdf",sep=""), width=3,height=2)
par(mar=c(0,0,0,0))
makeBlankPlot()
legend(legend=BPlegendText,pch=15,lwd=NA,lty=NA,col=legendCol,x="center",ncol=1,cex=1.4, pt.cex=3, bty="n", title="Biomass proportion")
dev.off()


## where is orange roughy getting the high effect from..?
f=grep("ORH", groupsDF$Code)
test<-allEndPoints[f,,3]


# what does just log of overall effect look like..?
plot(x=relImpact, y=log(relImpact),type="n")
text(seq(1,length(keyStoneAS)),x=relImpact[keystoneOrderNumber],y=log(relImpact[keystoneOrderNumber]),cex=1.2)

# does biomass prop have any effect..?
plot(x=biomProp, y=keyStoneAS, type="n",xlab="Biomass proportion", ylab="Keystoneness", cex.lab=thisCex, cex.axis=thisCex,ylim=c(-4,0.6))
text(seq(1,length(keyStoneAS)),x=biomProp[keystoneOrderNumber],y=keyStoneAS[keystoneOrderNumber],cex=1.4)


thisIndex<-asIndex; thisIndex[groupsDF$Code=="SB"]<-FALSE
## SB
test<-allEndPoints[grep("SB",groupsDF$Code),,]
plot(test[thisIndex,3][rev(order(abs(test[thisIndex,3])))], type="h")
