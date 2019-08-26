##looks at predators of given prey

source(paste(DIR$'General functions',"get_interaction_trophicLevel.R",sep=""))
source(paste(DIR$'General functions',"get_first_number.R",sep=""))
source(paste(DIR$'General functions',"myRounding.R",sep=""))
source(paste(DIR$'General functions',"getSNRN.R",sep=""))
source(paste(DIR$'General functions',"nonZeroMean.R",sep=""))
source(paste(DIR$'General functions',"getICNumbers.R",sep=""))

source(paste(DIR$'General functions',"get_interaction_spatial_bySpace.R",sep=""))
source(paste(DIR$'General functions',"get_interaction_gapeSize_byCohort.R",sep=""))

this_run<-"Base"

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
# plotPath<-paste(this_path,"Figures\\",this_run,"\\",sep="")
plotPath<-paste(DIR$'Figures',"Mortality\\testZoo", sep="")

this_out<-"BASE50yr_habDep"
outPath<-paste(this_path,"\\",this_run,"\\output",this_out,"\\",sep="")

diet_check<-read.csv(paste(outPath,"outputDietCheck.txt",sep=""),skip=0,sep=" ",header=TRUE)
predators<-unique(diet_check$Predator)

##we'll need the groups df 
groupsDF<-read.csv(paste(this_path,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
thisInitialConditionsFile<-paste(this_path,"CRAM_input_short.nc",sep="")
ThisIC.nc<-nc_open(thisInitialConditionsFile)
thisBiolFile<-paste(this_path,"CRAM_BH_hybrid_biol.prm",sep="")
biolLines<-readLines(thisBiolFile)
groupsTL<-read.csv(paste(this_path,"inputs\\biol_prm\\CRAM_trophicLevels_isotopes.csv",sep=""))

ThisNC.nc<-nc_open(paste(out_path,"output.nc", sep="")); thisVol<-ncvar_get(ThisNC.nc, "volume"); nts<-dim(thisVol)[3]

#get volume from ICs
volume<-ncvar_get(ThisIC.nc,"volume")[,,1]
nlayers<-dim(volume)[1]; nboxes<-dim(volume)[2]

load(paste(DIR$"Data","\\eating\\biomassIC",sep="")); #brings in biomassIC, dim= 55 10  6 30

##read in ageCohort linking
ageCohortLinking<-read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ageCohortLinking.csv",sep=""))

this_sampledAvails<-read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ExpectedValuesScaledForMax.csv",sep=""))[-1]

ppGroups<-as.character(read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ppGroups.csv",sep=""))[,1]);
npg<-length(ppGroups)

getBiomass<-function(code){
  thisPredName<-str_trim(groupsDF$Name[groupsDF$Code==code], side="both")
  thisTracer<-paste(thisPredName, "_N", sep=""); thisTracerData<-ncvar_get(ThisNC.nc, thisTracer); 
  if(length(dim(thisTracerData))==3){  
    thisBiomass<-apply(thisTracerData*thisVol, 3, sum)*X_CN * mg_2_tonne
  } else{
    thisBiomass<-apply(thisTracerData*thisVol[6,,], 2, sum)*X_CN * mg_2_tonne
  }
  return(thisBiomass)
}

#turn mortality (in diet_check) into approximately mg N eaten by multiply it by prey biomass
#first get all biomass'
predators<-as.character(groupsDF$Code[groupsDF$IsPredator==1]); npg<-length(predators)
storeBiomass<-array(NA, dim=c(ng, nts))
biomassEaten<-array(NA, dim=c(ng, npg, nts))
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; storeBiomass[g,]<-getBiomass(thisCode)
  thisEatenTemp<-diet_check[,c(1:3,grep(thisCode, colnames(diet_check)))]
  thisEaten<-tapply(as.double(thisEatenTemp[,c(as.character(thisCode))]), list(thisEatenTemp$Time, thisEatenTemp$Predator), sum, na.rm=TRUE)
  thisEatenBiomass<-array(NA, dim=c(npg, dim(thisEaten)[1]))
  for(i in 1:dim(thisEaten)[2]){
    thisEatenBiomass[i,]<-thisEaten[,i]*storeBiomass[g,1:(dim(thisEaten)[2])]
  }
  biomassEaten[g,,1:(dim(thisEaten)[1])]<-thisEatenBiomass
}
## biomassEaten is the actual tonnes eaten (approximately!)
#dimensions are  preys, predators, timesteps

## now plot for  prey which predators are eating them most
pdf(paste(plotPath,"ApproxTonnesEatenOfEachPrey.pdf", sep=""))
par(mfrow=c(3,1), mar=c(10,4,2,1))
for(p in 1:ng){
  thisPrey<-as.character(groupsDF$Code[p])
  preyIndex<-groupsDF$Code==thisPrey; thisPredatedBiomass<-biomassEaten[preyIndex,,]
  test<-rowSums(thisPredatedBiomass, na.rm=TRUE)
  index<-test>summary(test)[4]*0.5
  toPlot<-thisPredatedBiomass[index,]; this_npreds<-dim(toPlot)[1]
  if(sum(toPlot,na.rm=TRUE)>0){
    predCols<-colorRampPalette(colors=c(myGold,myOrange,"red", myPurple,myBlue,"midnightblue"))(this_npreds)
    ymax<-max(toPlot, na.rm=TRUE)
    if(length(dim(toPlot))==0){ #just one predator case
      plot(toPlot, type="l", col=predCols[1], lwd=2, ylim=c(0,ymax),xlab="Timestep (year)", ylab="Approximate tonnes eaten")
    }else{
      plot(toPlot[1,], type="l", col=predCols[i], lwd=2, ylim=c(0,ymax),xlab="Timestep (year)", ylab="Approximate tonnes eaten")
      for(i in 1:this_npreds){
        points(toPlot[i,], type="l", lwd=2, col=predCols[i])
      }
    }
    mtext(thisPrey,side=3,adj=0)
    # makeBlankPlot()
    par(xpd=NA)
    legend(legend=groupsDF$Code[groupsDF$IsPredator==1][index], col=predCols, lty=1, lwd=2, x="bottom", ncol=5, inset=-1.4)
  }
}
dev.off()

#############
#plot predator diets
pdf(paste(plotPath,"ApproxRealisedDietsZoo.pdf", sep=""), height=7)
par(mfrow=c(3,1), mar=c(10,4,2,1))
for(p in 1:npg){
  thisPred<-predators[p]
  myDiets<-biomassEaten[,p,]
  test<-apply(myDiets, 1, sum, na.rm=TRUE)
  test<-apply(myDiets, 1, max, na.rm=TRUE)
  index<-test>summary(test)[4]*0.01
  plotDiets<-myDiets[index,]
  if(length(dim(plotDiets))==0){
    index<-test>summary(test)[3]
    plotDiets<-myDiets[index,]
  }
  plotPreys<-as.character(groupsDF$Code[index]); npreys<-length(plotPreys)
  thisPreyColors<-colorRampPalette(colors=c(myGold,myOrange, "red", myPurple, myBlue,myAqua, myGreen))(npreys)
  ymax<-max(plotDiets, na.rm=TRUE)
  plot(plotDiets[1,], type="l", lwd=2, col=thisPreyColors[1], ylim=c(0,ymax), xlab="Timestep (years)", ylab="Approximate tonnes eaten")
  mtext(thisPred, side=3, adj=0)
  for(i in 2:npreys){
    points(plotDiets[i,], type="l", lwd=2, col=thisPreyColors[i])
  }
  par(xpd=NA)
  legend(legend=plotPreys, col=thisPreyColors, lty=1, lwd=2, x="bottom", ncol=5, inset=-1.4)
  
}
dev.off()




#################################################################################################################################################
#################################################################################################################################################

#######################
thisPrey<-"ASQjuv"
########################
temp<-this_sampledAvails[,grep(thisPrey,ppGroups)]; 
names(temp)<-ppGroups

preyAvail<-temp[temp>0]

plot(as.double(preyAvail),type="l",lty=2,col="red",xaxt="n",ylab="",xlab="")
par(las=2)
axis(at=seq(1,length(preyAvail)),labels=names(preyAvail),side=1)

thisSummary<-summary(as.double(preyAvail))

u3qIndex<-preyAvail>thisSummary[2]
plot(as.double(preyAvail)[u3qIndex],type="h",xaxt="n",xlab="",ylab="")
par(las=2)
axis(at=seq(1,length(preyAvail[u3qIndex])),labels=names(preyAvail)[u3qIndex],side=1)

#####################################
## go into more detail: for each predator, calculate the proportion of its diet this group makes up - 
# use realised diets from DIET_check. this doesn't split by ad/juv
thisPreyCode<-gsub("juv|ad","",thisPrey)

temp<-this_sampledAvails[,grep(thisPrey,ppGroups)]; 
names(temp)<-ppGroups

topAvails<-temp[temp>1e-5]; topPredators<-names(topAvails); ntp<-length(topPredators)
par(las=2)
plot(topAvails,type="h",lwd=5,lend=1,xaxt="n",xlab="",ylab="Proportion available")
axis(at=seq(1,length(topAvails)),labels=topPredators,side=1)

thisTime<-unique(diet_check$Time)
propInDiets<-array(NA,dim=c(length(thisTime),ntp))

for(p in 1:ntp){
  thisPred<-topPredators[p]
  thisCode<-gsub("juv|ad","",thisPred)
  thisAge<-gsub(thisCode,"",thisPred)
  #get this predators prey from diet_check
  if(thisAge==2){
    thisDietCheck<-diet_check[diet_check$Predator==thisCode & diet_check$Cohort==1,c(1,6:dim(diet_check)[2])]
  } else{
    thisDietCheck<-diet_check[diet_check$Predator==thisCode & diet_check$Cohort==0,c(1,6:dim(diet_check)[2])]
  }
  thisTotalDiets<-rowSums(thisDietCheck[,-1])
  propInDiets[1:length(thisTotalDiets),p]<-thisDietCheck[,c(thisPreyCode)]/thisTotalDiets
}
test<-colSums(propInDiets,na.rm=TRUE)
index<-test>=summary(test)[4]
topPropPreds<-topPredators[index]; ntpp<-length(topPropPreds)

ymax<-max(propInDiets[,index],na.rm=TRUE)

temp<-apply(propInDiets,2,max,na.rm=TRUE); maxPropInDiets<-temp[index]
plot(maxPropInDiets,type="h",lend=1,col=myBlue,lwd=3,xlab="",xaxt="n")
par(las=2)
axis(at=seq(1,length(maxPropInDiets), length.out=ntpp),topPropPreds,side=1)

par(las=1,mfrow=c(2,1),mar=c(3,4,0.5,1),oma=c(1,0,0,0))
colByPred<-colorRampPalette(colors=c(myYellow,myGreen,myAqua,myBlue,myPurple,"midnightblue"))(ntpp)
thisYMax<-max(propInDiets,na.rm=TRUE)
thisYMax<-0.0005
plot(propInDiets[,1],type="n",xlab="Year",ylab="Proportion of diet", ylim=c(0,ymax))

for(p in 1:ntpp){
  points(propInDiets[,index][,p],type="l",lwd=2,col=colByPred[p])
}
# par(mar=c(1,10,1,1))
# makeBlankPlot()
# legend(legend=topPropPreds[index],lwd=2,col=colByPred[index],x="center")

rowIndex<-ppGroups %in% topPropPreds
topPropPredAvails<-this_sampledAvails[rowIndex,grep(thisPrey,ppGroups)]; 

temp<-this_sampledAvails[,grep(thisPrey,ppGroups)]; 
names(temp)<-ppGroups

topAvails<-temp[temp>1e-5]; topPredators<-names(topAvails); ntp<-length(topPredators)
tppIndex<-topPredators %in% topPropPreds
par(las=2)
plot(topAvails,type="h",lwd=5,lend=1,xaxt="n",xlab="",ylab="Proportion available",col=myGrey_trans)
axis(at=seq(1,length(topAvails)),labels=topPredators,side=1)
points(x=seq(1,length(topAvails))[tppIndex],y=topAvails[tppIndex],type="h",lend=1,lwd=5,col=colByPred)




