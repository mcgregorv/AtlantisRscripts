##looks at prey of given predator

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
plotPath<-paste(this_path,"Figures\\",this_run,"\\",sep="")

##we'll need the groups df 
groupsDF<-read.csv(paste(this_path,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
thisInitialConditionsFile<-paste(this_path,"CRAM_input_short.nc",sep="")
ThisIC.nc<-nc_open(thisInitialConditionsFile)
thisBiolFile<-paste(this_path,"CRAM_BH_hybrid_biol.prm",sep="")
biolLines<-readLines(thisBiolFile)
groupsTL<-read.csv(paste(this_path,"inputs\\biol_prm\\CRAM_trophicLevels_isotopes.csv",sep=""))

#get volume from ICs
volume<-ncvar_get(ThisIC.nc,"volume")[,,1]
nlayers<-dim(volume)[1]; nboxes<-dim(volume)[2]

load(paste(DIR$"Data","\\eating\\biomassIC",sep="")); #brings in biomassIC, dim= 55 10  6 30

##read in ageCohort linking
ageCohortLinking<-read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ageCohortLinking.csv",sep=""))

this_sampledAvails<-read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ExpectedValuesScaledForMax.csv",sep=""))[-1]

ppGroups<-as.character(read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ppGroups.csv",sep=""))[,1]);
npg<-length(ppGroups)

thisPred<-"ZM"
temp<-this_sampledAvails[grep(thisPred,ppGroups),]; colnames(temp)<-ppGroups

preyAvail<-temp[colSums(temp)>0]
biomassAvail<-0*preyAvail

for(p in 1:length(preyAvail)){
  thisVar<-names(preyAvail)[p]
  thisCode<-gsub("juv|ad","",thisVar)
  rowIndex<-ageCohortLinking[,1]==thisCode; colIndex<-ageCohortLinking[rowIndex,-1]==thisVar
  preyCohorts<-seq(1,10)[colIndex]
  thisPPGroupIndex<-ageCohortLinking[c(thisCode),]
  thisBiomassIC<-sum(biomassIC[groupsDF$Code==thisCode,preyCohorts,,],na.rm=TRUE)
  if(thisCode==thisVar){
    thisBiomassIC<-sum(biomassIC[groupsDF$Code==thisCode,1,,]*volume,na.rm=TRUE)
  }
  biomassAvail[p]<-thisBiomassIC*preyAvail[p]
  
  # biomassAvail[p]<-thisBiomassIC
}
sum(biomassAvail)

par(mar=c(4,5,1,4),mfrow=c(3,1))
plot(as.double(biomassAvail),type="h",xaxt="n",xlab="",ylab="")
par(las=2)
axis(at=seq(1,length(biomassAvail)),labels=names(biomassAvail),side=1)

par(new=TRUE)
plot(as.double(preyAvail),type="l",lty=2,col="red",yaxt="n",xaxt="n",ylab="",xlab="")
axis(at=pretty(preyAvail),labels=pretty(preyAvail),side=4,col.axis="red")

u3qIndex<-biomassAvail>summary(as.double(biomassAvail))[2]
plot(as.double(biomassAvail)[u3qIndex],type="h",xaxt="n",xlab="",ylab="")
par(las=2)
axis(at=seq(1,length(biomassAvail[u3qIndex])),labels=names(biomassAvail)[u3qIndex],side=1)


u2qIndex<-biomassAvail>summary(as.double(biomassAvail))[3]
plot(as.double(biomassAvail)[u2qIndex],type="h",xaxt="n",xlab="",ylab="")
par(las=2)
axis(at=seq(1,length(biomassAvail[u2qIndex])),labels=names(biomassAvail)[u2qIndex],side=1)


par(new=TRUE)
plot(as.double(preyAvail)[u2qIndex],type="l",lty=2,col="red",yaxt="n",xaxt="n",ylab="",xlab="")
axis(at=pretty(preyAvail[u2qIndex]),labels=pretty(preyAvail[u2qIndex]),side=4,col.axis="red")






##################################################
##test alt availabilites
testAvail<-preyAvail
# testAvail[grep("BFF",names(testAvail))] <- 1e-4
testAvail[grep("ZL",names(testAvail))] <- 3e-4
testAvail[grep("ZM",names(testAvail))] <- 5e-4
testAvail[grep("HOK",names(testAvail))] <- 0
testAvail[grep("MACad",names(testAvail))] <- 1e-4
# # # testAvail[grep("HOKad",names(testAvail))] <- 1e-3
# # testAvail[grep("CEPad",names(testAvail))] <- 0.1
# # # testAvail[grep("PFSad",names(testAvail))] <- 8e-3
# testAvail$DL<-5e-5; testAvail$DR<-5e-4
# testAvail$EIDad<-0.05

# 
changePreys<-c("ZL","ZM","HOK","MACad");
changes<-c(3e-4, 5e-4, 0, 1e-4); nc<-length(changePreys)
rbind(changePreys,changes)


biomassAvail<-0*testAvail
for(p in 1:length(testAvail)){
  thisVar<-names(testAvail)[p]
  thisCode<-gsub("juv|ad","",thisVar)
  rowIndex<-ageCohortLinking[,1]==thisCode; colIndex<-ageCohortLinking[rowIndex,-1]==thisVar
  preyCohorts<-seq(1,10)[colIndex]
  thisPPGroupIndex<-ageCohortLinking[c(thisCode),]
  thisBiomassIC<-sum(biomassIC[groupsDF$Code==thisCode,preyCohorts,,],na.rm=TRUE)
  if(thisCode==thisVar){
    thisBiomassIC<-sum(biomassIC[groupsDF$Code==thisCode,1,,]*volume,na.rm=TRUE)
  }
  biomassAvail[p]<-thisBiomassIC*testAvail[p]
  
  # biomassAvail[p]<-thisBiomassIC
}

par(mfrow=c(1,1))
u2qIndex<-biomassAvail>summary(as.double(biomassAvail))[3]
plot(as.double(biomassAvail)[u2qIndex],type="h",xaxt="n",xlab="",ylab="")
par(las=2)
axis(at=seq(1,length(biomassAvail[u2qIndex])),labels=names(biomassAvail)[u2qIndex],side=1)


par(new=TRUE)
plot(as.double(testAvail)[u2qIndex],type="l",lty=2,col="red",yaxt="n",xaxt="n",ylab="")
axis(at=pretty(testAvail[u2qIndex]),labels=pretty(testAvail[u2qIndex]),side=4,col.axis="red")

