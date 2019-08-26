##looks at prey of given predator
# original version of this uses initial conditions - this one uses output.nc from model run, at a specified timestep

source(paste(DIR$'General functions',"get_interaction_trophicLevel.R",sep=""))
source(paste(DIR$'General functions',"get_first_number.R",sep=""))
source(paste(DIR$'General functions',"myRounding.R",sep=""))
source(paste(DIR$'General functions',"getSNRN.R",sep=""))
source(paste(DIR$'General functions',"nonZeroMean.R",sep=""))
source(paste(DIR$'General functions',"getICNumbers.R",sep=""))

source(paste(DIR$'General functions',"get_interaction_spatial_bySpace.R",sep=""))
source(paste(DIR$'General functions',"get_interaction_gapeSize_byCohort.R",sep=""))

this_run<-"..\\TBGB\\TBGB_JP"

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
plotPath<-paste(this_path,"Figures\\",this_run,"\\",sep="")

##we'll need the groups df 
groupsDF<-read.csv(paste(this_path,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

this_out<-"BaseLong3"
this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
volume<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
ntsteps<-dim(thisVol)[3]

thisBiolFile<-paste(this_path,"CRAM_BH_hybrid_biol.prm",sep="")
biolLines<-readLines(thisBiolFile)
groupsTL<-read.csv(paste(this_path,"inputs\\biol_prm\\CRAM_trophicLevels_isotopes.csv",sep=""))

nlayers<-dim(volume)[1]; nboxes<-dim(volume)[2]

##read in ageCohort linking
ageCohortLinking<-read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ageCohortLinking.csv",sep=""))
# ageCohortLinking[1:3,] ## this is a table with species codes in the first column followed by 10 columns that show if each age class (1-10) is juvenile or adult
# X     X1     X2    X3    X4    X5    X6    X7    X8    X9   X10
# 1 ASQ ASQjuv  ASQad  <NA>  <NA>  <NA>  <NA>  <NA>  <NA>  <NA>  <NA>
# 2 BAL BALjuv BALjuv BALad BALad BALad BALad BALad BALad BALad BALad
# 3  BB     BB   <NA>  <NA>  <NA>  <NA>  <NA>  <NA>  <NA>  <NA>  <NA>
  
this_sampledAvails<-read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ExpectedValuesScaledForMax.csv",sep=""))[-1]

## ppGroups is all groupCodes with ad and juv added if age-structured
ppGroups<-as.character(read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ppGroups.csv",sep=""))[,1]);
npg<-length(ppGroups)
# ppGroups
# [1] "ASQad"  "ASQjuv" "BALad"  "BALjuv" "BB"     "BC"     "BD"  ...

thisPred<-"ELI"
temp<-this_sampledAvails[grep(thisPred,ppGroups),]; colnames(temp)<-ppGroups
temp2<-apply(temp,2,sum)

preyAvail<-temp2[temp2>0]
biomassAvail<-0*preyAvail

# preyAvail <- ppGroups; biomassAvail <- rep(0, npg)


thisTimeStep<-50

for(p in 1:length(preyAvail)){
  # thisVar<-names(preyAvail)[p]
  thisVar<-preyAvail[p]
  thisCode<-gsub("juv|ad","",thisVar)
  thisPreyName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode], side="both")
  thisPreyNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisCode]
  if(thisPreyNumCohorts==1){
    thisTracerVar<-paste(thisPreyName,"_N", sep="")
    tracerData<-ncvar_get(ThisNC.nc, thisTracerVar)
    if(length(dim(tracerData))==3){
      thisPreyBiomass<-apply(tracerData*thisVol,3,sum)[thisTimeStep]
    } else{
      thisPreyBiomass<-apply(tracerData*thisVol[nlayers,,],2,sum)[thisTimeStep]
    }
    
  }else{
    #need to loop through cohorts and add them up
    
    rowIndex<-ageCohortLinking[,1]==thisCode; colIndex<-ageCohortLinking[rowIndex,-1]==thisVar
    preyCohorts<-seq(1,10)[colIndex]
    thisPreyBiomass<-0
    for(c in preyCohorts){
      if(!is.na(c)){
        # NUMBERS
        thisNumsVar<-paste(thisPreyName,c,"_Nums", sep="")
        xxNums<-ncvar_get(ThisNC.nc, thisNumsVar)
        thisNums<-apply(xxNums, 3, sum)[thisTimeStep]
        # RESN
        thisResNVar<-paste(thisPreyName,c,"_ResN", sep="")
        xxResN<-ncvar_get(ThisNC.nc, thisResNVar)
        thisResN<-apply(xxResN, 3, nonZeroMean)[thisTimeStep]
        # STRUCTN
        thisStructNVar<-paste(thisPreyName,c,"_StructN", sep="")
        xxStructN<-ncvar_get(ThisNC.nc, thisStructNVar)
        thisStructN<-apply(xxStructN, 3, sum)[thisTimeStep]
        #biomass
        thisPreyBiomass<-thisPreyBiomass+(thisStructN + thisResN) * thisNums
      }
    }
    
  }
  biomassAvail[p]<-thisPreyBiomass*preyAvail[p]
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

