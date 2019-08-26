#read in sampled availabilities from first version of this file, edit it then create pPREY lines for biol.prm file
this_sampledAvails<-read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ExpectedValuesScaledForMax.csv",sep=""))[,-1]
this_path<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
ppGroups<-as.character(read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ppGroups.csv",sep=""))[,1]); npg<-length(ppGroups)
groupsDF<-read.csv(paste(this_path,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
ndg<-ng+3 #number of diet groups, include sediment version of DL, DR, DC

adIndex<-grep("ad",ppGroups); juvIndex<-grep("juv",ppGroups)
naIndex<-grep("ad|juv",ppGroups,invert = TRUE)
not_adultIndex<-grep("ad",ppGroups,invert = TRUE)


scaled_sampledAvails<- this_sampledAvails

rowIndex<-ppGroups=="ZM"; colIndex<-ppGroups=="PS"
scaled_sampledAvails[rowIndex,colIndex]<-1e-4

rowIndex<-ppGroups=="ZL"; colIndex<-ppGroups=="ZG"
scaled_sampledAvails[rowIndex,colIndex]<-0

rowIndex<-grep("SSOad", ppGroups); colIndex<-ppGroups=="ZG"
scaled_sampledAvails[rowIndex,colIndex]<-0.2

rowIndex<-grep("SPE", ppGroups); colIndex<-grep("IVS",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-2*scaled_sampledAvails[rowIndex,colIndex]

rowIndex<-grep("SPDad", ppGroups); colIndex<-grep("HOK",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-2*scaled_sampledAvails[rowIndex,colIndex]

rowIndex<-grep("SB", ppGroups); colIndex<-grep("PFS",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-10*scaled_sampledAvails[rowIndex,colIndex]

rowIndex<-grep("PFL", ppGroups); colIndex<-grep("MJE|MAC",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-scaled_sampledAvails[rowIndex,colIndex]*0.2

rowIndex<-grep("LDO", ppGroups); colIndex<-grep("BIS",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-3*scaled_sampledAvails[rowIndex,colIndex]

rowIndex<-grep("LDO", ppGroups); colIndex<-grep("EID",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-0.01*scaled_sampledAvails[rowIndex,colIndex]

rowIndex<-grep("JAV", ppGroups); colIndex<-grep("PFSjuv",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-0.015
rowIndex<-grep("JAVad", ppGroups); colIndex<-grep("PFSad",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-0.0015

rowIndex<-grep("HAK", ppGroups); colIndex<-grep("JAV",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-2*scaled_sampledAvails[rowIndex,colIndex]
rowIndex<-grep("HAK", ppGroups); colIndex<-grep("MJE|BOE|SSO",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-0.5*scaled_sampledAvails[rowIndex,colIndex]

rowIndex<-grep("GSHjuv", ppGroups); colIndex<-grep("IVS",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-2* scaled_sampledAvails[rowIndex,colIndex]

rowIndex<-grep("ETB", ppGroups); colIndex<-grep("ORH",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-2*scaled_sampledAvails[rowIndex,colIndex]
colIndex<-grep("GSH",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-0.1*scaled_sampledAvails[rowIndex,colIndex]

rowIndex<-grep("ELP", ppGroups); colIndex<-grep("MJE",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-0.1*scaled_sampledAvails[rowIndex,colIndex]
colIndex<-grep("LIN",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-0.05*scaled_sampledAvails[rowIndex,colIndex]

rowIndex<-grep("ELI", ppGroups); colIndex<-grep("ELP|MJE|LIN|MAC|LDO",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-0.1*scaled_sampledAvails[rowIndex,colIndex]

rowIndex<-grep("DPI", ppGroups); colIndex<-grep("LIN",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-0.05*scaled_sampledAvails[rowIndex,colIndex]
colIndex<-grep("BIS|CEP|EIS",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-4*scaled_sampledAvails[rowIndex,colIndex]

rowIndex<-grep("CRA", ppGroups); colIndex<-grep("IVS",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-3*scaled_sampledAvails[rowIndex,colIndex]

rowIndex<-grep("CEP", ppGroups); colIndex<-grep("ZG",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-0.01*scaled_sampledAvails[rowIndex,colIndex]
colIndex<-grep("PFS|IVS|ASQ",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-3*scaled_sampledAvails[rowIndex,colIndex]
colIndex<-grep("EID|PFM",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-0.5*scaled_sampledAvails[rowIndex,colIndex]

rowIndex<-grep("BOE", ppGroups); colIndex<-grep("ZG",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-5*scaled_sampledAvails[rowIndex,colIndex]
colIndex<-grep("PFS",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-5*scaled_sampledAvails[rowIndex,colIndex]

rowIndex<-grep("BEE", ppGroups); colIndex<-grep("ASQ|BID",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-5*scaled_sampledAvails[rowIndex,colIndex]

rowIndex<-grep("ASQ", ppGroups); colIndex<-grep("PFS",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-5*scaled_sampledAvails[rowIndex,colIndex]

##check CET diets
rowIndex<-ppGroups=="CETad"
temp<-scaled_sampledAvails[rowIndex,]; index<-temp>0.1
plot(temp[index], xaxt="n"); par(las=2, mar=c(5,4,1,1))
axis(at=seq(1,length(temp[index])),labels=ppGroups[index], side=1)

#increase predation on BID adults
colIndex<-ppGroups=="EIDad"
temp<-scaled_sampledAvails[,colIndex]
scaled_sampledAvails[,colIndex]<-10*scaled_sampledAvails[,colIndex]

#############
##do some increasing of food avail for ELI
rowIndex<-grep("ELI", ppGroups); colIndex<-grep("IVSjuv|CEP|DL|DR|ASQ|BIS|SPE",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-scaled_sampledAvails[rowIndex,colIndex]+0.1

#HAK
rowIndex<-grep("HAK", ppGroups); colIndex<-grep("PFS|JAV|CEPjuv|ASQ",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-scaled_sampledAvails[rowIndex,colIndex]+0.1

#LIN
rowIndex<-grep("LINad", ppGroups); colIndex<-grep("IVS|BIS|BEE|BID|SPE|CEP|JAV",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-scaled_sampledAvails[rowIndex,colIndex]+0.1

rowIndex<-grep("LINjuv", ppGroups); colIndex<-grep("IVS|BIS|BEE|BID|SPEjuv|CEPjuv|JAV",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-scaled_sampledAvails[rowIndex,colIndex]+0.1

#ORH
rowIndex<-grep("ORHad", ppGroups); colIndex<-grep("PFS|CEPjuv|BID|EIDjuv",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-scaled_sampledAvails[rowIndex,colIndex]+0.1

rowIndex<-grep("ORHjuv", ppGroups); colIndex<-grep("PFS|CEPjuv|BIDjuv|EIDjuv",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-scaled_sampledAvails[rowIndex,colIndex]+0.1

#SND
rowIndex<-grep("SNDad", ppGroups); colIndex<-grep("PFSad|CEP|JAV|ASQ|BID|EID|BIS|EIS",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-scaled_sampledAvails[rowIndex,colIndex]+0.1

rowIndex<-grep("SNDjuv", ppGroups); colIndex<-grep("PFS|CEP|JAV|ASQ|BFF|BID|EIDjuv|BIS|EIS",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-scaled_sampledAvails[rowIndex,colIndex]+0.1

#SPD
rowIndex<-grep("SPDad", ppGroups); colIndex<-grep("CEP|DL|DR|IVS|BIS|JAV",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-scaled_sampledAvails[rowIndex,colIndex]+0.1

rowIndex<-grep("SPDjuv", ppGroups); colIndex<-grep("CEP|DL|DR|IVS|BD|BFF|BISjuv|JAVjuv",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-scaled_sampledAvails[rowIndex,colIndex]+0.1

##increase PFL pred on PFM and decrease on Z's
rowIndex<-grep("PFL", ppGroups); colIndex<-grep("^Z", ppGroups)
scaled_sampledAvails[rowIndex, colIndex]<-scaled_sampledAvails[rowIndex, colIndex]*0.1
colIndex<-grep("PFM", ppGroups)
scaled_sampledAvails[rowIndex, colIndex]<-scaled_sampledAvails[rowIndex, colIndex]+ 0.05

colIndex<-grep("ELP", ppGroups)
test<-scaled_sampledAvails[,colIndex]
index<-test[,1]>0.1
scaled_sampledAvails[!index,colIndex]<-2*scaled_sampledAvails[!index,colIndex]

colIndex<-grep("MACjuv", ppGroups)
scaled_sampledAvails[,colIndex]<-scaled_sampledAvails[,colIndex]*0.5
colIndex<-grep("MACad", ppGroups)
scaled_sampledAvails[,colIndex]<-scaled_sampledAvails[,colIndex]*1.5

colIndex<-grep("LDOjuv", ppGroups)
scaled_sampledAvails[,colIndex]<-scaled_sampledAvails[,colIndex]*0.5
colIndex<-grep("LDOad", ppGroups)
scaled_sampledAvails[,colIndex]<-scaled_sampledAvails[,colIndex]*1.5

##############
####
colIndex<-grep("RFIad", ppGroups)
test<-scaled_sampledAvails[,colIndex]
scaled_sampledAvails[,colIndex]<-5*scaled_sampledAvails[,colIndex]

### 
colIndex<-grep("EISad", ppGroups)
scaled_sampledAvails[,colIndex]<-2*scaled_sampledAvails[,colIndex]

##################
# BID
colIndex<-grep("BIDjuv", ppGroups)
test<-scaled_sampledAvails[,colIndex]
scaled_sampledAvails[,colIndex]<-scaled_sampledAvails[,colIndex]*0.5
##adults
colIndex<-grep("BIDad", ppGroups)
test<-scaled_sampledAvails[,colIndex]
scaled_sampledAvails[,colIndex]<-scaled_sampledAvails[,colIndex]*0.7


## PFS
colIndex<-grep("PFSad", ppGroups); rowIndex<-grep("HOK|JAV|PFM|PFL|ASQ|CEP",ppGroups)
test<-scaled_sampledAvails[rowIndex,colIndex]; 
scaled_sampledAvails[rowIndex, colIndex]<-2*scaled_sampledAvails[rowIndex, colIndex]
test<-scaled_sampledAvails[,colIndex]
scaled_sampledAvails[,colIndex]<-scaled_sampledAvails[,colIndex]*0.7
#juv
colIndex<-grep("PFSjuv", ppGroups)
test<-scaled_sampledAvails[,colIndex]
scaled_sampledAvails[,colIndex]<-scaled_sampledAvails[,colIndex]*0.5

#### ELI
colIndex<-grep("ELIad", ppGroups)
test<-scaled_sampledAvails[,colIndex]
index<-test>0.1 
scaled_sampledAvails[!index,colIndex]<-2*scaled_sampledAvails[!index,colIndex]
## juv
colIndex<-grep("ELIjuv", ppGroups)
test<-scaled_sampledAvails[,colIndex]
scaled_sampledAvails[,colIndex]<-0.7*scaled_sampledAvails[,colIndex]

## EID
colIndex<-grep("EID", ppGroups)
#take any values over 1 and divide by 10
fixTop<-function(x){
  y<-x
  if(x>1){
    y<-x/10
  }
  return(y)
}
scaled_sampledAvails[,colIndex]<-apply(scaled_sampledAvails[,colIndex], c(1,2), FUN=fixTop)
#now divide ad by 2
colIndex<-grep("EIDad", ppGroups)
scaled_sampledAvails[,colIndex]<-0.1*scaled_sampledAvails[,colIndex]
colIndex<-grep("EIDjuv", ppGroups)
scaled_sampledAvails[,colIndex]<-0.05*scaled_sampledAvails[,colIndex]

##
colIndex<-grep("LINjuv", ppGroups)
scaled_sampledAvails[,colIndex]<-scaled_sampledAvails[,colIndex]*0.2
colIndex<-grep("LINad", ppGroups)
scaled_sampledAvails[,colIndex]<-scaled_sampledAvails[,colIndex]*3

## DPI
colIndex<-grep("DPI", ppGroups)
scaled_sampledAvails[,colIndex]<-2*scaled_sampledAvails[,colIndex]

## CBO
colIndex<-grep("CBO", ppGroups)
scaled_sampledAvails[,colIndex]<-1.5*scaled_sampledAvails[,colIndex]

## SPE canabilism
colIndex<-grep("SPE", ppGroups)
scaled_sampledAvails[colIndex,colIndex]<-1e-2*scaled_sampledAvails[colIndex,colIndex]
# 
# colIndex<-grep("HAKjuv", ppGroups)
# test<-scaled_sampledAvails[,colIndex]
# scaled_sampledAvails[,colIndex]<-scaled_sampledAvails[,colIndex]/3
# 
# rowIndex<-grep("PFLad", ppGroups)
# test<-scaled_sampledAvails[rowIndex,]

#######
##BID #have taken off mL instead
# colIndex<-grep("BIDjuv", ppGroups)
# test<-scaled_sampledAvails[,colIndex]
# scaled_sampledAvails[,colIndex]<-0.5*scaled_sampledAvails[,colIndex]
# colIndex<-grep("BIDad", ppGroups)
# test<-scaled_sampledAvails[,colIndex]
# rowIndex<-test>0.08
# scaled_sampledAvails[rowIndex, colIndex]<-0.1*scaled_sampledAvails[rowIndex, colIndex]

#################################################################
##write them back out
this_sampledAvails<-scaled_sampledAvails

pPREYfile<-paste(DIR$'Base',"ATLANTISmodels\\inputs\\biol_prm\\pPREY\\sampledAvailsScaledForRealised.txt",sep="")
cat("## sampled availabilities take 1\n",file=pPREYfile,append=FALSE)

for(g in 1:npg){
  thisVar<-ppGroups[g]
  thisCode<-gsub('juv|ad',"",thisVar)
  thisAge<-gsub(thisCode,"",thisVar)
  thisVec<-this_sampledAvails[g,]
  if(thisAge==""){
    ppreyvar<-paste("pPREY",thisCode,"\t",ndg,"\n",sep="")
    cat(ppreyvar,file=pPREYfile,append=TRUE)
    temp<-c(signif(as.double(thisVec[not_adultIndex]),4),c(0,0,0))
    newLine<-paste(temp,collapse = " ")
    cat(newLine,file=pPREYfile,append=TRUE)
    cat("\n",file=pPREYfile,append=TRUE)
  } else{
    if(thisAge=="ad"){
      ageNumber<-2
    } else{
      ageNumber<-1
    }
    #first do not adult prey
    ppreyvar<-paste("pPREY",1,thisCode,ageNumber,"\t",ndg,"\n",sep="")
    cat(ppreyvar,file=pPREYfile,append=TRUE)
    temp<-c(signif(as.double(thisVec[not_adultIndex]),4),c(0,0,0))
    newLine<-paste(temp,collapse = " ")
    cat(newLine,file=pPREYfile,append=TRUE)
    cat("\n",file=pPREYfile,append=TRUE)
    #first do adult prey
    ppreyvar<-paste("pPREY",2,thisCode,ageNumber,"\t",ndg,"\n",sep="")
    cat(ppreyvar,file=pPREYfile,append=TRUE)
    temp<-c(signif(as.double(thisVec[adIndex]),4))
    filledTemp<-rep(0,ng); filledTemp[groupsDF$NumCohorts>1]<-temp; 
    temp<-c(filledTemp,c(0,0,0))
    temp[is.na(temp)]<-0
    newLine<-paste(temp,collapse = " ")
    cat(newLine,file=pPREYfile,append=TRUE)
    cat("\n",file=pPREYfile,append=TRUE)
  }
}

