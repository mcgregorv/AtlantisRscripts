#read in sampled availabilities and create pPREY lines for biol.prm file
this_sampledAvails<-read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ExpectedValues.csv",sep=""))
# this_sampledAvails<-read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ExpectedValuesScaledForMax.csv",sep=""))

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

ppGroups<-as.character(read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ppGroups.csv",sep=""))[,1]); npg<-length(ppGroups)

groupsDF<-read.csv(paste(this_path,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
ndg<-ng+3 #number of diet groups, include sediment version of DL, DR, DC


adIndex<-grep("ad",ppGroups); juvIndex<-grep("juv",ppGroups)
naIndex<-grep("ad|juv",ppGroups,invert = TRUE)
not_adultIndex<-grep("ad",ppGroups,invert = TRUE)

##check wehich prey have too much availbility on them
test<-colSums(this_sampledAvails)
index<-test>0.5; toPlot<-test[index]
plot(toPlot,type="h",xaxt="n",xlab="")
par(las=2)
axis(at=seq(1,length(toPlot)),labels=ppGroups[index],side=1)

names(test)<-ppGroups
test[index]

test2<-rowSums(this_sampledAvails)
index2<-test2>1; toPlot<-test2[index2]
plot(toPlot,type="h",xaxt="n",xlab="")
par(las=2)
axis(at=seq(1,length(toPlot)),labels=ppGroups[index2],side=1)

##make the most any entry can be 0.03
max(this_sampledAvails)
capAtMax<-function(x,capMax=0.05){
  y<-x
  if(x>capMax){y<-capMax}
  return(y)
}
scaled_sampledAvails<-apply(this_sampledAvails,c(1,2),FUN=capAtMax)

scaleToMax<-function(x,capMax,curMax){
  y<-(x/curMax)*capMax
  return(y)
}
curMax<-max(this_sampledAvails)
scaled_sampledAvails<-apply(this_sampledAvails,c(1,2),FUN=scaleToMax,capMax=0.01,curMax=curMax)

##how much predation on Z's?
index<-ppGroups %in% c("ZS","ZM","ZL","ZG")
test<-scaled_sampledAvails[,index]
colSums(test)

##all biomass pools?
index<-grep("juv|ad",ppGroups,invert = TRUE)
test<-scaled_sampledAvails[,index]
bpSums<-colSums(test)
names(bpSums)<-ppGroups[index]
plot(bpSums,type="h",xaxt="n",xlab="")
par(las=2)
axis(at=seq(1,length(bpSums)),labels=names(bpSums),side=1)

#ZL needs to be less; ZG more
index<-ppGroups=="ZL"
foodSum<-sum(scaled_sampledAvails[index,])
test<-scaled_sampledAvails[,index]
ppGroups[test>0]
curSum<-sum(test)
newSum<-0.2
scaled_sampledAvails[,index]<-(scaled_sampledAvails[,index]/curSum)*newSum

index<-ppGroups=="ZG"
foodSum<-sum(scaled_sampledAvails[index,])
test<-scaled_sampledAvails[,index]
ppGroups[test>0]
curSum<-sum(test)
newSum<-0.4
scaled_sampledAvails[,index]<-(scaled_sampledAvails[,index]/curSum)*newSum

#check predation on PFS
index<-ppGroups=="PFSjuv"
foodSum<-sum(scaled_sampledAvails[index,])
test<-scaled_sampledAvails[,index]
bpSums<-test[test>0.001]
names(bpSums)<-ppGroups[test>0.001]
plot(bpSums,type="h",xaxt="n",xlab="")
par(las=2)
axis(at=seq(1,length(bpSums)),labels=names(bpSums),side=1)
ppGroups[test>0]
curSum<-sum(test)
# newSum<-0.02
# scaled_sampledAvails[,index]<-(scaled_sampledAvails[,index]/curSum)*newSum

index<-ppGroups=="PFSad"
foodSum<-sum(scaled_sampledAvails[index,])
test<-scaled_sampledAvails[,index]
ppGroups[test>0]
curSum<-sum(test)
newSum<-0.3
scaled_sampledAvails[,index]<-(scaled_sampledAvails[,index]/curSum)*newSum

index<-ppGroups=="IVSad"
foodSum<-sum(scaled_sampledAvails[index,])
test<-scaled_sampledAvails[,index]
ppGroups[test>0]
curSum<-sum(test)

#cap prey Avail of Zs to verts to 5e-6
newMax<-5e-6
rowIndex<-grep("juv|ad",ppGroups); colIndex<-grep("^Z",ppGroups)
curMax<-max(scaled_sampledAvails[rowIndex,colIndex],na.rm=TRUE)
scaled_sampledAvails[rowIndex,colIndex]<-(scaled_sampledAvails[rowIndex,colIndex]/curMax)*newMax

##limit canabalism to 1e-4
ads<-ppGroups[grep("ad",ppGroups)]; na<-length(ads)
thisCurMax<-0
for(a in 1:na){
  this_ad<-ads[a]
  thisRow<-ppGroups==this_ad
  thisCol<-ppGroups==gsub("ad","juv",this_ad)
  if(scaled_sampledAvails[thisRow,thisCol]>thisCurMax){thisCurMax<-as.double(scaled_sampledAvails[thisRow,thisCol])}
}
newMax<-1e-4
for(a in 1:na){
  this_ad<-ads[a]
  thisRow<-ppGroups==this_ad
  thisCol<-ppGroups==gsub("ad","juv",this_ad)
  scaled_sampledAvails[thisRow,thisCol]<-(scaled_sampledAvails[thisRow,thisCol]/thisCurMax)*newMax
}

#take out canabalism in biomass pools
# bpGroups<-ppGroups[grep("juv|ad",ppGroups,invert=TRUE)]; nbp<-length(bpGroups)
# for(i in 1:nbp){
#   thisGroup<-bpGroups[i]
#   rowIndex<-ppGroups==thisGroup
#   scaled_sampledAvails[rowIndex,rowIndex]<-0
# }

# ##scale down predation on PFSjuv
# colIndex<-ppGroups=="PFSjuv"; newSum<-0
# curSum<-sum(scaled_sampledAvails[,colIndex])
# scaled_sampledAvails[,colIndex]<-(scaled_sampledAvails[,colIndex]/curSum)*newSum
# ##scale down predation on PFSad
# colIndex<-ppGroups=="PFSad"; newSum<-0.01
# curSum<-sum(scaled_sampledAvails[,colIndex])
# scaled_sampledAvails[,colIndex]<-(scaled_sampledAvails[,colIndex]/curSum)*newSum

#reduce baleen whale predation on PFS - by increasing ZM available
rowIndex<-ppGroups=="BALad"; colIndex<-ppGroups=="ZM"
test<-scaled_sampledAvails[rowIndex,]; index<-test>0
test<-test[index]; names(test)<-ppGroups[index]
scaled_sampledAvails[rowIndex,colIndex]<-1e-3
#and juveniles
rowIndex<-ppGroups=="BALjuv"; colIndex<-ppGroups=="ZM"
test<-scaled_sampledAvails[rowIndex,]; index<-test>0
test<-test[index]; names(test)<-ppGroups[index]
scaled_sampledAvails[rowIndex,colIndex]<-1e-3

#stop ZG eating ZL
colIndex<-ppGroups=="ZL"; rowIndex<-ppGroups=="ZG"
scaled_sampledAvails[rowIndex,colIndex]<-1e-5

#stop SB eating ZS
colIndex<-ppGroups=="ZS"; rowIndex<-ppGroups %in% c("SBjuv","SBad")
scaled_sampledAvails[rowIndex,colIndex]<-0

#increase PFS predation on ZM
rowIndex<-ppGroups=="PFSad"; colIndex<-ppGroups=="ZM"
scaled_sampledAvails[rowIndex,colIndex]<-1e-4
rowIndex<-ppGroups=="PFSjuv"; colIndex<-ppGroups=="ZM"
scaled_sampledAvails[rowIndex,colIndex]<-1e-4

#edit IVS diet
thisPreys<-c("BC", "BD", "BFF", "BO", "DC", "DL", "IVHad", "IVHjuv", "IVSad", "IVSjuv", "MB"); ntp<-length(thisPreys)
thisAvails<-c(0.2, 0.2, 0.01, 0.2, 1e-05, 1e-05, 0.0045, 0.01, 0.0045, 0.01, 0.2)
rowIndex<-ppGroups=="IVSad";
for(p in 1:ntp){
  colIndex<-ppGroups==thisPreys[p]
  scaled_sampledAvails[rowIndex,colIndex]<-thisAvails[p]
}
thisPreys<-c("BD", "BO", "DC", "DL", "MA", "MB", "PL", "PS", "ZM", "ZS"); ntp<-length(thisPreys)
thisAvails<-c(0.2, 0.2, 1e-05, 1e-05, 0.2, 0.2, 1e-06, 1e-06, 1e-06, 1e-06)
rowIndex<-ppGroups=="IVSjuv";
for(p in 1:ntp){
  colIndex<-ppGroups==thisPreys[p]
  scaled_sampledAvails[rowIndex,colIndex]<-thisAvails[p]
  # cat(scaled_sampledAvails[rowIndex,colIndex],"--",thisAvails[p],"\n")
}

##based on PFS run, scale up for preds with weight loss and prey with numbers up
#preds with weight down
rowIndex<-ppGroups %in% c("BEEad", "ETBad", "BISad", "CBOad", "DPI", "ELIad", "EISad", "HAKad", "HOKad", "IVSad", "LINad", "LDOad", 
                          "MACad", "MJEad", "ORHad", "PFLad", "PFMad", "PINad", "RFIad", "CRAad", "SBad", "SPEad", "SNDad", "SSOad", 
                          "SPDad") 
#prey with nums up
colIndex<-ppGroups %in% c("BEEad", "ETBad", "BIDad", "CBOad", "DPI", "ELIjuv", "ELIad", "GSHjuv", "HAKjuv", "HOKjuv", "IVHad", "IVHad", 
                          "JAVad", "MACjuv", "MJEjuv", "PFSad")
test<-scaled_sampledAvails[rowIndex,colIndex]
scaled_sampledAvails[rowIndex,colIndex]<-2*scaled_sampledAvails[rowIndex,colIndex]

### individual edits
##PFL diet
rowIndex<-ppGroups=="PFLad"; colIndex<-ppGroups=="MACad"
scaled_sampledAvails[rowIndex,colIndex]<-0.047
colIndex<-ppGroups=="MACjuv"
scaled_sampledAvails[rowIndex,colIndex]<-0.0036
## PFM diet
rowIndex<-ppGroups=="PFMad"; colIndex<-ppGroups=="ZL"
scaled_sampledAvails[rowIndex,colIndex]<-8.4e-5
colIndex<-ppGroups=="ZM"
scaled_sampledAvails[rowIndex,colIndex]<-3.6e-5

##BEEad prey
#EID
rowIndex<-ppGroups=="BEEad"; colIndex<-ppGroups=="EIDad"
scaled_sampledAvails[rowIndex,colIndex]<-0.023
colIndex<-ppGroups=="EIDjuv"
scaled_sampledAvails[rowIndex,colIndex]<-0.023
#CEP
colIndex<-ppGroups=="CEPad"
scaled_sampledAvails[rowIndex,colIndex]<-0.032
colIndex<-ppGroups=="CEPjuv"
scaled_sampledAvails[rowIndex,colIndex]<-0.032
#BFF
colIndex<-ppGroups=="BFF"
scaled_sampledAvails[rowIndex,colIndex]<-1e-5
#PFSad
colIndex<-ppGroups=="PFSad"
scaled_sampledAvails[rowIndex,colIndex]<-7e-3

## ETB
rowIndex<-ppGroups=="ETBad"; colIndex<-grep("ORH",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-7e-3
colIndex<-ppGroups=="CEPad"
scaled_sampledAvails[rowIndex,colIndex]<-0.2
colIndex<-ppGroups=="CEPjuv"
scaled_sampledAvails[rowIndex,colIndex]<-0.05
colIndex<-grep("HOK",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-2e-3
colIndex<-grep("PFSad",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-2e-3
colIndex<-grep("BEEad",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-0.01

##BISad
#testAvail[grep("PFSad",names(testAvail))] <- 1e-3
# testAvail[grep("BFF",names(testAvail))] <- 1e-6
# testAvail[grep("ZM",names(testAvail))]<-3e-5
# testAvail[grep("ZL",names(testAvail))]<-1e-5
#PFSad
rowIndex<-ppGroups=="BISad"; colIndex<-ppGroups=="PFSad"
scaled_sampledAvails[rowIndex,colIndex]<-1e-3
#BFF
colIndex<-ppGroups=="BFF"
scaled_sampledAvails[rowIndex,colIndex]<-1e-6
#ZM
colIndex<-ppGroups=="ZM"
scaled_sampledAvails[rowIndex,colIndex]<-3e-5
#ZL
colIndex<-ppGroups=="ZL"
scaled_sampledAvails[rowIndex,colIndex]<-1e-5
#

##CBOa
#BFF
rowIndex<-ppGroups=="CBOad"; colIndex<-ppGroups=="BFF"
scaled_sampledAvails[rowIndex,colIndex]<- 1e-6
# IVSad
colIndex<-ppGroups=="IVSad"
scaled_sampledAvails[rowIndex,colIndex]<-1e-2
# IVSjuv
colIndex<-ppGroups=="IVSjuv"
scaled_sampledAvails[rowIndex,colIndex]<-0.1
#BC
colIndex<-ppGroups=="BC"
scaled_sampledAvails[rowIndex,colIndex]<-1e-3
#BD
colIndex<-ppGroups=="BD"
scaled_sampledAvails[rowIndex,colIndex]<-1e-3
#ZM
colIndex<-ppGroups=="ZM"
scaled_sampledAvails[rowIndex,colIndex]<-1e-6
#ZG
colIndex<-ppGroups=="ZG"
scaled_sampledAvails[rowIndex,colIndex]<-1e-8
#ZL
colIndex<-ppGroups=="ZL"
scaled_sampledAvails[rowIndex,colIndex]<-1e-8
#PFSad
colIndex<-ppGroups=="PFSad"
scaled_sampledAvails[rowIndex,colIndex]<-1e-8

##DPIad
rowIndex<-ppGroups=="DPIad"; 
changePreys<-c("BFF","ASQad","ASQjuv","BISad","BISjuv","HOKjuv")
changes<-c(1e-6, 1e-1, 5e-3, 5e-2, 8e-3, 1e-3); nc<-length(changePreys)
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}

##ELIad
rowIndex<-ppGroups=="ELIad"; 
changePreys<-c("BFF","IVSad","HOKjuv","HOKad","BOEjuv"); changes<-c(1e-6, 1e-1, 1e-3, 1e-4, 1e-3); nc<-length(changePreys)
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}

rowIndex<-ppGroups=="ELPad"; 
colIndex<-grep("BEEad",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-0.02

##EISad
rowIndex<-ppGroups=="EISad"
changePreys<-c("BFF","PFSad","PFSjuv","ZM"); changes<-c(1e-6, 1e-4, 1e-4, 1e-4); nc<-length(changePreys)
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}

##HAKad
changePreys<-c("BFF","PFSad","HOKjuv","HOKad", "CEPad","JAVad"); changes<-c(1e-6, 1e-4, 1e-4, 3e-3, 1e-1, 3e-2); nc<-length(changePreys)
rowIndex<-ppGroups=="HAKad"
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}

##HOKad
changePreys<-c("BFF","ZG"); changes<-c(1e-6, 1e-10); nc<-length(changePreys)
rowIndex<-ppGroups=="HOKad"
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}

##IVSad
changePreys<-c("BFF","DR"); changes<-c(4e-3, 1e-5); nc<-length(changePreys)
rowIndex<-ppGroups=="IVSad"
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}

##LINad
changePreys<-c("HOKjuv","HOKad", "ZL","PFSad", "BISad", "IVSad", "BOEad"); changes<-c(1e-3, 1e-3, 1e-5, 3e-3, 1e-2, 1e-2, 5e-4); nc<-length(changePreys)
rowIndex<-ppGroups=="LINad"
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}

##LDOad
changePreys<-c("BFF","PFSad", "EIDad","EIDjuv", "ZL", "ZM"); changes<-c(1e-6, 1e-3, 1e-2, 1e-2, 1e-5, 1e-5); nc<-length(changePreys)
rowIndex<-ppGroups=="LDOad"
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}

##MACad
changePreys<-c("BFF","ZL", "ZM"); changes<-c(1e-6, 3e-5, 1e-4); nc<-length(changePreys)
rowIndex<-ppGroups=="MACad"
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}

##MJEad
changePreys<-c("BFF","PFSad", "ZG"); changes<-c(1e-6, 3e-5, 1e-4); nc<-length(changePreys)
rowIndex<-ppGroups=="MJEad"
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}

##ORHad
changePreys<-c("BFF","CEPad", "CEPjuv", "ZL"); changes<-c(1e-6, 1e-1, 1e-2, 5e-5); nc<-length(changePreys)
rowIndex<-ppGroups=="ORHad"
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}

##PFLad
changePreys<-c("BFF","HOKad", "CEPjuv", "HOKjuv", "PFSad", "MJEad","CETjuv"); changes<-c(1e-6, 1e-4, 1e-2, 1e-3, 2e-3, 1e-3,2e-1); nc<-length(changePreys)
rowIndex<-ppGroups=="PFLad"
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}

##PFMad
changePreys<-c("BFF","HOKjuv", "ZL", "ZM"); changes<-c(1e-6, 1e-4, 1e-4, 5e-5); nc<-length(changePreys)
rowIndex<-ppGroups=="PFMad"
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}

#PINad
changePreys<-c("HOKad","HOKjuv", "PFSad", "CEPad", "ASQad", "PFMad", "MJEad", "MJEjuv", "JAVad", "EISad", "SPDad"); 
changes<-c(1e-4, 2e-4, 1e-3, 0.2, 0.2, 1e-3, 5e-4, 1e-3, 1e-3, 1e-3, 1e-3); nc<-length(changePreys)
rowIndex<-ppGroups=="PINad"
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}

##RFIad
changePreys<-c("BFF","HOKjuv", "PFSad"); changes<-c(1e-4, 1e-3, 1e-3); nc<-length(changePreys)
rowIndex<-ppGroups=="RFIad"
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}

##CRAad
changePreys<-c("BFF","IVSad", "IVHad", "BD"); changes<-c(1e-4, 4e-2, 4e-2, 1e-3); nc<-length(changePreys)
rowIndex<-ppGroups=="CRAad"
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}

##SBad
changePreys<-c("BFF","ASQad", "ZL", "DL", "DR"); changes<-c(1e-6, 0.1, 1e-4, 1e-6, 1e-5); nc<-length(changePreys)
rowIndex<-ppGroups=="SBad"
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}

##SPEad
changePreys<-c("BFF","IVSad", "BISad", "ZG", "PFSad"); changes<-c(1e-6, 1e-2, 1e-2, 1e-5, 1e-3); nc<-length(changePreys)
rowIndex<-ppGroups=="SPEad"
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}

#SSOad
changePreys<-c("BFF","ZG"); changes<-c(1e-6, 8e-4); nc<-length(changePreys)
rowIndex<-ppGroups=="SSOad"
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}

##SPDad
changePreys<-c("BFF","HOKjuv", "PFSad", "CEPad", "ASQad", "DL", "DR"); changes<-c(1e-6, 4e-3, 0.1, 1e-3, 0.1, 1e-6, 1e-5); nc<-length(changePreys)
rowIndex<-ppGroups=="SPDad"
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}

##fix BFF predators - seem to want to be eaten by evertying! restrict here to what should actually be eating them
BFFpreds<-c("CBOad", "SPEad", "DPIad", "EIDad", "EIDjuv", "RFIad", "GSHad", "SPDjuv", "SNDjuv", "ELPjuv", 
            "CEPjuv", "CRAad", "CRAad", "IVSjuv", "IVSad")
rowIndex<-ppGroups %in% BFFpreds; colIndex<-ppGroups=="BFF"
scaled_sampledAvails[!rowIndex,colIndex]<-0

###### reduce food for ZG
rowIndex<-ppGroups=="ZG"
# test<-scaled_sampledAvails[rowIndex,]
# test[test>0]
scaled_sampledAvails[rowIndex,]<-scaled_sampledAvails[rowIndex,]/100

groupsStillHungry<-c("LDO", "MAC", "MJE", "PFL", "PFM", "RFI", "CRA", "BAL", "ELI", "HAK", "HOK", "LIN"); ngh<-length(groupsStillHungry)
for(i in 1:ngh){
  rowIndex<-grep(groupsStillHungry[i], ppGroups)
  scaled_sampledAvails[rowIndex,]<-2*scaled_sampledAvails[rowIndex,]
}

##increase avail of all age-structured
colIndex<-grep("juv|ad",ppGroups)
scaled_sampledAvails[,colIndex]<-2*scaled_sampledAvails[,colIndex]

## do even more! But skip a couple off the menu that are looking good for number
skipPreys<-c("BID", "EIS", "JAV", "PFS"); temp<-ppGroups[grep("juv|ad", ppGroups)]; eatMorePrey<-temp[!(temp %in% skipPreys)]
colIndex<-ppGroups %in% eatMorePrey
scaled_sampledAvails[,colIndex]<-2*scaled_sampledAvails[,colIndex]


##then scale for max over 0.8
curMax<-max(scaled_sampledAvails[,colIndex],na.rm=TRUE)
scaled_sampledAvails[,colIndex]<-apply(scaled_sampledAvails[,colIndex],c(1,2),FUN=function(x){min(x,0.8)})

#edit ivh food a little
rowIndex<-grep("IVH",ppGroups); colIndex<-ppGroups %in% c("PL", "MB", "MA")
scaled_sampledAvails[rowIndex,colIndex]<-scaled_sampledAvails[rowIndex,colIndex]*4

##increase avail of IVS to its age-structured preds
colIndex<-grep("IVS",ppGroups)
rowIndex<-grep("ad|juv",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-2*scaled_sampledAvails[rowIndex,colIndex]
scaled_sampledAvails[rowIndex,colIndex]<-apply(scaled_sampledAvails[rowIndex,colIndex],c(1,2),FUN=function(x){min(x,0.8)})

#increase PFL predation on age-structured prey
rowIndex<-grep("PFL",ppGroups)
colIndex<-grep("ad|juv",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-2*scaled_sampledAvails[rowIndex,colIndex]
scaled_sampledAvails[rowIndex,colIndex]<-apply(scaled_sampledAvails[rowIndex,colIndex],c(1,2),FUN=function(x){min(x,0.8)})

#increase PFM predation on age-structured prey, ZL and ZM
rowIndex<-grep("PFM",ppGroups)
colIndex<-grep("ad|juv|ZM|ZL",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-2*scaled_sampledAvails[rowIndex,colIndex]

#increase PIN predation on age-structured prey
rowIndex<-grep("PIN",ppGroups)
colIndex<-grep("ad|juv",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-4*scaled_sampledAvails[rowIndex,colIndex]
scaled_sampledAvails[rowIndex,colIndex]<-apply(scaled_sampledAvails[rowIndex,colIndex],c(1,2),FUN=function(x){min(x,0.8)})

#increase SSO predation on age-structured prey
rowIndex<-grep("SSO",ppGroups)
colIndex<-grep("ad|juv",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-4*scaled_sampledAvails[rowIndex,colIndex]

rowIndex<-grep("HAK|HOK|LIN|MAC|MJE|PFL|PFM|PIN|RFI",ppGroups)
scaled_sampledAvails[rowIndex,]<-2*scaled_sampledAvails[rowIndex,]
scaled_sampledAvails[rowIndex,]<-apply(scaled_sampledAvails[rowIndex,],c(1,2),FUN=function(x){min(x,0.8)})

tooMany<-"CRA|LIN|LDO|MAC|MJE|PFL"
colIndex<-grep(tooMany,ppGroups)
scaled_sampledAvails[,colIndex]<-4*scaled_sampledAvails[,colIndex]
scaled_sampledAvails[,colIndex]<-apply(scaled_sampledAvails[,colIndex],c(1,2),FUN=function(x){min(x,0.8)})

##increase food for CETad
rowIndex<-ppGroups=="CETad"
scaled_sampledAvails[rowIndex,]<-5*scaled_sampledAvails[rowIndex,]
changePreys<-c("ASQad","CEPad"); changes<-c(0.4,0.4); nc<-length(changePreys)
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}
#CET drink their mother's milk, but this isnt' set up in Atlantis, so they need to eat the same food as the adults instead. 
juvRowIndex<-ppGroups=="CETjuv"
scaled_sampledAvails[juvRowIndex,]<-scaled_sampledAvails[rowIndex,]

##same for PIN
#PIN drink their mother's milk, but this isnt' set up in Atlantis, so they need to eat the same food as the adults instead. 
rowIndex<-ppGroups=="PINad";  juvRowIndex<-ppGroups=="PINjuv"
scaled_sampledAvails[juvRowIndex,]<-scaled_sampledAvails[rowIndex,]

##same for BAL
#BAL drink their mother's milk, but this isnt' set up in Atlantis, so they need to eat the same food as the adults instead. 
rowIndex<-ppGroups=="BALad";  juvRowIndex<-ppGroups=="BALjuv"
scaled_sampledAvails[juvRowIndex,]<-scaled_sampledAvails[rowIndex,]


##more ZL for baleen whales
rowIndex<-ppGroups %in% c("BALad", "BALjuv"); colIndex<-ppGroups %in% c("ZL", "ZM")
scaled_sampledAvails[rowIndex,colIndex]<-0.1

##feed more to BIS
rowIndex<-ppGroups %in% c("BISad", "BISjuv"); 
scaled_sampledAvails[rowIndex,]<-5*scaled_sampledAvails[rowIndex,]

#more ZM for EISjuv
rowIndex<-ppGroups=="EISjuv"; colIndex<-ppGroups=="ZM"
scaled_sampledAvails[rowIndex,colIndex]<-1e-4

#More MA and MB for IVH
changePreys<-c("MA","MB"); changes<-c(0.4,0.4); nc<-length(changePreys)
rowIndex<-ppGroups %in% c("IVHad", "IVHjuv")
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}

#increase food for LIN
rowIndex<-grep("LINjuv", ppGroups)
scaled_sampledAvails[rowIndex,]<-2*scaled_sampledAvails[rowIndex,]

#increase food for PFLjuv
rowIndex<-ppGroups=="PFLjuv"; colIndex<-grep("PFS",ppGroups,invert = TRUE)
scaled_sampledAvails[rowIndex,colIndex]<-2*scaled_sampledAvails[rowIndex,colIndex]
scaled_sampledAvails[rowIndex,ppGroups=="ZM"]<-3e-3
scaled_sampledAvails[rowIndex,ppGroups=="ZL"]<-4e-3

##eat more CRA adults
colIndex<-ppGroups=="CRAad"
scaled_sampledAvails[,colIndex]<-2*scaled_sampledAvails[,colIndex]

# needMoreFood<-grep("BAL|LDO|LIN|MJE|PFS|PIN|RFI|CET|CRA|SB|SND|SSO|SPD")

needLessMort<-grep("BID|BIS|EIS|JAV|PFS",ppGroups)
test<-scaled_sampledAvails[,needLessMort]
scaled_sampledAvails[,needLessMort]<-scaled_sampledAvails[,needLessMort]*0.1

# more ZM and ZL for ZG
rowIndex<-ppGroups=="ZG"; colIndex<-ppGroups %in% c("ZM", "ZL")
scaled_sampledAvails[rowIndex,colIndex]<-1e-3
##increase other avails for ZG as well
scaled_sampledAvails[rowIndex,]<-5*scaled_sampledAvails[rowIndex,]


rowIndex<-ppGroups=="BEEad";
changePreys<-c("BFF","EIDad", "DL", "DR", "CEPad");
changes<-c(1e-6, 0.5, 1e-6, 1e-5, 0.1); nc<-length(changePreys)
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}

rowIndex<-ppGroups=="GSHad"
changePreys<-c("BFF","IVSad");
changes<-c(1e-4,0.2); nc<-length(changePreys)
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}

##decrease CET predation on PFS
rowIndex<-ppGroups=="CETad"
changePreys<-c("PFSad","PFSjuv"); changes<-c(1e-3, 1e-4); nc<-length(changePreys)
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}
#CET drink their mother's milk, but this isnt' set up in Atlantis, so they need to eat the same food as the adults instead. 
juvRowIndex<-ppGroups=="CETjuv"
scaled_sampledAvails[juvRowIndex,]<-scaled_sampledAvails[rowIndex,]

#reduce SPD predation on PFS
rowIndex<-ppGroups=="SPDad"; colIndex<-ppGroups=="PFSad"
scaled_sampledAvails[rowIndex,colIndex]<-1e-3

#increase hok juv pred on ZL, ZM
rowIndex<-ppGroups=="HOKjuv"
changePreys<-c("ZL","ZM");
changes<-c(3e-4,1e-3); nc<-length(changePreys)
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}
#same for Jav ad and juv
#increase hok juv pred on ZL, ZM
rowIndex<-ppGroups %in% c("JAVjuv", "JAVad")
changePreys<-c("ZL","ZM");
changes<-c(3e-4,1e-3); nc<-length(changePreys)
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}

#LIN ad to eat DL, DR
rowIndex<-ppGroups=="LINad"; colIndex<-ppGroups=="DR"
scaled_sampledAvails[rowIndex,colIndex]<- 5e-4
colIndex<-ppGroups=="DL"
scaled_sampledAvails[rowIndex,colIndex]<- 5e-5

#more ZL, ZM for LINjuv to eat
rowIndex<-ppGroups=="LINjuv"; colIndex<-ppGroups=="ZM"
scaled_sampledAvails[rowIndex,colIndex]<-1e-3
colIndex<-ppGroups=="ZL"
scaled_sampledAvails[rowIndex,colIndex]<-3e-4

#reduce CET and PIN pred on BIS and BAL to zero
rowIndex<-grep("CET|PIN",ppGroups)
colIndex<-grep("BIS",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-scaled_sampledAvails[rowIndex,colIndex]*1e-4
rowIndex<-grep("BAL",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-0

#edit LDO juv zoo diet
rowIndex<-ppGroups=="LDOjuv"
changePreys<-c("ZL","ZM","ZG");
changes<-c(1e-5, 1e-5, 1e-6); nc<-length(changePreys)
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}

#mac ad
rowIndex<-ppGroups=="MACad"
changePreys<-c("ZL","ZM","HOK","MACad");
changes<-c(3e-4, 5e-4, 0, 1e-4); nc<-length(changePreys)
for(c in 1:nc){
  scaled_sampledAvails[rowIndex,ppGroups==changePreys[c]]<-changes[c]
}

#stop BFF eating ZG
rowIndex<-ppGroups=="BFF"; colIndex<-ppGroups=="ZG"
scaled_sampledAvails[rowIndex,colIndex]<-0

#reduce CET pred on IVS and BID
rowIndex<-grep("CET",ppGroups); colIndex<-grep("IVS|BID",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-scaled_sampledAvails[rowIndex,colIndex]*1e-3

##increase pred on Hak and ling
colIndex<-grep("HAK|LIN",ppGroups)
scaled_sampledAvails[,colIndex]<-scaled_sampledAvails[,colIndex]*10

#MJE and SSO to eat lots more ZG
rowIndex<-ppGroups %in% c("SSOad", "MJEad")
colIndex<-ppGroups=="ZG"
scaled_sampledAvails[rowIndex,colIndex]<-0.1

##reduce predation on BIDad
colIndex<-ppGroups=="BIDad"
test<-scaled_sampledAvails[,colIndex]
scaled_sampledAvails[,colIndex]<-scaled_sampledAvails[,colIndex]/10

#reduce CET pred on LDO, LIN, MAC
colIndex<-ppGroups %in% c("LDOad", "LINad", "MACad", "RFIad")
rowIndex<-grep("CET",ppGroups); 
test<-scaled_sampledAvails[rowIndex,colIndex]
scaled_sampledAvails[rowIndex,colIndex]<-0.001
colIndex<-ppGroups %in% c("LINjuv")
scaled_sampledAvails[rowIndex,colIndex]<-0.001

colIndex<-ppGroups %in% c("LDOad", "LINad", "MACad", "RFIad")
rowIndex<-grep("PIN",ppGroups); 
test<-scaled_sampledAvails[rowIndex,colIndex]
scaled_sampledAvails[rowIndex,colIndex]<-0.02
colIndex<-ppGroups %in% c("LINjuv", "MACjuv")
test<-scaled_sampledAvails[rowIndex,colIndex]
scaled_sampledAvails[rowIndex,colIndex]<-0.02

colIndex<-ppGroups %in% c("LDOad", "LINad", "MACad", "RFIad")
rowIndex<-grep("PFL",ppGroups); 
test<-scaled_sampledAvails[rowIndex,colIndex]
scaled_sampledAvails[rowIndex,colIndex]<-0.02

#are RFI eating themselves?
rowIndex<-grep("RFIad", ppGroups); colIndex<-grep("RFI",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-1e-4



##increasae avail of PFL, HOK, HAK to PIN. note PFLad had high avail already but they are outside gapelimits
colIndex<-grep("HOKad|HAKad|PFLad",ppGroups)
rowIndex<-grep("CET",ppGroups)
test<-scaled_sampledAvails[rowIndex,colIndex]
scaled_sampledAvails[rowIndex,colIndex]<-.12
#add some juv too
colIndex<-grep("HOKjuv|HAKjuv|PFLjuv",ppGroups)
test<-scaled_sampledAvails[rowIndex,colIndex]
scaled_sampledAvails[rowIndex,colIndex]<-0.025

#increase CRA CEP, ELI and RFI pred on IVH
colIndex<-grep("IVH",ppGroups); rowIndex<-grep("CRA|CEP|IVS|RFI",ppGroups)
test<-scaled_sampledAvails[rowIndex,colIndex]
scaled_sampledAvails[rowIndex,colIndex]<-apply(scaled_sampledAvails[rowIndex,colIndex],c(1,2),FUN=function(x){max(x,0.05)})

#decrease pred on IVS
colIndex<-grep("IVS",ppGroups)
rowIndex<-grep("CBO|CRA|GSH|IVS|BID|SPE",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]
rowIndex<-ppGroups=="CBOad"; colIndex<-grep("IVS",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-0.04
rowIndex<-ppGroups=="CRAad"; colIndex<-grep("IVSad",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-0.04
rowIndex<-grep("IVSad",ppGroups); grep("IVS",ppGroups)
scaled_sampledAvails[rowIndex,colIndex]<-1e-5



#take canabalism off cEt
colIndex<-grep("CEP",ppGroups)
scaled_sampledAvails[colIndex,colIndex]<-0

#reduce PFS avails again
colIndex<-grep("PFS", ppGroups)
test<-scaled_sampledAvails[,colIndex]
scaled_sampledAvails[,colIndex]<-scaled_sampledAvails[,colIndex]*5e-2

###############################################

this_sampledAvails<-scaled_sampledAvails

pPREYfile<-paste(DIR$'Base',"ATLANTISmodels\\inputs\\biol_prm\\pPREY\\sampledAvailsScaledForMax.txt",sep="")
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

#write scaled avails out too
write.csv(this_sampledAvails,paste(DIR$'Tables',"SampledAvailabilities\\ExpectedValuesScaledForMax.csv",sep=""),row.names = TRUE)
          
















###
# 
# ##leave BC, BD, BFF, BO as are as they are eating and being eaten a lot and not tracked as individuals (they are biomass pools)
# ##ZG should be 0.02; ZL 0.35; ZM 1.4; ZS 2 (approx; based on daily growth rates)
# ##PL 1.5; PS 1.1; DF 0.5; MA 0.5; MB 0.35
# # 
# # # scaleDownVars<-c("BEEjuv", "BISjuv", "CBOjuv", "CRAjuv", "EISad",  "EISjuv", "IVSad",  "JAVad",  "JAVjuv", "LDOjuv", "MACjuv", "PFSad",  "PFSjuv", "RFIjuv",  "SPEjuv")
# # scaleDownVars<-ppGroups[index][!(ppGroups[index]%in%c("BO","BC","BD","BFF","BB","DC","DR","DL","PB"))]
# # scales<-rep(0.5,length(scaleDownVars))
# # scales[scaleDownVars=="ZG"]<-0.02; scales[scaleDownVars=="ZL"]<-0.35; scales[scaleDownVars=="ZM"]<-1.4; scales[scaleDownVars=="ZS"]<-0.02;  
# # scales[scaleDownVars=="MB"]<-0.35; scales[scaleDownVars=="PL"]<-1.5; scales[scaleDownVars=="PS"]<-1.1;   
# 
# ##some groups seem to cope with high availabilities - perhaps protected through gape size or spatilly or v. productive (like squid)
# #scale down these ones
# # scaled_sampledAvails<-this_sampledAvails
# # # newMax<-0.1
# # for(s in 1:length(scaleDownVars)){
# #   thisVar<-scaleDownVars[s]
# #   thisVec<-this_sampledAvails[,grep(thisVar,ppGroups)]; thisSum<-sum(thisVec,na.rm=TRUE)
# #   newVec<-(thisVec/thisSum)*scales[s]
# #   scaled_sampledAvails[,grep(thisVar,ppGroups)]<-newVec
# # }
# # 
# # #put zero predation on some that are still declining to test
# # scaleDownVars<-c("CRAjuv", "EISad",  "EISjuv", "IVSad",  "JAVad",  "MACjuv", "PFSad",  "PFSjuv")
# # newMax<-0
# # for(s in 1:length(scaleDownVars)){
# #   thisVar<-scaleDownVars[s]
# #   thisVec<-this_sampledAvails[,grep(thisVar,ppGroups)]; thisSum<-sum(thisVec,na.rm=TRUE)
# #   newVec<-(thisVec/thisSum)*newMax
# #   scaled_sampledAvails[,grep(thisVar,ppGroups)]<-newVec
# # }
# 
# scaled_sampledAvails<-sampledAvails
# ## any additional changes..?
# #reduce predation of BC and BFF on IVS (max 0.001)
# colIndex<-ppGroups=="IVSjuv"; 
# test<-scaled_sampledAvails[,colIndex]
# index<-test>0; toCheck<-test[index]; toCheckP<-ppGroups[index]; toCheckP[order(toCheck, decreasing=TRUE)]; sort(toCheck,decreasing = TRUE)
# rowIndex<-ppGroups %in% c("BC","BFF")
# scaled_sampledAvails[rowIndex,colIndex]<-0.0001
# 
# #reduce predation of BC and BFF on BEE (max 0.001)
# colIndex<-ppGroups=="BEEjuv"; rowIndex<-ppGroups %in% c("BC","BFF")
# scaled_sampledAvails[rowIndex,colIndex]<-0.001
# 
# colIndex<-ppGroups=="IVSad"; rowIndex<-ppGroups %in% c("BC","BFF")
# scaled_sampledAvails[rowIndex,colIndex]<-0.001
# 
# ##similar for BC and BFF on PFS juv
# colIndex<-ppGroups=="PFSjuv"; rowIndex<-ppGroups %in% c("BC","BFF")
# scaled_sampledAvails[rowIndex,colIndex]<-0.001
# 
# ##similar for BC and BFF on EIS juv
# colIndex<-ppGroups=="EISjuv"; rowIndex<-ppGroups %in% c("BC","BFF")
# scaled_sampledAvails[rowIndex,colIndex]<-0.001
# 
# #reduce predation of BC and BFF on IVH (max 0.001)
# colIndex<-ppGroups %in% c("IVHad","IVHjuv"); rowIndex<-ppGroups %in% c("BC","BD","BFF")
# scaled_sampledAvails[rowIndex,colIndex]<-0.001
# 
# #MAC juvs could do with some more predation, but not from biomass pools
# rowIndex<-grep("juv|ad",ppGroups)
# colIndex<-ppGroups=="MACjuv"
# curSum<-sum(scaled_sampledAvails[rowIndex,colIndex])
# scaled_sampledAvails[rowIndex,colIndex]<-(scaled_sampledAvails[rowIndex,colIndex]/curSum)*0.7
# 
# colsToIncrease<-c("SNDad","SSOad","SPDad","BALjuv","BEEjuv","BIDjuv","BOEad","CBOjuv","EIDad","GSHjuv",
#                   "CRAad","HAKjuv","HAKad","HOKad","HOKjuv","JAVad","LINad","LINad","LDOad","MACad","MACjuv","MJEad","ORHad")
# 
# # newSum<-0.7
# for(i in 1:length(colsToIncrease)){
#   thisSSA<-scaled_sampledAvails[,ppGroups==colsToIncrease[i]]
#   thisSum<-sum(thisSSA,na.rm=TRUE)
#   # scaled_sampledAvails[,ppGroups==colsToIncrease[i]]<-(scaled_sampledAvails[,ppGroups==colsToIncrease[i]]/thisSum)*newSum
#   cat(paste(colsToIncrease[i],signif(thisSum,2),"--"))
#   
# }
# ##mac adults only have 0.045 predation on them
# thisSSA<-scaled_sampledAvails[,ppGroups=="MACad"]
# thisSum<-sum(thisSSA,na.rm=TRUE)
# scaled_sampledAvails[,ppGroups=="MACad"]<-(scaled_sampledAvails[,ppGroups=="MACad"]/thisSum)*0.25
# 
# # SNDad 0.5 --SSOad 0.5 --SPDad 0.5 --BALjuv 0 --BEEjuv 0.5 --BIDjuv 0.5 --BOEad 0.5 --CBOjuv 0.5 --EIDad 0.5 --GSHjuv 0.35 
# # --CRAad 0.05 --HAKjuv 0.26 --HAKad 0.5 --HOKad 0.5 --HOKjuv 0.36 --JAVad 0.11 --LINad 0.5 --LINad 0.5 --LDOad 0.5 --MACad 0.2 
# # --MACjuv 0.7 --MJEad 0.35 --ORHad 0.5 --
# 
# #CRA adults could do with some more predation, but not from biomass pools
# #first, add to CEP, ELP and RFI, then scale up 
# rowIndex<-grep("CEP|ELP|RFI",ppGroups)
# colIndex<-ppGroups=="CRAad"
# scaled_sampledAvails[rowIndex,colIndex]<-scaled_sampledAvails[rowIndex,colIndex]+0.01
# rowIndex<-grep("juv|ad",ppGroups)
# curSum<-sum(scaled_sampledAvails[rowIndex,colIndex])
# scaled_sampledAvails[rowIndex,colIndex]<-(scaled_sampledAvails[rowIndex,colIndex]/curSum)*0.9
# 
# #zero predation on PFS, EIS, IVS juv from biomass pools
# colIndex<-ppGroups == "IVSad"
# index<-scaled_sampledAvails[,colIndex]>0
# ppGroups[index]
# plot(scaled_sampledAvails[index,colIndex],type="h",lwd=5,lend=1,xaxt="n")
# axis(at=seq(1,length(ppGroups[index])),labels=ppGroups[index],side=1)
# #reduce BO predation to 0.003
# rowIndex<-ppGroups=="BO"
# scaled_sampledAvails[rowIndex,colIndex]<-0.003
# #check total predation on adults
# curSum<-sum(scaled_sampledAvails[,colIndex])
# #make equal to 0.05
# scaled_sampledAvails[,colIndex]<-(scaled_sampledAvails[,colIndex]/curSum)*0.25
# 
# # colIndex<-ppGroups %in% c("PFSjuv","EISjuv","IVSjuv")
# # rowIndex<-grep("juv|ad",ppGroups,invert = TRUE)
# # scaled_sampledAvails[rowIndex,colIndex]<-0.01*scaled_sampledAvails[rowIndex,colIndex]
# # 
# # #reduce from other predators too
# # colIndex<-ppGroups %in% c("PFSjuv","EISjuv","IVSjuv")[1]
# # curSum<-sum(scaled_sampledAvails[,colIndex])
# # scaled_sampledAvails[,colIndex]<-(scaled_sampledAvails[,colIndex]/curSum)*0.1
# # 
# # colIndex<-ppGroups %in% c("PFSjuv","EISjuv","IVSjuv")[2]
# # curSum<-sum(scaled_sampledAvails[,colIndex])
# # scaled_sampledAvails[,colIndex]<-(scaled_sampledAvails[,colIndex]/curSum)*0.1
# # 
# # 
# # colIndex<-ppGroups %in% c("PFSjuv","EISjuv","IVSjuv")[3]
# # curSum<-sum(scaled_sampledAvails[,colIndex])
# # scaled_sampledAvails[,colIndex]<-(scaled_sampledAvails[,colIndex]/curSum)*0.1
