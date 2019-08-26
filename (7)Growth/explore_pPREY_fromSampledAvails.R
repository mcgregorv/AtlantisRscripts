#read in sampled availabilities and create pPREY lines for biol.prm file
this_sampledAvails<-read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ExpectedValuesScaledForMax.csv",sep=""))
# this_sampledAvails<-read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ExpectedValues_scaledDown.csv",sep=""))
# this_sampledAvails<-read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ExpectedValues_scaledUp.csv",sep=""))

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

ppGroups<-as.character(read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ppGroups.csv",sep=""))[,1]); npg<-length(ppGroups)

groupsDF<-read.csv(paste(this_path,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
ndg<-ng+3 #number of diet groups, include sediment version of DL, DR, DC


adIndex<-grep("ad",ppGroups); juvIndex<-grep("juv",ppGroups)
naIndex<-grep("ad|juv",ppGroups,invert = TRUE)
not_adultIndex<-grep("ad",ppGroups,invert = TRUE)

##who eats DPI
testCode<-"MJE"; thisMin<-0
test<-this_sampledAvails[,grep(testCode,ppGroups)]
index<-test[,1]>thisMin
juvPreds<-ppGroups[index]
juvAvails<-test[index,1]
rbind(juvPreds,juvAvails)
index<-test[,2]>thisMin
adPreds<-ppGroups[index]; adAvails<-test[index,2]
rbind(adPreds,adAvails)

par(mfrow=c(2,1),las=2)
plot(juvAvails,type="h",lwd=3,xaxt="n",xlab="")
axis(at=seq(1,length(juvPreds)),labels=juvPreds,side=1)
plot(adAvails,type="h",lwd=3,xaxt="n",xlab="")
axis(at=seq(1,length(adPreds)),labels=adPreds,side=1)



#check out who can eat ivs juveniles
test<-this_sampledAvails[,ppGroups %in% c("IVSjuv")]
index<-test>0.005
thisPreds<-ppGroups[index]
thisAvails<-test[index]
plot(sort(thisAvails,decreasing = TRUE),type="h",lwd=3,lend=1,xaxt="n",xlab="")
par(las=2)
axis(at=seq(1,length(thisAvails)),labels=thisPreds[order(thisAvails,decreasing = TRUE)],side=1)
##adults
test<-this_sampledAvails[,ppGroups %in% c("IVSad")]
index<-test>0.005
thisPreds<-ppGroups[index]
thisAvails<-test[index]
plot(sort(thisAvails,decreasing = TRUE),type="h",lwd=3,lend=1,xaxt="n",xlab="")
par(las=2)
axis(at=seq(1,length(thisAvails)),labels=thisPreds[order(thisAvails,decreasing = TRUE)],side=1)

##same for PFS
test<-this_sampledAvails[,ppGroups %in% c("PFSjuv")]
sum(test)
index<-test>0.01
thisPreds<-ppGroups[index]
thisAvails<-test[index]
plot(sort(thisAvails,decreasing = TRUE),type="h",lwd=3,lend=1,xaxt="n",xlab="")
par(las=2)
axis(at=seq(1,length(thisAvails)),labels=thisPreds[order(thisAvails,decreasing = TRUE)],side=1)

##same for EIS
test<-this_sampledAvails[,ppGroups %in% c("EISjuv")]
index<-test>0.01
thisPreds<-ppGroups[index]
thisAvails<-test[index]
plot(sort(thisAvails,decreasing = TRUE),type="h",lwd=3,lend=1,xaxt="n",xlab="")
par(las=2)
axis(at=seq(1,length(thisAvails)),labels=thisPreds[order(thisAvails,decreasing = TRUE)],side=1)

##same for IVH
test<-this_sampledAvails[,ppGroups %in% c("IVHjuv")]
index<-test>0.01
thisPreds<-ppGroups[index]
thisAvails<-test[index]
plot(sort(thisAvails,decreasing = TRUE),type="h",lwd=3,lend=1,xaxt="n",xlab="")
par(las=2)
axis(at=seq(1,length(thisAvails)),labels=thisPreds[order(thisAvails,decreasing = TRUE)],side=1)

##what are BFF, BC and BD eating?
test<-this_sampledAvails[ppGroups %in% c("BFF"),]
index<-test>0.01
thisPreys<-ppGroups[index]
thisAvails<-test[index]
thisPreys[order(thisAvails,decreasing = TRUE)]

#who eats MAC juvs?
test<-this_sampledAvails[,ppGroups %in% c("MACjuv")]
sum(test[grep("juv|ad",ppGroups)])
index<-test>0.01
thisPreds<-ppGroups[index]
thisAvails<-test[index]
plot(sort(thisAvails,decreasing = TRUE),type="h",lwd=3,lend=1,xaxt="n",xlab="")
par(las=2)
axis(at=seq(1,length(thisAvails)),labels=thisPreds[order(thisAvails,decreasing = TRUE)],side=1)

#adult rock lobsters need more predation
test<-this_sampledAvails[,ppGroups %in% c("CRAjuv")]
sum(test[grep("juv|ad",ppGroups)])
index<-test>0.01
thisPreds<-ppGroups[index]
thisAvails<-test[index]
plot(sort(thisAvails,decreasing = TRUE),type="h",lwd=3,lend=1,xaxt="n",xlab="")
par(las=2)
axis(at=seq(1,length(thisAvails)),labels=thisPreds[order(thisAvails,decreasing = TRUE)],side=1)



##check wehich prey have too much availbility on them
test<-colSums(this_sampledAvails)
index<-test>0.5; toPlot<-test[index]
plot(toPlot,type="h",xaxt="n",xlab="")
par(las=2)
axis(at=seq(1,length(toPlot)),labels=ppGroups[index],side=1)

names(test)<-ppGroups
test[index]

test2<-rowSums(this_sampledAvails)
index2<-test2>2; toPlot<-test2[index2]
plot(toPlot,type="h",xaxt="n",xlab="")
par(las=2)
axis(at=seq(1,length(toPlot)),labels=ppGroups[index2],side=1)
