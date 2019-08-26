# some of the groups have a high proportion or mortality from mL - how much of the total biomass?
plotPath <- paste(DIR$'Base',"Reports\\Chaos\\Figures\\", sep="")
basePath <-paste(DIR$'Base',"AtlantisModels\\chaos\\",sep="")

## biomass by mortality proportions
#read in B0 - this should be OK to use
ppGroups<-as.character(read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ppGroups.csv",sep=""))[,1]);
npg<-length(ppGroups)
# read in biomassByppGroup 
load(file=paste(DIR$'Base',"\\ATLANTISmodels\\ArchivedModels\\Base\\biomassByppGroup", sep=""))

#read in all the rankings - TL, lifespan, keystoneness,...
allTheRankings <- read.csv(paste(DIR$'Tables', "allTheRankings.csv", sep="")); nrankings <- dim(allTheRankings)[2]
groupsDF<- read.csv(paste(basePath,"..\\CRAM_groups.csv", sep="")); ng<-dim(groupsDF)[1]

## need numbers by group and ad/juv so can get weighted average for Mprops
load(paste(basePath,"storeNumbersByGroup",sep="")) #brings in storeNumbersByGroup, created in mortalitySummary.R in (10)Mortality\
# need to sum over adult and juvenile for each
# read in biol.prm file to get age at maturity for ad/juv split
biolLines <- readLines(paste(basePath, "CRAM_base_biol.prm", sep=""))
storeStageNumbers <- array(NA, dim=c(116, 55, 2))
for(g in 1:ng){
  thisCode <- as.character(groupsDF$Code[g]); thisNumCohorts <- groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    # if(thisNumCohorts > 2){
      #get age mature
      this_age_mat <- get_age_mat(thisCode)
      # juv
      temp<-storeNumbersByGroup[,g,1:this_age_mat]
      if(length(dim(temp))>1){
        thisNums <- apply(temp, 1, sum, na.rm=TRUE)
      } else{
        thisNums <- temp
      }
      storeStageNumbers[,g,1] <- thisNums
      # adults
      temp<-storeNumbersByGroup[,g,(this_age_mat+1):thisNumCohorts]
      if(length(dim(temp))>1){
        thisNums <- apply(temp, 1, sum, na.rm=TRUE)
      } else{
        thisNums <- temp
      }
      storeStageNumbers[,g,2] <- thisNums
      
    # }
  }
}

# set up data frame that can summarise - 2 columns, 1 for proportion mortality, 1 for biomass
# first will put biomass pool groups as NA - then will check on them too
biomassBymL <- data.frame(matrix(NA, ncol=3, nrow=npg)); colnames(biomassBymL) <- c("Code", "Mprop", "Biomass")
for(g in 1:npg){
  thisCode <- gsub("ad|juv","", ppGroups[g])
  thisAge <- gsub(thisCode , "", ppGroups[g])
  groupIndex <- groupsDF$Code==thisCode
  thisMprop <- NA
  if(thisAge=="ad"){
    thisMprop <- allTheRankings$propAdM[groupIndex]
  } else if(thisAge=="juv"){
    thisMprop <- allTheRankings$propJuvM[groupIndex]
  }
  biomassBymL$Mprop[g] <- thisMprop
  biomassBymL$Biomass[g] <- biomassByppGroup[g,1]*(mg_2_tonne*X_CN)
  biomassBymL$Code[g]<-ppGroups[g]
}
# also do m prop weighted average by numbers of adults/juveniles
MpropByGroup <- data.frame(matrix(NA, ncol=6, nrow=ng)); colnames(MpropByGroup)<-c("Code","NumJuv", "NumAd", "MpropJuv", "MpropAd", "MpropAve")
for(g in 1:ng){
  thisCode <- as.character(groupsDF$Code[g])
  thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts){
    juvNums <- median(storeStageNumbers[,g,1], na.rm=TRUE)
    adNums <- median(storeStageNumbers[,g,2], na.rm=TRUE)
    juvProp <- biomassBymL$Mprop[biomassBymL$Code==paste(thisCode,"juv", sep="")]
    adProp <- biomassBymL$Mprop[biomassBymL$Code==paste(thisCode,"ad", sep="")]
    aveMprop <- (juvNums * juvProp + adNums * adProp)/(adNums + juvNums)
    MpropByGroup[g,c("Code","NumJuv", "NumAd", "MpropJuv", "MpropAd", "MpropAve")]<-c(thisCode, juvNums, adNums, juvProp, adProp, aveMprop)
  }
  
}
as_MpropByGroup <- MpropByGroup[!is.na(MpropByGroup$MpropAve),]
orderMprop <- order(as.double(as_MpropByGroup$MpropAve), decreasing = TRUE)
as_MpropByGroup$Code[orderMprop]
thisSpeciesNames <- gsub("_", " ", str_trim(groupsDF$Name[match(as_MpropByGroup$Code[orderMprop], groupsDF$Code)], side="both"))

pdf(paste(plotPath, "MpropByGroup.pdf", sep=""), height=4, width=7)
par(lend=1, mar=c(10,4,1,1))
plot(as_MpropByGroup$MpropAve[orderMprop], type="h", lwd=5, xaxt="n", xlab="", col=myGrey, ylab="Proporion of M forced")
par(las=2)
axis(at=1:length(MpropByGroup$Code[orderMprop]), labels = thisSpeciesNames, side=1)
abline(h=seq(0.2,0.8,by=0.2), col=myGrey_trans, lwd=1, lty=2)
dev.off()


pdf(paste(plotPath, "MpropByGroup_withAge.pdf", sep=""), height=5, width=10)
par(lend=1, mar=c(10,4,0.1,8), las=1)
plot(as_MpropByGroup$MpropAve[orderMprop], type="h", lwd=5, xaxt="n", xlab="", col=myGrey, ylab="Proporion of M forced", ylim=c(0, 1.1))
par(las=2)
axis(at=1:length(MpropByGroup$Code[orderMprop]), labels = thisSpeciesNames, side=1)
abline(h=seq(0.2,1,by=0.2), col=myGrey_trans, lwd=1, lty=2)
points(as_MpropByGroup$MpropAd[orderMprop], pch=3, col=myBlue, lwd=2)
points(as_MpropByGroup$MpropJuv[orderMprop], pch=4,col=myOrange, lwd=2)
par(xpd=NA)
legend(legend = c("Combined", "Adult", "Juvenile"), col=c(myGrey, myBlue,myOrange), pch=c(NA, 3, 4), lty=c(1,NA, NA), lwd=c(5,2,2),x="right", inset=-0.2, bty="n")
dev.off()

## hist of adult propM and juv prop M
plot(x=c(1,2), y=c(0,1), type="n", xlim=c(0.5,2.5), ylim=c(0, 1.1))
boxplot(as.double(MpropByGroup$MpropJuv), at=1, col=myOrange_trans, border=myOrange, add=TRUE)
boxplot(as.double(MpropByGroup$MpropAd), at=2, col=myBlue_trans, border=myBlue, add=TRUE)

thisBreaks<-seq(0, 1.2, by=0.15)
h1 <- hist(as.double(MpropByGroup$MpropJuv), plot=FALSE, breaks=thisBreaks)
h2 <- hist(as.double(MpropByGroup$MpropAd), plot=FALSE, breaks=thisBreaks)
plot(x=h1$mids, y=h1$counts, type="h", lwd=5, col=myOrange_trans, xlim=c(0, 1.1), ylim=c(0, 17.5))
points(x=h2$mids, y=h2$counts, type="h", lwd=5, col=myBlue_trans)

biomassBymL$roundedMprop <- unlist(lapply(biomassBymL$Mprop, myRounding, fraction=0.1))
hist(biomassBymL$roundedMprop)

toPlot<-data.frame(cbind("mid"=c(h1$mids, h2$mids), "age"=c(rep("Juvenile",length(h1$mids)), rep("Adult", length(h2$mids))), "value"=c(h1$counts, h2$counts)))
toPlot$age <- factor(toPlot$age, levels=c("Juvenile", "Adult"))
bp<-ggplot(data = toPlot, aes(x = mid, fill = age, y = value)) + 
  geom_bar(stat = 'identity')
pdf(paste(plotPath,"PropForcedMbyAdJuv.pdf", sep=""), height=4, width=10)
bp + scale_fill_manual(values=c(myOrange, myBlue)) + labs(y="Number of groups", x="Proportion of forced mortality") + theme_igray() + 
  theme(axis.text=element_text(size=12, angle=0),axis.title=element_text(size=12)) + guides(fill=guide_legend(title="")) 
dev.off()


toPlot <- tapply(biomassBymL$Biomass, biomassBymL$roundedMprop, sum, na.rm=TRUE)

plot(x=names(toPlot), y=toPlot, type="h", lwd=10)

# do not include detritus and bacteria in total
thisIndex <- grep("DC|DR|DL|BB|PB", ppGroups, invert=TRUE)
thisIndex <- grep("DC", ppGroups, invert=TRUE)
allBtotal <- sum(biomassBymL$Biomass[thisIndex], na.rm=TRUE)

par(lend=1, las=0)
thisRoundings <- sort(unique(biomassBymL$roundedMprop))
nr <- length(thisRoundings); storeBiomPerc <- rep(NA, nr)
for(r in 1:nr){
  thisR <- thisRoundings[r]
  if(is.na(thisR)){
    thisIndex <- is.na(biomassBymL$roundedMprop)
  } else{
    thisIndex <- biomassBymL$roundedMprop==thisR & !is.na(biomassBymL$roundedMprop)
  }
  testB <- biomassBymL$Biomass[thisIndex]
  testCodes <- gsub("ad|juv", "", ppGroups[thisIndex])
  thisBtotal <- sum(testB); 
  storeBiomPerc[r]<-thisBtotal/allBtotal
  # plot(testB, type="h", lwd=10, col=myBlue, xaxt="n", xlab="", ylab="")
  # mtext("Biomass (tonnes)", side=2, line=3)
  # mtext(paste("Percentage of total biomass: ", signif(thisBtotal/allBtotal,2)*100, "%", sep=""), side=3, adj=0)
}
# get biomass pool proportion too
bpGroups <- ppGroups[is.na(biomassBymL$Mprop)]; bpGroups <- bpGroups[grep("BB|PB|DR|DL|DC", bpGroups, invert = TRUE)]
thisIndex <- ppGroups %in% bpGroups
bpProp <- sum(biomassBymL$Biomass[thisIndex])/allBtotal

plot(storeBiomPerc*100, type="h", lwd=3, ylab="Percentage of total biomass", xaxt="n", xlab="")
axis(at=1:nr, labels = thisRoundings, side=1)

test<-biomassByppGroup[thisIndex,] 
test2<-apply(test,2, sum, na.rm=TRUE)
plot(test, type="l")

bpGroups[order(test[,1], decreasing = TRUE)]

highMpropIndex <- biomassBymL>0.8

npg_plot <- npg-1 # going to drop DC as not using it so ~0
orderAll <- order(biomassByppGroup[,1], decreasing = TRUE)
ppOrdered <- ppGroups[orderAll]
BOrdered <- biomassByppGroup[orderAll,1] * mg_2_tonne * X_CN
HPorderedIndex <- ppOrdered %in% ppGroups[highMpropIndex]


thisRoundings <- c(0,0.2,0.4,0.6,0.8); nr <- length(thisRoundings)
biomassBymL$roundedMprop2 <- unlist(lapply(biomassBymL$Mprop, myRounding, fraction=0.2, direction="down"))
biomassBymL$roundedMprop2[biomassBymL$roundedMprop2>0.8]<-0.8

colByR <- rev(colorRampPalette(colors=c("midnightblue",myBlue,myAqua, myLightGrey))(nr))

yaxisLab <- c(1e+2, 1e+4, 1e+6, 1e+8); yaxisAt <- log(yaxisLab, base=10)
thisMax <- log(max(BOrdered[1:npg_plot], na.rm=TRUE),base=10)

thisRoundings_text <- c(expression("" >= 0), expression("" >= 0.2),expression("" >= 0.4),expression("" >= 0.6),expression("" >= 0.8))

pdf(paste(plotPath, "biomassByMprop.pdf", sep=""), width=15, height=8)
par(mar=c(4,4.5,1,1), lend=1)
plot(log(BOrdered[1:npg_plot], base=10), type="h", lwd=5, xaxt="n", xlab="", col=colByR[1], ylab="", yaxt="n", ylim=c(0, thisMax))
par(las=2)
axis(at=1:npg_plot, labels=ppOrdered[1:npg_plot], side=1)
axis(at=yaxisAt, labels = yaxisLab, side=2)
par(las=0)
mtext("Biomass (tonnes)", side=2, adj=0.5, line=3.5)
for( r in 1:nr){
  thisRgroups <- ppGroups[biomassBymL$roundedMprop2==thisRoundings[r]]
  thisCol <- colByR[r]
  RorderedIndex <- ppOrdered %in% thisRgroups
  points(x=(1:npg_plot)[RorderedIndex], y=log(BOrdered[RorderedIndex], base=10), col=thisCol, type="h", lwd=5)
}
legend(legend=thisRoundings_text, col=colByR, seg.len=3, lwd=3, x="topright", title= "mL/M")
dev.off()

biomassBymL$truncProp <- trunc((10*biomassBymL$Mprop))/10
# there is one set to 1, add to 0.9 and define as >=0.9
biomassBymL$truncProp[biomassBymL$truncProp==1]<-0.9
biomassBymL$truncProp[is.na(biomassBymL$truncProp)]<-0
biomassBymL$roundedMprop[is.na(biomassBymL$roundedMprop)]<-0
biomassBymL$roundedMprop2[is.na(biomassBymL$roundedMprop2)]<-0
sumBiom <- tapply(biomassBymL$Biomass, biomassBymL$truncProp, sum, na.rm=TRUE)
plot(x=names(sumBiom), y=log(sumBiom), type="h", xaxt="n", xlab="", ylab="log(Biomass (tonnes))")
sumGroups <- table(biomassBymL$truncProp)
par(new=TRUE)
plot(x=names(sumGroups), y=sumGroups, yaxt="n", ylab="")

# number of groups by mprop, and proportion of biomass by mprop
# for this, take out bacteria, detritus, and primary producers
dropGroups <- c("BB", "PB", "DR", "DL", "DF", "MA", "MB", "DC", "PL", "PS")
biomassByMprop <- tapply(biomassBymL$Biomass[!(ppGroups %in% dropGroups)], biomassBymL$roundedMprop2[!(ppGroups %in% dropGroups)], sum, na.rm=TRUE)
ngroupsByMprop <- table( biomassBymL$roundedMprop2[!(ppGroups %in% dropGroups)])
biomassByNAprop <- sum(biomassBymL$Biomass[!(ppGroups %in% dropGroups) & is.na(biomassBymL$roundedMprop2)])
# test<-tapply()
biomR0_notBP <- biomassByMprop[1]
biomassByMprop[1] <- biomassByNAprop + biomR0_notBP
ngroupsByMprop[1] <- ngroupsByMprop[1] + length(biomassBymL$Mprop[is.na(biomassBymL$roundedMprop2)])
plot(x=(1:nr), y=log(biomassByMprop, base=10), type="n", lwd=5, col=myGold)
points(x=(1:nr)-0.1, y=log(biomassByMprop, base=10), type="h", lwd=5, col=myGold)
# mark where biomass pool biomass comes to
points(x=1, y=log(biomR0_notBP, base=10), pch=8, col=myGold)

par(new=TRUE)
plot(x=(1:nr), y=ngroupsByMprop, type="n", lwd=5, col=myBlue)
points(x=(1:nr)+0.1, y=ngroupsByMprop, type="h", lwd=5, col=myBlue)

# perhaps plot the ordered one with groups dropped
dropIndex <- !(ppGroups %in% dropGroups)
thisbiomassByppGroup <- biomassByppGroup[dropIndex,]; thisppGroups <- ppGroups[dropIndex]; this_npg <- length(thisppGroups)
thisbiomassBymL <- biomassBymL[dropIndex,]
orderAll <- order(thisbiomassByppGroup[,1], decreasing = TRUE)
ppOrdered <- thisppGroups[orderAll]
BOrdered <- thisbiomassByppGroup[orderAll,1] * mg_2_tonne * X_CN

colByR <- rev(colorRampPalette(colors=c("midnightblue",myBlue,myAqua, myLightGrey))(nr))

yaxisLab <- c(1e+2, 1e+4, 1e+6, 1e+8); yaxisAt <- log(yaxisLab, base=10)
thisMax <- log(max(BOrdered, na.rm=TRUE),base=10)

pdf(paste(plotPath, "biomassByMprop_dropPPBactDet.pdf", sep=""), width=14, height=9)
par(mar=c(4,4.5,1,1), lend=1)
plot(log(BOrdered, base=10), type="h", lwd=5, xaxt="n", xlab="", col=colByR[1], ylab="", yaxt="n", ylim=c(0, thisMax))
par(las=2)
axis(at=1:this_npg, labels=ppOrdered, side=1)
axis(at=yaxisAt, labels = yaxisLab, side=2)
par(las=0)
mtext("Biomass (tonnes)", side=2, adj=0.5, line=3.5)
for( r in 1:nr){
  thisRgroups <- thisppGroups[thisbiomassBymL$roundedMprop2==thisRoundings[r]]
  thisCol <- colByR[r]
  RorderedIndex <- ppOrdered %in% thisRgroups
  points(x=(1:this_npg)[RorderedIndex], y=log(BOrdered[RorderedIndex], base=10), col=thisCol, type="h", lwd=5)
}
legend(legend=thisRoundings, col=colByR, seg.len=3, lwd=3, x="topright", title= "mL/M")
dev.off()









