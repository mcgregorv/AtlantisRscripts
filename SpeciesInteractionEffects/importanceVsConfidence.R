# plotPath<-paste(DIR$'Figures',"\\Sensitivities\\", sep="") ## overwrites figure in paper
this_run<-"base"
basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
outPath<-paste(basePath, "base\\IntSims\\", sep="")
plotPath<-paste(basePath,"Figures\\intEffects\\",sep="")
# plotPath<-paste(DIR$'Figures',"\\Sensitivities\\", sep="") ## overwrites figure in paper

dataOutPath<-outPath

groupsDF<-read.csv(paste(basePath,"\\CRAM_groups.csv", sep=""))
groupsDFPaper<-read.csv(paste(basePath,"\\CRAM_groupsPaper.csv", sep=""))
ageStructuredNames<-gsub("_", " ", groupsDFPaper$Name[groupsDF$NumCohorts>1])

load(paste(dataOutPath,"keyStoneAS",sep="")) #brings in keyStonesAS
load(paste(dataOutPath,"responsiveness",sep="")) # brings in responsiveness
load(file=paste(dataOutPath,"interactionData",sep=""))

perfAn<-read.csv(paste(DIR$'Tables',"speciesGroupPerformanceAnalyses.csv", sep=""))

thisX<-perfAn$SCORE_fit + perfAn$SCORE_knowledge

thisY<-keyStoneAS + responsiveness

plot(x=thisX, y=thisY, pch=20, col=myBlue, cex=1.5, xlab="Model performance + status of knowledge", ylab="Keystoneness + responsiveness")


keystoneOrderNumber<-rev(order(keyStoneAS + responsiveness)); 
keystoneOrderNumber<-rev(order(keyStoneAS )); 
groupsDF$Code[asIndex][keystoneOrderNumber]

nag<-length(ageStructuredNames)
legendText<-c()
for(f in 1:nag){
  thisNum<-f; thisName<-ageStructuredNames[keystoneOrderNumber][f]
  thisLegend<-paste(thisNum, thisName,sep=" ")
  legendText<-c(legendText,thisLegend)
}
colRamp<-colorRampPalette(colors=c("red",myOrange,myAqua,myBlue))(4)

uniqueGradings<-sort(unique(thisX))
getColor<-function(x){
  # y<-(x-17 + 0.5)*2
  thisCol<-colRamp[match(x,uniqueGradings)]
  return(thisCol)
}
colByPer<-unlist(lapply(thisX, getColor))

plot(x=seq(1,length(keyStoneAS)), y=(keyStoneAS + responsiveness)[keystoneOrderNumber],col=colByPer[keystoneOrderNumber], pch=20, cex=2)

# put together data to output as table
rankTable <- data.frame(matrix(NA, ncol=4, nrow=length(keyStoneAS))); colnames(rankTable)<-c("Group", "KeyRank","respRank", "informPerformRank")
rankTable$Group <- legendText
rankTable$KeyRank <- 1:nag
## get responsiveness order- the ?? way..
rorder <- order(responsiveness, decreasing = TRUE)
rankTable$respRank <- seq(1,nag)[match(seq(1,nag), rorder)]
rankTable$informPerformRank <- thisX[keystoneOrderNumber]
rankTable$Code <- groupsDF$Code[asIndex][keystoneOrderNumber]
#write this out
write.csv(rankTable, paste(DIR$'Table',"interactionEffectsRANKINGs.csv", sep=""),row.names = FALSE)


pdf(paste(plotPath,"importanceByPerformance1.pdf",sep=""), width=7,height=5)
par(mar=c(4.5,4.5,0.5,0.5))
plot(x=keyStoneAS, y=responsiveness,type="n", ylab="Responsiveness", xlab="Keystoneness", cex.lab=thisCex, cex.axis=thisCex)
text(x=keyStoneAS[keystoneOrderNumber],y=responsiveness[keystoneOrderNumber],cex=1.2,col=colByPer[keystoneOrderNumber])
dev.off()

pdf(paste(plotPath,"importanceByPerformance1LEGEND.pdf",sep=""), width=10,height=6)
par(mar=c(0,0,0,0))
makeBlankPlot()
legend(legendText,pch=15,lwd=NA,lty=NA,col=colByPer[keystoneOrderNumber],x="center",ncol=3,cex=1.4, bty="n")
dev.off()

perfUnique<-sort(unique(thisX))

perfText<-c("No data gaps, performed well, abundance index available",  "Slight data gaps and/or poor performance",  "Some data gaps and/or poor performance","Poorly specified")
pdf(paste(plotPath,"perfLEGEND.pdf",sep=""), width=7,height=3)
par(mar=c(0,1,0,0))
makeBlankPlot()
legend(legend=perfText,pch=15,lwd=NA,lty=NA,col=rev(colRamp[seq(1,length(colRamp))]),x="center",ncol=1,cex=1.4, pt.cex=3, bty="n", title="")
dev.off()



