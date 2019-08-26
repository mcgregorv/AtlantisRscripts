##intended for fishing scenarios initially, but will work for any multipl runs. 
## supply a base run, some/one alternative runs, some/one timesteps for which you want the comparison
this_run<-"base"
this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
plotPath<-paste(this_path,"..\\Figures\\CompareRuns_PreSENS_",sep="")

version<-"FS_C1"

baseOut<-"FS_C1All100catch"
scenarioScalars<-c(0,1,0.5,0.6,0.7,0.8,0.5,0.4,0.3,0.2,1.2,1.3,1.4,1.5,1.2,1.3,1.4,1.5); nScenarios<-length(scenarioScalars)
scenarioGroups<-c(rep("All", 6),rep("Hoki",8),rep("All",4))
scenarioCodes<-mapply(FUN=function(x,y){paste(x,y*100, "catch",sep="")},x=scenarioGroups,y=scenarioScalars)

# 
# scenarioScalars<-c(0,1,0.5,0.6,0.7,0.8,0.5,0.4,0.3,0.2,1.2,1.3,1.4,1.5,1.2,1.3,1.4,1.5); nScenarios<-length(scenarioScalars)
# scenarioCodes<-mapply(FUN=function(x,y){paste(x,y*100, "catch",sep="")},x=scenarioGroups,y=scenarioScalars)

altOuts<-unlist(lapply(scenarioCodes,FUN=function(x){paste(version,x,sep="")}))
# baseOut<-"PreSENS"
# altOuts<-c("XXX_mQB1", "XXX_mQB2", "XXX_mQB3", "XXX_mQB4", "XXX_mQA1")
nScenarios<-length(altOuts)

outPath<-paste(this_path,"output",baseOut,"\\",sep=""); baseNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(baseNC.nc,"volume")
nts<-dim(thisVol)[3]; nboxes<-dim(thisVol)[2]; nlayers<-dim(thisVol)[1]

#want the _N tracers (mg N per m^3)
allTracers<-sort(names(baseNC.nc$var))
temp<-allTracers[grep("_N",allTracers)]; tracers2plot<-temp[grep("Nums",temp,invert = TRUE)]; ntracers<-length(tracers2plot)

groups2plot<-c("HOK", "ORH", "SSO", "LIN", "BOE", "PFM", "PFS");
# groups2plot<-c("HOK", "ORH", "SSO", "LIN", "BOE", "PFM"); 
npg<-length(groups2plot)
##also want some others, like total biomass, total biomass by TL?
groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep=""))

storeBiomassArray<-array(NA,dim=c(npg,nScenarios+1,nts))

for(g in 1:npg){
  thisCode<-groups2plot[g]
  thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode],side="both")
  thisTracer<-paste(thisName,"_N",sep="")
  thisData<-ncvar_get(baseNC.nc,thisTracer)
  thisTimeSeries<-apply(thisData*thisVol,3,sum)*mg_2_tonne*X_CN
  storeBiomassArray[g,1,]<-thisTimeSeries
}

for(s in 1:nScenarios){
  thisOutPath<-paste(this_path,"output",altOuts[s],"\\",sep=""); tempNC.nc<-nc_open(paste(thisOutPath,"output.nc",sep=""))
  tempVol<-ncvar_get(tempNC.nc,"volume"); temp_nts<-dim(tempVol)[3]
  for(g in 1:npg){
    thisCode<-groups2plot[g]
    thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode],side="both")
    thisTracer<-paste(thisName,"_N",sep="")
    thisData<-ncvar_get(tempNC.nc,thisTracer)
    thisTimeSeries<-apply(thisData*tempVol,3,sum,na.rm=TRUE)*mg_2_tonne*X_CN
    storeBiomassArray[g,s+1,1:min(temp_nts,nts)]<-thisTimeSeries[1:min(temp_nts, nts)]
    cat(temp_nts,"--")
  }
}

propOfBase<-storeBiomassArray[,-1,]
diffsFromBase<-storeBiomassArray[,-1,]
propDiffFromBase<-storeBiomassArray[,-1,]

for(s in 1:nScenarios){
  propOfBase[,s,]<-storeBiomassArray[,s+1,]/storeBiomassArray[,1,]
  diffsFromBase[,s,]<-storeBiomassArray[,s+1,]-storeBiomassArray[,1,]
  propDiffFromBase[,s,]<-(storeBiomassArray[,s+1,]-storeBiomassArray[,1,])/storeBiomassArray[,1,]
}

############################################################################################################
############################################################################################################
cat(nts)

this_ts<-40
thisCatches<-data.frame(t(storeBiomassArray[,,this_ts]))
colnames(thisCatches)<-groups2plot
thisCatches$Scenario<-c("Base",altOuts)

catchColors<-c(colorRampPalette(colors=c(myGold,myGreen,myAqua,myBlue,myPurple,myRed,myGrey))(npg))


thisPlotData<-melt(thisCatches,id.var="Scenario")
bp<-ggplot(thisPlotData, aes(x = Scenario, y = value, fill = variable)) + 
  geom_bar(stat = "identity")

# pdf(paste(plotPath,"TopCatchByYear.pdf",sep=""),height=5,width=7)
bp + scale_fill_manual(values=catchColors) + labs(y="Biomass (tonnes)") + theme_igray() + 
  theme(axis.text=element_text(size=14, angle=90),axis.title=element_text(size=14)) + guides(fill=guide_legend(title=""))
# dev.off()

#do diffs from base
thisDiffs<-data.frame(t(diffsFromBase[,,this_ts]))
colnames(thisDiffs)<-groups2plot
thisMax<-max(thisDiffs); thisMin<-min(thisDiffs)

#do relative too
thisPropDiffs<-data.frame(t(propDiffFromBase[,,this_ts]))
colnames(thisPropDiffs)<-groups2plot


scenColors<-colorRampPalette(colors=c(myGold,myYellow,myGreen,myDarkGreen))(nScenarios)

plot(x=seq(1,npg),y=rep(1,npg),ylim=c(thisMin,thisMax),type="n", xlab="", ylab="",xaxt="n")
abline(h=1,col=myGrey,lwd=4)

shifts<-seq(-.2,.2,length.out=nScenarios)
for(s in 1:nScenarios){
  points(x=(seq(1,npg)+shifts[s]), y=thisDiffs[s,],type="h",lwd=5,col=scenColors[s],lend=1)
}
axis(at=seq(1,npg),labels=groups2plot,side=1)
mtext("Biomass relative to base (tonnes)",side=2,adj=0.5, line=2.5)

## do difference in biomass, seperating out the hoki only changes from the all changes
hokiRunIndex<-grep("Hoki",scenarioCodes)
hokiGroups<-c("HOK", "PFS"); nhg<-length(hokiGroups)
hokiGroupsIndex<-groups2plot %in% hokiGroups
hokiScens<-scenarioCodes[hokiRunIndex]; nhs<-length(hokiScens)
hokiScenColors<-colorRampPalette(colors=c(myOrange,"red",myRed,myPurple,"midnightblue"))(nhs)
hokiDiffs<-thisDiffs[hokiRunIndex,hokiGroupsIndex]
thisMax<-max(hokiDiffs); thisMin<-min(hokiDiffs)

plot(x=seq(1,nhg),y=rep(1,nhg),ylim=c(thisMin,thisMax),type="n", xlab="", ylab="",xaxt="n",xlim=c(0.5,nhg+0.5))
abline(h=1,col=myGrey,lwd=4)

shifts<-seq(-.2,.2,length.out=nhs)
for(s in 1:nhs){
  points(x=(seq(1,nhg)+shifts[s]), y=hokiDiffs[s,],type="h",lwd=5,col=hokiScenColors[s],lend=1)
}
axis(at=seq(1,nhg),labels=hokiGroups,side=1)
mtext("Biomass relative to base (tonnes)",side=2,adj=0.5, line=2.5)

ggHokDiffs<-data.frame(hokiDiffs); ggHokDiffs$Scenario<-hokiScens
plotHokiDiffs<-melt(ggHokDiffs, id.vars = "Scenario")

bp<-ggplot(plotHokiDiffs, aes(x = Scenario, y = value, fill = variable)) + 
  geom_bar(stat = "identity")

# pdf(paste(plotPath,"TopCatchByYear.pdf",sep=""),height=5,width=7)
bp + scale_fill_manual(values=hokiScenColors) + labs(y="Change in biomass (tonnes)") + theme_igray() + 
  theme(axis.text=element_text(size=14, angle=90),axis.title=element_text(size=14)) + guides(fill=guide_legend(title="")) + geom_hline(yintercept = 0) 
# dev.off()

##do relative to base biomass - either chose run with no fishing or with default (in this case average from last 5 years) fishing
hokiPropDiffs<-thisPropDiffs[hokiRunIndex,hokiGroupsIndex]
thisMax<-max(hokiPropDiffs); thisMin<-min(hokiPropDiffs)

ggHokDiffs<-data.frame(hokiPropDiffs); ggHokDiffs$Scenario<-hokiScens
plotHokiDiffs<-melt(ggHokDiffs, id.vars = "Scenario")

bp<-ggplot(plotHokiDiffs, aes(x = Scenario, y = value, fill = variable)) + 
  geom_bar(stat = "identity")

# pdf(paste(plotPath,"TopCatchByYear.pdf",sep=""),height=5,width=7)
bp + scale_fill_manual(values=hokiScenColors) + labs(y="Proportional change in biomass") + theme_igray() + 
  theme(axis.text=element_text(size=14, angle=90),axis.title=element_text(size=14)) + guides(fill=guide_legend(title="")) + geom_hline(yintercept = 0) 
# 

#############################################################################################
#########################################################
## do the all catch change version, but just plotting hoki and PFS still
hokiRunIndex<-grep("Hoki",scenarioCodes, invert = TRUE); hokiPropDiffs<-thisPropDiffs[hokiRunIndex,hokiGroupsIndex]
thisMax<-max(hokiPropDiffs); thisMin<-min(hokiPropDiffs)
hokiScens<-scenarioCodes[hokiRunIndex]; nhs<-length(hokiScens)
ggHokDiffs<-data.frame(hokiPropDiffs); ggHokDiffs$Scenario<-hokiScens
plotHokiDiffs<-melt(ggHokDiffs, id.vars = "Scenario")

bp<-ggplot(plotHokiDiffs, aes(x = Scenario, y = value, fill = variable)) + 
  geom_bar(stat = "identity")

# pdf(paste(plotPath,"TopCatchByYear.pdf",sep=""),height=5,width=7)
bp + scale_fill_manual(values=hokiScenColors) + labs(y="Proportional change in biomass") + theme_igray() + 
  theme(axis.text=element_text(size=14, angle=90),axis.title=element_text(size=14)) + guides(fill=guide_legend(title="")) + geom_hline(yintercept = 0) 
# 

#all groups, and plot them all
plotData<-t(thisPropDiffs)
colnames(plotData)<-hokiScens


thisMax<-max(thisPropDiffs); thisMin<-min(thisPropDiffs)

#do relative too
thisPropDiffs<-data.frame(t(propDiffFromBase[,,this_ts]))
colnames(thisPropDiffs)<-groups2plot

plot(x=seq(1,npg),y=rep(1,npg),ylim=c(thisMin,thisMax),type="n", xlab="", ylab="",xaxt="n")
abline(h=1,col=myGrey,lwd=4)
shifts<-seq(-.2,.2,length.out=nScenarios)
for(s in 1:nScenarios){
  points(x=(seq(1,npg)+shifts[s]), y=thisPropDiffs[s,],type="h",lwd=5,col=scenColors[s],lend=1)
}
axis(at=seq(1,npg),labels=groups2plot,side=1)
mtext("Biomass relative to base (tonnes)",side=2,adj=0.5, line=2.5)

##testing #############################################################################################################################
allScens<-scenarioCodes[grep("Hok",scenarioCodes,invert = TRUE)]
test<-data.frame(cbind("HOK"=thisPropDiffs$"HOK"[scenarioCodes %in% allScens],"Scenario"=allScens))




############################################################################################################
thisProps<-data.frame(t(propOfBase[,,this_ts]))
colnames(thisProps)<-groups2plot
thisProps$Scenario<-altOuts

propPlotData<-melt(thisProps,id.var="Scenario")
bp<-ggplot(propPlotData, aes(x = Scenario, y = value, fill = variable)) + 
  geom_bar(stat = "identity")

# pdf(paste(plotPath,"TopCatchByYear.pdf",sep=""),height=5,width=7)
bp + scale_fill_manual(values=catchColors) + labs(y="Biomass (tonnes)") + theme_igray() + 
  theme(axis.text=element_text(size=14, angle=90),axis.title=element_text(size=14)) + guides(fill=guide_legend(title=""))
# dev.off()

# 
# pp<-ggplot(propPlotData, aes(x = Scenario, y = value, fill = variable)) +
#   geom_line(data = propPlotData, aes(x = Scenario, y = value))
# 
