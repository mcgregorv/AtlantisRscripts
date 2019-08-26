## read in tracers from the 'repeat ROMS from 1 year for 50 years' runs (9 of them) and the 'bootstrap ROMS for each of 50 years' runs (50 of them)
## and compare each of the repeats within the range from the bootstrapped - if we picked any one year would it go outside bounds from the bootstraps..?
# 
## took the other (same-named) script and did some testing here - saved for the record
##
##
# tracers were read in and written out in exploringROMS_bootstrapCVs.R
basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
plotPath<-paste(DIR$'Figures',"ROMS\\BootVsRepeat_",sep="")
## repeat ROMS is Version D; bootstrap ROMS is Version B
VersionRepeat<-"D"; versionBoot<-"B"
outPath<-paste(basePath,"base\\ouputROMS",VersionRepeat,"\\",sep="")
dataOutPath<-paste(basePath,"BootstrapROMS\\modelOut",VersionRepeat,"_", sep="")
load(paste(dataOutPath,"modelTracers",sep="")); ##to bring "storeTracers", "storeTracersByCell", "storeTemperature", "store_nts"
## storeTracers and storeTracersByCell are already converted to tonnes
repeatTracers<-storeTracers; repeat_ntsList<-store_nts
repeatTemperature<-storeTemperature
##
dataOutPath<-paste(basePath,"BootstrapROMS\\modelOut",versionBoot,"_", sep="")
load(paste(dataOutPath,"modelTracers",sep="")); ##to bring "storeTracers", "storeTracersByCell", "storeTemperature", "store_nts"
## storeTracers and storeTracersByCell are already converted to tonnes
bootTracers<-storeTracers; boot_ntsList<-store_nts
bootTemperature<-storeTemperature

nBootRuns<-dim(bootTracers)[1]; nRepeatRuns<-dim(repeatTracers)[1]
bootYears<-seq(1996,2004)

groupsDF<-read.csv(paste(basePath,"\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

#############
# for each tracer, plot all the boots in grey_trans, then plot the repeats in black as a first check
## reading a base run in, just to get tracers names
#get _N tracers
BaseNC.nc<-nc_open(paste(basePath, "\\Base\\outputBase\\output.nc",sep=""))
allTracers<-names(BaseNC.nc$var); baseVol<-ncvar_get(BaseNC.nc, "volume")
x<-grep("_N",allTracers); temp<-allTracers[x]
y<-grep("Nums",temp, invert = TRUE); Ntracers<-temp[y]; ntracers<-length(Ntracers); ng<-ntracers


#for each tracer, get it's biomass range from all bootstrap sim and then from all repeat sims
bootstrapMins<-apply(bootTracers,2, min, na.rm=TRUE); bootstrapMaxs<-apply(bootTracers, 2, max, na.rm=TRUE)
repeatMins<-apply(repeatTracers, 2, min, na.rm=TRUE); repeatMaxs<-apply(repeatTracers, 2, max, na.rm=TRUE)
aboveMaxIndex<-repeatMaxs>bootstrapMaxs; belowMins<-repeatMins<bootstrapMins


par(lend=1)
plot((repeatMins/bootstrapMins), type="h", lwd=5)
plot((bootstrapMaxs/repeatMaxs), type="h", lwd=5)

###
belowMin<-array(NA, dim=c(ng, nRepeatRuns)); aboveMax<-belowMin
#from base model (which has ROMS repeated)
baseMins<-rep(NA, nRepeatRuns); baseMaxs<-rep(NA, nRepeatRuns)
for(g in 1:ng){
  thisY<-repeatTracers[,g,]
  thisMins<-apply(thisY,1,min, na.rm=TRUE); thisMaxs<-apply(thisY, 1, max, na.rm=TRUE)
  thisBootMin<-bootstrapMins[g]; thisBootMax<-bootstrapMaxs[g]
  indexMin<-thisMins<thisBootMin; indexMax<-thisMaxs>thisBootMax
  belowMin[g,]<-indexMin; aboveMax[g,]<-indexMax
}

baseMins<-rep(FALSE, ng); baseMaxs<-rep(FALSE, ng)
for(g in 1:ng){
  thisTracer<-Ntracers[g]
  temp<-ncvar_get(BaseNC.nc, thisTracer)
  if(length(dim(temp))==3){
    thisY<-apply(temp * baseVol, 3, sum) * X_CN * mg_2_tonne
  } else{
    thisY<-apply(temp * baseVol[6,,], 2, sum) * X_CN * mg_2_tonne
  }
  thisMins<-min(thisY[36:85], na.rm=TRUE); thisMaxs<-max(thisY[36:85], na.rm=TRUE)
  thisBootMin<-bootstrapMins[g]; thisBootMax<-bootstrapMaxs[g]
  if(thisMins < thisBootMin){ baseMins[g]<-TRUE } 
  if(thisMaxs > thisBootMax){ baseMaxs[g]<-TRUE } 
  
  plot(thisY[36:50], type="l", ylim=c(min(c(thisY, thisBootMin)), max(c(thisY, thisBootMax))))
  abline(h=c(thisBootMin, thisBootMax), col="red")
  mtext(thisTracer, side=3)
}
baseBelowMin<-sum(baseMins); baseAboveMax<-sum(baseMaxs)

belowMinByGroup<-apply(belowMin, 1, sum)
aboveMaxByGroup<-apply(aboveMax, 1, sum)

par(mfrow=c(2,1), mar=c(4,4,1,1))
plot(belowMinByGroup[order(belowMinByGroup)], type="h", lwd=5, xaxt="n", ylab="", xlab="")
par(las=2)
axis(at=seq(1,ng), labels=groupsDF$Code[order(belowMinByGroup)], side=1, cex=0.6)

plot(aboveMaxByGroup[order(aboveMaxByGroup)], type="h", lwd=5, xaxt="n", ylab="", xlab="")
par(las=2)
axis(at=seq(1,ng), labels=groupsDF$Code[order(aboveMaxByGroup)], side=1)
# warm years are 99, 00, 01; cooler years are 96, 97, 04
warmIndex<-bootYears %in% c(1999, 2000, 2001); coolIndex<-bootYears %in% c(1996, 1997, 2004)
aboveMaxCool<-apply(aboveMax[,coolIndex], 1, sum)  ; aboveMaxAny <- apply(aboveMax, 1, sum)
belowMinWarm <- apply(belowMin[, warmIndex], 1, sum) ; belowMinAny <- apply(belowMin, 1, sum)

groupsDF$Code[order(aboveMaxCool, decreasing = TRUE)]
groupsDF$Code[order(belowMinWarm, decreasing = TRUE)]

plot(aboveMaxCool[order(aboveMaxCool, decreasing = TRUE)], type="h", lwd=5, xaxt="n")
par(las=2)
axis(at=seq(1,ng), labels=groupsDF$Code[order(aboveMaxCool, decreasing = TRUE)], side=1)
par(las=1)
mtext("Above max when it was cooler", side=3, adj=0)

plot(belowMinWarm[order(belowMinWarm, decreasing = TRUE)], type="h", lwd=5, xaxt="n")
par(las=2)
axis(at=seq(1,ng), labels=groupsDF$Code[order(belowMinWarm, decreasing = TRUE)], side=1)
par(las=1)
mtext("Below min when it was warmer", side=3, adj=0)

## above or below in the other years
aboveMaxNOTCool<-apply(aboveMax[,!coolIndex], 1, sum)  
belowMinNOTWarm <- apply(belowMin[, !warmIndex], 1, sum) 
nx<-max(c(length(aboveMaxNOTCool[aboveMaxNOTCool>0]), length(aboveMaxCool[aboveMaxCool>0]), 
          length(belowMinNOTWarm[belowMinNOTWarm>0]), length(belowMinWarm[belowMinWarm>0])))+1

par(mfrow=c(2,1), mar=c(4,4,1,1), lend=1)
plot(aboveMaxNOTCool[order(aboveMaxNOTCool, decreasing = TRUE)], type="h", lwd=5, xaxt="n",xlim=c(0,nx))
par(las=2)
axis(at=seq(1,ng), labels=groupsDF$Code[order(aboveMaxNOTCool, decreasing = TRUE)], side=1)
par(las=1)
mtext("Above max when it was NOT cooler", side=3, adj=0)
plot(aboveMaxCool[order(aboveMaxCool, decreasing = TRUE)], type="h", lwd=5, xaxt="n",xlim=c(0,nx))
par(las=2)
axis(at=seq(1,ng), labels=groupsDF$Code[order(aboveMaxCool, decreasing = TRUE)], side=1)
par(las=1)
mtext("Above max when it was cooler", side=3, adj=0)


plot(belowMinNOTWarm[order(belowMinNOTWarm, decreasing = TRUE)], type="h", lwd=5, xaxt="n",xlim=c(0,nx))
par(las=2)
axis(at=seq(1,ng), labels=groupsDF$Code[order(belowMinNOTWarm, decreasing = TRUE)], side=1)
par(las=1)
mtext("Below min when it was NOT warmer", side=3, adj=0)
plot(belowMinWarm[order(belowMinWarm, decreasing = TRUE)], type="h", lwd=5, xaxt="n",xlim=c(0,nx))
par(las=2)
axis(at=seq(1,ng), labels=groupsDF$Code[order(belowMinWarm, decreasing = TRUE)], side=1)
par(las=1)
mtext("Below min when it was warmer", side=3, adj=0)

## do heat map of any above and below, then restrict to cooler/warmer years
getColor<-function(x,thisMax){
  thisCol<-"white"
  if(!is.na(x) & x>0){
    y<-round((x-thisMin)/(thisMax-thisMin),2)*100+1
    thisCol<-thisColRamp[y]
  }
  return(thisCol)
}
plotGrid<-function(x,y){
  thisX<-c(x-0.5,x-0.5,x+0.5,x+0.5); thisY<-c(y-0.5,y+0.5,y+0.5,y-0.5)
  thisCol<-plotColour[x,y]
  if(length(thisCol)>0){
    polygon(x=thisX,y=thisY,col=thisCol,border=NA)
  }
  return(NULL)
}
thisColRamp<-colorRampPalette(colors=c(myLightAqua,myAqua,"midnightblue"))(101)

totalOut<-aboveMaxAny + belowMinAny; totalOrder <- order(totalOut, decreasing = TRUE)

thisMax<-1; legendTextAll<-1:9; legendColAll<-unlist(lapply(legendTextAll/9, getColor, thisMax=thisMax))
legendText3<-1:3; legendCol3<-unlist(lapply(legendText3/3, getColor, thisMax=thisMax))

plotData<-cbind(belowMinAny[totalOrder], aboveMaxAny[totalOrder])/9; thisMax<-max(plotData, na.rm=TRUE); thisMin<-min(plotData, na.rm=TRUE)
plotColour<-apply(plotData,c(1,2),getColor,thisMax)
tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotData)[1]),dim(plotData)[2]),"y"=sort(rep(seq(1,dim(plotData)[2]),dim(plotData)[1]))))
# pdf(paste(plotPath,"TLSensSSRhistoric.pdf",sep=""),height=4,width=5)
par(mar=c(4,4,1.5,2.5))
plot(x=seq(1,dim(plotData)[1]),y=rep(dim(plotData)[2],dim(plotData)[1]),type="n",xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,dim(plotData)[2]))
axis(at= seq(1,ng),labels = groupsDF$Code[totalOrder],side=1,las=2)
axis(at=seq(1,2),labels=c("Below", "Above"),side=2,las=1)
temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
box()
legend(legend=c(legendText3, rep("",6), legendTextAll), col=c(legendCol3, rep(NA, 6), legendColAll), x="right", pch=15, pt.cex=1.5, cex=0.6, ncol=2)

# dev.off()


plotData<-cbind(belowMinWarm[totalOrder], aboveMaxCool[totalOrder])/3; ## 3 years in each (cooler/warm) set
thisMax<-max(plotData, na.rm=TRUE); thisMin<-min(plotData, na.rm=TRUE)
plotColour<-apply(plotData,c(1,2),getColor,thisMax)
tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotData)[1]),dim(plotData)[2]),"y"=sort(rep(seq(1,dim(plotData)[2]),dim(plotData)[1]))))
# pdf(paste(plotPath,"TLSensSSRhistoric.pdf",sep=""),height=4,width=5)
par(mar=c(4,4,1.5,2.5))
plot(x=seq(1,dim(plotData)[1]),y=rep(dim(plotData)[2],dim(plotData)[1]),type="n",xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,dim(plotData)[2]))
axis(at= seq(1,ng),labels = groupsDF$Code[totalOrder],side=1,las=2)
axis(at=seq(1,2),labels=c("Below", "Above"),side=2,las=1)
temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
box()
legend(legend=c(legendText3, rep("",6), legendTextAll), col=c(legendCol3, rep(NA, 6), legendColAll), x="right", pch=15, pt.cex=1.5, cex=0.6, ncol=2)

# dev.off()


plotData<-cbind( belowMinWarm[totalOrder]/3, belowMinAny[totalOrder]/9); ## 3 years in each (cooler/warm) set
thisMax<-max(plotData, na.rm=TRUE); thisMin<-min(plotData, na.rm=TRUE)
plotColour<-apply(plotData,c(1,2),getColor,thisMax)
tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotData)[1]),dim(plotData)[2]),"y"=sort(rep(seq(1,dim(plotData)[2]),dim(plotData)[1]))))
# pdf(paste(plotPath,"TLSensSSRhistoric.pdf",sep=""),height=4,width=5)
par(mar=c(4,5,1.5,2.5))
plot(x=seq(1,dim(plotData)[1]),y=rep(dim(plotData)[2],dim(plotData)[1]),type="n",xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0.5,(dim(plotData)[2]+0.5)))
axis(at= seq(1,ng),labels = groupsDF$Code[totalOrder],side=1,las=2)
axis(at=seq(1,2),labels=c("Below\nwarm", "Below\nall"),side=2,las=1)
temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
box()
legend(legend=c(legendText3, rep("",6), legendTextAll), col=c(legendCol3, rep(NA, 6), legendColAll), x="right", pch=15, pt.cex=1.5, cex=0.6, ncol=2)

# dev.off()


plotData<-cbind( aboveMaxCool[totalOrder]/3, aboveMaxAny[totalOrder]/9); ## 3 years in each (cooler/warm) set
thisMax<-max(plotData, na.rm=TRUE); thisMin<-min(plotData, na.rm=TRUE)
plotColour<-apply(plotData,c(1,2),getColor,thisMax)
tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotData)[1]),dim(plotData)[2]),"y"=sort(rep(seq(1,dim(plotData)[2]),dim(plotData)[1]))))
# pdf(paste(plotPath,"TLSensSSRhistoric.pdf",sep=""),height=4,width=5)
par(mar=c(4,5,1.5,2.5))
plot(x=seq(1,dim(plotData)[1]),y=rep(dim(plotData)[2],dim(plotData)[1]),type="n",xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0.5,(dim(plotData)[2]+0.5)))
axis(at= seq(1,ng),labels = groupsDF$Code[totalOrder],side=1,las=2)
axis(at=seq(1,2),labels=c("Above\ncool", "Above\nall"),side=2,las=1)
temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
box()
legend(legend=c(legendText3, rep("",6), legendTextAll), col=c(legendCol3, rep(NA, 6), legendColAll), x="right", pch=15, pt.cex=1.5, cex=0.6, ncol=2)
# dev.off()

groupsAboveWhenCooler<-groupsDF$Code[aboveMaxCool>0]
groupsAboveWhenNOTCooler<-groupsDF$Code[aboveMaxNOTCool>0]
groupsAboveWhenCooler[!(groupsAboveWhenCooler %in% groupsAboveWhenNOTCooler)]

groupsBelowWhenWarmer<-groupsDF$Code[belowMinWarm>0]
groupsBelowWhenNOTWarmer<-groupsDF$Code[belowMinNOTWarm>0]
groupsBelowWhenWarmer[!(groupsBelowWhenWarmer %in% groupsBelowWhenNOTWarmer)]

## zero those that are below (resp above) in at least 3 other years
belowNOTwarmIndex<-belowMinNOTWarm>0
thisBelowMin <- belowMinWarm[belowNOTwarmIndex]; thisBelowMinGroups<-groupsDF$Code[belowNOTwarmIndex][thisBelowMin>0]

aboveNOTcoolIndex<-aboveMaxNOTCool>0
thisAboveMax <- aboveMaxCool[aboveNOTcoolIndex]; thisAboveMaxGroups<-groupsDF$Code[aboveNOTcoolIndex][thisAboveMax>0]

## groups out in the other years
groupsOutOtherIndex<- belowMinNOTWarm>0 | aboveMaxNOTCool>0
thisBelowMin <- belowMinWarm[!groupsOutOtherIndex];  thisBelowMinGroups<-groupsDF$Code[!groupsOutOtherIndex][thisBelowMin>0]
thisAboveMax <- aboveMaxCool[!groupsOutOtherIndex];  thisAboveMaxGroups<-groupsDF$Code[!groupsOutOtherIndex][thisAboveMax>0]

combineAboveBelowGroups<-unique(sort(c(as.character(thisAboveMaxGroups), as.character(thisBelowMinGroups))))
ca_df<-data.frame(cbind(combineAboveBelowGroups, "Below"= as.double(thisBelowMin[thisBelowMin>0][match(combineAboveBelowGroups, thisBelowMinGroups)]), 
                        "Above"=as.double(thisAboveMax[thisAboveMax>0][match(combineAboveBelowGroups, thisAboveMaxGroups)])))

ca_df$Below<-as.double(ca_df$Below)
ca_df$Below[is.na(ca_df$Below)]<-0
ca_df$Total<-as.double(ca_df$Below) + as.double(ca_df$Above)

groupsAboveAll <- groupsDF$Code[aboveMaxAny>0]; groupsBelowAll <- groupsDF$Code[belowMinAny>0]
combineAboveBelowGroups<-unique(sort(c(as.character(groupsAboveAll), as.character(groupsBelowAll))))
ca_df<-data.frame(cbind("Code"=combineAboveBelowGroups, "Below"= as.double(belowMinAny[belowMinAny>0][match(combineAboveBelowGroups, groupsBelowAll)]), 
                        "Above"=as.double(aboveMaxAny[aboveMaxAny>0][match(combineAboveBelowGroups, groupsAboveAll)])))

ca_df$Below<-as.double(ca_df$Below); ca_df$Above <- as.double(ca_df$Above)
ca_df$Below[is.na(ca_df$Below)]<-0; ca_df$Above[is.na(ca_df$Above)]<-0
ca_df$Total<-as.double(ca_df$Below) + as.double(ca_df$Above)
ca_df[order(ca_df$Total, decreasing = TRUE),]

totalWarmCool<- belowMinWarm + aboveMaxCool; totalWarmCoolOrder <- order(totalWarmCool, decreasing = TRUE)
groupCodes<-groupsDF$Code
groupsInBothIndex<-groupCodes %in% groupsAboveWhenCooler & groupCodes %in% groupsBelowWhenWarmer
groupsInBothOrder<-order(totalWarmCool[groupsInBothIndex], decreasing = TRUE)
groupsInBothOrderd <- groupCodes[groupsInBothIndex][groupsInBothOrder]
groupsCountsInBothOrderd <- totalWarmCool[groupsInBothIndex][groupsInBothOrder]

par(mfrow=c(1,1), mar=c(1,4,2,10))
plotData<-t(cbind( aboveMaxCool[totalWarmCoolOrder]/3, belowMinWarm[totalWarmCoolOrder]/3)); ## 3 years in each (cooler/warm) set
temp<-apply(plotData, 2, sum)
plotData <- plotData[,temp>0]
thisMax<-max(plotData, na.rm=TRUE); thisMin<-min(plotData, na.rm=TRUE)
plotColour<-apply(plotData,c(1,2),getColor,thisMax)
tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotData)[1]),dim(plotData)[2]),"y"=sort(rep(seq(1,dim(plotData)[2]),dim(plotData)[1]))))
# pdf(paste(plotPath,"TLSensSSRhistoric.pdf",sep=""),height=4,width=5)
plot(x=seq(1,dim(plotData)[1]),y=rep(dim(plotData)[2],dim(plotData)[1]),type="n",xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0.5,(dim(plotData)[2]+0.5)), xlim=c(0.5,2.5))
axis(at= seq(1,ng)[temp>0],labels = groupsDF$Code[totalWarmCoolOrder][temp>0],side=2,las=2)
axis(at=seq(1,2),labels=c("Above cool", "Below warm"),side=3,las=1)
temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
box()
par(xpd=NA)
legend(legend=c(legendText3), inset=-0.1, col=c(legendCol3), x="right", pch=15, pt.cex=1.5, cex=1)
# dev.off()



temp<-groupsDF$Code[aboveMaxCool>2]
temp[!(temp %in% groupsAboveWhenNOTCooler)]


plot(x=belowMinByGroup, y=aboveMaxByGroup, pch=20, col=myGrey_trans)


belowMinByYear<-apply(belowMin, 2, sum)
aboveMaxByYear<-apply(aboveMax, 2, sum)

temp<-cbind(aboveMaxByYear, belowMinByYear); colnames(temp)<-c("above","below");  rownames(temp)<-as.character(bootYears)
temp<-rbind(temp,c(baseAboveMax, baseBelowMin)); rownames(temp)[dim(temp)[1]]<-"Base"
toPlot<-melt(temp);

pdf(paste(plotPath,"BootstrapOutsideLimitsByROMSyear.pdf", sep=""), height=3.5, width=7.5)
par(mar=c(4,4.5,0.5,4))
bp<-ggplot(data = toPlot, aes(x = as.character(X1), fill = X2, y = value)) + 
  geom_bar(stat = 'identity')
thisCols<-c(myLightBlue, "midnightblue")
bp   + scale_fill_manual(values=thisCols) + labs(y="Number of species groups", x="Year") + theme_igray() +
  theme(axis.text=element_text(size=12, angle=0),axis.title=element_text(size=12)) + guides(fill=guide_legend(title="")) 
dev.off()

## wrt group
temp<-cbind(aboveMaxByGroup, belowMinByGroup); colnames(temp)<-c("above","below");  rownames(temp)<-Ntracers
toPlot<-melt(temp);

bp<-ggplot(data = toPlot, aes(x = as.character(X1), fill = X2, y = value)) + 
  geom_bar(stat = 'identity')
thisCols<-c(myLightBlue, "midnightblue")
bp  + coord_flip()   + scale_fill_manual(values=thisCols) + labs(y="Number of models", x="Species group") + theme_igray() +
  theme(axis.text=element_text(size=12, angle=0),axis.title=element_text(size=12)) + guides(fill=guide_legend(title="")) 

###############
##group by trophic level  ## NEED TO REORDER IF GOING TO USE!!   groups are in order of tracer names, not DF codes as was assumed here
temp<-read.csv(paste(basePath,"inputs\\biol_prm\\CRAM_trophicLevels_isotopes.csv",sep=""))
groupsTL<-temp
groupsTL$Isotope[is.na(groupsTL$Isotope)]<-groupsTL$TrophicLevel2[is.na(groupsTL$Isotope)]
TLindex<-order(groupsTL$Isotope)

labelNames<-gsub("_","", groupsDF$LongName[TLindex])
labelNames<-gsub("_"," ", groupsDF$Name[TLindex])

temp<-cbind(aboveMaxByGroup[TLindex], belowMinByGroup[TLindex]); colnames(temp)<-c("above","below");  rownames(temp)<-as.character(labelNames)
skipGroup<-c("Carrion")
temp<-temp[!(rownames(temp) %in% skipGroup), ]

toPlot<-melt(temp);

pdf(paste(plotPath,"BootstrapOutsideLimitsByGroup.pdf", sep=""), height=8, width=7)
par(mar=c(4,9,0.5,4))
bp<-ggplot(data = toPlot, aes(x = as.character(X1), fill = X2, y = as.double(value))) + 
  geom_bar(stat = 'identity')
thisCols<-c(myLightBlue, "midnightblue")
bp  + coord_flip()   + scale_fill_manual(values=thisCols) + labs(y="Number of model simulations", x="") + theme_igray() +
  theme(axis.text=element_text(size=12, angle=0),axis.title=element_text(size=12)) + guides(fill=guide_legend(title="")) 
dev.off()


### just 2003 by group
belowMinByGroup<-belowMin[,bootYears==2003]
aboveMaxByGroup<-aboveMax[,bootYears==2003]


#########################

#calculate a version of RI wrt bounds
calc_RI_bounds<-function(P, Ol, Ou){
  A=0; n<-length(Ol)
  for(i in 1:n){
    if((P[i]<=Ou[i]) & (P[i]>=Ol[i])){
      thisP<-Ol[i]; 
      thisO<-thisP
    } else if(P[i] > Ou[i]){
      thisP <- P[i]; thisO <- Ou[i]
    } else if(P[i] < Ol[i]){
      thisP <- P[i]; thisO <- Ol[i]
    }
    C<-(log(thisO/thisP))^2
    A = A + C
  }
  RI<-exp(sqrt(A/n))
  return(RI)
}


storeRIs<-array(NA, dim=dim(repeatTracers)[c(1,2)])

for(t in 1:ntracers){
  thisTracer<-Ntracers[t]; thisTracerName<-gsub("_|_N", " ", thisTracer)
  
  thisOl<-rep(bootstrapMaxs[t], dim(repeatTracers)[3])
  thisOu<-rep(bootstrapMins[t], dim(repeatTracers)[3])
  thisRI<-apply(repeatTracers[,t,], 1, calc_RI_bounds, Ol=thisOl, Ou=thisOu)
  storeRIs[,t]<-thisRI
}
meanByTracers<-apply(storeRIs, 2, mean); minByTracers<-apply(storeRIs, 2, min); maxByTracers<-apply(storeRIs, 2, max)
focusIndex<-(maxByTracers-minByTracers)>0.1; focusTracers<-Ntracers[focusIndex]

colByRun<-colorRampPalette(colors=c(myGold, myOrange, "red",myPurple,myBlue))(dim(storeRIs)[1])

par(mar=c(9,4,1,1))
plot(storeRIs[1,focusIndex], type="n", ylim=c(0, max(storeRIs)), xaxt="n", ylab="RI", xlab="")
par(las=2)
axis(at=seq(1,length(focusTracers)), labels=gsub("_|_N", " ", focusTracers), side=1)
for(r in 1:(dim(storeRIs  )[1])){
  thisCol<-colByRun[r]
  points(storeRIs[r,focusIndex], pch=r, col=thisCol)
}
legend(legend=seq(1996, 2004), col=colByRun, pch=seq(1,length(colByRun)), x="topleft", bty="n", ncol=5)


meanByRun<-apply(storeRIs, 1, sum)


plot(meanByTracers, type="h", lwd=3)
abline(h=1, col="red", lty=2)
# 
# pdf(paste(plotPath,"AlltracersBounds.pdf", sep=""), height=8, width=8)
# par(mfrow=c(5,2), mar=c(4,4,1,1), oma=c(0,0,0,0))
# for(t in 1:ntracers){
#   thisTracer<-Ntracers[t]; thisTracerName<-gsub("_|_N", " ", thisTracer)
#   
#   thisOl<-rep(bootstrapMaxs[t], dim(repeatTracers)[3])
#   thisOu<-rep(bootstrapMins[t], dim(repeatTracers)[3])
#   thisRI<-apply(repeatTracers[,t,], 1, calc_RI_bounds, Ol=thisOl, Ou=thisOu)
#   
#   plot(thisOu, type='n', ylim=c(0, max(c(thisOu, thisMP))*1.2), ylab="Biomass (tonnes)", xlab="Timestep (years)")
#   polygon(x=c(seq(1,length(thisOu)), rev(seq(1, length(thisOu)))), y=c(thisOl, rev(thisOu)), col=myGrey_trans, border="NA")
#   points(repeatTracers[1,t,], type="l")
#   for(r in 1:9){
#     points(x=seq(1,length(thisOu)), y=repeatTracers[r,t,], type="l")
#   }
#   mtext(thisTracerName, side=3, adj=0)
# }
# dev.off()
# 
# 
# pdf(paste(plotPath,"AlltracersAllruns.pdf", sep=""), height=10, width=8)
# par(mfrow=c(5,1), mar=c(4,4,1,1), oma=c(0,0,0,0))
# for(t in 1:ntracers){
#   thisTracer<-Ntracers[t]; thisTracerName<-gsub("_|_N", " ", thisTracer)
#   thisMax<-max(c(bootTracers[,t,], repeatTracers[,t,]), na.rm=TRUE)
#   plot(bootTracers[,t,1], type="n", ylim=c(0, thisMax), ylab="Biomass (tonnes)", xlab="Timestep (years)")
#   for(r in 1:nBootRuns){
#     points(bootTracers[r,t,], type="l", lwd=2, col=myGrey_trans)
#   }
#   for(r in 1:nRepeatRuns){
#     points(repeatTracers[r,t,], type="l", lwd=1, col="black")
#   }
#   mtext(thisTracerName, side=3, adj=0)
# }
# dev.off() 
####################################

