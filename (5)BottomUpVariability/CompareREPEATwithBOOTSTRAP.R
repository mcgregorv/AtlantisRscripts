## read in tracers from the 'repeat ROMS from 1 year for 50 years' runs (9 of them) and the 'bootstrap ROMS for each of 50 years' runs (50 of them)
## and compare each of the repeats within the range from the bootstrapped - if we picked any one year would it go outside bounds from the bootstraps..?
# 
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
# 
belowMinByGroup<-apply(belowMin, 1, sum)
aboveMaxByGroup<-apply(aboveMax, 1, sum)

# warm years are 99, 00, 01; cooler years are 96, 97, 04
warmIndex<-bootYears %in% c(1999, 2000, 2001); coolIndex<-bootYears %in% c(1996, 1997, 2004)
aboveMaxCool<-apply(aboveMax[,coolIndex], 1, sum)  ; aboveMaxAny <- apply(aboveMax, 1, sum)
belowMinWarm <- apply(belowMin[, warmIndex], 1, sum) ; belowMinAny <- apply(belowMin, 1, sum)

groupsAboveAll <- groupsDF$Code[aboveMaxAny>0]; groupsBelowAll <- groupsDF$Code[belowMinAny>0]
combineAboveBelowGroups<-unique(sort(c(as.character(groupsAboveAll), as.character(groupsBelowAll))))
ca_df<-data.frame(cbind("Code"=combineAboveBelowGroups, "Below"= as.double(belowMinAny[belowMinAny>0][match(combineAboveBelowGroups, groupsBelowAll)]), 
                        "Above"=as.double(aboveMaxAny[aboveMaxAny>0][match(combineAboveBelowGroups, groupsAboveAll)])))

ca_df$Below<-as.double(ca_df$Below); ca_df$Above <- as.double(ca_df$Above)
ca_df$Below[is.na(ca_df$Below)]<-0; ca_df$Above[is.na(ca_df$Above)]<-0
ca_df$Total<-as.double(ca_df$Below) + as.double(ca_df$Above)
ca_df[order(ca_df$Total, decreasing = TRUE),]

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
thisIndex<-groupsDF$Code !="DC"
temp<-cbind(aboveMaxByGroup[thisIndex], belowMinByGroup[thisIndex]); colnames(temp)<-c("above","below");  rownames(temp)<-Ntracers[thisIndex]
toPlot<-melt(temp);

par(mar=c(4,4.5,0.5,4))
bp<-ggplot(data = toPlot, aes(x = as.character(X1), fill = X2, y = value)) + 
  geom_bar(stat = 'identity')
thisCols<-c(myLightBlue, "midnightblue")
bp  + coord_flip()   + scale_fill_manual(values=thisCols) + labs(y="Number of models (with one year repeated)", x="Species group") + theme_igray() +
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
