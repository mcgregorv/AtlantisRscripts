#reads in initial conditions and compares biomass by trophic level - over all, and by cell
this_path<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
source(paste(DIR$'General functions',"getSNRN.R",sep=""))
source(paste(DIR$'General functions',"nonZeroMean.R",sep=""))

plotPath<-paste(this_path,"Figures\\Validation\\TL_",sep="")

mg_2_tonne<-2e-8; X_CN<-5.7

##we'll need the groups df 
groupsDF<-read.csv(paste(this_path,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
thisInitialConditionsFile<-paste(this_path,"CRAM_input_short.nc",sep="")
ThisIC.nc<-nc_open(thisInitialConditionsFile)
thisBiolFile<-paste(this_path,"CRAM_base_biol.prm",sep="")
biolLines<-readLines(thisBiolFile)
groupsTL<-read.csv(paste(this_path,"inputs\\biol_prm\\CRAM_trophicLevels_isotopes.csv",sep=""))

ppGroups<-as.character(read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ppGroups.csv",sep=""))[,1]);
npg<-length(ppGroups)

#this was created in ExperimentingWithEatingfromDiets_array... handy, but if don't have it, can just index cohorts for each age by grabing age mature first
ageCohortLinking<-read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ageCohortLinking.csv",sep=""))

#get volume from ICs
volume<-ncvar_get(ThisIC.nc,"volume")[,,1]
nlayers<-dim(volume)[1]; nboxes<-dim(volume)[2]

biomassByAgedGroup<-array(0,dim=c(npg,dim(volume)))

for(g in 1:npg){
  thisVar<-ppGroups[g]
  thisCode<-gsub("juv|ad","",thisVar)
  thisAge<-gsub(thisCode,"",thisVar)
  thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode],side="both")
  thisNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisCode]
  if(thisNumCohorts>1){
    #need to know  cohorts for this age
    temp<-ageCohortLinking[ageCohortLinking[,1]==thisCode,-1]
    cohorts<-seq(1,thisNumCohorts)[temp==thisVar]; cohorts<-cohorts[!is.na(cohorts)]
    for(c in 1:length(cohorts)){
      thisCohort<-cohorts[c]
      #get numbers, resn, structn, then multiply for biomass
      tempvar<-paste(thisName,thisCohort,"_Nums",sep="")
      tempNums<-ncvar_get(ThisIC.nc,tempvar)[,,1]
      ResN<-getSNRN(name=thisName,cohort=thisCohort,whichN = "RN",ThisIC.nc)
      StructN<-getSNRN(name=thisName,cohort=thisCohort,whichN = "SN",ThisIC.nc)
      tempBiomass<-tempNums*(ResN+StructN)
      biomassByAgedGroup[g,,]<-biomassByAgedGroup[g,,]+tempBiomass
    }
  }else{
    tempvar<-paste(thisName,"_N",sep="")
    tempBiomass<-ncvar_get(ThisIC.nc,tempvar)
    if(length(dim(tempBiomass))==3){
      biomassByAgedGroup[g,,]<-tempBiomass[,,1]*volume
    } else{
      biomassByAgedGroup[g,nlayers,]<-tempBiomass[,1]*volume[nlayers,]
    }
  }
}

##sum up by TL
## first index ppGroups by TL
ppTLs<-rep(NA,npg)
for(g in 1:npg){
  thisVar<-ppGroups[g]
  thisCode<-gsub("juv|ad","",thisVar)
  thisAge<-gsub(thisCode,"",thisVar)
  thisTLrow<-groupsTL[groupsTL$Code==thisCode,]
  TL_A<-thisTLrow$Isotope
  if(thisAge=="juv"){
    thisTL<-max(thisTLrow$TrophicLevel1,TL_A-0.5,na.rm=TRUE)
  } else { #either adult or un-aged here
    thisTL<-TL_A
    if(is.na(TL_A)){thisTL<-thisTLrow$TrophicLevel2}
  }
  if(is.na(thisTL)){thisTL<-1}
  ppTLs[g]<-thisTL
}



ppTLs[ppGroups %in% c("PB","BB")]<-0
ppTLs[ppGroups %in% c("DL","DR")]<-0

biomassByppGroup<-apply(biomassByAgedGroup,1,sum,na.rm=TRUE)
biomasTLdf<-data.frame(cbind(ppGroups,biomassByppGroup,ppTLs)); colnames(biomasTLdf)<-c("Group","Biomass","TL")
#take ceiling of TL to use as factor
biomasTLdf$TLc<-unlist(lapply(biomasTLdf$TL,FUN=function(x){floor(as.double(as.character(x)))}))
biomasTLdf$TLc<-as.factor(as.character(biomasTLdf$TLc))
biomasTLdf$Biomass<-as.double(as.character(biomasTLdf$Biomass))
biomassByTL<-tapply(biomasTLdf$Biomass,biomasTLdf$TLc,sum,na.rm=TRUE)*(mg_2_tonne*X_CN)

##take out any zero biomass (DC in CRAM case)
biomasTLdf<-biomasTLdf[biomasTLdf$Biomass>0,]

bacteria<-sum(biomasTLdf$Biomass[biomasTLdf$Group %in% c("PB","BB")])*(mg_2_tonne*X_CN)
detritus<-sum(biomasTLdf$Biomass[biomasTLdf$Group %in% c("DL","DR","DC")])*(mg_2_tonne*X_CN)

#total volume
totalArea<-sum(volume[6,])*1e-6

biomassByKm2<-biomassByTL/totalArea

thisMin<-log10(0.001); thisMax<-log10(max(biomassByTL[-1],c(bacteria,detritus)))
logBiomassByTL<-log10(biomassByTL)
yaxis_at<-seq(0,thisMax,by=2); yaxis=10^(yaxis_at)

decline5perc<-logBiomassByTL[c("1")]-seq(0,4)*0.05
decline10perc<-logBiomassByTL[c("1")]-seq(0,4)*0.1

fitData<-data.frame(cbind("TL"=seq(1,5),"biomass"=logBiomassByTL[-1]))

model<-lm(TL~biomass,data=fitData)
thisInt<-summary(model)$'coefficients'[1,1]
thisSlope<-summary(model)$'coefficients'[2,1]
fitLine<-unlist(lapply(seq(1,5),FUN=function(x){thisSlope*x+thisInt}))

pdf(paste(DIR$'Figures',"Validation\\trophicLevels\\InitialConditions_BiomassByTL.pdf",sep=""),height=4,width=5)
par(las=1,mar=c(4.5,5,1,1))
plot(as.double(names(logBiomassByTL))[-1],type="n",xlab="Trophic level",ylab="",ylim=c(0,thisMax),yaxt="n",xlim=c(0.5,5))
points(x=as.double(names(logBiomassByTL))[-1],y=logBiomassByTL[-1],type="h",lwd=5,lend=1,col=myOrange)
axis(at=yaxis_at,labels=yaxis,side=2)
mtext("Biomass (tonnes)",side=2,adj=0.5,las=0,line=3.5)
points(x=seq(1,5),y=fitLine,type="l",lwd=2,col=myBlue)
points(x=c(0.5,0.7),y=c(log10(detritus),log10(bacteria)),type="h",lwd=5,lend=1,col=myGrey_trans)
mtext(paste("Slope: ",signif(thisSlope,2),sep=""),side=3,adj=0.98,line=-1.1,col=myBlue,font=2)
abline(v=0.90,col="red",lty=2,lwd=2)
par(las=2)
axis(at=c(0.5,0.7),labels=c("Detritus","Bacteria"),side=1)
dev.off()

##plot each functional group, ordered by trophic level
orderedBiomassByTL<-biomasTLdf[order(biomasTLdf$TL),]
index<-orderedBiomassByTL$TL!=0
toPlot<-orderedBiomassByTL[index,]
toPlot$BiomassTonnes<-toPlot$Biomass*mg_2_tonne*X_CN
toPlot$Group<-as.character(toPlot$Group)

thisSeq<-seq(1,dim(toPlot)[1])
TL1index<-thisSeq[toPlot$TLc == "1"]; TL2index<-thisSeq[toPlot$TLc == "2"]
TL3index<-thisSeq[toPlot$TLc == "3"]; TL4index<-thisSeq[toPlot$TLc == "4"];
TL5index<-thisSeq[toPlot$TLc == "5"]


thisMin<-log10(0.001); thisMax<-log10(max(toPlot$BiomassTonnes))
yaxis_at<-seq(0,thisMax,by=2); yaxis=10^(yaxis_at)
ypoly<-c(0,0,thisMax,thisMax)


pdf(paste(DIR$'Figures',"Validation\\trophicLevels\\InitialConditions_BiomassByGroupAndTL.pdf",sep=""),height=4,width=14)
par(las=0,mar=c(5,5,1,1))
plot(x=seq(1,dim(toPlot)[1]),y=log10(toPlot$BiomassTonnes),type="h",lwd=3,col=myBlue,xlab="",ylab="Biomass (tonnes)",xaxt="n",yaxt="n",lend=1)
polygon(x=c(min(TL1index)-0.5,max(TL1index)+0.5,max(TL1index)+0.5,min(TL1index)-0.5),y=ypoly,col=myGreen_trans,border=NA)
polygon(x=c(min(TL2index)-0.5,max(TL2index)+0.5,max(TL2index)+0.5,min(TL2index)-0.5),y=ypoly,col=myBlue_trans,border=NA)
polygon(x=c(min(TL3index)-0.5,max(TL3index)+0.5,max(TL3index)+0.5,min(TL3index)-0.5),y=ypoly,col=myPurple_trans,border=NA)
polygon(x=c(min(TL4index)-0.5,max(TL4index)+0.5,max(TL4index)+0.5,min(TL4index)-0.5),y=ypoly,col=myRed_trans,border=NA)
polygon(x=c(min(TL5index)-0.5,max(TL5index)+0.5,max(TL5index)+0.5,min(TL5index)-0.5),y=ypoly,col=myGrey_trans,border=NA)
axis(at=yaxis_at,labels=yaxis,side=2)
par(las=2)
axis(at=seq(1,dim(toPlot)[1]),labels=toPlot$Group,side=1,cex=0.15)
dev.off()

# 
# ########################################################
# ## ALTS ALTS ##################################################
# #set sediment and pelagic bacteria to TL 1 - based on Michelle's levels for TBGB
# ppTLs[ppGroups %in% c("PB","BB")]<-1
# #what if add DL and DR in too?
# ppTLs[ppGroups %in% c("DL","DR")]<-1
# 
# biomassByppGroup<-apply(biomassByAgedGroup,1,sum,na.rm=TRUE)
# biomasTLdf<-data.frame(cbind(ppGroups,biomassByppGroup,ppTLs)); colnames(biomasTLdf)<-c("Group","Biomass","TL")
# #take ceiling of TL to use as factor
# biomasTLdf$TLc<-unlist(lapply(biomasTLdf$TL,FUN=function(x){floor(as.double(as.character(x)))}))
# biomasTLdf$TLc<-as.factor(as.character(biomasTLdf$TLc))
# biomasTLdf$Biomass<-as.double(as.character(biomasTLdf$Biomass))
# biomassByTL<-tapply(biomasTLdf$Biomass,biomasTLdf$TLc,sum,na.rm=TRUE)
# 
# ##take out any zero biomass (DC in CRAM case)
# biomasTLdf<-biomasTLdf[biomasTLdf$Biomass>0,]
# 
# logBiomass<-log(biomassByTL[-1],base=exp(1)) #missing the zeros (DL, DR)
# fitData<-data.frame(cbind(seq(1,5),logBiomass,biomassByTL[-1])); colnames(fitData)<-c("TL","logBiomass","biomass")
# fitData<-fitData[1:3,]
# 
# model<-lm(logBiomass~TL,data=fitData)
# int<-as.double(model["coefficients"][[1]][1]); slope<-as.double(model["coefficients"][[1]][2])
# 
# fittedLine<-unlist(lapply(seq(1,5),FUN=function(x){int+slope*x}))
# 
# pdf(paste(DIR$'Figures',"Validation\\trophicLevels\\InitialConditions_detritusAndBacteriaLevel1.pdf",sep=""),height=5)
# plot(x=as.double(names(biomassByTL))[-1],y=biomassByTL[-1],type="h",lwd=5,lend=1,col=myOrange,xlab="Trophic level",ylab="Biomass (tonnes)")
# points(exp(fittedLine),type="l",lwd=2,col=myBlue)
# dev.off()
# 
# ########################
# ##without bacteria and detritus in TL1
# ppTLs[ppGroups %in% c("PB","BB")]<-0
# ppTLs[ppGroups %in% c("DL","DR")]<-0
# 
# biomassByppGroup<-apply(biomassByAgedGroup,1,sum,na.rm=TRUE)
# biomasTLdf<-data.frame(cbind(ppGroups,biomassByppGroup,ppTLs)); colnames(biomasTLdf)<-c("Group","Biomass","TL")
# #take ceiling of TL to use as factor
# biomasTLdf$TLc<-unlist(lapply(biomasTLdf$TL,FUN=function(x){floor(as.double(as.character(x)))}))
# biomasTLdf$TLc<-as.factor(as.character(biomasTLdf$TLc))
# biomasTLdf$Biomass<-as.double(as.character(biomasTLdf$Biomass))
# biomassByTL<-tapply(biomasTLdf$Biomass,biomasTLdf$TLc,sum,na.rm=TRUE)
# 
# ##take out any zero biomass (DC in CRAM case)
# biomasTLdf<-biomasTLdf[biomasTLdf$Biomass>0,]
# 
# logBiomass<-log(biomassByTL[-1],base=exp(1)) #missing the zeros (DL, DR)
# fitData<-data.frame(cbind(seq(1,5),logBiomass,biomassByTL[-1])); colnames(fitData)<-c("TL","logBiomass","biomass")
# fitData<-fitData[1:5,]
# 
# model<-lm(biomass~TL,data=fitData)
# int<-as.double(model["coefficients"][[1]][1]); slope<-as.double(model["coefficients"][[1]][2])
# 
# fittedLine<-unlist(lapply(seq(1,5),FUN=function(x){int+slope*x}))
# 
# pdf(paste(DIR$'Figures',"Validation\\trophicLevels\\InitialConditions_detritusAndBacteriaLevel0.pdf",sep=""),height=5)
# plot(x=as.double(names(biomassByTL))[-1],y=biomassByTL[-1],type="h",lwd=5,lend=1,col=myOrange,xlab="Trophic level",ylab="Biomass (tonnes)")
# points((fittedLine),type="l",lwd=2,col=myBlue)
# dev.off()
# 
# ########################
# ##without bacteria and detritus in TL1 - this is what Michelle has for TBGB validation
# ppTLs[ppGroups %in% c("PB","BB")]<-1
# ppTLs[ppGroups %in% c("DL","DR")]<-0
# 
# biomassByppGroup<-apply(biomassByAgedGroup,1,sum,na.rm=TRUE)
# biomasTLdf<-data.frame(cbind(ppGroups,biomassByppGroup,ppTLs)); colnames(biomasTLdf)<-c("Group","Biomass","TL")
# #take ceiling of TL to use as factor
# biomasTLdf$TLc<-unlist(lapply(biomasTLdf$TL,FUN=function(x){floor(as.double(as.character(x)))}))
# biomasTLdf$TLc<-as.factor(as.character(biomasTLdf$TLc))
# biomasTLdf$Biomass<-as.double(as.character(biomasTLdf$Biomass))
# biomassByTL<-tapply(biomasTLdf$Biomass,biomasTLdf$TLc,sum,na.rm=TRUE)
# 
# ##take out any zero biomass (DC in CRAM case)
# biomasTLdf<-biomasTLdf[biomasTLdf$Biomass>0,]
# 
# logBiomass<-log(biomassByTL[-1],base=exp(1)) #missing the zeros (DL, DR)
# fitData<-data.frame(cbind(seq(1,5),logBiomass,biomassByTL[-1])); colnames(fitData)<-c("TL","logBiomass","biomass")
# fitData<-fitData[1:5,]
# 
# model<-lm(biomass~TL,data=fitData)
# int<-as.double(model["coefficients"][[1]][1]); slope<-as.double(model["coefficients"][[1]][2])
# 
# fittedLine<-unlist(lapply(seq(1,5),FUN=function(x){int+slope*x}))
# 
# pdf(paste(DIR$'Figures',"Validation\\trophicLevels\\InitialConditions_detritusLevel0AndBacteriaLevel1.pdf",sep=""),height=5)
# plot(x=as.double(names(biomassByTL))[-1],y=biomassByTL[-1],type="h",lwd=5,lend=1,col=myOrange,xlab="Trophic level",ylab="Biomass (tonnes)")
# points((fittedLine),type="l",lwd=2,col=myRed)
# dev.off()
# 
