#reads in initial conditions and compares biomass by trophic level - over all, and by cell
## TL calculated in calculate_TL_fromDiet.R in R/Growth/
this_path<-paste(DIR$'Base',"TBGB\\",sep="")
source(paste(DIR$'General functions',"getSNRN.R",sep=""))
source(paste(DIR$'General functions',"nonZeroMean.R",sep=""))

plotPath<-paste(DIR$'Figures',"\\Validation\\TLcalc_",sep="")
this_out<-"BaseBurnin35"
outPath<-paste(this_path,"TBGB_JP2\\output",this_out,"\\",sep="")
TLpath<-paste(DIR$'Base',"\\TBGB\\EWE\\", sep="")
TLpath <- outPath

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1880
fishingStartYear<-1900
modelStartYear<-1900

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
ntsteps<-dim(thisVol)[3]

testTS<-75; lastYear <- 1960 # looking at pre-fishing
last_ts <- round((lastYear-year0)*numStepsPerYear)
timeIndex<-(modelStartYear-year0): last_ts


mg_2_tonne<-2e-8; X_CN<-5.7

getAgeMature <- function(Code){
  x <- grep(paste(Code,"_age_mat", sep=""), biolLines)
  y <- get_first_number(biolLines[x], n=1)
  return(y)
}

##we'll need the groups df 
groupsDF<-read.csv(paste(outPath,"..\\TBGB_groups.csv",sep="")); ng<-dim(groupsDF)[1]
thisBiolFile<-paste(outPath,"..\\TBGB_biol.prm",sep="")
biolLines<-readLines(thisBiolFile)
TLfile <- paste(TLpath,"TBGBGroupsTL.csv",sep="")
if(file.exists(TLfile)){
  groupsTL<-read.csv(TLfile)
} else{
  groupsTL<-read.csv(paste(this_path,"EwE\\TBGBGroupsTL.csv",sep=""))
}

ppGroups <- c()
for(g in 1:ng){
  thisNumCohorts <- groupsDF$NumCohorts[g]; thisCode <- groupsDF$Code[g]
  if(thisNumCohorts>1){
    thisAdd <- paste(thisCode,c("ad", "juv"), sep="")
  } else{
    thisAdd <- as.character(thisCode)
  }
  ppGroups <- c(ppGroups, thisAdd)
}
npg <- length(ppGroups)

ageCohortLinking <- array(NA, dim=c(ng, 11))
for(g in 1:ng){
  thisCode <- as.character(groupsDF$Code[g]); thisNumCohorts <- groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    thisAgeMat <- getAgeMature(thisCode)
    ageCohortLinking[g,2:(thisAgeMat+1)]<-paste(thisCode, "juv", sep="")
    ageCohortLinking[g,(thisAgeMat+2):(thisNumCohorts+1)]<-paste(thisCode, "ad", sep="")
  } else{
    ageCohortLinking[g,2]<-as.character(thisCode)
  }
}
ageCohortLinking[,1]<- as.character(groupsDF$Code)

#get volume from tracers
volume<-thisVol[,,timeIndex]
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
      tempNums<-ncvar_get(ThisNC.nc,tempvar)[,,timeIndex]
      #ResN
      tempvar<-paste(thisName,thisCohort,"_ResN",sep="")
      ResN<-ncvar_get(ThisNC.nc,tempvar)[,,timeIndex]
      #StructN
      tempvar<-paste(thisName,thisCohort,"_StructN",sep="")
      StructN<-ncvar_get(ThisNC.nc,tempvar)[,,timeIndex]
      tempBiomass<-tempNums*(ResN+StructN)
      biomassByAgedGroup[g,,,]<-biomassByAgedGroup[g,,,]+tempBiomass
    }
  }else{
    tempvar<-paste(thisName,"_N",sep="")
    tempBiomass<-ncvar_get(ThisNC.nc,tempvar)
    if(length(dim(tempBiomass))==3){
      biomassByAgedGroup[g,,,]<-tempBiomass[,,timeIndex]*volume
    } else{
      biomassByAgedGroup[g,nlayers,,]<-tempBiomass[,timeIndex]*volume[nlayers,,]
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
  TL_A<-thisTLrow$TL
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

biomassByppGroup<-apply(biomassByAgedGroup,c(1,4),sum,na.rm=TRUE)
biomasTLdf<-data.frame(cbind(ppGroups,ppTLs)); colnames(biomasTLdf)<-c("Group","TL")
#take ceiling of TL to use as factor
biomasTLdf$TLc<-unlist(lapply(biomasTLdf$TL,FUN=function(x){round(as.double(as.character(x)))}))
biomasTLdf$TLc<-as.factor(as.character(biomasTLdf$TLc))
# nearest half
biomasTLdf$TL <- as.double(as.character(biomasTLdf$TL))
biomasTLdf$TLpoint5 <- unlist(lapply(biomasTLdf$TL, myRounding, fraction=0.5))


# biomasTLdf$TLc <- biomasTLdf$TLpoint5
allBiomassByTL <- array(NA, dim=c(length(unique(biomasTLdf$TLc)), dim(biomassByAgedGroup)[4]))

for(t in 1:length(timeIndex)){
  thisBioms<-biomassByppGroup[,t]
  thisBiomsByTL <- tapply(thisBioms,biomasTLdf$TLc,sum,na.rm=TRUE)*(mg_2_tonne*X_CN)
  allBiomassByTL[,t]<-thisBiomsByTL
}
thisMax<-max(allBiomassByTL, na.rm=TRUE); thisMin <- min(allBiomassByTL, na.rm=TRUE)
nTLs <- length(unique(biomasTLdf$TLc))
plot(1, type="n", ylim=c(0,thisMax), xlim=c(0,nTLs+0.5))
for(T in 1:nTLs){
  boxplot(allBiomassByTL[T,], at=T, add=TRUE)
}

# write out biomassByppGroup so can use elsewhere
# save(list=c("biomassByppGroup"), file=paste(DIR$'Base',"\\ATLANTISmodels\\ArchivedModels\\Base\\biomassByppGroup", sep=""))

test<-apply(allBiomassByTL, 2, max, na.rm=TRUE)

biomasTLdf$Biomass <- biomassByppGroup[,116]

##take out any zero biomass (DC in CRAM case) - NAs should be OK
takeLog<-function(x){
  y<-NA
  if(x>0){y<-log10(x)}
  return(y)
}
logAllBiomassByTL <- apply(allBiomassByTL, c(1,2), takeLog)
colByTime<-colorRampPalette(colors=c(myLightBlue,myBlue,"midnightblue"))(dim(logAllBiomassByTL)[2])
plot(logAllBiomassByTL[,117], type="l", col=myBlue)
for(t in 1:dim(logAllBiomassByTL)[2]){points(logAllBiomassByTL[,t], type="l", col=colByTime[t])}


bacteriaIndex<-ppGroups %in% c("BB", "PB"); detritusIndex <- ppGroups %in% c("DL", "DR")
bacteria<-apply(biomassByppGroup[bacteriaIndex,], 2, sum, na.rm=TRUE)*(mg_2_tonne*X_CN); detritus <- apply(biomassByppGroup[detritusIndex,], 2, sum, na.rm=TRUE)*(mg_2_tonne*X_CN)
bacteriaLog<-log10(bacteria); detritusLog <- log10(detritus)

test<-sort(bacteriaLog); nb <- length(test); bacteriaMedian <- test[round(nb*0.5)]; bacteriaLC<-test[round(nb*0.025)]; bacteriaUC <- test[round(nb*0.975)]
test<-sort(detritusLog); nb <- length(test); detritusMedian <- test[round(nb*0.5)]; detritusLC<-test[round(nb*0.025)]; detritusUC <- test[round(nb*0.975)]

## fit lines to start and end of timeseries

lm_fn<-function(x,a,b){
  y<-a*x+b
  return(y)
}
ll_lm<-function(pars){
  a<-pars[1]; b<-pars[2]
  est<-unlist(lapply(seq(1,5), lm_fn, a=a, b=b))
  ss<-0
  for(i in 1:length(obs)){
    ss<-ss+ (obs[i]-est[i])^2
  }
  return(ss)
}

fitData<-data.frame(cbind("TL"=seq(1,5),"biomass"=logAllBiomassByTL[-1,1]))
obs<-fitData[,2]
thisOptim<-optim(par=c(10,-5), fn=ll_lm)
optim_a<-thisOptim$par[1]; optim_b<-thisOptim$par[2]
thisInt<-optim_b
thisSlope<-optim_a
fitLine1<-unlist(lapply(seq(1,5),FUN=function(x){thisSlope*x+thisInt}))
thisSlope1<-thisSlope

fitData<-data.frame(cbind("TL"=seq(1,5),"biomass"=logAllBiomassByTL[-1,117]))
obs<-fitData[,2]
thisOptim<-optim(par=c(10,-5), fn=ll_lm)
optim_a<-thisOptim$par[1]; optim_b<-thisOptim$par[2]
thisInt<-optim_b
thisSlope<-optim_a
fitLine2<-unlist(lapply(seq(1,5),FUN=function(x){thisSlope*x+thisInt}))
thisSlope2<-thisSlope

thisMax<-max(logAllBiomassByTL, na.rm=TRUE)
par(las=1,mar=c(4.5,5,1,1))
plot(as.double(1:5),type="n",xlab="Trophic level",ylab="",ylim=c(0,thisMax),yaxt="n",xlim=c(0.5,5))
for(t in 1:dim(logAllBiomassByTL)[2]){points(logAllBiomassByTL[-1,t], pch=15, col=colByTime[t])}
axis(at=yaxis_at,labels=yaxis,side=2)
mtext("Biomass (tonnes)",side=2,adj=0.5,las=0,line=3.5)
points(x=seq(1,5),y=fitLine1,type="l",lwd=2,col=myLightBlue); points(x=seq(1,5),y=fitLine2,type="l",lwd=2,col="midnightblue")
mtext(paste("Slope 1: ",signif(thisSlope1,2),sep=""),side=3,adj=0.98,line=-1.1,col=myLightBlue,font=2)
mtext(paste("Slope 2: ",signif(thisSlope2,2),sep=""),side=3,adj=0.98,line=-2.1,col="midnightblue",font=2)

abline(v=0.90,col="red",lty=2,lwd=2)
par(las=2)
axis(at=c(0.5,0.7),labels=c("Detritus","Bacteria"),side=1)

# grab the median biomass at each trophic level
medBiomByTL<-apply(logAllBiomassByTL, 1, median)
calcQuartile<-function(x, q){
  y<-sort(x); ny<-length(y)
  z <- y[round(ny*q)]
  return(z)
}
# using 95% CIs in place of upper and lower quartiles
LQbiomByTL<-apply(logAllBiomassByTL, 1, calcQuartile, q=0.025); UQbiomByTL <- apply(logAllBiomassByTL, 1, calcQuartile, q=0.975)

fitData<-data.frame(cbind("TL"=seq(1,5),"biomass"=medBiomByTL[-1]))
obs<-fitData[,2]
thisOptim<-optim(par=c(10,-5), fn=ll_lm)
optim_a<-thisOptim$par[1]; optim_b<-thisOptim$par[2]
thisInt<-optim_b
thisSlope<-optim_a
fitLineMed<-unlist(lapply(seq(1,5),FUN=function(x){thisSlope*x+thisInt}))
thisSlopeMed<-thisSlope

## fit slope to LC and UC
fitData<-data.frame(cbind("TL"=seq(1,5),"biomass"=LQbiomByTL[-1]))
obs<-fitData[,2]
thisOptim<-optim(par=c(10,-5), fn=ll_lm)
optim_a<-thisOptim$par[1]; optim_b<-thisOptim$par[2]
thisInt<-optim_b
thisSlope<-optim_a
fitLineMed<-unlist(lapply(seq(1,5),FUN=function(x){thisSlope*x+thisInt}))
thisSlopeLQ<-thisSlope
#UC
fitData<-data.frame(cbind("TL"=seq(1,5),"biomass"=UQbiomByTL[-1]))
obs<-fitData[,2]
thisOptim<-optim(par=c(10,-5), fn=ll_lm)
optim_a<-thisOptim$par[1]; optim_b<-thisOptim$par[2]
thisInt<-optim_b
thisSlope<-optim_a
fitLineMed<-unlist(lapply(seq(1,5),FUN=function(x){thisSlope*x+thisInt}))
thisSlopeUQ<-thisSlope


pdf(paste(DIR$'Figures',"Validation\\trophicLevels\\TS_BiomassByTL_TBGB",this_out,".pdf",sep=""),height=4,width=5)
par(las=1,mar=c(4.5,5,1,1))
plot(seq(1,5),type="n",xlab="Trophic level",ylab="",ylim=c(0,thisMax),yaxt="n",xlim=c(0.5,5))
points(x=seq(1,5),y=medBiomByTL[-1],type="h",lwd=5,lend=1,col=myOrange)
axis(at=yaxis_at,labels=yaxis,side=2)
mtext("Biomass (tonnes)",side=2,adj=0.5,las=0,line=3.5)
points(x=seq(1,5),y=fitLineMed,type="l",lwd=2,col=myBlue)
for(T in 1:5){
  test<-sort(logAllBiomassByTL[(T+1),]); nb<-length(test); LC<-test[(round(nb*0.025))]; UC<-test[(round(nb*0.975))]
  y0=LC; y1=UC
  points(x=c(T,T), y=c(y0,y1), type="l", lty=1, lwd=1.2)
  segments(x0=T-0.05, y0=y0, x1=T+0.05, y1=y0, lwd=1.2); segments(x0=T-0.05, y0=y1, x1=T+0.05, y1=y1, lwd=1.2)
}
points(x=c(0.5,0.7),y=c(detritusMedian,bacteriaMedian),type="h",lwd=5,lend=1,col=myGrey_trans)
mtext(paste("Slope: ",signif(thisSlopeMed,2),sep=""),side=3,adj=0.98,line=-1.1,col=myBlue,font=2)
mtext(paste("(",signif(thisSlopeLQ,2),", ", signif(thisSlopeUQ,2),")",sep=""),side=3,adj=0.98,line=-2,col="black",font=2)
#detritus CIs

dev.off()

