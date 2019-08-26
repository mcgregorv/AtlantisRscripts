## read in the initial conditions for each species group, and temperature from the ROMS
## to estimate the Tcorr affect on mL for each species group - so know how much background M is being applied
inputsPath <- paste(DIR$'Base', "ATLANTISmodels\\inputs\\", sep="")
inputsNC.nc <- nc_open(paste(inputsPath,"..\\CRAM_input_fromBase50yr.nc", sep=""))
tempNC.nc <- nc_open(paste(inputsPath, "ROMS\\Chatham30_tempAll.nc", sep=""));
groupsDF <- read.csv(paste(inputsPath,"..\\CRAM_Groups.csv", sep="")); ng <- dim(groupsDF)[1]
tempData <- ncvar_get(tempNC.nc, "temperature"); nts <- dim(tempData)[3]
biolLines <- readLines( paste(inputsPath,"..\\CRAM_BH_hybrid_biol.prm", sep=""))
nlayers <- dim(tempData)[1]; nboxes <- dim(tempData)[2]

calc_tempScalar <- function(temperature){
  y <- 2^((temperature - 15)/10)
  return(y)
}
get_age_mat<- function(thisCode){
  thisVar<-paste(thisCode,"_age_mat", sep="")
  y<-grep(thisVar,biolLines); z<-get_first_number(biolLines[y])
  return(z)
}
get_mL <- function(thisCode){
  thisVar<-paste(thisCode,"_mL", sep="")
  y<-grep(thisVar,biolLines); z<-get_first_number(biolLines[y+1], n="all")[1:2]
  return(z)
}
getNightOrDay <- function(thisCode){
  x=grep(paste("flag",thisCode, "day", sep=""), biolLines)
  y <- get_first_number(biolLines[x])
  pref<- c("night", "day", "both")[y+1]
  return(pref)
}
zero2NA <- function(x){
  y<-x
  if(x <= 1e-16){y <- NA}
  return(y)
}
calcTempFromTcorr <- function(Tcorr){
  T <- 10 *( (log(Tcorr))/(log(2)) ) + 15
  return(T)
}
getValuesWeighted <- function(scalars, values, k=2){
  scalarSum<-sum(scalars, na.rm=TRUE)
  scalarProps <- scalars/scalarSum
  scalarWholes <- round(scalarProps, k) * 10^k; I<-dim(scalarWholes)[1]; J <- dim(scalarWholes)[2]; 
  thisString <- c()
  for(i in 1:I){
    for(j in 1:J){
      thisMult <- scalarWholes[i,j]; 
      if(!is.na(thisMult)){
        thisTemp <- values[i,j,t]
        if(sum(c(thisMult, thisTemp), na.rm=TRUE)>0){
          thisString <- c(thisString, rep(thisTemp, thisMult))
        }
      }
    }
  }
  return(thisString)
}
getValuesWeighted_nts <- function(scalars, values, k=2, delta_t = 15){
    scalarSum<-sum(scalars, na.rm=TRUE)
    scalarProps <- scalars/scalarSum
    nts <- dim(values)[3]
    timeIndex <- seq(1,nts, by=delta_t)
    thisValues <- values[,,timeIndex]
    scalarWholes <- round(scalarProps, k) * 10^k; I<-dim(scalarWholes)[1]; J <- dim(scalarWholes)[2]; 
    thisString <- c()
    for(i in 1:I){
      for(j in 1:J){
        thisMult <- scalarWholes[i,j]; 
        if(!is.na(thisMult)){
          for(t in 1:dim(thisValues)[3]){
            thisTemp <- thisValues[i,j,t]
            if(sum(c(thisMult, thisTemp), na.rm=TRUE)>0){
              thisString <- c(thisString, rep(thisTemp, thisMult))
            }
          }
        }
      }
    }
  return(thisString)
}
# calcQuartile<-function(scalars, values, k=2, q=0.25, timesteps=10){
#   #scalars are light weightings or number of repetitions - in this case, the density of the species group
#   # values are what we want the quartiles of - in this case, the temperatures 
#   # k is the number of significant figures for the scalars
#   #turn the scalars into proportions then whole numbers
#   thisString <- getValuesWeighted(scalars = scalars, values = values, k=k)
#   sortedString <- sort(thisString); ns <- length(thisString)
#   lq_x <- round(ns*q); uq_x <- (1-q)* ns
#   lq <- sortedString[lq_x]; uq <- sortedString[uq_x]
#   thisMed <- sortedString[0.5*ns]
#   return(c(lq, thisMed, uq))
# }

calcQuartile_1timestep<-function(scalars, values, k=2, q=0.25){
  #scalars are light weightings or number of repetitions - in this case, the density of the species group
  # values are what we want the quartiles of - in this case, the temperatures 
  # k is the number of significant figures for the scalars
  #turn the scalars into proportions then whole numbers
  thisString <- getValuesWeighted(scalars = scalars, values = values, k=k)
  sortedString <- sort(thisString); ns <- length(thisString)
  lq_x <- round(ns*q); uq_x <- (1-q)* ns
  lq <- sortedString[lq_x]; uq <- sortedString[uq_x]
  thisMed <- sortedString[0.5*ns]
  return(c(lq, thisMed, uq))
}

g = grep("HAK", groupsDF$Code)
nts<-dim(tempData)[3]
storeApplied_mLStats<-array(NA, dim=c(ng, 10, nts, 6))
store_scaledTcorr_ALL<-array(NA, dim=c(ng, 10, nlayers*nboxes* nts))
storeApplied_mL_Summary <- array(NA, dim=c(ng, 10,6))
statDesc<-c("mean", "median","LQ", "UQ", "LCI", "UCI")

for(g in 1:ng){
  thisCode<-groupsDF$Code[g];      cat(as.character(thisCode),"--")
  thisNumCohorts <- groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    thisName <- str_trim(groupsDF$Name[g], side="both")
    thisPref<-getNightOrDay(thisCode)
    this_age_mat <- get_age_mat(thisCode)
    this_mLs <- get_mL(thisCode)
    for(c in 1:thisNumCohorts){
      this_mL <- this_mLs[1]
      if(c > this_age_mat){
        this_mL <- this_mLs[2]
      }
      if(this_mL >0){
       thisTracer <- paste(thisName, c, "_Nums", sep="")
        thisData <- ncvar_get(inputsNC.nc, thisTracer)[,,1]
        thisDataProps <- thisData / sum(thisData, na.rm=TRUE)
        thisMeanByTime <- rep(NA, nts)
        nonZeroData <- apply(thisData, c(1,2), zero2NA)
        # thisTempQs <- calcQuartile(scalars = nonZeroData, values=tempData, k=2, q=0.025, timesteps = 10)
         ## 95 percentailes
        thisLQ<-rep(NA, nts); thisMed <- rep(NA, nts); thisUQ <- rep(NA, nts)
        # for(t in 1:nts){
        #   thisTempQs <- calcQuartile_1timestep(scalars = nonZeroData, values=tempData[,,t], k=2, q=0.025)
        #   thisLQ[t]<- thisTempQs[1]; thisMed[t] <- thisTempQs[2]; thisUQ[t]<-thisTempQs[3]
        #   thisMeanByTime[t] <- sum(nonZeroData * tempData[,,t], na.rm=TRUE)/sum(nonZeroData, na.rm=TRUE)
        # }
        thisTempsWeighted <- getValuesWeighted_nts(scalars = nonZeroData, values= tempData, k=2, delta_t=15)
        thisTempScalars <- unlist(lapply(thisTempsWeighted, calc_tempScalar))
        store_scaledTcorr_ALL[g,c,1:length(thisTempScalars)]<- thisTempScalars
        
        # thisQs <- calcQuartile(scalars=nonZeroData, values = tempData, k=2, q=0.25, timesteps=10)
        # thisQs_tempScalar <- unlist(lapply(thisQs, calc_tempScalar))
        # storeApplied_mL_Summary[g,c,grep("median",statDesc)]<- thisQs_tempScalar[2]*this_mL
        # storeApplied_mL_Summary[g,c,grep("LQ",statDesc)]<- thisQs_tempScalar[1]*this_mL
        # storeApplied_mL_Summary[g,c,grep("UQ",statDesc)]<- thisQs_tempScalar[3]*this_mL
        
        # thisTcorrByTime_mean <- unlist(lapply(thisMeanByTime, calc_tempScalar))
        # thisTcorrByTime_med <- unlist(lapply(thisMed, calc_tempScalar))
        # thisTcorrByTime_LQ <- unlist(lapply(thisLQ, calc_tempScalar))
        # thisTcorrByTime_UQ <- unlist(lapply(thisUQ, calc_tempScalar))
        # # store them
        # storeApplied_mLStats[g,c, ,(grep("median",statDesc))]<- thisTcorrByTime_med
        # storeApplied_mLStats[g,c, ,(grep("mean",statDesc))]<- thisTcorrByTime_mean
        # storeApplied_mLStats[g,c, ,(grep("LQ",statDesc))]<- thisTcorrByTime_LQ
        # storeApplied_mLStats[g,c, ,(grep("UQ",statDesc))]<- thisTcorrByTime_UQ
        # # thisYmax <- max(thisTcorrByTime_UQ, na.rm=TRUE); thisYmin <- min(thisTcorrByTime_LQ, na.rm=TRUE)
        # # plot(thisTcorrByTime_med, type="l", ylim=c(thisYmin, thisYmax))
        # # polygon(y=c(thisTcorrByTime_LQ, rev(thisTcorrByTime_UQ)), x=c(1:nts, rev(1:nts)), col=myBlue_trans, border=NA)
        # # do again for quartiles
        # thisLQ<-rep(NA, nts); thisMed <- rep(NA, nts); thisUQ <- rep(NA, nts)
        # for(t in 1:nts){
        #   thisTempQs <- calcQuartile_1timestep(scalars = nonZeroData, values=tempData[,,t], k=2, q=0.25)
        #   thisLQ[t]<- thisTempQs[1];  thisUQ[t]<-thisTempQs[3]
        # }
        # 
        # thisTcorrByTime_LQ <- unlist(lapply(thisLQ, calc_tempScalar))
        # thisTcorrByTime_UQ <- unlist(lapply(thisUQ, calc_tempScalar))
        # # store these
        # storeApplied_mLStats[g,c, ,(grep("LCI",statDesc))]<- thisTcorrByTime_LQ
        # storeApplied_mLStats[g,c, ,(grep("UCI",statDesc))]<- thisTcorrByTime_UQ
        # polygon(y=c(thisTcorrByTime_LQ, rev(thisTcorrByTime_UQ)), x=c(1:nts, rev(1:nts)), col=myBlue_trans, border=NA)
      }
    }
  }
}

### UP TO HERE UP TO HERE
## TURN store_scaledTcorr_ALL INTO M ## for all groups juv and adult - by scaled, I mean weighted

## turn corrections into mL applied and M
calcMfrom_mL<-function(mL){
  M <- (-365)* log(1-mL)
  return(M)
}
storeMboxplots_juv<-NULL; storeMboxplots_ad<-NULL; storemLboxplots_juv<-NULL; storemLboxplots_ad<-NULL
for(g in 1:ng){
  thisCode<-groupsDF$Code[g];      cat(as.character(thisCode),"--")
  thisNumCohorts <- groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    thisName <- str_trim(groupsDF$Name[g], side="both")
    thisPref<-getNightOrDay(thisCode)
    this_age_mat <- get_age_mat(thisCode)
    this_mLs <- get_mL(thisCode)
    for(c in c(1,(this_age_mat+1))){
      this_mL <- this_mLs[1]
      if(c > this_age_mat){
        this_mL <- this_mLs[2]
      }
      if(this_mL >0){
        thisTracer <- paste(thisName, c, "_Nums", sep="")
        # thisTcorrs <- storeApplied_mLStats[g,c,,]
        thisTcorrs <- store_scaledTcorr_ALL[g,c,]
        thisScaled_mL <- this_mL * thisTcorrs
        if(thisPref %in% c("night","day")){thisScaled_mL <- 0.5*this_mL * thisTcorrs}
        thisM <- unlist(lapply(thisScaled_mL, calcMfrom_mL))
        thisMstats <- boxplot.stats(thisM)$stats
        thisScaledmLstats <- boxplot.stats(thisScaled_mL)$stats
        if(c > this_age_mat){
          storeMboxplots_ad[[as.character(thisCode)]]<-thisMstats
          storemLboxplots_ad[[as.character(thisCode)]]<-thisScaledmLstats
        } else{
          storeMboxplots_juv[[as.character(thisCode)]]<-thisMstats
          storemLboxplots_juv[[as.character(thisCode)]]<-thisScaledmLstats
        }
      }else{
        if(c > this_age_mat){
          storeMboxplots_ad[[as.character(thisCode)]]<-rep(0,5)
          storemLboxplots_ad[[as.character(thisCode)]]<-rep(0,5)
        } else{
          storeMboxplots_juv[[as.character(thisCode)]]<-rep(0,5)
          storemLboxplots_juv[[as.character(thisCode)]]<-rep(0,5)
        }
        
      }
    }
  }
}
## dump store_scaledTcorr_ALL as it's large, so don't need to keep it here
# save("store_scaledTcorr_ALL", file=paste(inputsPath,"..\\storeChaosTemp\\store_scaledTcorr_ALL", sep=""))
# rm(store_scaledTcorr_ALL)

makeBoxPlotFromStats <- function(stats,at,col,border, width=1){
  x0<-at-0.5*width;  x1<-at+0.5*width
  y0<-stats[2]; y1<-stats[4]
  polygon(x=c(x0,x0,x1,x1), y=c(y0,y1,y1,y0), col=col, border=border)
  segments(x0=at, x1=at, y0=stats[1], y1=stats[2], col=border)
  segments(x0=at, x1=at, y0=stats[4], y1=stats[5], col=border)
  segments(x0=x0, x1=x1, y0=stats[3], y1=stats[3], col=border, lwd=2.5)
  x2 <- at-0.25*width;  x3<-at+0.25*width
  segments(x0=x2, x1=x3, y0=stats[5], y1=stats[5], col=border)
  segments(x0=x2, x1=x3, y0=stats[1], y1=stats[1], col=border)
}

##also read in M from literature to compare
plotCodes <- groupsDF$Code[groupsDF$NumCohorts>1]; nplot<-length(plotCodes)
groupBioPars<-read.csv(paste(inputsPath,"..\\CRAM_B0.csv", sep=""))
groupsDFPaper<-read.csv(paste(inputsPath,"..\\CRAM_groupsPaper.csv", sep=""))
groupsDFOrderPaper <- groupsDFPaper[order(groupsDFPaper$Name),]
paperIndex <- groupsDFOrderPaper$Code %in% plotCodes
par(mar=c(10,4,1,1))
plot(x=1:nplot, y=rep(1, nplot), type="n", ylim=c(0,3.5), ylab="Applied M / Literature M", xlab="", xaxt="n")
abline(h=1, col="red", lty=3)
par(las=2)
axis(labels=gsub("_", " ",groupsDFOrderPaper$Name[paperIndex]), at=1:nplot, side=1)
for(g in 1:nplot){
  thisCode<-as.character(groupsDFOrderPaper$Code[paperIndex][g])
  thisGindex <- grep(thisCode, groupsDF$Code)
  thisLitM <-  groupBioPars$M[groupBioPars$Code==thisCode]
  thisAdstats <- storeMboxplots_ad[[thisCode]]/thisLitM; thisJuvstats <- storeMboxplots_juv[[thisCode]]/thisLitM
  makeBoxPlotFromStats(stats=thisAdstats, at=g+0.1, col=myBlue_trans, border=myBlue, width=0.5)
  makeBoxPlotFromStats(stats=thisJuvstats, at=g-0.1, col=myOrange_trans, border=myOrange, width=0.5)
} 
 
  
#     thisMads <- as.vector(storeAll_appliedMad[thisGindex,adIndex,,k])
#   thisMjuvs <- as.vector(storeAll_appliedMjuv[thisGindex,juvIndex,,k])
#   thisLitM <-  groupBioPars$M[groupBioPars$Code==thisCode]                 
#   
#   if(sum(thisMads, na.rm=TRUE)>0){
#     boxplot(thisMads/thisLitM, outline=FALSE, col=myBlue_trans, border=myBlue, add=TRUE, at=(g-0.1))
#   }
#   if(sum(thisMjuvs, na.rm=TRUE)>0){
#     boxplot(thisMjuvs/thisLitM, outline=FALSE, col=myOrange_trans, border=myOrange, add=TRUE, at=(g+0.1))
#   }
# }

## version using median, LQ, ect
# storeAll_appliedmLjuv<- 0 * storeApplied_mLStats; storeAll_appliedMjuv <- 0 * storeApplied_mLStats
# storeAll_appliedmLad<- 0 * storeApplied_mLStats; storeAll_appliedMad <- 0 * storeApplied_mLStats
# for(g in 1:ng){
#   thisCode<-groupsDF$Code[g];      cat(as.character(thisCode),"--")
#   thisNumCohorts <- groupsDF$NumCohorts[g]
#   if(thisNumCohorts>1){
#     thisName <- str_trim(groupsDF$Name[g], side="both")
#     thisPref<-getNightOrDay(thisCode)
#     this_age_mat <- get_age_mat(thisCode)
#     this_mLs <- get_mL(thisCode)
#     for(c in 1:thisNumCohorts){
#       this_mL <- this_mLs[1]
#       if(c > this_age_mat){
#         this_mL <- this_mLs[2]
#       }
#       if(this_mL >0){
#         thisTracer <- paste(thisName, c, "_Nums", sep="")
#         thisTcorrs <- storeApplied_mLStats[g,c,,]
#         thisScaled_mL <- this_mL * thisTcorrs
#         if(thisPref %in% c("night","day")){thisScaled_mL <- 0.5*this_mL * thisTcorrs}
#         thisM <- apply(thisScaled_mL,c(1,2), calcMfrom_mL)
#         if(c > this_age_mat){
#           storeAll_appliedmLad[g,c,,] <- thisScaled_mL 
#           storeAll_appliedMad[g,c,,] <- thisM
#         } else{
#           storeAll_appliedmLjuv[g,c,,] <- thisScaled_mL 
#           storeAll_appliedMjuv[g,c,,] <- thisM
#         }
#       }
#     }
#   }
# }

#Index those to plot
# testSumAd <- apply(storeAll_appliedMad, 1, sum, na.rm=TRUE)
# testSumJuv <- apply(storeAll_appliedMjuv, 1, sum, na.rm=TRUE)
# plotIndex <- (testSumAd + testSumJuv) >0; 
# plotCodes <- as.character(groupsDF$Code[plotIndex]); nplot <- length(plotCodes)

# ##also read in M from literature to compare
# groupBioPars<-read.csv(paste(inputsPath,"..\\CRAM_B0.csv", sep=""))
# groupsDFPaper<-read.csv(paste(inputsPath,"..\\CRAM_groupsPaper.csv", sep=""))
# groupsDFOrderPaper <- groupsDFPaper[order(groupsDFPaper$Name),]
# paperIndex <- groupsDFOrderPaper$Code %in% plotCodes
# k <- grep("median|LQ|UQ", statDesc)
# par(mar=c(10,4,1,1))
# plot(x=1:nplot, y=rep(1, nplot), type="n", ylim=c(0,3), ylab="Applied M/Literature M", xlab="", xaxt="n")
# abline(h=1, col="red", lty=3)
# par(las=2)
# axis(labels=gsub("_", " ",groupsDFOrderPaper$Name[paperIndex]), at=1:nplot, side=1)
# for(g in 1:nplot){
#   thisCode<-groupsDFOrderPaper$Code[paperIndex][g]
#   thisGindex <- grep(thisCode, groupsDF$Code)
#   thisNumCohorts <- groupsDF$NumCohorts[thisGindex]
#   this_age_mat <- get_age_mat(thisCode)
#   adIndex <- (this_age_mat+1):thisNumCohorts; juvIndex <- 1:this_age_mat
#   # alt - just use first cohort as they are the same anyway
#   adIndex <- this_age_mat+1; juvIndex <- 1
#   thisMads <- as.vector(storeAll_appliedMad[thisGindex,adIndex,,k])
#   thisMjuvs <- as.vector(storeAll_appliedMjuv[thisGindex,juvIndex,,k])
#   thisLitM <-  groupBioPars$M[groupBioPars$Code==thisCode]                 
#   
#   if(sum(thisMads, na.rm=TRUE)>0){
#     boxplot(thisMads/thisLitM, outline=FALSE, col=myBlue_trans, border=myBlue, add=TRUE, at=(g-0.1))
#   }
#   if(sum(thisMjuvs, na.rm=TRUE)>0){
#     boxplot(thisMjuvs/thisLitM, outline=FALSE, col=myOrange_trans, border=myOrange, add=TRUE, at=(g+0.1))
#   }
# }
############################################
## need realised M to compare
########################################
##Read in tracers and summarise how much mortality each group has (estimated from numbers at age) and compare this with M from mL
mg_2_tonne<-2e-8; X_CN<-5.7

#################################### TRACERS
#plot all tracers for a given box and layer
this_run<-"base"

this_out<-"BASE"
burnin<-35 #number of years to skip in plot
# burnin<-1
this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")

doingMedian<-TRUE
doingMedian<-FALSE

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1900
fishingStartYear<-1900
modelStartYear<-1900

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
ntsteps<-dim(thisVol)[3]

nts<-dim(thisVol)[3]-burnin #number of timesteps
cat(paste(nts,"\n"))
xLabsTemp<-seq(0,(nts*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0

#get all tracer names
allTracers<-sort(names(ThisNC.nc$var))

#grab numbers tracers and store
storeNumbersByGroup<-array(NA, dim=c(nts, ng, 10))
for(g in 1:ng){
  thisName<-str_trim(groupsDF$Name[g], side="both"); thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    for(c in 1:thisNumCohorts){
      thisTracer<-paste(thisName,c,"_Nums", sep="")
      temp<-ncvar_get(ThisNC.nc, thisTracer)
      thisNumbers<-apply(temp,3, sum)
      storeNumbersByGroup[,g,c]<-thisNumbers[(burnin+1):length(thisNumbers)]
    }
  }
}
##get proportions at age for each timestep and group, so can get 95% CIs
storeAllPropsAtAge<-array(NA,dim=c(nts,ng,10))
for(g in 1:ng){
  thisNumCohorts<-groupsDF$NumCohorts[g]
  thisNumByTS<-apply(storeNumbersByGroup[,g,],1,sum,na.rm=TRUE)
  if(thisNumCohorts>2){
    for(t in 1:nts){
      storeAllPropsAtAge[t,g,1:thisNumCohorts]<-storeNumbersByGroup[t,g,1:thisNumCohorts]/thisNumByTS[t]
    }
  }
}

numsICFromM_fn<-function(M){
  thisNumFromM<-rep(NA,length(thisNumbers)); thisNumFromM[1]<-thisNumbers[1]
  for(i in 2:length(thisNumbers)){
    thisNumFromM[i]<-(thisNumFromM[i-1])*exp(-M*ageClassSize)
  }
  return(thisNumFromM)
}
ll_IC_fn<-function(pars){
  optimM<-pars[1]
  estNums<-numsICFromM_fn(optimM)
  ss<-0
  for(c in 1:length(thisNumbers)){
    ss<-ss+abs(thisNumbers[c]-estNums[c])
  }
  return(ss)
}

calc_P1<-function(M, thisNumCohorts, ageClassSize){
  #what should the first proportion be for all to add to 1, for a given number of age classes and a given M?
  xx=0
  for(i in 1:thisNumCohorts){
    xx = xx + exp(-(i-1) * M * ageClassSize)
  }
  p1<-1/xx
  return(p1)
}
getProportionsAtAge<-function(M, thisNumCohorts, ageClassSize){
  thisP1<-calc_P1(M, thisNumCohorts, ageClassSize)
  propsAtAge<-rep(NA, thisNumCohorts); propsAtAge[1]<-thisP1
  for(i in 2:thisNumCohorts){
    propsAtAge[i]<-propsAtAge[i-1]* exp(-(M * ageClassSize))
  }
  return(propsAtAge)
}
#now loop through and estimate M 
storeJuvMs<-rep(NA,ng); storeAdMs<-rep(NA,ng)
for(g in 1:ng){
  thisName<-str_trim(groupsDF$Name[g], side="both"); thisNumCohorts<-groupsDF$NumCohorts[g]; thisCode<-as.character(groupsDF$Code[g])
  if(thisNumCohorts>2){
    ## define props at age using M ('obs')
    thisM<-groupBioPars$M[groupBioPars$Code==thisCode]
    if(!is.na(thisM)){
      thisVar<-paste(thisCode, "_AgeClassSize", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageClassSize<-get_first_number(temp)
      thisVar<-paste(thisCode, "_age_mat", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageMature<-get_first_number(temp)
      # 'observed' props at age using M from lit
      obsPropsAtAge<-getProportionsAtAge(thisM, thisNumCohorts, ageClassSize)
     
      thisLongName<-gsub("_"," ", groupsDFPaper$Name[g])
      # thisLongName<-gsub("-","\n",thisLongName)
      # thisLongName<-gsub("\\(","\n",thisLongName); thisLongName<-gsub(")", "", thisLongName)
      
      #get props at age CIs
      propLowerCI<-rep(NA,thisNumCohorts); propUpperCI<-rep(NA,thisNumCohorts); propMedian<-rep(NA,thisNumCohorts)
      for(c in 1:thisNumCohorts){
        test<-sort(storeAllPropsAtAge[,g,c])
        propLowerCI[c]<-test[round(0.025*length(test))]; propUpperCI[c]<-test[round(0.975*length(test))]; propMedian[c]<-test[round(0.5*length(test))]
      }
      thisData<-storeNumbersByGroup[,g,]
      totalNumsTracer<-apply(thisData,1,sum, na.rm=TRUE)
      propsByCohortTime<- 0 * thisData
      for(t in 1:dim(thisData)[1]){propsByCohortTime[t,]<-thisData[t,]/totalNumsTracer[t]}
      
      thisProps<-thisData[t,]/totalNumsTracer[t]
      # if(doingMedian==TRUE){thisProps<-propMedian}
      #do juveniles first
      thisNumbers<-thisProps[1:(ageMature+1)]
      juvFittedM<-optim(par=0.1,fn=ll_IC_fn, method="Brent", upper=2, lower=0)$par
      juvFittedNums<-numsICFromM_fn(juvFittedM)
      ## do adults
      thisNumbers<-thisProps[(ageMature+1):thisNumCohorts]
      adFittedM<-optim(par=0.1,fn=ll_IC_fn, method="Brent", upper=2, lower=0)$par
      adFittedNums<-numsICFromM_fn(adFittedM)
      storeAdMs[g]<-adFittedM; storeJuvMs[g]<-juvFittedM
    }
  }
  
}

### now plot boxplots relative to realised M

par(mar=c(10,4,1,1))
plot(x=1:nplot, y=rep(1, nplot), type="n", ylim=c(0,1.5), ylab="Applied M / Model M", xlab="", xaxt="n")
abline(v=seq(0,nplot), col=myGrey_trans)
abline(h=seq(0,2,by=0.125), col=myGrey_trans)
abline(h=1, col="red", lty=2, lwd=1.5)
par(las=2)
axis(labels=gsub("_", " ",groupsDFOrderPaper$Name[paperIndex]), at=1:nplot, side=1)
for(g in 1:nplot){
  thisCode<-as.character(groupsDFOrderPaper$Code[paperIndex][g])
  thisGindex <- grep(thisCode, groupsDF$Code)
  thisModelMjuv <-  storeJuvMs[grep(thisCode, groupsDF$Code)]; thisModelMad <-  storeAdMs[grep(thisCode, groupsDF$Code)]      
  thisAdstats <- storeMboxplots_ad[[thisCode]]/thisModelMad; thisJuvstats <- storeMboxplots_juv[[thisCode]]/thisModelMjuv
  if(sum(storeMboxplots_ad[[thisCode]])==0){thisAdstats<-rep(0,5)}
  if(sum(storeMboxplots_juv[[thisCode]])==0){thisJuvstats<-rep(0,5)}
  makeBoxPlotFromStats(stats=thisJuvstats, at=g, col=myOrange_trans, border=myOrange, width=0.5)
  makeBoxPlotFromStats(stats=thisAdstats, at=g, col=myBlue_trans, border=myBlue, width=0.5)
}  

# do hist of median proportions for adults, and for females
medianProps_ad <- rep(0, nplot); medianProps_juv <- medianProps_ad
for(g in 1:nplot){
  thisCode<-as.character(groupsDFOrderPaper$Code[paperIndex][g])
  thisGindex <- grep(thisCode, groupsDF$Code)
  thisModelMjuv <-  storeJuvMs[grep(thisCode, groupsDF$Code)]; thisModelMad <-  storeAdMs[grep(thisCode, groupsDF$Code)]      
  thisAdstats <- storeMboxplots_ad[[thisCode]]/thisModelMad; thisJuvstats <- storeMboxplots_juv[[thisCode]]/thisModelMjuv
  if(sum(storeMboxplots_ad[[thisCode]])==0){thisAdstats<-rep(0,5)}
  if(sum(storeMboxplots_juv[[thisCode]])==0){thisJuvstats<-rep(0,5)}
  medianProps_ad[g]<-thisAdstats[3]; medianProps_juv[g]<-thisJuvstats[3]
}
medianProps_ad[is.na(medianProps_ad)]<-0; medianProps_juv[is.na(medianProps_juv)]<-0
hist(medianProps_ad, col=myBlue_trans, border=myBlue, xlab="Applied M / Model M", main="")
hist(medianProps_juv, col=myOrange_trans,  border=myOrange, xlab="Applied M / Model M", main="")

## store values for these plots, so can re-read them to re-do the plots as/when needed without re-doing the analyses
plotValues <- c("medianProps_ad", "medianProps_juv", "storeMboxplots_ad", "storeMboxplots_juv", 
                "storeMboxplots_ad", "storeMboxplots_juv", "storeAdMs", "storeJuvMs")

# save(list=plotValues, file=paste(inputsPath, "..\\storeChaosTemp\\storeValuesToPlot",sep=""))
# load(paste(inputsPath, "..\\storeChaosTemp\\storeValuesToPlot",sep=""))
# 
# k <- grep("median|LQ|UQ", statDesc)
# k <- grep("median", statDesc)
# par(mar=c(10,4,1,1))
# plot(x=1:nplot, y=rep(1, nplot), type="n", ylim=c(0,1.25), ylab="Applied M / Model M", xlab="", xaxt="n")
# abline(v=seq(0,nplot), col=myGrey_trans)
# abline(h=seq(0,2,by=0.125), col=myGrey_trans)
# abline(h=1, col="red", lty=2, lwd=1.5)
# par(las=2)
# axis(labels=gsub("_", " ",groupsDFOrderPaper$Name[paperIndex]), at=1:nplot, side=1)
# for(g in 1:nplot){
#   thisCode<-groupsDFOrderPaper$Code[paperIndex][g]
#   thisGindex <- grep(thisCode, groupsDF$Code)
#   thisNumCohorts <- groupsDF$NumCohorts[thisGindex]
#   this_age_mat <- get_age_mat(thisCode)
#   adIndex <- (this_age_mat+1):thisNumCohorts; juvIndex <- 1:this_age_mat
#   # alt - just use first cohort as they are the same anyway
#   adIndex <- this_age_mat+1; juvIndex <- 1
#   thisMads <- as.vector(storeAll_appliedMad[thisGindex,adIndex,,k])
#   thisMjuvs <- as.vector(storeAll_appliedMjuv[thisGindex,juvIndex,,k])
#   thisLitMjuv <-  storeJuvMs[grep(thisCode, groupsDF$Code)]                 
#   thisLitMad <-  storeAdMs[grep(thisCode, groupsDF$Code)]                 
#   
#   if(sum(thisMads, na.rm=TRUE)>0){
#     boxplot(thisMads/thisLitMad, outline=FALSE, col=myBlue_trans, border=myBlue, add=TRUE, at=(g-0.1))
#   }else{
#     boxplot(0, outline=FALSE, col=myBlue_trans, border=myBlue, add=TRUE, at=(g-0.1))
#   }
#   if(sum(thisMjuvs, na.rm=TRUE)>0){
#     boxplot(thisMjuvs/thisLitMjuv, outline=FALSE, col=myOrange_trans, border=myOrange, add=TRUE, at=(g+0.1))
#   } else{
#     boxplot(0, outline=FALSE, col=myOrange_trans, border=myOrange, add=TRUE, at=(g+0.1))
#   }
# }





