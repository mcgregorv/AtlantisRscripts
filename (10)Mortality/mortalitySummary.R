##Read in tracers and summarise how much mortality each group has (estimated from numbers at age) and compare this with M from mL
mg_2_tonne<-2e-8; X_CN<-5.7

#################################### TRACERS
#plot all tracers for a given box and layer
this_run<-"chaos"

this_out<-"DChaosBASE"
burnin<-35 #number of years to skip in plot
# burnin<-1
this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")
groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

# plotPath<-paste(DIR$'Figures',"Validation\\",sep="") ## use this to overwrite the figure used in the report
# plotPath<-paste(DIR$'Figures',"Validation\\Med",sep="") ## doing version using median prop at age rather than final (year 2016)
plotPath<-paste(DIR$'Figures',"\\Validation\\testing\\",this_out,sep="")

doingMedian<-TRUE
# doingMedian<-FALSE


thisBiolFile<-paste(this_path,"..\\CRAM_BH_hybrid_biol.prm",sep="")
biolLines<-readLines(thisBiolFile)


##also read in M from literature to compare
groupBioPars<-read.csv(paste(this_path,"..\\CRAM_B0.csv", sep=""))

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

#grab _N tracers and store as tonnes for each
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
## dump storeNumbersByGroup so can read in elsewhere
# save(list=c("storeNumbersByGroup"), file=paste(this_path,"storeNumbersByGroup",sep=""))

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
t=nts
#now loop through and estimate M and create numbers at age plots
storeJuvMs<-rep(NA,ng); storeAdMs<-rep(NA,ng)
for(g in 1:ng){
  thisName<-str_trim(groupsDF$Name[g], side="both"); thisNumCohorts<-groupsDF$NumCohorts[g]; thisCode<-as.character(groupsDF$Code[g])
  if(thisNumCohorts>2){
    #get props at age CIs
    propLowerCI<-rep(NA,thisNumCohorts); propUpperCI<-rep(NA,thisNumCohorts); propMedian<-rep(NA,thisNumCohorts)
    for(c in 1:thisNumCohorts){
      test<-sort(storeAllPropsAtAge[,g,c])
      propLowerCI[c]<-test[round(0.025*length(test))]; propUpperCI[c]<-test[round(0.975*length(test))]; propMedian[c]<-test[round(0.5*length(test))]
    }
    thisData<-storeNumbersByGroup[,g,]
    totalNumsTracer<-apply(thisData,1,sum, na.rm=TRUE)
    
    thisVar<-paste(thisCode, "_AgeClassSize", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageClassSize<-get_first_number(temp)
    thisVar<-paste(thisCode, "_age_mat", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageMature<-get_first_number(temp)
  
    thisProps<-thisData[t,]/totalNumsTracer[t]
    if(doingMedian==TRUE){thisProps<-propMedian}
    #do juveniles first
    thisNumbers<-thisProps[1:(ageMature+1)]
    juvFittedM<-optim(par=0.1,fn=ll_IC_fn, method="Brent", upper=2, lower=0)$par
    juvFittedNums<-numsICFromM_fn(juvFittedM)
    ## do adults
    thisNumbers<-thisProps[(ageMature+1):thisNumCohorts]
    adFittedM<-optim(par=0.1,fn=ll_IC_fn, method="Brent", upper=2, lower=0)$par
    adFittedNums<-numsICFromM_fn(adFittedM)
    
    
    plotFile<-paste(plotPath,"FittedM",thisCode,".pdf",sep="")
    pdf(plotFile,height=4,width=5)
    par(mar=c(4.5,4.5,1,1), oma=c(0,0,0,0))
    plot(thisProps[!is.na(thisProps)], pch=20,col=myBlue,xlab="Age (years)",ylab="Proportion at age", xaxt="n",cex=thisCex, cex.axis=thisCex, cex.lab=thisCex,ylim=c(0,max(propUpperCI)))
    segments(x0=seq(1,thisNumCohorts), y0=propLowerCI, x1=seq(1,thisNumCohorts), y1=propUpperCI, col=myBlue, lwd=1.5)
    points(juvFittedNums,type="l",lwd=2,col=myOrange)
    points(x=seq((ageMature+1),thisNumCohorts), y=adFittedNums[1:length(adFittedNums)],type="l",lwd=2,col=myGreen)
    axis(at=seq(1,thisNumCohorts), labels=seq(1,thisNumCohorts)*ageClassSize, side=1, cex.axis=thisCex, cex.lab=thisCex)
    mtext(thisCode,side=3,adj=0)
    abline(v=(ageMature+1), col=myGrey_trans, lwd=3)
    dev.off()
    storeAdMs[g]<-adFittedM; storeJuvMs[g]<-juvFittedM
  }
  
}


## compare the fitted to the mL values
mL2M <- function(mL){
  M <- (-365)*log(1-mL)
  return(M)
}
get_mL<-function(x){
  xx<-grep(paste(x, "_mL", sep=""), biolLines)
  this_mLs<-get_first_number(biolLines[(xx+1)], n="all")
  return(this_mLs)
}
storeMLad<-rep(NA, ng); storeMLjuv <- rep(NA, ng); storeMLad_alt <- storeMLad; storeMLjuv_alt <- storeMLjuv
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    this_mLs<-get_mL(as.character(thisCode))
    thisMfrom_mL<-unlist(lapply(this_mLs, mL2M))
    storeMLad[g]<-thisMfrom_mL[2]; storeMLjuv[g]<-thisMfrom_mL[1]
    thisMfrom_mL_alt<-unlist(lapply(this_mLs/2, mL2M))
    storeMLad_alt[g]<-thisMfrom_mL_alt[2]; storeMLjuv_alt[g]<-thisMfrom_mL_alt[1]
  }
}

par(lend=1)
plot(storeMLjuv, type="h", lwd=5, col=myOrange)
points(storeMLad, type="h", lwd=3, col=myBlue)

plot(storeAdMs, type="h", lwd=5, col=myOrange)
points(storeMLad, type="h", lwd=3, col=myBlue)

plot(storeJuvMs, type="h", lwd=5, col=myOrange)
points(storeMLjuv, type="h",  lwd=3, col=myBlue)
points(storeMLjuv_alt, type="h", lwd=1.5, col=myRed)

plot(storeMLjuv/storeJuvMs, pch=20, col=myBlue)
points(storeMLjuv_alt/storeJuvMs, pch=8, col=myRed)
abline(h=1, col="red", lty=2)

plot(storeMLad/storeAdMs, pch=20, col=myOrange)
points(storeMLad_alt/storeAdMs, pch=8, col=myRed)
abline(h=1, col="red", lty=2)

adIndex<-!is.na(storeMLad_alt) & ((storeMLad_alt/storeAdMs)>1)

calcPropsAtAge <-function(M, ageClassSize=1, numCohorts=10){
  x<-rep(NA, numCohorts*ageClassSize)
  x[1]<-1
  for(a in 2:(numCohorts*ageClassSize)){
    x[a]<-x[(a-1)] * exp(-M)
  }
  xx <- x/sum(x)
  xIndex<- sort(rep((1:numCohorts), ageClassSize))
  y <- tapply(xx, xIndex, sum)
  return(y)
}

calcPropsAtSPLITAge <-function(M1, M2, ageClassSize=1, numCohorts=10, ageMature=1){
  x<-rep(NA, numCohorts*ageClassSize)
  x[1]<-1
  for(a in 2:(numCohorts*ageClassSize)){
    if(a <= (ageMature+1)){
      x[a] <- x[(a-1)] * exp(-M1)
    } else{
      x[a]<-x[(a-1)] * exp(-M2)
    }
  }
  xx <- x/sum(x)
  xIndex<- sort(rep((1:numCohorts), ageClassSize))
  y <- tapply(xx, xIndex, sum)
  return(y)
}

for(g in 1:ng){
  # g=grep("LIN", groupsDF$Code)
  thisCode=groupsDF$Code[g]
  Mlit <- groupBioPars$M[groupBioPars$Code==thisCode]
  Mrad <- storeAdMs[g]; Mlad <- storeMLad[g]
  Mrjuv<-storeJuvMs[g]; Mljuv <- storeMLjuv[g]
  thisVar<-paste(thisCode, "_AgeClassSize", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageClassSize<-get_first_number(temp)
  thisVar<-paste(thisCode, "_age_mat", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageMature<-get_first_number(temp)
  
  thisName<-str_trim(groupsDF$Name[g], side="both"); thisNumCohorts<-groupsDF$NumCohorts[g]; thisCode<-as.character(groupsDF$Code[g])
  if(thisNumCohorts>2){
    #get props at age CIs
    propLowerCI<-rep(NA,thisNumCohorts); propUpperCI<-rep(NA,thisNumCohorts); propMedian<-rep(NA,thisNumCohorts)
    for(c in 1:thisNumCohorts){
      test<-sort(storeAllPropsAtAge[,g,c])
      propLowerCI[c]<-test[round(0.025*length(test))]; propUpperCI[c]<-test[round(0.975*length(test))]; propMedian[c]<-test[round(0.5*length(test))]
    }
    plot(propMedian, type="l", ylim=c(0, 0.5), lty=1, col=myGrey)
    for(t in 1:dim(storeNumbersByGroup)[1]){
      thisNumsByAgeC <-storeNumbersByGroup[t,g,]; thisAgeProps<-thisNumsByAgeC/sum(thisNumsByAgeC)
      points(thisAgeProps, type="l", col=myGrey_trans)
    }
 
    
    adLitNums<-calcPropsAtAge(Mlit, ageClassSize=ageClassSize, numCohorts=thisNumCohorts)
    points(adLitNums, type="l", col=myBlue)
    
    MLNums <- calcPropsAtSPLITAge(storeMLjuv[g], storeMLad[g], ageClassSize=ageClassSize, numCohorts = thisNumCohorts, ageMature= ageMature)
    points(MLNums, type="l", col=myGold)

    MLNums_alt <- calcPropsAtSPLITAge(storeMLjuv_alt[g], storeMLad_alt[g], ageClassSize=ageClassSize, numCohorts = thisNumCohorts, ageMature = ageMature)
    points(MLNums_alt, type="l", col=myAqua)
    mtext(thisCode, side=3, adj=0)
  }
}
makeBlankPlot()
legend(legend=c("realised", "M lit", "mL", "mL half"), col=c(myGrey, myBlue, myGold, myAqua), lwd=2, lty=1, x="center")

GroupsAboveMlit <- c(paste(groupsDF$Code[storeMLad > groupBioPars$M],"_ad",  sep=""), paste(groupsDF$Code[storeMLjuv > groupBioPars$M],"_juv",  sep=""))
GroupsAboveMlit<-GroupsAboveMlit[grep("NA_", GroupsAboveMlit, invert = TRUE)]
GroupsZero <- c(paste(groupsDF$Code[storeMLad ==0],"_ad",  sep=""), paste(groupsDF$Code[storeMLjuv ==0],"_juv",  sep=""))
GroupsZero<-GroupsZero[grep("NA_", GroupsZero, invert = TRUE)]
GroupsLessThanHalf <- c(paste(groupsDF$Code[storeMLad < (0.5*groupBioPars$M)],"_ad",  sep=""), paste(groupsDF$Code[storeMLjuv < (0.5*groupBioPars$M)],"_juv",  sep=""))
GroupsLessThanHalf<-GroupsLessThanHalf[grep("NA_", GroupsLessThanHalf, invert = TRUE)]
GroupsLessThanHalf <- GroupsLessThanHalf[grep(paste(GroupsZero,collapse="|"), GroupsLessThanHalf, invert = TRUE)]
#
GroupsGreaterThanHalf <- c(paste(groupsDF$Code[storeMLad >= (0.5*groupBioPars$M)],"_ad",  sep=""), paste(groupsDF$Code[storeMLjuv >= (0.5*groupBioPars$M)],"_juv",  sep=""))
GroupsGreaterThanHalf<-GroupsGreaterThanHalf[grep("NA_", GroupsGreaterThanHalf, invert = TRUE)]
GroupsGreaterThanHalf <- GroupsGreaterThanHalf[grep(paste(GroupsAboveMlit,collapse="|"), GroupsGreaterThanHalf, invert = TRUE)]

nAgeGroups <- length(unique(c(GroupsZero, GroupsAboveMlit, GroupsGreaterThanHalf, GroupsLessThanHalf)))


# nAgeGroups
# [1] 71
# > length(GroupsZero)/nAgeGroups
# [1] 0.2535211
# > length(GroupsAboveMlit)/nAgeGroups
# [1] 0.3943662
# > length(GroupsGreaterThanHalf)/nAgeGroups
# [1] 0.1126761
# > length(GroupsLessThanHalf)/nAgeGroups
# [1] 0.2394366
# > length(GroupsLessThanHalf)/nAgeGroups + length(GroupsZero)/nAgeGroups
# [1] 0.4929577
# > length(GroupsGreaterThanHalf)/nAgeGroups + length(GroupsAboveMlit)/nAgeGroups
# [1] 0.5070423



## of those greater than Mlit, are they greater than M realised..?

calcNumsAtSPLITAge<-function(M1, M2, numsAtAge, ageClassSize=1, numCohorts=10, ageMature=1){
  x<-rep(NA, numCohorts)
  x[1]<-numsAtAge[1]
  for(a in 2:(numCohorts)){
    if(a <= (ageMature+1)){
      x[a] <- numsAtAge[(a-1)] * exp(-M1*ageClassSize)
    } else{
      x[a]<-numsAtAge[(a-1)] * exp(-M2*ageClassSize)
    }
  }
  y <- x
  return(y)
}
## test the mL
par(mfrow=c(2,2), mar=c(4,4,1,1))
g=grep("DPI", groupsDF$Code)
for(g in 1:ng){
  thisCode <- groupsDF$Code[g]
  thisNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisCode]
  if(thisNumCohorts>1){
    Mlit <- groupBioPars$M[groupBioPars$Code==thisCode]
    newNumsByGroup_mL <- 0* storeNumbersByGroup
    newNumsByGroup_mL[1,g,] <- storeNumbersByGroup[1,g,]
    thisVar<-paste(thisCode, "_AgeClassSize", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageClassSize<-get_first_number(temp)
    thisVar<-paste(thisCode, "_age_mat", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageMature<-get_first_number(temp)
    numCohorts <- groupsDF$NumCohorts[g]
    
    plot(newNumsByGroup_mL[1,g,1:thisNumCohorts], type="l", lwd=1.5)
    for(t in 2:nts){
      thisStartNums <- newNumsByGroup_mL[(t-1), g, ]  
      
      newNums <- calcNumsAtSPLITAge(M1=(0.5*storeMLjuv[g]), M2=(0.5*storeMLad[g]), numsAtAge = thisStartNums, ageClassSize = ageClassSize, numCohorts = numCohorts, ageMature= ageMature)
      newNumsByGroup_mL[t,g,(1:length(newNums))] <- newNums
      points(newNums, type="l", col=myGrey_trans)
    }
    points(newNumsByGroup_mL[nts,g,], type="l", lty=2, col=myGreen, lwd=2)
    points(storeNumbersByGroup[nts, g, ], type="l", col="red", lwd=1.5)
    legend(legend=c(paste("mLad",signif(0.5*storeMLad[g],2)), paste("mLjuv", signif(0.5*storeMLjuv[g],2)), paste("Mlit", Mlit), "Mrealised"), col=c(myGreen, myGreen, "black", "red"), lwd=c(2,2,2,2), lty=c(2,2,1,1), x="topright")
    mtext(thisCode, side=3, adj=0)
  }
}
######################################################################

## create tex insert
count<-1; nperpage<-15
skipGroups<-as.character(groupsDF$Code[groupsDF$NumCohorts<3] )## not enough age classes to est M
texFile<-paste(DIR$'Reports',"(01)BaseReport\\EstM_figures.tex", sep="")
cat("", file=texFile, append=FALSE)
thisCaption<-paste("Proportions at age from no-fishing model at 2016 (blue dots) with 95\\% confidence intervals from 1900\\textendash 2016 model outputs (blue lines) with proportions at age from estimated M from juvenile age classes (orange line) and adult age classes (green line).",sep="")
thisLab<-"EstM"
thisFigText<-paste("\\begin{figure}[H]
                   \\centering", sep="")
cat(thisFigText, file=texFile, append=TRUE)
for(g in 1:ncg){
  thisCode<-catchCodes[g]; thisFig<-paste("FittedM", thisCode,".pdf", sep="")
  cat(as.character(thisCode), count, "\n")
  if(!(thisCode %in% skipGroups)){
    if(count==(nperpage+1)){
      #start a new figure
      texFile<-paste(DIR$'Reports',"(01)BaseReport\\EstM_figures2.tex", sep="")
      cat("", file=texFile, append=FALSE)
      thisCaption<-paste("Proportions at age from no-fishing model (blue dots) with proportions at age from estimated M (orange line).",sep="")
      thisLab<-"EstMP2"
      thisFigText<-paste("\\begin{figure}[H]
                         \\centering", sep="")
      cat(thisFigText, file=texFile, append=TRUE)
    }
    thisFigText<-paste("\\includegraphics[width=4.5cm]{", thisFig,"}\n", sep="")
    cat(thisFigText, file=texFile, append=TRUE)
    if(count==nperpage){
      #close off previous figure
      thisFigText<-paste("\\caption{",thisCaption,"}\\label{gif:", thisLab,"}
                         \\end{figure}\n",sep="")
      cat(thisFigText, file=texFile, append=TRUE)
    }
    count<-count+1
  }
    }
thisFigText<-paste("\\caption{",thisCaption,"}\\label{gif:", thisLab,"}
                   \\end{figure}\n",sep="")
cat(thisFigText, file=texFile, append=TRUE)

Mfrom_mL<-function(mL){
  M<-mL
  if(M>0){
   M<-(-1/365)*(log(mL))
  }
  return(M)
}
#get mL for each and convert it into M
storeAdMfrom_mL<-rep(NA, ng); storeJuvMfrom_mL<-rep(NA,ng)
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>2){
    thisVar<-paste(thisCode, "_mL", sep=""); temp<-biolLines[grep(thisVar,biolLines)+1]; mL<-get_first_number(temp, n="all")
    thisM<-unlist(lapply(mL,Mfrom_mL))
    storeAdMfrom_mL[g]<-thisM[2]; storeJuvMfrom_mL[g]<-thisM[1]
  }
}

plotFile<-paste(plotPath,"Msummary.pdf",sep="")
pdf(plotFile,height=5,width=6)
index<-groupsDF$NumCohorts>2
thisX<-seq(1,length(storeJuvMs[index]))
par(lend=1, mar=c(0,4.5,1,1), oma=c(5,0,0,0),mfrow=c(2,1),las=1)
plot(x=thisX,y=storeJuvMs[index], type="h", lwd=5,xaxt="n",ylab="")
abline(v=seq(2,ng,by=2),col=myGrey_trans)
mtext("M",side=2,adj=0.5,line=3.5)
points(x=thisX,y=storeJuvMfrom_mL[index],type="h", lwd=5,col=myOrange)
mtext("Juveniles",side=3,adj=0,line=-1)
plot(x=thisX,y=storeAdMs[index], type="h", lwd=5,xaxt="n",ylab="")
abline(v=seq(2,ng,by=2),col=myGrey_trans)
mtext("M",side=2,adj=0.5,line=3.5)
points(x=thisX,y=storeAdMfrom_mL[index],type="h", lwd=5,col=myGreen)
mtext("Adults",side=3,adj=0,line=-1)
par(las=2)
axis(at=thisX,labels = groupsDF$Code[index], side=1)
dev.off()


plotFile<-paste(plotPath,"Msummary_withLitM.pdf",sep="")
pdf(plotFile,height=5,width=6)
index<-groupsDF$NumCohorts>2
thisX<-seq(1,length(storeJuvMs[index]))
par(lend=1, mar=c(0,4.5,1,1), oma=c(5,0,0,0),mfrow=c(2,1),las=1)
plot(x=thisX,y=storeJuvMs[index], type="h", lwd=5,xaxt="n",ylab="")
abline(v=seq(2,ng,by=2),col=myGrey_trans)
mtext("M",side=2,adj=0.5,line=3.5)
points(x=thisX,y=storeJuvMfrom_mL[index],type="h", lwd=5,col=myOrange)
points(x=thisX, y=groupBioPars$M[index], pch=8,col=myBlue)
points(x=thisX, y=groupBioPars$M[index], pch=8,col=myBlue)
mtext("Juveniles",side=3,adj=0,line=-1)
plot(x=thisX,y=storeAdMs[index], type="h", lwd=5,xaxt="n",ylab="")
abline(v=seq(2,ng,by=2),col=myGrey_trans)
mtext("M",side=2,adj=0.5,line=3.5)
points(x=thisX,y=storeAdMfrom_mL[index],type="h", lwd=5,col=myGreen)
points(x=thisX, y=groupBioPars$M[index], pch=8,col=myBlue)
mtext("Adults",side=3,adj=0,line=-1)
par(las=2)
axis(at=thisX,labels = groupsDF$Code[index], side=1)
dev.off()


thisCode<-"CET"
g=grep(thisCode,groupsDF$Code); 
thisVar<-paste(thisCode, "_AgeClassSize", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageClassSize<-get_first_number(temp)
thisNumsByAgeClass<-storeNumbersByGroup[,g,]
thisNumbers<-thisNumsByAgeClass[1,]
thisMoptm<-optim(par=0.1,fn=ll_IC_fn, method="Brent", upper=2, lower=0)$par
fittedNums<-numsICFromM_fn(signif(thisMoptm,2))
plot(thisNumbers, type="l")
points(fittedNums, type="l", col="green")
signif(thisMoptm,2)


