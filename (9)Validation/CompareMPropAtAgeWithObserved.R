##Read in tracers and summarise how much mortality each group has (estimated from numbers at age) and compare this with M from mL
mg_2_tonne<-2e-8; X_CN<-5.7

#################################### TRACERS
#plot all tracers for a given box and layer
this_run<-"base"

this_out<-"base"
burnin<-35 #number of years to skip in plot
# burnin<-1
this_path<-paste(DIR$'Base',"ATLANTISmodels\\ArchivedModels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")
groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

plotPath<-paste(DIR$'Figures',"Validation\\",sep="") ## use this to overwrite the figure used in the report
plotPath<-paste(DIR$'Figures',"\\Validation\\testing\\",this_out,sep="")

doingMedian<-TRUE
doingMedian<-FALSE

##also read in M from literature to compare
groupBioPars<-read.csv(paste(this_path,"..\\..\\CRAM_B0.csv", sep=""))
groupsDFPaper<-read.csv(paste(this_path,"..\\CRAM_groupsPaper.csv", sep=""))


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
#now loop through and estimate M and create numbers at age plots
storeJuvMs<-rep(NA,ng); storeAdMs<-rep(NA,ng)
toPlot<-c()
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
      lowerByC<-rep(NA, thisNumCohorts); upperByC<-rep(NA, thisNumCohorts)
      for(c in 1:thisNumCohorts){
        thisAmean<-obsPropsAtAge[c]; thisAsd<-0.10*thisAmean;
        thisLowerCI<-thisAmean-1.96*thisAsd; thisUpperCI<-thisAmean + 1.96*thisAsd
        lowerByC[c]<-thisLowerCI; upperByC[c]<-thisUpperCI
      }
      
      thisLongName<-gsub("_"," ", groupsDFPaper$Name[g])
      thisLongName<-gsub("-","\n",thisLongName)
      thisLongName<-gsub("\\(","\n",thisLongName); thisLongName<-gsub(")", "", thisLongName)
      
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

      thisYmax<-max(c(upperByC, max(thisProps, na.rm=TRUE)))*1.2
  
      plotFile<-paste(plotPath,"FittedM",thisCode,".pdf",sep="")
      pdf(plotFile,height=4.5,width=5)
      toPlot<-c(toPlot, thisCode)
      par(mar=c(4.5,4.5,3.5,0.5))
      bpdata<-melt(propsByCohortTime)
      boxplot(value ~ X2, data=bpdata, xlab="Age (years)", ylab="Proportion at age", xaxt="n",cex=thisCex, cex.axis=thisCex, cex.lab=thisCex,ylim=c(0,thisYmax), outline=FALSE)
      yy<-c(lowerByC, rev(upperByC)); xx<-c(seq(1,thisNumCohorts), rev(seq(1, thisNumCohorts)))
      polygon(x=xx, y=yy, col=myOrange_trans, border=NA)
      # points(x=seq(1,thisNumCohorts), y=obsPropsAtAge[1:thisNumCohorts], type="l", lwd=2, col="red", lty=1)
      axis(at=seq(1,thisNumCohorts), labels=seq(1,thisNumCohorts)*ageClassSize, side=1, cex.axis=thisCex, cex.lab=thisCex)
      mtext(thisLongName,side=3, adj=0, cex=thisCex)
      abline(v=(ageMature+1), col=myGrey_trans, lwd=3)
     
      # plot(thisProps[!is.na(thisProps)], pch=20,col=myBlue,xlab="Age (years)",ylab="Proportion at age", xaxt="n",cex=thisCex, cex.axis=thisCex, cex.lab=thisCex,ylim=c(0,thisYmax))
      # points(juvFittedNums,type="l",lwd=2,col=myOrange)
      # points(x=seq((ageMature+1),thisNumCohorts), y=adFittedNums[1:length(adFittedNums)],type="l",lwd=2,col=myGreen)
      dev.off()
      storeAdMs[g]<-adFittedM; storeJuvMs[g]<-juvFittedM
    }
  }
  
}

## for poster
g<-grep("SPD", groupsDF$Code)
thisName<-str_trim(groupsDF$Name[g], side="both"); thisNumCohorts<-groupsDF$NumCohorts[g]; thisCode<-as.character(groupsDF$Code[g])
  ## define props at age using M ('obs')
  thisM<-groupBioPars$M[groupBioPars$Code==thisCode]
    thisVar<-paste(thisCode, "_AgeClassSize", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageClassSize<-get_first_number(temp)
    thisVar<-paste(thisCode, "_age_mat", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageMature<-get_first_number(temp)
    # 'observed' props at age using M from lit
    obsPropsAtAge<-getProportionsAtAge(thisM, thisNumCohorts, ageClassSize)
    lowerByC<-rep(NA, thisNumCohorts); upperByC<-rep(NA, thisNumCohorts)
    for(c in 1:thisNumCohorts){
      thisAmean<-obsPropsAtAge[c]; thisAsd<-0.10*thisAmean;
      thisLowerCI<-thisAmean-1.96*thisAsd; thisUpperCI<-thisAmean + 1.96*thisAsd
      lowerByC[c]<-thisLowerCI; upperByC[c]<-thisUpperCI
    }
    
    thisLongName<-gsub("_"," ", groupsDFPaper$Name[g])
    thisLongName<-gsub("-","\n",thisLongName)
    thisLongName<-gsub("\\(","\n",thisLongName); thisLongName<-gsub(")", "", thisLongName)
    
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
    
    thisYmax<-max(c(upperByC, max(thisProps, na.rm=TRUE)))*1.2
    
    thisCex<-8
    plotFile<-paste(DIR$'Figures',"FittedMforPoster",thisCode,".png",sep="")
    png(plotFile, width=3500, height=2000, bg="transparent", type="cairo")
    par(mar=c(15,30,10,0.5), las=1)
    bpdata<-melt(propsByCohortTime)
    boxplot(value ~ X2, data=bpdata, xlab="", ylab="", yaxt="n", xaxt="n",cex=thisCex, cex.axis=thisCex, cex.lab=thisCex,ylim=c(0,thisYmax), outline=FALSE)
    yy<-c(lowerByC, rev(upperByC)); xx<-c(seq(1,thisNumCohorts), rev(seq(1, thisNumCohorts)))
    polygon(x=xx, y=yy, col=myOrange_trans, border=NA)
    axis(at=seq(1,thisNumCohorts), labels=seq(1,thisNumCohorts)*ageClassSize,lwd=0, line=2, side=1,cex=thisCex, cex.axis=8, col.axis="white", col.ticks = "white")
    axis(at=seq(0,1, by=0.1), labels=seq(0,1, by=0.1), side=2, col.axis="white", col.ticks = "white", cex.axis=8, cex=thisCex)
    par(las=0)
    mtext(paste(thisLongName,", mortality rate check", sep=""), font=3, col=myOrange,side=3, adj=0, cex=thisCex*1.2)
    mtext("Proportion at age", side=2, line=20, adj=0.5, cex=thisCex, col="white")
    mtext("Age (years)", side=1, line=8, adj=0.5, cex=thisCex, col="white")
    dev.off()

    

###
## create tex insert
orderByLName<-match(sort(groupsDFPaper$Name), groupsDFPaper$Name)
codesByLName<-groupsDF$Code[orderByLName]

count<-1; nperpage<-24
skipGroups<-as.character(groupsDF$Code[groupsDF$NumCohorts<3] )## not enough age classes to est M
texFile<-paste(DIR$'Reports',"(01)BaseReport\\PropAtAge_figures.tex", sep="")
cat("", file=texFile, append=FALSE)
thisCaption<-paste("Proportions at age using M based on literature (Table \\ref{tab:bioPars}) where available (orange shaded shows 95\\% confidence intervals using CV 10\\%) and from CRAM simulated years 1900\\textendash 2015 (boxplots).",sep="")
thisLab<-"EstM"
thisFigText<-paste("\\begin{figure}[H]
                   \\centering", sep="")
cat(thisFigText, file=texFile, append=TRUE)
for(g in 1:ng){
  thisCode<-codesByLName[g]; thisFig<-paste("FittedM", thisCode,".pdf", sep="")
  if(thisCode %in% toPlot){
    if(count==(nperpage+1)){
      #start a new figure
      texFile<-paste(DIR$'Reports',"(01)BaseReport\\PropAtAge_figures2.tex", sep="")
      cat("", file=texFile, append=FALSE)
      # thisCaption<-paste("Proportions at age from no-fishing model 1900\\textendash 2016 with proportions at age from M in the literature (see Table \ref{tab:bioPars} for references and values) with 95\\% confidence intervals based on CV=10\\% (gold).",sep="")
      thisLab<-"EstMP2"
      thisFigText<-paste("\\begin{figure}[H]
                         \\centering", sep="")
      cat(thisFigText, file=texFile, append=TRUE)
    }
    thisFigText<-paste("\\includegraphics[width=3.5cm]{", thisFig,"}\n", sep="")
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


