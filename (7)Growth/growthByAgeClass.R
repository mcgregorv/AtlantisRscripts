## code written by Vidette McGregor, NIWA, 2018 ##
###############################################

#read in outputs from a model, reads in length to weight conversion pars, then plots length at age (with CIs) with intended curve
this_run<-"base"
source(paste(DIR$'General functions',"nonZeroMean.R", sep="")); source(paste(DIR$'General functions',"makeBlankPlot.R", sep=""))
source(paste(DIR$'General functions',"calc_length.R", sep=""));
nlayers<-6

this_out<-c("BASE") #the name of the folder the Atlantis outputs are in
thisLineCol<-"black"; ptCol<-myBlue_trans ## i dont think i use this, but just change any colours to what you wnt

# plotPath<-paste(DIR$'Figures',"Validation\\",sep="") ## use this to overwrite the figure used in the report
plotPath<-paste(this_path,"..\\Figures\\growth\\",this_out,sep="")

thisCex=2
burnin<-35 #this is just if you have a burnin period for the model. set as number of timesteps

nruns<-length(this_out)

mg_2_tonne<-2e-8; X_CN<-5.7; mg_2_grams<-2e-5
## get_KWRR and SR only work if biolLines is global - use readLines() to bring in biol.prm lines (line 51 below)
get_KWRR<-function(Code){
  thisPar<-paste("KWRR_",Code,sep="")
  thisOut<-get_first_number(biolLines[grep(thisPar,biolLines)])
  return(thisOut)
}
get_KWSR<-function(Code){
  thisPar<-paste("KWSR_",Code,sep="")
  thisOut<-get_first_number(biolLines[grep(thisPar,biolLines)])
  return(thisOut)
}

this_path<-paste(DIR$'Base',"ATLANTISmodels\\ArchivedModels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
this_nts<-dim(thisVol)[3]
daysTimeStep<-365; numStepsPerYear<-365/daysTimeStep
year0<-1865

xLabsTemp<-seq(0,(this_nts*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0

groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
groupsDFPaper<-read.csv(paste(this_path,"..\\CRAM_groupsPaper.csv", sep=""))

biolLines<-readLines(paste(this_path,"..\\CRAM_BH_hybrid_biol.prm",sep=""))

#you might not have these - just create a .csv file with these values where you have them from other research and such
## a and b are the length to weight conversion pars, the other 3 are the 3-parameter VB growth pars
lw_pars<-read.csv(paste(this_path,"..\\..\\inputs\\supporting\\length2weights.csv", sep="")); 
colnames(lw_pars)<-c("Code", "a", "b","Linf", "k", "t0")

# these store from Atlantis outputs
storeWeightsAtAge<-array(NA,dim=c(ng, 10, this_nts)); storeNumbersAtAge<-array(NA,dim=c(ng, 10, this_nts))
storeBiomassByGroup<-array(NA,dim=c(ng, this_nts))

lengthsAtAge<-0*storeWeightsAtAge

for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode], side="both")
  thisNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisCode]
  if(thisNumCohorts>1){
    #then it's age structured and need numbers and individual weights
    thisBiomass<-rep(0,dim(thisVol)[3])
    for(c in 1:thisNumCohorts){
      thisVar<-paste(thisName,c,"_Nums", sep=""); temp<-ncvar_get(ThisNC.nc, thisVar)
      thisNums<-apply(temp, 3, sum)
      storeNumbersAtAge[g,c,]<-thisNums
      #resN
      thisVar<-paste(thisName,c,"_ResN", sep=""); temp<-ncvar_get(ThisNC.nc, thisVar)
      thisResN<-apply(temp, 3, nonZeroMean)
      #structN
      thisVar<-paste(thisName,c,"_StructN", sep=""); temp<-ncvar_get(ThisNC.nc, thisVar)
      thisStructN<-apply(temp, 3, nonZeroMean)
      thisWeights<-thisResN + thisStructN
      storeWeightsAtAge[g, c, ]<-thisWeights * mg_2_tonne * X_CN
      #biomass
      thisBiomass<-thisBiomass+thisWeights * thisNums
    }
    storeBiomassByGroup[g,]<-thisBiomass
  } else{
    thisVar<-paste(thisName, "_N", sep=""); temp<-ncvar_get(ThisNC.nc, thisVar)
    if(length(dim(temp))==3){
      thisBiomass<-apply(temp*thisVol,3,sum)
    } else{
      thisBiomass<-apply(temp*thisVol[nlayers,,],2,sum)
    }
    storeBiomassByGroup[g,]<-thisBiomass
  }
}

#turn storeWeightsAtAge into lengths at age
for(g in 1:ng){
  thisCode<-as.character(groupsDF$Code[g]); thisNumCohorts<-groupsDF$NumCohorts[g]
  this_a<-lw_pars$a[lw_pars$Code==thisCode]; this_b<-lw_pars$b[lw_pars$Code==thisCode]
  if(length(this_a)>0 & length(this_b)>0){
    for(c in 1:thisNumCohorts){
      thisWeights<-storeWeightsAtAge[g,c,]
      
      thisLengths<-unlist(lapply(thisWeights, calc_length, this_a=this_a, this_b=this_b))
      lengthsAtAge[g,c,]<-thisLengths
    }
  }
}

###########################################
## up to here up to here
## do the plots, with lengths rather than weights, and with CIs

## also, check out the ling mum sensitivity runs and the FSP steepness (base fish) runs
lw_pars[,c("estLinf", "estt0", "estk")]<-NA
storePvalues<-rep(NA, ng); storePbycohort<-array(NA, dim=c(ng, 10)); storeCVrequiredByC<-array(NA, dim=c(ng, 10))
thisCex<-0.8
orderByLName<-match(sort(groupsDFPaper$Name), groupsDFPaper$Name)
codesByLName<-groupsDF$Code[orderByLName]
pdf(paste(plotPath,"AppendixFigB_SizeAtAgeALLgroups.pdf", sep=""), width=6, height=9)
par(mar=c(4.5,4,2,0.1), oma=c(0,0,0,0), mfrow=c(6,3))
for(g in 1:ng){
  thisCode<-codesByLName[g];
  gIndex<-(1:ng)[as.character(groupsDF$Code)==thisCode]
  thisName<-str_trim(groupsDF$Name[gIndex],side="both")
  thisLongName<-gsub("_"," ", groupsDFPaper$Name[gIndex])
  thisLongName<-gsub("-","\n",thisLongName)
  thisLongName<-gsub("\\(","\n",thisLongName); thisLongName<-gsub(")", "", thisLongName)
  cat(as.character(thisCode),"--")
  thisNumCohorts<-groupsDF$NumCohorts[gIndex]
  if(thisNumCohorts>1){
    #get ageclass size
    thisVar<-paste(thisCode, "_AgeClassSize", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageClassSize<-get_first_number(temp)
    lw_index<-lw_pars$Code==as.character(thisCode)
    thisLengths<-lengthsAtAge[gIndex,,burnin:this_nts]
    lengthsMax<-max(thisLengths, na.rm=TRUE)
    thisLinf<-lw_pars$Linf[lw_index]; this_t0<-lw_pars$t0[lw_index]; this_k<-lw_pars$k[lw_index]
    
    if(lengthsMax>0){
      length_fn<-function(t,Linf,t0,k){
        l<-Linf*(1-exp((-k)*(t-t0)))
        return(l)
      }
      test_t<-seq(0.5,9.5, by=1)*ageClassSize
      test_l<-unlist(lapply(test_t, length_fn, Linf=thisLinf, t0=this_t0, k=this_k))
      
      if(sum(test_l, na.rm=TRUE)>0){
        slength<-dim(thisLengths)[2]
        lowerByC<-rep(NA, thisNumCohorts); upperByC<-rep(NA, thisNumCohorts)
        for(c in 1:thisNumCohorts){
          thisAmean<-test_l[c]; thisBmean<-mean(thisLengths[c,])
          thisAsd<-0.10*thisAmean; thisBsd<-sqrt(var(thisLengths[c,])); this_n<-slength
          
          thisLowerCI<-thisAmean-1.96*thisAsd; thisUpperCI<-thisAmean + 1.96*thisAsd
          lowerByC[c]<-thisLowerCI; upperByC[c]<-thisUpperCI
        }
        
        toPlot<-c(toPlot, as.character(thisCode))
        
        # pdf(paste(plotPath,"sizeAtAge_", thisCode,".pdf", sep=""), width=5, height=4.5)
        # par(mar=c(4.5,4.5,3.5,0.5))
        thisYmax<-max(c(test_l, thisLengths),  na.rm=TRUE)
        
        bpdata<-melt(thisLengths)
        boxplot(value ~ X1, data=bpdata, xlab="Age (years)", ylab="Length (cm)", xaxt="n", ylim=c(0,thisYmax), cex.axis=thisCex, cex.lab=thisCex)
        axis(at=seq(1,thisNumCohorts), labels=seq(1,thisNumCohorts)*ageClassSize, side=1, cex.axis=thisCex, cex.lab=thisCex)
        yy<-c(lowerByC, rev(upperByC)); xx<-c(seq(1,thisNumCohorts), rev(seq(1, thisNumCohorts)))
        polygon(x=xx, y=yy, col=myOrange_trans, border=NA)
        mtext(thisLongName,side=3, adj=0, cex=thisCex)
      }
    }
  }
}
dev.off()

## for poster
thisCode<-"HOK";
gIndex<-(1:ng)[as.character(groupsDF$Code)==thisCode]
thisName<-str_trim(groupsDF$Name[gIndex],side="both")
thisLongName<-gsub("_"," ", groupsDFPaper$Name[gIndex])
thisLongName<-gsub("-","\n",thisLongName)
thisLongName<-gsub("\\(","\n",thisLongName); thisLongName<-gsub(")", "", thisLongName)
cat(as.character(thisCode),"--")
thisNumCohorts<-groupsDF$NumCohorts[gIndex]
  #get ageclass size
  thisVar<-paste(thisCode, "_AgeClassSize", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageClassSize<-get_first_number(temp)
  lw_index<-lw_pars$Code==as.character(thisCode)
  thisLengths<-lengthsAtAge[gIndex,,burnin:this_nts]
  lengthsMax<-max(thisLengths, na.rm=TRUE)
  thisLinf<-lw_pars$Linf[lw_index]; this_t0<-lw_pars$t0[lw_index]; this_k<-lw_pars$k[lw_index]

    length_fn<-function(t,Linf,t0,k){
      l<-Linf*(1-exp((-k)*(t-t0)))
      return(l)
    }
    test_t<-seq(0.5,9.5, by=1)*ageClassSize
    test_l<-unlist(lapply(test_t, length_fn, Linf=thisLinf, t0=this_t0, k=this_k))
    
      slength<-dim(thisLengths)[2]
      lowerByC<-rep(NA, thisNumCohorts); upperByC<-rep(NA, thisNumCohorts)
      for(c in 1:thisNumCohorts){
        thisAmean<-test_l[c]; thisBmean<-mean(thisLengths[c,])
        thisAsd<-0.10*thisAmean; thisBsd<-sqrt(var(thisLengths[c,])); this_n<-slength
        
        thisLowerCI<-thisAmean-1.96*thisAsd; thisUpperCI<-thisAmean + 1.96*thisAsd
        lowerByC[c]<-thisLowerCI; upperByC[c]<-thisUpperCI
      }
      
      # thisCex<-3
      # png(paste(DIR$Figures,"sizeAtAge_", thisCode,".png", sep=""), width=1000, height=500, bg="transparent")
      # par(mar=c(7,8,3.5,0.5), las=1)
      # thisYmax<-max(c(test_l, thisLengths),  na.rm=TRUE)
      # 
      # bpdata<-melt(thisLengths)
      # boxplot(value ~ X1, data=bpdata, xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(0,thisYmax*1.2), cex.axis=thisCex, cex.lab=thisCex, lwd=1)
      # axis(at=seq(1,thisNumCohorts), labels=seq(1,thisNumCohorts)*ageClassSize, side=1, cex.axis=thisCex, cex.lab=thisCex, col.axis="white", col.ticks = "white")
      # axis(at=seq(20,100, by=20), labels=seq(20,100, by=20), side=2, col.axis="white", col.ticks = "white", cex.axis=thisCex, cex.lab=thisCex)
      # yy<-c(lowerByC, rev(upperByC)); xx<-c(seq(1,thisNumCohorts), rev(seq(1, thisNumCohorts)))
      # polygon(x=xx, y=yy, col=myOrange_trans, border=NA)
      # par(las=0)
      # mtext(paste(thisLongName,", growth check", sep=""), font=3, col="white",side=3, adj=0, cex=thisCex*1.2)
      # mtext("Age (years)", side=2, line=5, adj=0.5, cex=thisCex, col="white")
      # mtext("Length (cm)", side=1, line=4.5, adj=0.5, cex=thisCex, col="white")
      # dev.off()
      
      
      thisCex<-8
      plotFile<-paste(DIR$'Figures',"sizeAtAge_",thisCode,".png",sep="")
      png(plotFile, width=3500, height=2000, bg="transparent", type="cairo")
      par(mar=c(15,30,10,0.5), las=1)
      thisYmax<-max(c(test_l, thisLengths),  na.rm=TRUE)
      
      bpdata<-melt(thisLengths)
      boxplot(value ~ X1, data=bpdata, xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(0,thisYmax*1.2), cex.axis=thisCex, cex.lab=thisCex, lwd=1)
      axis(at=seq(1,thisNumCohorts), labels=seq(1,thisNumCohorts)*ageClassSize, side=1, cex.axis=thisCex, cex.lab=thisCex, col.axis="white", col.ticks = "white")
      axis(at=seq(20,100, by=20), labels=seq(20,100, by=20), side=2, col.axis="white", col.ticks = "white", cex.axis=thisCex, cex.lab=thisCex)
      yy<-c(lowerByC, rev(upperByC)); xx<-c(seq(1,thisNumCohorts), rev(seq(1, thisNumCohorts)))
      polygon(x=xx, y=yy, col=myOrange_trans, border=NA)
      par(las=0)

            mtext(paste(thisLongName,", growth check", sep=""), font=3, col=myOrange,side=3, adj=0, cex=thisCex*1.2)
      mtext("Length (cm)", side=2, line=20, adj=0.5, cex=thisCex, col="white")
      mtext("Age (years)", side=1, line=8, adj=0.5, cex=thisCex, col="white")
      dev.off()

## BEE
# Paired t-test
# data:  length by group
# t = 35.873, df = 1169, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   9.10292 10.15625
# sample estimates:
#   mean of the differences 
# 9.629587 

##plot them individually for paper, then create latex input file too
thisCex <- 2
toPlot<-c()
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[g],side="both")
  thisLongName<-gsub("_"," ", groupsDFPaper$Name[g])
  thisLongName<-gsub("-","\n",thisLongName)
  thisLongName<-gsub("\\(","\n",thisLongName); thisLongName<-gsub(")", "", thisLongName)
  cat(as.character(thisCode),"--")
  thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    #get ageclass size
    thisVar<-paste(thisCode, "_AgeClassSize", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageClassSize<-get_first_number(temp)
    lw_index<-lw_pars$Code==as.character(thisCode)
    thisLengths<-lengthsAtAge[g,,burnin:this_nts]
    lengthsMax<-max(thisLengths, na.rm=TRUE)
    thisLinf<-lw_pars$Linf[lw_index]; this_t0<-lw_pars$t0[lw_index]; this_k<-lw_pars$k[lw_index]
    
    if(lengthsMax>0){
      length_fn<-function(t,Linf,t0,k){
        l<-Linf*(1-exp((-k)*(t-t0)))
        return(l)
      }
      test_t<-seq(0.5,9.5, by=1)*ageClassSize
      test_l<-unlist(lapply(test_t, length_fn, Linf=thisLinf, t0=this_t0, k=this_k))
      
      if(sum(test_l, na.rm=TRUE)>0){
        slength<-dim(thisLengths)[2]
        lowerByC<-rep(NA, thisNumCohorts); upperByC<-rep(NA, thisNumCohorts)
        for(c in 1:thisNumCohorts){
          thisAmean<-test_l[c]; thisBmean<-mean(thisLengths[c,])
          thisAsd<-0.10*thisAmean; thisBsd<-sqrt(var(thisLengths[c,])); this_n<-slength
          
          thisLowerCI<-thisAmean-1.96*thisAsd; thisUpperCI<-thisAmean + 1.96*thisAsd
          lowerByC[c]<-thisLowerCI; upperByC[c]<-thisUpperCI
        }
        
        toPlot<-c(toPlot, as.character(thisCode))
        
        pdf(paste(plotPath,"sizeAtAge_", thisCode,".pdf", sep=""), width=5, height=4.5)
        par(mar=c(4.5,4.5,3.5,0.5))
        thisYmax<-max(c(test_l, thisLengths),  na.rm=TRUE)
        
        bpdata<-melt(thisLengths)
        boxplot(value ~ X1, data=bpdata, xlab="Age (years)", ylab="Length (cm)", xaxt="n", ylim=c(0,thisYmax), cex.axis=thisCex, cex.lab=thisCex)
        axis(at=seq(1,thisNumCohorts), labels=seq(1,thisNumCohorts)*ageClassSize, side=1, cex.axis=thisCex, cex.lab=thisCex)
        # 
        # plot(x=seq(1,thisNumCohorts), y=rep(lengthsMax, thisNumCohorts), type="n", xlab="Age Class", ylab="Length (cm)", xaxt="n", ylim=c(0,thisYmax), cex.axis=thisCex, cex.lab=thisCex)
        # axis(at=seq(1,thisNumCohorts), labels=seq(1,thisNumCohorts), side=1, cex.axis=thisCex)
        # for(c in 1:thisNumCohorts){
        #   points(x=rep(c,dim(thisLengths)[2]), y=thisLengths[c,], col=myGrey, pch=20, cex=2)
        # }
        # points(x=seq(1,thisNumCohorts), y=test_l[1:thisNumCohorts], type="l", lwd=2, col="red", lty=1)
        yy<-c(lowerByC, rev(upperByC)); xx<-c(seq(1,thisNumCohorts), rev(seq(1, thisNumCohorts)))
        polygon(x=xx, y=yy, col=myOrange_trans, border=NA)
        mtext(thisLongName,side=3, adj=0, cex=thisCex)
        dev.off()
      }
    }
  }
}

### create latex input file. make sure in alphabetical order based on long name

## create tex insert
orderByLName<-match(sort(groupsDFPaper$Name), groupsDFPaper$Name)
codesByLName<-groupsDF$Code[orderByLName]
count<-1; nperpage<-24
skipGroups<-c()
texFile<-paste(DIR$'Reports',"(01)BaseReport\\sizeAtAge_figures.tex", sep="")
cat("", file=texFile, append=FALSE)
thisCaption<-paste("Size-at-age using values based on literature (Table \\ref{tab:bioPars}) where available (orange shaded shows 95\\% confidence intervals using CV 10\\%) and from CRAM simulated years 1900\\textendash 2015 (boxplots).",sep="")
thisLab<-"sizeAtAge"
thisFigText<-paste("\\begin{figure}[H]
                   \\centering", sep="")
cat(thisFigText, file=texFile, append=TRUE)
for(g in 1:ng){
  thisCode<-codesByLName[g]; thisFig<-paste("sizeAtAge_", thisCode,".pdf", sep="")
  if(thisCode %in% toPlot){
    thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode],side="both")
    if(!(thisCode %in% skipGroups)){
      if(count==(nperpage+1)){
        #start a new figure
        texFile<-paste(DIR$'Reports',"(01)BaseReport\\sizeAtAge_figuresPart2.tex", sep="")
        cat("", file=texFile, append=FALSE)
        thisCaption<-paste("Size-at-age using values based on literature (Table \\ref{tab:bioPars}) where available (orange line) and from CRAM simulated years 1900\\textendash 2015 (blue dots).",sep="")
        thisLab<-"sizeAtAge2"
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
}
thisFigText<-paste("\\caption{",thisCaption,"}\\label{gif:", thisLab,"}
                   \\end{figure}\n",sep="")
cat(thisFigText, file=texFile, append=TRUE)


## write cvs out
colnames(storeCVrequiredByC)<-seq(1,10)
tableFile<-paste(DIR$'Tables',"sizeAtAgeCVs.csv", sep="")
write.csv(storeCVrequiredByC, tableFile, row.names = as.character(groupsDF$Code))

#write p values
colnames(storePbycohort)<-seq(1,10)
tableFile<-paste(DIR$'Tables',"sizeAtAge_pvalues.csv", sep="")
write.csv(storePbycohort, tableFile, row.names = as.character(groupsDF$Code))


