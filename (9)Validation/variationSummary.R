##Read in tracers and summarise how much each groups' biomass varied (goal is within 20% of mean)
source(paste(DIR$'General functions',"getCIfromCV.R", sep=""))
mg_2_tonne<-2e-8; X_CN<-5.7

#################################### TRACERS
#plot all tracers for a given box and layer
this_run<-"base"
this_out<-"BASE"
# this_out<-"TEST150yrbase"
burnin<-35 #number of years to skip in plot

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")
groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
groupsDFPaper<-read.csv(paste(this_path,"..\\CRAM_groupsPaper.csv", sep=""))

thisCex<-1.5

# plotPath<-paste(DIR$'Figures',"Validation\\",sep="") ## use this to overwrite the figure used in the report
plotPath<-paste(DIR$'Figures',"\\Validation\\testing\\",this_out,sep="")

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1900
fishingStartYear<-1900
modelStartYear<-1900

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
ntsteps<-dim(thisVol)[3]; nlayers<-dim(thisVol)[1]

nts<-dim(thisVol)[3]-burnin #number of timesteps
cat(paste(nts,"\n"))
xLabsTemp<-seq(0,(nts*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0

#get all tracer names
allTracers<-sort(names(ThisNC.nc$var))

#grab _N tracers and store as tonnes for each
storeBiomassByGroup<-array(NA, dim=c(nts, ng))
for(g in 1:ng){
  thisName<-str_trim(groupsDF$Name[g], side="both")
  thisTracer<-paste(thisName, "_N", sep="")
  temp<-ncvar_get(ThisNC.nc, thisTracer)
  if(length(dim(temp))==3){
    xx<-apply(temp*thisVol,3,sum)*mg_2_tonne*X_CN
  } else{
    xx<-apply(temp*thisVol[nlayers,,],2,sum)*mg_2_tonne*X_CN
  }
  thisBiomass<-xx[(burnin+1):length(xx)]
  storeBiomassByGroup[,g]<-thisBiomass
}

##plot them with median and mean as ref point, to see which makes most sense
plotFile<-paste(plotPath,"testingRefPoints.pdf", sep="")
pdf(plotFile)
par(mfrow=c(5,2), mar=c(4,4,1,1), oma=c(1,1,1,1))
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]
  thisBiomass<-storeBiomassByGroup[,g]
  plot(thisBiomass, type="l", lwd=2, ylim=c(0,max(thisBiomass)*1.2), xaxt="n")
  axis(at=xLabsAt, labels=xLabs, side=1)
  thisMean<-mean(thisBiomass); thisMed<-median(thisBiomass)
  abline(h=thisMean, col=myOrange, lwd=2, lty=2)
  abline(h=c(0.8,1.2)*thisMean, col=myOrange_trans, lwd=2, lty=1)
  mtext(thisCode, side=3,adj=0)
}
dev.off()

##plot them with median and mean as ref point, to see which makes most sense
orderByLName<-match(sort(groupsDFPaper$Name), groupsDFPaper$Name)
codesByLName<-groupsDF$Code[orderByLName]

storeCVs<-rep(NA,ng)
for(g in 1:ng){
  thisCode<-codesByLName[g]; thisNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisCode]
  if(thisNumCohorts>1){
    thisLongName<-gsub("_"," ", groupsDF$LongName[groupsDF$Code==thisCode])
    thisLongName<-gsub("-","\n",thisLongName)
    thisLongName<-gsub("\\(","\n",thisLongName); thisLongName<-gsub(")", "", thisLongName)
    
    plotFile<-paste(plotPath,"testingRefPoints_CV_",thisCode,".pdf", sep="")
    pdf(plotFile, width=6, height=3.2)
    par(mar=c(3.5,4.5,3.5,0.1))
    thisBiomass<-storeBiomassByGroup[,groupsDF$Code==thisCode]
    thisVariance<-var(thisBiomass); thisMean<-mean(thisBiomass); thisCV<-sqrt(thisVariance)/thisMean
    storeCVs[g]<-thisCV
    thisLimits<-getCIfromCV(indices = thisMean, cv=0.2)
    plot(thisBiomass, type="l", lwd=2, ylim=c(0,max(thisBiomass)*1.2), xaxt="n", ylab="Biomass (tonnes)", xlab="", cex.axis=thisCex, cex.lab=thisCex)
    axis(at=xLabsAt, labels=xLabs, side=1, cex.axis=thisCex)
    polygon(x=c(0, 0, length(thisBiomass),  length(thisBiomass)), y=c(thisLimits$LowerCI, thisLimits$UpperCI, thisLimits$UpperCI, thisLimits$LowerCI), col=myOrange_trans, border=NA)
    mtext(thisLongName, side=3,adj=0, cex=thisCex)
    dev.off()
  }
}
# dev.off()

## create tex insert
orderByLName<-match(sort(groupsDFPaper$Name), groupsDFPaper$Name)
codesByLName<-groupsDF$Code[orderByLName]
count<-1; nperpage<-24
skipGroups<-groupsDF$Code[groupsDF$NumCohorts==1]
texFile<-paste(DIR$'Reports',"(01)BaseReport\\biomassCVs_figures.tex", sep="")
cat("", file=texFile, append=FALSE)
thisCaption<-paste("Simulated biomass from the un-fished model (black line) with 95\\% confidence intervals based on 20\\% CVs (Coefficient of Variation) shaded orange by species group.",sep="")
thisLab<-"biomassVC"
thisFigText<-paste("\\begin{figure}[H]
                   \\centering", sep="")
cat(thisFigText, file=texFile, append=TRUE)
for(g in 1:ng){
  thisCode<-codesByLName[g]; thisFig<-paste("testingRefPoints_CV_", thisCode,".pdf", sep="")
  # if(thisCode %in% toPlot){
    thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode],side="both")
    if(!(thisCode %in% skipGroups)){
      if(count==(nperpage+1)){
        #start a new figure
        texFile<-paste(DIR$'Reports',"(01)BaseReport\\biomassCVs_figuresPart2.tex", sep="")
        cat("", file=texFile, append=FALSE)
        thisCaption<-paste("Simulated biomass from the un-fished model (black line) with 95\\% confidence intervals based on 20\\% CVs (Coefficient of Variation) shaded orange by species group.",sep="")
        thisLab<-"sizeAtAge2"
        thisFigText<-paste("\\begin{figure}[H]
                           \\centering", sep="")
        cat(thisFigText, file=texFile, append=TRUE)
      }
      thisFigText<-paste("\\includegraphics[width=4.7cm]{", thisFig,"}\n", sep="")
      cat(thisFigText, file=texFile, append=TRUE)
      if(count==nperpage){
        #close off previous figure
        thisFigText<-paste("\\caption{",thisCaption,"}\\label{gif:", thisLab,"}
                           \\end{figure}\n",sep="")
        cat(thisFigText, file=texFile, append=TRUE)
      }
      count<-count+1
    }
      # }
}
thisFigText<-paste("\\caption{",thisCaption,"}\\label{gif:", thisLab,"}
                   \\end{figure}\n",sep="")
cat(thisFigText, file=texFile, append=TRUE)



##cv summary
index<-!is.na(storeCVs) #cvs are already in order of long names

plotFile<-paste(plotPath,"testingRefPoints_CVsummary.pdf", sep="")
pdf(plotFile,width=10, height=3)
par(mar=c(4,4,1,1), oma=c(1,1,1,1), lend=1)
maxCV<-max(storeCVs, na.rm=TRUE)
plot(rep(maxCV,ng)[index], type="n", ylab="CV", xlab="",xaxt="n",ylim=c(0,maxCV))
points(storeCVs[index], type="h", col=myBlue, lwd=5)
abline(h=0.2,col=myOrange,lwd=2, lty=2)
par(las=2)
axis(at=seq(1,length(codesByLName[index])), labels=codesByLName[index], side=1)
dev.off()

