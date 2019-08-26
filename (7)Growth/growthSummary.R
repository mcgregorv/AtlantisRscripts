##Read in tracers and summarise size at age for each group and compare with growth curve
mg_2_tonne<-2e-8; X_CN<-5.7
mg_2_grams<-mg_2_tonne*1e+6

#################################### TRACERS
#plot all tracers for a given box and layer
this_run<-"base"

this_out<-"Base"
burnin<-35 #number of years to skip in plot
# burnin<-1

this_path<-paste(DIR$'Base',"ATLANTISmodels\\ArchivedModels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")
groupsDF<-read.csv(paste(this_path,"..\\..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

# plotPath<-paste(DIR$'Figures',"Validation\\",sep="") ## use this to overwrite the figure used in the report
plotPath<-paste(DIR$'Figures',"\\Validation\\testing\\",this_out,sep="")

doingMedian<-TRUE

##also read in M from literature to compare
groupGrowthPars<-read.csv(paste(this_path,"..\\..\\inputs\\supporting\\length2weights.csv", sep=""))

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
storeWeightsByGroup<-array(NA, dim=c(nts, ng, 10))
for(g in 1:ng){
  thisName<-str_trim(groupsDF$Name[g], side="both"); thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    for(c in 1:thisNumCohorts){
      thisTracer<-paste(thisName,c,"_ResN",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer);   thisRN<-apply(thisTemp,3,nonZeroMean)
      thisTracer<-paste(thisName,c,"_StructN",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer); thisSN<-apply(thisTemp,3,nonZeroMean)
      thisWeights<-(thisRN+thisSN)*mg_2_grams*X_CN
      storeWeightsByGroup[,g,c]<-thisWeights[(burnin+1):length(thisWeights)]
    }
  }
}
##now loop through and calculate length
lengthFromWeight<-function(w,a,b){
  l<-(w/a)^(1/b)
  return(l)
}
VBgrowth_fn<-function(t, Linf, k, t0){
  l<-Linf * (1- exp(-k * (t-t0)))
  return(l)
}
# groupGrowthPars
# Code         a      b   Linf      K      t0
# 1   ASQ 2.900e-02 3.0000  35.00 2.4000   0.000
pdf(paste(plotPath, "GrowthCurves.pdf", sep=""), height=10)
par(mar=c(4,4,1,1), oma=c(1,1,1,1), mfrow=c(5,3))
for(g in 1:ng){
  thisCode<-as.character(groupsDF$Code[g]); thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>2 & thisCode %in% groupGrowthPars$Code){
    index<-groupGrowthPars$Code==thisCode
    this_a<-groupGrowthPars$a[index]; this_b<-groupGrowthPars$b[index]
    this_Linf<-groupGrowthPars$Linf[index]; this_K<-groupGrowthPars$K[index]; this_t0<-groupGrowthPars$t0[index]
    if(!is.na(sum(c(this_a, this_b, this_Linf, this_K, this_t0)))){
      #need years per age class
      thisVar<-paste(thisCode, "_AgeClassSize", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageClassSize<-get_first_number(temp)
      thisAges<-seq(1,thisNumCohorts)*ageClassSize; 
      thisAges<-seq(1,thisNumCohorts)*ageClassSize - ageClassSize/2
      thisAgesFilledIn<-seq(0,(thisNumCohorts*ageClassSize))
      ## get weight at age, and convert to length at age from tracers
      thisWeights<-storeWeightsByGroup[,g,]
      thisLengths<-0*thisWeights
      for(c in 1:thisNumCohorts){
        cweights<-thisWeights[,c]
        thisLengths[,c]<-unlist(lapply(cweights, lengthFromWeight, a=this_a, b=this_b))
      }
      ## fit VB growth curve
      intendedLengths<-unlist(lapply(thisAges, VBgrowth_fn, Linf=this_Linf, k=this_K, t0=this_t0))
      intendedLengthsFilledIn<-unlist(lapply(thisAgesFilledIn, VBgrowth_fn, Linf=this_Linf, k=this_K, t0=this_t0))
      plot(x=thisAges, y=intendedLengths, type="n",ylim=c(0,max(max(thisLengths, na.rm=TRUE), max(intendedLengths))), xlab="Age (years)", ylab="Length (cm)")
      for(c in 1:thisNumCohorts){
        points(x=rep(thisAges[c], dim(thisLengths)[1]), y=thisLengths[,c], col=myBlue_trans, pch=20)
      }
      points(x=thisAgesFilledIn, y=intendedLengthsFilledIn, type="l", col="red", lty=2, lwd=2)
      mtext(thisCode,side=3, adj=0)
    }
  }
}
dev.off()

g=grep("HAK", groupsDF$Code)
       thisCode<-as.character(groupsDF$Code[g]); thisNumCohorts<-groupsDF$NumCohorts[g]
            index<-groupGrowthPars$Code==thisCode
         this_a<-groupGrowthPars$a[index]; this_b<-groupGrowthPars$b[index]
         this_Linf<-groupGrowthPars$Linf[index]; this_K<-groupGrowthPars$K[index]; this_t0<-groupGrowthPars$t0[index]
             #need years per age class
           thisVar<-paste(thisCode, "_AgeClassSize", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageClassSize<-get_first_number(temp)
           thisAges<-seq(1,thisNumCohorts)*ageClassSize; 
           thisAges<-seq(1,thisNumCohorts)*ageClassSize - ageClassSize/2
           thisAgesFilledIn<-seq(0,(thisNumCohorts*ageClassSize))
           ## get weight at age, and convert to length at age from tracers
           thisWeights<-storeWeightsByGroup[,g,]
           thisLengths<-0*thisWeights
           for(c in 1:thisNumCohorts){
             cweights<-thisWeights[,c]
             thisLengths[,c]<-unlist(lapply(cweights, lengthFromWeight, a=this_a, b=this_b))
           }
           ## fit VB growth curve
           intendedLengths<-unlist(lapply(thisAges, VBgrowth_fn, Linf=this_Linf, k=this_K, t0=this_t0))
           intendedLengthsFilledIn<-unlist(lapply(thisAgesFilledIn, VBgrowth_fn, Linf=this_Linf, k=this_K, t0=this_t0))
           plot(x=thisAges, y=intendedLengths, type="n",ylim=c(0,max(max(thisLengths, na.rm=TRUE), max(intendedLengths))), xlab="Age (years)", ylab="Length (cm)")
           for(c in 1:thisNumCohorts){
             points(x=rep(thisAges[c], dim(thisLengths)[1]), y=thisLengths[,c], col=myBlue_trans, pch=20)
           }
           points(x=thisAgesFilledIn, y=intendedLengthsFilledIn, type="l", col="red", lty=2, lwd=2)
           mtext(thisCode,side=3, adj=0)
 
