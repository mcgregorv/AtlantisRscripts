#read in tracers and plot overlap between given prey and predator, either as a snapshot or a timeseries

this_run<-"base"

this_out<-"EatMoreAS4"

mg_2_tonne<-2e-8; X_CN<-5.7

burnin<-0 #number of years to skip in plot

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")
r<-""

plotPath<-paste(this_path,"..\\Figures\\",this_run,"\\",this_out,"",sep="")

groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1900
fishingStartYear<-1900
modelStartYear<-1900

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
ntsteps<-dim(thisVol)[3]

nlayers<-dim(thisVol)[1]; nboxes<-dim(thisVol)[2]

nts<-dim(thisVol)[3]-burnin #number of timesteps
cat(paste(nts,"\n"))
xLabsTemp<-seq(0,(nts*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0

#get all tracer names
allTracers<-sort(names(ThisNC.nc$var))

cohortCols<-colorRampPalette(colors=c(myYellow,myGreen,myAqua,myBlue,myPurple,myRed,"red"))(10)

thisPrey<-"PFS"; thisPred<-"HOK"
for(thisPred in c("HOK", "ORH","BOE", "PFM", "JAV")){
  cat("\n",thisPred)
  thisOverlap<-getCondPredOverlap(thisPrey,thisPred)
  cat(thisOverlap$min,thisOverlap$mean, thisOverlap$max,sep="--")
}

getCondPredOverlap<-function(thisPrey, thisPred){
  
  preyName<-str_trim(groupsDF$Name[groupsDF$Code==thisPrey],side="both")
  predName<-str_trim(groupsDF$Name[groupsDF$Code==thisPred],side="both")
  
  #just do by total biomass
  
  temp<-ncvar_get(ThisNC.nc,paste(preyName,"_N",sep=""))
  preyBiomass<-(temp*thisVol)
  temp2<-ncvar_get(ThisNC.nc,paste(predName,"_N",sep=""))
  predBiomass<-(temp2*thisVol)
  
  preyScore<-rep(0,nts) #this is a count of cells with the prey in them
  condPredScore<-rep(0,nts) #This is the count of cells with prey 
  for(t in 1:nts){
    for(l in 1:nlayers){
      for(b in 1:nboxes){
        if(preyBiomass[l,b,t]>0){
          preyScore[t]<-preyScore[t]+1
          if(predBiomass[l,b,t]>0){
            condPredScore[t]<-condPredScore[t]+1
          }
        }
      }
    }
  }
  
  #this is a proportion of the cells with prey in them that also have the predator in them
  ratio<-condPredScore/preyScore
  thisMean<-signif(mean(ratio,na.rm=TRUE),2)
  thisMin<-signif(min(ratio,na.rm=TRUE),2)
  thisMax<-signif(max(ratio,na.rm=TRUE),2)
  thisOut<-list("mean"=thisMean, "min"=thisMin, "max"=thisMax)
  
  return(thisOut)
}





