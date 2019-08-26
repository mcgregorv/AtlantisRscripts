#plot all tracers for a given box and layer
this_run<-"base"

this_out<-c("TestExtMortOFF","TestExtMortONall2")  #this is the outputs went into - 
# I did two matching runs, one with forced external mort on and one with off to compare

nruns<-length(this_out)
burnin<-1 #number of years to skip in plot

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")

plotPath<-paste(this_path,"..\\figures\\test2_",sep="")

#will need groups.csv file so can lookup group name and number of cohorts - but can enter these manually if you prefer
groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep=""))

#this is just used for the time-axis
daysTimeStep<-1
numStepsPerYear<-365/daysTimeStep
year0<-1900

#reads in the output.nc files and stores them in a list; also stores number of timesteps for each run
nc_list<-NULL; nts_list<-NULL; min_nts<-1e+12
for(r in 1:nruns){
  outPath<-paste(this_path,"output",this_out[r],"\\",sep="")
  
  nc_list[[r]]<-nc_open(paste(outPath,"output.nc",sep=""))
  thisVol<-ncvar_get(nc_list[[r]],"volume")
  thisDz<-ncvar_get(nc_list[[r]],"dz")
  nts_list[[r]]<-dim(thisVol)[3]-burnin #number of timesteps
  if(nts_list[[r]]<min_nts){min_nts<-nts_list[[r]]}
}

# xLabsTemp<-seq(0,(min_nts*daysTimeStep),by=365)/365
# xLabsAt<-xLabsTemp*numStepsPerYear
# xLabs<-xLabsTemp+year0

#just getting numLayers and numBoxes - can replace and do manually if prefer
numLayers<-dim(thisVol)[1]; numBoxes<-dim(thisVol)[2]

#this is the file was created in SetUpExternalForcingMortality.R
#reading it back in so can mark change timesteps on the plot
mortForcingFile<-paste(this_path,"mortForcingTEST.txt",sep="")

#actually, I'm just going to copy these manually from the set up R script..
#define timesteps for mortality changes
timeSteps<-c(5,10,50,150,300); nts<-length(timeSteps) 

#set scalars for mortality. most likely will want first value to be 1 since 
mortScalars<-c(1,5,0.2,1,10)  ##needs to be same length as timesteps - will be applied in same order as timesteps

#replace this with "SCA"
thisGroupCode<-"HOK"
#get name and NumCohorts from groupsDF - replace with manual if prefer
thisGroupName<-str_trim(groupsDF$Name[groupsDF$Code==thisGroupCode],side="both") #str_trim just trims white space
thisNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisGroupCode]

#want to check out numbers for each of the cohorts - this defines the tracer variables for each\
tracers2plot<-paste(thisGroupName,seq(1,thisNumCohorts),"_Nums",sep=""); ntracers<-length(tracers2plot)


runCols<-c("blue","black")

plotsFile<-paste(plotPath,"TesingExternalScaling.pdf",sep="")
pdf(plotsFile)
par(mfrow=c(4,1),mar=c(4,4,2,4),oma=c(1,0,0,0))
for(t in 1:ntracers){
  thisTracer<-tracers2plot[t]
  temp<-ncvar_get(nc_list[[1]],thisTracer)
  thisVol<-ncvar_get(nc_list[[1]],"volume")
  if(length(dim(temp))==3){
    xx<-apply(temp*thisVol,3,sum)
    thisymax<-max(xx)*1.1
    thisymin<-min(0,min(xx)*1.1)
    plot(xx,type="l",col=runCols[1],lwd=2.5,ylim=c(thisymin,thisymax),ylab="",xlab="Day")
    mtext(thisTracer,side=3,adj=0,font=2)

    for(r in 2:nruns){
      temp<-ncvar_get(nc_list[[r]],thisTracer)
      thisVol<-ncvar_get(nc_list[[r]],"volume")
      xx<-apply(temp*thisVol,3,sum)
      points(xx,type="l",col=runCols[r],lwd=2.5,lty=r)
      
    }
    ##add mortality changes
    par(new=TRUE)
    plot(x=xx,ylim=c(0,max(mortScalars)),type="n",xlab="",ylab="",xaxt="n", yaxt="n")
    points(x=timeSteps, y=mortScalars,type="l",col="yellow",lwd=3,lend=1)
    axis(at=mortScalars,labels=mortScalars,side=4)
    mtext("Mortality scalar",side=4,adj=0.5,line=2.5)
    mtext("Numbers",side=2,adj=0.5,line=2.5)
  }
}
dev.off()


