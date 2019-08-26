#plot all tracers for a given box and layer
this_run<-"base"

plotPath<-paste(this_path,"..\\Figures\\BALnumbersA_",sep="")

this_out<-c("BAL", "BAL2", "BAL3","BAL4")
# this_out<-c("ShortLonger","ShortLongerCompare"); 
this_out<-c("XXX_mLB1", "XXX_mLB2", "XXX_mLB3","XXX_mLB4","XXX_mLB5","XXX_mLA5")
this_out<-c("XXX_mLA1", "XXX_mLA2", "XXX_mLA3","XXX_mLA4","XXX_mLA5","SENselALL1")


nruns<-length(this_out)
burnin<-1 #number of years to skip in plot

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1900
fishingStartYear<-1900
modelStartYear<-1900

nc_list<-NULL; nts_list<-NULL; min_nts<-1e+12
for(r in 1:nruns){
  outPath<-paste(this_path,"output",this_out[r],"\\",sep="")
  
  nc_list[[r]]<-nc_open(paste(outPath,"output.nc",sep=""))
  thisVol<-ncvar_get(nc_list[[r]],"volume")
  thisDz<-ncvar_get(nc_list[[r]],"dz")
  nts_list[[r]]<-dim(thisVol)[3]-burnin #number of timesteps
  if(nts_list[[r]]<min_nts){min_nts<-nts_list[[r]]}
}

xLabsTemp<-seq(0,(min_nts*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0

#get all tracer names
allTracers<-sort(names(nc_list[[r]]$var))
tracers2plot<-allTracers[grep("Nums",allTracers,invert = FALSE)]; ntracers<-length(tracers2plot)

runCols<-colorRampPalette(colors=c(myGreen,myAqua,myBlue,"midnightblue"))(nruns)

plotsFile<-paste(plotPath,"ALL_Nums_mL.pdf",sep="")
pdf(plotsFile)
par(mfrow=c(4,1),mar=c(3,4,2,0),oma=c(1,0,0,0))
for(t in 1:ntracers){
  thisTracer<-tracers2plot[t]
  temp<-ncvar_get(nc_list[[2]],thisTracer)
  thisVol<-ncvar_get(nc_list[[2]],"volume")
  if(length(dim(temp))==3){
    xx<-apply(temp*thisVol,3,sum)
    thisymax<-max(xx)*1.1
    thisymin<-min(0,min(xx)*1.1)
    plot(xx,type="l",col=runCols[1],lwd=2.5,ylim=c(thisymin,thisymax),ylab="",xlab="Day",xaxt="n")
    mtext(thisTracer,side=3,adj=0,font=2)
    abline(h=1,col="red",lty=2,lwd=1.5)
    axis(at=xLabsAt,labels=xLabs,side=1)
    
    for(r in 1:nruns){
      temp<-ncvar_get(nc_list[[r]],thisTracer)
      thisVol<-ncvar_get(nc_list[[r]],"volume")
      xx<-apply(temp*thisVol,3,sum)
      points(xx,type="l",col=runCols[r],lwd=2.5,lty=r)
      
    }
    legend(legend=this_out,col=runCols,lty=seq(1,nruns),x="bottomright")
  }
}
dev.off()


