#plot fluxes
this_run<-"base"

this_out<-""

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")

plotPath<-paste(this_path,"..\\Figures\\",this_run,"\\",this_out,"",sep="")

daysTimeStep<-1
numStepsPerYear<-365/daysTimeStep
year0<-1880
fishingStartYear<-1900
modelStartYear<-1880

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
ntsteps<-dim(thisVol)[3]
# xLabs<-seq(1,(ntsteps*daysTimeStep),by=daysTimeStep)

nts<-dim(thisVol)[3] #number of timesteps

xLabsTemp<-seq(0,(nts*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0

nrows=4
naxisLabs<-10

theseTracers<-c("eflux","vflux","Temp","salt","Stress","volume")
for(thisTracer in theseTracers){
  # thisTracer<-"eflux"
  # thisTracer<-"vflux"
  thisData<-ncvar_get(ThisNC.nc,thisTracer)
  axisIndex<-unique(round(pretty(seq(1,dim(thisData)[3]),n=naxisLabs)))+1
  thisAxis<-(axisIndex)/numStepsPerYear+year0-1
  
  nl<-dim(thisData)[1]
  layerCols<-colorRampPalette(colors=c(myRed_trans,myPurple_trans,myBlue_trans,myGreen_trans))(nl)
  
  
  pdf(paste(thisPath,"Figures\\",thisTracer,".pdf",sep=""))
  par(mfrow=c(nrows,1),mar=c(0,4,1,1),oma=c(10,1,1,1),las=1)
  for(b in 1:dim(thisData)[2]){
    thisMax<-max(thisData[,b,])*1.1
    thisMin<-min(thisData[,b,])*0.9
    plot(thisData[1,b,],type="n",xaxt="n",xlab="",ylab="",ylim=c(thisMin,thisMax))
    for(l in 1:nl){
      xx<-thisData[l,b,]
      points(xx,lty=l,lwd=2,col=layerCols[l],type="l")
    }
    mtext(paste("Box ",b,sep=""),side=3,adj=0.01,line=-1.01,font=2,cex=0.8)
    if(round(b/nrows)==b/nrows | b==dim(thisData)[2]){
      # axis(at=axisIndex,labels=thisAxis,side=1,outer=TRUE)
      axis(at=xLabsAt,labels=xLabs,side=1,outer=TRUE)
      par(xpd=NA)
      legend(legend=seq(1,nl),lwd=2,lty=seq(1,nl),col=layerCols,x="bottom",inset = -0.8,horiz=TRUE,title="Layer")
      mtext("Year",side=1,line=3)
      par(las=0)
      mtext(thisTracer,side=2,outer=TRUE,line=-1)
      par(xpd=TRUE,las=1)
    }
  }
  dev.off()
}
