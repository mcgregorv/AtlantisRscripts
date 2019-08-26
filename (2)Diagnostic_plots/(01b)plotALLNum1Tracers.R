#plot all tracers for a given box and layer
this_run<-"base"

burnin<-2 #number of years to skip in plot
burnin<-0

  this_out<-"Test1"
  # 
  this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
  outPath<-paste(this_path,"output",this_out,"\\",sep="")
  
  plotPath<-paste(this_path,"..\\Figures\\",this_run,"\\",this_out,"",sep="")

  daysTimeStep<-1
  numStepsPerYear<-365/daysTimeStep
  year0<-1900
  fishingStartYear<-1900
  modelStartYear<-1900
  
  ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
  thisVol<-ncvar_get(ThisNC.nc,"volume")
  thisDz<-ncvar_get(ThisNC.nc,"dz")
  ntsteps<-dim(thisVol)[3]
  
  nts<-dim(thisVol)[3]-burnin #number of timesteps
  
  xLabsTemp<-seq(0,(nts*daysTimeStep),by=365)/365
  xLabsAt<-xLabsTemp*numStepsPerYear
  xLabs<-xLabsTemp+year0
  
  #get all tracer names
  allTracers<-sort(names(ThisNC.nc$var))
  num1Tracers<-allTracers[grep("1_Nums",allTracers)]
  
  cohortCols<-colorRampPalette(colors=c(myYellow,myGreen,myAqua,myBlue,myPurple,myRed,"red"))(10)
  
  thisBox<-"all"
  thisLayer<-"all"
  #get layer index. 0 is closest to sediment, it is deepest. up through water column to layer 4
  skip<-c("nominal_dz")
  plotsFile<-paste(plotPath,"ALL_Num1_TRACERS_box",thisBox,"_l",thisLayer,"_compressed.pdf",sep="")
  pdf(plotsFile)
  par(mfrow=c(4,1),mar=c(3,4,2,0),oma=c(1,0,0,0))
  for(t in 1:(length(num1Tracers))){
    thisTracer<-num1Tracers[t]
    temp<-ncvar_get(ThisNC.nc,thisTracer)
      if(length(dim(temp))==3){
        xx<-apply(temp*thisVol,3,sum)
      } else{
        xx<-apply(temp*thisVol[6,,],2,sum)
      }
      #skip any burn in years
      
      if(length(dim(xx))==0){
        if(xx[1]==0){
          thisY<-xx
        } else{
          thisY<-xx/xx[1]
        }
        thisymax<-max(thisY)*1.1
        thisymin<-min(0,min(thisY)*1.1)
        plot(x=seq(1,nts),y=thisY,type="l",col=myGreen,lwd=2.5,ylim=c(thisymin,thisymax),ylab="",xlab="Day",xaxt="n")
        mtext(thisTracer,side=3,adj=0,font=2)
        abline(h=1,col="red",lty=2,lwd=1.5)
        axis(at=xLabsAt,labels=xLabs,side=1)
      } 

    
  }
  dev.off()
