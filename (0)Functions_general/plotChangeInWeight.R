plotChangeInWeight<-function(changeInWeightData,plotPath){
  nrows<-3
  thisData<-changeInWeightData[[1]]
  xx<-unique(round(pretty(seq(1,length(thisData)),n=10)))+1
  axisIndex<-xx[xx<length(thisData)]
  thisAxis<-names(thisData)[axisIndex]
  pdf(paste(plotPath,"ChangeInWeight.pdf",sep=""))
  par(mfrow=c(nrows,1),mar=c(0,4,0,1),oma=c(5,1,1,1))
  for(i in 1:(length(changeInWeightData))){
    thisName<-names(changeInWeightData)[i]
    thisData<-changeInWeightData[[i]]
    plot(thisData,xaxt="n",xlab="",ylab="",ylim=c(0,5),type="n")
    abline(h=1,col="red",lty=2,lwd=1.5)
    points(thisData,col="cornflowerblue",lwd=2,type="l")
    mtext(thisName,side=3,line=-1.1,adj=0.01,font=2,cex=0.8)
    if(round(i/nrows)==i/nrows){
      axis(at=axisIndex,labels=thisAxis,side=1,outer=TRUE)
    }
  }
  dev.off()
  
  #those that go up a lot
  upLim<-5
  pdf(paste(plotPath,"ChangeInWeightUP.pdf",sep=""))
  par(mfrow=c(nrows,1),mar=c(0,4,0,1),oma=c(5,1,1,1))
  axis_count<-0
  for(i in 1:(length(changeInWeightData))){
    thisName<-names(changeInWeightData)[i]
    thisData<-changeInWeightData[[i]]
    thisMax<-max(thisData)
    if(thisMax>upLim){
      axis_count<-axis_count+1
      plot(thisData,xaxt="n",xlab="",ylab="",ylim=c(0,10),type="n")
      abline(h=1,col="red",lty=2,lwd=1.5)
      points(thisData,col="cornflowerblue",lwd=2,type="l")
      mtext(thisName,side=3,line=-1.1,adj=0.01,font=2,cex=0.8)
      if(round(i/nrows)==i/nrows){
        axis(at=axisIndex,labels=thisAxis,side=1,outer=TRUE)
      }
    }
    
  }
  dev.off()
  
  
  
  #those that go down a lot
  lowLim<-0.1
  pdf(paste(plotPath,"ChangeInWeightDOWN.pdf",sep=""))
  par(mfrow=c(nrows,1),mar=c(0,4,0,1),oma=c(5,1,1,1))
  axis_count<-0
  for(i in 1:(length(changeInWeightData))){
    thisName<-names(changeInWeightData)[i]
    thisData<-changeInWeightData[[i]]
    thisMin<-min(thisData)
    if(thisMin<lowLim){
      axis_count<-axis_count+1
      plot(thisData,xaxt="n",xlab="",ylab="",ylim=c(0,2),type="n")
      abline(h=1,col="red",lty=2,lwd=1.5)
      points(thisData,col="cornflowerblue",lwd=2,type="l")
      mtext(thisName,side=3,line=-1.1,adj=0.01,font=2,cex=0.8)
      if(round(i/nrows)==i/nrows){
        axis(at=axisIndex,labels=thisAxis,side=1,outer=TRUE)
      }
    }
    
  }
  dev.off()
  
}
