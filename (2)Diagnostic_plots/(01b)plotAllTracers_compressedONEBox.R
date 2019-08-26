#plot all tracers for a given box and layer
this_run<-"base"
this_run<-"chaos"

this_out<-"CChaosdown05"; runFolder=""

thisDesc <- paste(runFolder, this_out,sep="")

## pick a box
boxIndex<-2

mg_2_tonne<-2e-8; X_CN<-5.7

thisRun<-"output"

burnin<-0 #number of years to skip in plot

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")

# outPath<-paste(this_path,"ZmQpermute\\",this_out,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")

r<-""

#read in B0's 
groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

plotPath<-paste(this_path,"..\\Figures\\",this_run,"\\",this_out,"",sep="")

  ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
  daysTimeStep<-365
  numStepsPerYear<-365/daysTimeStep
  year0<-1865
  fishingStartYear<-1900
  modelStartYear<-1900


thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
ntsteps<-dim(thisVol)[3]
nboxes<-dim(thisVol)[2]

nts<-dim(thisVol)[3]-burnin #number of timesteps
cat(paste(nts,"\n"))
xLabsTemp<-seq(0,(nts*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0

#get all tracer names
allTracers<-sort(names(ThisNC.nc$var))

cohortCols<-colorRampPalette(colors=c(myYellow,myGreen,myAqua,myBlue,myPurple,myRed,"red"))(10)

thisBox<-paste(boxIndex, collapse="_")
thisLayer<-"all"
#get layer index. 0 is closest to sediment, it is deepest. up through water column to layer 4
skip<-c("nominal_dz", "MicroPB_Cover")
plotsFile<-paste(plotPath,"ALL_TRACERS_box",thisDesc,thisBox,"_l",thisLayer,"_compressed",r,".pdf",sep="")
pdf(plotsFile)
par(mfrow=c(4,1),mar=c(3,4,2,0),oma=c(1,0,0,0))
for(t in 1:(length(allTracers))){
  thisTracer<-allTracers[t]
  temp<-ncvar_get(ThisNC.nc,thisTracer)
  if(!thisTracer %in% skip){
    if(length(dim(temp))==3){
      #check if it is part of a group that should go on one plot
      x<-grep("1_",thisTracer)
      if(length(x)==0){
        biomTemp<-temp[,boxIndex,]*thisVol[,boxIndex,]
        xx<-apply(biomTemp,length(dim(biomTemp)),sum)*(mg_2_tonne*X_CN)
        if(burnin>0){
          nx<-length(xx); xx<-xx[(burnin+1):nx]
        }
      } else{
        thisVar1<-unlist(str_split(thisTracer,"1_"))[1]; thisVar2<-unlist(str_split(thisTracer,"1_"))[2]
        yy<-allTracers[grep(thisVar1,allTracers)]
        theseVars<-yy[grep(thisVar2,yy)]
        if(length(theseVars)==10){
          theseVars<-theseVars[order(c(1,10,seq(2,9)))]
        }
        skip<-c(skip,theseVars) #so not plotted twice
        xx<-data.frame(matrix(NA,ncol=length(theseVars),nrow=dim(thisVol)[3]))
        colnames(xx)<-theseVars
        for(i in 1:(length(theseVars))){
          thisTracer<-theseVars[i]
          temp<-ncvar_get(ThisNC.nc,thisTracer)
          biomTemp<-temp[,boxIndex,]
          xx[,i]<-apply(biomTemp,length(dim(biomTemp)),sum)
        }
        if(burnin>0){
          nx<-dim(xx)[1]
          xx<-xx[(burnin+1):nx,]
        }
      }
      
    } else{
      biomTemp<-temp[boxIndex,]*thisVol[6,boxIndex,]
      if(length(dim(biomTemp))>1){
        xx<-apply(biomTemp,length(dim(biomTemp)),sum)
      } else{
        xx<-biomTemp
      }
      if(burnin>0){
        nx<-length(xx); xx<-xx[(burnin+1):nx]
      }
    }
 
    if(length(dim(xx))==0){
      thisY<-xx
      thisymax<-max(thisY)*1.1
      thisymin<-min(0,min(thisY)*1.1)
      #it we're in here, and it is age-structure, grab B0 too
      thisB0<-NA; thisNumCohorts<-1
      test<-grep("_N",thisTracer)
      if(length(test)>0){
        thisName<-gsub("_N","",thisTracer)
        thisCode<-groupsDF$Code[grep(thisName,groupsDF$Name)]; thisNumCohorts<-groupsDF$NumCohorts[grep(thisName,groupsDF$Name)]
        if(thisNumCohorts[1]>1){
          # thisB0<-thisB0df$B0[thisB0df$Code==thisCode]
          thisymax<-max(thisB0,thisymax,na.rm=TRUE)
        }
      }
      plot(x=seq(1,nts),y=thisY,type="l",col=myGreen,lwd=2.5,ylim=c(thisymin,thisymax),ylab="",xlab="Day",xaxt="n")
      mtext(thisTracer,side=3,adj=0,font=2)
      abline(h=1,col="red",lty=2,lwd=1.5)
      abline(h=thisB0,col=myGrey,lwd=2,lty=4)
      axis(at=xLabsAt,labels=xLabs,side=1)
    } else{
      thisLab<-paste(thisVar1,thisVar2,sep=" ")
      thisRel<-xx
      for(j in 1:dim(thisRel)[2]){
        thisRel[,j]<-thisRel[,j]/thisRel[1,j]
      }
      thisMax<-max(thisRel,na.rm=TRUE)*1.1
      plot(x=seq(1,nts),y=thisRel[,1],type="l",col=cohortCols[1],lwd=2.5,ylim=c(0,thisMax),ylab="",xlab="Day",xaxt="n")
      mtext(thisLab,side=3,adj=0,font=2)
      for(j in 2:length(theseVars)){
        points(x=seq(1,nts),y=thisRel[,j],type="l",col=cohortCols[j],lwd=2.5)
      }
      abline(h=1,col=myGrey,lty=2,lwd=2)
      abline(h=1.5,col=myGrey_trans,lty=1,lwd=2)
      abline(h=0.5,col=myGrey_trans,lty=1,lwd=2)
      axis(at=xLabsAt,labels=xLabs,side=1)
      legend(legend=seq(1,length(theseVars)),lwd=2,lty=1,x="topleft",col=cohortCols[1:(length(theseVars))],ncol=5,bty="n")
    }
    
  }
  
}
dev.off()
