#plot selected tracers by box (not layer) comparing multiple runs
library(ggmap)
library(rgdal)
library(gtable)
library(maps)
library(mapdata)
library(rgeos)
source(paste(DIR$'General functions',"\\get_first_number.R",sep=""))
source(paste(DIR$'General functions',"formatShape.R",sep=""))
this_run<-"TBGB_JP2"
this_path = paste(DIR$'Base', "TBGB\\",sep="")
runDesc <- "jPbase"
this_out<-c(paste(this_run,"\\outputST1",sep=""), paste(this_run,"\\outputBase", sep=""), paste(this_run,"\\outputBaseBurnin35", "",sep=""));
plotPath<-paste(this_path,"\\Figures\\Testing\\Spatial\\",runDesc, sep="")

nruns<-length(this_out)
burnin<-rep(0,nruns) #number of years to skip in plot

runCols<-colorRampPalette(colors=c(myBlue,myAqua,myGold,  myOrange, "red"))(nruns)

year0<-1899

daysTimeStep<-rep(5,nruns)
nc_list<-NULL; nts_list<-NULL; min_nts<-1e+12; daysTimesteps<-list()
for(r in 1:nruns){
  outPath<-paste(this_path,"\\",this_out[r],"\\output",sep="")
  
  nc_list[[r]]<-nc_open(paste(outPath,".nc",sep=""))
  thisVol<-ncvar_get(nc_list[[r]],"volume")
  thisDz<-ncvar_get(nc_list[[r]],"dz")
  nts_list[[r]]<-dim(thisVol)[3]-burnin[r] #number of timesteps
  thisNdays<-nts_list[[r]]*daysTimeStep[r]
  if(nts_list[[r]]<min_nts){min_nts<-nts_list[[r]]}
  daysTimesteps[[r]]<-seq(1,thisNdays, by=daysTimeStep[r])
}

max_nts<-max(nts_list)
daysList<-nts_list*daysTimeStep; 
max_ndays<-max(daysList); maxDaysIndex<-(daysList==max_ndays); 
nboxes<-dim(thisVol)[2]
xLabsTemp<-seq(0,(max_ndays),by=365)/365
xLabsAt<-xLabsTemp
xLabs<-xLabsTemp+year0

## get allTracers, then select ones to plot
allTracers<-sort(names(nc_list[[1]]$var))
temp <- allTracers[grep("_N", allTracers)]; Ntracers<-temp[grep("Nums", temp, invert = TRUE)]

allPlotTracers<-c(Ntracers,"Oxygen", "NO3", "NH3")
nPT<-length(allPlotTracers)

#populate array of tracers
tracersArray<-array(NA, dim=c((nruns), nPT,nboxes, max_nts))
for(r in 1:nruns){
  thisVol<-ncvar_get(nc_list[[r]], "volume"); nlayers <- dim(thisVol)[1]
  for(t in 1:nPT){
    thisTracer<-allPlotTracers[t]
    thisData<-ncvar_get(nc_list[[r]], thisTracer)
    if(length(dim(thisData))==2){
      tracersArray[(r),t,,1:nts_list[[r]]]<-thisData * thisVol[nlayers,,] * X_CN * mg_2_tonne
    } else{
      tracersArray[(r),t,,1:nts_list[[r]]]<-apply(thisData * thisVol,c(2,3), sum, na.rm=TRUE) * mg_2_tonne * X_CN
    }
  }
}

for(t in 1:nPT){
  
  # thisTracer<-"Carniv_Zoo_N"
  thisTracer<-allPlotTracers[t]
  thisData<-tracersArray[,t,,]; 
  
  pdf(paste(plotPath,"Tracer2DwithMap_",thisTracer,runDesc,".pdf",sep=""))
  par(mfrow=c(2,2))
  for(b in 1:nboxes){
    # b=3
    boxData<-thisData[,b,]
    thisMax<-max(max(boxData,na.rm=TRUE)*1.2, 0.1)
    if(thisMax>0){
      plot(boxData[1,1:nts_list[[1]]],type="n",ylim=c(0,thisMax),ylab="Biomas (tonnes)",xlab="Timestep", xaxt="n")
      axis(at=xLabsAt, labels = xLabs, side=1)
      for(r in 1:(nruns)){
        points(boxData[r,1:nts_list[[r]]], type="l", col=runCols[r], lwd=2, lty=r)
      }
      mtext(paste(thisTracer,", box ",b-1,sep=""),side=3,adj=0.5,outer=FALSE) #b-1 starts at zero (which is boundary box)
      # plot(shape)
      # # map('nzHires',add=TRUE,col="black",lwd=2)
      # # map.axes()
      # for(plotB in 1:dim(labeldf)[1]){
      #   if(b== plotB){thisCol=myGreen}else{thisCol="white"}
      #   polygon(sdata$shp$shp[[plotB]]$points,col=thisCol)
      # }
    }
  }
  dev.off()
}

pdf(paste(plotPath, "SpatialCompareLEGEND.pdf"), height=5, width=5)
par(mar=c(0,0,0,0), mfrow=c(1,1))
makeBlankPlot()
legend(legend=this_out, col=runCols, lty=1, lwd=2, x="center")
dev.off()

thisTracer<-"NH3"; thisData<-ncvar_get(nc_list[[3]],thisTracer)
thisTracer<-"NO3"; thisDataNO<-ncvar_get(nc_list[[3]], thisTracer)
temp<-ncvar_get(nc_list[[3]], "Temp")
light<-ncvar_get(nc_list[[3]],"Light")

testSurfaceNO<-apply(thisDataNO[5,,], 2, nonZeroMean)

