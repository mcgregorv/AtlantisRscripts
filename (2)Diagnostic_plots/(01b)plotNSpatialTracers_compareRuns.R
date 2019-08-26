#plot selected tracers by box (not layer) comparing multiple runs
library(ggmap)
library(rgdal)
library(gtable)
library(maps)
library(mapdata)
library(rgeos)
source(paste(DIR$'General functions',"\\get_first_number.R",sep=""))
source(paste(DIR$'General functions',"formatShape.R",sep=""))

this_run<-"base"
this_run<-"Chaos"
this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")

this_out <- c("ChaosUp05", "ChaosUp02", "ChaosUp01", "ChaosUp005", "ChaosBASE", "ChaosDown005", "ChaosDown01", "ChaosDown02", "ChaosDown05");
plotDescrip<-"ChaosUpDownBASE"


this_out <- paste(c("ChaosUp05", "ChaosUp02", "ChaosUp01", "ChaosBASE", "ChaosDown005", "ChaosDown01", "ChaosDown02", "ChaosDown05"), "short", sep="");
plotDescrip<-"ChaosUpDownBASEshort"

this_out <- c(rep(paste("B",c("Chaosdown005", "Chaosup01", "Chaosup02", "Chaosup005", "Chaosup05"), sep=""),1),"BOChaosup05");
plotDescrip<-"ChaosBOBASE"


this_out <- c(paste("C",c("ChaosUp05", "ChaosUp02", "ChaosUp01", "ChaosUp005","ChaosBASE", "ChaosDown005", "ChaosDown01", "ChaosDown02", "ChaosDown05"), sep=""));
plotDescrip<-"ChaosCBASE"

## D is missing up05 for now..
this_out <- c(paste("D",c("ChaosUp02", "ChaosUp01", "ChaosUp005","ChaosBASE", "ChaosDown005", "ChaosDown01", "ChaosDown02", "ChaosDown05"), sep=""));
plotDescrip<-"ChaosDBASE"



# this_out <- c(paste("BO",c( "ChaosUp02","ChaosBASE","ChaosDown02"), sep=""));
# plotDescrip<-"ChaosBOBASE"


# this_out <- c("CChaosBase", paste("MB0", 4, sep=""),"MBmum", "MBLight","MBLYSIS"); plotDescrip <-"testMB"


plotPath<-paste(this_path,"..\\Figures\\",plotDescrip,sep="")

nruns<-length(this_out)
burnin<-rep(0,nruns) #number of years to skip in plot

runCols<-colorRampPalette(colors=c(myBlue,myAqua,myGold,  myOrange, "red"))(nruns)
runCols <- c(colorRampPalette(colors=c(myBlue,myAqua,myGreen))(4), "black", colorRampPalette(colors=c(myYellow, myOrange,"red"))(4))
# runCols <- c(colorRampPalette(colors=c(myBlue,myAqua,myGreen))(3), "black", colorRampPalette(colors=c(myYellow, myOrange,"red"))(4))

# runCols <- c(colorRampPalette(colors=c(myBlue,myAqua,myGreen))(4), "black", colorRampPalette(colors=c(myYellow, myOrange,"red"))(4))

# # #read in shape file
# pdf("testShape.pdf")
# shapeFile<-paste(DIR$'Base',"ATLANTISmodels\\inputs\\bgm\\CHAT30_LL",sep="")
# sdata<-read.shapefile(shapeFile)
# shape<-formatShape(shapeFile=shapeFile)
# ns<-length(shape)
# SpDF <- SpatialPolygonsDataFrame(shape,data.frame( z=1:ns, row.names=paste("P",seq(1,(ns)),sep="")))
# labels<-seq(1,(ns))
# plot(shape)
# LABELpOS<-polygonsLabel(shape, labels = labels, cex=.1,doPlot=FALSE)
# labeldf<-data.frame(cbind("x"=LABELpOS[1:ns],"y"=LABELpOS[(ns+1):(2*ns)]))
# dev.off()

# this_path<-paste(DIR$'Base',"ATLANTISmodels\\projections\\",sep="")
this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")

year0<-1865
timeList<-NULL; timeMin <- 30000; timeMax <- 0
daysTimeStep<-rep(365,nruns)
nc_list<-NULL; nts_list<-NULL; min_nts<-1e+12; daysTimesteps<-list()
max_nts<-0
for(r in 1:nruns){
  outPath<-paste(this_path,"output",this_out[r],"\\",sep="")
  nc_list[[r]]<-nc_open(paste(outPath,"output.nc",sep=""))
  thisVol<-ncvar_get(nc_list[[r]],"volume")
  thisDz<-ncvar_get(nc_list[[r]],"dz")
  nts_list[[r]]<-dim(thisVol)[3]-burnin[r] #number of timesteps
  thisNdays<-nts_list[[r]]*daysTimeStep[r]
  if(nts_list[[r]]<min_nts){min_nts<-nts_list[[r]]}
  daysTimesteps[[r]]<-seq(1,thisNdays, by=daysTimeStep[r])
  this_time <- year0 : (year0 + nts_list[[r]] -1)
  timeList[[r]]<-this_time
  if(max(this_time) > timeMax){timeMax<-max(this_time)}
  if(min(this_time) < timeMin){timeMin <- min(this_time)}
  if(nts_list[[r]]>max_nts){max_nts<-nts_list[[r]]}
  
}

numStepsPerYear<-365/daysTimeStep[1]
xLabsTemp<-seq(0,(max_nts*daysTimeStep[1]),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0+burnin[1]


allPlotTracers<-c("PicoPhytopl_N","DinoFlag_N", "Oxygen", "Diatom_N", "Pelag_Bact_N", "MicroPB_N", "Lab_Det_N","NH3", "NO3", "Nitrification","Macroalgae_N")
allPlotTracers<-c("Invert_comm_Herb_N", "Nitrification","Rock_lobster_N", "Invert_comm_Scav_N", "Arrow_squid_N", "Filter_Other_N", "Carniv_Zoo_N", "Gelat_Zoo_N",  "MicroZoo_N",   "Zoo_N" , "PicoPhytopl_N", "Deposit_Feeder_N","DinoFlag_N", "MicroPB_N", "dz", "Macroalgae_N","Diatom_N", "Light", "NH3", "NO3", "Meiobenth_N", "Pelag_Bact_N", "Sed_Bact_N", "Lab_Det_N", "Ref_Det_N", "Meiobenth_N", "Benthic_Carniv_N")
# allPlotTracers<-c( "Invert_comm_Herb_N")
# allPlotTracers <-c("Rock_lobster_N", "Invert_comm_Scav_N", "Arrow_squid_N", "Filter_Other_N")
# allPlotTracers <- c("soft", "flat", "reef")
nPT<-length(allPlotTracers)


# populate array of tracers
nboxes<-30
tracersArray<-array(NA, dim=c((nruns), nPT,nboxes, max_nts))
for(r in 1:nruns){
  for(t in 1:nPT){
    thisTracer<-allPlotTracers[t]
    thisData<-ncvar_get(nc_list[[r]], thisTracer)
    if(length(dim(thisData))==2){
      tracersArray[(r),t,,1:nts_list[[r]]]<-thisData
    } else{
      tracersArray[(r),t,,1:nts_list[[r]]]<-apply(thisData,c(2,3), sum, na.rm=TRUE)
    }
  }
}

for(t in 1:nPT){
  
  # thisTracer<-"Carniv_Zoo_N"
  thisTracer<-allPlotTracers[t]
  thisData<-tracersArray[,t,,]; 
    
  pdf(paste(plotPath,"Tracer2DwithMap_",thisTracer,this_out,".pdf",sep=""))
  par(mfrow=c(2,2))
  for(b in 1:nboxes){
  # b=3
    boxData<-thisData[,b,]
    thisMax<-max(boxData,na.rm=TRUE)
    if(thisMax>0){
      plot(x=timeList[[1]],y=boxData[1,1:nts_list[[1]]],type="n",ylim=c(0,thisMax),ylab="Tracer abundance",xlab="Timestep",  xlim=c(timeMin, timeMax))
       for(r in 1:(nruns)){
        points(x=timeList[[r]],y=boxData[r,1:nts_list[[r]]], type="l", col=runCols[r], lwd=2, lty=r)
      }
      mtext(paste(thisTracer,", box ",b,sep=""),side=3,adj=0.5,outer=FALSE)
    plot(shape)
    # map('nzHires',add=TRUE,col="black",lwd=2)
    # map.axes()
    for(plotB in 1:dim(labeldf)[1]){
      if(b== plotB){thisCol=myGreen}else{thisCol="white"}
      polygon(sdata$shp$shp[[plotB]]$points,col=thisCol)
    }
  }
  }
  dev.off()
}

makeBlankPlot()
legend(legend=this_out, col=runCols, lty=1:nruns, lwd=2, x="center")


1560
1650
1680


