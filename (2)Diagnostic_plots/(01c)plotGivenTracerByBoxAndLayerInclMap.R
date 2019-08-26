library(ggmap)
library(rgdal)
library(gtable)
library(maps)
library(mapdata)
library(rgeos)
source(paste(DIR$'General functions',"\\get_first_number.R",sep=""))
source(paste(DIR$'General functions',"formatShape.R",sep=""))

this_run<-"base"
this_out<-paste("TestPL2",sep="")
# 
basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
outPath<-paste(basePath,this_run,"\\","output",this_out,"\\",sep="")

plotPath<-paste(basePath,"\\Figures\\",this_run,"Zoo1\\TestPL2",sep="")

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1900
fishingStartYear<-1900
modelStartYear<-1900

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
allTracers<-sort(names(ThisNC.nc$var))

nts<-dim(thisVol)[3] #number of timesteps
nbox<-dim(thisVol)[2]; nlayer<-dim(thisVol)[1]

# altNC.nc<-nc_open(paste(outPath,"..\\outputShort\\output.nc",sep=""))
# altDF<-ncvar_get(altNC.nc,"DinoFlag_N")

# df<-ncvar_get(ThisNC.nc,"DinoFlag_N")
# par(mfrow=c(2,2),mar=c(3,4,0,0))
# plot(df[5,2,],type="l"); mtext(2,side=3,line=-1); points(altDF[5,2,],type="l",col=myRed,lty=2,lwd=2)
# plot(df[5,9,],type="l"); mtext(9,side=3,line=-1); points(altDF[5,9,],type="l",col=myRed,lty=2,lwd=2)
# plot(df[5,20,],type="l"); mtext(20,side=3,line=-1); points(altDF[5,20,],type="l",col=myRed,lty=2,lwd=2)
# plot(df[5,12,],type="l"); mtext(12,side=3,line=-1); points(altDF[5,12,],type="l",col=myRed,lty=2,lwd=2)

# df<-ncvar_get(ThisNC.nc,"MicroPB_N")
# plot(df[6,2,],type="l"); mtext(2,side=3,line=-1); plot(df[6,9,],type="l"); mtext(9,side=3,line=-1); plot(df[6,20,],type="l"); mtext(20,side=3,line=-1); plot(df[6,12,],type="l"); mtext(12,side=3,line=-1); 

# shortO2<-ncvar_get(ThisNC.nc,"Oxygen")
# testO2<-ncvar_get(ThisNC.nc,"Oxygen")
# testDepth<-ncvar_get(ThisNC.nc,"dz")
# 
# shortO2vec<-as.double(as.vector(shortO2)); testO2vec<-as.double(as.vector(testO2)); testDepthVec<-as.double(as.vector(testDepth))
# plot(shortO2vec[1:length(testO2vec)],testO2vec)
# index<-testDepth>1
# plot(testDepth[index],testO2vec[index])
# 
# #read in shape file
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
# wcColors<-colorRampPalette(colors=c(myAqua,myBlue,"midnightblue"))((nlayer-1))
# waterColumnColors<-data.frame(cbind("layer"=seq(0,(nlayer-2)),"color"=wcColors))
# sedCol<-myGreen
# #################
## Tracer to plot ##
###################
# thisTracer<-"NH3"
# thisTracer<-"NO3"
# thisTracer<-"Nitrification"
# thisTracer<-"Sed_Bact_N"
# thisTracer<-"Oxygen"
# thisTracer<-"Pelag_Bact_N"
# thisTracer<-"dz"
# thisTracer<-"Pelagic_fish_lge1_Nums"
# thisTracer<-"Arrow_squid1_Nums"
# thisTracer<-"Arrow_squid2_Nums"
# thisTracer<-"DinoFlag_N"
# allPlotTracers<-sort(unique(c("Meiobenth_N","Pelagic_fish_lge1_Nums","Nitrification","NH3","NO3","Oxygen","Pelag_Bact_N","Sed_Bact_N","dz","Lab_Det_N","Ref_Det_N","sedoxdepth","Si","vflux","volume","Macroalgae_N", "Macroalgae_Cover", "MicroPB_N", "MicroPB_Cover", "MicroPB_S","Diatom_N","Diatom_S","PicoPhytopl_N","DinoFlag_N","Chl_a")))
# allPlotTracers<-sort(unique(c("Light","Si","Nitrification","NH3","NO3","Oxygen","Pelag_Bact_N","Sed_Bact_N","Lab_Det_N","Ref_Det_N","PicoPhytopl_N","DinoFlag_N","Chl_a","Macroalgae_N","Diatom_N","Diatom_S","MicroPB_N","Zoo_N")))
# allPlotTracers<-c("PicoPhytopl_N", "Diatom_N")
# 
# allPlotTracers<-c("Chl_a","MicroPB_N","Meiobenth_N","Lab_Det_N")
# 
# allPlotTracers<-paste("Baleen_whales",seq(1,10),"_Nums",sep="")
# 
# allPlotTracers<-paste("Invert_comm_Scav",seq(1,10),"_Nums",sep="")
# 
allPlotTracers<-allTracers[grep("zoo",allTracers,ignore.case = TRUE)]
# 
# allPlotTracers<-paste("Pelagic_fish_sml",seq(1,4),"_Nums",sep="")
# 
# allPlotTracers<-"Chl_a"
# 
# allPlotTracers<-c("Mesopel_fish_Invert_N", "Smooth_oreo_N")
# 
# allPlotTracers<-c("DinoFlag_N","Chl_a", "Ref_Det_N", "Lab_Det_N",allTracers[grep("zoo",allTracers,ignore.case = TRUE)])
# 
# allPlotTracers<-c(allTracers[grep("zoo",allTracers,ignore.case = TRUE)])
# 
# allPlotTracers<-c("DinoFlag_N", "MicroZoo_N", "Diatom_N", "Pelag_Bact_N", "MicroPB_N", "Lab_Det_N")
# 
# allPlotTracers<-c("Meiobenth_N")
# 
# allPlotTracers<-c("Light", "Light_Adaptn_DF")
# allPlotTracers<-c("NH3", "NO3")
# allPlotTracers<-c("Sed_Bact_N")


allPlotTracers<-c("DinoFlag_N", "Filter_Other_N", "Invert_comm_Herb_N" )
allPlotTracers<-c("Det_Si", "Diatom_N", "Si", "PicoPhytopl_N")


for(thisTracer in allPlotTracers){
  
  # thisTracer<-"Carniv_Zoo_N"
  
  thisData<-ncvar_get(ThisNC.nc,thisTracer)

  pdf(paste(plotPath,"Tracer3DwithMap_",thisTracer,this_out,".pdf",sep=""))
  par(mfrow=c(2,2))
  for(b in 1:nbox){
    if(length(dim(thisData))==2){
      temp<-matrix(0,ncol=nts,nrow=nlayer)
      temp[nlayer,]<-thisData[b,]
    }else{
      temp<-thisData[,b,]
    }
    #get depth for this box
    thisDepth<-round(sum(thisDz[,b,1]-1))
    thisMax<-max(temp)
    if(thisMax>0){
      plot(x=0*as.double(temp[1,]),type="n",ylim=c(0,thisMax),ylab="Tracer abundance",xlab="Timestep")
      layerIndex<-rowSums(temp)>0;  layerIndex[nlayer]<-TRUE; layerIndex[is.na(layerIndex)]<-FALSE
      temp<-temp[layerIndex,]
      #last layer is sediment. Other layers are wter column and need to be reversed
      thisWCLayers<-rev(seq((nlayer-2),0)[layerIndex[1:(nlayer-1)]])
      thisNL<-length(thisWCLayers); thisWCcolors<-(waterColumnColors$color[match(thisWCLayers,waterColumnColors$layer)])
      #get last value
      if(length(dim(temp))==0){
        fvalue<-signif(temp[nts],2)
        points(temp,col=sedCol,lwd=2,type="l",lty=2)
        
      }else{
        fvalue<-signif(mean(temp[,nts]),2)
        for(l in 1:thisNL){
          thisCol<-as.character(thisWCcolors[l])
          if(is.na(thisCol)){thisCol<-myGrey}
          points(temp[l,],col=thisCol,lwd=2,type="l")
        }
        points(temp[(thisNL+1),],col=sedCol,lwd=2,type="l",lty=2)
        
      }
      mtext(paste(thisTracer,", box ",b,sep=""),side=3,adj=0.5,outer=FALSE)
      axis(at=fvalue,labels=fvalue,side=4)
      plot(shape)
      # map('nzHires',add=TRUE,col="black",lwd=2)
      # map.axes()
      for(plotB in 1:dim(labeldf)[1]){
        if(b== plotB){thisCol=myGreen}else{thisCol="white"}
        polygon(sdata$shp$shp[[plotB]]$points,col=thisCol)
      }
      mtext(paste("Depth= ",thisDepth," m",sep=""),side=3,adj=1,outer=FALSE)
    }
  }
  dev.off()
}



