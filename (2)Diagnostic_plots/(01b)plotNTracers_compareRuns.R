#plot all tracers for a given box and layer
this_run<-"SS3"
# this_run<-"ClimateChange"
this_run<-"TBGB_JP2"
TBGB=FALSE
TBGB=TRUE
this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
if(TBGB==TRUE){
  this_path = paste(DIR$'Base', "TBGB\\",this_run,"\\",sep="")
}

this_out <- c(paste("TestsSCA",c(1:5), sep="")); plotDescrip <-"SCA"

nlayers<-6

if(TBGB==TRUE){
  groupsDF<-read.csv(paste(this_path,"\\TBGB_Groups.csv",sep="")); ng<-dim(groupsDF)[1]
} else{
  groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
}
# thisB0df<-read.csv(paste(this_path,"..\\CRAM_B0.csv",sep=""))

plotPath<-paste(this_path,"..\\Figures\\", plotDescrip,sep="")

showall<-TRUE

nruns<-length(this_out)
burnin<-rep(1,nruns) #number of years to skip in plot
# burnin<-c(36,1)

runCols<-c(  colorRampPalette(colors=c("midnightblue",myBlue,myAqua,myGold,  myOrange, "red"))(nruns))
# runCols<-c(colorRampPalette(colors=c("midnightblue",myBlue,myAqua,myGreen))(nruns-1), "red")
# 
# runCols <- c(colorRampPalette(colors=c(myBlue,myAqua,myGreen))(4), "black", colorRampPalette(colors=c(myYellow, myOrange,"red"))(4))
# # runCols <- c( "black",colorRampPalette(colors=c(myBlue,myAqua,myGreen))(4))

daysTimeStep<-73
numStepsPerYear<-365/daysTimeStep
year0<-1865
fishingStartYear<-1865
modelStartYear<-1865

if(TBGB==TRUE){
  year0<-1899
  fishingStartYear<-1899
  modelStartYear<-1899
  
}

nc_list<-NULL; nts_list<-NULL; min_nts<-1e+12
for(r in 1:nruns){
  outPath<-paste(this_path,"output",this_out[r],"\\",sep="")
  # if(TBGB==TRUE){
  #   thisRun<-thisRuns[r]
  #   nc_list[[r]]<-nc_open(paste(outPath,thisRun,".nc",sep=""))
  #   
  # } else{
    nc_list[[r]]<-nc_open(paste(outPath,"output.nc",sep=""))
    
  # }
  thisVol<-ncvar_get(nc_list[[r]],"volume")
  thisDz<-ncvar_get(nc_list[[r]],"dz")
  nts_list[[r]]<-dim(thisVol)[3]-burnin[r] #number of timesteps
  if(showall==TRUE){nts_list[[r]]<-dim(thisVol)[3]}
  if(nts_list[[r]]<min_nts){min_nts<-nts_list[[r]]}
}


nts_list
max_nts<-max(nts_list, na.rm=TRUE)



timeList<-NULL; timeMin <- 30000; timeMax <- 0
for(r in 1:nruns){
  this_nts<-nts_list[[r]]; this_burnin <- burnin[r]
  thisYear0<-1865 - this_burnin + 1
  thisSeq <- seq(1, (this_nts-this_burnin +1)*daysTimeStep, by=daysTimeStep)/365
  this_time <-thisYear0 + thisSeq
  timeList[[r]]<-this_time
  if(max(this_time) > timeMax){timeMax<-max(this_time)}
  if(min(this_time) < timeMin){timeMin <- min(this_time)}
}

xLabsTemp<-seq(0,(max_nts*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0+burnin[1]

#get all tracer names
allTracers<-sort(names(nc_list[[r]]$var))
temp<-allTracers[grep("_N",allTracers)]; tracers2plot<-temp[grep("Nums",temp,invert = TRUE)]; 
# tracers2plot<-c(tracers2plot,"Oxygen","Temp","Si", "NO3")
ntracers<-length(tracers2plot)

dynBoxes<-2:24
# dynBoxes<-2:3
storeTracers<-array(NA, dim=c(nruns, length(tracers2plot), max(nts_list)+1))


plotsFile<-paste(plotPath,"ALL_N.pdf",sep="")
pdf(plotsFile)
par(mfrow=c(4,1),mar=c(3,4,2,0),oma=c(1,0,0,0))
for(t in 1:ntracers){
  thisTracer<-tracers2plot[t]
  temp<-ncvar_get(nc_list[[1]],thisTracer)
  thisVol<-ncvar_get(nc_list[[1]],"volume")
    if(length(dim(temp))==3){
      yy<-apply(temp[,dynBoxes,]*thisVol[,dynBoxes,],3,sum) * mg_2_tonne * X_CN
    } else{
      yy<-apply(temp[dynBoxes,]*thisVol[nlayers,dynBoxes,],2,sum) * mg_2_tonne * X_CN
    }
      xx<-yy[burnin[1]:length(yy)]
      if(showall==TRUE){
        xx <- yy
      }
      # storeTracers[1, t, burnin[r]:length(yy)]<- xx
      # storeTracers[1, t, ]<- xx
      thisymax<-max(xx)*1.1
      thisymin<-min(0,min(xx)*1.1)
      if(showall==TRUE){
        plot(x=timeList[[1]], y=xx,type="l",col=runCols[1],lwd=2,ylim=c(thisymin,thisymax*1.5),ylab="Biomass (tonnes)",xlab="Day", xlim=c(timeMin, timeMax))
        mtext(thisTracer,side=3,adj=0,font=2)
        for(r in 2:nruns){
          temp<-ncvar_get(nc_list[[r]],thisTracer)
          thisVol<-ncvar_get(nc_list[[r]],"volume")
          if(length(dim(temp))==3){
            yy<-apply(temp[,dynBoxes,]*thisVol[,dynBoxes,],3,sum) * mg_2_tonne * X_CN
          } else{
            yy<-apply(temp[dynBoxes,]*thisVol[nlayers,dynBoxes,],2,sum) * mg_2_tonne * X_CN
          }
          xx<-yy
          points(x=timeList[[r]], y=xx,type="l",col=runCols[r],lwd=1.5,lty=r)
          
          # storeTracers[r, t, burnin[r]:length(yy)]<- xx
          # storeTracers[r, t, ]<- xx
          
          # legend(legend=this_out,col=runCols,lty=seq(1,nruns),x="bottomleft")
        }
      } else{
        plot(xx,type="l",col=runCols[1],lwd=2,ylim=c(thisymin,thisymax*1.5),ylab="Biomass (tonnes)",xlab="Day",xaxt="n")
        mtext(thisTracer,side=3,adj=0,font=2)
        # abline(h=1,col="red",lty=2,lwd=1.5)
        axis(at=xLabsAt,labels=xLabs,side=1)
        
        for(r in 2:nruns){
          temp<-ncvar_get(nc_list[[r]],thisTracer)
          thisVol<-ncvar_get(nc_list[[r]],"volume")
          if(length(dim(temp))==3){
            yy<-apply(temp[,dynBoxes,]*thisVol[,dynBoxes,],3,sum) * mg_2_tonne * X_CN
          } else{
            yy<-apply(temp[dynBoxes,]*thisVol[nlayers,dynBoxes,],2,sum) * mg_2_tonne * X_CN
          }
          xx<-yy[burnin[r]:length(yy)]
          points(xx,type="l",col=runCols[r],lwd=1.5,lty=r)
          
          # storeTracers[r, t, burnin[r]:length(yy)]<- xx
          # storeTracers[r, t, ]<- xx
          
          # legend(legend=this_out,col=runCols,lty=seq(1,nruns),x="bottomleft")
        }
      }
}

dev.off()


# pdf(paste(plotPath,"_LEGEND.pdf", sep=""), height=7, width=5)
makeBlankPlot()
legend(legend=this_out,col=runCols,lty=seq(1,nruns),x="center", seg.len=3, lwd=3)
# dev.off()

# 
# ## do the high keystoneness ones on their own - add MB and BO too, as they are spectacularly variable so far!
# ks_codes <-c("HOK", "SPD", "PFS", "ORH", "BIS", "SB", "PFM", "CET", "HAK", "LIN", "SND", "MJE"); 
# nks<-length(ks_codes)
# for(k in 1:nks){
#   thisCode <- ks_codes[k]
#   thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode])
#   thisTracer<-paste(thisName,"_N", sep="")
#   temp<-ncvar_get(nc_list[[1]],thisTracer)
#   thisVol<-ncvar_get(nc_list[[1]],"volume")
#   if(length(dim(temp))==3){
#     yy<-apply(temp[,dynBoxes,]*thisVol[,dynBoxes,],3,sum) * mg_2_tonne * X_CN
#   } else{
#     yy<-apply(temp[dynBoxes,]*thisVol[nlayers,dynBoxes,],2,sum) * mg_2_tonne * X_CN
#   }
#   xx <- yy
#   thisymax<-max(xx)*1.1
#   thisymin<-min(0,min(xx)*1.1)
#   thisPlotFile<-paste(plotPath, "fullBiomassTracers",thisCode,sep="")
#   jpeg(paste(thisPlotFile,".jpg", sep=""), quality=3000)
#     plot(x=timeList[[1]], y=xx,type="l",col=runCols[1],lwd=2,ylim=c(thisymin,thisymax*1.5),ylab="Biomass (tonnes)",xlab="Day", xlim=c(timeMin, timeMax))
#     mtext(thisTracer,side=3,adj=0,font=2)
#     for(r in 2:nruns){
#       temp<-ncvar_get(nc_list[[r]],thisTracer)
#       thisVol<-ncvar_get(nc_list[[r]],"volume")
#       if(length(dim(temp))==3){
#         yy<-apply(temp[,dynBoxes,]*thisVol[,dynBoxes,],3,sum) * mg_2_tonne * X_CN
#       } else{
#         yy<-apply(temp[dynBoxes,]*thisVol[nlayers,dynBoxes,],2,sum) * mg_2_tonne * X_CN
#       }
#       xx<-yy
#       points(x=timeList[[r]], y=xx,type="l",col=runCols[r],lwd=1.5,lty=r)
#     }
#     dev.off()
# }
# 
# 
