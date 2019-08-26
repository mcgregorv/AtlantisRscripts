#plot all tracers for a given box and layer

this_run<-"ClimateChange"
this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")

this_out<-c("TestOutIN50yr", "ForceNFromBurnin","ForceSiFromBurnin", "ForceTemperatureFromBurnin", "ForceClimateFromBurnin") # have the climate run going through now (N, Si, temperature)
# this_out<-c("TestOutIN50yr", "TestNFromBurnin","TestSiFromBurnin", "TestTemperatureFromBurnin") 
plotDescrip<-"ForceChangesFromBurnin"; thisRuns<-""

groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]


plotPath<-paste(this_path,"..\\Figures\\", plotDescrip,sep="")

nruns<-length(this_out)
burnin<-rep(1,nruns) #number of years to skip in plot

runCols<-colorRampPalette(colors=c("midnightblue",myBlue,myAqua,myGold,  myOrange, "red"))(nruns)

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1865
fishingStartYear<-1865
modelStartYear<-1865

nc_list<-NULL; nts_list<-NULL; min_nts<-1e+12
for(r in 1:nruns){
  outPath<-paste(this_path,"output",this_out[r],"\\",sep="")
  nc_list[[r]]<-nc_open(paste(outPath,"output.nc",sep=""))
  thisVol<-ncvar_get(nc_list[[r]],"volume")
  thisDz<-ncvar_get(nc_list[[r]],"dz")
  nts_list[[r]]<-dim(thisVol)[3]-burnin[r] #number of timesteps
  if(nts_list[[r]]<min_nts){min_nts<-nts_list[[r]]}
}
nts_list
max_nts<-max(nts_list, na.rm=TRUE)
daysTimeStep=365; numStepsPerYear=365/daysTimeStep
xLabsTemp<-seq(0,(max_nts*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0+burnin[1]

#get all tracer names
allTracers<-sort(names(nc_list[[r]]$var))
temp<-allTracers[grep("_N",allTracers)]; tracers2plot<-temp[grep("Nums",temp,invert = TRUE)]; 
tracers2plot<-c(tracers2plot,"Oxygen","Temp","Si", "NO3")
ntracers<-length(tracers2plot)

dynBoxes<-2:26
## compare nut
temp1a<-ncvar_get(nc_list[[2]],"NO3"); temp2a <- ncvar_get(nc_list[[5]], "NO3")
temp1b<-ncvar_get(nc_list[[2]],"NH3"); temp2b <- ncvar_get(nc_list[[5]], "NH3")
temp1 <- temp1a + temp1b; temp2 <- temp2a + temp2b

temp1mean<-apply(temp1[,dynBoxes,], 3, nonZeroMean); temp2mean<-apply(temp2[,dynBoxes,], 3, nonZeroMean); thisMax<-max(c(max(temp1mean, na.rm=TRUE), max(temp2mean, na.rm=TRUE)))
plot(temp2mean, type="l", xaxt="n", ylim=c(8,thisMax)); axis(at=xLabsAt,labels=xLabs,side=1)
points(temp1mean,type="l",col="red")
legend(legend=c(this_out[1], this_out[2]), col=c("black","red"), lty=1, x="topright")


temp3<-ncvar_get(nc_list[[4]],"Si"); tempReducedSi <- ncvar_get(nc_list[[3]], "Si")
temp3mean<-apply(temp3[,dynBoxes,], 3, sum, na.rm=TRUE); temp4mean<-apply(temp4[,dynBoxes,], 3, sum, na.rm=TRUE); thisMax<-max(c(max(temp3mean, na.rm=TRUE), max(temp4mean, na.rm=TRUE)))
plot(temp3mean, type="l", xaxt="n", ylim=c(8,thisMax)); axis(at=xLabsAt,labels=xLabs,side=1)
points(temp4mean,type="l",col="red")
legend(legend=c(this_out[1], this_out[3]), col=c("black","red"), lty=1, x="bottomright")

climateSi <- ncvar_get(nc_list[[4]], "Si")


dynBoxes<-2:24

testDiff <- apply(temp1[-nlayers,dynBoxes,]-temp2[-nlayers, dynBoxes,],3,mean, na.rm=TRUE)
plot(testDiff, type="l")
testDiff <- apply(temp1[(nlayers-1),dynBoxes,]-temp2[(nlayers-1), dynBoxes,],2,mean, na.rm=TRUE)
plot(testDiff, type="l")


SST1<-apply(temp1[(nlayers-1),dynBoxes,],2,mean, na.rm=TRUE); SST2<-apply(temp2[(nlayers-1),dynBoxes, ], 2, mean, na.rm=TRUE)
timeIndex<-675:825
plot(SST1[timeIndex], type="l"); points(SST2[timeIndex], type="l", col=myBlue)


storeTracers<-array(NA, dim=c(nruns, length(tracers2plot), max(nts_list)+1))

plotsFile<-paste(plotPath,"ALL_N.pdf",sep="")
pdf(plotsFile)
par(mfrow=c(4,1),mar=c(3,4,2,0),oma=c(1,0,0,0))
for(t in 1:ntracers){
  thisTracer<-tracers2plot[t]
  temp<-ncvar_get(nc_list[[1]],thisTracer)
  thisVol<-ncvar_get(nc_list[[1]],"volume"); nlayers<-dim(thisVol)[1]
  if(length(dim(temp))==3){
    yy<-apply(temp*thisVol,3,sum) * mg_2_tonne * X_CN
  } else{
    yy<-apply(temp*thisVol[nlayers,,],2,sum) * mg_2_tonne * X_CN
  }
  xx<-yy[burnin[r]:length(yy)]
  # storeTracers[1, t, burnin[r]:length(yy)]<- xx
  # storeTracers[1, t, ]<- xx
  thisymax<-max(xx)*1.1
  thisymin<-min(0,min(xx)*1.1)
  plot(xx,type="l",col=runCols[1],lwd=2,ylim=c(thisymin,thisymax*1.5),ylab="Biomass (tonnes)",xlab="Day",xaxt="n")
  mtext(thisTracer,side=3,adj=0,font=2)
  abline(h=1,col="red",lty=2,lwd=1.5)
  axis(at=xLabsAt,labels=xLabs,side=1)
  
  for(r in 2:nruns){
    temp<-ncvar_get(nc_list[[r]],thisTracer)
    thisVol<-ncvar_get(nc_list[[r]],"volume")
    if(length(dim(temp))==3){
      yy<-apply(temp*thisVol,3,sum) * mg_2_tonne * X_CN
    } else{
      yy<-apply(temp*thisVol[nlayers,,],2,sum) * mg_2_tonne * X_CN
    }
    xx<-yy[burnin[r]:length(yy)]
    points(xx,type="l",col=runCols[r],lwd=1.5,lty=r)
    
    # storeTracers[r, t, burnin[r]:length(yy)]<- xx
    # storeTracers[r, t, ]<- xx
    
    # legend(legend=this_out,col=runCols,lty=seq(1,nruns),x="bottomleft")
  }
}

dev.off()

pdf(paste(plotPath,"_LEGEND.pdf", sep=""), height=5, width=10)
makeBlankPlot()
legend(legend=this_out,col=runCols,lty=seq(1,nruns),x="center", seg.len=3, lwd=3)
dev.off()

## summarise difference
ts_index<-(75*9):dim(storeTracers)[3]
t<-grep("Si", tracers2plot)
diffProp<-(storeTracers[3,,ts_index] - storeTracers[4,,ts_index])/ storeTracers[4,,ts_index]
meanByTracer<-apply(diffProp, 1, median, na.rm=TRUE)
plot(diffProp, type="l")
maxByTracer<-apply(diffProp, 1, max, na.rm=TRUE)
plot(maxByTracer)

tracersByMean <- tracers2plot[order(meanByTracer)]
par(lend=1)
plot(meanByTracer[order(meanByTracer)]*100, type="h", lwd=5, xaxt="n", xlab="", ylab="Percentage change")
par(las=2)
axis(at=seq(1,length(tracersByMean)), labels = gsub("_N|_"," ",tracersByMean), side=1, cex.axis=0.8)
yy<-grep("Temp",tracersByMean)
abline(v=yy, col=myOrange_trans, lwd=5)

#box
pdf(paste(plotPath,"TemperatureTesingPercChange_boxplots.pdf", sep=""), height=5, width=8)
par(mar=c(8,4,1,1))
thisMax<-150; thisMin <- -50
plot(meanByTracer[order(meanByTracer)]*100, type="n", xaxt="n", xlab="", ylab="Percentage change", ylim=c(thisMin, thisMax))
abline(v=grep("Temp", tracersByMean), col=myOrange_trans, lwd=5)
abline(h=meanByTracer[order(meanByTracer)][grep("Temp", tracersByMean)]*100, col=myOrange_trans, lwd=5)
abline(h=0, col="red", lwd=1.5, lty=3)
for(t in 1:length(tracers2plot)){
  thisData<-diffProp[order(meanByTracer)[t],]*100
  boxplot(thisData, at=t, add=TRUE, outline=FALSE, lwd=1.2, yaxt="n")
}
par(las=2)
axis(at=seq(1,length(tracersByMean)), labels = gsub("_N|_"," ",tracersByMean), side=1, cex.axis=0.8)
yy<-grep("Temp",tracersByMean)
dev.off()

# this_out<-c("TestTempNO3Si5day","TestNO3Temp5day", "TestTemp5day", "TestTempBase")

t <- grep("Benthic_Carn", tracers2plot)
t<-grep("NO3", tracers2plot)
thisData1<-ncvar_get(nc_list[[1]],"Si"); thisData2<-ncvar_get(nc_list[[2]], "Si")
xx1<-apply(thisData1[,dynBoxes,], 3, nonZeroMean)
xx2<-apply(thisData2[,dynBoxes,], 3, nonZeroMean)

plot(xx1[100:length(xx2)],type="l", ylim=c(0,max(c(xx1,xx2))))
points(xx2[100:length(xx2)],type="l", lty=2, col="red")

thisData<-storeTracers[,t,]
plot(thisData[1,], type="l"); points(thisData[2,], type="l", lty=2, col=myRed)
plot(diffProp[t,], type="l")


testNC.nc<-nc_open(paste(this_path,"outputTestNO3Time\\output.nc",sep=""))
testNO3time<-ncvar_get(testNC.nc, "NO3")
dim(testNO3time)
yy <- apply(testNO3time[,dynBoxes,], 3, mean, na.rm=TRUE)

## read in TLs so can check out these by TL


# 
# index<-grep("Seaperch", tracers2plot)
# test<-storeTracers[,index,]
# 
# }

# 
# temp<-ncvar_get(nc_list[[1]],thisTracer); temp2<-ncvar_get(nc_list[[2]], thisTracer)
# test<-apply(temp,c(2,3), sum)
# test2<-apply(temp2,c(2,3), sum)
# 
# x_rs_fn<-function(r){
#   x=r/(1-r)
#   return(x)
# }
# test_r<-seq(0.6,0.8,length.out = 20)
# test_x<-unlist(lapply(test_r, x_rs_fn)); mean(test_x)
# plot(x=test_r, y=test_x, pch=20, col=myBlue_trans)
# abline(h=x_rs_fn(0.7), col=myGrey_trans, lwd=5); abline(h=2.75, col=myGrey_trans, lwd=5)
# abline(h=2.5, col=myGrey)
# # 
# # 
# # 
# # ##testing
# # testVar<-"Ben_fish_shal1_Nums"
# # temp1<-ncvar_get(nc_list[[1]],testVar)
# # xx1<-apply(temp1,3,sum)
# # yy1<-xx1/xx1[1]
# # temp2<-ncvar_get(nc_list[[2]],testVar)
# # 
# # xx2<-apply(temp2,3,sum)
# # yy2<-xx2/xx2[1]
# # 
# #                  