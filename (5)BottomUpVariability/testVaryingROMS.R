#read in tracers and check out temperature to see if correct ROMS being used
this_run<-"base"
this_out<-paste("BASE50yrTestROMS5day",sep="")

basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
ROMSpath<-paste(basePath,"inputs\\ROMS\\",sep="")

ThisNC.nc<-nc_open(paste(basePath,"base\\output",this_out,"\\output.nc", sep=""))
thisTemp<-ncvar_get(ThisNC.nc, "Temp")

testTemp<-apply(thisTemp,3,mean); test_nts<-10*length(testTemp)

#test roms nc files for temperature
# par(mfrow=c(3,1))
plotTempByYear<-NULL
plot(plotTemp,type="n",ylim=c(5,12))
for(y in 1:3){
  thisYear<-c("2003", "1996","1999")[y]
  # thisYear<-"1996"
  thisRoms<-nc_open(paste(ROMSpath,"Chatham_temp_year",thisYear,".nc",sep=""))
  thisTemp<-ncvar_get(thisRoms,"temperature")
  plotTempByYear[[y]]<-apply(thisTemp,3,mean,na.rm=TRUE)
  # 
  # plotTemp<-apply(thisTemp,3,mean, na.rm=TRUE)
  # points(plotTemp,type="l",col=c(myGold,myRed,myBlue)[y])
  # mtext(thisYear,side=3)
}
legend(legend=c("2003","1996","1999"), col=c(myGold,myRed,myBlue), lwd=2,lty=1,x="bottomleft",bty="n")

thisPlotTemp<-c(plotTempByYear[[1]],plotTempByYear[[2]],plotTempByYear[[2]],plotTempByYear[[3]],plotTempByYear[[3]])
this5dayTemp<-thisPlotTemp[seq(1,trunc(length(thisPlotTemp)/5),by=5)]
par(mfrow=c(2,1))
plot(this5dayTemp)
plot(testTemp)

plot(thisPlotTemp[1:test_nts])
