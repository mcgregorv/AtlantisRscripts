## what happens with BO?? - first check out what they eat most of

this_run<-"Chaos"
this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")

thisOut <- "outputCChaosdown05"; 
thisDietCheck <- read.csv(paste(this_path, thisOut,"\\outputDietCheck.txt", sep=""), sep=" ")

BOdiets <- thisDietCheck[thisDietCheck$Predator=="BO",]
BObyprey <- apply(BOdiets[,c(6:dim(BOdiets)[2])],2, sum, na.rm=TRUE)
sort(BObyprey, decreasing = TRUE)

# DR          MA          BB          BD          MB          DL          ZS          ZM         ASQ         BAL          BC         BEE         BFF 
# 60.31567034 50.29447631 19.60274262 12.48157225  6.16249611  0.84280152  0.27593408  0.02430709  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000 

MBdiets <- thisDietCheck[thisDietCheck$Predator=="MB",]
BObyprey <- apply(BOdiets[,c(6:dim(BOdiets)[2])],2, sum, na.rm=TRUE)
sort(BObyprey, decreasing = TRUE)

## MB predators
xx<-(1:dim(thisDietCheck)[2])[colnames(thisDietCheck)=="MB"]
temp<-thisDietCheck[,c(1:4,xx)]
index<-temp[,c("MB")]>0
MBpredData<-temp[index,]

predByYear<-tapply(MBpredData$MB, MBpredData$Time, sum, na.rm=TRUE)

toPlot<-MBpredData[MBpredData$Time < 10000,]
toPlot<-MBpredData

bp<-ggplot(data = toPlot, aes(x = Time, fill = Predator, y = MB)) + 
  geom_bar(stat = 'identity')
bp

#################
## check out depth of sediment in box 3
thisNC.nc <- nc_open(paste(this_path, thisOut,"\\output.nc", sep=""))
thisData<-ncvar_get( thisNC.nc, "dz")
thisData<-ncvar_get(thisNC.nc, "Light")

thisBoxes<-dynBoxes
thisBoxes<-c(2,3,9)
test<-apply(thisData[6,thisBoxes,], 2,mean, na.rm=TRUE)

boxCol<-colorRampPalette(colors=c("red",myYellow,myBlue))(length(dynBoxes))
plot(test[1:65], type="l", ylim=c(0,max(thisData[6,dynBoxes,1:65])))
points(thisData[6,3,1:65], type="l", col=myOrange, lwd=2)

for(b in dynBoxes){
  par(new=TRUE)
  plot(thisData[6,b,1:65], type="l", lty=b, col=boxCol[b])
}
