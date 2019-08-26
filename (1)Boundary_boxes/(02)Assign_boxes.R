source(paste(DIR$'General functions',"findMyPolygon.R",sep=""))

colRamp<-colorRampPalette(colors=c(myYellow,myOrange,myRed,myPurple,myBlue,myGreen))

#assigning boxes to fisheries data
thisFile<-paste(DIR$'Data',"\\Fisheries\\CHAT_Lat_Long.txt",sep="")
thisData<-read.csv(thisFile,sep="\t")
thisData$box<-NA
thisData$start_latitude<-as.double(thisData$start_latitude)
thisData$start_longitude<-as.double(thisData$start_longitude)

#get box boundaries. this brings in dynamicBoxes
thisFile<-paste(DIR$'Base',"ATLANTISmodels\\inputs\\Chatham_Boxes.R",sep="")
source(thisFile)

nboxes<-length(dynamicBoxes)
cols<-colRamp(nboxes)

names(dynamicBoxes[[1]])
# "x"        "y"        "minDepth" "maxDepth"

for(i in 1:(nrow(thisData))){
  xx<-findMyPolygon(x=thisData$start_longitude[i],y=thisData$start_latitude[i],boxes=dynamicBoxes,nboxes)
  if(length(xx)==1){
    thisData$box[i]<-xx
    cat(xx,"--")
  }

}
#write it out
thisFile<-paste(DIR$'Data',"\\Fisheries\\CHAT_Lat_Long_boxes.txt",sep="")
write.csv(thisData,thisFile)

#plot if
new.graph(0.5)
nz(xlim = c(172,187),ylim=c(-42,-46.2),fill.col="DarkOliveGreen3",mgp=c(3,0.5,0))
for(b in 1:nboxes){
  thisY<-dynamicBoxes[[b]]$y
  thisX<-dynamicBoxes[[b]]$x
  nz.polygon(x=thisX,y=thisY,border="black",lwd=2)
  thisIndex<-thisData$box==b
  thisCol<-cols[b]
  nz.points(x=as.double(thisData$start_longitude[thisIndex]),y=as.double(thisData$start_latitude[thisIndex]),col=thisCol,pch=18)
}
thisIndex<-is.na(thisData$box)
nz.points(x=as.double(thisData$start_longitude[thisIndex]),y=as.double(thisData$start_latitude[thisIndex]),col=myGrey,pch=18)
SavePlot(paste(DIR$'Figures',"CHAT_Lat_Long",sep=""))

############################################################################################
############################################################################################
## CHAT_lat_long_FSU ##
############################################################################################
#assigning boxes to fisheries data
thisFile<-paste(DIR$'Data',"\\Fisheries\\CHAT_Lat_Long_FSU.txt",sep="")
thisData<-read.csv(thisFile,sep="\t")
thisData$box<-NA
thisData$start_latitude<-as.double(thisData$lat_s)
thisData$start_longitude<-as.double(thisData$long_s)

for(i in 1:(nrow(thisData))){
  xx<-findMyPolygon(x=thisData$start_longitude[i],y=thisData$start_latitude[i],boxes=dynamicBoxes,nboxes)
  if(length(xx)==1){
    thisData$box[i]<-xx
    cat(xx,"--")
  }
  
}
#write it out
thisFile<-paste(DIR$'Data',"\\Fisheries\\CHAT_Lat_Long_FSU_boxes.txt",sep="")
write.csv(thisData,thisFile)

#plot if
new.graph(0.5)
nz(xlim = c(172,187),ylim=c(-42,-46.2),fill.col="DarkOliveGreen3",mgp=c(3,0.5,0))
thisIndex<-is.na(thisData$box)
nz.points(x=as.double(thisData$start_longitude[thisIndex]),y=as.double(thisData$start_latitude[thisIndex]),col=myGrey,pch=18)
for(b in 1:nboxes){
  thisIndex<-thisData$box==b
  thisCol<-cols[b]
  nz.points(x=as.double(thisData$start_longitude[thisIndex]),y=as.double(thisData$start_latitude[thisIndex]),col=thisCol,pch=18)
  thisY<-dynamicBoxes[[b]]$y
  thisX<-dynamicBoxes[[b]]$x
  nz.polygon(x=thisX,y=thisY,border="black",lwd=2)
}
SavePlot(paste(DIR$'Figures',"CHAT_Lat_Long_FSU",sep=""))












