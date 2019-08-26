#read in xlsx file from Peter and plot the polygons
thisFile<-paste(DIR$'Base',"ATLANTISmodels\\inputs\\Chatham_Boxes.xlsx",sep="")
nboxes<-24 ##dynamic boxes
nboundaries<-6 ##non-dynamic, boundary boxes

dynamicBoxes<-NULL

for(b in 1:nboxes){
  thisBox<-seq(1,nboxes)[b]
  thisSheetName<-paste("Box",thisBox,sep="")
  temp<-read.xlsx2(thisFile,sheetName=thisSheetName,header=FALSE)
  thisDepth<-temp[1,3]
  if(substr(thisDepth,start=1,stop=1)=="<"){
    thisMinDepth<-0
    thisMaxDepth<-as.double(unlist(str_split(thisDepth,"<"))[2])
  } else{
    temp<-unlist(str_split(thisDepth,"-"))
    thisMinDepth<-as.double(temp[1])
    thisMaxDepth<-as.double(temp[2])
  }
  thisBoxData<-read.xlsx(thisFile,sheetName=thisSheetName,header=TRUE,startRow=4)
  dynamicBoxes[[b]]<-list()
  dynamicBoxes[[b]]$x<-thisBoxData$"lon.dec"
  dynamicBoxes[[b]]$y<-thisBoxData$'lat.dec'
  dynamicBoxes[[b]]$minDepth<-thisMinDepth
  dynamicBoxes[[b]]$maxDepth<-thisMaxDepth
}

#dump the boxes so can read in easily again
thisFile<-paste(DIR$'Base',"ATLANTISmodels\\inputs\\Chatham_Boxes.R",sep="")
dump("dynamicBoxes",file=thisFile)

new.graph(0.5)
nz(xlim = c(172,187),ylim=c(-42,-46.2),fill.col="DarkOliveGreen3",mgp=c(3,0.5,0))
for(b in 1:nboxes){
  thisY<-dynamicBoxes[[b]]$y
  thisX<-dynamicBoxes[[b]]$x
  nz.polygon(x=thisX,y=thisY,border=myOrange,lwd=2)
}
SavePlot(paste(DIR$'Figures',"Chatham_Rise_Boxes",sep=""))

