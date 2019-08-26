#First do overall biomass, only of verts
thisPath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
groupsDF<-read.csv(paste(thisPath,"CRAM_groups.csv",sep=""))
vertsDF<-groupsDF[groupsDF$NumCohorts>=2,]
nv<-dim(vertsDF)[1]

this_file<-paste(thisPath,"inputs\\biol_initial\\CRAM_Initial_biomass_distributions.xlsx",sep="")
wb <- loadWorkbook(this_file)
sheetNames<-getSheets(wb)

nb<-24 #number of dynamic boxes

vertB0<-NULL
vertB0byBox<-NULL
for(v in 1:nv){
  thisGroup<-vertsDF$Code[v]
  sheetIndex<-grep(thisGroup,sheetNames)
  thisBiolData<-readWorksheet(wb,sheet=sheetIndex) 
  thisB0<-as.double(thisBiolData[2,2])
  vertB0[[as.character(thisGroup)]]<-thisB0
  thisByBox<-as.double(thisBiolData[5:(5+nb-1),2])*thisB0
  vertB0byBox[[as.character(thisGroup)]]<-thisByBox
}

B0sum<-sum(vertB0,na.rm=TRUE)
B0Max<-max(vertB0,na.rm=TRUE)

cols<-colorRampPalette(colors=c(myYellow,myGreen,"green",myAqua,myBlue,"black",myGrey,myPurple,myRed,"red","brown",myOrange,"lemon"))(nv)
vert2col<-data.frame(cbind(names(vertB0),cols))
colnames(vert2col)<-c("vert","col")

pdf(paste(DIR$'Figures',"biomass\\B0_vertebrates.pdf",sep=""))
plot(x=1,y=1,type="n",xlab="",ylab="Biomass (tonnes)",xaxt="n",xaxt="n",xlim=c(0.5,(nv+0.5)),ylim=c(0,B0Max))

par(lend=1)
for(v in 1:nv){
  thisB0<-vertB0[v]
  points(x=v,y=thisB0,type="h",lwd=5,col=as.character(vert2col$col[vert2col$vert==names(thisB0)]))
}
par(las=2)
axis(at=seq(1,nv),labels = names(vertB0),side=1,cex=0.2)
dev.off()


#plot by box, as proportions, then split out individually
#first get sum by box, and max by box
B0sumByBox<-rep(0,nb)
B0maxByBox<-rep(0,nb)
for(b in 1:nb){
  thisSum<-0
  thisMax<-0
  for(v in 1:nv){
    thisSum<-thisSum+vertB0byBox[[v]][b]
    if(vertB0byBox[[v]][b]>thisMax){thisMax<-vertB0byBox[[v]][b]}
  }
  B0sumByBox[b]<-thisSum
  B0maxByBox[[b]]<-thisMax
}

pdf(paste(DIR$'Figures',"biomass\\B0_vertebrates_proportionsByBox.pdf",sep=""))
par(oma=c(7,1,0,1))
plot(x=1,y=1,type="n",xlab="Polygon",ylab="Proportion of biomass",xlim=c(0.5,(nb+0.5)),ylim=c(0,1))
for(b in 1:nb){
  thisX<-c(b-0.5,b-0.5,b+0.5,b+0.5)
  prevProp<-0
  for(v in 1:nv){
    thisVertCode<-names(vertB0)[[v]]
    thisCol<-as.character(vert2col$col[vert2col$vert==thisVertCode])
    thisProp<-prevProp+vertB0byBox[[v]][b]/B0sumByBox[b]
    thisY<-c(prevProp,thisProp,thisProp,prevProp)
    polygon(x=thisX,y=thisY,col=thisCol,border=NA)
    prevProp<-thisProp
  }
}
par(xpd=NA)
legend(legend=vert2col$vert,col=as.character(vert2col$col),lwd=4,seg.len=1,x="bottom",inset=-0.6,ncol=7)
dev.off()

library(maps)
library(mapdata)
library(rgeos)
source(paste(generalFunctionsPath,"formatShape.R",sep=""))
#read the shape file
graphics.off()
shapeFile<-paste(DIR$'Base',"ATLANTISmodels\\inputs\\bgm\\CHAT30_LL",sep="")
sdata<-read.shapefile(shapeFile)
shape<-formatShape(shapeFile=shapeFile)
ns<-length(shape)
SpDF <- SpatialPolygonsDataFrame(shape,data.frame( z=1:ns, row.names=paste("P",seq(1,(ns)),sep="")))
labels<-seq(1,(ns))
plot(shape)
LABELpOS<-polygonsLabel(shape, labels = labels, cex=.1,doPlot=FALSE)
labeldf<-data.frame(cbind("x"=LABELpOS[1:ns],"y"=LABELpOS[(ns+1):(2*ns)]))

#a plot for each box
pdf(paste(DIR$'Figures',"biomass\\B0_vertebratesForEachBox.pdf",sep=""))
par(mfrow=c(3,2),lend=1,mar=c(4,4,0,0))
for(b in 1:nb){
  
  plot(x=1,y=1,type="n",xlab="",ylab="Biomass (tonnes)",xaxt="n",xaxt="n",xlim=c(0.5,(nv+0.5)),ylim=c(0,B0maxByBox[b]))
  
  for(v in 1:nv){
    thisB0<-vertB0byBox[[v]][b]
    points(x=v,y=thisB0,type="h",lwd=5,col=as.character(vert2col$col[vert2col$vert==names(vertB0byBox)[v]]))
  }
  par(las=2)
  axis(at=seq(1,nv),labels = names(vertB0),side=1,cex=0.2,font=2)
  
  #plot map
  plot(shape)
  map('nzHires',add=TRUE,col="black",lwd=2)
  map.axes()
  for(plotB in 1:dim(labeldf)[1]){
    if(plotB==(b+1)){
      polygon(sdata$shp$shp[[plotB]]$points,col=myGreen)
    }
    text((plotB-1),x=labeldf$x[plotB],y=labeldf$y[plotB],col="black",cex=0.8,font=2)
  }
  
}
dev.off()

###########################
## how does this one look if I read in from atlantis initial conditions file instead?
## first read in and store values for each vert
mg2tonnes<-2e-8
X_CN<-5.7

getBoxTotal<-function(x){
  xx<-gsub(";","",x)
  thisVec<-as.double(unlist(str_split(xx,",")))
  thisBoxSum<-sum(thisVec,na.rm=TRUE)
  return(thisBoxSum)
}

biolInLines<-readLines(paste(thisPath,"inputs\\biol_initial\\CRAM_input.txt",sep=""))
#get box volumes first
xx<-grep("volume =",biolInLines)
yy<-grep(";",biolInLines)
this_yy<-min(yy[yy>xx])
thisInitialData<-biolInLines[(xx+1):(this_yy-1)]
volByBox<-unlist(lapply(thisInitialData,getBoxTotal))

BI_vertsByBox<-NULL
for(v in 1:nv){
  thisVertCode<-as.character(vert2col$vert[v])
  thisVertName<-str_trim(vertsDF$Name[vertsDF$Code==thisVertCode],side="both")
  thisVertVar<-paste(thisVertName,"_N =",sep="")
  xx<-grep(thisVertVar,biolInLines)
  this_yy<-min(yy[yy>xx])
  thisInitialData<-biolInLines[(xx+1):(this_yy-1)]
  thisNbBox<-unlist(lapply(thisInitialData,getBoxTotal))
  thisTonnesByBox<-thisNbBox*volByBox*mg2tonnes*X_CN
  BI_vertsByBox[[v]]<-thisTonnesByBox
}

#first get sum by box, and max by box
BIB0sumByBox<-rep(0,nb)
BIB0maxByBox<-rep(0,nb)
for(b in 1:nb){
  thisSum<-0
  thisMax<-0
  for(v in 1:nv){
    thisSum<-thisSum+BI_vertsByBox[[v]][b]
    if(BI_vertsByBox[[v]][b]>thisMax){thisMax<-BI_vertsByBox[[v]][b]}
  }
  BIB0sumByBox[b]<-thisSum
  BIB0maxByBox[[b]]<-thisMax
}

pdf(paste(DIR$'Figures',"biomass\\B0_vertebratesForEachBox_fromInput.pdf",sep=""))
par(mfrow=c(3,2),lend=1,mar=c(4,4,0,0))
for(b in 1:nb){
  
  plot(x=1,y=1,type="n",xlab="",ylab="Biomass (tonnes)",xaxt="n",xaxt="n",xlim=c(0.5,(nv+0.5)),ylim=c(0,BIB0maxByBox[b]))
  
  for(v in 1:nv){
    thisB0<-BI_vertsByBox[[v]][b]
    points(x=v,y=thisB0,type="h",lwd=5,col=as.character(vert2col$col[vert2col$vert==names(vertB0byBox)[v]]))
  }
  par(las=2)
  axis(at=seq(1,nv),labels = names(vertB0),side=1,cex=0.2,font=2)
  
  #plot map
  plot(shape)
  map('nzHires',add=TRUE,col="black",lwd=2)
  map.axes()
  for(plotB in 1:dim(labeldf)[1]){
    if(plotB==b){
      polygon(sdata$shp$shp[[plotB]]$points,col=myGreen)
    }
    text((plotB-1),x=labeldf$x[plotB],y=labeldf$y[plotB],col="black",cex=0.8,font=2)
  }
  
}
dev.off()


