bgmFileName<-"NGOM"

source(paste(DIR$'General functions',"read_boxes.R",sep=""))

bgm<-read_boxes(fnm=paste(DIR$'Base',"\\ATLANTISmodels\\inputs\\bgm\\",bgmFileName,".bgm",sep=""))
numpolys<-bgm$nbox

yMin<-1e+20; yMax<-0; xMin<-0; xMax<--1e+20
##get lims for plotting
for(b in 1:numpolys){
  this_box<-bgm$verts[[b]]
  x<-unlist(str_split(this_box," "))
  this_x<-as.double(x[seq(1,length(x),by=2)])
  this_y<-as.double(x[seq(2,length(x),by=2)])
  thisXmin<-min(this_x); thisXmax<-max(this_x); thisYmin<-min(this_y); thisYmax<-max(this_y)
  if(thisXmin<xMin){xMin<-thisXmin}
  if(thisYmin<yMin){yMin<-thisYmin}
  if(thisYmax>yMax){yMax<-thisYmax}
  if(thisXmax>xMax){xMax<-thisXmax}
  
}

thisScale<-seq(0,150,by=10)

colRamp<-colorRampPalette(colors=c(myOrange,myYellow,myGreen,myAqua,myBlue,myPurple))
thisColors<-colRamp(length(thisScale))

depthScale_byPoly<-abs(as.numeric(unlist(bgm$b_botz)))  
color_byPoly<-thisColors[match(depthScale_byPoly,thisScale)]



#do one that is transparent with no background
plotType="png"
saveName<-paste(DIR$'Figures',"\\DepthByPolygon.",plotType,sep="")
new.graph(0.5,quiet=TRUE,filename=saveName,type=plotType)
par(mar=c(12,5,0,0),oma=c(0,1,0,0),bg = "transparent")
plot(x=NULL,y=NULL,xlim=c(xMin,xMax),ylim=c(yMin,yMax),type="n",xaxt="n",yaxt="n",asp=TRUE,xlab="",ylab="",bty="n")
for(b in 1:numpolys){
  this_box<-bgm$verts[[b]]
  x<-unlist(str_split(this_box," "))
  this_x<-as.double(x[seq(1,length(x),by=2)])
  this_y<-as.double(x[seq(2,length(x),by=2)])
  thisCol<-color_byPoly[b]
  polygon(x=this_x,y=this_y,col=thisCol)  
} 
par(xpd=TRUE)
legend(legend=thisScale,title="Depth (m)",pch=15,cex=0.65,pt.cex=7,bty="n",col=thisColors,x="bottom",adj=c(1,1),inset=-0.1,horiz=TRUE) 
par(xpd=FALSE)
dev.off()

thisProjection<-"+proj=utm +a=6378137.0 +es=0.006694380022900787 +lon_0=-81d00 +lat_0=0d00 +x_0=500000.0 +y_0=0.0 +k=0.9996 +zone=17"

TBGBprojection<-"+proj=tmerc +a=6378137.0 +es=0.006694380022900787 +lon_0=173d00 +lat_0=0d00 +x_0=1600000.0 +y_0=1.0E7 +k=0.9996"
         EAproj="+proj=tmerc +lat_0=0 +lon_0=173 +k=0.9996 +x_0=1600000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs " #EPSG:2193

"wgs84"


for(b in 1:numpolys){
  thisPoly<-bgm$verts[[b]]
  thisNumPoints<-length(thisPoly)
  for(p in 1:thisNumPoints){
    xx<-unlist(str_split(thisPoly[p]," "))
    thisPoint<-c(as.double(xx[1]),as.double(xx[2]))
    this_llPoint<-unlist(project(xy=thisPoint,thisProjection,inverse=TRUE))
    projShape$shp$shp[[s]]$points[p,]<-this_aeaPoint
  }
  #replace x and y min and max too
  thisYmin<-min(projShape$shp$shp[[s]]$points$Y)
  thisYmax<-max(projShape$shp$shp[[s]]$points$Y)
  thisXmin<-min(projShape$shp$shp[[s]]$points$X)
  thisXmax<-max(projShape$shp$shp[[s]]$points$X) 
  projShape$shp$shp[[s]]$box<-c("xmin"=thisXmin,"ymin"=thisYmin,"xmax"=thisXmax,"ymax"=thisYmax)  
}
