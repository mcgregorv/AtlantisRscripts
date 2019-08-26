library(ggmap)
library(rgdal)
library(gtable)
library(maps)
library(mapdata)
library(rgeos)
source(paste(DIR$'General functions',"\\get_first_number.R",sep=""))
source(paste(DIR$'General functions',"formatShape.R",sep=""))

this_run<-"base"
this_out<-"TEST150yrfish"

burnin<-35 #number of years to skip in plot

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")
baseOutPath<-paste(this_path,"output",base_out,"\\", sep="")

################################
#read in satelite data
months<-c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
nm<-length(months)
chl_a_all<-data.frame(matrix(NA,ncol=4,nrow=0))
colnames(chl_a_all)<-c("lon","lat","chl_a","month")
for(m in 1:nm){
  thisMonth<-months[m]
  thisMonthNC<-nc_open(paste(DIR$'Data',"chlorophyll\\",thisMonth,"_chl.nc",sep=""))
  monthTracers<-sort(names(thisMonthNC$var))
  this_chlor_a<-ncvar_get(thisMonthNC,"chlor_a")
  thisLon<-ncvar_get(thisMonthNC,"lon")
  thisLat<-ncvar_get(thisMonthNC,"lat")
  rownames(this_chlor_a)<-thisLon; colnames(this_chlor_a)<-thisLat
  chl_a_df<-adply(this_chlor_a,c(1,2))
  colnames(chl_a_df)<-c("lon","lat","chl_a")
  chl_a_df$month<-thisMonth  
  chl_a_all<-rbind(chl_a_all,chl_a_df)
}

thisMax<-max(chl_a_all$chl_a,na.rm=TRUE)

##assign to polygon
source(paste(DIR$'General functions',"findMyPolygon.R",sep=""))

colRamp<-colorRampPalette(colors=c(myYellow,myOrange,myRed,myPurple,myBlue,myGreen))

#get box boundaries. this brings in dynamicBoxes
thisFile<-paste(this_path,"..\\inputs\\bgm\\Chatham_Boxes.R",sep="")
source(thisFile)

nboxes<-length(dynamicBoxes)

names(dynamicBoxes[[1]])
# "x"        "y"        "minDepth" "maxDepth"

chl_a_all$box<-NA

# test_chl_a_all<-chl_a_all[1:100,]
# test_chl<-chl_a_all
# test_chl$lon<-as.double(as.character(test_chl$lon))
# test_chl$lat<-as.double(as.character(test_chl$lat))
# test<-mapply(FUN=findMyPolygon,x=test_chl$lon,y=test_chl$lat,boxes=dynamicBoxes,nboxes)

for(i in 1:(nrow(chl_a_all))){
  xx<-findMyPolygon(x=as.double(as.character(chl_a_all$lon[i])),y=as.double(as.character(chl_a_all$lat[i])),boxes=dynamicBoxes,nboxes)
  if(length(xx)==1){
    chl_a_all$box[i]<-xx
  }
}
index<-chl_a_all$chl_a==-9999
chl_a_all$chl_a[index]<-NA
monthsOrder<-cbind(months,seq(1,12))
chl_by_month<-tapply(chl_a_all$chl_a,chl_a_all[,c("month")],mean,na.rm=TRUE)
chl_by_box<-tapply(chl_a_all$chl_a,chl_a_all[,c("box")],mean,na.rm=TRUE)
chl_by_boxMonth<-tapply(chl_a_all$chl_a,chl_a_all[,c("month","box")],mean,na.rm=TRUE)
plotByMonths<-as.double(chl_by_month[match(months,names(chl_by_month))])
plot(plotByMonths,type="h",xaxt="n",lend=1,lwd=5,col=myGrey_trans,xlab="",ylab="Chl_a")
axis(at=seq(1,12),labels=months,side=1,las=2)

plot(chl_by_box,type="h",lend=1,lwd=5,col=myGrey_trans,xlab="Box number",ylab="Chl_a")


#######################
#read in tracers file
thisBaseOut<-paste("Short",sep="")
# thisFishOut<-paste("FishShort",sep="")

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
baseOutPath<-paste(this_path,"output",thisBaseOut,"\\",sep="")

plotPath<-paste(this_path,"..\\Figures\\",this_run,"\\",this_out,"",sep="")

plotBurnIn<-FALSE

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1900 #this is the first year of the burn-in part
fishingStartYear<-1900 #first year historical catches are removed
modelStartYear<-1920 

#read in the base tracers file
BaseNC.nc<-nc_open(paste(baseOutPath,"output.nc",sep=""))
baseVol<-ncvar_get(BaseNC.nc,"volume")
baseDz<-ncvar_get(BaseNC.nc,"dz")
base_nts<-dim(baseVol)[3]

chl_a_tracer<-ncvar_get(BaseNC.nc,"Chl_a")
chl_tracer_byBox<-apply(chl_a_tracer,2,mean,na.rm=TRUE)
chl_tracer_BoxIC<-apply(chl_a_tracer[,,1],2,mean,na.rm=TRUE)
chl_tracer_BoxEnd<-apply(chl_a_tracer[,,(base_nts-1)],2,mean,na.rm=TRUE)


par(mfrow=c(4,3),mar=c(3,3,0,0))
for(m in 1:nm){
  thisMonth<-months[m]
  thisData<-chl_by_boxMonth[thisMonth,]
  plot(thisData,type="h",lwd=5,lend=1,col=myRed)
  mtext(thisMonth,side=3,adj=0,line=-1)
  par(new=TRUE)
  plot(chl_tracer_byBox[2:25],type="l",lwd=2)
}
for(m in 1:nm){
  plot(x=chl_tracer_byBox[2:25],y=chl_by_box,ylim=c(0,1),xlim=c(0,1))
  points(x=chl_tracer_byBox[2:25],chl_by_boxMonth[m,],col=myOrange,pch=8)
  mtext(rownames(chl_by_boxMonth)[m])
  points(x=c(0,1),y=c(0,1),col="red",lty=2,type="l")
}

for(m in 1:nm){
  plot(x=chl_tracer_BoxEnd[2:25],y=chl_by_box,ylim=c(0,1),xlim=c(0,1))
  points(x=chl_tracer_BoxEnd[2:25],chl_by_boxMonth[m,],col=myGreen,pch=8)
  mtext(rownames(chl_by_boxMonth)[m])
  points(x=c(0,1),y=c(0,1),col="red",lty=2,type="l")
}

#################### plot spatially #################
chl_colors<-rev(rainbow(n=11,start=0.05,end=0.8))

chl_max<-max(chl_tracer_byBox,na.rm=TRUE)

getCol<-function(x){
  y<-round((log(10*x))/(log(10*chl_max)),1)*10
  # y<-min(1,y); 
  y<-max(0,y)
  ycol<-chl_colors[y]
  return(ycol)
}
chl_col_tracers<-unlist(lapply(chl_tracer_byBox,FUN=getCol))

plot(seq(1,11),col=rev(chl_colors),pch=8)
#read in shape file
shapeFile<-paste(DIR$'Base',"ATLANTISmodels\\inputs\\bgm\\CHAT30_LL",sep="")
sdata<-read.shapefile(shapeFile)
shape<-formatShape(shapeFile=shapeFile)
ns<-length(shape)
SpDF <- SpatialPolygonsDataFrame(shape,data.frame( z=1:ns, row.names=paste("P",seq(1,(ns)),sep="")))
labels<-seq(1,(ns))
pdf("test.pdf")
plot(shape)
LABELpOS<-polygonsLabel(shape, labels = labels, cex=.1,doPlot=FALSE)
dev.off()
labeldf<-data.frame(cbind("x"=LABELpOS[1:ns],"y"=LABELpOS[(ns+1):(2*ns)]))

# pdf(paste(plotPath,"plotChl_a.pdf",sep=""))
# par(mfrow=c(3,2),mar=c(0,0,0,0),oma=c(0,0,0,0))
# for(t in 1:nts){
#   timeData<-thisData[,t]
#   timeColors<-unlist(lapply(timeData,getCol))
plot(shape)
map('nzHires',add=TRUE,col="black",lwd=2)
map.axes()

  for(plotB in 1:dim(labeldf)[1]){
    polygon(sdata$shp$shp[[plotB]]$points,col=chl_col_tracers[plotB-1])
  }


plot(chl_tracer_BoxIC[2:25],type="h",lwd=3,col=myBlue,lend=1)
points(chl_by_box,type="h",lwd=7,col=myGrey_trans)

plot(x=seq(1,24,by=1),y=chl_by_box,type="h",lend=1,lwd=5,col=myGrey_trans,xlab="Box number",ylab="Chl_a")
points(x=seq(1.5,24.5,by=1),y=chl_tracer_byBox[2:25],type="h",lend=1,lwd=5,col=myBlue)

  
  
  

################################
# 
# 
# this_nbox<-length(sdata$shp$shp)
# boxes<-seq(1,this_nbox)
# 
# groupsDF<-read.csv(paste(this_path,"CRAM_groups.csv",sep=""))
# 
# if(thisGroup %in% groupsDF$Code){
#   thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisGroup],side="both")
# } else{
#   thisName<-thisGroup; thisVar<-thisGroup
# }
# 
# #read in nc file
# ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
# thisVol<-ncvar_get(ThisNC.nc,"volume")
# thisDz<-ncvar_get(ThisNC.nc,"dz")
# 
# nts<-dim(thisVol)[3] #number of timesteps
