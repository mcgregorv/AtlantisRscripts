#plot gapesize for given predator with size of given prey
library(ggmap)
library(rgdal)
library(gtable)
library(maps)
library(mapdata)
library(rgeos)
source(paste(DIR$'General functions',"\\get_first_number.R",sep=""))
source(paste(DIR$'General functions',"formatShape.R",sep=""))
source(paste(DIR$'General functions',"nonZeroMean.R",sep=""))
this_run<-"base"

this_out<-"FISH"

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")
plotPath<-paste(this_path,"..\\Figures\\",this_run,"\\",this_out,"",sep="")

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1900
fishingStartYear<-1900
modelStartYear<-1900

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
ntsteps<-dim(thisVol)[3]

burnin=35

nts<-dim(thisVol)[3]-burnin #number of timesteps
cat(paste(nts,"\n"))
xLabsTemp<-seq(0,(nts*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0

thisBiolFile<-paste(this_path,"..\\CRAM_BH_hybrid_biol.prm",sep="")
biolLines<-readLines(thisBiolFile)

groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

#get all tracer names
allTracers<-sort(names(ThisNC.nc$var))

cohortCols<-colorRampPalette(colors=c(myYellow,myGreen,myAqua,myBlue,myPurple,myRed,"red"))(10)
cohortCols_trans<-paste(cohortCols,"88",sep="")
# # legends 
# pdf(paste(plotPath,"GapeSizePreyLegend.pdf", sep=""), height=3, width=3)
# par(mar=c(0,1,0,1), lend=1)
# makeBlankPlot()
# legend(legend=seq(1,10), col=cohortCols, lwd=3, lty=2, seg.len=4, x="left", bty="n", title="Prey cohort")
# legend(legend=seq(1,10), col=cohortCols_trans, lwd=7, lty=1, seg.len=4, x="right", bty="n", title="Predator cohort")
# dev.off()


# can run as a function, but I usually just run for one species - then can edit KLP or KUP here and re-run to check effect
# plotGapeLimits<-function(thisPred,thisPrey){
  thisPred<-"HOK"; 
  # for(thisPred in c("SND", "ELI", "SPE", "LDO", "LIN", "PFL", "CEP", "MAC")){
  
  thisPrey<-"LIN"
  # pdf(paste(plotPath,"Pred",thisPred,"_prey",thisPrey,"_gapeSize.pdf",sep=""))
  #get KLP and KUP for predator
  thisVar<-paste("KLP_",thisPred,sep="")
  thisKLP<-get_first_number(biolLines[grep(thisVar,biolLines)])
  # thisKLP<-2e-5 # overwrite to test alt value
  thisVar<-paste("KUP_",thisPred,sep="")
  thisKUP<-get_first_number(biolLines[grep(thisVar,biolLines)])
  # thisKUP<-0.5  # overwrite to test alt value
  #get SN for each cohort wrt time
  thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisPred],side="both"); thisNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisPred]
  predSNarray<-array(NA,dim=c(nts,thisNumCohorts))
  for(c in 1:thisNumCohorts){
    thisVar<-paste(thisName,c,"_StructN",sep="")
    temp<-ncvar_get(ThisNC.nc,thisVar)
    predSNarray[,c]<-apply(temp,3,nonZeroMean)[(burnin+1):(nts+burnin)]
  }
  #get prey SN
  preyName<-str_trim(groupsDF$Name[groupsDF$Code==thisPrey],side="both"); preyNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisPrey]
  preySNarray<-array(NA,dim=c(nts,preyNumCohorts))
  for(c in 1:preyNumCohorts){
    thisVar<-paste(preyName,c,"_StructN",sep="")
    temp<-ncvar_get(ThisNC.nc,thisVar)
    preySNarray[,c]<-apply(temp,3,nonZeroMean)[(burnin+1):(nts+burnin)]
  }
  maxPredLimit<-max(predSNarray*thisKUP)
  SNmax<-max(max(preySNarray,na.rm=TRUE),max(predSNarray,na.rm=TRUE)); SNmin<-min(min(preySNarray,na.rm=TRUE),min(predSNarray,na.rm=TRUE))
  par(mfrow=c(5,2),mar=c(1,4,1,1),oma=c(4,1,1,1))
  for(c in 1:thisNumCohorts){
    thisMax<-max(predSNarray,na.rm=TRUE)*thisKUP
    thisMax<-max(preySNarray,na.rm=TRUE)
    plot(x=seq(1,nts),y=rep(thisMax,nts),type="n",ylim=c(0,thisMax),xlab="",xaxt="n",ylab="SN")
    if(c %in% c(thisNumCohorts,(thisNumCohorts-1))){
      axis(at=xLabsAt,labels=xLabs,side=1)
    }
    thisLower<-thisKLP*predSNarray[,c]; thisUpper<-thisKUP*predSNarray[,c]
    thisY<-c(thisLower,rev(thisUpper)); thisX<-c(seq(1,nts),seq(nts,1))
    polygon(x=thisX,y=thisY,col=cohortCols_trans[c],border=NA)
    for(k in 1:preyNumCohorts){
      points(x=seq(1,nts),y=preySNarray[,k],col=cohortCols[k],lwd=2,lty=2,type="l")
    }
  }
  par(las=1)
  mtext(paste("Predator: ",thisPred,", prey: ",thisPrey,sep=""),side=3,adj=0.5,outer=TRUE,line=-0.5)
 
  
   # }
  
  # dev.off()

# }
# 
# for(thisPred in c("HOK", "HAK", "PIN", "PFL", "MAC")){
#   plotGapeLimits(thisPred=thisPred, thisPrey = "RFI")
# }
# 
# 
ThisIC.NC<-nc_open(paste(this_path,"..\\CRAM_input_short.nc", sep=""))
# testIC<-ncvar_get(ThisIC.NC, thisTracer)
# 
# thisTracer<-paste(preyName, 1,"_Nums", sep=""); temp<-ncvar_get(ThisNC.nc, thisTracer)
# test<-apply(temp[,,1], 2, sum)

##################
## spatial overlap
thisPred<-"HOK"; thisPrey<-"RFI"
predName<-str_trim(groupsDF$Name[groupsDF$Code==thisPred], side="both"); preyName<-str_trim(groupsDF$Name[groupsDF$Code==thisPrey])
predTracer<-paste(predName,"_N", sep=""); preyTracer<-paste(preyName,"_N", sep="")
predData<-ncvar_get(ThisIC.NC, predTracer)[,,1]; preyData<-ncvar_get(ThisIC.NC,preyTracer)[,,1]

# #read in shape file
shapeFile<-paste(DIR$'Base',"ATLANTISmodels\\inputs\\bgm\\CHAT30_LL",sep="")
sdata<-read.shapefile(shapeFile)
shape<-formatShape(shapeFile=shapeFile)
ns<-length(shape)
SpDF <- SpatialPolygonsDataFrame(shape,data.frame( z=1:ns, row.names=paste("P",seq(1,(ns)),sep="")))
labels<-seq(1,(ns))
pdf("test.pdf")
plot(shape)
LABELpOS<-polygonsLabel(shape, labels = labels, cex=.1,doPlot=FALSE)
labeldf<-data.frame(cbind("x"=LABELpOS[1:ns],"y"=LABELpOS[(ns+1):(2*ns)]))
dev.off()

plotColour <- colorRampPalette(colors=c("snow", "midnightblue"))(101)
# functions needed
getPropColor <- function(x, log=FALSE, base=10){
  thisCol<-"white"
  if(!is.na(x)){
    thisIndex <- round((x-thisMin)/(thisMax-thisMin),2)*100 +1
    if(log==TRUE){
      thisIndex <- round((log(x, base=base)/log(thisMax, base=base)),2)*100 +1
    }
    if(length(thisIndex)>0){
      if(thisIndex > 101){
        thisIndex <- 101
      }
      if(sum(x, na.rm=TRUE)>0 & thisIndex>0){
        thisCol<-plotColour[thisIndex]
      }
    }
  }
  return(thisCol)
}

thisMax<-max(predData); thisMin <-0
colByYearPred<-apply(predData,c(1,2), getPropColor)

pdf(paste(plotPath,"SpatialDistrib",thisPred,"_", thisPrey,".pdf", sep=""), height=5, width=5)
par(mar=c(2,3,1,1), las=1, mfrow=c(2,1))
plot(shape)
map('nzHires',add=TRUE,col="black",lwd=2)
map.axes()
for(b in 1:dim(labeldf)[1]){
  thisCol<-colByYearPred[1,b]
  polygon(sdata$shp$shp[[b]]$points,col=thisCol)
}  
mtext(thisPred,side=3,adj=0)

thisMax<-max(preyData)
colByYearPrey<-apply(preyData,c(1,2), getPropColor)

par(mar=c(2,3,1,1), las=1)
plot(shape)
map('nzHires',add=TRUE,col="black",lwd=2)
map.axes()
for(b in 1:dim(labeldf)[1]){
  thisCol<-colByYearPrey[1,b]
  polygon(sdata$shp$shp[[b]]$points,col=thisCol)
}  
mtext(thisPrey,side=3,adj=0)
dev.off()

