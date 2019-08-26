source(paste(DIR$'General functions',"make_grid.r",sep=""))

#read in outputDietCheck.txt
options(stringsAsFactors = FALSE)

WC_COL<-myBlue
SD_COL<-myRed
EP_COL<-myGreen

habitat_layers<-c("WC","SED","EPIBENTHIC")
habitat_cols<-c(WC_COL,SD_COL,EP_COL)

this_run<-"base"
this_out<-paste("FISH",sep="")

# 
basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
outPath<-paste(basePath,this_run,"\\","output",this_out,"\\",sep="")
burnin<-35
burnin<-1

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1900
fishingStartYear<-1900
modelStartYear<-1900

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
ntsteps<-dim(thisVol)[3]

groupsDF<-read.csv(paste(basePath, "CRAM_groups.csv", sep="")); 
ng<-dim(groupsDF)[1]

nts<-dim(thisVol)[3]-burnin #number of timesteps
cat(paste(nts,"\n"))

xLabsTemp<-seq(0,(nts*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0

plotPath<-paste(basePath,"\\Figures\\",this_run,"\\DIET_check\\", this_out,sep="")

diet_check<-read.csv(paste(outPath,"outputDietCheck.txt",sep=""),skip=0,sep=" ",header=TRUE)
predators<-unique(diet_check$Predator); npreds<-length(predators)
all_preys<-colnames(diet_check)[6:dim(diet_check)[2]]
thisColorRamp<-colorRampPalette(colors=c(myOrange,myRed,myPurple,"midnightblue",myBlue,myAqua,myGreen))

colorByPredProp<-colorRampPalette(colors=c(myLightAqua, myBlue,"midnightblue"))(101)
getPropColor<-function(x){
  thisCol<-"white"
  if(!is.na(x) & x>0){
    y<-round(x/maxProp,2)*100 + 1
    thisCol<-colorByPredProp[y]
    if(x>maxProp){thisCol=colorByPredProp[101]}
  }
  return(thisCol)
}
plotGrid<-function(x,y){
  thisX<-c(x-0.5,x-0.5,x+0.5,x+0.5); thisY<-c(y-0.5,y+0.5,y+0.5,y-0.5)
  thisCol<-colByYearPred[x,y]
  if(length(thisCol)>0){
    polygon(x=thisX,y=thisY,col=thisCol,border=NA)
  }
  return(NULL)
}

p=47 # for RFI
p=27 #GSH
p=25 #ELP
p=26 # ETB
p<-grep("HOK",all_preys)
for(p in 1:length(all_preys)){
  thisPrey<-all_preys[p]
  thisIndex<-diet_check[,c(thisPrey)]>0
  thisPredData<-diet_check[thisIndex,c(1:3,grep(thisPrey,colnames(diet_check)))]
  myPredators<-sort(unique(thisPredData$Predator)); nmyPreds<-length(myPredators)
  if(length(myPredators)>1){
    timeIndex<-thisPredData$Time/365 > burnin
    sumByPredatorYear<-tapply(thisPredData[timeIndex,c(thisPrey)],thisPredData[timeIndex,c("Time","Predator")],sum,na.rm=TRUE)
    #turn into props within each year
    yearSums<-rowSums(sumByPredatorYear, na.rm=TRUE)
    propsByYear<-sumByPredatorYear/yearSums
    propsByYear<-sumByPredatorYear
    maxProp<-max(propsByYear, na.rm=TRUE)
    test<-sort(as.vector(propsByYear)); 
    minProp<-5e-3 # only keep those with at least 1% at least once
    minProp<-0
    testMin<-apply(propsByYear, 2, max, na.rm=TRUE)
    minIndex<-testMin>minProp
    colByYearPred<-apply(propsByYear[,minIndex],c(1,2), getPropColor)
    tempDF<-data.frame(cbind("x"=rep(seq(1,dim(propsByYear[,minIndex])[1]),dim(propsByYear[,minIndex])[2]),"y"=sort(rep(seq(1,dim(propsByYear[,minIndex])[2]),dim(propsByYear[,minIndex])[1]))))
    
    # pdf(paste(plotPath, thisPrey,"predation.pdf", sep=""), height=10, width=6)
    jpeg(paste(plotPath, thisPrey,"predation.jpg", sep=""), height=700, width=600)
    par(mar=c(4,4,1.5,1), mfrow=c(2,1), fig=c(0,1,0.25,1))
    plot(x=seq(1,dim(propsByYear[,minIndex])[1]),y=rep(dim(propsByYear[,minIndex])[2],dim(propsByYear[,minIndex])[1]),type="n",xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,dim(propsByYear[,minIndex])[2]))
    axis(at=xLabsAt,labels=xLabs,side=1)
    axis(at=seq(1,dim(propsByYear[,minIndex])[2]),labels=colnames(propsByYear)[minIndex],side=2,las=1)
    temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
    box()
    thisPreyName<-str_trim(groupsDF$Name[groupsDF$Code==thisPrey]); thisVar<-paste(thisPreyName,"_N", sep="")
    xx<-ncvar_get(ThisNC.nc, thisVar); 
    if(length(dim(xx))==3){
      xxx<-apply(xx*thisVol,3,sum)* X_CN * mg_2_tonne; thisBiomass<-xxx[burnin:(nts+burnin-1)]
    } else{
      xxx<-apply(xx*thisVol[nlayers,,],2,sum)* X_CN * mg_2_tonne; thisBiomass<-xxx[burnin:(nts+burnin-1)]
    }
    mtext(paste("Predators of ", thisPrey,sep=""), side=3, adj=0)
    par(mar=c(4,4,1.5,1), mfrow=c(2,1), fig=c(0,1,0,0.3))
    par(new=TRUE)
    plot(thisBiomass,type="l",col=myBlue,lwd=3, xaxt="n",xlab="", ylim=c(0,max(thisBiomass)), ylab="Biomass (tonnes)")
    
    axis(at=xLabsAt,labels=xLabs,side=1)
    abline(v=seq(0,nts,by=10),col=myGrey_trans)
    dev.off()
    
    # 
    # #do legend
    # pdf(paste(plotPath, "RFIpredationLEGEND.pdf", sep=""), height=3, width=3)
    # legend<-pretty(seq(minProp,maxProp, length.out = 5)); legend[1]<-signif(minProp,2); legend[length(legend)]<-signif(maxProp,2); legendCols<-unlist(lapply(legend, getPropColor))
    # par(mar=c(0,0,0,0))
    # makeBlankPlot()
    # legend(legend=legend, col=legendCols, pch=15, x="center", bty="n", pt.cex=1.5, title="Proportion of predation\n within each year")
    # dev.off()
    
  } 
}
  
realisedDiets<-diet_check
realisedDiets$year_ts<-realisedDiets$Time/365
yearIndex<-realisedDiets$year_ts %in% dietTimeSteps
yearRelDiets<-realisedDiets[yearIndex, ]

for(p in 1:npreds){
  thisPred<-predators[p]
  temp<-realisedDiets[realisedDiets$Predator==thisPred,]
  test<-temp[,6:(dim(temp)[2]-1)]
  posPreyIndex<-colSums(test, na.rm=TRUE)>0
  thisDiet<-cbind("year"=temp$year_ts, test[,posPreyIndex])
  thisPreys<-colnames(test)[posPreyIndex]; nprey<-length(thisPreys)
  #get proportions of each prey
  totalPrey<-sum(thisDiet[,-1], na.rm=TRUE); preySums<-apply(thisDiet[,-1],2,sum,na.rm=TRUE); preyProps<-preySums/totalPrey
  if(totalPrey>0){
    timeSums<-array(NA, dim=c(length(unique(thisDiet$year)),(dim(thisDiet)[2]-1)))
    for(i in 1:length(thisPreys)){
      xx<-thisDiet[,c(1,(i+1))]; xxx<-tapply(xx[,2], xx[,1], sum, na.rm=TRUE)
      timeSums[,i]<-xxx
    }
    minProp<-0
    testMin<-apply(timeSums, 2, max, na.rm=TRUE)
    minIndex<-testMin>minProp; maxProp<-max(timeSums)
    colByYearPred<-apply(timeSums[,minIndex],c(1,2), getPropColor)
    tempDF<-data.frame(cbind("x"=rep(seq(1,dim(timeSums[,minIndex])[1]),dim(timeSums[,minIndex])[2]),"y"=sort(rep(seq(1,dim(timeSums[,minIndex])[2]),dim(timeSums[,minIndex])[1]))))
    
    jpeg(paste(plotPath, thisPred,"prey.jpg", sep=""), height=700, width=600)
    par(mar=c(4,4,1.5,1), mfrow=c(2,1), fig=c(0,1,0.25,1))
    plot(x=seq(1,dim(timeSums[,minIndex])[1]),y=rep(dim(timeSums[,minIndex])[2],dim(timeSums[,minIndex])[1]),type="n",xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,dim(timeSums[,minIndex])[2]))
    axis(at=xLabsAt,labels=xLabs,side=1)
    axis(at=seq(1,dim(timeSums[,minIndex])[2]),labels=thisPreys[minIndex],side=2,las=1)
    temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
    box()
    thisPreyName<-str_trim(groupsDF$Name[groupsDF$Code==thisPred]); thisVar<-paste(thisPreyName,"_N", sep="")
    xx<-ncvar_get(ThisNC.nc, thisVar); 
    if(length(dim(xx))==3){
      xxx<-apply(xx*thisVol,3,sum)* X_CN * mg_2_tonne; thisBiomass<-xxx[burnin:(nts+burnin-1)]
    } else{
      xxx<-apply(xx*thisVol[nlayers,,],2,sum)* X_CN * mg_2_tonne; thisBiomass<-xxx[burnin:(nts+burnin-1)]
    }
    mtext(paste("Prey of ", thisPred,sep=""), side=3, adj=0)
    par(mar=c(4,4,1.5,1), mfrow=c(2,1), fig=c(0,1,0,0.3))
    par(new=TRUE)
    plot(thisBiomass,type="l",col=myBlue,lwd=3, xaxt="n",xlab="", ylim=c(0,max(thisBiomass)), ylab="Biomass (tonnes)")
    
    axis(at=xLabsAt,labels=xLabs,side=1)
    abline(v=seq(0,nts,by=10),col=myGrey_trans)
    dev.off()
  }
} 


