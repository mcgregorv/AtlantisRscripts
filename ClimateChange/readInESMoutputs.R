thisESM<-"gggfdlext"
# thisESM<-"ggcesmext"
thisRCP<-"85"
thisRCP<-"45"

thisTestStat<-"mean"
thisTestStat <-"min"
thisTestStat<-"max"
thisTestStat<-"median"

for(E in 1:2){
  thisESM<-c("gggfdlext", "ggcesmext")[E]
  
  dataPath<-paste(DIR$'Base',"data\\ClimateModels\\", sep="")
  ppData<-nc_open(paste(dataPath, thisESM, "_1861_2005_intpp.nc", sep=""))
  
  names(ppData$var)
  ppLon<-ncvar_get(ppData, "bound_longitude")
  ppLat<-ncvar_get(ppData, "bound_latitude")
  intPP<-ncvar_get(ppData, "data")
  
  modelYears<-1861:2005; nyears<-length(modelYears)
  
  dataMonths<-rep(seq(1,12), nyears)
  dataYears<-sort(rep(modelYears,12))
  
  if(thisESM=="ggcesmext"){
    modelYears<-1860:2005; 
    nyears<-length(modelYears)
    
    dataMonths<-rep(seq(1,12), nyears)[-1]
    dataYears<-sort(rep(modelYears,12))[-1]
  }
  
  test<-tapply(apply(intPP,3,thisTestStat, na.rm=TRUE), dataYears, thisTestStat, na.rm=TRUE)
  ppBTime_hist<-apply(intPP,3,thisTestStat, na.rm=TRUE)
  
  futureYears<-2006:2100; nfyears<-length(futureYears)
  
  datafMonths<-rep(seq(1,12), nfyears)
  datafYears<-sort(rep(futureYears,12))
  
  futurePPdata<-nc_open(paste(dataPath, thisESM,  "_rcp",thisRCP,"_2006_2100_intpp.nc", sep=""))
  future_intPP<-ncvar_get(futurePPdata, "data")
  
  testf<-tapply(apply(future_intPP,3,thisTestStat, na.rm=TRUE), datafYears, thisTestStat, na.rm=TRUE)
  ppBTime_future<-apply(future_intPP,3,thisTestStat, na.rm=TRUE)
  
  xRange=c(min(dataYears), max(datafYears))
  yRange<-c(min(c(ppBTime_future, ppBTime_hist), na.rm=TRUE), max(c(ppBTime_future, ppBTime_hist), na.rm=TRUE))
  
  # par(mfrow=c(2,1))
  plot(x=dataYears, y=ppBTime_hist,pch=20, xlim=xRange, ylim=yRange)
  points(x=datafYears, y=ppBTime_future, pch=20, col=myOrange)
  points(x=modelYears, y=test, type="l", col=myAqua, lwd=3)
  points(x=futureYears, y=testf, type="l", col="red", lwd=3)
  mtext(paste(thisESM,", RCP", as.double(thisRCP)/10,": ", thisTestStat, sep=""), side=3, adj=0)
}

###########################
## spatially
E=1; thisRCP<-"85"
storeMonthAverages<-array(NA, dim=c(2,nyearSets, 12))
storeMonthQAverages<- array(NA, dim=c(2,nyearSets, 4, 12))
storeMonthQMins<- 0 * storeMonthQAverages; storeMonthQMaxs <- 0* storeMonthQAverages
storeALL <- NULL
storeModelYears<-NULL
storeEAST<-NULL; storeWEST <- NULL

for(E in 1:2){
  thisESM<-c("gggfdlext", "ggcesmext")[E]

  dataPath<-paste(DIR$'Base',"data\\ClimateModels\\", sep="")
  ppData<-nc_open(paste(dataPath, thisESM, "_1861_2005_intpp.nc", sep=""))
  
  names(ppData$var)
  ppLon<-ncvar_get(ppData, "bound_longitude")
  ppLat<-ncvar_get(ppData, "bound_latitude")
  intPP<-ncvar_get(ppData, "data")
  
  # So for model CESM, the months go from (2,1860) to (12,2005) and then {1,2006} to (12,2100).
  # For model GFDL the months are (1,1861) to (12,2005) and then (1,2006) to (12,2100).
  modelYears<-1861:2005; nyears<-length(modelYears)
  dataMonths<-rep(seq(1,12), nyears)
  dataYears<-sort(rep(modelYears,12))
  if(thisESM=="ggcesmext"){
    modelYears<-1860:2005; 
    nyears<-length(modelYears)
    
    dataMonths<-rep(seq(1,12), nyears)[-1]
    dataYears<-sort(rep(modelYears,12))[-1]
  }
  
  
  futureYears<-2006:2100; nfyears<-length(futureYears)
  datafMonths<-rep(seq(1,12), nfyears)
  datafYears<-sort(rep(futureYears,12))
  
  futurePPdata<-nc_open(paste(dataPath, thisESM,  "_rcp",thisRCP,"_2006_2100_intpp.nc", sep=""))
  future_intPP<-ncvar_get(futurePPdata, "data")
  ppfLon<-ncvar_get(futurePPdata, "bound_longitude")
  ppfLat<-ncvar_get(futurePPdata, "bound_latitude")

  ## summarise spatially the 1995:2005, 1945:55, 2090:2100 year-sets
  currentYears<-1995:2005; midCenturyYears<-2045:2055; endCenturyYears<-2090:2100
  yearSets<-list("CYs"=currentYears, "MCYs"=midCenturyYears, "ECYs"=endCenturyYears)
  nyearSets<-length(yearSets)
  
  yearSets<-list(); midYears<-seq(2000,2090, by=10); nyearSets<-length(midYears)
  for(i in 1:nyearSets){
    thisYearSet<-(midYears[i]-4) : (midYears[i]+5)
    yearSets[[i]]<-thisYearSet
  }
  
  # define quadrant index 
  Q1lat<-ppLat[1,]> -43.5 & ppLat[1,] < -42; Q1lon <- ppLon[1,] <180 & ppLon[1,] > 174
  Q2lat<-ppLat[1,] > -43.5 & ppLat[1,] < -42; Q2lon <- ppLon[1,] >180
  
  Q3lat<-ppLat[1,]<= -43.5 & ppLat[1,] > -45.5; Q3lon <- ppLon[1,] <180 & ppLon[1,] > 172
  Q4lat<-ppLat[1,] <= -43.5 & ppLat[1,] > -45.5; Q4lon <- ppLon[1,] >180  & ppLon[1,] < 186
  
  QlatIndexes <- list(Q1lat, Q2lat, Q3lat, Q4lat); QlonIndexes <- list(Q1lon, Q2lon, Q3lon, Q4lon)
  
  yy<-Q1lat + Q2lat + Q3lat + Q4lat;  QcombinedLatIndex<-yy>0
  xx <- Q1lon + Q2lon + Q3lon + Q4lon; QcombinedLonIndex <- xx>0
  
  xx <- Q1lon + Q3lon; Qwest<-xx>0
  xx <- Q2lon + Q4lon; Qeast<-xx>0
  
  for(ys in 1:nyearSets){
    thisYears<-yearSets[[ys]]; yearSet<-names(yearSets)[ys]
    if(sum(thisYears %in% dataYears)==length(thisYears)){
      thisIntPP<-intPP[,,dataYears %in% thisYears]
    } else{
      thisIntPP<-future_intPP[,,datafYears %in% thisYears]
    }
    # over these years, get the average by month, and the average by cell
    thisMonths<-rep(1:12, length(thisYears))
    monthAve<-tapply(apply(thisIntPP[QcombinedLonIndex, QcombinedLatIndex , ], 3, mean, na.rm=TRUE), thisMonths, mean, na.rm=TRUE)
    storeMonthAverages[E,ys,]<-monthAve
    ## spatially, split by quadrant
    for(q in 1:4){
      thisQdata<-thisIntPP[ QlonIndexes[[q]], QlatIndexes[[q]], ]
      monthQAve<-tapply(apply(thisQdata, 3, mean, na.rm=TRUE), thisMonths, mean, na.rm=TRUE)
      storeMonthQAverages[E,ys, q,]<- monthQAve
      
      monthQAve<-tapply(apply(thisQdata, 3, min, na.rm=TRUE), thisMonths, min, na.rm=TRUE)
      storeMonthQMins[E,ys, q,]<- monthQAve
      
      
      monthQAve<-tapply(apply(thisQdata, 3, max, na.rm=TRUE), thisMonths, max, na.rm=TRUE)
      storeMonthQMaxs[E,ys, q,]<- monthQAve
      
    }
    
  }
  ## populate the ALL array
  westPP <- intPP[Qwest,,]; futureWestPP<-future_intPP[Qwest,,]
  eastPP <- intPP[Qeast,,]; futureEastPP<-future_intPP[Qeast,,]
  
  
  intQPP <- intPP[QcombinedLonIndex, QcombinedLatIndex,]
  future_intQPP <- future_intPP[QcombinedLonIndex, QcombinedLatIndex,]
  ntotalYears<-nyears + nfyears;   storeModelYears[[E]]<-c(modelYears, futureYears)
  thisALL <- array( NA, dim=c(dim(intQPP)[1:2], ntotalYears, 12))
  thisEast <- array(NA, dim=c(dim(eastPP)[1:2], ntotalYears, 12)); thisWest <- array(NA, dim=c(dim(westPP)[1:2], ntotalYears, 12))
  for(y in 1:ntotalYears){
    if(y <= nyears){
      thisYear<-modelYears[y]
      thisData<-intQPP[,,dataYears ==thisYear]; thisMonths<-dataMonths[dataYears==thisYear]
      thisWestData <- westPP[,,dataYears==thisYear]; thisEastData <- eastPP[,,dataYears==thisYear]
    } else{
      thisYear <- futureYears[(y-nyears)]
      thisData<-future_intQPP[,,datafYears ==thisYear]; thisMonths <- datafMonths[datafYears==thisYear]
      thisWestData <- futureWestPP[,,datafYears==thisYear]; thisEastData <- futureEastPP[,,datafYears==thisYear]
      
    }
    for(M in 1:length(thisMonths)){
      thisM<-thisMonths[M]
      thisALL[,,y,thisM]<- thisData[,,M]
      thisWest[,,y,thisM]<- thisWestData[,,M]; thisEast[,,y,thisM]<- thisEastData[,,M]
    }
  }
  storeALL[[E]]<-thisALL
  storeWEST[[E]]<-thisWest
  storeEAST[[E]] <- thisEast  
}

calcCV<-function(x){
  thisMean<-mean(x, na.rm=TRUE); thisVar<-var(x, na.rm=TRUE)
  thisCV<- sqrt(thisVar)/ thisMean
  return(thisCV)
}
CVs1<-apply(storeALL[[1]], 3, calcCV)
CVs2<-apply(storeALL[[2]], 3, calcCV)

yRange<-c(min(c(CVs1, CVs2)), max(c(CVs1, CVs2)))
plot(CVs1, pch=20, col=myBlue, ylim=yRange)
points(CVs2 , col=myOrange, pch=20)

thisTestStat<-"median"

par(mfrow=c(2,2), mar=c(4,4,1,1))
for(M in 1:12){
  mean1 <- apply(storeALL[[1]][,,,M], 3, thisTestStat, na.rm=TRUE)
  mean2<-apply(storeALL[[2]][,,,M], 3, thisTestStat, na.rm=TRUE)
  yRange<-c(max(0,min(c(mean1, mean2), na.rm=TRUE)), max(c(mean1, mean2), na.rm=TRUE))
  plot(mean1, pch=20, col=myBlue, ylim=yRange)
  points(mean2, pch=20, col=myOrange)
  mtext(M, side=3, adj=0)

}

for(M in 1:12){
  mean1 <- apply(storeEAST[[1]][,,,M], 3, thisTestStat, na.rm=TRUE)
  mean2<-apply(storeEAST[[2]][,,,M], 3, thisTestStat, na.rm=TRUE)
  yRange<-c(max(0,min(c(mean1, mean2), na.rm=TRUE)), max(c(mean1, mean2), na.rm=TRUE))
  plot(mean1, pch=20, col=myBlue, ylim=yRange)
  points(mean2, pch=20, col=myOrange)
  mtext(M, side=3, adj=0); mtext("EAST", side=3, adj=1)
  
}

for(M in 1:12){
  mean1 <- apply(storeWEST[[1]][,,,M], 3, thisTestStat, na.rm=TRUE)
  mean2<-apply(storeWEST[[2]][,,,M], 3, thisTestStat, na.rm=TRUE)
  yRange<-c(max(0,min(c(mean1, mean2), na.rm=TRUE)), max(c(mean1, mean2), na.rm=TRUE))
  plot(mean1, pch=20, col=myBlue, ylim=yRange)
  points(mean2, pch=20, col=myOrange)
  mtext(M, side=3, adj=0); mtext("WEST", side=3, adj=1)
  
}
  ## CVs by month and year
CVs1<-apply(storeALL[[1]], 4, calcCV)
CVs2<-apply(storeALL[[2]], 4, calcCV)

plotYears<-2000:2050; nplotYears<-length(plotYears)
colByTime<-colorRampPalette(colors=c("red", myYellow,myGreen, myLightBlue, myBlue,"midnightblue"))(nplotYears)

yRange <- c(min(c(min(storeALL[[1]], na.rm=TRUE), min(storeALL[[2]], na.rm=TRUE))), 
            max(c(max(storeALL[[1]], na.rm=TRUE), max(storeALL[[2]], na.rm=TRUE))))
par(mfrow=c(2,1), mar=c(4,4,1,1))
for(E in 1:2){
  plot(x=1:12, y=rep(0,12), type="n", ylim=yRange, ylab="IntPP", xlab="Month")
  for(y in 1:nplotYears){
    thisY<-plotYears[y]; thisYindex<-storeModelYears[[E]]==thisY
    thisData<-apply(storeALL[[E]][,,thisYindex,], 3, mean, na.rm=TRUE)
    points(thisData, type="l", col=colByTime[y])
  }
}

colByYearSet<-colorRampPalette(colors=c(myLightBlue, myBlue,  "midnightblue"))(nyearSets)
E=1
xx<-(storeMonthQAverages[E,1,1,]-storeMonthQAverages[E,nyearSets,1,])/ storeMonthQAverages[E,1,1,]
yRange<- c(myRounding(min(xx),0.1), myRounding(max(xx),0.1))*2.5
par(mfrow=c(2,2), mar=c(4,4,1,1))
for(q in 1:4){
  plot((storeMonthQAverages[E,1,q,]-storeMonthQAverages[E,nyearSets,q,])/ storeMonthQAverages[E,1,q,], 
       xlab="Month", type="n", ylab="Proportional intPP decrease", ylim=yRange)
  for(ys in 1:nyearSets){
    points((storeMonthQAverages[E,1, q,]-storeMonthQAverages[E,ys, q,])/ storeMonthQAverages[E,1, q,], pch=20, col=colByYearSet[ys], cex=1.5)
    
  }
  abline(h=0, col="red", lty=2, lwd=1.5)
  
}

par(mfrow=c(2,2), mar=c(4,4,1,1))
yRange <- c(min(storeMonthQAverages)*0.9, max(storeMonthQAverages)*1.1)
for(q in 1:4){
  plot((storeMonthQAverages[E,1,q,]), 
       xlab="Month", type="n", ylab="IntPP", ylim=yRange)
  for(ys in 1:nyearSets){
    points((storeMonthQAverages[E,ys, q,]), pch=20, col=colByYearSet[ys], cex=1.5)
    
  }

}
yRange <- c(min(storeMonthQMins)*0.9, max(storeMonthQMaxs)*1.1)

E=2
for(q in 1:4){
  y2<-storeMonthQAverages[E,1,q,]; y1<-storeMonthQMins[E,1,q,]; y3<-storeMonthQMaxs[E,1,q,];
  x1<-seq(1, length(y2))
  plot((storeMonthQAverages[E,1,q,]), 
       xlab="Month", type="n", ylab="IntPP", ylim=yRange)
  polygon(x= c(x1, rev(x1)), y=c(y1, rev(y3)), col=myBlue_trans, border=NA)
  y2<-storeMonthQAverages[E,nyearSets,q,]; y1<-storeMonthQMins[E,nyearSets,q,]; y3<-storeMonthQMaxs[E,nyearSets,q,];
  polygon(x= c(x1, rev(x1)), y=c(y1, rev(y3)), col=myOrange_trans, border=NA)
}

colByQ<-c(myGold, myBlue,myRed, myGreen)

plot(shape,ylim=c(-48,-40))
map('nzHires',add=TRUE,col="black",lwd=2)
map.axes()
for(q in 1:4){
  points(x=rep(ppLon[1,QlonIndexes[[q]]], length(ppLat[1,QlatIndexes[[q]]])), 
             y=sort(rep(ppLat[1,QlatIndexes[[q]]], length(ppLon[1,QlonIndexes[[q]]]))), col=colByQ[q], pch=20)

}

yRange<-c(-0.1, 0.2)

plot((storeMonthAverages[1,]-storeMonthAverages[nyearSets,])/ storeMonthAverages[1,], 
     xlab="Month", type="n", ylab="Proportional intPP decrease", ylim=yRange)
for(ys in 1:nyearSets){
  points((storeMonthAverages[1,]-storeMonthAverages[ys,])/ storeMonthAverages[1,], pch=20, col=colByYearSet[ys], cex=1.5)
  
}
abline(h=0, col="red", lty=2, lwd=1.5)

## look at December- for example
M=12

