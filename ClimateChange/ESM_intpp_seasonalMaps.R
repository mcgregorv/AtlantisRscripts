plotPath<-paste(DIR$'Figures',"ClimateChange\\PP",sep="")

thisESM<-"gggfdlext"

ESMdesc<-c("GFDL", "CESM")
# thisESM<-"ggcesmext"
thisRCP<-"85"
thisRCP<-"45"

RCPs<-c("45","85")

thisTestStat<-"mean"
thisTestStat <-"min"
thisTestStat<-"max"
thisTestStat<-"median"

thisColRamp<-colorRampPalette(colors=c(myLightAqua, myAqua, myBlue, "midnightblue"))(101)
thisColRamp<-colorRampPalette(colors=c( "midnightblue", myBlue, myAqua, myGreen, myYellow, myGold, myOrange,"red"))(101)

getColor<-function(x){
  thisCol<-"white"
  y<-x; 
   if(!is.na(x)){
    if(x>thisMax){y<-thisMax}
    if(x<thisMin){y<-thisMin}
    z<-round((y-thisMin)/ (thisMax-thisMin),2) *100 +1
    thisCol<-thisColRamp[z]
  }
  return(thisCol)
}
## log scale alt
getColor<-function(x){
  thisCol<-"white"
  y<-x; 
  if(!is.na(x)){
    if(x>thisMax){y<-thisMax}
    if(x<thisMin){y<-thisMin}
    z<-round((log(y)-log(thisMin))/ (log(thisMax)-log(thisMin)),2) *100 +1
    thisCol<-thisColRamp[z]
  }
  return(thisCol)
}

yearBlocks<-NULL; yearBlocks[[1]]<-1986:2005; yearBlocks[[2]]<-2036:2055; yearBlocks[[3]]<-2081:2100
yearDesc<-c("Present","Mid-century","End-century")

yearBlocks<-NULL; yearBlocks[[1]]<-1866:1885; yearBlocks[[2]]<-1986:2005; yearBlocks[[3]]<-2036:2055; yearBlocks[[4]]<-2081:2100; 
yearDesc<-c("Past","Present","Mid-century","End-century")

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
  
  ## limit intPP data to only points within our region
  # define quadrant index 
  Q1lat<-ppLat[1,]> -43.5 & ppLat[1,] < -42; Q1lon <- ppLon[1,] <180 & ppLon[1,] > 174
  Q2lat<-ppLat[1,] > -43.5 & ppLat[1,] < -42; Q2lon <- ppLon[1,] >180
  Q3lat<-ppLat[1,]<= -43.5 & ppLat[1,] > -45.5; Q3lon <- ppLon[1,] <180 & ppLon[1,] > 172
  Q4lat<-ppLat[1,] <= -43.5 & ppLat[1,] > -45.5; Q4lon <- ppLon[1,] >180  & ppLon[1,] < 186
  
  QlatIndexes <- list(Q1lat, Q2lat, Q3lat, Q4lat); QlonIndexes <- list(Q1lon, Q2lon, Q3lon, Q4lon)
  
  yy <- Q1lat + Q2lat + Q3lat + Q4lat; xx <- Q1lon +  Q2lon +  Q3lon +  Q4lon
  QallLat <- yy>0; QallLon <- xx>0
  
  
  test<-tapply(apply(intPP[QallLon, QallLat,],3,thisTestStat, na.rm=TRUE), dataYears, thisTestStat, na.rm=TRUE)
  ppBTime_hist<-apply(intPP[QallLon, QallLat,],3,thisTestStat, na.rm=TRUE)
  
  futureYears<-2006:2100; nfyears<-length(futureYears)
  
  datafMonths<-rep(seq(1,12), nfyears)
  datafYears<-sort(rep(futureYears,12))
  for(R in 1:2){
    thisRCP<-RCPs[R]
      
    futurePPdata<-nc_open(paste(dataPath, thisESM,  "_rcp",thisRCP,"_2006_2100_intpp.nc", sep=""))
    future_intPP<-ncvar_get(futurePPdata, "data")
    ppfLat<-ncvar_get(ppData, "bound_latitude"); ppfLon<-ncvar_get(ppData, "bound_longitude")
    ## limit intPP data to only points within our region
    # define quadrant index 
    Q1flat<-ppfLat[1,]> -43.5 & ppfLat[1,] < -42; Q1flon <- ppfLon[1,] <180 & ppfLon[1,] > 174
    Q2flat<-ppfLat[1,] > -43.5 & ppfLat[1,] < -42; Q2flon <- ppfLon[1,] >180
    Q3flat<-ppfLat[1,]<= -43.5 & ppfLat[1,] > -45.5; Q3flon <- ppfLon[1,] <180 & ppfLon[1,] > 172
    Q4flat<-ppfLat[1,] <= -43.5 & ppfLat[1,] > -45.5; Q4flon <- ppfLon[1,] >180  & ppfLon[1,] < 186
    QflatIndexes <- list(Q1flat, Q2flat, Q3flat, Q4flat); QflonIndexes <- list(Q1flon, Q2flon, Q3flon, Q4flon)
    yy <- Q1flat + Q2flat + Q3flat + Q4flat; xx <- Q1flon +  Q2flon +  Q3flon +  Q4flon
    QallfLat <- yy>0; QallfLon <- xx>0
    
    testf<-tapply(apply(future_intPP[QallfLon, QallfLat,],3,thisTestStat, na.rm=TRUE), datafYears, thisTestStat, na.rm=TRUE)
    ppBTime_future<-apply(future_intPP[QallfLon, QallfLat,],3,thisTestStat, na.rm=TRUE)
    
    xRange=c(min(dataYears), max(datafYears))
    yRange<-c(min(c(ppBTime_future, ppBTime_hist), na.rm=TRUE), max(c(ppBTime_future, ppBTime_hist), na.rm=TRUE))
    
    ## first set up the range of values for the colours
    thisMax<-max(max(intPP[QallLon, QallLat,], na.rm=TRUE), max(future_intPP[QallfLon, QallfLat,], na.rm=TRUE)); thisMin<-min(min(intPP[QallLon, QallLat,], na.rm=TRUE), min(future_intPP[QallfLon, QallfLat,], na.rm=TRUE))
    
    # cut down the data for the year range interested in for this plot
    thisYears<-1986:2005; 
    thisYears<-2036:2055
    thisYears<-2081:2100
    
    for(Y in 1:length(yearBlocks)){
      thisYears<-yearBlocks[[Y]]
      for(M in 1:12){
        # months2plot<-c(11,12)
        months2plot<-M
        if(sum(thisYears %in% dataYears)==length(thisYears)){
         yearIndex<-dataYears %in% thisYears 
         year_intPP<- intPP[,,yearIndex]; yearMonths<-dataMonths[yearIndex]
         yearLats<-ppLat; yearLons<-ppLon
        } else{
          yearIndex<-datafYears %in% thisYears 
          year_intPP<- future_intPP[,,yearIndex]; yearfMonths<-dataMonths[yearIndex]
          yearLats<-ppfLat; yearLons<-ppfLon
        }
        
        ncells <- dim(yearLats)[2] * dim(yearLons)[2]
        thisGrid <- data.frame(matrix(NA, ncol=3, nrow=ncells))
        colnames(thisGrid) <- c("Cell","Lat","Lon")
        thisGrid$Cell <- 1:ncells
        thisGrid$Lat <- rep(yearLats[1,], dim(yearLons)[2]); thisGrid$Lon <- sort(rep(yearLons[1,], dim(yearLats)[2]))
        
        png(paste(plotPath,"Map_ESM",E, "RCP",thisRCP,"Years",min(thisYears),"to",max(thisYears),months2plot,".png",sep=""))
        plot(shape,ylim=c(-48,-40))
        for(p in 1:ncells){
          thisLon<- thisGrid$Lon[p]; thisLat <- thisGrid$Lat[p]
          thisLon2<-min(thisGrid$Lon[thisGrid$Lon>thisLon]); thisLat2<-min(thisGrid$Lat[thisGrid$Lat>thisLat])
          thisX<- c(thisLon, thisLon2, thisLon2, thisLon); thisY<-c(thisLat, thisLat, thisLat2, thisLat2)
          thisLatIndex<-yearLats[1,]==thisLat; thisLonIndex<-yearLons[1,]==thisLon; thisPPs<-year_intPP[thisLonIndex, thisLatIndex,]
          thisAverage<-mean(thisPPs[yearMonths %in% months2plot]); thisColor<-getColor(thisAverage)
          polygon(x=thisX, y=thisY, col=thisColor, border=NA)
        }
        map('nzHires',add=TRUE,col="black",lwd=2)
        map.axes()
        mtext(paste(ESMdesc[E],", RCP", as.double(thisRCP)/10,"\n", min(thisYears),":", max(thisYears), ", Months ",min(months2plot),":",max(months2plot), sep=""), side=3, adj=0, cex=2)
        dev.off()
      }
    }
  }
}

## create latex insert
reportPath<-paste(DIR$'Reports',"(02)ClimateChange\\", sep="")
Mblocks<-NULL; Mblocks[[1]]<-1:5; Mblocks[[2]]<-6:10; Mblocks[[3]]<-11:12
Mblocks<-NULL; Mblocks[[1]]<-1:8; Mblocks[[2]]<-9:12;

for(E in 1:2){
  for(R in 1:2){
    for(M in 1:length(Mblocks)){
  thisMs<-Mblocks[[M]]
  texfile<-paste(reportPath,"ESMppFiguresInsert",E,R,M,".txt", sep="")
  cat("\\begin{figure}[H]
  	\\centering
      ", file=texfile,append=FALSE)
  for(m in thisMs){
    for(Y in 1:length(yearBlocks)){
        
      thisFigFile<-paste("PPMap_ESM",E,"RCP",RCPs[R],"Years",min(yearBlocks[[Y]]),"to", max(yearBlocks[[Y]]),m, ".png", sep="")
      cat(paste("\\includegraphics[width=3cm]{",thisFigFile,"}
                ", sep=""), file=texfile, append=TRUE)
    }
  }
  cat(paste("\\caption{PP",E,R,Y,".}\\label{PP",E,R,Y,"}
              ", sep=""), file=texfile, append=TRUE)
  cat("\\end{figure}", file=texfile, append=TRUE)
    }
  }
}
