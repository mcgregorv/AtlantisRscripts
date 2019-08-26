# nFutureYears<-50
basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

plotPath<-paste(DIR$'Figures',"ClimateChange\\", sep="")

dynBoxes<-2:25
nboxes<-30; nlayers<-6
load(paste(basePath,"ClimateChange\\depthByCell", sep=""))## brings in depyByCell


## read in the input .nc file for temperature
thisTracer<-"temperature"
tempNC.nc<- nc_open(paste(basePath, "inputs\\ROMS\\Chatham30_tempAll.nc", sep=""))
tempInputs<-ncvar_get(tempNC.nc, "temperature")
inputYears <- sort(rep(1996:2004, 365*2));  uniqueInputYears <-unique(round(inputYears)); nROMSyears<-length(uniqueInputYears)
monthDays<-c(31,28,31,30,31,30,31,31,30,31,30,31)
monthYears <- c(rep(1, monthDays[1]), rep(2, monthDays[2]), rep(3, monthDays[3]), rep(4, monthDays[4]), rep(5, monthDays[5]), rep(6, monthDays[6])
                , rep(7, monthDays[7]), rep(8, monthDays[8]), rep(9, monthDays[9]), rep(10, monthDays[10]), rep(11, monthDays[11]), rep(12, monthDays[12]))
inputMonths <- rep(sort(rep(monthYears, 2)), 9)


## GET MEAN temp for each layer, box, timestep(within year) over the ROMS years
# first need to define using 4th dimension for year (to seperate from the 12 hour timesteps)
tempInputsByYears<-array(NA, dim=c(nlayers, nboxes, 730, nROMSyears)) # 730 is timesteps within year; nROMSyears is number of ROMS years
for(y in 1:nROMSyears){
  thisYear <- unique(inputYears)[y]
  thisYearIndex <- round(inputYears)==thisYear
  tempInputsByYears[,,,y]<- tempInputs[,,thisYearIndex]
}

# get mean for each cell and timestep over the ROMS years
meanInputTemperatures <- apply(tempInputsByYears, c(1,2,3), mean, na.rm=TRUE)
SSTInputs <- apply(tempInputsByYears, c(2,3,4), max, na.rm=TRUE) # for each box, timestep, year, get the SST (as the max temp for this cell)
meanSSTs<-apply(SSTInputs, c(1,2), mean, na.rm=TRUE) ## the mean SST over all cells from the ROMS years

meanSSTByTS<-apply(SSTInputs, 2, mean, na.rm=TRUE)

## the 1999 ROMS year is consistently warmer than the other years, so use this to learn from
## how do the warmer temps look at depths and does this change by box?
baseYear <-1999; baseYearIndex <- round(inputYears)==baseYear
baseTemps <- tempInputs[,,baseYearIndex]
baseSST <- apply(baseTemps, c(2,3), max, na.rm=TRUE)
BaseSSTByTS <- apply(baseSST, 2, mean, na.rm=TRUE)
baseRelSST <- baseSST / meanSSTs

plot(meanSSTByTS, type="l", col=myGrey)
points(BaseSSTByTS, type="l", col="midnightblue")
# for each depth, how does 1999 look compared to the average of all ROMS..?
Depths<-sort(unique(round(as.vector(depthByCell),-1)))
## restructure depths by cell to match the input temperatures (matches the output tracers as it is)
inputDepthsByCell <- 0 * baseTemps[,,1]
for(b in 1:nboxes){
  thisTemplate <-!is.na(baseTemps[1:(nlayers-1),b,1]); thisDepths<-round(depthByCell[1:(nlayers-1),b])[rev(thisTemplate)]; 
  inputDepthsByCell[1:(nlayers-1),b][thisTemplate]<-thisDepths; inputDepthsByCell[nlayers, b]<- max(thisDepths, na.rm=TRUE)+1
}

save(inputDepthsByCell, file=paste(basePath,"ClimateChange\\inputDepthsByCell", sep=""))

meanTempByDepth<-array(NA, dim=c(length(Depths), 730)); meanBaseTempByDepth<-meanTempByDepth
for(d in 1:length(Depths)){
  thisDepth<-Depths[d]
  thisBaseTemps<-data.frame(matrix(NA, ncol=730, nrow=0)); thisAllTemps<- thisBaseTemps
  for(b in 1:nboxes){
    thisDepthIndex<-round(inputDepthsByCell[,b])==thisDepth
    if(sum(thisDepthIndex, na.rm=TRUE)>0){
      thisBaseTemp <- baseTemps[thisDepthIndex, b,]
      thisAllTemp <- tempInputsByYears[thisDepthIndex, b, ,]
      thisAllTemp <- apply(thisAllTemp, length(dim(thisAllTemp))-1, mean, na.rm=TRUE)
      
      if(length(dim(thisBaseTemp))>1){
        thisBaseTemp <- apply(thisBaseTemp, 2, mean, na.rm=TRUE)
       }
      thisBaseTemps <- rbind(thisBaseTemps, thisBaseTemp)
      thisAllTemps <- rbind(thisAllTemps, thisAllTemp)
    }
  }
  meanBaseTempByDepth[d,]<- apply(thisBaseTemps, 2, mean, na.rm=TRUE)
  meanTempByDepth[d,] <- apply(thisAllTemps, 2, mean, na.rm=TRUE)
}
thisMax<-max(max(meanBaseTempByDepth, na.rm=TRUE), max(meanTempByDepth, na.rm=TRUE))
thisMin<-min(min(meanBaseTempByDepth, na.rm=TRUE), min(meanTempByDepth, na.rm=TRUE))

par(mfrow=c(2,2))
for(d in 1:length(Depths)){
  plot(meanTempByDepth[d,], type="l", ylim=c(thisMin,thisMax), col=myGrey)
  points(meanBaseTempByDepth[d,], type="l", col="midnightblue")
  mtext(paste(Depths[d]," m", sep=""), side=3, adj=0)
}

deltaSST <- mean( BaseSSTByTS- meanSSTByTS, na.rm=TRUE)
## create DF with diff in temps (base - means)

diffTempsDF <- melt(meanBaseTempByDepth - meanTempByDepth); colnames(diffTempsDF)<-c("Depth", "TS","Diff")
diffTempsDF$Depth <- Depths[diffTempsDF$Depth]
## fit glm
# diffTempsDF$Depth <- as.factor(diffTempsDF$Depth)
modelDepth <- glm(Diff ~ Depth, data=diffTempsDF)

Xdepth<-round(cbind(anova(modelDepth),"Rsqr" = ((modelDepth$null.deviance - anova(modelDepth)[,4])/modelDepth$null.deviance)),2)

testData<-data.frame("Depth"=Depths); testData$Depth<-as.factor(testData$Depth)
testPRed<-predict(modelDepth, newdata = testData)
plot(x=Depths, y=testPRed, type="h")

## use the depth model to predict the 1999 data at depth based on the meanRoms temps and plot with actual
par(mfrow=c(2,2), mar=c(4,4,2,2))
for(d in 1:length(Depths)){
  thisPredDiff <- testPRed[d]
  plot(meanTempByDepth[d,], type="l", ylim=c(thisMin,thisMax), col=myGrey)
  points(meanBaseTempByDepth[d,], type="l", col="midnightblue")
  points(meanTempByDepth[d,]+ thisPredDiff, type="l", col="red", lwd=2, lty=2)
  mtext(paste(Depths[d]," m", sep=""), side=3, adj=0)
}
colByDepth <- colorRampPalette(colors=c(myLightBlue, myBlue, "midnightblue"))(length(Depths))
thisMin<-min(meanBaseTempByDepth, na.rm=TRUE); thisMax<-max(meanBaseTempByDepth, na.rm=TRUE)
par(mfrow=c(1,1), mar=c(4,4,2,2))
plot(1,xlim=c(thisMin, thisMax), ylim=c(thisMin, thisMax), type="n", xlab="Observerd",ylab="Fitted")
for(d in 1:length(Depths)){
  thisPredDiff <- testPRed[d]
  points(x=meanBaseTempByDepth[d,], y=meanTempByDepth[d,]+ thisPredDiff, pch=20, col=colByDepth[d])
}
points(x=seq(1,20), y=seq(1,20), type="l", col="red",lty=2, lwd=2)

save(modelDepth, file=paste(basePath,"ClimateChange\\modelDepth", sep=""))

modelDayDepth <- glm(Diff ~ Depth + TS, data=diffTempsDF)
Xdaydepth<-round(cbind(anova(modelDayDepth),"Rsqr" = ((modelDayDepth$null.deviance - anova(modelDayDepth)[,4])/modelDayDepth$null.deviance)),2)

baseRelTemps <- baseTemps/meanInputTemperatures  
## turn to df so can analyse by depth and box
baseRelTempsDF<-data.frame(matrix(NA, ncol=6, nrow=0))
colnames(baseRelTempsDF)<-c("Depth","Box","TS","RelTemp", "RelDeltaTemp", "BoxDepth")
for(l in 1:(nlayers-1)){
  for(b in 1:nboxes){
    thisDepth<-round(depthByCell[l,b])
    thisSurfaceDepth <- min(depthByCell[,b])
    if(thisDepth != thisSurfaceDepth){
      thisData<-baseRelTemps[l,b,]; thisDeltaSST <- baseRelSST[b,]
      thisRelDelta <- thisDeltaSST/thisData
      boxDepth <- round(sum(depthByCell[,b]))
      thisDF<-data.frame(cbind("Depth"=rep(thisDepth, length(thisData)), "Box"=rep(b,length(thisData)), "TS"=seq(1,length(thisData))
                               , "RelTemp"=thisData, "RelDeltaTemp" = thisRelDelta, "BoxDepth"= boxDepth))
      baseRelTempsDF<-rbind(baseRelTempsDF, thisDF)
    }
  }
}
for(d in 1:length(unique(baseRelTempsDF$Depth))){
  
}


baseRelTempsDF$Box <- as.factor(baseRelTempsDF$Box)
baseRelTempsDF$Depth <- as.factor(baseRelTempsDF$Depth)
baseRelTempsDF$BoxDepth <- as.factor(baseRelTempsDF$BoxDepth)
## fit models
model1 <- glm(RelDeltaTemp ~ Depth, data=baseRelTempsDF)
model2 <- glm(RelDeltaTemp ~ BoxDepth, data=baseRelTempsDF)
model3 <- glm(RelDeltaTemp ~ Box + Depth, data=baseRelTempsDF)
model4 <- glm(RelDeltaTemp ~ Box + Box:Depth, data=baseRelTempsDF)


modelAICs<- c(AIC(model1), AIC(model2), AIC(model3), AIC(model4))
AICs_axis<-pretty(seq(0,max(abs(modelAICs)), length.out = 10))

X1<-round(cbind(anova(model1),"Rsqr" = ((model1$null.deviance - anova(model1)[,4])/model1$null.deviance)),2)
X2<-round(cbind(anova(model2),"Rsqr" = ((model2$null.deviance - anova(model2)[,4])/model2$null.deviance)),2)
X3<-round(cbind(anova(model3),"Rsqr" = ((model3$null.deviance - anova(model3)[,4])/model3$null.deviance)),2)
X4<-round(cbind(anova(model4),"Rsqr" = ((model4$null.deviance - anova(model4)[,4])/model4$null.deviance)),2)

modelR2s <- c(X1$Rsqr[dim(X1)[1]], X2$Rsqr[dim(X2)[1]], X3$Rsqr[dim(X3)[1]], X4$Rsqr[dim(X4)[1]])



byDepth<-tapply(baseRelTempsDF$RelDeltaTemp, baseRelTempsDF$Depth, mean, na.rm=TRUE)
byBox <- tapply(baseRelTempsDF$RelDeltaTemp, baseRelTempsDF$Box, mean, na.rm=TRUE)
byBoxSST <- apply(baseRelSST, 1, mean, na.rm=TRUE)

baseRelTemp <- byBoxSST/byBox

tempOverSST<- 0 * tempInputs
SSTinput<-tempInputs[(nlayers-1),,]
for(b in dynBoxes){
  xx<-tempInputs[,b,]
  index<-!is.na(xx[,1]); index[nlayers]<-FALSE;  xxx <- seq(1,(nlayers))[index]
  thisSST <- xx[max(xxx),]
  thisDepths<-depthByCell[,b]
  thisLayerIndex <- sort(nlayers-xxx)
  for(l in 1:length(xxx)){
    thisLI <- thisLayerIndex[l]
    tempOverSST[thisLI,b,] <- xx[l,] / thisSST
  }
}
meanByLayer<-apply(tempOverSST[,dynBoxes,], 1, nonZeroMean)
## take out the sediment and surface layers
minByts<-apply(tempInputs, 3, sum, na.rm=TRUE); baseTS<-seq(1,length(minByts))[minByts==min(minByts)]
baseTemperature <- tempInputs[,,baseTS]
baseSST <- apply(baseTemperature, 2, max, na.rm=TRUE)
relTemperature <- 0 * tempInputs
for(t in 1:dim(tempInputs)[3]){
  thisSST <- apply(tempInputs[,,t],2,max, na.rm=TRUE) / baseSST ## actually this delta SST
  test<- (tempInputs[,,t]/baseTemperature)
  thisSSTfill <- t(0*t(test) + (thisSST))
  relTemperature[,,t] <-thisSSTfill / test
}
colByLayer <- colorRampPalette(colors=c(myLightBlue,myBlue,"midnightblue"))(nlayers-1)

thisMin<-min(relTemperature, na.rm=TRUE); thisMax <- max(relTemperature, na.rm=TRUE)

pdf(paste(plotPath,"relativeSeaTemperatureByBoxLayer.pdf", sep=""))
par(mfrow=c(5,2), mar=c(4,4,2,1))
for(b in 1:nboxes){
  if(b %in% dynBoxes){
    plot(relTemperature[1,1,], type="n", ylim=c(thisMin, thisMax), xlab="Timestep", ylab="Relative temperature")
    thisDepths<- depthByCell[,b]
    thisNl<-length(unique(thisDepths))-1 # one is sediment, not water column
    for(l in 1:thisNl){
      lIndex<-thisNl-l+1
      points(relTemperature[l,b,], type="l", col=colByLayer[lIndex], lwd=2)
    }
    points(relTemperature[nlayers, b, ], type="l", col=myGreen, lwd=2, lty=2)
    mtext(paste("Box: ",b), side=3, adj=0)
  }
}
dev.off()

# set up df with year, month, box, depth, reltemp
inputYears <- sort(rep(1996:2004, 365*2)); 
monthDays<-c(31,28,31,30,31,30,31,31,30,31,30,31)
monthYears <- c(rep(1, monthDays[1]), rep(2, monthDays[2]), rep(3, monthDays[3]), rep(4, monthDays[4]), rep(5, monthDays[5]), rep(6, monthDays[6])
                , rep(7, monthDays[7]), rep(8, monthDays[8]), rep(9, monthDays[9]), rep(10, monthDays[10]), rep(11, monthDays[11]), rep(12, monthDays[12]))
inputMonths <- rep(sort(rep(monthYears, 2)), 9)
fitData<-data.frame(matrix(NA, ncol=5, nrow=0))
colnames(fitData)<- c("Year","Month","Box","Depth","relTemperature")
for(b in 1:nboxes){
  if(b %in% dynBoxes){
    thisDepths<- depthByCell[,b]
    thisNl<-length(unique(thisDepths))-1 # one is sediment, not water column
    
    for(l in 1:thisNl){
      thisDepth <- unique(thisDepths)[l]
      thisData<-data.frame(cbind("Year"=inputYears, "Month"=inputMonths, "Box"=rep(b, length(inputMonths)), 
                                 "Depth"=rep(thisDepth, length(inputMonths)), "relTemperature"=relTemperature[l,b,] ))
      fitData <- rbind(fitData, thisData)
    }
  }
}
fitData$Month <- as.factor(fitData$Month)
model1 <- glm( relTemperature ~ Month + Depth, data=fitData); modelDescrip <-"MonthPlusDepth"
model2 <- glm( relTemperature ~ Depth, data=fitData); modelDescrip <-"Depth"
model3 <- glm( relTemperature ~ Month, data=fitData); modelDescrip<-"Month"
model4 <- glm( relTemperature ~ Month + Depth + Month:Depth , data=fitData); modelDescrip<-"MonthIntDepth"

modelAICs<- c(AIC(model1), AIC(model2), AIC(model3), AIC(model4))
AICs_axis<-pretty(seq(0,max(abs(modelAICs)), length.out = 10))

X1<-round(cbind(anova(model1),"Rsqr" = ((model1$null.deviance - anova(model1)[,4])/model1$null.deviance)),2)
X2<-round(cbind(anova(model2),"Rsqr" = ((model2$null.deviance - anova(model2)[,4])/model2$null.deviance)),2)
X3<-round(cbind(anova(model3),"Rsqr" = ((model3$null.deviance - anova(model3)[,4])/model3$null.deviance)),2)
X4<-round(cbind(anova(model4),"Rsqr" = ((model4$null.deviance - anova(model4)[,4])/model4$null.deviance)),2)

modelR2s <- c(X1$Rsqr[dim(X1)[1]], X2$Rsqr[dim(X2)[1]], X3$Rsqr[dim(X3)[1]], X4$Rsqr[dim(X4)[1]])
plotModelLabels<-modelLabels; plotModelLabels[4]<-"Month + Depth +\nMonth:Depth"

pdf(paste(plotPath, "TemperatureSSTmodelCompare.pdf", sep=""), height=5, width=7)
par(lend=1, yaxs="i", las=1, mar=c(12,5,2,7))
plot(x=seq(1,length(modelR2s))-0.1, y=modelR2s, type="h", lwd=10, col=myBlue, yaxt="n", xaxt="n",  xlab="",ylab="", ylim=c(0, max(modelR2s)*1.1), xlim=c(0.5,4.5))
axis(at=seq(0.1,1,by=0.1), labels=seq(0.1,1,by=0.1), side=2, col.axis=myBlue, col=myBlue, cex.axis=thisCex)
par(las=0)
mtext(expression(R^2), col=myBlue, side=2, adj=0.5, font=2, line=3, cex=thisCex)
mtext("|AIC|", col=myOrange, side=4, adj=0.5, font=2, line=5, cex=thisCex)
par(las=2)
axis(at=seq(1,length(plotModelLabels)), labels = plotModelLabels, side=1, cex.axis=thisCex)
par(new=TRUE)
plot(x=seq(1,length(modelR2s))+0.1, y=abs(modelAICs), type="h", lwd=10, col=myOrange, xaxt="n", yaxt="n",ylab="", xlab="", xlim=c(0.5,4.5), ylim=c(0,max(abs(modelAICs))*1.05))
axis(at=AICs_axis, labels=AICs_axis, side=4, col.axis=myOrange, col=myOrange, cex.axis=thisCex)
dev.off()


for(m in 1:3){
  model <- list(model1, model2, model3)[[m]]
  modelDescrip <- c("MonthPlusDepth", "Depth", "Month")[m]
## use the relationship between SST and base temp to turn it back into actual temperature - for the predicted as well
# we can assume we 'know' SST and deltaSST, and need to estimate temperature at depth
### NEXT NEXT
## turn this into loop so can see each cell - could do pred vs obs by box too
## check model diags, and check if depth or box should be added - or perhaps quadrant or long or lat
  fittedRelTemp <- predict(model, data=fitData)
  pdf(paste(plotPath, "CompareTemperatureFitted_MODEL",modelDescrip,".pdf", sep=""))
  par(mfrow=c(5,2), mar=c(4,4,2,1))
  for(b in dynBoxes){
    thisSST <- apply(tempInputs[,b,], 2, max, na.rm=TRUE)
    thisDeltaSST <- thisSST / baseSST[b]
    thisDepths <- unique(depthByCell[,b])
    this_nlayers<-length(thisDepths)-2 # -1 takes off sediment, -2 also takes off surface layer
    if(this_nlayers>0){
      for(l in 1:this_nlayers){
        thisDepth<-thisDepths[l]
        index<-fitData$Box == b & round(fitData$Depth)==round(thisDepth)
        thisData <- fitData[index,]; predData<-fittedRelTemp[index]
        thisTemperatures <- tempInputs[l,b,]
        predDeltaTemperatures <- thisDeltaSST/predData 
        predTemperatures <- predDeltaTemperatures * baseTemperature[l,b]
        thisYmax<-max(max(predTemperatures, na.rm=TRUE), max(thisTemperatures, na.rm=TRUE))
        thisYmin<-min(min(predTemperatures, na.rm=TRUE), min(thisTemperatures, na.rm=TRUE))
        plot(thisTemperatures, pch=20, ylim=c(thisYmin, thisYmax), ylab="Temperatue", xlab="Timestep")
        points(predTemperatures, type="l", col="red", lwd=1.2); mtext(paste("Box: ",b,", depth: ", round(thisDepth)," m", sep=""), side=3, adj=0)
      }
    }
  }
  dev.off()
}
# layerIndex<-round(thisDepths)==300
# thisInputTemperature <- tempInput[,b,]
# plot(thisData$relTemperature,pch=20, lwd=2)
# points(predData, type="l", col="red", lwd=3)

## save the model out
save(model4,file=paste(basePath,"ClimateChange\\TemperatureSSTmodel", sep=""))
save(model2,file=paste(basePath,"ClimateChange\\TemperatureSSTmodelDepth", sep=""))

write.csv(fitDataExtended, file=paste(basePath, "ClimateChange\\TemperatureSSTmodelData.csv", sep=""), row.names = FALSE)

#########################################################################
calcCV<-function(x){
  thisMean<-mean(x,na.rm=TRUE); thisVar<-var(x,na.rm=TRUE)
  thisCV<-sqrt(thisVar)/thisMean
  return(thisCV)
}
calcPearsonResids <- function(O,F,CV){
  pr<-(O-F)/(F*CV)
  return(pr)
}

## add actual temperatures, and fitted temperatures for each model to fitData
modelDescriptions <- c("MonthPlusDepth", "Depth", "Month", "MonthIntDepth")
modelLabels <- c("Month + Depth", "Depth","Month","Month + Depth + Month:Depth")
fitDataExtended <- fitData
for(m in 1:4){
  model <- list(model1, model2, model3, model4)[[m]]
  modelDescrip <- modelDescriptions[m]
  fitDataExtended[,c(paste("M",m, sep=""))]<-NA
  if(m==1){fitDataExtended[,c("obsTemperature")]<-NA}
  ## use the relationship between SST and base temp to turn it back into actual temperature - for the predicted as well
  # we can assume we 'know' SST and deltaSST, and need to estimate temperature at depth
  ### NEXT NEXT
  ## turn this into loop so can see each cell - could do pred vs obs by box too
  ## check model diags, and check if depth or box should be added - or perhaps quadrant or long or lat
  fittedRelTemp <- predict(model, data=fitData)
  for(b in dynBoxes){
    thisSST <- apply(tempInputs[,b,], 2, max, na.rm=TRUE)
    thisDeltaSST <- thisSST / baseSST[b]
    thisDepths <- unique(depthByCell[,b])
    this_nlayers<-length(thisDepths)-2 # -1 takes off sediment, -2 also takes off surface layer
    if(this_nlayers>0){
      for(l in 1:this_nlayers){
        thisDepth<-thisDepths[l]
        index<-fitData$Box == b & round(fitData$Depth)==round(thisDepth)
        thisData <- fitData[index,]; predData<-fittedRelTemp[index]
        thisTemperatures <- tempInputs[l,b,]
        predDeltaTemperatures <- thisDeltaSST/predData 
        predTemperatures <- predDeltaTemperatures * baseTemperature[l,b]
        if(m==1){ ## these don't change by model, so only need to populate once
          fitDataExtended[index,c("obsTemperature")] <- thisTemperatures
        }
        fitDataExtended[index, c(paste("M",m, sep=""))] <- predTemperatures
      }
    }
  }
}

obsCV<-calcCV(x=fitDataExtended$obsTemperature)

for(m in 1:4){
  fitDataExtended[,c(paste("resid_M",m,sep=""))] <- mapply(FUN=calcPearsonResids, "O"=fitDataExtended$obsTemperature, "F"=fitDataExtended[,c(paste("M",m, sep=""))], CV=obsCV)
}

thisYmax <- max(fitDataExtended[,grep("resid",colnames(fitDataExtended))], na.rm=TRUE)
thisYmin <- min(fitDataExtended[,grep("resid", colnames(fitDataExtended))], na.rm=TRUE)

# par(mfrow=c(4,1), mar=c(3,4,2,1))
for(m in 1:4){
  pdf(paste(plotPath, "TemperatureResidualsByMonthMODEL",m,".pdf", sep=""), height=3.5, width=5)
  par(mar=c(4,4.5,2,1), las=1)
  plot(1:12, ylim=c(thisYmin, thisYmax), type="n", ylab="Pearson's residuals", xlab="Month", cex.lab=thisCex, cex.axis=thisCex)
  abline(h=0, col="red",lty=2, lwd=2)
  mtext(modelLabels[m], side=3, adj=0, cex=thisCex)
  for(mnth in 1:12){
     boxplot( fitDataExtended[fitDataExtended$Month==mnth, paste("resid_M",m,sep="")], add=TRUE, at=mnth, pch=20, cex=0.6, yaxt="n")
  }
  dev.off()
}

residSumByModel <- apply(fitDataExtended[,grep("resid", colnames(fitDataExtended))],2, sum, na.rm=TRUE)

depthsTable <-table(round(fitDataExtended$Depth[!is.na(fitDataExtended$obsTemperature)]))
fittedDepths <- as.double(names(depthsTable)); ndepths<-length(fittedDepths)

# par(mfrow=c(4,1), mar=c(3,4,2,1))
for(m in 1:4){
  pdf(paste(plotPath, "TemperatureResidualsByDepthMODEL",m,".pdf", sep=""), height=3.5, width=5)
  par(mar=c(4,4.5,2,1), las=1)
  plot(fittedDepths, type="n", ylim=c(thisYmin, thisYmax), ylab="Pearson's residuals", xlab="Depth (m)", xaxt="n", cex.lab=thisCex, cex.axis=thisCex)
  abline(h=0, col="red",lty=2, lwd=2)
  axis(at=1:ndepths, labels=fittedDepths, side=1, cex.lab=thisCex, cex.axis=thisCex)
  mtext(modelLabels[m], side=3, adj=0, cex=thisCex)
  for(d in 1:ndepths){
    thisDepth <- fittedDepths[d]
    boxplot(fitDataExtended[round(fitDataExtended$Depth)==thisDepth, paste("resid_M",m,sep="")], add=TRUE, at=d, yaxt="n", pch=20, cex=0.6)
  }
  dev.off()
}


# par(mfrow=c(4,1), mar=c(3,4,2,1))
for(m in 1:4){
  pdf(paste(plotPath, "TemperatureResidualsByBoxMODEL",m,".pdf", sep=""), height=3.5, width=5)
  par(mar=c(4,4.5,2,1), las=1)
  plot(dynBoxes, type="n", ylim=c(thisYmin, thisYmax), ylab="Pearson's residuals", xlab="Box", xlim=c(min(dynBoxes), max(dynBoxes)), cex.lab=thisCex, cex.axis=thisCex)
  abline(h=0, col="red",lty=2, lwd=2)
  mtext(modelLabels[m], side=3, adj=0)
  for(b in dynBoxes){
    boxplot(fitDataExtended[(fitDataExtended$Box)==b, paste("resid_M",m,sep="")], add=TRUE, at=b, yaxt="n", pch=20, cex=0.6)
  }
  dev.off()
}


residSumByMonth <- tapply(fitDataExtended[,grep("resid",colnames(fitDataExtended))], fitDataExtended$Month, sum, na.rm=TRUE)

# > meanByLayer
# [1] 0.4608811 0.6019230 0.7011310 0.8301316 1.0000000       NaN
par(lend=1)
plot(meanByLayer[1:(nlayers-1)], type="h", lwd=5, xlab="Depth layer", ylab="Proportion of SST")

fitData<-data.frame(cbind("depth"=rep(as.vector(depthByCell[1:(nlayers-2),dynBoxes]), dim(tempOverSST)[3]), 
                          "relTemp"=as.vector(tempOverSST[1:(nlayers-2),dynBoxes,])))
thismod<-lm(relTemp ~ depth, data=fitData)
thisInt<-summary(thismod )$coefficient[1]; thisSlope <- summary(thismod )$coefficient[2]
# want to rescale the non zero mean, then re-aporiton the new mean by the previous distrution wrt depth layer and box
scaleYearIndex <- allYears %in% scaleYears

## compare delta SST to delta temperature in all cells and at depth
# take one timestep as base, then the others can be difference relative to this ts



newMean <- diatomMean[scaleYearIndex] * (1- testPPtimeseriesScalars/100)

sumByScaleYear<-apply(diatomData[-nlayers, dynBoxes, ], 3, sum)[scaleYearIndex]
scaledDiatomData<-0 * diatomData[,,scaleYearIndex]
for(y in 1:dim(scaledDiatomData)[3]){
  temp<-diatomData[,,scaleYearIndex][,,y]/diatomMean[scaleYearIndex][y]
  scaledDiatomData[,,y]<- temp * newMean[y]
}

test<-apply(scaledDiatomData[-nlayers, dynBoxes,] * baseVol[-nlayers, dynBoxes,scaleYearIndex], 3, sum)*mg_2_tonne *X_CN
testBase<-diatomBiomass[scaleYearIndex]

plot(testBase,type="l", ylim=c(0,max(testBase)))
points(test,type="l", col="red")

## now add the scaled PP data to the force file
thisForceFile<-paste(basePath, "\\inputs\\PPforcing_",thisTracer,".txt", sep="")
thisHeader<-paste("netcdf Chatham_PPtest {
                  dimensions:
                  t = UNLIMITED ; // (4 currently)
                  b = 30 ;
                  z = 6 ;
                  variables:
                  double t(t) ;
                  t:units = \"seconds since ",startYear,"-01-01 00:00:00 +10\" ;
                  t:dt = 43200. ;
                  double Diatom_N(t, b, z) ;
                  Diatom_N:_FillValue = -999.;
                  Diatom_N:missing_value = -999.;
                  Diatom_N:valid_min = 0.;
                  Diatom_N:valid_max = 300.;
                  Diatom_N:units = \"mg N m-3\" ;
                  
                  // global attributes:
                  :title = \"trivial\" ;
                  :geometry = \"CHAT30_aea.bgm\" ;
                  :parameters = \"\" ;
                  data:\n\n",sep="")
cat(thisHeader, file=thisForceFile, append=FALSE)

#need the timesteps in seconds
scaleSeconds <- (scaleYears - startYear) * 60*60*24*365

cat("t = ",paste(scaleSeconds, collapse=", ", sep=""), ";\n", file=thisForceFile, append=TRUE)

## now add the forced data
cat(paste(thisTracer, " = \n", sep=""), file=thisForceFile, append=TRUE)
for(y in 1:nFutureYears){
  thisData<-scaledDiatomData[,,y]
  for(b in 1:dim(thisData)[2]){
    xx<-signif(thisData[,b],4)
    ## check no negatives
    xx[xx<0]<-1e-8
    cat(paste(xx,collapse=", "), file=thisForceFile, append=TRUE); 
    if(b==dim(thisData)[2] & y==nFutureYears){
      cat(";\n", file=thisForceFile, append=TRUE)
    } else{
      cat(",\n", file=thisForceFile, append=TRUE)
    }
    
    
  }
}
cat("}", file=thisForceFile, append=TRUE)


# 
# testNC.nc<-nc_open(paste(basePath,"TBGBMM\\output\\output.nc", sep=""))
# 
# testDino<-ncvar_get(testNC.nc, "Diatom_N")
# testVol<-ncvar_get(testNC.nc, "volume")
# testNlayers<-dim(testVol)[1]; testDynBoxes<-2:dim(testVol)[2]
# 
# testDinoTonnes <- apply(testDino[-testNlayers, testDynBoxes,] * testVol[-testNlayers, testDynBoxes,], 3, sum) * mg_2_tonne * X_CN
# 
# modelArea<-sum(testVol[nlayers,,1])



