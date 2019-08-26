### DELETE DELETE DELETE #####


## similar to the original script, but takes each ROMS year in turn as the base year and fits the model to all relative temperatures

# nFutureYears<-50
basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

plotPath<-paste(DIR$'Figures',"ClimateChange\\", sep="")

dynBoxes<-2:25
nboxes<-30; nlayers<-6
load(paste(basePath,"ClimateChange\\inputDepthsByCell", sep=""))## brings in inputDepysByCell
depthsByCell<-inputDepthsByCell

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
# set up df with year, month, box, depth, reltemp
fitData<-data.frame(matrix(NA, ncol=5, nrow=0))
colnames(fitData)<- c("Year","Month","Box","Depth","relTemperature")
for(y in 1:nROMSyears){
  thisYear<-inputYears[y]; 
  thisTempInput <- tempInputsByYears[,,,y]
  ## what is the mean difference between these SSTs and the mean SSTs?
  thisSSTs<-apply(thisTempInput, c(2,3), max, na.rm=TRUE)
  thisRelSST <- thisSSTs/meanSSTs
  ## want to comapre to overal mean difference in SST
  deltaSST <- mean(thisRelSST)
  thisRelTemps <- thisTempInput / meanInputTemperatures
  thisRatioTemp <- deltaSST / thisRelTemps
  
  # thisRatioTemp <- 0 * thisRelTemps
  # for(b in 1:nboxes){
  #   for(t in 1:730){
  #     thisRatioTemp[,b,t]<- thisRelSST[b,t]  / thisRelTemps[,b,t]
  #   }
  # }
  # populate df to fit
  for(b in 1:nboxes){
    if(b %in% dynBoxes){
      thisDepths<- depthByCell[,b]
      thisNl<-length(unique(thisDepths))-1 # one is sediment, not water column
      
      for(l in 1:thisNl){
        thisDepth <- unique(thisDepths)[l]
        thisData<-data.frame(cbind("Year"=inputYears, "Month"=inputMonths, "Box"=rep(b, length(inputMonths)), 
                                   "Depth"=rep(thisDepth, length(inputMonths)), "relTemperature"=thisRatioTemp[l,b,] ))
        fitData <- rbind(fitData, thisData)
      }
    }
  }

}




fitData$Month <- as.factor(fitData$Month)
model1 <- glm( relTemperature ~ Month + Depth, data=fitData); modelDescrip <-"MonthPlusDepth"
model2 <- glm( relTemperature ~ Depth, data=fitData); modelDescrip <-"Depth"
model3 <- glm( relTemperature ~ Month, data=fitData); modelDescrip<-"Month"
model4 <- glm( relTemperature ~ Month + Depth + Month:Depth , data=fitData); modelDescrip<-"MonthAndIntDepth"
model5 <- glm( relTemperature ~ Month +  Month:Depth , data=fitData); modelDescrip<-"MonthIntDepth"


modelAICs<- c(AIC(model1), AIC(model2), AIC(model3), AIC(model4), AIC(model5))
AICs_axis<-pretty(seq(0,max(abs(modelAICs)), length.out = 10))

X1<-round(cbind(anova(model1),"Rsqr" = ((model1$null.deviance - anova(model1)[,4])/model1$null.deviance)),2)
X2<-round(cbind(anova(model2),"Rsqr" = ((model2$null.deviance - anova(model2)[,4])/model2$null.deviance)),2)
X3<-round(cbind(anova(model3),"Rsqr" = ((model3$null.deviance - anova(model3)[,4])/model3$null.deviance)),2)
X4<-round(cbind(anova(model4),"Rsqr" = ((model4$null.deviance - anova(model4)[,4])/model4$null.deviance)),2)
X5<-round(cbind(anova(model5),"Rsqr" = ((model5$null.deviance - anova(model5)[,4])/model5$null.deviance)),2)

modelR2s <- c(X1$Rsqr[dim(X1)[1]], X2$Rsqr[dim(X2)[1]], X3$Rsqr[dim(X3)[1]], X4$Rsqr[dim(X4)[1]], X5$Rsqr[dim(X5)[1]])
modelDescriptions <- c("MonthPlusDepth", "Depth", "Month", "MonthAndIntDepth", "MonthIntDepth")
modelLabels <- c("Month + Depth", "Depth","Month","Month + Depth + Month:Depth","Month + Month:Depth")

plotModelLabels<-modelLabels; plotModelLabels[4]<-"Month + Depth +\nMonth:Depth"; plotModelLabels[5]<-"Month +\nMonth:Depth"


pdf(paste(plotPath, "TemperatureSSTmodelCompareFromSamples.pdf", sep=""), height=5, width=7)
par(lend=1, yaxs="i", las=1, mar=c(12,5,2,7))
plot(x=seq(1,length(modelR2s))-0.1, y=modelR2s, type="h", lwd=10, col=myBlue, yaxt="n", xaxt="n",  xlab="",ylab="", ylim=c(0, max(modelR2s)*1.1), xlim=c(0.5,5.5))
axis(at=seq(0.1,1,by=0.1), labels=seq(0.1,1,by=0.1), side=2, col.axis=myBlue, col=myBlue, cex.axis=thisCex)
par(las=0)
mtext(expression(R^2), col=myBlue, side=2, adj=0.5, font=2, line=3, cex=thisCex)
mtext("|AIC|", col=myOrange, side=4, adj=0.5, font=2, line=5, cex=thisCex)
par(las=2)
axis(at=seq(1,length(plotModelLabels)), labels = plotModelLabels, side=1, cex.axis=thisCex)
par(new=TRUE)
plot(x=seq(1,length(modelR2s))+0.1, y=abs(modelAICs), type="h", lwd=10, col=myOrange, xaxt="n", yaxt="n",ylab="", xlab="", xlim=c(0.5,5.5), ylim=c(0,max(abs(modelAICs))*1.05))
axis(at=AICs_axis, labels=AICs_axis, side=4, col.axis=myOrange, col=myOrange, cex.axis=thisCex)
dev.off()


for(m in 1:3){
  model <- list(model1, model3, model5)[[m]]
  modelDescrip <- c("MonthPlusDepth",  "Month", "MonthIntDepth")[m]
  ## use the relationship between SST and base temp to turn it back into actual temperature - for the predicted as well
  # we can assume we 'know' SST and deltaSST, and need to estimate temperature at depth
  ### NEXT NEXT
  ## turn this into loop so can see each cell - could do pred vs obs by box too
  ## check model diags, and check if depth or box should be added - or perhaps quadrant or long or lat
  fittedRelTemp <- predict(model, data=fitData)
  pdf(paste(plotPath, "CompareTemperatureFitted_MODEL",modelDescrip,"FromSamples.pdf", sep=""))
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
save(model4,file=paste(basePath,"ClimateChange\\TemperatureSSTmodelFromSamples", sep=""))
save(model2,file=paste(basePath,"ClimateChange\\TemperatureSSTmodelDepthFromSamples", sep=""))

write.csv(fitDataExtended, file=paste(basePath, "ClimateChange\\TemperatureSSTmodelDataFromSamples.csv", sep=""), row.names = FALSE)

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
# modelDescriptions <- c("MonthPlusDepth", "Depth", "Month", "MonthIntDepth")
# modelLabels <- c("Month + Depth", "Depth","Month","Month + Depth + Month:Depth")
fitDataExtended <- fitData
thisModelIndex<-c(1,3,5)
model=model5; modelDescrip<-modelDescriptions[5]  
fitDataExtended[,c(paste("M",m, sep=""))]<-NA
fitDataExtended[,c("obsTemperature")]<-NA
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

for(m in 1:3){
  fitDataExtended[,c(paste("resid_M",m,sep=""))] <- mapply(FUN=calcPearsonResids, "O"=fitDataExtended$obsTemperature, "F"=fitDataExtended[,c(paste("M",m, sep=""))], CV=obsCV)
}

thisYmax <- max(fitDataExtended[,grep("resid",colnames(fitDataExtended))], na.rm=TRUE)
thisYmin <- min(fitDataExtended[,grep("resid", colnames(fitDataExtended))], na.rm=TRUE)

# par(mfrow=c(4,1), mar=c(3,4,2,1))
for(m in 1:3){
  pdf(paste(plotPath, "TemperatureResidualsByMonthMODEL",m,"FromSamples.pdf", sep=""), height=3.5, width=5)
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
  pdf(paste(plotPath, "TemperatureResidualsByDepthMODEL",m,"FromSamples.pdf", sep=""), height=3.5, width=5)
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
  pdf(paste(plotPath, "TemperatureResidualsByBoxMODEL",m,"FromSamples.pdf", sep=""), height=3.5, width=5)
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
# 
# 
# # select some timesteps as random to use as the base temperature
# sampleTS <- sample(1:dim(tempInputs)[3], size=1000, replace=FALSE)
# ## check on month coverage
# sampledMonths <- inputMonths[sampleTS]
# hist(sampledMonths)
# 
# for(T in 1:length(uniqueInputYears)){
#   baseTS<-sampleTS[T]
#   baseTemperature <- tempInputs[,,baseTS]
#   baseSST <- apply(baseTemperature, 2, max, na.rm=TRUE)
#   relTemperature <- 0 * tempInputs
#   for(t in 1:dim(tempInputs)[3]){
#     thisSST <- apply(tempInputs[,,t],2,max, na.rm=TRUE) / baseSST ## actually this delta SST
#     test<- (tempInputs[,,t]/baseTemperature)
#     thisSSTfill <- t(0*t(test) + (thisSST))
#     relTemperature[,,t] <-thisSSTfill / test
#   }
#   for(b in 1:nboxes){
#     if(b %in% dynBoxes){
#       thisDepths<- depthByCell[,b]
#       thisNl<-length(unique(thisDepths))-1 # one is sediment, not water column
#       
#       for(l in 1:thisNl){
#         thisDepth <- unique(thisDepths)[l]
#         thisData<-data.frame(cbind("Year"=inputYears, "Month"=inputMonths, "Box"=rep(b, length(inputMonths)), 
#                                    "Depth"=rep(thisDepth, length(inputMonths)), "relTemperature"=relTemperature[l,b,] ))
#         fitData <- rbind(fitData, thisData)
#       }
#     }
#   }
# }
# 
