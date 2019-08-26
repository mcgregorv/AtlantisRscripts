
getPredVars <- function(this_df){
  temp<- this_df$Factor
  temp2 <- sort(unique(unlist(lapply(temp, FUN=function(x){unlist(str_split(x,":"))}))))
  return(temp2)
}
absColorRamp <- colorRampPalette(colors=c(myLightAqua, myAqua, myBlue,"midnightblue"))(101)

getVarValues <- function(thisVar, thisData){
  thisTab <- table(thisData[,c(thisVar)])
  thisValues<-signif(as.double(names(thisTab)),2)
  if(length(thisValues)>10 | length(unique(thisValues))!=length(thisValues)){
    thisValues <- pretty(as.double(thisValues), n=10)
  }
  if(sum(thisValues, na.rm=TRUE)==0){
    thisValues <- names(thisTab)
  }
  return(thisValues)
}
#calc cv for each timestep - split fished from not-fished runs
calc_cv <- function(x){
  thisVar <- var(as.double(x), na.rm=TRUE); thisMean <- mean(as.double(x), na.rm=TRUE)
  thisCV <- (sqrt(thisVar))/thisMean
  return(thisCV)
}
# only keep values within fitted data ranges
getPlausiblePairs <- function(var1, var2, thisData){
  # var1 <- thisIntVars[1]; var2 <- thisIntVars[2]
  values1 <- unlist(thisPredVarValues[grep(var1, thisPredVars)]); values2 <- unlist(thisPredVarValues[grep(var2, thisPredVars)])
  thisPairs <- data.frame(cbind("v1"=sort(rep(values1, length(values2))), "v2"= rep(values2, length(values1))))
  thisPairs <- data.frame(matrix(NA, ncol=length(values2), nrow=length(values1))); rownames(thisPairs)<-values1; colnames(thisPairs)<- values2
  for(i in 1:length(values1)){
    thisLims <- values1[i]*c(0.5,1.5)
    thisIndex <- thisData[, c(var1)]<= max(thisLims) & thisData[, c(var1)]>= min(thisLims) 
    var2range <- c(min(thisData[thisIndex, c(var2)]), max(thisData[thisIndex, c(var2)]))
    thisPairIndex <- colnames(thisPairs)>=var2range[1] & colnames(thisPairs)<=var2range[2]
    thisPairs[i,thisPairIndex]<-1
  }
  return(thisPairs)
  
}
testPlausiblePair <- function(value1, value2, var1, var2, thisData){
  thisPairs <- getPlausiblePairs(var1, var2, thisData)
  test<-thisPairs[as.double(rownames(thisPairs))==value1, as.double(colnames(thisPairs))==value2]
  return(test)
}



plotIntEffect_noChaosRun <- function(thisIntVar, this_model, thisMin=0, inputMax=1, updateMax=TRUE, thisData=thisData){
  lab_cex=1.5
  # get predictor variable values
  thisPredVars <- unique(unlist(str_split(names(this_model$coefficients[-1]),":"))); 
  # fix any that have ChaosAlt as these are a factor
  thisPredVars[grep("^ChaosAlt", thisPredVars)]<-"ChaosAlt"
  thisPredVars <- sort(unique(thisPredVars))
  npvars <- length(thisPredVars)
  thisPredVarLengths <- rep(NA, npvars)
  thisPredVarValues <- NULL
  for(v in 1:npvars){
    thisVar <- thisPredVars[v]
       tempVarValues<-getVarValues(thisVar, thisData = this_model$data)
       if(is.factor(thisData[,thisVar])){
         tempVarValues <- sort(unique(thisData[,thisVar]))
       }
 
    ## fix the odd things - don't bother with zero biomass, lifespan, Linf, 
    # and keystonennes and responsiveness should just have unique values
    if(thisVar %in% c("B0", "Lifespan", "Linf", "NumL1cons", "PropByTopPrey","TL")){
      tempVarValues <- tempVarValues[tempVarValues!=0]
      # check not outside data limits
      dataMax <- max(this_model$data[,c(thisVar)], na.rm=TRUE); dataMin <- min(this_model$data[,c(thisVar)], na.rm=TRUE)
      tempVarValues <- tempVarValues[tempVarValues<=dataMax & tempVarValues>=dataMin]
      # } else if(thisVar %in% c("Keystone", "Informance")){
    } else if(thisVar %in% c("Informance", "Keystone","Response")){
      thisUnique <- sort(unique(this_model$data[,c(thisVar)]))
      tempVarValues <- thisUnique
    }
    # if(thisVar=="ChaosRun" & grep("ChaosRun", names(thisPredVarLengths)))
    thisPredVarValues[[v]] <- tempVarValues
    thisPredVarLengths[v]<-length(thisPredVarValues[[v]])
  }
  
  thisIntVars <- (unlist(str_split(thisIntVar,":")))
  thisValues1 <- as.double(thisPredVarValues[[grep(thisIntVars[1], thisPredVars)[1]]])
  thisValues2 <- as.double(thisPredVarValues[[grep(thisIntVars[2], thisPredVars)[1]]])
  if(is.factor(this_model$data[,c(thisIntVars[1])])){
    thisValues1<-thisPredVarValues[[grep(thisIntVars[1], thisPredVars)]]
  }
  if(is.factor(this_model$data[,c(thisIntVars[2])])){
    thisValues2<-unlist(thisPredVarValues[grep(thisIntVars[2], thisPredVars)])
  }
  test_df <- data.frame(matrix(0, ncol=dim(this_model$data)[2], nrow=(length(thisValues1) * length(thisValues2))))
  colnames(test_df)<- colnames(this_model$data)
  for(s in 1:npvars){
    fillVar <- thisPredVars[s]
    if(length(grep(fillVar, colnames(thisData)))>0){
      if( is.factor(thisData[,fillVar])==FALSE){
        temp<- boxplot(this_model$data[,c(fillVar)], plot=FALSE)$stats
        fillValueMed <- temp[3]; 
        test_df[,c(fillVar)]<- fillValueMed; 
      } else{
        test_df[,c(fillVar)]<-names(sort(table(thisData[,fillVar]),decreasing = TRUE)[1])
        test_df[,c(fillVar)] <- factor(test_df[,c(fillVar)], levels=sort(unique(thisData[,c(fillVar)])))
      }
    }
  }
  test_df[,c(thisIntVars[1])] <- sort(rep(thisValues1, length(thisValues2)));
  test_df[,c(thisIntVars[2])] <-  rep(thisValues2, length(thisValues1))
  if(is.factor(this_model$data[,c(thisIntVars[1])])){
    test_df[,c(thisIntVars[1])] <- factor(test_df[,c(thisIntVars[1])], levels=levels(this_model$data[,c(thisIntVars[1])]))
  }
  if(is.factor(this_model$data[,c(thisIntVars[2])])){
    test_df[,c(thisIntVars[2])] <- factor(test_df[,c(thisIntVars[2])], levels=levels(this_model$data[,c(thisIntVars[2])]))
  }
  
  test_df$Predicted <- ((predict.glm(this_model, newdata = test_df)))^3
  if(updateMax==TRUE){
    thisMax <<- max(test_df$Predicted, na.rm=TRUE); 
    # thisMin <- min(test_df$Predicted, na.rm=TRUE)
  } else{
    thisMax <<- inputMax
  }
  cat(paste("-- ",thisMax, thisMin, sep=", "))
  plotColour <<- unlist(lapply(test_df$Predicted, getColor, log=FALSE, thisMax=thisMax))
  # tempDF<-data.frame(cbind("x"=test_df[,c(thisIntVars[1])]),"y"=test_df[,c(thisIntVars[2])])
  x1<-sort(rep(seq(1,length(thisValues1)), length(thisValues2))); y1 <- rep(seq(1,length(thisValues2)), length(thisValues1))
  testDF <- cbind(x1,y1)

  tempDF<<-data.frame(testDF); colnames(tempDF)<-c("x", "y")
  # 
  #   # pdf(paste(plotPath,"KemptonsQ_SSRsens.pdf",sep=""),height=4,width=10)
  plot(x=tempDF$x,y=tempDF$y,type="n",xlab="",ylab="",yaxt="n", xaxt="n", xlim=c(min(tempDF$x)-0.5, max(tempDF$x)+0.5))
  temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
  par(las=1)
  xx <- thisValues1; xxx <- seq(1, length(thisValues1))
  if(length(xx)>6){
    aIndex<-seq(1,length(xx), by=ceiling(length(xx)/5))
    axis1Labs  <- xx[aIndex]; axis1At <- xxx[aIndex]
  }else{
    axis1Labs  <- xx; axis1At <- xxx
  }
  xx <- thisValues2; xxx <- seq(1, length(thisValues2))
  if(length(xx)>6){
    aIndex<-seq(1,length(xx), by=ceiling(length(xx)/5))
    axis2Labs  <- xx[aIndex]; axis2At <- xxx[aIndex]
  } else{
    axis2Labs  <- xx; axis2At <- xxx
  }
  axis(at=axis1At, labels=axis1Labs, side=1)
  axis(at=axis2At, labels=axis2Labs, side=2)
  par(las=0)
  thisLabel1 <- varDesc$Decr[varDesc$Variable == thisIntVars[1]]
  if(thisLabel1==""){
    thisLabel1 <- expression(B[0]~"('000 tonnes)")
  }
  thisLabel2 <- varDesc$Decr[varDesc$Variable == thisIntVars[2]]
  if(thisLabel2==""){
    thisLabel2 <- expression(B[0]~"('000 tonnes)")
  }
  mtext(thisLabel1, side=1, adj=0.5, line=2.5, cex=lab_cex); mtext(thisLabel2, side=2, adj=0.5, line=4, cex=lab_cex)
  box()
  if(updateMax==TRUE){
    mtext(paste("Max. CV: ", 100*round(thisMax, 2), "%", sep=""), side=3, adj=1, cex=lab_cex)
    # rm(thisMax)
  }
}


plotIntEffect <- function(thisIntVar, this_model, thisMin=0, inputMax=1, updateMax=TRUE, thisData=thisData){
  # get predictor variable values
  thisPredVars <- unique(unlist(str_split(names(this_model$coefficients[-1]),":"))); 
  # check for 'chaosA, B,..' and change to just 'chaos'
  x <- grep("Chaos", thisPredVars)
  if(length(x)>0){
    thisPredVars[x]<-"ChaosRun"
  }
  thisPredVars <- sort(unique(thisPredVars))
  npvars <- length(thisPredVars)
  thisPredVarLengths <- rep(NA, npvars)
  thisPredVarValues <- NULL
  for(v in 1:npvars){
    thisVar <- thisPredVars[v]
    if(length(grep("ChaosRun", thisVar))>0){
      # thisLevel <- gsub("ChaosRun","", thisVar)
      # thisVar <- "ChaosRun"
      # tempVarValues <- thisLevel
      tempVarValues<-thisChaosRuns
    } else{
      tempVarValues<-getVarValues(thisVar, thisData = this_model$data)
    }
    ## fix the odd things - don't bother with zero biomass, lifespan, Linf, 
    # and keystonennes and responsiveness should just have unique values
    if(thisVar %in% c("B0", "Lifespan", "Linf", "NumL1cons", "PropByTopPrey","TL")){
      tempVarValues <- tempVarValues[tempVarValues!=0]
      # check not outside data limits
      dataMax <- max(this_model$data[,c(thisVar)], na.rm=TRUE); dataMin <- min(this_model$data[,c(thisVar)], na.rm=TRUE)
      tempVarValues <- tempVarValues[tempVarValues<=dataMax & tempVarValues>=dataMin]
      # } else if(thisVar %in% c("Keystone", "Informance")){
    } else if(thisVar %in% c("Informance", "Keystone","Response")){
      thisUnique <- sort(unique(this_model$data[,c(thisVar)]))
      tempVarValues <- thisUnique
    }
    # if(thisVar=="ChaosRun" & grep("ChaosRun", names(thisPredVarLengths)))
    thisPredVarValues[[v]] <- tempVarValues
    thisPredVarLengths[v]<-length(thisPredVarValues[[v]])
  }
  
  thisIntVars <- (unlist(str_split(thisIntVar,":")))
  thisValues1 <- as.double(thisPredVarValues[[grep(thisIntVars[1], thisPredVars)[1]]])
  thisValues2 <- as.double(thisPredVarValues[[grep(thisIntVars[2], thisPredVars)[1]]])
  if(is.factor(this_model$data[,c(thisIntVars[1])])){
    thisValues1<-thisPredVarValues[[grep(thisIntVars[1], thisPredVars)]]
  }
  if(is.factor(this_model$data[,c(thisIntVars[2])])){
    thisValues2<-unlist(thisPredVarValues[grep(thisIntVars[2], thisPredVars)])
  }
  test_df <- data.frame(matrix(0, ncol=dim(this_model$data)[2], nrow=(length(thisValues1) * length(thisValues2))))
  colnames(test_df)<- colnames(this_model$data)
  for(s in 1:npvars){
    fillVar <- thisPredVars[s]
    if(length(grep(fillVar, colnames(thisData)))>0){
      if( is.factor(thisData[,fillVar])==FALSE){
        temp<- boxplot(this_model$data[,c(fillVar)], plot=FALSE)$stats
        fillValueMed <- temp[3]; 
        test_df[,c(fillVar)]<- fillValueMed; 
      } else{
        test_df[,c(fillVar)]<-names(sort(table(thisData[,fillVar]),decreasing = TRUE)[1])
        test_df[,c(fillVar)] <- factor(test_df[,c(fillVar)], levels=sort(unique(thisData[,c(fillVar)])))
      }
    }
  }
  test_df[,c(thisIntVars[1])] <- sort(rep(thisValues1, length(thisValues2)));
  test_df[,c(thisIntVars[2])] <-  rep(thisValues2, length(thisValues1))
  if(is.factor(this_model$data[,c(thisIntVars[1])])){
    test_df[,c(thisIntVars[1])] <- factor(test_df[,c(thisIntVars[1])], levels=levels(this_model$data[,c(thisIntVars[1])]))
  }
  if(is.factor(this_model$data[,c(thisIntVars[2])])){
    test_df[,c(thisIntVars[2])] <- factor(test_df[,c(thisIntVars[2])], levels=levels(this_model$data[,c(thisIntVars[2])]))
  }
  test_df$ChaosRun <- factor(test_df$ChaosRun, levels=thisChaosRuns)
  
   # for(a in 1:dim(test_df)[1]){
  #   this_v1 <- test_df[a,c(thisIntVars[1])]; this_v2 <- test_df[a, c(thisIntVars[2])]
  #   test<- testPlausiblePair(value1=this_v1, value2=this_v2, var1=thisIntVars[1], var2 = thisIntVars[2], thisData=thisData)
  #   if(is.na(test)){
  #      test_df[a,c(thisIntVars)] <- NA
  #   }
  # }
  
  test_df$Predicted <- ((predict.glm(this_model, newdata = test_df)))^5
  if(updateMax==TRUE){
    thisMax <<- max(test_df$Predicted, na.rm=TRUE); 
    # thisMin <- min(test_df$Predicted, na.rm=TRUE)
  } else{
    thisMax <<- inputMax
  }
  cat(paste("-- ",thisMax, thisMin, sep=", "))
  plotColour <<- unlist(lapply(test_df$Predicted, getColor, log=FALSE, thisMax=thisMax))
  # tempDF<-data.frame(cbind("x"=test_df[,c(thisIntVars[1])]),"y"=test_df[,c(thisIntVars[2])])
  x1<-sort(rep(seq(1,length(thisValues1)), length(thisValues2))); y1 <- rep(seq(1,length(thisValues2)), length(thisValues1))
  # cat("x1 = ", length(x1))
  # cat(x1,"\n")
  # cat("y1 = ", length(y1))
  # cat(y1, "\n")
  # # tempDF <- data.frame(cbind("x"=x1,"y"=y1))
  testDF <- cbind(x1,y1)
  # cat(dim(testDF),"\n")
  
  tempDF<<-data.frame(testDF); colnames(tempDF)<-c("x", "y")
# 
#   # pdf(paste(plotPath,"KemptonsQ_SSRsens.pdf",sep=""),height=4,width=10)
  plot(x=tempDF$x,y=tempDF$y,type="n",xlab="",ylab="",yaxt="n", xaxt="n")
  temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
  par(las=1)
  xx <- thisValues1; xxx <- seq(1, length(thisValues1))
  if(length(xx)>5){
    aIndex<-seq(1,length(xx), by=ceiling(length(xx)/5))
    axis1Labs  <- xx[aIndex]; axis1At <- xxx[aIndex]
  }else{
    axis1Labs  <- xx; axis1At <- xxx
  }
  xx <- thisValues2; xxx <- seq(1, length(thisValues2))
  if(length(xx)>6){
    aIndex<-seq(1,length(xx), by=ceiling(length(xx)/5))
    axis2Labs  <- xx[aIndex]; axis2At <- xxx[aIndex]
  } else{
    axis2Labs  <- xx; axis2At <- xxx
  }
  axis(at=axis1At, labels=axis1Labs, side=1)
  axis(at=axis2At, labels=axis2Labs, side=2)
  par(las=0)
  thisLabel1 <- varDesc$Decr[varDesc$Variable == thisIntVars[1]]
  if(thisLabel1==""){
    thisLabel1 <- expression(B[0]~"('000 tonnes)")
  }
  thisLabel2 <- varDesc$Decr[varDesc$Variable == thisIntVars[2]]
  if(thisLabel2==""){
    thisLabel2 <- expression(B[0]~"('000 tonnes)")
  }
  mtext(thisLabel1, side=1, adj=0.5, line=2.5, cex=lab_cex); mtext(thisLabel2, side=2, adj=0.5, line=4, cex=lab_cex)
  box()
  if(updateMax==TRUE){
     mtext(paste("Max. relative change: ", 100*round(thisMax, 2), "%", sep=""), side=3, adj=0)
    # rm(thisMax)
  }
 # plot(x=test_df[,c(thisIntVars[1])], y=test_df$Predicted*100, xlab=thisIntVars[1], ylab="Predicted relative difference (%)", pch=20, col=myOrange); 
 # plot(x=test_df[,c(thisIntVars[2])], y=test_df$Predicted*100, xlab=thisIntVars[2], ylab="Predicted relative difference (%)", pch=20, col=myOrange)
}

plotSingleEffect <- function(thisVar, this_model){
  thisVarValues <- sort(unique(this_model$data[,c(thisVar)])); nvv <- length(thisVarValues)
  testMed<-data.frame(matrix(0, ncol=dim(this_model$data)[2], nrow=nvv)); colnames(testMed)<-colnames(this_model$data); 
  temp <- as.character(names(this_model$coefficients)[-1]); thisModelPredictors <-unique(unlist(str_split(temp,":")))
  for(i in 1:length(thisModelPredictors)){
    fillVar <- thisModelPredictors[i]
    if( is.factor(this_model$data[,fillVar])==FALSE){
      temp<- boxplot(this_model$data[,c(fillVar)], plot=FALSE)$stats
      fillValueMed <- temp[3]; 
      testMed[,c(fillVar)]<- fillValueMed; 
    } else{
      testMed[,c(fillVar)]<-names(sort(table(this_model$data[,fillVar]),decreasing = TRUE)[1])
      testMed[,c(fillVar)] <- factor(testMed[,c(fillVar)], levels=sort(unique(this_model$data[,c(fillVar)])))
    }
  }
  testMed[,c(thisVar)]<-thisVarValues;
  thisPred <- predict.glm(this_model,newdata = testMed)
  transPred <- unlist(lapply(thisPred, FUN= function(x){x^5}))*100
  thisYmin <- min(c(exp(c( thisPred)),0), na.rm=TRUE); thisYmax <- max(exp(c(thisPred)), na.rm=TRUE)
  par(las=1)
  plot(x=thisVarValues, y=transPred, type="l", xlab="", ylab="", cex.axis=lab_cex, cex.lab=lab_cex) 
  thisLabel1 <- varDesc$Decr[varDesc$Variable == thisVar]
  if(thisLabel1==""){
    thisLabel1 <- expression(B[0]~"('000 tonnes)")
  }
  par(las=0)
  mtext(thisLabel1, side=1, adj=0.5, line=2.5, cex=lab_cex); mtext("Relative change in biomass (%)", side=2, adj=0.5, line=4, cex=lab_cex)
  
}

getPercentile<-function(x, k){
  y <- sort(x[!is.na(x)]); ny <- length(y)
  i1<-round(ny*(k/2)); i2 <- round(ny*(1-(k/2)))
  return(y[c(i1, i2)])
}
getColor <- function(x, log=FALSE, base=10, thisMax=1){
  thisCol<-"white"
  if(!is.na(x)){
    thisIndex <- round((x-thisMin)/(thisMax-thisMin),2)*100 +1
    if(log==TRUE){
      if(x<1){x<-1}
      thisIndex <- round((log(x, base=base)/log(thisMax, base=base)),2)*100 +1
    }
    if(length(thisIndex)>0){
      if(thisIndex > 101){
        thisIndex <- 101
      }
      if(sum(x, na.rm=TRUE)>0 & thisIndex>0){
        thisCol<-absColorRamp[thisIndex]
      }
    }
  }
  return(thisCol)
}
plotGrid<-function(x,y){
  thisX<-c(x-0.5,x-0.5,x+0.5,x+0.5); thisY<-c(y-0.5,y+0.5,y+0.5,y-0.5)
  index <- tempDF$x==x & tempDF$y==y
  thisCol<-plotColour[index]
  if(length(thisCol)>0){
    polygon(x=thisX,y=thisY,col=thisCol,border=NA)
  }
  return(NULL)
}

getBestLModel <- function(thisVar, thisFactors, cutoff, data2fit){
  this_df <- data.frame(matrix(NA, ncol=2, nrow=length(thisFactors))); colnames(this_df)<-c("Factor", "Rsq")
  existingRsq <- 0
  modelFactors <- c()
  leftFactors<-thisFactors
  keepGoing <- TRUE
  i <- 1
  while(keepGoing==TRUE){
    testRsq <- c(); keep4nowFactors<-c()
    for(f in 1:length(leftFactors)){
      testFactors <- c(modelFactors, leftFactors[f])
      testFormula <- as.formula(paste(thisVar,"~", paste(testFactors, collapse="+")))  
      testModel <- glm(testFormula, data=data2fit)
      thisRsq <- 1-(summary(testModel)$"deviance"/summary(testModel)$"null.deviance")
      if((thisRsq-existingRsq) >cutoff){
        keep4nowFactors <- c(keep4nowFactors, leftFactors[f])
        testRsq <- c(testRsq, thisRsq)
      }
    }
    if(length(keep4nowFactors)>0){
      if((max(testRsq)-existingRsq)>cutoff){
        # then keep the top one
        index <- testRsq == max(testRsq)
        thisFactor <- keep4nowFactors[index]; thisRsq <- testRsq[index]
        thisRsq_improvement <- thisRsq - existingRsq
        this_df[i,] <- c(thisFactor, thisRsq)
        leftFactors <- leftFactors[leftFactors != thisFactor]
        existingRsq <- thisRsq
        modelFactors <- c(modelFactors, thisFactor)
        i <- i +1
      }else{
        keepGoing <- FALSE
      }
    } else{
      keepGoing=FALSE
    }
  }
  index<-!is.na(this_df$Factor)
  return(this_df[index,])
}
getIntVars <- function(vars){
  nvars <- length(vars)
  allInts <- c()
  for(i in 1:nvars){
    for(j in (i):nvars){
      if(i != j){ 
        thisInt <- paste(vars[i],":", vars[j], sep="")
        checkIntRev <- paste(vars[j],":", allInts[i], vars="")
        if(!(checkIntRev %in% allInts)){
          allInts <- c(allInts, thisInt)
        }
      }
    }
  }
  return(allInts)
}
