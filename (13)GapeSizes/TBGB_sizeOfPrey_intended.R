##based on intended diets, defined for adults and juveniles, what is the approximate mean size and variance of size of prey for each predator
## will need to put approx values for biomass pool prey

source(paste(DIR$'General functions',"get_interaction_trophicLevel.R",sep=""))
source(paste(DIR$'General functions',"get_first_number.R",sep=""))
source(paste(DIR$'General functions',"myRounding.R",sep=""))
source(paste(DIR$'General functions',"getSNRN.R",sep=""))
source(paste(DIR$'General functions',"nonZeroMean.R",sep=""))
source(paste(DIR$'General functions',"getICNumbers.R",sep=""))


thisPath <-  paste(DIR$'Base', "TBGB\\",sep="")
plotPath<-paste(thisPath,"\\Figures\\Testing\\DIET\\",thisDesc, sep="")

# bring in bp_size and storeSizeAtAge - created in TBGB_initialConditionsDump.R - in setting_up/Biology_initial
load(paste(thisPath, "\\SizeAtAge_intialConditions", sep=""))
biolLines  <- readLines( paste(thisPath, "TBGB_JP\\TBGB_biol.prm", sep=""))

getAgeMat <- function(Code){
  x <- paste(Code,"_age_mat", sep="")
  y <- get_first_number(biolLines[grep(x, biolLines)])
  return(y)
}
groupsDF<-read.csv(paste(thisPath,"\\TBGB_JP\\TBGB_Groups.csv",sep="")); ng<-dim(groupsDF)[1]

year0<-1899; #this is when the model starts - includes 35 year burn-in that takes it up to the real start of 1900

dietYears<-seq(1900,1980) #just replace with 1970 if only want one year, or give a range of years - will give average over them
dietTimeSteps<-dietYears-year0 #to index dietYears

preyGroupsDF<-read.csv(paste(thisPath,"\\TBGB_PreyGroups.csv",sep=""))
preyGroups<-sort(unique(preyGroupsDF$PreyGroup)); nPreyGroups<-length(preyGroups)

## first plot intended diets
predators<-as.character(groupsDF$Code[groupsDF$IsPredator==1]); npreds<-length(predators)

getAgeColor<-function(x){
  thisCol<-myGrey
  if(is.na(x)){
    thisCol<-myGrey
  } else{
    if(x==2){thisCol="black"}
  }
  return(thisCol)
}
getAgeText<-function(x){
  thisText<-"Juv"
  if(is.na(x)){
    thisText<-""
  }else{
    if(x==2){thisText="Adult"}
  }
  return(thisText)
}

preyGroupColours<-colorRampPalette(colors=c("black", myGrey,"midnightblue",myBlue,myLightBlue,myLightAqua, myAqua,myDarkGreen,myGreen,myGold,myOrange,"red","brown",myRed))(nPreyGroups)
#####################################

## do intended diets
## colnames will have predators, with pPREY, prey age, PredCode, pred age
predCode<-groupsDF$Code[groupsDF$IsPredator==1]; npreds <- length(predCode)
ageIndex <- groupsDF$NumCohorts[groupsDF$IsPredator==1]>1
agedVrs <- paste("pPREY",c(1,1,2,2), sort(rep(predCode[ageIndex],4)),c(1,2,1,2), sep="")
noageVrs <- paste("pPREY",predCode[!ageIndex], sep="")

predVars<-sort(c(agedVrs, noageVrs)); npredVrs<-length(predVars)
nPrey <- dim(groupsDF)[1]

dietDF <- data.frame(matrix(NA, ncol=npredVrs, nrow=nPrey))
colnames(dietDF)<- predVars; rownames(dietDF)<- as.character(groupsDF$Code)
dietsOriginal <- read.csv(paste(basePath, "..\\TBBOriginalDiet.csv",sep=""))

preyCodesOrig<-colnames(dietsOriginal)[7:(7+ng-1)]; preyCodeIndex<-match(groupsDF$Code, preyCodesOrig)

for(j in 1:npredVrs){
  thisPredVar <- predVars[j]
  preyAge<-get_first_number(thisPredVar, n=1); predAge<-get_first_number(thisPredVar, n=2); thisCode<-gsub(paste("pPREY", predAge,preyAge,sep="|"), "", thisPredVar)
  xx<-grep(thisCode, dietsOriginal$Code); 
  if(!is.na(preyAge) & !is.na(predAge)){
    xxx<-dietsOriginal[xx:(xx+3),]
    predAgeText<-c("Juvenile", "Adult")[predAge]; preyAgeText<-c("Juvenile", "Adult")[preyAge]
  } else{
    xxx<-dietsOriginal[xx,]
    predAgeText<-"All"; preyAgeText<-"All"
  }  
  index<-xxx$Predator.age==predAgeText & xxx$Prey.age==preyAgeText
  thisDiets<-xxx[index,]
  newColumn<-thisDiets[7:dim(dietsOriginal)[2]]
  dietDF[1:ng,j]<-as.double(newColumn[preyCodeIndex])
  
}
preys<-dietDF[,1]; npreys<-length(preys)
diets<-dietDF

getSizeRange <- function(Code, mature=NA){
  thisNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==Code]
  if(thisNumCohorts>1){
    thisAgeMat <- getAgeMat(Code)
     thisJuvCohorts <- seq(1:thisAgeMat); thisAdCohorts <- seq((thisAgeMat+1):thisNumCohorts)
    test <- grep("ad", mature, ignore.case = TRUE)
    if(length(test)==0){
      thisAges <- storeSizeAtAge[grep(Code, groupsDF$Code),thisJuvCohorts]
    } else{
      thisAges <- storeSizeAtAge[grep(Code, groupsDF$Code), thisAdCohorts]
    }
    this_out <- c(min(thisAges, na.rm=TRUE), median(thisAges, na.rm=TRUE), max(thisAges, na.rm=TRUE))
  }else{
    this_out <- bp_size[bp_size$Code==Code,c("Min", "Median", "Max")]
  }
  return(as.double(this_out))
}
## with age structure - only do for predators with age structure
dietsByPreyGroup<-array(NA, dim=c(npreds, nPreyGroups)); colnames(dietsByPreyGroup)<-preyGroups; rownames(dietsByPreyGroup)<-predators
storePreySizeRange <- array(NA, dim=c(npreds, 2, 3))
for(p in 1:npreds){
  thisPred<-predators[p]; thisNumCohorts <- groupsDF$NumCohorts[groupsDF$Code==thisPred]
  if(thisNumCohorts>1){
    temp<-diets[,grep(thisPred,colnames(diets))]; thisDiet<-temp
    predAge<-unlist(lapply(colnames(temp), get_first_number, n=2))
    preyAge <- unlist(lapply(colnames(temp), get_first_number, n=1))
    adultDiet <- temp[,predAge==2]; adSums <- apply(adultDiet,1, sum, na.rm=TRUE)
    adultDiet <- adultDiet[adSums>0,]
    # any prey that are not age-structured need to have the same (use highest) values in both adult prey and juvenile prey ages
    fixRows <- rownames(adultDiet) %in% groupsDF$Code[groupsDF$NumCohorts==1]
    maxByRow <- apply(adultDiet, 1, max)
    adultDiet[fixRows,]<- maxByRow[fixRows]
    juvDiet <- temp[,predAge==1]; juvSums <- apply(juvDiet, 1, sum, na.rm=TRUE)
    juvDiet <- juvDiet[juvSums>0,]
    # any prey that are not age-structured need to have the same (use highest) values in both adult prey and juvenile prey ages
    fixRows <- rownames(juvDiet) %in% groupsDF$Code[groupsDF$NumCohorts==1]
    maxByRow <- apply(juvDiet, 1, max)
    juvDiet[fixRows,]<- maxByRow[fixRows]
    # prey sizes
    # juveniles
    xx <- unlist(lapply(rownames(juvDiet), getSizeRange, mature="juvenile"))
    juvSizes1 <- matrix(xx, ncol=3, nrow=length(xx)/3, byrow=TRUE)
    xx <- unlist(lapply(rownames(juvDiet), getSizeRange, mature="adult"))
    juvSizes2 <- matrix(xx, ncol=3, nrow=length(xx)/3, byrow=TRUE)
    thisMean <- sum(juvSizes1[,2] * juvDiet[,1]) + sum(juvSizes2[,2] * juvDiet[,2])
    thisMin <- min(c(juvSizes1[juvDiet[,1]>0,1]), c(juvSizes2[juvDiet[,1]>0,1] ))
    thisMax <- max(c(juvSizes1[juvDiet[,1]>0,3]), c(juvSizes2[juvDiet[,2]>0,3] ))
    storePreySizeRange[p,1,]<- c(thisMin, thisMean, thisMax)
    # adults
    xx <- unlist(lapply(rownames(adultDiet), getSizeRange, mature="juvenile"))
    adSizes1 <- matrix(xx, ncol=3, nrow=length(xx)/3, byrow=TRUE)
    xx <- unlist(lapply(rownames(adultDiet), getSizeRange, mature="adult"))
    adSizes2 <- matrix(xx, ncol=3, nrow=length(xx)/3, byrow=TRUE)
    thisMean <- sum(adSizes1[,2] * adultDiet[,1]) + sum(adSizes2[,2] * adultDiet[,2])
    thisMin <- min(c(adSizes1[adultDiet[,1]>0,1]), c(adSizes2[adultDiet[,2]>0,1] ))
    thisMax <- max(c(adSizes1[adultDiet[,1]>0,3]), c(adSizes2[adultDiet[,2]>0,3] ))
    storePreySizeRange[p,2,]<- c(thisMin, thisMean, thisMax)
  } else{
    temp<-diets[,grep(thisPred,colnames(diets))]; thisDiet<-temp[temp>0]; 
    thisPreys <- rownames(diets)[temp>0]
    xx <- unlist(lapply(thisPreys, getSizeRange, mature=NA))
    thisPreySizes <- matrix(xx, ncol=3, nrow=length(xx)/3, byrow=TRUE)
    thisMean <- sum(thisPreySizes[,2] * thisDiet)
    thisMin <- min(c(thisPreySizes[thisDiet>0]))
    thisMax <- max(c(thisPreySizes[thisDiet>0]))
    storePreySizeRange[p,1,]<- c(thisMin, thisMean, thisMax); storePreySizeRange[p,2,]<- c(thisMin, thisMean, thisMax)
  }
}
    

# for each predator, plot mean and lims of prey size
thismin <- min(storePreySizeRange, na.rm=TRUE); thismax <- max(storePreySizeRange, na.rm=TRUE)
par(las=1, mar=c(4,5,1,1))
plot(y=1:npreds, x=rep(1, npreds), type="n", ylab="",  xlab="Prey size (cm)", xlim=c(thismin, thismax), yaxt="n")
axis(at=1:npreds, labels = predators, side=2)
for(p in 1:npreds){
  thisCode <- predators[p]; 
  groupIndex <- groupsDF$Code==thisCode
  thisNumCohorts<-groupsDF$NumCohorts[groupIndex]
  thisY<- c(p-0.1, p-0.1, p+0.1, p+0.1)
  if(thisNumCohorts>1){
    thisJ <- storePreySizeRange[p, 1,]; thisA <- storePreySizeRange[p, 2, ]
    polygon(x=thisJ[c(1,3,3,1)], y=thisY, col=myGreen_trans, border=myGreen)
    points(x=thisJ[2], y=p, pch=8, col=myGreen)
    polygon(x=thisA[c(1,3,3,1)], y=thisY, col=myBlue_trans, border=myBlue)
    points(x=thisA[2], y=p, pch=8, col=myBlue)
  } else{
    thisA <- storePreySizeRange[p, 1, ]
    polygon(x=thisA[c(1,3,3,1)], y=thisY, col=myBlue_trans, border=myBlue)
    points(x=thisA[2], y=p, pch=8, col=myBlue)
  }
}


thismin <- min(storePreySizeRange, na.rm=TRUE); thismax <- max(storePreySizeRange, na.rm=TRUE)
thisXaxis <- c(1e-6, 1e-4, 1e-2, 1, 100); thisXaxis_at <- log(thisXaxis, base=10)
par(las=1, mar=c(4,5,1,1))
plot(y=1:npreds, x=rep(1, npreds), type="n", ylab="",  xlab="Prey size (cm)", xlim=log(c(thismin, thismax), base=10), yaxt="n", xaxt="n")
axis(at=1:npreds, labels = predators, side=2)
axis(at=thisXaxis_at, labels=thisXaxis, side=1)
for(p in 1:npreds){
  thisCode <- predators[p]; 
  groupIndex <- groupsDF$Code==thisCode
  thisNumCohorts<-groupsDF$NumCohorts[groupIndex]
  thisYj<- c(p-0.2, p-0.2, p, p); thisYa <- c(p, p, p+0.2, p+0.2)
  thisY <- c(p-0.2, p-0.2, p+0.2, p+0.2)
  if(thisNumCohorts>1){
    thisJ <- log(storePreySizeRange[p, 1,],base=10); thisA <- log(storePreySizeRange[p, 2, ], base=10)
    polygon(x=thisJ[c(1,3,3,1)], y=thisYj, col=myGreen_trans, border=myGreen)
    points(x=thisJ[2], y=p, pch=8, col=myGreen)
    polygon(x=thisA[c(1,3,3,1)], y=thisYa, col=myBlue_trans, border=myBlue)
    points(x=thisA[2], y=p, pch=8, col=myBlue)
  } else{
    thisA <- log(storePreySizeRange[p, 1, ], base=10)
    polygon(x=thisA[c(1,3,3,1)], y=thisY, col=myBlue_trans, border=myBlue)
    points(x=thisA[2], y=p, pch=8, col=myBlue)
    
  }
}


# how about the max's and means - figure they can all eat pretty small
thismin <- min(storePreySizeRange[,,2], na.rm=TRUE); thismax <- max(storePreySizeRange[,,2], na.rm=TRUE)
thisXaxis <- c(1e-6, 1e-4, 1e-2, 1, 100); thisXaxis_at <- log(thisXaxis, base=10)
par(las=1, mar=c(4,5,1,1))
plot(y=1:npreds, x=rep(1, npreds), type="n", ylab="",  xlab="Prey size (cm)", xlim=log(c(thismin, thismax), base=10), yaxt="n", xaxt="n")
axis(at=1:npreds, labels = predators, side=2)
axis(at=thisXaxis_at, labels=thisXaxis, side=1)
for(p in 1:npreds){
  thisCode <- predators[p]; 
  groupIndex <- groupsDF$Code==thisCode
  thisNumCohorts<-groupsDF$NumCohorts[groupIndex]
  thisYj<- c(p-0.2, p-0.2, p, p); thisYa <- c(p, p, p+0.2, p+0.2)
  thisY <- c(p-0.2, p-0.2, p+0.2, p+0.2)
  if(thisNumCohorts>1){
    thisJ <- log(storePreySizeRange[p, 1,],base=10); thisA <- log(storePreySizeRange[p, 2, ], base=10)
    # polygon(x=thisJ[c(1,3,3,1)], y=thisYj, col=myGreen_trans, border=myGreen)
    points(x=thisJ[2], y=p, pch=8, col=myGreen)
    # polygon(x=thisA[c(1,3,3,1)], y=thisYa, col=myBlue_trans, border=myBlue)
    points(x=thisA[2], y=p, pch=8, col=myBlue)
  } else{
    thisA <- log(storePreySizeRange[p, 1, ], base=10)
    # polygon(x=thisA[c(1,3,3,1)], y=thisY, col=myBlue_trans, border=myBlue)
    points(x=thisA[2], y=p, pch=8, col=myBlue)
    
  }
}

this_df <- data.frame(rbind(storePreySizeRange[,1,], storePreySizeRange[,2,])); 
this_df$Code <- rep(predators, 2); this_df$Age <- c(rep("juv",npreds), rep("ad", npreds))
this_df <- this_df[order(this_df$Code),]
colnames(this_df)[1:3]<- c("min", "mean", "max")
write.csv(this_df, paste(thisPath, "\\SizeBased\\PreySizeRanges.csv", sep=""), row.names = FALSE )

## max only
# how about the max's and means - figure they can all eat pretty small
predOrder <- order(predators)
thismin <- min(storePreySizeRange[,,2], na.rm=TRUE); thismax <- max(storePreySizeRange[,,2], na.rm=TRUE)
thisXaxis <- c(1e-6, 1e-4, 1e-2, 1, 100); thisXaxis_at <- log(thisXaxis, base=10)
par(las=1, mar=c(4,5,1,1))
plot(y=1:npreds, x=rep(1, npreds), type="n", ylab="",  xlab="Prey size (cm)", xlim=log(c(thismin, thismax), base=10), yaxt="n", xaxt="n")
abline(h=seq(0,npreds, by=2), col=myGrey_trans)
axis(at=1:npreds, labels = predators[predOrder], side=2)
axis(at=thisXaxis_at, labels=thisXaxis, side=1)
for(p in 1:npreds){
  thisCode <- predators[predOrder][p]; 
  groupIndex <- groupsDF$Code==thisCode
  thisNumCohorts<-groupsDF$NumCohorts[groupIndex]
  thisYj<- c(p-0.2, p-0.2, p, p); thisYa <- c(p, p, p+0.2, p+0.2)
  thisY <- c(p-0.2, p-0.2, p+0.2, p+0.2)
  if(thisNumCohorts>1){
    thisJ <- log(storePreySizeRange[p, 1,],base=10); thisA <- log(storePreySizeRange[p, 2, ], base=10)
    # polygon(x=thisJ[c(1,3,3,1)], y=thisYj, col=myGreen_trans, border=myGreen)
    points(x=thisJ[2], y=p, pch=8, col=myGreen)
    # polygon(x=thisA[c(1,3,3,1)], y=thisYa, col=myBlue_trans, border=myBlue)
    points(x=thisA[2], y=p, pch=8, col=myBlue)
  } else{
    thisA <- log(storePreySizeRange[p, 1, ], base=10)
    # polygon(x=thisA[c(1,3,3,1)], y=thisY, col=myBlue_trans, border=myBlue)
    points(x=thisA[2], y=p, pch=8, col=myOrange)
  }
}


# suppose the limits are the 95% CIs, then approx. standard deviation, and then CVs
getCV <- function(max, mean){
  x <- max-mean
  sd <- x/1.96
  cv <- sd/mean
  return(cv)
}
CVbyPred <- matrix(mapply(getCV, max=storePreySizeRange[,,3], mean=storePreySizeRange[,,2]), ncol=2, nrow=dim(storePreySizeRange)[1], byrow = TRUE)
rownames(CVbyPred)<-predators













