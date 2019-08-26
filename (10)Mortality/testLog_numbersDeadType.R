thisModel <- "Base"; thisRun <- "TestMORT2"
thisRun <- "TestMORTmL"

thisPath <- paste(DIR$'Base',"ATLANTISmodels\\" ,thisModel, "\\output", thisRun,"\\", sep="")

thisCode <- "HAK"

thisLogFile <- paste(thisPath, "log.txt", sep="")
logLines <- readLines(thisLogFile)
x <- grep("Vert_MortalityTEST", logLines)
mortLines <- logLines[x]


getCohort <- function(x){
  y<-unlist(str_split(x, paste(thisCode,":", sep="")))[2]
  yy<-get_first_number(y)
  return(yy)
}
getNumTotalDead <- function(x){
  y<- unlist(str_split(x, ":|,"))
  z <- y[grep("Number dead are", y)]
  thisNum <- get_first_number(z)
  return(thisNum)
}
getApplied_mL <- function(x){
  y<- unlist(str_split(x, ":|,"))
  z <- y[grep("mL", y)+1]
  thisNum <- get_first_number(z)
  return(thisNum)
}
getNums <- function(x){
  y<- unlist(str_split(x, ":|,"))
  z <- y[grep("Nums", y)+1]
  thisNum <- get_first_number(z)
  return(thisNum)
}

numsecs <- 12*60*60
input_mL <- c(5e-5, 8.5e-4)/numsecs

allApplied_mL <- unlist(lapply(mortLines, getApplied_mL))
allNatDead <- unlist(lapply(mortLines, getNumTotalDead))
allNums <- unlist(lapply(mortLines, getNums))
allProps <- allNatDead/allNums

# storeApplied_mL <- allApplied_mL

test <- allApplied_mL/allProps

hist(allApplied_mL)
abline(v=input_mL, col="red")

max(allApplied_mL, na.rm=TRUE) * numsecs

## temperature
inputsPath <- paste(DIR$'Base',"ATLANTISmodels\\inputs\\ROMS\\Chatham30_tempAll.nc", sep="")
ThisIC.nc <- nc_open(inputsPath)
thisTracer <- "temperature"
thisData <- ncvar_get(ThisIC.nc, thisTracer)
ndays <- round(dim(thisData)[3]/2)
dayIndex <- seq(1, ndays, by=2); nightIndex <- seq(2,ndays, by=2)
dayTemps <- thisData[,,dayIndex]
nightTemps <- thisData[,,nightIndex]
maxByTS<-apply(thisData, 3, max, na.rm=TRUE)

calc_tempScalar <- function(temperature){
  y <- 2^((temperature - 15)/10)
  return(y)
}

tempScalars <- apply(thisData, c(1,2,3), calc_tempScalar)
