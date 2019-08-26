thisModel <- "Base"; thisRun <- "TestH2Otemp"

thisPath <- paste(DIR$'Base',"ATLANTISmodels\\" ,thisModel, "\\output", thisRun,"\\", sep="")

thisCode <- "HAK"

thisLogFile <- paste(thisPath, "log.txt", sep="")
logLines <- readLines(thisLogFile)
x <- grep("H2OtemperatureTest", logLines)
H2OLines <- logLines[x]


getTemp <- function(x){
  y<-unlist(str_split(x, paste(thisCode,":| ", sep="")))
  z <- grep("Temperature", y)+1
  yy<-get_first_number(y[z])
  return(yy)
}
getBox <- function(x){
  y<-unlist(str_split(x, paste(thisCode,":| ", sep="")))
  z <- grep("Box", y)+1
  yy<-get_first_number(y[z])
  return(yy)
}
getLayer <- function(x){
  y<-unlist(str_split(x, paste(thisCode,":| ", sep="")))
  z <- grep("Layer", y)+1
  yy<-get_first_number(y[z])
  return(yy)
}
calc_tempScalar <- function(temperature){
  y <- 2^((temperature - 15)/10)
  return(y)
}

allTemps <- unlist(lapply(H2OLines, getTemp))
allTcorrs <- unlist(lapply(allTemps, calc_tempScalar))

hist(allTcorrs)

######################################
## TESTING TCORR
thisModel <- "Base"; thisRun <- "TestH2Otemp"

thisPath <- paste(DIR$'Base',"ATLANTISmodels\\" ,thisModel, "\\output", thisRun,"\\", sep="")

thisLogFile <- paste(thisPath, "log.txt", sep="")
logLines <- readLines(thisLogFile)
x <- grep("H2OtemperatureTest", logLines)
H2OLines <- logLines[x]

getTcorr <- function(x){
  y<-unlist(str_split(x, paste(thisCode,":| ", sep="")))
  z <- grep("Tcorr", y)+1
  yy<-get_first_number(y[z])
  return(yy)
}

allTcorrs <- unlist(lapply(H2OLines, getTcorr))
hist(allTcorrs)
thisALlTemps <- unlist(lapply(H2OLines, getTemp))
thisCalcTCorrs <- unlist(lapply(thisALlTemps, calc_tempScalar))
