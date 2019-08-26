thisModel <- "Base"; thisRun <- "TestMORT2"
thisRun <- "TestH2Otemp"

thisPath <- paste(DIR$'Base',"ATLANTISmodels\\" ,thisModel, "\\output", thisRun,"\\", sep="")

thisCode <- "HAK"

thisLogFile <- paste(thisPath, "log.txt", sep="")
logLines <- readLines(thisLogFile)
x <- grep("Number dead are", logLines)
deadLines <- logLines[x]
getCohort <- function(x){
  y<-as.double(unlist(str_split(x, ":| "))[3])
  return(y)
}

wchab <- deadLines[grep("WChabitat", deadLines)]
wchab <- wchab[1:(length(wchab)-1)] # take off the last line in case incomplete
otherHab <- deadLines[grep("WChabitat", deadLines, invert = TRUE)]
otherHab <- otherHab[1:(length(otherHab)-1)] # take off the last line in case incomplete


getNumDead <- function(x){
  y<- unlist(str_split(x, ":|,"))
  z <- y[grep("Number dead", y)]
  thisNum <- get_first_number(z)
  return(thisNum)
}
getNumDead_Epi <- function(x){
  y<- unlist(str_split(x, ":|,"))
  z <- y[grep("Number dead", y)]
  thisNum <- get_first_number(z, n=2)
  return(thisNum)
}
getNumEaten <- function(x){
  y <- unlist(str_split(x, ":|,"))
  z <- y[grep("Eaten", y)]
  thisNum <- get_first_number(z)
  return(thisNum)
}
getTotalNumFlux <- function(x){
  y <- unlist(str_split(x, ":|,"))
  z <- y[(1+grep("tot num changed", y))]
  thisNum <- get_first_number(z)
  return(thisNum)
}
getDayT <- function(x){
  y <- unlist(str_split(x, ":|,"))
  z <- y[(1+grep("Day", y))]
  thisNum <- get_first_number(z)
  return(thisNum)
}
getCohorts <- function(x){
  y <- unlist(str_split(x, ":|,"))
  z <- y[(1+grep(paste("Group ",thisCode, sep=""), y))]
  thisNum <- get_first_number(z)
  return(thisNum)
}
getColByNight <- function(x){
  thisCol <- "midnightblue"
  if(round(x)==x){
    thisCol <- myGold
  }
  return(thisCol)
}
getTotalNumbers <- function(x){
  y <- unlist(str_split(x, ":|,"))
  z <- y[(1+grep("total number", y))]
  thisNum <- get_first_number(z)
  return(thisNum)
}
WC_totalNumbers <- unlist(lapply(wchab, getTotalNumbers))
WC_numsDead <- unlist(lapply(wchab, getNumDead))
WC_numsEaten <- unlist(lapply(wchab, getNumEaten))
WC_numsFlux <- unlist(lapply(wchab, getTotalNumFlux))
WC_day_ts<- unlist(lapply(wchab, getDayT))
WC_cohorts <- unlist(lapply(wchab, getCohorts))
colByNight <- unlist(lapply(WC_day_ts, getColByNight))
cohortColors<-colorRampPalette(colors=c(myGold, myGreen,myAqua,myBlue, "midnightblue", myPurple, myRed))(10)
WC_cohortColors <- cohortColors[(WC_cohorts +1)]
plot(WC_numsFlux[WC_day_ts==0], pch=20, col=WC_cohortColors[WC_day_ts==0])
plot(WC_numsEaten[WC_day_ts==0], pch=20, col=WC_cohortColors[WC_day_ts==0])
plot(WC_numsDead[WC_day_ts==0], pch=20, col=WC_cohortColors[WC_day_ts==0])

index<-WC_cohorts==0 & WC_day_ts<1
plot(WC_numsFlux[index], pch=20, col=(WC_day_ts[index]+1)*2)

plot(WC_numsFlux, pch=20, col=colByNight)

index<-round(WC_day_ts)==WC_day_ts
plot(WC_numsFlux[index], pch=20, col=myGrey_trans)

plot(WC_numsDead[index], pch=20, col=myGrey_trans)
points(WC_numsDead[!index], pch=8, col=myGold)

plot(WC_numsEaten[index], pch=20, col=myGrey_trans)
points(WC_numsEaten[!index], pch=8, col=myGold)

############
oh_numsDead <- unlist(lapply(otherHab, getNumDead_Epi))
oh_numsEaten <- unlist(lapply(otherHab, getNumEaten))
oh_numsFlux <- unlist(lapply(otherHab, getTotalNumFlux))
oh_day_ts<- unlist(lapply(otherHab, getDayT))
oh_cohorts <- unlist(lapply(otherHab, getCohorts))
oh_totalNumbers <- unlist(lapply(otherHab, getTotalNumbers))


index<-round(oh_day_ts)==oh_day_ts
plot(oh_numsFlux[index], pch=20, col=myGrey_trans)
points(oh_numsDead[!index], pch=8, col=myGold)

plot(oh_numsDead[index], pch=20, col=myGrey_trans)
points(oh_numsDead[!index], pch=8, col=myGold)

plot(oh_numsEaten[index], pch=20, col=myGrey_trans)
points(oh_numsEaten[!index], pch=8, col=myGold)

numsecs <- 12*60*60
est_mL <- (WC_numsDead * numsecs)/WC_totalNumbers
hist(est_mL)

input_mL <- c(5e-5, 8.5e-4)
age_mat <- 3

#what is the Tcorr? to get it, est_mL = Tcorr * input_mL/numsecs
c=age_mat
thisNightIndex <- round(WC_day_ts)!=WC_day_ts & WC_cohorts==(c-1)
thisDayIndex <- round(WC_day_ts)==WC_day_ts & WC_cohorts==(c-1)
night_est_mL <- (WC_numsDead[thisNightIndex] * numsecs)/WC_totalNumbers[thisNightIndex]
this_dayt <- WC_day_ts[thisNightIndex]

day_est_mL <- (WC_numsDead[thisDayIndex] * numsecs)/WC_totalNumbers[thisDayIndex]
dayTcorr <- day_est_mL/input_mL[2]
hist(dayTcorr)


for(c in 1:10){
  this_mL <- 0.5*input_mL[2]
  if(c<age_mat){this_mL <- 0.5*input_mL[1]}
  thisNightIndex <- round(WC_day_ts)!=WC_day_ts & WC_cohorts==(c-1)
  thisDayIndex <- round(WC_day_ts)==WC_day_ts & WC_cohorts==(c-1)
  night_est_mL <- (WC_numsDead[thisNightIndex] * numsecs)/WC_totalNumbers[thisNightIndex]
  this_dayt <- WC_day_ts[thisNightIndex]
  
  day_est_mL <- (WC_numsDead[thisDayIndex] * numsecs)/WC_totalNumbers[thisDayIndex]
  test<-night_est_mL[this_dayt==2]
  
  par(mfrow=c(2,1))
  hist(night_est_mL)
  abline(v=this_mL, col=myRed_trans, lwd=3)
  mtext(paste("Age class ", c, sep=""), side=3, adj=0)
  hist(day_est_mL)
  abline(v=this_mL, col=myRed_trans, lwd=3)
  abline(v=mean(day_est_mL, na.rm=TRUE), col=myBlue_trans)
}

# storeNumsDead_nopref <- WC_numsDead
# storeNumbers_nopref <- WC_totalNumbers

# storeNumsDead_day <- WC_numsDead
# storeNumbers_day <- WC_totalNumbers
# 
# storeNumsDead_night <- WC_numsDead
# storeNumbers_night <- WC_totalNumbers
# 

ThisNC.nc <- nc_open(paste(thisPath, "output.nc", sep=""))
thisName <- "Hake"
thisData <- ncvar_get(ThisNC.nc, paste(thisName,10,"_Nums", sep=""))
thisNumbers <- apply(thisData, 3, sum)

# storeC10Numbers_noPref<-thisNumbers
# storeC10Numbers_day<-thisNumbers
# storeC10Numbers_night<-thisNumbers

thisYmin<-min(c(storeC10Numbers_day, storeC10Numbers_night, storeC10Numbers_noPref), na.rm=TRUE)
thisYmax<-max(c(storeC10Numbers_day, storeC10Numbers_night, storeC10Numbers_noPref), na.rm=TRUE)

par(mfrow=c(1,1))
plot(storeC10Numbers_noPref, type="l", ylim=c(thisYmin, thisYmax))
points(storeC10Numbers_day, type="l", col=myBlue, lwd=2)
points(storeC10Numbers_night, type="l", lty=2, col=myOrange)

# makeBlankPlot()
legend(legend=c("No pref", "Day", "Night"), col=c("black", myBlue, myOrange), lwd=c(1,2,1), lty=c(1,1,2), x="bottomright")

this_mL <- 8.5e-04
plot(est_mL, pch=20, col=colByNight, ylim=c(0, this_mL))
abline(h=this_mL, col="red")
abline(h=0.5*this_mL, col="red", lty=3)

oh_colByNight <- unlist(lapply(oh_day_ts, getColByNight))
oh_est_mL <- (oh_numsDead * numsecs)/oh_totalNumbers
plot(oh_est_mL, pch=20, col=oh_colByNight, ylim=c(0, this_mL))
abline(h=this_mL, col="red")
abline(h=0.5*this_mL, col="red", lty=3)


plot(oh_numsDead, pch=20)
plot(oh_totalNumbers, pch=20)
points(WC_totalNumbers, pch=8, col=myOrange)

index<-WC_cohorts==1 & WC_day_ts==0
oh_index<-oh_cohorts==1 & oh_day_ts==0

oh_totalNumbers[oh_index]
WC_totalNumbers[index]

oh_numsDead[oh_index]
WC_numsDead[index]

test<- tapply(oh_numsDead, list(oh_cohorts, oh_day_ts), nonZeroMean)

## check numbers tracer from nc file
ThisNC.nc <- nc_open(paste(thisPath, "output.nc", sep=""))
thisTracer <- "Hake1_Nums"
thisData <- ncvar_get(ThisNC.nc, thisTracer)

nboxes <-dim(thisData)[2]
thisEpiData <- rep(NA, nboxes)
for(b in 1:nboxes){
  x <- thisData[,b,1]; xx <- x>0
  y <- x[xx]
  if(length(y)>0){
    thisEpiData[b]<-y[1]
  }
}
