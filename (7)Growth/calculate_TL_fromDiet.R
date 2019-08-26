##calculate TL from diet for comparisons
## diet summarised in SummaryDiets_AtlantisOutput.R in EWE\diets\

this_path<-paste(DIR$'Base',"ATLANTISmodels\\base\\EWEbase\\", sep="")
outPath<-this_path

preyGroupsDF<-read.csv(paste(this_path,"CRAM_groups.csv",sep=""))

dietMat <- read.csv(paste0(outPath,"realiseDietSnapshot1900_2015.csv"))
row.names(dietMat) <- dietMat$X
colnames(dietMat)[1] <- "Code"
dietMat$Code<-as.character(dietMat$Code)

preyGroupsDF$Predator <- preyGroupsDF$IsPredator
predators<-as.character(preyGroupsDF$Code[preyGroupsDF$Predator==1])
prey<-as.character(preyGroupsDF$Code[preyGroupsDF$Predator==0])
npreds<-sum(preyGroupsDF$Predator)

preyGroupsDF$Code<-as.character(preyGroupsDF$Code)

#Assign TL 1 to primary producers 
  preyGroupsDF$TL <- ifelse(preyGroupsDF$Predator==1,0,1)
  ## not all of these are primary producers - includes bacteria and detritus
  ## read in pre-defined trophic levels and use values from here for those that are not predators
  si_tls<-read.csv(paste(this_path,"CRAM_trophicLevels_isotopes.csv", sep=""))
  preyGroupsDF$TL[preyGroupsDF$Predator==0]<-si_tls$TrophicLevel2[ match(preyGroupsDF$Code[preyGroupsDF$Predator==0], si_tls$Code)]
  ## oh, doesn't get off the ground
  ## set 0's to 0.1..?
  preyGroupsDF$TL[preyGroupsDF$Predator==0 & preyGroupsDF$TL==0]<-0.1
  
#Assign TL 2 to macrozoobenthos (ZM) ## ZM is mesozoo - and why '2'?
  # preyGroupsDF[preyGroupsDF$Code=="ZM",'TL'] <- 1
##---select just predators
  dietMatrix <- as.matrix(dietMat[,predators])  
###---Remove cannibalism and reallocate proportions
  diag(dietMatrix) <- NA
##----Reallocate diets <0.5% of total. Probably need to find a better way than this....
  dietMatrix[dietMatrix<0.05] <- NA
##--Make matrix back into data frame
  dietMatrix <- data.frame((dietMatrix)) 
  dietMatrix$Code <- row.names(dietMatrix)
  dietMatrix <- left_join (dietMatrix, dietMat[,c("Code",prey)],by = "Code")
##---Sum the diets and rescale to 1
  dietMatrix$sums <- rowSums(select(dietMatrix,-Code),na.rm = T)
  dietMatrix <- data.frame(select(dietMatrix,Code),apply(select(dietMatrix,-Code), 2, function(x) x/dietMatrix$sums))
  dietMatrix <- select(dietMatrix, -sums)
  row.names(dietMatrix) <- dietMatrix$Code
  


###-----Function to calculate   
  trophic.levels <- function() {
    success <- FALSE
    r=0
    while (!success) {    # do something over and over until success
      
      #---Order diet by TL
      newdietMat <- data.frame(t(select(dietMatrix, -Code)))
      newdietMat$Code <- rownames(newdietMat)
      newdietMat <- left_join(newdietMat, preyGroupsDF[,c("Code","TL")], by = "Code")
      
      ##Create a matrix layer to calculate proportions
      ##--Add the updated TL's each time its run
      layer1 <- data.frame(apply(newdietMat %>% select(-Code), 2, function(x) newdietMat$TL*x))
      TL1 <- data.frame(cbind(apply(layer1[,c(1:npreds)], 2, function(x) sum(x, na.rm = T)+1),
                              colSums(!is.na(layer1[,1:npreds]))-colSums(layer1[,1:npreds]>0, na.rm = T)))
      TL1$Code <- row.names(TL1)
      TL1 <- TL1[TL1$X2==0,] 
      # TL1[TL1$Code %in% c("ZM"),1]<-2
      # TL1[TL1$Code %in% c("ZL", "ZG"),1]<-2.5
      TL1[TL1$Code %in% c("DF"),1]<-1
      #this means it only keeps the ones that were calculated based on all non-zero (or NA) entries in the diet matrix 
      # - ie, those for which TLs already exist for all of their prey
      if("ZL" %in% TL1$Code){cat(r," has ZL TL ", TL1[TL1$Code=="ZL",1], "\n")}
      
      ##Add the newly calculated TL's to the main data
      for(i in 1:nrow(TL1)){
        preyGroupsDF[preyGroupsDF$Code==TL1[i,c("Code")],'TL'] <- TL1[i,1]       
      }    
      # check for success - if true then stop loop
      success <- sum(preyGroupsDF$TL == 0)==0
      if(r>200){success=1}
      r=r+1
      cat(paste("run ", r, " has ", sum(preyGroupsDF$TL == 0), " to go\n"))
    }
    return(preyGroupsDF)
  }

  testCode<-"ZM"
  test<-newdietMat[!is.na(newdietMat[,c(testCode)]),c(testCode, "Code", "TL")]
  test
  sum(test[,1] * test[,3])  +1
  
preyGroupsDF <- trophic.levels()

preyGroupsDF$Code[preyGroupsDF$TL==0]

outGroups<-preyGroupsDF[,c("Code", "Name", "LongName", "NumCohorts", "IsPredator", "TL")]

write.csv(outGroups, paste0(outPath,"CRAMGroupsTL.csv"))

## use those from stable isotopes and compare
si_codes<- si_tls$Code[!is.na(si_tls$Isotope)]

compareTLs<-data.frame(matrix(NA, ncol=3, nrow=length(si_codes))); colnames(compareTLs)<-c("Code", "TLcalc", "TLiso")
compareTLs$Code<-si_codes; compareTLs$TLiso<-si_tls$Isotope[!is.na(si_tls$Isotope)]
compareTLs$TLcalc<-preyGroupsDF$TL[match(si_codes, preyGroupsDF$Code)]

plot(x=compareTLs$TLiso, y=compareTLs$TLcalc, ylim=c(2,5.2), xlim=c(2,5.5))

compareTLs<-compareTLs[rev(order(compareTLs$TLiso)),]
si_codes<-compareTLs$Code

par(lend=1)
plot(x=seq(1,length(si_codes)), y=rep(0,length(si_codes)), ylim=c(0,5.5), type="n", xaxt="n", xlab="", ylab="TL")
for(g in 1:length(si_codes)){
  points(x=g-0.1, y=compareTLs$TLiso[g], type="h", lwd=3, col=myOrange)
  points(x=g+0.1, y=compareTLs$TLcalc[g], type="h", lwd=3, col="midnightblue")
}
par(las=2)
axis(at=seq(1,length(si_codes)), labels = si_codes, side=1)

## compare to prev estimates where isotopes not avail
compare2<-outGroups
compare2$prevGuess<-si_tls$TrophicLevel2[match(si_tls$Code, compare2$Code)]
compare2$prevGuessJuv<-si_tls$TrophicLevel1[match(si_tls$Code, compare2$Code)]


par(lend=1)
plot(x=seq(1,dim(compare2)[1]), y=rep(0,dim(compare2)[1]), ylim=c(0,5.5), type="n", xaxt="n", xlab="", ylab="TL")
for(g in 1:dim(compare2)[1]){
  points(x=g-0.1, y=compare2$prevGuess[g], type="h", lwd=3, col=myOrange)
  points(x=g+0.1, y=compare2$TL[g], type="h", lwd=3, col="midnightblue")
}
par(las=2)
axis(at=seq(1,dim(compare2)[1]), labels = compare2$Code, side=1)



par(lend=1)
plot(x=seq(1,dim(compare2)[1]), y=rep(0,dim(compare2)[1]), ylim=c(0,5.5), type="n", xaxt="n", xlab="", ylab="TL")
for(g in 1:dim(compare2)[1]){
  points(x=g-0.1, y=compare2$prevGuessJuv[g], type="h", lwd=3, col=myOrange)
  points(x=g+0.1, y=compare2$TL[g], type="h", lwd=3, col="midnightblue")
}
par(las=2)
axis(at=seq(1,dim(compare2)[1]), labels = compare2$Code, side=1)

plot(x=compare2$prevGuess, y=compare2$TL, pch=20)
points(x=compare2$prevGuessJuv, y=compare2$TL, pch=4, col=myOrange)

