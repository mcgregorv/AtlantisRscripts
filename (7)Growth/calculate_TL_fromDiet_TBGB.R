##calculate TL from diet for comparisons
## diet summarised in SummaryDiets_AtlantisOutput.R in R\\Growth\\

this_out<-"TestZG"; runFolder="TBGB_JP2"
thisDesc <- paste(runFolder, this_out,sep="")

this_path = paste(DIR$'Base', "TBGB\\",runFolder,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")

plotPath<-paste(DIR$'Base',"TBGB\\Figures\\Testing\\FIshing\\",sep="")


preyGroupsDF<-read.csv(paste(this_path,"TBGB_groups.csv",sep=""))


dietMat <- read.csv(paste0(outPath,"realiseDietSnapshot1880_1960.csv"))
row.names(dietMat) <- dietMat$X
colnames(dietMat)[1] <- "Code"
dietMat$Code<-as.character(dietMat$Code)

preyGroupsDF$Predator <- preyGroupsDF$IsPredator
predators<-as.character(preyGroupsDF$Code[preyGroupsDF$Predator==1])
prey<-as.character(preyGroupsDF$Code[preyGroupsDF$Predator==0])
npreds<-sum(preyGroupsDF$Predator)

preyGroupsDF$Code<-as.character(preyGroupsDF$Code)

#Assign TL 1 to primary producers 
preyGroupsDF$TL <- ifelse(preyGroupsDF$Predator==1,0,0.01)

#Assign TL 2 to macrozoobenthos (ZM) ## ZM is mesozoo - and why '2'?
# preyGroupsDF[preyGroupsDF$Code=="ZM",'TL'] <- 1
preyGroupsDF[preyGroupsDF$Code=="DF","TL"]<-0.01
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
dietMatrix$sums <- rowSums(dplyr::select(dietMatrix,-Code),na.rm = T)
dietMatrix <- data.frame(dplyr::select(dietMatrix,Code),apply(dplyr::select(dietMatrix,-Code), 2, function(x) x/dietMatrix$sums))
dietMatrix <- dplyr::select(dietMatrix, -sums)
row.names(dietMatrix) <- dietMatrix$Code



###-----Function to calculate   
trophic.levels <- function() {
  success <- FALSE
  r=0
  while (!success) {    # do something over and over until success
    
    #---Order diet by TL
    newdietMat <- data.frame(t(dplyr::select(dietMatrix, -Code)))
    newdietMat$Code <- rownames(newdietMat)
    newdietMat <- left_join(newdietMat, preyGroupsDF[,c("Code","TL")], by = "Code")
    
    ##Create a matrix layer to calculate proportions
    ##--Add the updated TL's each time its run
    layer1 <- data.frame(apply(newdietMat %>% dplyr::select(-Code), 2, function(x) newdietMat$TL*x))
    TL1 <- data.frame(cbind(apply(layer1[,c(1:npreds)], 2, function(x) sum(x, na.rm = T)+1),
                            colSums(!is.na(layer1[,1:npreds]))-colSums(layer1[,1:npreds]>0, na.rm = T)))
    TL1$Code <- row.names(TL1)
    TL1 <- TL1[TL1$X2==0,] 
    # TL1[TL1$Code %in% c("ZM"),1]<-2
    # TL1[TL1$Code %in% c("ZL", "ZG"),1]<-2.5
    TL1[TL1$Code %in% c("DF"),1]<-0.01
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
preyGroupsDF <- trophic.levels()

preyGroupsDF$Code[preyGroupsDF$TL==0]

outGroups<-preyGroupsDF[,c("Code", "Name", "LongName", "NumCohorts", "IsPredator", "TL")]

write.csv(outGroups, paste0(outPath,"TBGBGroupsTL.csv"))
