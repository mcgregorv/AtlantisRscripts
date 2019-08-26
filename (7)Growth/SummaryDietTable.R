## reads in diets and creates latex table for appendix in paper - at the species group level rather than prey type
year0<-1865
dietYears<-seq(1900,2016) #just replace with 1970 if only want one year, or give a range of years - will give average over them
dietTimeSteps<-dietYears-year0+1 #to index dietYears

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
inputsPath<-paste(this_path,"inputs\\",sep="")
# plotPath<-paste(DIR$'Reports',"Figures\\",sep="") #overwrites the one in the report
plotPath<-paste(DIR$'Figures',"Mortality\\test", sep="")
outPath<-paste(DIR$'Reports',"CRAMpaper1_PeerJ\\DataAndCode\\Data\\",sep="")
tablePath <- paste(outPath,"..\\..\\table_", sep="")

preyGroupsDF<-read.csv(paste(this_path,"inputs\\supporting\\preyGroups.csv",sep=""))
preyGroups<-sort(unique(preyGroupsDF$PreyGroup)); nPreyGroups<-length(preyGroups)

## first plot intended diets
groupsDF<-read.csv(paste(this_path,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
groupsDFPaper<-read.csv(paste(this_path,"CRAM_groupsPaper.csv", sep=""))

predators<-as.character(groupsDF$Code[groupsDF$IsPredator==1]); npreds<-length(predators)

realisedDiets<-read.csv(paste(outPath,"Historic_withFishing_outputDietCheck.txt", sep=""), sep=" ")
realisedDiets$year_ts<-realisedDiets$Time/365
yearIndex<-realisedDiets$year_ts %in% dietTimeSteps
yearRelDiets<-realisedDiets[yearIndex, ]

summaryWRTtime <- data.frame(matrix( NA, nrow=ng, ncol=ng))
colnames(summaryWRTtime)<- colnames(yearRelDiets)[6:(dim(yearRelDiets)[2]-1)]
for(p in 1:ng){
  thisCode<-groupsDFPaper$Code[p]
  if(thisCode %in% predators){
    thisDiets <- yearRelDiets[yearRelDiets$Predator==as.character(thisCode),]
    dietSums <- apply(thisDiets[,c(6:(dim(yearRelDiets)[2]-1))],2, sum, na.rm=TRUE)
    dietProps <- round(dietSums/sum(dietSums),2)
  } else{
    dietProps <- rep(0,ng)
  }
  thisName <- gsub("_"," ",groupsDFPaper$Name[groupsDFPaper$Code==thisCode])
  rownames(summaryWRTtime)[p]<- thisName
  summaryWRTtime[p,]<-dietProps
}

## use numbers for prey, and make sure predators are in the same order (and include all groups)
# writeTable <-summaryWRTtime

preyNames <- gsub("_", " ", groupsDFPaper$Name[match(colnames(summaryWRTtime), groupsDFPaper$Code)])
colnames(summaryWRTtime)<- preyNames

## write as .csv file for supp. material
write.csv(summaryWRTtime, paste(DIR$'Reports',"CRAMpaper1_PeerJ\\DataAndCode\\Data\\Historic_withFishing_outputDietSummary.csv",sep=""))


## need to do in 2 parts
colIndex<-1:28
## write table to latex
latexFile <- paste(tablePath,"dietSummaries1.txt", sep="")
tableHead <- "\\fontsize{8}{8}\\selectfont
\\begin{longtable}{lrrrrrrrrrrrrrrrrrrrrrrrrrrrr}	
	\\centering
"
cat(tableHead, file=latexFile, append=FALSE)

temp<-paste(colnames(summaryWRTtime)[colIndex],collapse="} & \\textbf{")
# columnHeadings <-paste("\\textbf{Predator} &	 \\textbf{", temp, "}\\\\", collapse="", sep="")
columnHeadings<-paste((1:ng)[colIndex],collapse="} & \\textbf{")
cat("\\textbf{Predator} & \\textbf{", file=latexFile, append=TRUE)
cat(columnHeadings, file=latexFile, append=TRUE)
cat("}\\\\
    ", file=latexFile, append=TRUE)
for(p in 1:ng){
  thisLine<-paste(summaryWRTtime[p,colIndex],collapse=" & ")
  cat(paste(rownames(summaryWRTtime)[colIndex][p]," & ", thisLine, "\\\\
            ", sep=""), file=latexFile, append=TRUE)
}
tableEnd <- "\\end{longtable} 
\\normalsize"
cat(tableEnd, file=latexFile, append=TRUE)

