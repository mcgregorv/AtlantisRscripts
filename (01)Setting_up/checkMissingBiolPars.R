biolFile <- "C:\\Projects\\2017\\ATLANTISChathamRiseRepository\\TBGB\\TBGB_JP2\\TBGB_biol.prm"
biolLines <- readLines(biolFile)

groupsDF<-read.csv("C:\\Projects\\2017\\ATLANTISChathamRiseRepository\\TBGB\\TBGB_JP2\\TBGB_groups.csv")
codes <- as.character(groupsDF$Code)

testLines <- biolLines[grep("^#", biolLines, invert = TRUE)]
codesCheck <- paste(codes, collapse="|")
testLines <- testLines[grep(codesCheck, testLines)]

missingPars <-c()
nlines<-length(testLines)
for(i in 1:nlines){
  thisLine <- unlist(str_split(testLines[i], " |\t"))[1]
  thisPattern <- gsub(codesCheck,"XXX", thisLine)
  for(g in 1:ng){
    testPar <- gsub("XXX", as.character(groupsDF$Code[g]), thisPattern)
    test4par <- grep(testPar, biolLines)
    if(length(test4par)==0){
      missingPars <- c(missingPars, testPar)
    }
  }
}


