##gets R0 and B0 from CASAL MPD files

casalPath<-paste(DIR$'Base',"Biology\\CASAL\\",sep="")
modelPath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
groupsDF<-read.csv(paste(modelPath,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

groupsDF$B0<-NA; groupsDF$R0<-NA

casalFileNames<-list.files(casalPath)

for(g in 1:ng){
  thisCode<-groupsDF$Code[g]
  x<-grep(thisCode,casalFileNames)
  thisFileName<-casalFileNames[x]
  if(length(thisFileName)==1){
    cat(as.character(thisCode),"--")
    if(file.exists(paste(casalPath,thisFileName,"\\MPD.txt",sep=""))){
      thisCasalLines<-readLines(paste(casalPath,thisFileName,"\\MPD.txt",sep=""))
      x<-grep("initialization.B0",thisCasalLines)[1]
      thisLine<-thisCasalLines[x]
      thisB0<-get_first_number(thisLine)
      groupsDF$B0[g]<-thisB0
      ##R0
      x<-grep("* R0",thisCasalLines)[1]
      thisLine<-thisCasalLines[x+1]
      thisR0<-get_first_number(thisLine)
      groupsDF$R0[g]<-thisR0
      
    }
  }
}

#write it out so can be used again
write.csv(groupsDF[,c("Code","B0","R0")],paste(modelPath,"CRAM_B0.csv",sep=""),row.names = FALSE)
