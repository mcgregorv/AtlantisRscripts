dataPath<-paste(DIR$'Base',"Biology\\CASAL\\", sep="")

fnames<-list.files(dataPath); nfiles<-length(fnames)
years<-seq(1900,2014); nyears<-length(years)

storeSSBs<-array(NA, dim=c(nfiles, nyears)); storeSpeciesCodes<-rep(NA, nfiles)
for(f in 1:nfiles){
  thisFile<-paste(dataPath, fnames[f],"\\MPD.txt", sep="")
  if(file.exists(thisFile)){
    thisLines<-readLines(thisFile)
    x<-grep("SSB", thisLines)
    y<-get_first_number(thisLines[x], n="all")
    storeSSBs[f,]<-y
    x<-unlist(str_split(fnames[f]," "))[2]; storeSpeciesCodes[f]<-x
  }
}

storeSSBs<-storeSSBs[!is.na(storeSpeciesCodes),]; storeSpeciesCodes<-storeSpeciesCodes[!is.na(storeSpeciesCodes)]
outSSBs<-data.frame(storeSSBs); colnames(outSSBs)<-years
write.csv(outSSBs, file=paste(dataPath,"AllCASALSSBs.csv", sep=""), row.names=storeSpeciesCodes)
