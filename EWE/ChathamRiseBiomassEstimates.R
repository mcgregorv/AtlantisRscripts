## read in the atlantis estimated biomasses
# also grab the trawl survey absolute biomass estiamtes (from CASAL file), the trawl survey rel estiamtes, and the estimated B0
## approximate q
## combine with atlantis B(1900) and perhaps B(1970)

load(file=paste(DIR$'Base',"\\ATLANTISmodels\\base\\EWEbase\\CRAM_biomasses", sep=""))

groupsDF<-read.csv(paste(DIR$'Base',"\\ATLANTISmodels\\CRAM_groups.csv", sep=""))

dataPath<-paste(DIR$'Base',"Biology\\CASAL\\", sep="")

fnames<-list.files(dataPath); nfiles<-length(fnames)
years<-seq(1900,2014); nyears<-length(years)

storeB0s<-rep(NA, nfiles); storeSpeciesCodes<-rep(NA, nfiles)
storeAbsTS <- rep(NA, nfiles)
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


