#read in sampled availabilities and create pPREY lines for biol.prm file
this_sampledAvails<-read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ExpectedValues.csv",sep=""))
# this_sampledAvails<-read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ExpectedValues_scaledDown.csv",sep=""))
# this_sampledAvails<-read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ExpectedValues_scaledUp.csv",sep=""))

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

ppGroups<-as.character(read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ppGroups.csv",sep=""))[,1]); npg<-length(ppGroups)

groupsDF<-read.csv(paste(this_path,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
ndg<-ng+3 #number of diet groups, include sediment version of DL, DR, DC

pPREYfile<-paste(DIR$'Base',"ATLANTISmodels\\inputs\\biol_prm\\pPREY\\sampledAvails.txt",sep="")
cat("## sampled availabilities take 1\n",file=pPREYfile,append=FALSE)

adIndex<-grep("ad",ppGroups); juvIndex<-grep("juv",ppGroups)
naIndex<-grep("ad|juv",ppGroups,invert = TRUE)
not_adultIndex<-grep("ad",ppGroups,invert = TRUE)

for(g in 1:npg){
  thisVar<-ppGroups[g]
  thisCode<-gsub('juv|ad',"",thisVar)
  thisAge<-gsub(thisCode,"",thisVar)
  thisVec<-this_sampledAvails[g,]
  if(thisAge==""){
    ppreyvar<-paste("pPREY",thisCode,"\t",ndg,"\n",sep="")
    cat(ppreyvar,file=pPREYfile,append=TRUE)
    temp<-c(signif(as.double(thisVec[not_adultIndex]),4),c(0,0,0))
    newLine<-paste(temp,collapse = " ")
    cat(newLine,file=pPREYfile,append=TRUE)
    cat("\n",file=pPREYfile,append=TRUE)
  } else{
    if(thisAge=="ad"){
      ageNumber<-2
    } else{
      ageNumber<-1
    }
    #first do not adult prey
    ppreyvar<-paste("pPREY",1,thisCode,ageNumber,"\t",ndg,"\n",sep="")
    cat(ppreyvar,file=pPREYfile,append=TRUE)
    temp<-c(signif(as.double(thisVec[not_adultIndex]),4),c(0,0,0))
    newLine<-paste(temp,collapse = " ")
    cat(newLine,file=pPREYfile,append=TRUE)
    cat("\n",file=pPREYfile,append=TRUE)
    #first do adult prey
    ppreyvar<-paste("pPREY",2,thisCode,ageNumber,"\t",ndg,"\n",sep="")
    cat(ppreyvar,file=pPREYfile,append=TRUE)
    temp<-c(signif(as.double(thisVec[adIndex]),4))
    filledTemp<-rep(0,ng); filledTemp[groupsDF$NumCohorts>1]<-temp; 
    temp<-c(filledTemp,c(0,0,0))
    temp[is.na(temp)]<-0
    newLine<-paste(temp,collapse = " ")
    cat(newLine,file=pPREYfile,append=TRUE)
    cat("\n",file=pPREYfile,append=TRUE)
  }
}
