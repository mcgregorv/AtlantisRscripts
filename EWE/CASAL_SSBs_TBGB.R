dataPath<-paste(DIR$'Base',"Biology\\CASAL_TBGB\\", sep="")
TBGBgroups<-read.csv(paste(dataPath,"TBGB_groups.csv", sep=""))

fnames<-list.files(dataPath); nfiles<-length(fnames)
years<-seq(1900,2012); nyears<-length(years)

storeSSBs<-array(NA, dim=c(nfiles, nyears)); storeSpeciesCodes<-rep(NA, nfiles)
storeCatches<-0*storeSSBs
for(f in 1:nfiles){
  thisFile<-paste(dataPath, fnames[f],"\\MPD.txt", sep="")
  if(file.exists(thisFile)){
    thisLines<-readLines(thisFile)
    x<-grep("SSB", thisLines)
    y<-get_first_number(thisLines[x], n="all")
    storeSSBs[f,]<-y
    x<-unlist(str_split(fnames[f]," "))[2]; storeSpeciesCodes[f]<-x
    ## grab the catches too
    x<-grep("^catch", thisLines); xx<-get_first_number(thisLines[x], n="all")
    storeCatches[f,]<-xx
  }
}

storeSSBs<-storeSSBs[!is.na(storeSpeciesCodes),]; storeSpeciesCodes<-storeSpeciesCodes[!is.na(storeSpeciesCodes)]
outSSBs<-data.frame(storeSSBs); colnames(outSSBs)<-years
write.csv(outSSBs, file=paste(dataPath,"AllCASALSSBs_TBGB.csv", sep=""), row.names=storeSpeciesCodes)

thisCode<-"BAR"; 
pdf(paste(dataPath,"TBGB_SSBs.pdf", sep=""))
par(mfrow=c(3,2), lend=1, mar=c(4,4,1,4))
for(f in 1:nfiles){
  thisCode<-storeSpeciesCodes[f]
  thisIndex<-grep(thisCode, fnames)
  thisSSB<-storeSSBs[thisIndex,]
  if(sum(thisSSB, na.rm=TRUE)){
    thisCatches<-storeCatches[f,]; thisCatchLabs<-pretty(thisCatches)
    plot(thisSSB, type="l", ylim=c(0,max(thisSSB)), xlab="Time (years from 1900)", ylab="SSB (tonnes)")
    mtext(thisCode,side=3, adj=0)
    par(new=TRUE)
    plot(thisCatches, type="h", lwd=5, col=myGrey_trans, xaxt="n",ylab="", xlab="", yaxt="n")
    axis(at=thisCatchLabs, labels = thisCatchLabs, side=4, col=myGrey)
    mtext("Catches (tonnes)", side=4, adj=0.5, line=2.5, col=myGrey, cex=0.8)
  }
}
dev.off()


#read in tracers to get model area volume
TBGB_NC.nc<-nc_open(paste(dataPath, "TBGB_output.nc", sep=""))
TBGBvol<-ncvar_get(TBGB_NC.nc,"volume")
# depth layer 6 is 1 m deep, so used for area
# box 1 is boundary, so taken out of model area
# timestep 1 is the initial conditions (not that it matters much which timestep used)
vol_m2<-sum(TBGBvol[6,-1,1], na.rm=TRUE); vol_km2<-vol_m2 * 1e-6 

thisB0_km2<-thisSSB[1]/vol_km2

BAR_tkm2<-0.64
BAR_t<-BAR_tkm2 * vol_km2
# 3149.58



