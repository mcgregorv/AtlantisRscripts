## catch histories were created in setUpScenario_ts_files_versionC1.R - they are full historic + 50 year future runs
## biol.prm files for eddy_scale are  created first

thisPath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

version<-"C"
runPath<-paste(thisPath,"eddySens\\",sep="")

## set up DF with runs specified
## set up DF with runs specified
## set up DF with runs specified
scenarios<-c("All0catch", "All50catch", "All80catch",  "All100catch",  "All120catch",  "All150catch") ## C
nscenarios<-length(scenarios)
eddysens<-c("EddyNOchange", "Eddy")
neddys<-length(eddysens); eddyIndex<-seq(1,neddys)
nruns<-nscenarios * neddys

lookup_df<-data.frame(array(NA, dim=c(nruns,3)))
colnames(lookup_df)<-c("Run", "Scenario", "Eddy")
lookup_df$Run<-seq(1,nruns)
lookup_df$Scenario<-rep( sort(rep(scenarios, neddys)))
lookup_df$Eddy<-rep(eddyIndex, (nscenarios))

lookup_df<-lookup_df

nruns<-dim(lookup_df)[1]

baseModelPath<-paste(thisPath,"base\\outputBase\\",sep="")
ThisNC.nc<-nc_open(paste(baseModelPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
nts<-dim(thisVol)[3]; nboxes<-dim(thisVol)[2]; nlayers<-dim(thisVol)[1]

nts<-202
allTracers<-names(ThisNC.nc$var); 
expTracers<-allTracers[grep("Nums|StructN|ResN", allTracers, invert = TRUE)] # take out numbers and individual weights
ntracers<-length(expTracers)

storeTracers<-array(NA, dim=c(nruns, ntracers,  nlayers, nboxes, nts))
nts_list<-rep(NA, nruns)

thisRuns<-1:12
# for(r in 1:nruns){
for(r in thisRuns){
  cat(r,"--")
  thisRunPath<-paste(runPath,"outputEDDY",version,"_",r,"\\",sep="")
  runNC.nc<-nc_open(paste(thisRunPath, "output.nc",sep=""))
  runVol<-ncvar_get(runNC.nc,"volume"); run_nts<-dim(runVol)[3]
  nts_list[r]<-run_nts
}
nts_list

storeBiomass<-array(NA, dim=c(nruns, ntracers, nts))

# for(r in 1:nruns){
for(r in thisRuns){
  cat(r,"--")
  thisRunPath<-paste(runPath,"outputEDDY",version,"_",r,"\\",sep="")
  runNC.nc<-nc_open(paste(thisRunPath, "output.nc",sep=""))
  runVol<-ncvar_get(runNC.nc,"volume"); run_nts<-dim(runVol)[3]
  for(t in 1:ntracers){
    thisTracer<-expTracers[t]
    thisData<-ncvar_get(runNC.nc,thisTracer)
    test<-grep("_N", thisTracer)
    if(length(dim(thisData))==3){
      storeTracers[r,t,,,1:min(run_nts, nts)]<-thisData[,,1:(min(run_nts, nts))]
    } else if(length(dim(thisData))==2 & dim(thisData)[1]==nboxes){
      storeTracers[r,t,nlayers,,1:min(run_nts,nts)]<-thisData[,1:(min(run_nts, nts))]
    }
    ##do biomasses
    if(length(test)>0){
      if(length(dim(thisData))==3){
        yy<-apply(thisData[-nlayers,,]* runVol[-nlayers,,], 3, sum) *mg_2_tonne * X_CN 
       } else if(length(dim(thisData))==2 & dim(thisData)[1]==nboxes){
         yy<-apply(thisData* runVol[nlayers,,], 2, sum) *mg_2_tonne * X_CN 
         yy<-rep(NA, length(yy)) ## just analyse water column biomasses
       }
      storeBiomass[r,t,1:min(run_nts, nts)]<- yy[1:(min(run_nts, nts))]
    }
  }
}

## prob only need the future period (can cut off the historic)
## and take 2015 as the base, and compare all to this
# useTracers<-storeTracers[,,,,150:dim(storeTracers)[5]]
# 
# baseTracers<-useTracers[,,,,1]
# relTracers<-0*useTracers
# 
# for(i in 1:dim(useTracers)[5]){
#   relTracers[,,,,i]<-useTracers[,,,,i]/baseTracers
# }
# 
# cleanNAs<-function(x){
#   y<-x
#   if(length(x)==0){y<-NA}
#   else if(is.na(x)){y<-NA}
#   else if(x==""){y<-NA}
#   else if(x=="NaN"){y<-NA}
#   else if(x=="Inf"){y<-NA}
#   return(y)
# }
# relTracers<-apply(relTracers,seq(1,5), cleanNAs)


# storeTracers<-useTracers
####################################
#write out the tracers and rels so can read in elsewhere
save(list=c( "nts_list","storeBiomass", "storeTracers","lookup_df", "eddysens"), file=paste(thisPath,"eddySens\\modelTracers_V",version,sep=""))
####################################

####################################

test_ts<-10
testTracerIndex<-grep("_N", expTracers); testTracers<-expTracers[testTracerIndex]
test<-apply(relTracers[,testTracerIndex,,,test_ts],1, mean, na.rm=TRUE)

timeIndex<-seq(1,53)
## take ZM
testTracerIndex<-grep("^Zoo_N",expTracers)
testTracerIndex<-grep("Zoo_N",expTracers); expTracers[testTracerIndex]

testData<-relTracers[,testTracerIndex,,,timeIndex]

testMean<-apply(testData,c(1,2), mean, na.rm=TRUE); ymax<-max(testMean, na.rm=TRUE)
colByZ<-c(myGold, myGreen,myAqua,myBlue)
plot(testMean[1:3,1], ylim=c(0,1), type="n", xlab="Run", ylab="Biomass relative to base")
for(z in 1:4){
  points(testMean[3:4,z], pch=20, col=colByZ[z],cex=1.5)
}
legend(legend=expTracers[testTracerIndex], pch=20, pt.cex=1.5, col=colByZ, x="bottomright", bty="n")

ppcodes<-c("DF", "MA", "MB", "PL", "PS")
primaryProducers<-c(paste(groupsDF$Name[groupsDF$Code %in% ppcodes],"_N", sep=""), paste(groupsDF$Name[groupsDF$Code %in% ppcodes],"_S", sep=""))
testTracerIndex<-grep(paste(primaryProducers,collapse="|", sep=""), expTracers)
expTracers[testTracerIndex]


testData<-relTracers[,testTracerIndex,,,timeIndex]

testMean<-apply(testData,c(1,2), mean, na.rm=TRUE)
colByZ<-colorRampPalette(colors=c(myGold, myGreen,myAqua,myBlue,myRed))(length(testTracerIndex))
plot(testMean[1:5,1], ylim=c(0,1.2), type="n", xlab="Run", ylab="Biomass relative to base")
for(z in 1:length(testTracerIndex)){
  points(testMean[1:5,z], pch=20, col=colByZ[z],cex=1.5)
}
legend(legend=expTracers[testTracerIndex], pch=20, pt.cex=1.5, col=colByZ, x="bottomright", bty="n")

## read in trophic levels so can plot by these
trophicLevels<-read.csv(paste(thisPath, "base\\EWEbase\\CRAMGroupsTL.csv", sep=""))

pdf(paste(plotPath,"relBiomByTL.pdf", sep=""), height=10, width=10)
par(mfrow=c(5,2), mar=c(4,4,1,1))
for(l in 1:5){
  levelNames<-trophicLevels$Name[trunc(trophicLevels$TL)==l]
  testTracerIndex<-grep(paste(levelNames,collapse="|", sep=""), expTracers)
  expTracers[testTracerIndex]
  
  testData<-relTracers[,testTracerIndex,,,timeIndex]
  
  testMean<-apply(testData,c(1,2), mean, na.rm=TRUE)
  thisMax<-max(testMean[1:5,], na.rm=TRUE)*1.05
  thisMax<-1.2
  colByZ<-colorRampPalette(colors=c(myGold, myGreen,myAqua,myBlue,myRed))(length(testTracerIndex))
  plot(testMean[1:5,1], ylim=c(0,thisMax), type="n", xlab="Run", ylab="Biomass relative to base")
  for(z in 1:length(testTracerIndex)){
    points(testMean[1:5,z], pch=20, col=colByZ[z],cex=1.5)
  }
  makeBlankPlot()
  legend(legend=expTracers[testTracerIndex], pch=20, pt.cex=1.5, col=colByZ, x="center", bty="n", ncol=3, cex=0.8)
}
dev.off()

i<-grep("Macroalgae_N", expTracers)
testMA<-relTracers[1:5,i,,,]
testMA1<-apply(testMA, c(1,4), mean, na.rm=TRUE)
for(j in 1:5){
plot(testMA1[j,])
}
test2<-apply(testMA[1,,,],2,mean, na.rm=TRUE)

test2<-storeTracers[3,i,,5,]
test2<-apply()

i<-grep("Meio", expTracers)
test<-storeTracers[1,i,,,]
testBase<-storeTracers[5,i,,,]

reltest<-relTracers[1,i,,,36]






