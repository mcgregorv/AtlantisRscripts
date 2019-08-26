## read in bootstrap tracers
## calc CVs based on biomass tracers
## rank groups by CV

Version<-"B"

dataPath<-paste(DIR$'Base',"ATLANTISmodels\\BootstrapROMS\\", sep="")

load(paste(dataPath,"modelOut",Version,"_modelTracers", sep=""))
## storeTracers and storeTracersByCell are already converted to tonnes

#get _N tracers
outPath<-paste(DIR$'Base',"ATLANTISmodels\\base\\ouputROMS",Version,"\\",sep="")
BaseNC.nc<-nc_open(paste(outPath, "outputROMSBootstrap",Version,"1\\output.nc",sep=""))
allTracers<-names(BaseNC.nc$var)
thisVol<-ncvar_get(BaseNC.nc, "volume"); nlayers<-dim(thisVol)[1]; nts<-dim(thisVol)[3]
x<-grep("_N",allTracers); temp<-allTracers[x]
y<-grep("Nums",temp, invert = TRUE); Ntracers<-temp[y]; ntracers<-length(Ntracers)


calcCV<-function(x){
  thisCV<-NA
  thisMean<-mean(x, na.rm=TRUE); thisVar<-var(x, na.rm=TRUE)
  if(thisVar>0){
    thisCV<-sqrt(thisVar)/thisMean
  }
  return(thisCV)
}

cvByGroupTime<-apply(storeTracers,c(2,3), calcCV)

maxCVByGroup<-apply(cvByGroupTime, 1, max, na.rm=TRUE)
index<-maxCVByGroup>1e-2
plot(maxCVByGroup[index])

topCVs<-maxCVByGroup[index]; topGroups<-Ntracers[index]
topGroups[rev(order(as.double(as.character(topCVs))))]


cbind(Ntracers[rev(order(maxCVByGroup))], signif(maxCVByGroup[rev(order(maxCVByGroup))],2))



