## give this as many model runs as you like - it will test for correlation between species. 
## best if they are all the same length, but can cope if not. #set a max number of timesteps first up to set up arrays
## first just using ROMS runs, but will add other sensitivities later

#and set a base model for reading in box and layer dimensions and tracers. Only considers the _N tracers, converted to tonnes

X_CN<-5.7; mg_2_tonne<-2e-8

basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
plotPath<-paste(DIR$'Figures',"Sensitivities\\CORR_", sep="")
path1<-paste(basePath,"base\\ouputROMS",Version,"\\",sep="")
models1<-paste(path1,paste("outputROMSBootstrapB",seq(6,47),"\\output.nc", sep=""),sep="")

BaseNC.nc<-nc_open(paste(path1, "outputROMSBootstrapB6\\output.nc",sep=""))
allTracers<-names(BaseNC.nc$var); thisVol<-ncvar_get(BaseNC.nc, "volume"); nlayers<-dim(thisVol)[1]; nboxes<-dim(thisVol)[2];
x<-grep("_N",allTracers); temp<-allTracers[x]
y<-grep("Nums",temp, invert = TRUE); Ntracers<-temp[y]; ntracers<-length(Ntracers)


max_nts<-151

ALLmodels<-c(models1); nModels<-length(ALLmodels)

storeALLsensTracers<-array(NA, dim=c(nModels, ntracers, max_nts))

# storeALLsensTracers<-storeTracersByCell
##summarise by timestep, tracer and model (cut out spatial)
for(r in 1:nModels){
  cat(r," -- ")
  thisOutFile<-ALLmodels[r]
  ThisNC.nc<-nc_open(thisOutFile)
  thisVol<-ncvar_get(ThisNC.nc,"volume"); this_nts<-dim(thisVol)[3]
  for(t in 1:ntracers){
    thisTracer<-Ntracers[t]; thisData<-ncvar_get(ThisNC.nc, thisTracer)
    if(length(dim(thisData))==3){
      yy<-apply(thisData*thisVol,3,sum)*mg_2_tonne*X_CN
      storeALLsensTracers[r,t,1:min(max_nts,this_nts)]<-yy[1:min(max_nts,this_nts)]
    } else{
      yy<-apply(thisData*thisVol[nl,,],2,sum)*mg_2_tonne*X_CN
      storeALLsensTracers[r,t,1:min(max_nts,this_nts)]<-yy[1:min(max_nts,this_nts)]
    }
  }
}




storeCorrelations<-array(NA, dim=c(ntracers, ntracers))

for(t in 1:ntracers){
  for(u in 1:ntracers){
    thisTdata<-as.vector(storeALLsensTracers[,t,]); thisUdata<-as.vector(storeALLsensTracers[,u,])
    xTdata<-thisTdata[!is.na(thisTdata)]; xUdata<-thisUdata[!is.na(thisUdata)]  
    thisCorr<-rcorr(xTdata, xUdata)
    storeCorrelations[t,u]<-thisCorr[[1]][1,2]
  }
}
posColRamp<-colorRampPalette(colors=c(myLightBlue, myBlue, "midnightblue"))(101)
negColRamp<-colorRampPalette(colors=c(myYellow,myGold,myOrange,"red"))(101)
maxNeg<-max(-1,myRounding(min(storeCorrelations, na.rm=TRUE), 0.01,direction='down')*2)
getCorrColor<-function(x){
  thisColor<-myGrey_trans
  if(length(x)>0){
    if(!is.na(x)){
      if(x==0){
        thisColor<-"white"
      } else if(x>0){
        thisColor<-posColRamp[myRounding(x,0.01)*100+1]
      } else if(x<0){
        thisColor<-negColRamp[myRounding(x/maxNeg,0.01)*100+1]
      }
    }
  }
  return(thisColor)
}

plotColour<-apply(storeCorrelations, c(1,2), getCorrColor)
plotData<-storeCorrelations

tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotData)[1]),dim(plotData)[2]),"y"=sort(rep(seq(1,dim(plotData)[2]),dim(plotData)[1]))))

# pdf(paste(plotPath,"CorrbyTracer.pdf",sep=""),height=4,width=7)
par(mar=c(4,7,1.5,1),las=1)
plot(x=seq(1,dim(plotData)[1]),y=rep(dim(plotData)[2],dim(plotData)[1]),type="n",ylab="",xlab="Timestep (years)",xaxt="n",yaxt="n",ylim=c(0.5,dim(plotData)[2]+0.5))
axis(at=seq(1,dim(plotData)[2]),labels = Ntracers,side=2,las=1)
axis(at=seq(1,dim(plotData)[1]),labels=Ntracers,side=1,las=2)
temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)

## sort the tracers by trophic level
trophicLevelDF<-read.csv(paste(basePath,"inputs\\biol_prm\\CRAM_trophicLevels_isotopes.csv",sep=""))
trophicLevelDF$TL<-trophicLevelDF$Isotope
trophicLevelDF$TL[is.na(trophicLevelDF$Isotope)]<-trophicLevelDF$TrophicLevel2[is.na(trophicLevelDF$Isotope)]

xx<-trophicLevelDF$Name[order(trophicLevelDF$TL)]; orderedNames<-unlist(lapply(xx, FUN=function(x){str_trim(x, side="both")}))
orderedTracers<-paste(orderedNames,"_N", sep="")

orderdCorr<-storeCorrelations[match(orderedTracers, Ntracers), match(orderedTracers, Ntracers)]

plotColour<-apply(orderdCorr, c(1,2), getCorrColor)
plotData<-orderdCorr

tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotData)[1]),dim(plotData)[2]),"y"=sort(rep(seq(1,dim(plotData)[2]),dim(plotData)[1]))))

pdf(paste(plotPath,"CorrbyTracer.pdf",sep=""),height=10,width=10)
par(mar=c(9,9,1.5,1),las=1)
plot(x=seq(1,dim(plotData)[1]),y=rep(dim(plotData)[2],dim(plotData)[1]),type="n",ylab="",xlab="",xaxt="n",yaxt="n",ylim=c(0.5,dim(plotData)[2]+0.5))
axis(at=seq(1,dim(plotData)[2]),labels = gsub("_"," ",orderedNames),side=2,las=1)
axis(at=seq(1,dim(plotData)[1]),labels=gsub("_", " ", orderedNames),side=1,las=2)
temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
dev.off()


#check out a couple
# MB and EIS - negatively correlated
t1<-grep("MicroPB_N", Ntracers); t2<-grep("Epiben_fish_shal_N", Ntracers)
xt1<-storeALLsensTracers[,t1,]; xt2<-storeALLsensTracers[,t2,]
plot(x=as.vector(xt1), y=as.vector(xt2))
plot(xt1[1,], type="l",col=myBlue)
par(new=TRUE)
plot(xt2[1,], type="l",col=myOrange)

