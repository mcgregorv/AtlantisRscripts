#read in outputs from a model, then plot together with re-estimated growth curves
this_run<-"base"
this_run<-"chaos"
source(paste(DIR$'General functions',"nonZeroMean.R", sep="")); source(paste(DIR$'General functions',"makeBlankPlot.R", sep=""))
nlayers<-6
this_path<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")


this_out<-c("FISH20","FISHH1", "FISHH2", "FISHH3")
plotPath<-paste(this_path,"\\Figures\\growth\\FISHH2",sep="")

## run C is testing clearance
this_out <-c(paste("C_XXXA",1:5, sep=""), "DChaosBASE");
plotDescrip<-"ClearanceA"
plotPath<-paste(this_path,"\\Figures\\growth\\", plotDescrip  ,sep="")

# 
# this_out<-c("BASE", "BASEdistrib", "BASE50yr2");
# plotPath<-paste(this_path,"\\Figures\\growth\\TestBase50yr",sep="")
nruns<-length(this_out)
colByRun<-c(colorRampPalette(colors=c(myGold,myBlue,"midnightblue"))(nruns-1),"black")

# this_out<-c( "TEST150yrbase", "TEST150yrfish")
# plotPath<-paste(this_path,"..\\Figures\\growth\\Base_",sep="")

# this_out<-c("PFS_BHbase_R05", "PFS_BHbase_R06", "PFS_BHbase_R07", "PFS_BHbase_R08", "PFS_BHbase_R09")
# plotPath<-paste(this_path,"..\\Figures\\growth\\SSR_R0_",sep="")

nruns<-length(this_out)

mg_2_tonne<-2e-8; X_CN<-5.7; mg_2_grams<-2e-5
get_KWRR<-function(Code){
  thisPar<-paste("KWRR_",Code,sep="")
  thisOut<-get_first_number(biolLines[grep(thisPar,biolLines)])
  return(thisOut)
}
get_KWSR<-function(Code){
  thisPar<-paste("KWSR_",Code,sep="")
  thisOut<-get_first_number(biolLines[grep(thisPar,biolLines)])
  return(thisOut)
}

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
# this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\ROMS1yearRepeat\\",sep="")

nc_list<-NULL; nts_list<-NULL; min_nts<-1e+12
for(r in 1:nruns){
  outPath<-paste(this_path,"output",this_out[r],"\\",sep="")
  
  nc_list[[r]]<-nc_open(paste(outPath,"output.nc",sep=""))
  thisVol<-ncvar_get(nc_list[[r]],"volume")
  thisDz<-ncvar_get(nc_list[[r]],"dz")
  nts_list[[r]]<-dim(thisVol)[3] #number of timesteps
  if(nts_list[[r]]<min_nts){min_nts<-nts_list[[r]]}
}

xLabsTemp<-seq(0,(min_nts*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0

groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

biolLines<-readLines(paste(this_path,"..\\CRAM_BH_hybrid_biol.prm",sep=""))

storeAllGs<-array(NA,dim=c(min_nts,ng,10)); storeAllMs<-array(NA,dim=c(min_nts,ng))
storeWeightsAtAge<-array(NA,dim=c(nruns, ng, 10, min_nts)); storeNumbersAtAge<-array(NA,dim=c(nruns, ng, 10, min_nts))
storeBiomassByGroup<-array(NA,dim=c(nruns, ng, min_nts))

for(r in 1:nruns){
  thisVol<-ncvar_get(nc_list[[r]], "volume")
  for(g in 1:ng){
    thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode], side="both")
    thisNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisCode]
    if(thisNumCohorts>1){
      #then it's age structured and need numbers and individual weights
      thisBiomass<-rep(0,dim(thisVol)[3])
      for(c in 1:thisNumCohorts){
        thisVar<-paste(thisName,c,"_Nums", sep=""); temp<-ncvar_get(nc_list[[r]], thisVar)
        thisNums<-apply(temp, 3, sum)
        storeNumbersAtAge[r,g,c,]<-thisNums[1:min_nts]
        #resN
        thisVar<-paste(thisName,c,"_ResN", sep=""); temp<-ncvar_get(nc_list[[r]], thisVar)
        thisResN<-apply(temp, 3, nonZeroMean)
        #structN
        thisVar<-paste(thisName,c,"_StructN", sep=""); temp<-ncvar_get(nc_list[[r]], thisVar)
        thisStructN<-apply(temp, 3, nonZeroMean)
        thisWeights<-thisResN + thisStructN
        storeWeightsAtAge[r, g, c, ]<-thisWeights[1:min_nts]
        #biomass
        thisBiomass<-thisBiomass+thisWeights * thisNums
        
      }
      storeBiomassByGroup[r,g,]<-thisBiomass[1:min_nts]
    } else{
      thisVar<-paste(thisName, "_N", sep=""); temp<-ncvar_get(nc_list[[r]], thisVar)
      if(length(dim(temp))==3){
        thisBiomass<-apply(temp*thisVol,3,sum)
      } else{
        thisBiomass<-apply(temp*thisVol[nlayers,,],2,sum)
      }
      storeBiomassByGroup[r,g,]<-thisBiomass[1:min_nts]
    }
  }
}
nts_list

testTimeStep<-min(150, min_nts-1)

for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[g],side="both")
  cat(as.character(thisCode),"--")
  thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    #get ageclass size
    thisVar<-paste(thisCode, "_AgeClassSize", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageClassSize<-get_first_number(temp)
    
    NumsMax<-max(storeNumbersAtAge[,g,,1], na.rm=TRUE); WeightsMax<-max(storeWeightsAtAge[,g,,1], na.rm=TRUE)*mg_2_grams*X_CN
    for(r in 1:nruns){
      thisNumbers<-storeNumbersAtAge[r,g,,testTimeStep]
      thisWeights<-storeWeightsAtAge[r,g,,testTimeStep]*mg_2_grams*X_CN
      if(max(thisNumbers, na.rm=TRUE)>NumsMax){NumsMax<-max(thisNumbers, na.rm=TRUE)}
      if(max(thisWeights, na.rm=TRUE)>WeightsMax){WeightsMax<-max(thisWeights, na.rm=TRUE)}
    }
    jpeg(paste(plotPath,"CompareWeightAndNumbersAtAge_", thisCode,"_ts", testTimeStep,".jpg", sep=""), quality=5000, width=700, height=300)
    par(mfrow=c(1,2), mar=c(4,4,1,0.1), oma=c(0,0,0,0))
    par(las=1)
    plot(thisNumbers, type="n", ylim=c(0,NumsMax), xlab="Age Class", ylab="Numbers", xaxt="n", xlim=c(1,thisNumCohorts))
    axis(at=seq(1,thisNumCohorts), labels=seq(1,thisNumCohorts), side=1)
    for(r in 1:nruns){
      points(storeNumbersAtAge[r,g,,testTimeStep], type="l", lwd=2, lty=r, col=colByRun[r])
    }
    points(storeNumbersAtAge[1,g,,1], col="red", pch=8, cex=thisCex)
    mtext(gsub("_", " ", thisName), side=3, adj=0)
    #weights
    plot(thisWeights, type="n", ylim=c(0,WeightsMax), xlab="Age Class", ylab="Weights of individuals (kg)", xaxt="n", xlim=c(1,thisNumCohorts))
    axis(at=seq(1,thisNumCohorts), labels=seq(1,thisNumCohorts), side=1)
    for(r in 1:nruns){
      points(storeWeightsAtAge[r,g,,testTimeStep]*mg_2_grams*X_CN, type="l", lwd=2, lty=r, col=colByRun[r])
    }
    points(storeWeightsAtAge[1,g,,1]*mg_2_grams*X_CN, col="red", pch=8, cex=thisCex)
    mtext(gsub("_", " ", thisName), side=3, adj=0)
    legend(legend=1:5, col=colByRun, lty=1:5, x="bottomright", lwd=1.2, bty="n")
    dev.off()
  }
}

makeBlankPlot()
legend(legend=this_out, col=colByRun, lty=seq(1,nruns), lwd=2, seg.len = 5, x="center")
