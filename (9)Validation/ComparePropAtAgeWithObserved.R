#bring in the run, the observed data, the catch history and plot together
source(paste(DIR$'General functions',"getCIfromCV.R", sep=""))
source(paste(DIR$'General functions',"getCIfromCVs.R", sep=""))
nboxes<-30

mg_2_tonne<-2e-8; X_CN<-5.7
catchYears<-seq(1900,2014); ny<-length(catchYears)

#################################### TRACERS
#plot all tracers for a given box and layer
this_run<-"base"
base_out<-"BASEEddy"
this_out<-"FISH"

this_run <- "Archivedmodels"

burnin<-35 #number of years to skip in plot

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")
baseOutPath<-paste(this_path,"output",base_out,"\\", sep="")
r<-""
biolLines<-readLines(paste(this_path,"..\\CRAM_BH_hybrid_biol.prm", sep=""))

#read in groups csv file to get codes for fished groups
groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
catchCodes<-groupsDF$Code[groupsDF$IsFished==1]; ncg<-length(catchCodes)

catchPath<-paste(this_path,"inputs\\catch_history\\",sep="")
#read in catch ARRAY
catch_array<-read.csv(paste(DIR$'Tables',"ALLcatchHistories.csv",sep=""))

# plotPath<-paste(DIR$'Figures',"Validation\\",sep="") ## use this to overwrite the figure used in the report
plotPath<-paste(DIR$'Figures',"\\Validation\\testing\\",this_out,sep="")

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1900
fishingStartYear<-1900
modelStartYear<-1900

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
ntsteps<-dim(thisVol)[3]

BaseNC.nc<-nc_open(paste(baseOutPath,"output.nc",sep=""))
baseVol<-ncvar_get(BaseNC.nc, "volume")

harvestLines<-readLines(paste(this_path,"..\\CRAM_harvest_short.prm", sep=""))

nts<-dim(thisVol)[3]-burnin #number of timesteps
cat(paste(nts,"\n"))
xLabsTemp<-seq(0,(nts*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0

#get all tracer names
allTracers<-sort(names(ThisNC.nc$var))

B0data<-read.csv(paste(this_path,"..\\CRAM_B0.csv",sep=""))
## store SSB tracers for each group -only fill with age-structured
storeNUMSbase<-array(NA, dim=c(ng,10,nts)); storeNUMSfish<-storeNUMSbase
storeWEIGHTSbase<-storeNUMSbase; storeWEIGHTSfish<-storeWEIGHTSbase
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    cat(as.character(thisCode))
    thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode],side="both")
    #if only 2 cohorts, just the adults are SSB
    for(c in 1:thisNumCohorts){
      thisVar<-paste(thisName,c,"_Nums", sep=""); tempNums<-ncvar_get(BaseNC.nc,thisVar); fishNums<-ncvar_get(ThisNC.nc,thisVar);
      thisNUMS<-apply(tempNums, 3, sum, na.rm=TRUE); thisFishNUMS<-apply(fishNums, 3, sum, na.rm=TRUE)
      storeNUMSbase[g,c,]<-thisNUMS[burnin:(burnin+nts-1)]; storeNUMSfish[g,c,]<-thisFishNUMS[burnin:(burnin+nts-1)]
      #get weights
      thisVar<-paste(thisName,c,"_ResN", sep=""); tempResN<-ncvar_get(BaseNC.nc,thisVar); fishResN<-ncvar_get(ThisNC.nc,thisVar);
      thisVar<-paste(thisName,c,"_StructN", sep=""); tempStructN<-ncvar_get(BaseNC.nc,thisVar); fishStructN<-ncvar_get(ThisNC.nc,thisVar);
      thisWeightsBase<-apply(tempResN+tempStructN,3,nonZeroMean)*mg_2_tonne*X_CN*1e+3
      thisWeightsFish<-apply(fishResN+fishStructN,3,nonZeroMean)*mg_2_tonne*X_CN*1e+3
      storeWEIGHTSbase[g,c,]<-thisWeightsBase[burnin:(burnin+nts-1)]; storeWEIGHTSfish[g,c,]<-thisWeightsFish[burnin:(burnin+nts-1)]; 
    }
  }
}

# 
# 
# g=38 # ORH
# thisCode<-groupsDF$Code[g]
# baseNums<-storeNUMSbase[g,,116]
# fishNums<-storeNUMSfish[g,,116]
# thisVar<-paste("CatchTS_agedistrib", thisCode, sep=""); temp<-harvestLines[grep(thisVar, harvestLines)+1]; thisCatchPropByAge<-get_first_number(temp,n="all")
# # get age class sizes
# thisVar<-paste(thisCode,"_AgeClassSize", sep=""); temp<-biolLines[grep(thisVar, biolLines)]; thisAgeClassSize<-get_first_number(temp)
# xx<-thisCatchPropByAge * baseNums; propByAge<-xx/sum(xx)
# xx<-thisCatchPropByAge * fishNums; fishPropByAge<-xx/sum(xx)
# 
# plot(x=seq(1,length(propByAge))-0.1,y=propByAge, type="h", lwd=5, ylim=c(0,max(c(propByAge, fishPropByAge))), xaxt="n", xlab="Age (years)", ylab="Proportion at age")
# points(x=seq(1,length(propByAge))+0.1,y=fishPropByAge, type="h", lwd=5, col=myOrange)
# axis(at=seq(1,10), labels=seq(1,10)*thisAgeClassSize, side=1)
# par(new=TRUE); plot(thisCatchPropByAge, type="l", lwd=2,col=myGreen, xlab="", ylab="", xaxt="n",yaxt="n")
# 
# thisWeightsBase<-storeWEIGHTSbase[g,,nts]
# thisWeightsFish<-storeWEIGHTSfish[g,,nts]
# a<-0.0921; b<-2.71
# calc_length<-function(w,a,b){
#   l<-((w*1e+3)/a)^(1/b)
#   
#   return(l)
# }
# lengthsBase<-unlist(lapply(thisWeightsBase,calc_length,a=a,b=b))
# lengthsFish<-unlist(lapply(thisWeightsFish, calc_length,a=a,b=b))
# 
# plot(x=lengthsBase, y=propByAge, type="l", ylim=c(0,max(c(propByAge, fishPropByAge))))
# points(x=lengthsFish, y=fishPropByAge, type="l",lty=2,col=myOrange)
# 
# 
# xx<-unlist(str_split(c("0 0.0007814836 0 0 0.006726165 0.007721561 0.005457647
#                        0.005730768 0.01529701 0.01899051 0.02372589 0.02183696 0.03957551
#                        0.0432098 0.04079977 0.03448648 0.05753911 0.05046981 0.05509916
#                        0.04208103 0.05553336 0.03145382 0.03317786 0.04592065 0.02275182
#                        0.02986571 0.02657233 0.01854391 0.01539028 0.01794228 0.01548356
#                        0.01304583 0.02524875 0.01134949 0.02380642 0.01000043 0.009284529
#                        0.007677311 0.002632531 0.01363473 0.007956802 0.007742889
#                        0.01036683 0.005898211 0.00355376 0.008738286 0.006887238
#                        0.003488182 0.003274269 0.002920611 0.002558364 0.003488182
#                        0.005398438 0.0007814836 0.002780866 0.001209309 0.004335244 0
#                        0.0006417383 0 0.001637135 0 0.001423222 0.001562967 0 0 0
#                        0.0008556511 0.004335244 0 0.001562967 0.001851048 0 0.0008556511 0
#                        0.001069564 0.0007814836 0 0 0 0.003200102")," |\n|\t"))
# obsProp<-as.double(xx[xx!=""])
# obsAges<-seq(20,100)
# 
# modelAges<-seq(1,10)*thisAgeClassSize
# 
# plot(x=obsAges, y=obsProp, type="h",lwd=2, col=myGrey)
# points(x=modelAges, y=(propByAge/max(propByAge))*max(obsProp), type="l",lwd=2,col="black", lty=2)
# points(x=modelAges, y=(fishPropByAge/max(fishPropByAge))*max(obsProp), type="l",lwd=2,col=myOrange)
# par(new=TRUE); plot(thisCatchPropByAge, type="l", lwd=2,col=myGreen, xlab="", ylab="", xaxt="n",yaxt="n")
# 
# obsProps2016<-as.double(unlist(str_split("0 0 0.007056693 0.00347577 0.007161846 0.004409342 0.00419527 0.01585127 0.005152948 0.01261368 0.01737885 0.02373639 0.03635007 0.03963963 0.02694114 0.05080678 0.0372633 0.02846495 0.0159477 0.04707249 0.03701262 0.04157157 0.03643111 0.02442305 0.04627195 0.02556316 0.02332362 0.0212665 0.02657282 0.03195233 0.02040149 0.02184425 0.01575396 0.02127524 0.004538602 0.01575773 0.01140533 0.01418194 0.01363991 0.003487382 0.009239299 0.01655828 0.005448069 0.005848342 0.006324699 0.01257331 0.006749079 0.005896555 0.005035301 0.007953665 0.0009576786 0.003463275 0.00621578 0 0.001252799 0.002352231 0.004704462 0.009869023 0.004987088 0.00360503 0.0009576786 0 0.004914768 0 0.003119944 0.001228692 0.004324529 0.00390015 0.004938876 0 0.0006384524 0.001276905 0.002186371 0.001252799 0.001276905 0 0.0009576786 0.001596132 0 0 0.008203457"," ")))
# plot(x=obsAges, y=obsProps2016, type="h",lwd=2, col=myGrey)
# points(x=modelAges, y=(propByAge/max(propByAge))*max(obsProps2016), type="l",lwd=2,col="black", lty=2)
# points(x=modelAges, y=(fishPropByAge/max(fishPropByAge))*max(obsProps2016), type="l",lwd=2,col=myOrange)
# 
# 
# thisCode<-"HOK"; codeIndex<-groupsDF$Code==thisCode
# 
# thisWeights<-storeWEIGHTSbase[codeIndex,,]
# a<-0.00479; b<-2.89
# thisLengths<-apply(thisWeights, c(1,2), calc_length,a=a,b=b)
# Linf<-100.8; k<-	0.164; t0<-(-2.16)
# thisAgeClassSize<-2
# 
# VBgrowth_fn<-function(t, Linf, k, t0){
#   l<-Linf * (1- exp(-k * (t-t0)))
#   return(l)
# }
# compSize<-unlist(lapply(seq(1,10)*2, VBgrowth_fn, Linf=Linf, k=k, t0=t0))
# plot(1,type="n", ylim=c(0,max(thisLengths)), xlim=c(1,10), xaxt="n",xlab="",ylab="Length (cm)")
# for(c in 1:10){
#   xx<-sort(thisLengths[c,])
#   boxplot(xx,add=TRUE, at=c)
#   # points(x=c, y=mean(xx), pch=20,col=myBlue, cex=thisCex)
#   # points(x=rep(c,length(xx)), y=xx,pch=20, col=myBlue_trans)
# }
# points(x=seq(1,10), y=compSize, type="l", lwd=2, col=myOrange)
# axis(at=seq(1,10), labels=seq(1,10)*2, side=1)
# mtext("Age (years)", side=1, adj=0.5, line=3)
# 
# 
# thisCode<-"SPD"; ageClassSize<-3
# xx<-splitString("104.8	0.093	-3.17		0.0013	3.2639", split="\t"); pars<-xx[!is.na(xx)]
# a<-pars[4]; b<-pars[5]
# Linf<-pars[1]; k<-pars[2]; t0<-pars[3]
# 
# codeIndex<-groupsDF$Code==thisCode
# thisWeights<-storeWEIGHTSbase[codeIndex,,]
# thisLengths<-apply(thisWeights, c(1,2), calc_length,a=a,b=b)
# 
# compSize<-unlist(lapply(seq(1,10)*ageClassSize, VBgrowth_fn, Linf=Linf, k=k, t0=t0))
# pdf(paste(plotPath,thisCode,"SizeAtAge.pdf",sep=""), height=4, width=5)
# par(mar=c(4.5,4.5,1,1))
# plot(1,type="n", ylim=c(0,max(thisLengths)), xlim=c(1,10), xaxt="n",xlab="",ylab="Length (cm)",cex.axis=thisCex, cex.lab=thisCex)
# for(c in 1:10){
#   xx<-sort(thisLengths[c,])
#   boxplot(xx,add=TRUE, at=c, outline = FALSE,yaxt="n")
# }
# points(x=seq(1,10), y=compSize, type="l", lwd=2, col=myOrange)
# axis(at=seq(1,10), labels=seq(1,10)*ageClassSize, side=1, cex.axis=thisCex)
# mtext("Age (years)", side=1, adj=0.5, line=3, cex=thisCex)
# dev.off()
# 
# 
# 
# makeBlankPlot()
# legend(legend=c("Model","Observed"),col=c("black",myOrange),lwd=c(NA,2), pch=c(0,NA),x="center",bty="n")
