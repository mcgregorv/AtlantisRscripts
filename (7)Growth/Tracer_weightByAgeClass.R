#plot all tracers for a given box and layer
this_run<-"base"

this_out<-"PreSENS2"

mg_2_tonne<-2e-8; X_CN<-5.7
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

burnin<-0 #number of years to skip in plot

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")
r<-""

plotPath<-paste(this_path,"..\\Figures\\growth\\",this_out,"",sep="")

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
nts<-dim(thisVol)[3]

groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

biolLines<-readLines(paste(this_path,"..\\CRAM_BH_hybrid_biol.prm",sep=""))


thisCode<-"HOK"; 

for(g in 1:ng){
  thisCode<-groupsDF$Code[g]
  thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode],side="both")
  thisNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisCode]
  if(thisNumCohorts>1){
    weight_array<-array(NA,dim=c(thisNumCohorts,nts))
    
    for(c in 1:thisNumCohorts){
      thisTracer<-paste(thisName,c,"_ResN",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer);   thisRN<-apply(thisTemp,3,nonZeroMean)
      thisTracer<-paste(thisName,c,"_StructN",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer); thisSN<-apply(thisTemp,3,nonZeroMean)
      weight_array[c,]<-(thisRN+thisSN)*mg_2_tonne*X_CN
    }
    tempArray<-data.frame(weight_array); tempArray$Cohort<-seq(1,thisNumCohorts)
    plotWeights<-melt(tempArray, id.var="Cohort")
    
    #get KWRR and KWSR
    kwrr<-get_KWRR(thisCode); kwsr<-get_KWSR(thisCode)
    recWeight<-(kwsr+kwrr)*mg_2_tonne*X_CN
    
    jpeg(paste(plotPath,"GrowthCurve_",thisCode,".jpg",sep=""),quality=300)
    plot(x=plotWeights$Cohort,y=plotWeights$value,col=myBlue_trans,pch=20,cex=1.5)
    points(x=seq(1,thisNumCohorts),y=weight_array[,1],col=myOrange,type="l",lwd=2)
    points(x=seq(1,thisNumCohorts),y=weight_array[,nts],col="midnightblue",type="l",lwd=2)
    mtext(thisCode,side=3,adj=0)
    points(x=1,y=recWeight,col="red",pch=8)
    dev.off()
  }
}

# #BOE
# kwsr<-267; kwrr<-707
# recWeight<-(kwsr+kwrr)*mg_2_tonne*X_CN
# 

colByTime<-colorRampPalette((colors=c(myYellow,myGreen,myAqua,myBlue,"midnightblue")))(nts)
formatCol<-function(x){
  y<-as.character(gsub("X\\.","#",x))
  return(y)
}

thisCode<-"BOE"
thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode],side="both")
thisNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisCode]

  weight_array<-array(NA,dim=c(thisNumCohorts,nts))
  
  for(c in 1:thisNumCohorts){
    thisTracer<-paste(thisName,c,"_ResN",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer);   thisRN<-apply(thisTemp,3,nonZeroMean)
    thisTracer<-paste(thisName,c,"_StructN",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer); thisSN<-apply(thisTemp,3,nonZeroMean)
    weight_array[c,]<-(thisRN+thisSN)*mg_2_tonne*X_CN
  }
  colnames(weight_array)<-colByTime
  tempArray<-data.frame(weight_array); tempArray$Cohort<-seq(1,thisNumCohorts)
  plotWeights<-melt(tempArray, id.var="Cohort")
  plotWeights$Color<-unlist(lapply(plotWeights$variable, FUN=formatCol))
  plot(x=plotWeights$Cohort,y=plotWeights$value,col=plotWeights$Color,pch=20,cex=1.5)
  
  
weightC2<-weight_array[2,1]/(mg_2_tonne*X_CN); weightC1<-weight_array[1,1]/(mg_2_tonne*X_CN); ageClassSize<-12
  G1<-(log(weightC2) - log(weightC1))/(ageClassSize*365)
  
weightC3<-weight_array[3,1]/(mg_2_tonne*X_CN)
G2<-(log(weightC3) - log(weightC2))/(ageClassSize*365)

temp<-ncvar_get(ThisNC.nc,paste(thisName,"1_Nums",sep="")); thisNums1<-sum(temp[,,1])
temp<-ncvar_get(ThisNC.nc,paste(thisName,"2_Nums",sep="")); thisNums2<-sum(temp[,,1])
nyears<-150
numsVec1<-rep(NA,nyears); numsVec2<-numsVec1
weightVec1<-numsVec1; weightVec1[1]<-weightC1; numsVec1[1]<-thisNums1
weightVec2<-numsVec1; weightVec2[1]<-weightC2; numsVec2[1]<-thisNums2

storeG2s<-rep(NA,nyears); 

testM<-0.12
plot(x=plotWeights$Cohort,y=plotWeights$value,col=plotWeights$Color,pch=20,cex=1.5)

for(t in 2:nyears){
  newWeight<-weightVec1[t-1]*exp(G1*365)
  newNums<-(numsVec1[t-1]*(ageClassSize-1)/ageClassSize)*exp(-testM)
  recWeight<-weightVec1[1]; recNums<-numsVec1[1]
  totalNewWeight<-(recWeight*recNums) + (newWeight*newNums)
  totalWeightPerInd<-totalNewWeight/(newNums+recNums)
  weightVec1[t]<-totalWeightPerInd; numsVec1[t]<-(newNums+recNums)
  #second cohort
  thisG<-calcG(Wex=weightVec2[t-1], Nex=numsVec2[t-1], Wrec=weightVec1[t], Nrec=numsVec1[t]/ageClassSize, A=ageClassSize, M=testM)
  storeG2s[t]<-thisG
  
  newWeight2<-weightVec2[t-1]*exp(thisG*365)
  newNums2<-(numsVec2[t-1]*(ageClassSize-1)/ageClassSize)*exp(-testM)
  recWeight<-weightVec1[t]; recNums<-numsVec1[t]*(1/ageClassSize)
  
  totalNewWeight<-(recWeight*recNums) + (newWeight2*newNums2)
  totalWeightPerInd<-totalNewWeight/(newNums2+recNums)
  weightVec2[t]<-totalWeightPerInd; numsVec2[t]<-(newNums2+recNums)
  
}
weight1t<-weightVec1*(mg_2_tonne*X_CN); weight2t<-weightVec2*(mg_2_tonne*X_CN)
thisCols<-colorRampPalette(colors=c(myGrey_trans,"black"))(nyears)
points(x=rep(1,nyears),y=weight1t,pch=4,col=thisCols)

points(x=rep(2,nyears),y=weight2t,pch=4,col=thisCols)


##suppose M, then calculate G
t=2

newWeight_existing<-weight_existing*exp(G1*365)

newWeight<-weightVec1[t-1]*exp(G1*365)
newNums<-(numsVec1[t-1]*(ageClassSize-1)/ageClassSize)*exp(-testM)
recWeight<-weightVec1[1]; recNums<-numsVec1[1]
totalNewWeight<-(recWeight*recNums) + (newWeight*newNums)
totalWeightPerInd<-totalNewWeight/(newNums+recNums)
weightVec1[t]<-totalWeightPerInd; numsVec1[t]<-(newNums+recNums)



calcG<-function(Wi, Wex, Nex, Wrec, Nrec, A, M){
  Nnew<-Nex*((A-1)/A)*exp(-M)
  G<-(1/365)*log( (Wi*(Nnew + Nrec) - Wrec * Nrec) / (Wex * Nnew))
  return(G)
}

thisG1<-calcG(Wex=weightVec1[t-1], Nex=numsVec1[t-1], Wrec=weightVec1[1], Nrec=numsVec1[1], A=ageClassSize, M=testM)

thisG2<-calcG(Wex=weightVec2[t-1], Nex=numsVec2[t-1], Wrec=weightVec1[t], Nrec=numsVec1[t]/ageClassSize, A=ageClassSize, M=testM)


## do array version
storeGs<-array(NA, dim=c(thisNumCohorts,nyears))
storeNumbers<-array(NA,dim=c(thisNumCohorts, nyears)); storeWeights<-storeNumbers 
calcM<-rep(NA,thisNumCohorts)
#fill weights and numbers from initial conditions (or first timestep from a run)
for(c in 1:thisNumCohorts){
  thisTracer<-paste(thisName,c,"_ResN",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer);   thisRN<-apply(thisTemp,3,nonZeroMean)[1]
  thisTracer<-paste(thisName,c,"_StructN",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer); thisSN<-apply(thisTemp,3,nonZeroMean)[1]
  storeWeights[c,1]<-(thisRN+thisSN)
  ## numbers
  thisTracer<-paste(thisName,c,"_Nums",sep="");   thisTemp<-ncvar_get(ThisNC.nc,thisTracer);
  storeNumbers[c,1]<-sum(thisTemp[,,1],na.rm=TRUE)
  if(c>1){
    calcM[c]<-(1/ageClassSize)*(-1)*log(storeNumbers[c,1]/storeNumbers[c-1,1])
  }
}

testM<-mean(calcM,na.rm=TRUE)
# testM<-0.12
##need to calculate growth for first cohort. 
kwrr<-get_KWRR(thisCode); kwsr<-get_KWSR(thisCode)
recWeight<-(kwsr+kwrr)
## starting G for cohort 1 doesn't take into account M, but OK place to start
G1<-(log(storeWeights[2,1]) - log(storeWeights[1,1]))/(ageClassSize*365/2)
G1=0
# G1<-1e-4

for(t in 2:nyears){
  for(c in 1:(thisNumCohorts-1)){
    if(c>1){
      #second cohort and higher
      thisG<-calcG(Wi=storeWeights[c,1]  ,Wex=storeWeights[c,t-1], Nex=storeNumbers[c,t-1], Wrec=storeWeights[c-1,t], Nrec=storeNumbers[c-1,t]/ageClassSize, A=ageClassSize, M=testM)
      storeGs[c,t]<-thisG
      recWeight<-storeWeights[c-1,t]; recNums<-storeNumbers[c-1,t]
    } else{
      thisG<-G1
      recWeight<-storeWeights[1,1]; recNums<-storeNumbers[1,1]
    }
    newWeight<-storeWeights[c,t-1]*exp(thisG*365)
    newNums<-(storeNumbers[c,t-1]*(ageClassSize-1)/ageClassSize)*exp(-testM)
    totalNewWeight<-(recWeight*recNums) + (newWeight*newNums)
    totalWeightPerInd<-totalNewWeight/(newNums+recNums)
    storeWeights[c,t]<-totalWeightPerInd; storeNumbers[c,t]<-(newNums+recNums)
 
  }
}
  
thisCols<-colorRampPalette(colors=c(myGrey_trans,"black"))(nyears)

plot(x=plotWeights$Cohort,y=plotWeights$value,col=plotWeights$Color,pch=20,cex=1.5,ylim=c(0,max(plotWeights$value))
for(c in 1:thisNumCohorts){
  points(x=rep(c,nyears), y=storeWeights[c,]*(mg_2_tonne*X_CN), pch=4, col=thisCols)
  points(x=c,y=storeWeights[c,nyears]*(mg_2_tonne*X_CN),pch=8,col="red")
}

plot(x=seq(1,10),y=storeGs[,nyears],type="h",lwd=5,lend=1)
plot(x=seq(1,10),y=exp(storeGs[,nyears]),type="h",lwd=5,lend=1)

exp(storeGs[,nyears])*storeWeights[,1]

exp(storeGs[,nyears])

tesetNC.nc<-nc_open(paste(this_path,"\\outputMUMtest\\output.nc",sep=""))
testVol<-ncvar_get(tesetNC.nc,"volume"); test_nts<-dim(testVol)[3]
thisCode<-"BOE"
thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode],side="both")
thisNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisCode]
  weight_array<-array(NA,dim=c(thisNumCohorts,test_nts))
  
  for(c in 1:thisNumCohorts){
    thisTracer<-paste(thisName,c,"_ResN",sep="");   thisTemp<-ncvar_get(tesetNC.nc,thisTracer);   thisRN<-apply(thisTemp,3,nonZeroMean)
    thisTracer<-paste(thisName,c,"_StructN",sep="");   thisTemp<-ncvar_get(tesetNC.nc,thisTracer); thisSN<-apply(thisTemp,3,nonZeroMean)
    weight_array[c,]<-(thisRN+thisSN)*mg_2_tonne*X_CN
  }
  tempArray<-data.frame(weight_array); tempArray$Cohort<-seq(1,thisNumCohorts)
  plotWeights<-melt(tempArray, id.var="Cohort")
  
  #get KWRR and KWSR
  kwrr<-get_KWRR(thisCode); kwsr<-get_KWSR(thisCode)
  recWeight<-(kwsr+kwrr)*mg_2_tonne*X_CN
  
 plot(x=plotWeights$Cohort,y=plotWeights$value,col=myBlue_trans,pch=20,cex=1.5,ylim=c(0,max(plotWeights$value)))
  points(x=seq(1,thisNumCohorts),y=weight_array[,1],col=myOrange,type="l",lwd=2)
  points(x=seq(1,thisNumCohorts),y=weight_array[,test_nts],col="midnightblue",type="l",lwd=2)
  mtext(thisCode,side=3,adj=0)
  points(x=1,y=recWeight,col="red",pch=8)
  
