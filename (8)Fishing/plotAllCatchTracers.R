#plot all tracers for a given box and layer
this_run<-"base"

this_out<-"FS_A1All100catch"
this_out<-"Fish"

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")
r<-""

mg_2_tonne<- 0.00000002; X_CN<-5.7  

plotPath<-paste(this_path,"..\\Figures\\",this_run,"\\",this_out,"",sep="")

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1900
fishingStartYear<-1900
modelStartYear<-1900

CatchNC.nc<-nc_open(paste(outPath,"outputCATCH.nc",sep=""))
temp<-names(CatchNC.nc$var); catchTracers<-temp[grep("_Catch",temp)]

thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
ntsteps<-dim(thisVol)[3]

testCatchNC<-nc_open(paste(outPath,"output.nc",sep=""))
testSquidCatch<-ncvar_get(testCatchNC,"Arrow_squid2_Nums")
testSquidCatchNums<-apply(testSquidCatch,3,sum,na.rm=TRUE)
test<-ncvar_get(CatchNC.nc,"Arrow_squid2_Catch")
testCatch<-apply(test,2,sum)

testCatch<-ncvar_get(testCatchNC,"Lookdown_dory6_Nums")
testCatchNums<-apply(testCatch,3,sum,na.rm=TRUE)
test<-ncvar_get(CatchNC.nc,"Lookdown_dory6_Catch")
testCatch<-apply(test,2,sum)*X_CN*mg_2_tonne

LDOcatch<-array(0,dim=c(10,length(testCatch)))
for(c in 1:10){
  test<-ncvar_get(CatchNC.nc,paste("Lookdown_dory",c,"_Catch",sep=""))
  testCatch<-apply(test,2,sum)*X_CN*mg_2_tonne
  LDOcatch[c,]<-testCatch
}

testName<-"Hoki"
testCatch<-ncvar_get(testCatchNC,paste(testName,"6_Nums",sep=""))
testCatchNums<-apply(testCatch,3,sum,na.rm=TRUE)
test<-ncvar_get(CatchNC.nc,paste(testName,"6_Catch",sep=""))
testCatch<-apply(test,2,sum)*X_CN*mg_2_tonne

thiscatch<-array(0,dim=c(10,length(testCatch)))
for(c in 1:10){
  test<-ncvar_get(CatchNC.nc,paste(testName,c,"_Catch",sep=""))
  testCatch<-apply(test,2,sum)
  thiscatch[c,]<-testCatch
}
catchByYear<-apply(thiscatch,2,sum)

####################
TotCatchNC.nc<-nc_open(paste(outPath,"outputTOTCATCH.nc",sep=""))
temp<-names(TotCatchNC.nc$var); catchTracers<-temp[grep("_Catch",temp)]

testHok<-ncvar_get(TotCatchNC.nc,"Tot_HOK_Catch")
testHokCatch<-apply(testHok,2,sum)

modelYears<-seq(1900,2014)[1:length(testHokCatch)]

catchObs<-c(320.9, 884, 1407.9,  1429.5,  5191.2,  2298.5,   7270.6, 3471.2, 5784.6,  9317.9, 5068.2)
yearObs<-seq(1976,1986)
plot(x=yearObs,y=catchObs,type="h",lend=1,lwd=5)
points(x=modelYears,y=testHokCatch/12,col=myOrange,pch=8)

testAsq<-ncvar_get(TotCatchNC.nc,"Tot_ASQ_Catch")
testASQCatch<-apply(testAsq,2,sum)


HOKtest<-catch_array[,((grep("HOK",catchGroupsDF$Code))+1)]

HokLast5<-c(525992,  682864,  238593,  265682,  277086)
testHok<-721105

ASQlast5<-c(1044507,  199204,  124623,   39427,    9284)
testAsq<-140969



