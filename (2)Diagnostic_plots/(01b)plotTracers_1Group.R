#plot all tracers for a given box and layer
this_run<-"base"

this_out<-"Short2"
this_out<-"TestFish"
# this_out<-"ZeroSomePred"
# this_out<-"ShortLonger"
this_out<-"NUMS4"
# this_out<-"Chla"
this_out<-"WeightedDiets2"
this_out<-"MultiSpawn"
# this_out<-"Zoo"

mg_2_tonne<-2e-8; X_CN<-5.7

burnin<-0 #number of years to skip in plot

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")
r<-""


plotPath<-paste(this_path,"..\\Figures\\",this_run,"\\",this_out,"",sep="")

daysTimeStep<-5
numStepsPerYear<-365/daysTimeStep
year0<-1900
fishingStartYear<-1900
modelStartYear<-1900

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
ntsteps<-dim(thisVol)[3]


nts<-dim(thisVol)[3]-burnin #number of timesteps
cat(paste(nts,"\n"))
xLabsTemp<-seq(0,(nts*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0

#get all tracer names
allTracers<-sort(names(ThisNC.nc$var))

cohortCols<-colorRampPalette(colors=c(myYellow,myGreen,myAqua,myBlue,myPurple,myRed,"red"))(10)

testVar<-"Chl_a"
temp<-ncvar_get(ThisNC.nc,testVar)
xx<-apply(temp*thisVol,3,sum,na.rm=TRUE)
plot(xx/xx[1],type="l",ylim=c(0,1.))


testVar<-"NH3"
nh3<-ncvar_get(ThisNC.nc,testVar)
xx<-apply(nh3*thisVol,3,sum,na.rm=TRUE)
plot(xx/xx[1],type="l",ylim=c(0,1.))

testVar<-"NO3"
no3<-ncvar_get(ThisNC.nc,testVar)
xx<-apply(no3*thisVol,3,sum,na.rm=TRUE)
plot(xx/xx[1],type="l",ylim=c(0,1.))

testVar<-"Oxygen"
o2<-ncvar_get(ThisNC.nc,testVar)
xx<-apply(o2*thisVol,3,sum,na.rm=TRUE)
plot(xx/xx[1],type="l",ylim=c(0,2))

testVar<-"Pelagic_fish_sml_N"
temp<-ncvar_get(ThisNC.nc,testVar)
xx<-apply(temp*thisVol,3,sum,na.rm=TRUE)
plot(xx,type="l",ylim=c(0,max(xx)))

xxTonnes<-xx*mg_2_tonne*X_CN

par(mfrow=c(2,2),mar=c(4,4,1,1))
for(c in 1:4){
  testVar<-paste("Pelagic_fish_sml",c,"_Nums",sep="")
  temp<-ncvar_get(ThisNC.nc,testVar)
  xx<-apply(temp,3,sum,na.rm=TRUE)
  plot(xx,type="l",xaxt="n",xlab="")
  mtext(c,side=3,adj=1,line=-1.3)
  axis(at=seq(1,length(xx)),labels=seq(daysTimeStep,daysTimeStep*length(xx),by=daysTimeStep),side=1)
  
}

# SSBlines<-read.csv(paste(outPath,"outputSSB.txt",sep=""),sep=" ")
# YOYlines<-read.csv(paste(outPath,"outputYOY.txt",sep=""),sep=" ")

par(mfrow=c(5,2),mar=c(1,4,1,1))
for(c in 1:4){
  testVar<-paste("Pelagic_fish_sml",c,"_StructN",sep="")
  temp<-ncvar_get(ThisNC.nc,testVar)
  xx<-apply(temp,3,sum,na.rm=TRUE)
  testVar<-paste("Pelagic_fish_sml",c,"_ResN",sep="")
  temp<-ncvar_get(ThisNC.nc,testVar)
  yy<-apply(temp,3,sum,na.rm=TRUE)
  plot(xx+yy,type="l")
  mtext(c,side=3,adj=1,line=-1.3)
  
}

par(mfrow=c(5,2),mar=c(1,4,1,1))
for(c in 1:4){
  testVar<-paste("Pelagic_fish_sml",c,"_StructN",sep="")
  temp<-ncvar_get(ThisNC.nc,testVar)
  xx<-apply(temp,3,sum,na.rm=TRUE)
  testVar<-paste("Pelagic_fish_sml",c,"_ResN",sep="")
  temp<-ncvar_get(ThisNC.nc,testVar)
  yy<-apply(temp,3,sum,na.rm=TRUE)
  plot(yy/xx,type="l")
  mtext(c,side=3,adj=1,line=-1.3)
  
}

totalSSB<-rep(0,nts)

par(mfrow=c(2,2),mar=c(4,4,1,1))
for(c in 1:4){
  testVar<-paste("Pelagic_fish_sml",c,"_StructN",sep="")
  temp<-ncvar_get(ThisNC.nc,testVar)
  xx<-apply(temp,3,sum,na.rm=TRUE)
  testVar<-paste("Pelagic_fish_sml",c,"_ResN",sep="")
  temp<-ncvar_get(ThisNC.nc,testVar)
  yy<-apply(temp,3,sum,na.rm=TRUE)
  testVar<-paste("Pelagic_fish_sml",c,"_Nums",sep="")
  temp<-ncvar_get(ThisNC.nc,testVar)
  zz<-apply(temp,3,sum,na.rm=TRUE)
  thisSSB<-((xx+yy)*zz)*mg_2_tonne*X_CN
  totalSSB<-totalSSB+thisSSB
  plot(thisSSB,type="l",xaxt="n",xlab="")
  axis(at=seq(1,length(yy)),labels=seq(daysTimeStep,daysTimeStep*length(yy),by=daysTimeStep),side=1)
  mtext(c,side=3,adj=1,line=-1.3)
  
}
x=totalSSB/(mg_2_tonne*X_CN)
a<-1744191200; b<-398724082935
calcR<-function(x){
  R<-(a*x)/(b+x)
  return(R)
}
testR<-unlist(lapply(x,calcR))
plot(x=x,y=testR,type="l",ylim=c(0,a))

c=1
testVar<-paste("Pelagic_fish_sml",c,"_Nums",sep="")
temp<-ncvar_get(ThisNC.nc,testVar)
xx<-apply(temp,3,sum,na.rm=TRUE)
xx

#need R(5.97e+14) = 4.8e+11
R0<-6e+14; B0<-4.8e+11; h<-0.7
alpha<-(4*h*R0)/(5*h-1); beta<-(B0*(1-h))/(5*h-1)

##read in recruit_hdistrib
thisBiolFile<-paste(this_path,"..\\CRAM_BH_hybrid_biol.prm",sep="")
biolLines<-readLines(thisBiolFile)

xx<-grep("PFS_recruit_hdis",biolLines)
thisVec<-get_first_number(biolLines[xx+1],n="all")

# 
# par(mfrow=c(5,2),mar=c(1,2,1,1))
# for(c in 1:10){
#   testVar<-paste("Invert_comm_Scav",c,"_Nums",sep="")
#   temp<-ncvar_get(ThisNC.nc,testVar)
#   xx<-apply(temp,3,sum,na.rm=TRUE)
#   plot(xx/xx[1],type="l")
#   mtext(c,side=3,adj=1,line=-1.3)
#   
# }
# 
# par(mfrow=c(5,2),mar=c(1,2,2,1),las=1)
# for(c in 1:4){
#   testVar<-paste("Invert_comm_Scav",c,"_StructN",sep="")
#   temp<-ncvar_get(ThisNC.nc,testVar)
#   xx<-apply(temp,3,sum,na.rm=TRUE)
#   testVar<-paste("Invert_comm_Scav",c,"_ResN",sep="")
#   temp<-ncvar_get(ThisNC.nc,testVar)
#   yy<-apply(temp,3,sum,na.rm=TRUE)
#   plot(yy+xx,type="l")
#   mtext(c,side=3,adj=1,line=-1.3)
#   axis(at=seq(1,length(yy)),labels=seq(daysTimeStep,daysTimeStep*length(yy),by=daysTimeStep),side=1)
# }
# mtext("Weight of individuals",side=3,outer=TRUE,line=-1)
# 
# par(mfrow=c(5,2),mar=c(1,2,2,1),las=1)
# for(c in 1:4){
#   testVar<-paste("Invert_comm_Scav",c,"_StructN",sep="")
#   temp<-ncvar_get(ThisNC.nc,testVar)
#   xx<-apply(temp,3,sum,na.rm=TRUE)
#   testVar<-paste("Invert_comm_Scav",c,"_ResN",sep="")
#   temp<-ncvar_get(ThisNC.nc,testVar)
#   yy<-apply(temp,3,sum,na.rm=TRUE)
#   plot(yy/xx,type="l")
#   mtext(c,side=3,adj=1,line=-1.3)
#   axis(at=seq(1,length(yy)),labels=seq(daysTimeStep,daysTimeStep*length(yy),by=daysTimeStep),side=1)
# }
# mtext("RN/SN",side=3,outer=TRUE,line=-1)
# 
# 
# temperature<-ncvar_get(ThisNC.nc,"Temp")
# salt<-ncvar_get(ThisNC.nc,"salt")
