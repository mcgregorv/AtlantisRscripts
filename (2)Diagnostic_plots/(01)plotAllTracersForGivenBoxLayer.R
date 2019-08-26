#plot all tracers for a given box and layer
this_run<-"base"
this_out<-"FISH21"

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")

plotPath<-paste(this_path,"..\\Figures\\",this_run,"\\",this_out,"",sep="")

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1880
fishingStartYear<-1900
modelStartYear<-1880

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
ntsteps<-dim(thisVol)[3]

nts<-dim(thisVol)[3] #number of timesteps

xLabsTemp<-seq(0,(nts*daysTimeStep),by=365)/365
xLabsAt<-xLabsTemp*numStepsPerYear
xLabs<-xLabsTemp+year0

#get all tracer names
allTracers<-sort(names(ThisNC.nc$var))

thisBox<-"all"
thisLayer<-"all"
#get layer index. 0 is closest to sediment, it is deepest. up through water column to layer 4
skip<-c("nominal_dz")
plotsFile<-paste(plotPath,"ALL_TRACERS_box",thisBox,"_l",thisLayer,".pdf",sep="")
pdf(plotsFile)
par(mfrow=c(4,1),mar=c(3,4,2,0),oma=c(1,0,0,0))
for(t in 1:(length(allTracers))){
  thisTracer<-allTracers[t]
  temp<-ncvar_get(ThisNC.nc,thisTracer)
  if(!thisTracer %in% skip){
    if(length(dim(temp))==3){
      xx<-apply(temp*thisVol,3,sum)
    } else{
      xx<-apply(temp*thisVol[6,,],2,sum)
    }
    if(xx[1]==0){
      thisY<-xx
    } else{
      thisY<-xx/xx[1]
    }
    thisymax<-max(thisY)*1.1
    thisymin<-min(0,min(thisY)*1.1)
    plot(x=seq(1,nts),y=thisY,type="l",col=myGreen,lwd=2.5,ylim=c(thisymin,thisymax),ylab="",xlab="Day",xaxt="n")
    mtext(thisTracer,side=3,adj=0,font=2)
    abline(h=1,col="red",lty=2,lwd=1.5)
    axis(at=xLabsAt,labels=xLabs,side=1)
  }
  
}
dev.off()

############
# sort(names(ThisNC.nc$var))
# 
# xx<-grep("Arrow",names(ThisNC.nc$var))
# names(ThisNC.nc$var)[xx]
# ASQ_N<-ncvar_get(ThisNC.nc,"Arrow_squid_N")
# volume<-ncvar_get(ThisNC.nc,"volume")
# 
# layerCols<-rainbow(dim(volume)[1])
# 
# par(mfrow=c(4,2))
# for(b in 1:dim(volume)[2]){
#   thisMax<-max(ASQ_N[,b,])
#   plot(x=seq(1,dim(volume)[3]),y=rep(0,dim(volume)[3]),type="n",ylim=c(0,thisMax))
#   for(l in 1:dim(volume)[1]){
#     points(x=seq(1,dim(volume)[3]),y=ASQ_N[l,b,],type="l",lwd=2,col=layerCols[l])
#   }
#   mtext(b,side=3,adj=0,line=-1)
# }
# 
# xx<-grep("^Macro",names(ThisNC.nc$var))
# test<-ncvar_get(ThisNC.nc,"Macroalgae_N")
# boxCols<-colorRampPalette(colors=c(myYellow,myOrange,"red",myRed,myPurple,myBlue,myAqua,myGreen))(dim(test[1]))
# plot(x=seq(1,dim(test)[2]),y=rep(0,dim(test)[2]),type="n",ylim=c(0,max(test)))
# par(mfrow=c(5,2),mar=c(4,4,1,0))
# for(b in 1:(dim(test)[1])){
#   plot(test[b,],type="l",lwd=2,col=boxCols[b],lty=(b%%4+1))
#   mtext(b,side=3,adj=0)
#   # points(test[b,],type="l",lwd=2,col=boxCols[b],lty=(b%%4+1))
# }
# 
# test<-ncvar_get(ThisNC.nc,"MicroPB_N")
# volume<-ncvar_get(ThisNC.nc,"volume")
# # test<-test*volume
# boxCols<-colorRampPalette(colors=c(myYellow,myOrange,"red",myRed,myPurple,myBlue,myAqua,myGreen))(dim(test)[1])
# plot(x=seq(1,dim(test)[2]),y=rep(0,dim(test)[2]),type="n",ylim=c(0,max(test)))
# par(mfrow=c(5,2),mar=c(4,4,1,0))
# par(mfrow=c(5,2),mar=c(4,4,1,0))
# for(b in 1:(dim(test)[2])){
#   plot(test[1,b,],type="n",ylim=c(0,max(test[,b,])))
#   mtext(b,side=3,adj=0)
#   for(l in 1:(dim(test)[1])){
#     points(test[l,b,],type="l",lwd=2,lty=2,col=boxCols[l])
#   }
#   legend(legend=seq(1,(dim(test)[1])),col=boxCols[1:(dim(test)[1])],x="bottom",lwd=2)
#   # points(test[b,],type="l",lwd=2,col=boxCols[b],lty=(b%%4+1))
# }
# 
# 
# # [1] "Cephalopod_other_N"        "Cephalopod_other1_Nums"    "Cephalopod_other1_ResN"    "Cephalopod_other1_StructN" "Cephalopod_other2_Nums"
# # [6] "Cephalopod_other2_ResN"    "Cephalopod_other2_StructN"
# 
# # #copmare Cep1nums and Cep2nums
# # cep1nums<-ncvar_get(ThisNC.nc,"Cephalopod_other1_Nums")
# # cep2nums<-ncvar_get(ThisNC.nc,"Cephalopod_other2_Nums")
# # C1N<-apply(cep1nums,3,sum)
# # C2N<-apply(cep2nums,3,sum)
# # plot(C2N,col=myOrange,type="l",lwd=3)
# # par(new=TRUE)
# # plot(C1N,col=myGreen,type="l",lwd=3)
# # 
# # xx<-grep("Arrow",names(ThisNC.nc$var))
# # # > names(ThisNC.nc$var)[xx]
# # # "Arrow_squid_N"        "Arrow_squid1_Nums"    "Arrow_squid1_ResN"    "Arrow_squid1_StructN" "Arrow_squid2_Nums"    "Arrow_squid2_ResN"    "Arrow_squid2_StructN"
# # 
# # #copmare Cep1nums and Cep2nums
# # cep1nums<-ncvar_get(ThisNC.nc,"Arrow_squid1_Nums")
# # cep2nums<-ncvar_get(ThisNC.nc,"Arrow_squid2_Nums")
# # C1N<-apply(cep1nums,3,sum)
# # C2N<-apply(cep2nums,3,sum)
# # plot(C2N,col=myOrange,type="l",lwd=3)
# # par(new=TRUE)
# # plot(C1N,col=myGreen,type="l",lwd=3)
# # 
