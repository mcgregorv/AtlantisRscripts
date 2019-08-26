#what to do once all the predators for a group have eaten what they should and numbers still increasing..?
#should be that the group is fairly high in the foodchain, in which case add in background mortality
#so long as it is small, mL is proba OK and it's more straight forward
#this calculates missing mortality for a group in order for numbers to remain constant
this_out<-"DietsBFF"
thisTracer<-"Basketwork_eel10_Nums"

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")

daysTimeStep<-365

ThisNC.nc<-nc_open(paste(outPath,"output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
thisDz<-ncvar_get(ThisNC.nc,"dz")
ntsteps<-dim(thisVol)[3]

temp<-ncvar_get(ThisNC.nc,thisTracer)
xx<-apply(temp,3,sum,na.rm=TRUE)

xx[2]/(xx[1])

addM<-0*xx
for(i in 2:length(xx)){
  thisFrac<-xx[i]/xx[i-1]
  addM[i]<--log(thisFrac)
}

#shoul dbe 0.14
Mest<-0.14

#set difference to recruitment
thisRec<-0*xx
for(i in 2:length(xx)){
  thisRec[i]<-xx[i]-xx[i-1]
}
N0<-xx[1]

(N0+thisRec[2])*exp(-0.14)-N0
