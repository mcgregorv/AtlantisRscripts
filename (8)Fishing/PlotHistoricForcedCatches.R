#these are summary plots of historic catch (not spatially)
this_run<-"Base"

nboxes<-30

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
plotPath<-paste(DIR$"Figures","CatchHistories\\",sep="")

catchYears<-seq(1900,2014); ny<-length(catchYears)

#read in groups csv file to get codes for fished groups
groupsDF<-read.csv(paste(this_path,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
catchCodes<-groupsDF$Code[groupsDF$IsFished==1]; ncg<-length(catchCodes)

catchPath<-paste(this_path,"inputs\\catch_history\\",sep="")
#read in one catch xlsx file so can get timesteps to set up array
tempCatch<-thisBiolData<-read.xlsx(paste(catchPath,"ASQ_catch.xlsx",sep=""))
thisCatch<-tempCatch[3:dim(tempCatch)[1],c(1,3:dim(tempCatch)[2])]
colnames(thisCatch)<-c("Year",seq(1,dim(thisCatch)[2]-1))

catch_array<-array(NA,dim=c(ny,ncg))

for(g in 1:ncg){
  thisCode<-catchCodes[g]
  tempCatch<-thisBiolData<-read.xlsx(paste(catchPath,thisCode,"_catch.xlsx",sep=""))
  thisCatch<-tempCatch[3:dim(tempCatch)[1],c(1,3:dim(tempCatch)[2])]
  if(dim(thisCatch)[2]>2){
    temp<-apply(thisCatch[,-1],1,sum)
    sumCatch<-tapply(temp,thisCatch[,1],sum)
  } else {
    sumCatch<-tapply(thisCatch[,2],thisCatch[,1],sum)
  }
  catch_array[match(names(sumCatch),catchYears),g]<-sumCatch
}

##add col and row names then write out catch array
colnames(catch_array)<-catchCodes
rownames(catch_array)<-catchYears
write.csv(catch_array,paste(DIR$'Tables',"ALLcatchHistories.csv",sep=""))


test<-colSums(catch_array,na.rm=TRUE)

pdf(paste(plotPath,"TotalCatchALL.pdf",sep=""),height=4)
par(mar=c(4,4.5,1,1))
plot(sort(test,decreasing=TRUE)/1000,type="h",lwd=5,lend=1,xaxt="n",xlab="",ylab="Total catch (tonnes)",col=myBlue)
par(las=2)
axis(at=seq(1,length(test)),labels=catchCodes[order(test,decreasing = TRUE)],side=1)
dev.off()

##total catch by year
test<-rowSums(catch_array,na.rm=TRUE)

pdf(paste(plotPath,"TotalCatchALLbyYear.pdf",sep=""),height=4)
par(mar=c(4,4.5,1,1))
plot(x=catchYears,y=test/1000,type="h",lwd=5,lend=1,xlab="",ylab="Total catch (tonnes)",col=myBlue)
dev.off()

##color by top 5 (HOK ORH SSO LIN BOE) and 'other' 
topCatchGroups<-c("HOK","ORH","SSO","LIN","BOE"); ntcg<-length(topCatchGroups)
topCatchYears<-seq(1965,2014); ncy<-length(topCatchYears)

topCatchs<-test[catchYears %in% topCatchYears]
catchCols<-colorRampPalette(colors=c(myYellow,myOrange,"red",myPurple,myBlue))(ncg)


pdf(paste(plotPath,"TotalCatchALLbyYearTop5Groups.pdf",sep=""),height=7)
par(mar=c(3,4.5,0,1),mfrow=c(ncg+1,1),oma=c(0,0,1,1))
for(c in 1:ntcg){
  thisGroup<-topCatchGroups[c]
  thisData<-catch_array[catchYears %in% topCatchYears,catchCodes==thisGroup]
  thisData[is.na(thisData)]<-0
  plot(x=topCatchYears,y=thisData/1000,type="h",lwd=5,lend=1,xlab="",ylab="",col=catchCols[c])
  mtext(thisGroup,side=3,adj=0,line=-1.2,cex=0.8)
  
}
thisData<-apply(catch_array[catchYears %in% topCatchYears,!(catchCodes %in% topCatchGroups)],1,sum,na.rm=TRUE)
thisData[is.na(thisData)]<-0
plot(x=topCatchYears,y=thisData/1000,type="h",lwd=5,lend=1,xlab="",ylab="",col=myGrey)
mtext("Other",side=3,adj=0,line=-1.2,cex=0.8)
mtext("Total catch (tonnes)",side=2,adj=0.5,outer=TRUE,line=-1.5)
dev.off()

###########################################
##do them all seperately too
###################################
topCatchYears<-seq(1900,2014)
for(g in 1:ncg){
  thisCode<-catchCodes[g]
  pdf(paste(plotPath,"TotalCatchbyYear_",thisCode,".pdf",sep=""),height=4)
  par(mar=c(3,4.5,2,1))
  thisData<-catch_array[catchYears %in% topCatchYears,catchCodes==thisCode]
  thisData[is.na(thisData)]<-0
  plot(x=topCatchYears,y=thisData/1000,type="h",lwd=5,lend=1,xlab="",ylab="",col=myBlue,cex.axis=thisCex)
  mtext(thisCode,side=3,adj=1,line=0.2,cex=thisCex)
  mtext("Catch (tonnes)",side=2,cex=thisCex,line=3)
  dev.off()
}
