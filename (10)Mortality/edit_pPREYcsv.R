this_path<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
preyFile<-paste(this_path,"pPREY_OUT.csv", sep="")

preycsv<-read.csv(preyFile)
colIndex<-c(1:3,grep("CEP|ASQ",colnames(preycsv)))

plot(preycsv$"ASQ")
points(preycsv$"ASQ"[preycsv$PreyAge==1],col=myBlue,pch=20)
points(preycsv$"ASQ"[preycsv$PreyAge==2],col=myOrange,pch=20)

unique(preycsv$PredatorCode[preycsv$"ASQ">0.005])

allMolluscPreds<-c("SSO", "SPE", "SPD", "SND", "SB", "RFI", "PIN", "PFM", "PFL", "ORH", "MJE", "MAC", "LIN", "LDO", "JAV", "HOK", "HAK", "ETB", "ELP", "ELI", "EIS", "EID",
                   "DPI", "CET", "CEP", "CBO", "BOE", "BIS", "BID", "BEE", "BAL", "ASQ")
allData<-preycsv[preycsv$PredatorCode %in% allMolluscPreds,colIndex]
checkZeros<-preycsv[!(preycsv$PredatorCode %in% allMolluscPreds),colIndex]
topMplluscPreds<-c("RFI", "SPE", "PFM", "HOK", "ELP", "BEE", "SPD", "SND", "PFL", "HAK", "ELI", "ASQ", "SB", "PIN", "ORH", "ETB", "DPI", "CEP", "CET")
thisData<-preycsv[preycsv$PredatorCode %in% topMplluscPreds,colIndex]
thisData$PredatorCode<-factor(as.character(thisData$PredatorCode))

points(thisData$ASQ[thisData$PreyAge==1], pch=20, col=myGreen)
points(thisData$ASQ[thisData$PreyAge==2], pch=20, col=myRed)

notThisData<-preycsv[!(preycsv$PredatorCode %in% topMplluscPreds),colIndex]
points(notThisData$ASQ, pch=20, col=myGrey)
reduce<-notThisData[notThisData$CEP>0.001,]


#increase avail of cep and asq to top preds
# index<-preycsv$PredatorCode %in% topMplluscPreds
# newPreyCsv<-preycsv
# newPreyCsv$ASQ[index]<-5*preycsv$ASQ[index]
# newPreyCsv$CEP[index]<-5*preycsv$CEP[index]

index<-preycsv$PredatorCode %in% allMolluscPreds & preycsv$PreyAge==2
newPreyCsv<-preycsv
newPreyCsv$ASQ[index]<-unlist(lapply(preycsv$ASQ[index], FUN=function(x){min(2*x, max(x,0.2))}))
newPreyCsv$CEP[index]<-unlist(lapply(preycsv$CEP[index], FUN=function(x){min(2*x, max(x,0.2))}))

## PFS predators
newPreyCsv<-preycsv
index<-newPreyCsv[,c("PFS")]>0; thisData<-newPreyCsv[index,c(1,1,2,grep("PFS",colnames(newPreyCsv)))]
newPreyCsv$PFS[index]<-newPreyCsv$PFS[index]*1.4


## write it out
# write.csv(newPreyCsv,file=preyFile,row.names = FALSE)

sumASQByGroup<-tapply(newPreyCsv$ASQ,newPreyCsv$PredatorCode,sum, na.rm=TRUE)
sumCEPByGroup<-tapply(newPreyCsv$CEP,newPreyCsv$PredatorCode,sum, na.rm=TRUE)


sumASQByGroup<-tapply(thisData$ASQ,thisData$PredatorCode,sum, na.rm=TRUE)
sumCEPByGroup<-tapply(thisData$CEP,thisData$PredatorCode,sum, na.rm=TRUE)
sumByGroup<-sumASQByGroup + sumCEPByGroup
par(lend=1)
plot(sumByGroup,type="h",lwd=5, xaxt="n"); axis(at=seq(1,length(sumByGroup)),labels = names(sumByGroup), side=1)

points(x=seq(1,length(sumByGroup))-0.1, y=sumASQByGroup,type="h",lwd=5, col=myAqua); 
points(x=seq(1,length(sumByGroup))+0.1, y=sumCEPByGroup,type="h",lwd=5, col=myOrange); 

## CEP adults
index<-thisData$PreyAge==2
sumCEPAdByGroup<-tapply(thisData$CEP[index],thisData$PredatorCode[index],sum, na.rm=TRUE)
points(x=seq(1,length(sumByGroup))+0.1, y=sumCEPAdByGroup,type="h",lwd=3, col=myRed,lty=2); 
## ASQ adults
index<-thisData$PreyAge==2
sumASQAdByGroup<-tapply(thisData$ASQ[index],thisData$PredatorCode[index],sum, na.rm=TRUE)
points(x=seq(1,length(sumByGroup))-0.1, y=sumASQAdByGroup,type="h",lwd=3, col=myDarkGreen,lty=2); 


# 
increase<-thisData[thisData$ASQ<=0.001 | thisData$CEP<=0.001,]

pred<-"RFI"
thisData<-preycsv[preycsv$PredatorCode==pred,colIndex]


# increase pred on DPI adults

index<- !is.na(preycsv$PreyAge) & preycsv$PreyAge!=1
newPreyCsv<-preycsv
newPreyCsv$DPI[index]<-unlist(lapply(newPreyCsv$DPI[index], FUN=function(x){min(2*x, 0.2)}))

# and PFS adults
newPreyCsv$PFS[index]<-unlist(lapply(preycsv$PFS[index], FUN=function(x){min(2*x,0.2)}))

## make BID similar to BIS
newPreyCsv<-preycsv
BID<-newPreyCsv[newPreyCsv$PredatorCode=='BID',]
BIS<-newPreyCsv[newPreyCsv$PredatorCode=='BIS',]
newPreyCsv[newPreyCsv$PredatorCode=='BID',-1]<-newPreyCsv[newPreyCsv$PredatorCode=='BIS',-1]

## write it out
# write.csv(newPreyCsv,file=preyFile,row.names = FALSE)




