#read in the tracers (numbers), calculate spawning stock abundance, calculate rectruited abundance, then calculate steepness
this_out<-"outputBaseLong5"; runFolder<-"TBGB_JP"; 
# this_out<-"output"; runFolder<-"TBGBReportBase";

# if(runFolder=="TBGBReportBase"){
#   minCohorts<-3
# } else{
  minCohorts <-2
# }
# this_out<-"MyRun_Fish1899_better1_codeupdate_rewriteDmatrix1"; runFolder<-"TBGBFish";
thisDesc <- paste(runFolder, this_out,sep="")

basePath<-  paste(DIR$'Base', "TBGB\\",runFolder,"\\",sep="")
outPath<-paste(basePath,"\\",this_out,"\\",sep="")
plotPath<-paste(basePath,"..\\Figures\\Testing\\SSR\\SSR",thisDesc, sep="")

groupsDF<-read.csv(paste(basePath,"\\TBGB_Groups.csv",sep="")); ng<-dim(groupsDF)[1]


ThisNC.nc<-nc_open(paste(outPath,"\\output.nc",sep=""))
thisVol<-ncvar_get(ThisNC.nc,"volume")
nts<-dim(thisVol)[3]; nboxes<-dim(thisVol)[2]; nlayers<-dim(thisVol)[1]

allTracers<-sort(unique(names(ThisNC.nc$var)))
numTracers<-allTracers[grep("_Nums",allTracers)]

YOY_df<-read.csv(paste(outPath,"outputYOY.txt",sep=""),sep=" ")

biolLines <- readLines(paste(outPath, "..\\TBGB_biol.prm",sep=""))

B0df<-read.csv(paste(basePath,"..\\TBGB_B0.csv", sep=""))

thisCode<-"ELI";

calcAlpha<-function(h,R0){
  alpha<-(4*h*R0)/ (5*h-1)
  return(alpha)
}
calcBeta<-function(h, B0){
  beta <- (B0*(1-h))/(5*h-1)
  return(beta)
}
calcR<-function(B,a,b){
  R<-(a*B)/(b+B)
  return(R)
}
thisCode<-"DEM";

# pdf(paste(plotPath,"SSRallGroups.pdf", sep=""))
# par(mfrow=c(3,2), mar=c(4,4,2,4), oma=c(2,2,2,2))
for(g in 1:ng){
  if(groupsDF$NumCohorts[g]>=minCohorts){
    thisCode <- as.character(groupsDF$Code[g])

    thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode]); thisNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisCode]
    temp<-ncvar_get(ThisNC.nc,paste(thisName,"_N",sep=""))
    
    mg_2_tonne<- 0.00000002; X_CN<-5.7
    
    yoyTimeIndex<-YOY_df$Time %in% seq(0,(nts*365),by=365)
    thisYOY<-YOY_df[,grep(thisCode,colnames(YOY_df))][yoyTimeIndex]
    # cat(thisCode,thisYOY,"\n")
    # #get KWRR and KWSR from biol.prm
    thisVar<-paste("KWRR_", thisCode, sep=""); temp<-biolLines[grep(thisVar,biolLines)]; thisKWRR<-get_first_number(temp)
    thisVar<-paste("KWSR_", thisCode, sep=""); temp<-biolLines[grep(thisVar,biolLines)]; thisKWSR<-get_first_number(temp)
    juvWeight<-thisKWSR+thisKWRR
    thisRecruitNums<-thisYOY/(juvWeight*X_CN * mg_2_tonne)
    # thisRecruitNums<-thisYOY/(juvWeight)
    #ageclass size
    thisVar<-paste(thisCode,"_AgeClassSize", sep=""); temp<-biolLines[grep(thisVar,biolLines)]; ageClassSize<-get_first_number(temp)
    
    
    temp<-ncvar_get(ThisNC.nc, paste(thisName,"1_Nums",sep=""))
    thisR0<-as.double(as.character(B0df$R0[B0df$Code==thisCode]))
    this_h<-as.double(as.character(B0df$h[B0df$Code==thisCode]))
    
    temp<-ncvar_get(ThisNC.nc, paste(thisName,"_N",sep="")); 
    tempBiomass <- (apply(temp*thisVol,3,sum)*mg_2_tonne*X_CN)
    thisBiomass<-tempBiomass[1:length(thisYOY)]
    thisB0<-as.double(as.character(B0df$B0[B0df$Code==thisCode]))
    
    thisYmax<-max(thisR0*1.2, max(thisRecruitNums[-1]))
    thisYmax<-max(thisR0*1.5)
    thisXmax<-max(c(thisB0*1.5, max(thisBiomass)))
    # thisYmax<-2*thisR0
    thisAlpha<-calcAlpha(h=this_h, R0=thisR0); thisBeta<-calcBeta(h=this_h, B0=thisB0)
    testSSB<- seq(0,thisXmax,length.out=100); testR<-unlist(lapply(testSSB, calcR, a=thisAlpha, b=thisBeta))
    
    # jpeg(paste(plotPath,"BHCompareRealised_h_",this_out,"_",thisCode,".jpg",sep=""),quality=300,height=300,width=340)
    # par(mar=c(4,4,2,2), las=0)
   par(las=0)
   plot(x=testSSB, y=testR, type="l", ylim=c(0,thisYmax),xlab="", ylab="")
     points(x=thisBiomass, y=thisRecruitNums,pch=20)
    abline(v=thisB0,col="red",lty=2)
    abline(h=thisR0,col="red",lty=2)
    axis(at=thisB0,labels=expression(B[0]),side=3,col="red",col.axis="red")
    mtext("Recruitment (numbers)",side=2,line=3,adj=0.5); mtext("SSB (tonnes)", side=1, adj=0.5, line=4)
    par(las=1)
    axis(at=thisR0,labels=expression(R[0]),side=4,col="red",col.axis="red")
    mtext(thisCode,side=3, adj=1)
     # dev.off()
    # 
    # jpeg(paste(plotPath,"BHCompareRealised_h_",this_out,"_",thisCode,"ColorByTime.jpg",sep=""),quality=300,height=300,width=340)
    # par(mar=c(4,4,2,2), las=0)
    # colByTime<-colorRampPalette(colors=c(myYellow,myGreen, myDarkGreen, "black"))(length(thisBiomass))
    # plot(x=thisBiomass, y=thisRecruitNums,pch=20,xlab="SSB (tonnes)", ylab="",col=colByTime, ylim=c(0,thisYmax))
    # abline(v=thisB0,col="red",lty=2)
    # abline(h=thisR0,col="red",lty=2)
    # axis(at=thisB0,labels=expression(B[0]),side=3,col="red",col.axis="red")
    # mtext("Recruitment (numbers)",side=2,line=3,adj=0.5)
    # par(las=1)
    # axis(at=thisR0,labels=expression(R[0]),side=4,col="red",col.axis="red")
    # dev.off()
  }
}
dev.off()
# 
# getLarvae <- function(x){
#   y<-unlist(str_split(gsub(".\\(|)", "", x), ":"))[11]
#   return(y)
# }
# getIndSpawn <- function(x){
#   y<-unlist(str_split(gsub(".\\[|]", "", x), ":|,"))[8]
#   return(y)
# }
# 
# x1 <- grep("recruit_hdistrib", loglines)
# length(x1)
# loglines[x1]
# 
# 
# # loglines2 <- readLines(paste(outPath, "log.txt", sep=""))
# x2 <- grep("recruit_hdistrib", loglines2)
# length(x2)
# loglines2[x2]
# 
# larvae1 <- unlist(lapply(loglines[x1], getLarvae))
# larvae2 <- unlist(lapply(loglines2[x2], getLarvae))
# 
# # loglines3 <- readLines(paste(outPath, "log.txt", sep=""))
# x3 <- grep("IndSpawn", loglines3)
# length(x3)
# 
# loglines4 <- readLines(paste(outPath, "log.txt", sep=""))
# x4 <- grep("Larvae", loglines4)
# length(x4)
# 
# #check out ratio or RN:SN 
par(mfrow=c(2,2))
for(g in 1:ng){
  if(groupsDF$NumCohorts[g]>2){
    thisCode <- groupsDF$Code[g]
    thisName <- str_trim(groupsDF$Name[groupsDF$Code==thisCode], side="both")
    thisRN <- ncvar_get(ThisNC.nc, paste(thisName,5,"_ResN", sep=""))
    thisSN <- ncvar_get(ThisNC.nc, paste(thisName,5,"_StructN", sep=""))
    
    rn <- apply(thisRN[-6,,], 3, nonZeroMean)
    sn <- apply(thisSN[-6,,], 3, nonZeroMean)
    
    thisRatio <- (rn+sn)/sn
    
    plot(thisRatio); abline(h=3.65, col="red")
    mtext(thisCode, side=3, adj=0)
    
     
  }
}

thisCode <- "PFL"
thisName <- str_trim(groupsDF$Name[groupsDF$Code==thisCode], side="both")
thisRN <- ncvar_get(ThisNC.nc, paste(thisName,5,"_ResN", sep=""))
thisSN <- ncvar_get(ThisNC.nc, paste(thisName,5,"_StructN", sep=""))
thisNum1s <- apply(ncvar_get(ThisNC.nc, paste(thisName, 1, "_Nums", sep="")), 3, sum, na.rm=TRUE)

rn <- apply(thisRN[-6,,], 3, nonZeroMean)
sn <- apply(thisSN[-6,,], 3, nonZeroMean)

thisRatio <- (rn+sn)/sn

par(mfrow=c(3,1), mar=c(2,4,1,1))
plot(rn, type="l")
plot(sn, type="l")
plot(thisRatio, type="l"); abline(h=3.65, col="red")

plot(thisNum1s, type="l")

w_for_s <- 3.65*sn
FSP<-0.8; KSPA<-1
step1 <- FSP * w_for_s - KSPA
SPWN <- FSP * (sn + 2.65*sn -1) + rn + sn - KSPA












