#summarise catch by fishing fleet, also which groups are caught in which fishery
# source(paste(DIR$'General functions',"read_boxes.R",sep=""))
# source(paste(DIR$'General functions',"\\pad_cols.R",sep=""))

this_run<-"TBGB_JP"
this_out<-"FishON"

this_path<-paste(DIR$'Base',"TBGB\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")
inputsPath<-this_path
# bgmf<-"CHAT30_aea.bgm"
# bgmFile<-paste(inputsPath,"bgm\\",bgmf,sep="")

plotPath<-paste(DIR$'Base',"TBGB\\Figures\\CatchHistories\\",sep="")

X_CN<-5.7
mg_2_tonne<-0.00000002
kg_2_mg<-1e-3/mg_2_tonne

#read in box file
this_bgm<-read_boxes(bgmFile)
ThisNC_nc<-nc_open(paste(outPath,"output.nc", sep=""))
thisVol<-ncvar_get(ThisNC_nc, "volume")
nbox<-dim(thisVol)[2] #this is the number of  boxes

#set catch history path
CH_path<-paste(inputsPath,"catch_history\\",sep="")

fleets<-read.csv(paste(inputsPath,"TBGB_Fisheries.csv", sep="")); nf<-dim(fleets)[1]

groupsDF<-read.csv(paste(this_path,"TBGB_groups.csv",sep=""))
fishedGroups<-groupsDF[groupsDF$IsFished==1,]
nfg<-dim(fishedGroups)[1]
#get time limits
stime<-as.Date("2020-01-01",format="%Y-%m-%d"); ftime<-as.Date("1900-01-01",format="%Y-%m-%d")
#get times
for(g in 1:nfg){
  #####################
  thisCode<-as.character(fishedGroups$Code[g])
  xlsxFN<-paste(thisCode,".xlsx",sep="")
  this_file<-paste(inputsPath,"catch_history\\",xlsxFN,sep="")
  wb <- loadWorkbook(this_file)
  sheetNames<-names(wb)
  sheetIndex<-grep("month-fishery",sheetNames)
  temp <- read.xlsx(xlsxFile = this_file, sheet = sheetIndex, skipEmptyRows = TRUE)
  
  thisCatchHistory<-temp[3:dim(temp)[1],]
  colnames(thisCatchHistory)<-temp[2,]
  year1<-min(as.double(thisCatchHistory$Yr))
  yearn<-max(as.double(thisCatchHistory$Yr))
  mnth1<-min(as.double(thisCatchHistory$Month[thisCatchHistory$Yr==year1]))
  mnthn<-max(as.double(thisCatchHistory$Month[thisCatchHistory$Yr==yearn]))
  if(mnth1<10){month1<-paste("0",mnth1,sep="")}else{month1=mnth1}
  this_stime<-as.Date(paste(year1,"-",month1,"-01",sep=""),format="%Y-%m-%d")
  if(mnthn<10){monthn<-paste("0",mnthn,sep="")}else{monthn=mnthn}
  this_ftime<-as.Date(paste(yearn,"-",monthn,"-01",sep=""),format="%Y-%m-%d")
  if(this_stime<stime){stime<-this_stime}
  if(this_ftime>ftime){ftime<-this_ftime}
}

all_dates<-seq.Date(from=stime,to=ftime,by="month")
nts<-length(all_dates)
#function for subtracting time
subtract_time<-function(time1,time2){
  j1<-julian(time1)
  j2<-julian(time2)
  xx<-j1-j2
  return(xx)
}
# dfs<-subtract_time(time1=all_dates,time2=all_dates[1])  #days from start
dfs<-subtract_time(time1=all_dates,time2=all_dates[1])*60*60*24   #seconds from start
nts<-length(dfs)
catch_years<-seq(1900,2012); ncatchYears<-length(catch_years)
getDate<-function(year,month){
  if(month<10){month1<-paste("0",month,sep="")}else{month1=month}
  thisDate<-as.Date(paste(year,"-",month1,"-01",sep=""),format="%Y-%m-%d")
  return(thisDate)
}
#set up array for storing catch histories by fleet
storeCatchByFleet<-array(NA, dim=c(nf, nfg, ncatchYears))

for(g in 1:nfg){
  thisCode<-as.character(fishedGroups$Code[g])
  xlsxFN<-paste(thisCode,".xlsx",sep="")
  this_file<-paste(inputsPath,"catch_history\\",xlsxFN,sep="")
  wb <- loadWorkbook(this_file)
  sheetNames<-names(wb)
  sheetIndex<-grep("month-fishery",sheetNames)
  rawData <- read.xlsx(xlsxFile = this_file, sheet = sheetIndex, skipEmptyRows = TRUE)
  
  justDates<-rawData[3:dim(rawData)[1],c(1,2)]
  colnames(justDates)<-c("year","month")
  
  thisDates<-mapply(getDate,year=justDates$year,month=justDates$month,SIMPLIFY = FALSE,USE.NAMES = FALSE)
  timeIndex<- all_dates %in% thisDates
  
  this_data<-rawData[3:(dim(rawData)[1]),]
  
  #get fleets for this group
  thisFleets<-colnames(rawData[grep("^X", colnames(rawData), invert = TRUE)])
  xx<-unique(as.double(rawData[1,c(3:dim(rawData)[2])])); thisFleetNumbers<-xx[!is.na(xx)]
  nfleets<-length(thisFleetNumbers); 
  
  for(f in 1:nfleets){
    thisFleetNum<-thisFleetNumbers[f]; fleetIndex<-as.double(rawData[1,])==thisFleetNum; fleetIndex[is.na(fleetIndex)]<-FALSE
    temp<-cbind(rawData[,1], rawData[,fleetIndex])
    if(dim(temp)[2]==2){
      rowSums<-as.double(temp[,2])
      thisFleetCatchByYear <- tapply(rowSums, temp[,1], sum, na.rm=TRUE)
    }else{
      thisFleetData<-temp[c(3:dim(temp)[1]),]; colnames(thisFleetData)[1]<-"year"
      rowSums<-apply(thisFleetData[,-1],1,sum,na.rm=TRUE)
      thisFleetCatchByYear <- tapply(rowSums, thisFleetData$year, sum, na.rm=TRUE)
    }
    thisCatchByYear <- thisFleetCatchByYear[!is.na(as.double(names(thisFleetCatchByYear)))]
    yearIndex<-match(names(thisCatchByYear), catch_years)
    storeCatchByFleet[thisFleetNum, g,yearIndex]<-as.double(thisCatchByYear)/1000 # these were in kg - turn to tonnes
  }
}

totalCatchByFleet<-apply(storeCatchByFleet,1,sum,na.rm=TRUE)

totalCatchByFleetYear<-apply(storeCatchByFleet,c(1,3), sum, na.rm=TRUE)

fleetColors<-colorRampPalette(colors=c("red",myOrange,myGold,myGreen,myAqua,myBlue))(nf)

prePlot<-data.frame(totalCatchByFleetYear)
colnames(prePlot)<-catch_years; prePlot$Fleet<-fleets$Code

toPlot<-melt(prePlot, id.var="Fleet")
bp<-ggplot(data = toPlot, aes(fill = Fleet, x = variable, y = value)) + 
  geom_bar(stat = 'identity')

breakYears<-seq(1900,2014,by=10)

pdf(paste(plotPath,"CatchesByFleetAndYear.pdf", sep=""), height=6,width=9)
par(mar=c(4,4,1,1))
bp + labs(y="Catch (tonnes)", x="Year")+ scale_fill_manual(values=fleetColors) + theme_igray() + scale_x_discrete(breaks = breakYears) + 
  theme(axis.text=element_text(size=12, angle=0),axis.title=element_text(size=12)) + guides(fill=guide_legend(title="")) 
dev.off()

##number of fish groups caught by each fleet
numGroupsByFleet<-rep(0,nf)
for(f in 1:nf){
  temp<-storeCatchByFleet[f,,]; xx<-apply(temp,1,sum,na.rm=TRUE); index<-xx>0
  countGroups<-sum(rep(1,nfg)[index])
  numGroupsByFleet[f]<-countGroups
}

## index groups in each fishery
test<-apply(storeCatchByFleet,c(1,2), sum, na.rm=TRUE)
for(f in 1:nf){
  fishedIndex<-test[f,]>0
  thisCodes<-as.character(fishedGroups$Code[fishedIndex])
  groupsIndex<-groupsDF$Code %in% thisCodes
  cat("\n",f,"\n")
  cat(paste(as.double(groupsIndex),collapse = " "))
}
rownames(test)<-fleets$Code
#table catch by species group and fishing fleet
outFile<-paste(DIR$'Tables',"TBGBCatchBySpeciesGroupAndFleet.csv",sep="")
write.csv(t(test),file=outFile,row.names = fishedGroups$Code)

## plot by species code
totalCatchByStockYear<-apply(storeCatchByFleet,c(2,3), sum, na.rm=TRUE)

stockColors<-colorRampPalette(colors=c("red",myOrange,myGold,myGreen,myAqua,myBlue))(nfg)

prePlot<-data.frame(totalCatchByStockYear)
colnames(prePlot)<-catch_years; prePlot$Stock<-as.character(fishedGroups$Code)

toPlot<-melt(prePlot, id.var="Stock")
bp<-ggplot(data = toPlot, aes(fill = Stock, x = variable, y = value)) + 
  geom_bar(stat = 'identity')

breakYears<-seq(1900,2014,by=10)

pdf(paste(plotPath,"CatchesByStockAndYear.pdf", sep=""), height=6,width=9)
par(mar=c(4,4,1,1))
bp + labs(y="Catch (tonnes)", x="Year")+ scale_fill_manual(values=stockColors) + theme_igray() + scale_x_discrete(breaks = breakYears) + 
  theme(axis.text=element_text(size=12, angle=0),axis.title=element_text(size=12)) + guides(fill=guide_legend(title="")) 
dev.off()

############## by top species
totalCatchBySpecies <- apply(storeCatchByFleet, 2, sum, na.rm=TRUE)
fishedCodes <- as.character(fishedGroups$Code)
fishedCodes[order(totalCatchBySpecies, decreasing = TRUE)]
totalCatchBySpecies[order(totalCatchBySpecies, decreasing = TRUE)]

topCodes <- fishedCodes[order(totalCatchBySpecies, decreasing = TRUE)][1:7]
topCodesPrePlot <- prePlot
index <- topCodesPrePlot$Stock %in% topCodes
topCodesPrePlot[!index,c("Stock")]<-"Other"
topCodesPrePlot$Stock <- factor(topCodesPrePlot$Stock, levels=c(topCodes, "Other"))
topCodes2plot <- melt(topCodesPrePlot, id.var="Stock")

topStockColours <-c( colorRampPalette(colors=c("red",myOrange,myGold,myGreen,myAqua,myBlue))(length(topCodes)), myGrey)

bp<-ggplot(data = topCodes2plot, aes(fill = Stock, x = variable, y = value)) + 
  geom_bar(stat = 'identity')

breakYears<-seq(1900,2014,by=10)

pdf(paste(plotPath,"CatchesByStockAndYear_top.pdf", sep=""), height=6,width=9)
par(mar=c(4,4,1,1))
bp + labs(y="Catch (tonnes)", x="Year")+ scale_fill_manual(values=topStockColours) + theme_igray() + scale_x_discrete(breaks = breakYears) + 
  theme(axis.text=element_text(size=12, angle=0),axis.title=element_text(size=12)) + guides(fill=guide_legend(title="")) 
dev.off()


pdf(paste(plotPath,"CatchesByStock_ordered.pdf", sep=""), height=6,width=9)
par(mar=c(4,4,1,1))
par(lend=1)
plot(totalCatchBySpecies[order(totalCatchBySpecies, decreasing = TRUE)]/1000, type="h", lwd=5, xaxt="n", ylab="Total biomass caught (1000 tonnes)")
abline(v=7.5, col="red", lty=2, lwd=2)
par(las=2)
axis(at=1:nfg, labels=fishedCodes[order(totalCatchBySpecies, decreasing = TRUE)], side=1)
dev.off()


storeCatchByStock <- apply(storeCatchByFleet,c(2,3), sum, na.rm=TRUE)
catchYears <- 1900:2012

# grab B0 from tracers - use initial condition
B0byfishedCode <- rep(NA, nfg)
for(g in 1:nfg){
  thisCode <- fishedGroups$Code[g]; thisName <- str_trim(fishedGroups$Name[g], side="both")
  thisTracer <- paste(thisName, "_N", sep="")
  thisData <- ncvar_get(ThisNC_nc, thisTracer)
  thisB <- apply(thisData * thisVol, 3, sum, na.rm=TRUE)*mg_2_tonne*X_CN; thisB0 <- thisB[1]
  B0byfishedCode[g]<-thisB0
}

# for those with casal files, grab SSB from these
SSBbyfishedCode <- data.frame(matrix(NA,nrow= nfg, ncol=length(catchYears)))
CASALpath <- paste(this_path, "..\\Data\\CASAL\\", sep="")
casalFiles <- list.files(CASALpath)
for(g in 1:nfg){
  thisCode <- fishedGroups$Code[g]; thisName <- str_trim(fishedGroups$Name[g], side="both")
  x <- grep(thisCode, casalFiles)
  if(length(x)>0){
    thisFolder <- casalFiles[x]
    thisFile <- paste(CASALpath, thisFolder,"\\MPD.txt", sep="")
    if(file.exists(thisFile)){
      thisLines <- readLines(thisFile)
      x <- grep("SSB", thisLines); 
      thisSSB <- get_first_number(thisLines[x[2]],n="all")
      thisSSByears <- get_first_number(thisLines[x[2]+1], n="all")
      SSBbyfishedCode[g,match(thisSSByears, catchYears)]<-as.double(thisSSB)
    }
    
  }
  #any gaps in SSBbyfishedCode, fill in with B0
  if(sum(SSBbyfishedCode[g,], na.rm=TRUE)==0){
    thisB0 <- B0byfishedCode[g]
    SSBbyfishedCode[g,]<- thisB0
  }
}

# get F, using SSBbyfishedCode
TBGB_F <- storeCatchByStock/SSBbyfishedCode 

# output catch by species and year, and also dump in 3D (species, fleet, year)

save(list=c("storeCatchByFleet","fleets","fishedCodes", "fishedGroups", "catchYears", "storeCatchByStock", "TBGB_F", "SSBbyfishedCode"), file=paste(this_path,"Catch_history\\storeCatchByFleet", sep=""))
# load(paste(this_path,"Catch_history\\storeCatchByFleet", sep=""))

# plot them to check look right
pdf(paste(plotPath, "TBGB_CatchHistoriesAndF.pdf", sep=""), height=10, width=10)
par(mfrow=c(5,2), mar=c(4,4.5,1.5,8), lend=1,  par(las=0))
for(g in 1:nfg){
  thisName <- gsub("_", " ",fishedGroups$Name[g])
  thisCatchHistory <- storeCatchByStock[g,]; thisF <- TBGB_F[g,]; thisSSB <- SSBbyfishedCode[g,]
  thisFaxis <- pretty(seq(0,max(thisF), length.out=10)); thisCatchAxis <- pretty(0:max(thisCatchHistory))
  
  plot(x=catchYears, y=thisCatchHistory, type="h", col=myGrey_trans, lwd=3, yaxt="n", xaxt="n", ylab="", xlab="")
  axis(at=thisCatchAxis, labels=thisCatchAxis, side=4, line=1, col.axis=myGrey, col.ticks=myGrey, col=myGrey)
  par(new=TRUE)
  plot(x=catchYears, y=thisSSB,type="l", lty=2, lwd=1.5, ylab="Biomass (tonnes)", xlab="", ylim=c(0, max(thisSSB, na.rm=TRUE)), cex.axis=thisCex, cex.lab=thisCex)
  par(las=0)
  mtext("Catch (tonnes)", side=4, line=3, col=myGrey, font=2)
  
  
  par(new=TRUE)
  plot(x=catchYears, y=thisF, type="l", col=myOrange, lwd=2, yaxt="n", xaxt="n", ylab="", xlab="")
  axis(at=thisFaxis, labels=thisFaxis, side=4, line=5, col.axis=myOrange, col.ticks=myOrange, col=myOrange)
  mtext("F", side=4, line=7, col=myOrange, font=2)
  
  mtext(thisName, side=3, adj=0, font=2)
}
dev.off()





