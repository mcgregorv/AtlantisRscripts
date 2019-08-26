#summarise catch by fishing fleet, also which groups are caught in which fishery
library(xlsx)
source(paste(DIR$'General functions',"read_boxes.R",sep=""))
source(paste(DIR$'General functions',"\\pad_cols.R",sep=""))

this_run<-"base"
this_out<-"Fish"

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")
inputsPath<-paste(DIR$'Base',"ATLANTISmodels\\inputs\\",sep="")
bgmf<-"CHAT30_aea.bgm"
bgmFile<-paste(inputsPath,"bgm\\",bgmf,sep="")

plotPath<-paste(DIR$'Figures',"CatchHistories\\",sep="")

X_CN<-5.7
mg_2_tonne<-0.00000002
kg_2_mg<-1e-3/mg_2_tonne

#read in box file
this_bgm<-read_boxes(bgmFile)
nbox<-24 #this is the number of dynamic boxes

#set catch history path
CH_path<-paste(inputsPath,"catch_history\\",sep="")

fleets<-read.csv(paste(inputsPath,"CRAM_Fisheries_full.csv", sep="")); nf<-dim(fleets)[1]

groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep=""))
fishedGroups<-groupsDF[groupsDF$IsFished==1,]
nfg<-dim(fishedGroups)[1]
#get time limits
stime<-as.Date("2020-01-01",format="%Y-%m-%d"); ftime<-as.Date("1900-01-01",format="%Y-%m-%d")
#get times
for(g in 1:nfg){
  #####################
  thisCode<-as.character(fishedGroups$Code[g])
  xlsxFN<-paste(thisCode,"_catch.xlsx",sep="")
  this_file<-paste(inputsPath,"catch_history\\",xlsxFN,sep="")
  wb <- loadWorkbook(this_file)
  sheetNames<-names(wb)
  sheetIndex<-grep(thisCode,sheetNames)
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
catch_years<-seq(1900,2014); ncatchYears<-length(catch_years)
getDate<-function(year,month){
  if(month<10){month1<-paste("0",month,sep="")}else{month1=month}
  thisDate<-as.Date(paste(year,"-",month1,"-01",sep=""),format="%Y-%m-%d")
  return(thisDate)
}
#set up array for storing catch histories by fleet
storeCatchByFleet<-array(NA, dim=c(nf, nfg, ncatchYears))
storeAllCatch <- array(NA, dim=c(nf, ng, ncatchYears, nboxes))
for(g in 1:nfg){
  thisCode<-as.character(fishedGroups$Code[g])
  xlsxFN<-paste(thisCode,"_catch.xlsx",sep="")
  this_file<-paste(inputsPath,"catch_history\\",xlsxFN,sep="")
  wb <- loadWorkbook(this_file)
  sheetNames<-names(wb)
  sheetIndex<-grep(thisCode,sheetNames)
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
      tempDF<-data.frame(cbind("year"=as.character(temp[3:(dim(temp)[1]),1]), "catch"=temp[3:(dim(temp)[1]),2])); tempDF$year<-as.factor(tempDF$year)
    }else{
      thisFleetData<-temp[c(3:dim(temp)[1]),]; colnames(thisFleetData)[1]<-"year"
      rowSums<-apply(thisFleetData[,-1],1,sum,na.rm=TRUE)
      tempDF<-data.frame(cbind("year"=as.character(thisFleetData$year), "catch"=rowSums)); tempDF$year<-as.factor(tempDF$year)
    }
    thisCatchByYear<-tapply(as.double(tempDF$catch),tempDF$year, sum, na.rm=TRUE)
    yearIndex<-match(names(thisCatchByYear), catch_years)
    storeCatchByFleet[thisFleetNum, g,yearIndex]<-as.double(thisCatchByYear)
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

#table catch by species group and fishing fleet
outFile<-paste(DIR$'Tables',"CatchBySpeciesGroupAndFleet.csv",sep="")
write.csv(t(test),file=outFile,row.names = fishedGroups$Code)



