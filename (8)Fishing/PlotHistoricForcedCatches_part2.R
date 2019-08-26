#these are summary plots of historic catch (not spatially)
this_run<-"Base"

nboxes<-30

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
plotPath<-paste(DIR$"Figures","CatchHistories\\",sep="")

catchYears<-seq(1900,2014); ny<-length(catchYears)

#read in groups csv file to get codes for fished groups
groupsDF<-read.csv(paste(this_path,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
catchCodes<-groupsDF$Code[groupsDF$IsFished==1]; ncg<-length(catchCodes)

catch_array<-read.csv(paste(DIR$'Tables',"ALLcatchHistories.csv",sep=""))
colnames(catch_array)[1]<-"Year"

totalByGroup<-colSums(catch_array[,-1],na.rm=TRUE)

nTopGroups<-6

#keep out top groups by catch and lump the rest into 'other'
topGroups<-names(sort(totalByGroup,decreasing = TRUE))[1:nTopGroups]
otherIndex<-!(colnames(catch_array) %in% c("Year",topGroups))
otherByYear<-rowSums(catch_array[,otherIndex],na.rm=TRUE)/1000 #turn to tonnes

top_catch_array<-cbind("Year"=catch_array$Year,catch_array[,colnames(catch_array) %in% topGroups]/1000,"Other"=otherByYear)
colnames(top_catch_array)<-c("Year", "Black oreo", "Hoki", "Ling", "Orange roughy", "Pelagic fish medium", "Smooth oreo", "Other")
# colnames(top_catch_array)<-c("Year", "Black oreo", "Hoki", "Ling", "Orange roughy", "Pelagic fish medium", "Smooth oreo", "Javelinfish", "Mesopelagic jelly-eaters","Other")

catchColors<-c(colorRampPalette(colors=c(myGold,myGreen,myAqua,myBlue,myPurple,myRed))(nTopGroups),myGrey)

catchDF<-melt(top_catch_array,id.var="Year")
bp<-ggplot(catchDF, aes(x = Year, y = value, fill = variable)) + 
  geom_bar(stat = "identity")

pdf(paste(plotPath,"TopCatchByYear.pdf",sep=""),height=5,width=7)
bp + scale_fill_manual(values=catchColors) + labs(y="Catch (tonnes)") + theme_igray() + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14)) + guides(fill=guide_legend(title=""))
dev.off()

options(scipen=999999)
png(paste(DIR$'Figures',"TopCatchByYear.png",sep=""), bg="transparent", height=2000, width=3500, type="cairo")
bp + scale_fill_manual(values=catchColors) + labs(y="Catch (tonnes)") + ggpubr::theme_transparent() + theme(legend.text=element_text(size=80, color=myLightAqua), legend.key.size = unit(10, 'lines'))+ 
  theme(axis.text=element_text(size=80, colour =myLightAqua), legend.spacing = unit(1.5, "cm"),axis.title=element_text(size=80, colour =myLightAqua), axis.line = element_line(colour = 'white', size = 2)) + guides(fill=guide_legend(title=""))
  
dev.off()


##future catch scenarios
#get last 5 year average yearly catch for each group
yearIndex<-catch_array$Year %in% seq(2010,2014)
last5YearAve<-apply(catch_array[yearIndex,-1],2,mean,na.rm=TRUE)

scenario_array<-data.frame(array(NA,dim=c(5,dim(catch_array)[2])))
colnames(scenario_array)<-colnames(catch_array)
colnames(scenario_array)[1]<-"Scenario"
scenario_array$Scenario<-c("Base", "Increase all", "Decrease all", "Increase hoki", "Decrease hoki")

scenario_array[1,-1]<-last5YearAve
scenario_array[2,-1]<-last5YearAve*1.2
scenario_array[3,-1]<-last5YearAve*0.8
scenario_array[4,-1]<-last5YearAve
scenario_array[4,colnames(scenario_array)=="HOK"]<-scenario_array[1,colnames(scenario_array)=="HOK"]*1.2
scenario_array[5,-1]<-last5YearAve
scenario_array[5,colnames(scenario_array)=="HOK"]<-scenario_array[1,colnames(scenario_array)=="HOK"]*0.8

otherScenByYear<-rowSums(scenario_array[,otherIndex],na.rm=TRUE)/1000 #turn to tonnes

top_scenario_array<-cbind("Scenario"=scenario_array$Scenario,scenario_array[,colnames(scenario_array) %in% topGroups]/1000,"Other"=otherScenByYear)
colnames(top_scenario_array)<-c("Scenario", "Black oreo", "Hoki", "Ling", "Orange roughy", "Pelagic fish medium", "Smooth oreo", "Other")


scenarioDF<-melt(top_scenario_array,id.vars = "Scenario")
bp<-ggplot(scenarioDF, aes(x = Scenario, y = value, fill = variable)) + 
  geom_bar(stat = "identity")
pdf(paste(plotPath,"ScenarioCatchByYear.pdf",sep=""),height=5,width=7)
bp + scale_fill_manual(values=catchColors) + labs(y="Catch (tonnes)") + theme_igray() + 
  theme(axis.text=element_text(size=14, angle=90),axis.title=element_text(size=14)) + guides(fill=guide_legend(title=""))
dev.off()


#output 5 year averages so can read in elsewhere to create catch forcing files
write.csv(last5YearAve,paste(DIR$'Tables',"CatchHist_last5YearAverage.csv",sep=""),row.names = FALSE)



