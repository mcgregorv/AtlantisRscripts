## read in .ts files and add lines for burnin


this_path<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
catchPath<-paste(this_path,"inputs\\catch_history\\catchts\\",sep="")

groupsDF<-read.csv(paste(this_path,"CRAM_groups.csv",sep=""))
catchGroupsDF<-groupsDF[groupsDF$IsFished==1,]; ncg<-dim(catchGroupsDF)[1]

version<-"Burnin1865"

X_CN<-5.7
mg_2_tonne<-0.00000002
kg_2_mg<-1e-3/mg_2_tonne

newStartYear<-1865; prevStartYear<-1900 #prevstartyear is the year in the historic catch ts files

allBoxes<-0:29; dynBoxes<-1:24

##set up zero catches for the burnin period

burnin_nyears<-prevStartYear-newStartYear
burnin_years<-seq(newStartYear,newStartYear+burnin_nyears)
burnin_months<-seq(0,burnin_nyears*12)
#set up timestep values in seconds from start day. 
numSecsPerMonth<-60*60*24*30.5
burnin_seconds<-burnin_months*numSecsPerMonth; n_burninSecs<-length(burnin_seconds)
#get seconds to add to the original catch ts times
secs2add<-max(burnin_seconds)+numSecsPerMonth

baseFolder<-paste(catchPath,"..\\catchts_",version,"\\",sep="")
dir.create(baseFolder)
  
  #read in one of the historic catch ts files to edit for the header lines in the new ones
  thisTSfile<-paste(catchPath,"catch1.ts",sep="")
  thisTempLines<-readLines(thisTSfile)
  thisTSlines<-thisTempLines[grep("#", thisTempLines, invert = TRUE)]
  
  newTSlines<-thisTempLines
  #replace start year
  x<-grep("seconds since",newTSlines)
  newTSlines[x]<-gsub(prevStartYear,newStartYear,newTSlines[x])
  #only want to keep the bits that start with # as replacing the other lines
  newTSlines<-newTSlines[grep("^#",newTSlines)]
  
  ##add in the new burnin lines
  for(s in 1:n_burninSecs){
    newTSlines[length(newTSlines)+1]<-paste(c(burnin_seconds[s], rep(0,33)), collapse=" ")
  }
  
  #now loop through the boxes, grab the catch value lines and add these with new time values to the burnin lines
  #then write them out to the new folder
  for(b in allBoxes){
    prevBoxFileName<-paste(catchPath, "catch",b,".ts", sep="")
    this_tsFileName<-paste(baseFolder,"catch",b,".ts",sep="")
    if(b==0){
      this_tsFileName<-paste(baseFolder,"boundary.ts",sep="")
      prevBoxFileName<-paste(catchPath, "boundary.ts", sep="")
    }
    prevBoxLines<-readLines(prevBoxFileName)
    keepBoxLines<-prevBoxLines[grep("#", prevBoxLines, invert = TRUE)]; nkl<-length(keepBoxLines)
    newBoxLines<-newTSlines
    for(k in 1:nkl){
      tempLine<-keepBoxLines[k]; tempNumbers<-get_first_number(tempLine, n="all")
      tempNumbers[1]<-tempNumbers[1]+secs2add
      newLine<-paste(tempNumbers, collapse=" ")
      newBoxLines[length(newBoxLines)+1]<-newLine
    }
    
    
    writeLines(newBoxLines,this_tsFileName)
  }

