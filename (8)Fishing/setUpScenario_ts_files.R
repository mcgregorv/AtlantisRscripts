#read in last 5 year average (total per year) catches - these were created in PlotHistoricForcedCatches_part2
#scale the value from the last year (up or down depending on scenario) and proportion it spatially based on last 5 year average proportions
this_path<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
catchPath<-paste(this_path,"inputs\\catch_history\\catchts\\",sep="")

groupsDF<-read.csv(paste(this_path,"CRAM_groups.csv",sep=""))
catchGroupsDF<-groupsDF[groupsDF$IsFished==1,]; ncg<-dim(catchGroupsDF)[1]

version<-"FS_A1"

X_CN<-5.7
mg_2_tonne<-0.00000002
kg_2_mg<-1e-3/mg_2_tonne

newStartYear<-1900; prevStartYear<-1900 #prevstartyear is the year in the historic catch ts files

allBoxes<-0:29; dynBoxes<-1:24

last5YearAve<-read.csv(paste(DIR$'Tables',"CatchHist_last5YearAverage.csv",sep=""))[,1]

aveBoxProp<-read.csv(paste(DIR$'Tables',"CatchHist_last5YearPropByBox.csv",sep=""))

# thisScenario<-"Base"; thisScenarioName<-"Base"; thisScalar<-1; thisGroups<-catchGroupsDF$Code
scenarioScalars<-c(0,1,0.5,0.6,0.7,0.8,0.5,0.4,0.3,0.2,1.2,1.3,1.4,1.5,1.2,1.3,1.4,1.5); nScenarios<-length(scenarioScalars)
scenarioGroups<-c(rep("All", 6),rep("Hoki",8),rep("All",4))
scenarioCodes<-mapply(FUN=function(x,y){paste(x,y*100, "catch",sep="")},x=scenarioGroups,y=scenarioScalars)
scenarioNames<-mapply(FUN=function(x,y){paste(x," ",y*100, "% catch",sep="")},x=scenarioGroups,y=scenarioScalars)
scenarioDF<-data.frame(cbind(scenarioCodes,scenarioNames,scenarioScalars,scenarioGroups))
nScenarios<-dim(scenarioDF)[1]

scenario_nyears<-10
scenario_years<-seq(newStartYear,newStartYear+scenario_nyears)
scenario_months<-seq(0,scenario_nyears*12)
#set up timestep values in seconds from start day. 
numSecsPerMonth<-60*60*24*30.5
scenario_seconds<-scenario_months*numSecsPerMonth

for(s in 1:nScenarios){
  thisScenario<-scenarioDF$scenarioCodes[s]
  cat(thisScenario," -- ")
  thisScalar<-as.double(scenarioDF$scenarioScalars[s])
  thisGroups<-catchGroupsDF$Code
  if(scenarioDF$scenarioGroups[s] != "All"){thisGroups<-"HOK"}

  #first set up base case with 5 year average catch and 5 year average proportions
  newCatchArray<-0*aveBoxProp
  for(g in 1:ncg){
    thisCode<-catchGroupsDF$Code[g]
    newCatchArray[,g]<-aveBoxProp[,g]*last5YearAve[g]
    if(thisCode %in% thisGroups){
      newCatchArray[,g]<-aveBoxProp[,g]*last5YearAve[g]*thisScalar
    }
  }
  
  baseFolder<-paste(catchPath,"..\\catchts_",thisScenario,"\\",sep="")
  dir.create(baseFolder)
  
  #read in one of the historic catch ts files to edit for the new ones
  thisTSfile<-paste(catchPath,"catch1.ts",sep="")
  thisTempLines<-readLines(thisTSfile)
  thisTSlines<-thisTempLines[grep("#", thisTempLines, invert = TRUE)]
  
  newTSlines<-thisTempLines
  #replace start year
  x<-grep("seconds since",newTSlines)
  newTSlines[x]<-gsub(prevStartYear,newStartYear,newTSlines[x])
  #only want to keep the bits that start with # as replacing the other lines
  newTSlines<-newTSlines[grep("^#",newTSlines)]
  
  
  for(b in allBoxes){
    this_tsFileName<-paste(baseFolder,"catch",b,".ts",sep="")
    if(b==0){
      this_tsFileName<-paste(baseFolder,"boundary.ts",sep="")
    }
    if(b %in% dynBoxes){
      thisData<-newCatchArray[b,]
    } else{
      thisData<-0*newCatchArray[1,]
    }
    #turn thisData into mg N caught per second
    x<-as.double(thisData)/12
    convertedData<-(x*kg_2_mg/(X_CN))/numSecsPerMonth
    
    writeLines(newTSlines,this_tsFileName)
    for(t in 1:length(scenario_seconds)){
      cat(scenario_seconds[t],"\t",file=this_tsFileName,append=TRUE)
      cat(as.double(convertedData),file=this_tsFileName,append=TRUE)
      cat("\n",file=this_tsFileName,append=TRUE)
    }
  
  }
}

#create the run file
runFile<-paste(this_path,"base\\RunFishScenarios_",version,"",sep="")

runText<-paste("#Run multiple fish scenarios ",version,sep="")
cat(runText,file=runFile,append=FALSE)
cat(paste("\n#",changes,sep=""),file=runFile,append=TRUE)

existingForceFile<-"inputs/CRAM_force.prm"

thisInitialConditionsFile<-"CRAM_input_short.nc" #can change this to one creates from outputs of another run

for(s in 1:nScenarios){
  thisCode<-scenarioDF$scenarioCodes[s]; thisScenario<-scenarioDF$scenarioCodes[s]
  #create new forcing file, which points to the appropriate catch ts files
  thisForceFile<-paste("inputs/CRAM_force_",thisCode,".prm",sep="")
  file.copy(paste(this_path,existingForceFile,sep=""), paste(this_path,thisForceFile,sep=""), overwrite = TRUE)
  thisForceLines<-readLines(paste(this_path,thisForceFile,sep=""))
  thisCatchFolder<-paste("catchts_",thisScenario,sep="")
  newForceLines<-gsub("catchts", thisCatchFolder,thisForceLines)
  writeLines(newForceLines,paste(this_path,thisForceFile,sep=""))
  
  thisOutFolder<-paste(version,thisCode,sep="")
  
  runText<-paste("\nWD=\"$(pwd)\"\nRUN=\"../../bin/bin/atlantisMerged -i ", thisInitialConditionsFile," 0 -o output.nc -r CRAM_baseFish_run_short.prm -f ", thisForceFile, " -p inputs/CRAM_physics.prm -b CRAM_BH_hybrid_biol.prm -h CRAM_harvest_short.prm  -s CRAM_Groups.csv -q CRAM_Fisheries.csv -d base/output",thisOutFolder,"\"\necho $RUN > RUN\nCMD=\"msub -l nodes=1 -l walltime=50:00:00 -l partition=slurm -l qos=standby -p -1000 -q large -o CRAM",editVar,p,".log.%j -e CRAMfishscen",thisScenario,version,".err.%j -S /bin/bash RUN\"\necho \"Running Atlantis ",editVar,version,p," for CRAM on MOAB in directory:\" $WD\necho -n \"Job started at: \" ; date\necho $RUN\nCOMMAND=\"cd $WD ; $CMD\"\nssh turbine $COMMAND\nsleep 0.5",sep="")
  cat(runText,file=runFile,append=TRUE)
}


## this should turn all files in a folder to unix files
# find . -type f -print0 | xargs -0 dos2unix

#create run file to turn all files created here to unix
thisFile<-paste(this_path,"Turn2dos_",version,sep="")
cat("",file=thisFile,append=FALSE)
for(s in 1:nScenarios){
  thisCode<-scenarioDF$scenarioCodes[s]; 
  thisForceFile<-paste("CRAM_force_",thisCode,".prm",sep="")
  cat("dos2unix",thisForceFile,"\n",file=thisFile, append = TRUE)
}

# 
# ###### testing testing ##########################
# #read in the original ts file for box 1
# origFile<-paste(this_path,"\\inputs\\catch_history\\catchts\\catch1.ts",sep="")
# origLines<-readLines(origFile); onl<-length(origLines)
# origTest<-origLines[(onl-4):onl]
# oTest<-unlist(lapply(origTest,FUN=function(x){get_first_number(x, n="all")[-1]}))
# # oTest<-get_first_number(origLines[onl],n="all")
# 
# newFile<-paste(this_path,"inputs\\catch_history\\catchts_All100catch\\catch1.ts",sep="")
# newLines<-readLines(newFile); nnl<-length(newLines)
# newTest<-newLines[(nnl-4):nnl]
# nTest<-unlist(lapply(newTest,FUN=function(x){get_first_number(x, n="all")[-1]}))
# # nTest<-get_first_number(newLines[nnl],n="all")
