createSIMInteractionSensitivityRuns<-function(thisRun,runLength="", runFile, thisPath,version="",editGroup, runPrmFile="CRAM_base_run_short", biolFileName, ICfile="CRAM_input.nc", forceFile="CRAM_force"){
  
  source(paste(DIR$'General functions',"get_first_number.R",sep=""))
  source(paste(DIR$'General functions',"get_biol_vector.R",sep=""))

  ## NOTE this function only takes in one group
  if(length(editGroup)>1){cat("Warning: only using first group"); editGroup<-editGroup[1]}
  
    #read in the biol file
  biolLines<-readLines(paste(thisPath,biolFileName,".prm",sep=""))
  #get recommended and current mLa
  thisVar<-paste(editGroup,"_mL", sep=""); x<-grep(thisVar, biolLines)
  hashlines<-grep("#", biolLines); thisHashLine<-min(hashlines[hashlines>x])
  this_mLa<-get_first_number(biolLines[(x+1)], n="all")[2]
  this_recom_mLa<-get_first_number(biolLines[thisHashLine], n="all")[2]
  if(this_mLa>0){
    thisMin<-max(0,(this_mLa - this_recom_mLa/2)); thisMax<-min(0.5, this_mLa + this_recom_mLa/2)
    changes<-c(thisMin, this_mLa, thisMax)
  } else{
    changes<-c(0, this_recom_mLa/2, this_recom_mLa)
  }
  
  for(p in 1:(length(changes))){
    thisV<-changes[p]
    #take a copy of the original biolLines
    thisBiolLines<-biolLines
    thisGroup<-editGroup
    # thisVar<-str_replace(editVar,"XXX",thisGroup)
    thisVar<-gsub("XXX",thisGroup,editVar); thisVar<-paste("^", thisVar, sep="")
    thisLineIndex<-grep(thisVar,thisBiolLines)
    thisLine<-thisBiolLines[thisLineIndex+1]
    thisBaseValue<-get_first_number(thisLine,n="all")
    vector_length<-length(thisBaseValue)
    thisNewValue<-thisBaseValue
    thisNewValue[2]<-thisV
    thisNewLine<-paste(c(thisNewValue,"\n"),collapse = " ")
    #write over the old line
    thisBiolLines[thisLineIndex+1]<-thisNewLine  
    #write out to a new version of biol.prm file
    thisBiolFile<-paste(thisPath,biolFileName,editVar,version,editGroup,p,".prm",sep="")
    writeLines(thisBiolLines,thisBiolFile)
    runText<-paste("\nWD=\"$(pwd)\"\nRUN=\"../../bin/bin/atlantisMerged -i ", ICfile," 0 -o output.nc -r ", runPrmFile, ".prm -f inputs/",forceFile,".prm -p inputs/CRAM_physics.prm -b CRAM_BH_hybrid_biol",editVar,version,editGroup,p,".prm -s CRAM_Groups.csv -q inputs/CRAM_Fisheries.csv -d output",editVar,version,editGroup,p,"\"\necho $RUN > RUN\nCMD=\"msub -l nodes=1 -l walltime=50:00:00 -l partition=slurm -l qos=standby -p -1000 -q large -o CRAM",editVar,p,".log.%j -e CRAM",editVar,version,p,".err.%j -S /bin/bash RUN\"\necho \"Running Atlantis ",editVar,version,p," for CRAM on MOAB in directory:\" $WD\necho -n \"Job started at: \" ; date\necho $RUN\nCOMMAND=\"cd $WD ; $CMD\"\nssh turbine $COMMAND\nsleep 0.5",sep="")
    cat(runText,file=runFile,append=TRUE)
  }
}


