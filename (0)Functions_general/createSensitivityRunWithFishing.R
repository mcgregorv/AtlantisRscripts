createSensitivityRunWithFishing<-function(thisRun,runLength="",thisPath,editVar,version="",editGroups,changes,proportionalChange,nextLine,isVector=FALSE,changeVectorIndex=NULL,groupsFileName,biolFileName,varMax,ICfile="CRAM_input.nc", forceFile="CRAM_force"){
  
  source(paste(DIR$'General functions',"get_first_number.R",sep=""))
  source(paste(DIR$'General functions',"get_biol_vector.R",sep=""))
  
  #these might be handy
  groups<-read.csv(paste(thisPath,groupsFileName,".csv",sep=""))
  vertsIndex<-groups$NumCohorts>=2
  verts<-groups$Code[vertsIndex]
  inverts<-groups$Code[!vertsIndex]
  # ageInverts<-groups$Code[groups$Code %in% inverts & groups$NumCohorts>1]
  
  if(!exists("editGroups")){editGroups<-groups$Code}
  
  #read in the biol file
  # biolFileName<-"TBGB_biol"
  biolLines<-readLines(paste(thisPath,biolFileName,".prm",sep=""))
  
  runFile<-paste(thisPath,"RunMulti_",editVar,version,"",sep="")
  runText<-paste("#Run multiple runs with varying ",editVar," by ",ifelse(proportionalChange,"proportional changes","changes"),sep="")
  cat(runText,file=runFile,append=FALSE)
  cat(paste("\n#",changes,sep=""),file=runFile,append=TRUE)
  
  for(p in 1:(length(changes))){
    thisV<-changes[p]
    #take a copy of the original biolLines
    thisBiolLines<-biolLines
    for(g in 1:(length(editGroups))){
      thisGroup<-editGroups[g]
      # thisVar<-str_replace(editVar,"XXX",thisGroup)
      thisVar<-gsub("XXX",thisGroup,editVar); thisVar<-paste("^", thisVar, sep="")
      thisLineIndex<-grep(thisVar,thisBiolLines)
      if(nextLine){thisLine<-thisBiolLines[thisLineIndex+1]}else{thisLine<-thisBiolLines[thisLineIndex]}   
      if(isVector){
        # if(is.null(changeVectorIndex)){
        thisBaseValue<-get_first_number(thisLine,n="all")
        vector_length<-length(thisBaseValue)
        # }else{
        #   vector_length<-length(changeVectorIndex)
        #   thisBaseValue<-get_first_number(thisLine,n="all")[changeVectorIndex]
        # }
      }else{
        thisBaseValue<-get_first_number(thisLine)
      }
      if(proportionalChange){
        if(isVector){
          thisNewValue<-thisBaseValue
          if(is.null(changeVectorIndex)){
            thisNewValue<-unlist(lapply(thisNewValue,FUN=function(x){min(x*thisV,varMax)}))
          }else{
            thisNewValue[changeVectorIndex]<-unlist(lapply(thisNewValue[changeVectorIndex]*thisV,FUN=function(x){min(varMax,x)}))
          }
        }else{
          thisNewValue<-min(thisBaseValue*thisV,varMax)
        }
      } else{
        if(isVector){
          thisNewValue<-thisBaseValue
          thisNewValue[changeVectorIndex]<-thisV
        }else{
          thisNewValue<-thisV
        }
      }
      if(nextLine){
        thisNewLine<-paste(c(thisNewValue,"\n"),collapse = " ")
        #write over the old line
        thisBiolLines[thisLineIndex+1]<-thisNewLine  
      }else{
        thisNewLine<-paste(thisVar,"\t",thisNewValue,"\n",sep="")
        #write over the old line
        thisBiolLines[thisLineIndex]<-thisNewLine  
      }
      
    }
    #write out to a new version of biol.prm file
    thisBiolFile<-paste(thisPath,biolFileName,editVar,version,p,".prm",sep="")
    writeLines(thisBiolLines,thisBiolFile)
    #add to run file
    runText<-paste("\nWD=\"$(pwd)\"\nRUN=\"../../bin/bin/atlantisMerged -i ", ICfile," 0 -o output.nc -r CRAM_base_run_short.prm -f inputs/",forceFile,".prm -p inputs/CRAM_physics.prm -b CRAM_BH_hybrid_biol",editVar,version,p,".prm -s CRAM_Groups.csv -q CRAM_Fisheries.csv -d output",editVar,version,p,"\"\necho $RUN > RUN\nCMD=\"msub -l nodes=1 -l walltime=50:00:00 -l partition=slurm -l qos=standby -p -1000 -q large -o CRAM",editVar,p,".log.%j -e CRAM",editVar,version,p,".err.%j -S /bin/bash RUN\"\necho \"Running Atlantis ",editVar,version,p," for CRAM on MOAB in directory:\" $WD\necho -n \"Job started at: \" ; date\necho $RUN\nCOMMAND=\"cd $WD ; $CMD\"\nssh turbine $COMMAND\nsleep 0.5",sep="")
    cat(runText,file=runFile,append=TRUE)
    #
    ## add fishing run command
    runText<-paste("\nWD=\"$(pwd)\"\nRUN=\"../../bin/bin/atlantisMerged -i ", ICfile," 0 -o output.nc -r CRAM_baseFish_run_long.prm -f inputs/",forceFile,".prm -p inputs/CRAM_physics.prm -b CRAM_BH_hybrid_biol",editVar,version,p,".prm  -h CRAM_harvest_short.prm -s CRAM_Groups.csv -q CRAM_Fisheries.csv -d outputFISH",editVar,version,p,"\"\necho $RUN > RUN\nCMD=\"msub -l nodes=1 -l walltime=50:00:00 -l partition=slurm -l qos=standby -p -1000 -q large -o CRAM",editVar,p,".log.%j -e CRAM",editVar,version,p,".err.%j -S /bin/bash RUN\"\necho \"Running Atlantis ",editVar,version,p," for CRAM on MOAB in directory:\" $WD\necho -n \"Job started at: \" ; date\necho $RUN\nCOMMAND=\"cd $WD ; $CMD\"\nssh turbine $COMMAND\nsleep 0.5",sep="")
    cat(runText,file=runFile,append=TRUE)
  }
}