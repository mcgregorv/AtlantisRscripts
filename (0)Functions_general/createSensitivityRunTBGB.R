createSensitivityRunTBGB<-function(thisRun,runLength="",thisPath,editVar,version="",editGroups, runPrmFile="TBGB_base_run_short",changes,proportionalChange,nextLine,isVector=FALSE,changeVectorIndex=NULL,groupsFileName,biolFileName,varMax,ICfile="CRAM_input.nc", forceFile="CRAM_force"){

  source(paste(DIR$'General functions',"get_first_number.R",sep=""))
  source(paste(DIR$'General functions',"get_biol_vector.R",sep=""))
  
  # thisRun<-"base"
  # thisPath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
  # outPath<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
  # groupsFileName<-"CRAM_Groups"
  # biolFileName<-"CRAM_base_biol"
  # editVar<-"XXX_mL"
  # nextLine<-TRUE
  # version=""
  # proportionalChange=TRUE
  # changes<-c(0,1e-2,1e-1,1,1e+1,1e+2)
  # isVector=TRUE
  # chageVectorIndex<-c(1,1)
  # runLength=""
  # varMax<-0.5

  #read in outputDietCheck.txt
  # options(stringsAsFactors = FALSE)
  # 
  # thisRun<-"TBGBLinuxBase"
  # runLength<-"_short" #if empty, selects default run (which is long and has 20 year burn in). if _short selects short run which is 20 years and no burn in
  # runLength<-""
  # 
  # thisPath<-paste(DIR$'Base',"\\ATLANTIS\\Model_runs\\TBGB\\",thisRun,"\\",sep="")
  # # plotPath<-paste(DIR$'Base',"\\ATLANTIS\\Model_runs\\testing\\figures",this_out,"\\",sep="")
  # 
  # #variables to edit
  # editVar<-"KDENR_XXX"
  # version<-"D"  #no version had the original linear natural mortality. B had zero linear M for all SCA and co, C has 0 for most, some for IVS and SCL
  # #version D has CEP KDENR fixed at 5, SCL, IVS KDENR also fixed at 5, the rest with changes higher than for prevous runs
  # #Version E is the same as D, but has a 20 year burn-in then does the full run (if it can!)
  # #groups to edit
  # # editGroups<-ageInverts #for versions "", "B", "C"
  # editGroups<-sort(c("IVH","MUS","SCA","OYS"))
  # #changes them all at the same time, to this value
  # # changes<-c(1e-1,5e-1,1,2e+0,5e+0,1e+1) #these were used for versions "", "B", "C"
  # changes<-c(1e+1,1.5e+1,2.0e+1,3.0e+1,5.0e+1) ## used for version "D"
  # proportionalChange<-FALSE #set to true if want the existing values to be multiplied by the change values, FALSE if the change values replace the existing ones
  # nextLine<-TRUE #some paramters are on the same line as the parameter name (in which case set this to FALSE), some are the next line
  # 
  
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
        thisNewLine<-paste(gsub("\\^","",thisVar),"\t",thisNewValue,"\n",sep="")
        #write over the old line
        thisBiolLines[thisLineIndex]<-thisNewLine  
      }
      
    }
    #write out to a new version of biol.prm file
    thisBiolFile<-paste(thisPath,biolFileName,editVar,version,p,".prm",sep="")
    writeLines(thisBiolLines,thisBiolFile)
    #add to run file
    # 
    # "atlantismain -i CRAM_input.nc 0 -o output.nc -r CRAM_base_run.prm -f inputs\CRAM_force.prm -p inputs\CRAM_physics.prm -b CRAM_base_biol.prm -s CRAM_Groups.csv -q inputs\CRAM_Fisheries.csv -d base\output"
    # 
    # RUN="../bin/bin/atlantisMerged -i create_biol_input/TBGB_input.nc 0 -o output.nc -r TBGB_run.prm -f TBGB_force.prm -p TBGB_physics.prm -b TBGB_biolXXX_mL1.prm -h TBGB_harvest.prm -s TBGB_Groups.csv -q TBGB_Fisheries.csv -d outputXXX_mL1"
    # 
    runText<-paste("\nWD=\"$(pwd)\"\nRUN=\"../bin/bin/atlantisMerged -i ", ICfile," 0 -o output.nc -r ", runPrmFile, ".prm -f ",forceFile,".prm -p TBGB_physics.prm -b TBGB_biol",editVar,version,p,".prm -s TBGB_Groups.csv -q inputs/TBGB_Fisheries.csv -d output",editVar,version,p,"\"\necho $RUN > RUN\nCMD=\"msub -l nodes=1 -l walltime=50:00:00 -l partition=slurm -l qos=standby -p -1000 -q large -o TBGB",editVar,p,".log.%j -e TBGB",editVar,version,p,".err.%j -S /bin/bash RUN\"\necho \"Running Atlantis ",editVar,version,p," for TBGB on MOAB in directory:\" $WD\necho -n \"Job started at: \" ; date\necho $RUN\nCOMMAND=\"cd $WD ; $CMD\"\nssh turbine $COMMAND\nsleep 0.5",sep="")
    
    # 
    # runText<-paste("\nWD=\"$(pwd)\"\nRUN=\"../bin/bin/atlantisMerged -i create_biol_input/TBGB_input.nc 0 -o output.nc -r TBGB_run",runLength,".prm -f TBGB_force.prm -p TBGB_physics.prm -b TBGB_biol",editVar,version,p,".prm -h TBGB_harvest.prm -s TBGB_Groups.csv -q TBGB_Fisheries.csv -d output",editVar,version,p,"\"\necho $RUN > RUN\nCMD=\"msub -l nodes=1 -l walltime=50:00:00 -l partition=slurm -l qos=standby -p -1000 -q large -o TBGB",editVar,p,".log.%j -e TBGB",editVar,version,p,".err.%j -S /bin/bash RUN\"\necho \"Running Atlantis ",editVar,version,p," for TBGB on MOAB in directory:\" $WD\necho -n \"Job started at: \" ; date\necho $RUN\nCOMMAND=\"cd $WD ; $CMD\"\nssh turbine $COMMAND\nsleep 0.5",sep="")
    cat(runText,file=runFile,append=TRUE)
    #   cat("\n\n",file=runFile,append=TRUE)
  }
}