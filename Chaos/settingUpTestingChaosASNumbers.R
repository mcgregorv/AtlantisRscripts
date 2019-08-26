#if we take the initial conditions as pre-burnt (35 years), how sensitive is our model to starting it slightly off this? Start with shifts of 5%
#read in the new initial conditions file that is 35 years into the base model
## shift these values all up (and then all down) by 5% and create new initial conditions file
#### NOTE: for age-structured groups, vary numbers, but no weights
### also, leave FLAT, REEF, SOFT unchanged

source(paste(DIR$`General functions`,"getVolDepth.R",sep=""))
source(paste(DIR$`General functions`,"read_boxes.R",sep=""))

version<- ""
basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

existingICfile<-"CRAM_input_fromBase50yr"; 

chaosDirections <- c("up", "down"); chaosShifts <- c(0.05, 0.1, 0.2, 0.5); nd <- length(chaosDirections); ns <- length(chaosShifts)

for(d in 1:nd){
  chaosDirection <- chaosDirections[d]
  for(s in 1:ns){
    chaosShift<-chaosShifts[s]
    
    if(chaosDirection=="up"){
      thisScale <- 1 + chaosShift
    } else{
      thisScale <- 1- chaosShift
    }
    newICFile<-paste(existingICfile,version, chaosDirection,gsub("\\.", "", chaosShift),sep="")
   
    if(!file.exists(paste(basePath,newICFile,sep=""))){
      file.copy(from=paste(basePath,existingICfile,".nc",sep=""),to=paste(basePath,newICFile,".nc",sep=""))
    }
    
    init <- nc_open(paste(basePath,newICFile,".nc",sep=""), write = T) 
    nlayers=6
    nboxes<-30
 
    Fix_negs<-function(x){
      y=x
      if(x<0){
        y<-0
      }
      return(y)
    }
    
    # Get info from init 
    var_names_init <- names(init$var) 
    
    ## what other vars are not species tracers..?
    # storeTestVars <-c()
    # for(i in seq_along(vars)){
    #   thisVar <- vars[i]
    #   ## if there is Nums, ResN, StructN, or _N then it is a species tracer - ignore for this test
    #   test<-grep("Nums|ResN|StructN|_N", thisVar)
    #   if(length(test)==0){storeTestVars <- c(storeTestVars,thisVar)}
    # }
    
    vars<-unique(sort(var_names_init))
    # skipVars<-c("nominal_dz", "volume", "dz", "numlayers", "topk")
    skipVars <- c("canyon", "Chl_a", "CO2", "Denitrifiction",  "DiagNFlux", "DiagNGain", "DiagNLoss", "DON", "eddy", "eflux", "erosion_rate", "Filter_Other_Cover", 
                  "flat", "hdsink", "hdsource", "Light", "Light_Adaptn_DF", "Light_Adaptn_MB", "Light_Adaptn_PL", "Light_Adaptn_PS", "Macroalgae_Cover", "MicroNut", 
                  "MicroPB_Cover",   "Nitrification",  "porosity", "reef", "salt", "sedbiodens", "sedbiodepth", "seddetdepth", "sedirrigenh", "sedoxdepth", "sedturbenh", 
                  "soft", "Stress", "Temp", "vflux", "water", "nominal_dz", "volume", "dz", "numlayers", "topk", "Oxygen", "NH3", "NO3", "Si", "Det_Si")
    vars<-vars[!(vars %in% skipVars)]
    
    for(i in seq_along(vars)){ 
      thisVar <- vars[i]
      ## if thisVar has ResN or StructN, skip as don't want to change weights
      test<-grep("ResN|StructN|soft|reef|flat", thisVar)
      if(length(test)==0){
        dataTemp<-ncvar_get(init,vars[i])
        
        newData <- dataTemp * thisScale
        if(length(dim(dataTemp))==3){
          dataTemp[,,1]<-newData[,,1]; dataTemp[,,2]<-"_"
        } else{
          dataTemp[,1]<-newData[,1]; dataTemp[,2]<-"_"
        }
        ncvar_put(init,varid = vars[i],vals = dataTemp)
      }
    } 
    # close the file. 
    nc_close(init)
  }
}

## create a .bat file to dump all inton text files so can fix the dimensions - run this file from command prompt once have created it
batFile <- paste(basePath,version, "dumpChaosInput2text.bat", sep="")
cat("##\n", file=batFile, append=FALSE)
for(d in 1:nd){
  chaosDirection <- chaosDirections[d]
  for(s in 1:ns){
    chaosShift<-chaosShifts[s]
    thisNCfile <- paste(existingICfile, version, chaosDirection,gsub("\\.", "", chaosShift),sep="")
    
    thisLine <- paste("ncdump ", thisNCfile, ".nc > ", thisNCfile, ".txt \n", sep="")
    cat(thisLine, file=batFile, append=TRUE)
  }
}

## fix the dimensions
for(d in 1:nd){
  chaosDirection <- chaosDirections[d]
  for(s in 1:ns){
    chaosShift<-chaosShifts[s]
    thisfile<-thisNCfile <- paste(basePath,existingICfile, chaosDirection,gsub("\\.", "", chaosShift), ".txt",sep="")
    
    newfile<-thisfile
    thisLines<-readLines(thisfile)
    
    x <- thisLines[grep("_,",thisLines)[1]]
    
    ## only edit lines after 'data:'
    lineStartIndex<-grep("data:", thisLines)
    
    newLines<-thisLines
    newLines[lineStartIndex:length(newLines)] <- unlist(lapply(thisLines[lineStartIndex:length(newLines)], FUN=function(x){str_trim(gsub("_,|;|NaN,|NaN","",x), side="both")}))
    newLines[grep("_,",thisLines)[1]]
    
    ## add in the end ;
    index<-grep(",",newLines); ni<-length(index)
    for(i in 1:ni){
      this_i <- index[i]
      nextLine<-newLines[(this_i+1)]
      if(nextLine==""){
        thisLine<-newLines[this_i]
        temp <- gsub(", ", " TEMP ", thisLine); temp2<-gsub(",",";", temp);
        if(length(grep(";", temp2))==0){ temp2 <- paste(temp, ";", collapse="", sep="")} #if there is no , at the end of the line - then just add the ; to the end
        thisNewLine<-gsub(" TEMP", ",", temp2)
        newLines[this_i]<-thisNewLine
      }
    }
    
    ## check for lone underscores left
    index <- newLines=="_"
    newLines <- newLines[!index]
    ## might as well take out white space
    index <- newLines==""
    newLines <- newLines[!index]
    
    ## a couple of ad-hoc fixes
    x <- grep("^t =", newLines)
    newLines[x]<- "t = 0;"
    
    writeLines(newLines, newfile)
    
  }
}

## file to turn back into .nc files
## create a .bat file to dump all inton text files so can fix the dimensions - run this file from command prompt once have created it
batFile <- paste(basePath, "dumpChaosText2input.bat", sep="")
cat("##\n", file=batFile, append=FALSE)
for(d in 1:nd){
  chaosDirection <- chaosDirections[d]
  for(s in 1:ns){
    chaosShift<-chaosShifts[s]
    thisNCfile <- paste(existingICfile, chaosDirection,gsub("\\.", "", chaosShift),sep="")
    
    thisLine <- paste("ncgen -o ", thisNCfile, ".nc  ", thisNCfile, ".txt \n", sep="")
    cat(thisLine, file=batFile, append=TRUE)
  }
}


## set up the run file
baseICfile<- "CRAM_input"
baseRunCommand <- "../../bin/bin/atlantisMerged -i CRAM_input.nc 0 -o output.nc -r CRAM_base_run.prm -f inputs/CRAM_forceBURNIN1865.prm -p inputs/CRAM_physics.prm -b CRAM_BH_hybrid_biol.prm -s CRAM_Groups.csv -q CRAM_Fisheries.csv -d outputFolder"
fishRunCommand <- "../../bin/bin/atlantisMerged -i CRAM_input.nc 0 -o output.nc -r CRAM_baseFish_run.prm -f inputs/CRAM_forceBURNIN1865.prm -p inputs/CRAM_physics.prm -b CRAM_BH_hybrid_biol.prm  -h CRAM_harvest_short.prm -s CRAM_Groups.csv -q CRAM_Fisheries.csv -d outputFolder"
runFile<-paste(basePath, "RunChaos", sep="")
cat("#Run base run and historical catches removed run with scaled ICs", file=runFile, append=FALSE)
for(d in 1:nd){
  chaosDirection <- chaosDirections[d]
  for(s in 1:ns){
    chaosShift<-chaosShifts[s]
    thisNCfile <- paste(existingICfile, chaosDirection,gsub("\\.", "", chaosShift),sep="")
    thisRunCommand <- gsub(baseICfile, thisNCfile, baseRunCommand)
    thisRunCommand <- gsub("outputFolder", paste("outputChaos",chaosDirection,gsub("\\.", "", chaosShift), sep="" ), thisRunCommand)
    cat(paste("WD=\"$(pwd)\"
              RUN=\"", thisRunCommand ,"\"
              echo $RUN > RUN
              CMD=\"msub -l nodes=1 -l walltime=50:00:00 -l partition=slurm -l qos=standby -p -1000 -q large -o CRAMBase1.log.%j -e CRAMBase.err.%j -S /bin/bash RUN\"
              echo \"Running Atlantis Base for CRAMBase on MOAB in directory:\" $WD
              echo -n \"Job started at: \" ; date
              echo $RUN
              COMMAND=\"cd $WD ; $CMD\"
              ssh turbine $COMMAND
              sleep 1"), file=runFile, append=TRUE)
    ## this fishing version
    thisRunCommand <- gsub(baseICfile, thisNCfile, fishRunCommand)
    thisRunCommand <- gsub("outputFolder", paste("outputChaosFISH",chaosDirection,gsub("\\.", "", chaosShift), sep="" ), thisRunCommand)
    cat(paste("
              WD=\"$(pwd)\"
              RUN=\"", thisRunCommand ,"\"
              echo $RUN > RUN
              CMD=\"msub -l nodes=1 -l walltime=50:00:00 -l partition=slurm -l qos=standby -p -1000 -q large -o CRAMBase1.log.%j -e CRAMBase.err.%j -S /bin/bash RUN\"
              echo \"Running Atlantis Base for CRAMBase on MOAB in directory:\" $WD
              echo -n \"Job started at: \" ; date
              echo $RUN
              COMMAND=\"cd $WD ; $CMD\"
              ssh turbine $COMMAND
              sleep 1
              "), file=runFile, append=TRUE)
    
  }
}
## also add in the base runs - this uses the new base IC file
thisNCfile<- "CRAM_input_fromBase50yr"
thisRunCommand <- gsub(baseICfile, thisNCfile, baseRunCommand)
thisRunCommand <- gsub("outputFolder", "outputChaosBASE", thisRunCommand)
cat(paste("WD=\"$(pwd)\"
          RUN=\"", thisRunCommand ,"\"
          echo $RUN > RUN
          CMD=\"msub -l nodes=1 -l walltime=50:00:00 -l partition=slurm -l qos=standby -p -1000 -q large -o CRAMBase1.log.%j -e CRAMBase.err.%j -S /bin/bash RUN\"
          echo \"Running Atlantis Base for CRAMBase on MOAB in directory:\" $WD
          echo -n \"Job started at: \" ; date
          echo $RUN
          COMMAND=\"cd $WD ; $CMD\"
          ssh turbine $COMMAND
          sleep 1"), file=runFile, append=TRUE)
## this fishing version
thisRunCommand <- gsub(baseICfile, thisNCfile, fishRunCommand)
thisRunCommand <- gsub("outputFolder", "outputChaosFISH", thisRunCommand)
cat(paste("
          WD=\"$(pwd)\"
          RUN=\"", thisRunCommand ,"\"
          echo $RUN > RUN
          CMD=\"msub -l nodes=1 -l walltime=50:00:00 -l partition=slurm -l qos=standby -p -1000 -q large -o CRAMBase1.log.%j -e CRAMBase.err.%j -S /bin/bash RUN\"
          echo \"Running Atlantis Base for CRAMBase on MOAB in directory:\" $WD
          echo -n \"Job started at: \" ; date
          echo $RUN
          COMMAND=\"cd $WD ; $CMD\"
          ssh turbine $COMMAND
          sleep 1
          "), file=runFile, append=TRUE)




