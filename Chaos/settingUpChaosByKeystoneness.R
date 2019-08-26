####### didn't really end up using this as shifting all by the same amount doesn't really test much..
##
## similar to settingUpTestingChaos.R, but uses a scalar on the initial conditions based on how well informed we defined the group in the first paper
## with the option of only changes the top XX groups by keystoneness, or all of them
source(paste(DIR$`General functions`,"getVolDepth.R",sep=""))
source(paste(DIR$`General functions`,"read_boxes.R",sep=""))

## fix keystonenss csv
temp<- read.csv(paste(DIR$'Tables', "Keystoneness_baseModel.csv", sep=""),header=FALSE); nk<-dim(temp)[1]
key_df<-array(NA, dim=c(nk, 2))
for(k in 1:nk){
  if(k==1){
    thisNum <-1
    thisGroup <- unlist(str_split(temp[k,]," "))[2]
  } else{
    thisNum<-get_first_number(temp[k,]); thisGroup<-str_trim(gsub(thisNum, "", temp[k,]), side="both")
  }
  xx<-unlist(str_split(thisGroup," ")); 
  xxx<-paste(xx, collapse="_")
  xxxx<-grep(xxx,groupsDF$Name); thisCode<-as.character(groupsDF$Code[xxxx])
  key_df[k,]<-c(thisNum, thisCode)
}



basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

##need ##
existingICfile<-"CRAM_input_fromBase50yr"; 

chaosDirections <- c("up", "down"); chaosShifts <- c(0.05, 0.1, 0.2, 0.5); nd <- length(chaosDirections); ns <- length(chaosShifts)

for(d in 1:nd){
  chaosDirection <- chaosDirections[d]
  for(s in 1:ns){
    chaosShift<-chaosShifts[s]
    
    # 
    # chaosDirection<-"up"; 
    # chaosDirection<-"down"
    # chaosShift <- 0.05
    # chaosShift <- 0.1
    # chaosShift <- 0.2
    # chaosShift <- 0.5
    
    if(chaosDirection=="up"){
      thisScale <- 1 + chaosShift
    } else{
      thisScale <- 1- chaosShift
    }
    
    ## need ## - call it what you like
    newICFile<-paste(existingICfile, chaosDirection,gsub("\\.", "", chaosShift),sep="")
    
    ### testing testing ###
    # read this in to check how it went - prob with 'up'?
    # newICdata <- nc_open(paste(basePath,newICFile, ".nc", sep=""))
    # thisCode<-"ASQ"; thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode], side="both")
    # thisTracer<-paste(thisName,"1_Nums", sep="")
    # thisICdata <- ncvar_get(newICdata, thisTracer)
    # testExistingData <- nc_open(paste(basePath, existingICfile, ".nc", sep=""))
    # thisExistingData <- ncvar_get(testExistingData, thisTracer)
    
    
    
    if(!file.exists(paste(basePath,newICFile,sep=""))){
      ## need this
      file.copy(from=paste(basePath,existingICfile,".nc",sep=""),to=paste(basePath,newICFile,".nc",sep=""))
    }
    ## need this
    init <- nc_open(paste(basePath,newICFile,".nc",sep=""), write = T) 
    # 
    # bgmf<-"CHAT30_aea.bgm"
    # bgmFile<-paste(basePath,"..\\",bgmf,sep="")
    # 
    # depthLayers<-c(1050,700,500,300,100)
    # num_wc<-length(depthLayers)
    # 
    # volAndSuch<-getVolDepth(bgmFile,depthLayers)
    
    nlayers=6
    
    # ndynBoxes<-24
    nboxes<-30
    # boundBoxes<-c(1,seq(26,30)) #labels are 1 less than these (0,25,...,29)
    # dynBoxes<-seq(2,25) #labels are 1 less than these 
    
    Fix_negs<-function(x){
      y=x
      if(x<0){
        y<-0
      }
      return(y)
    }
    
    # Get info from init 
    var_names_init <- names(init$var) 
    
    # vars <- var_names_init[is.element(var_names_init, var_names_out)] 
    vars<-unique(sort(var_names_init))
    skipVars<-c("nominal_dz", "volume", "dz", "numlayers", "topk")
    vars<-vars[!(vars %in% skipVars)]
    
     thisNames <- groupsDF$Name[groupsDF$GroupType=="Vert"]
     thisNumTracers <- paste(thisNames, 1,"_Nums", sep="")
     
     
    for(i in seq_along(vars)){ 
      dataTemp<-ncvar_get(init,vars[i])
      
      newData <- dataTemp * thisScale
      if(length(dim(dataTemp))==3){
        dataTemp[,,1]<-newData[,,1]; dataTemp[,,2]<-"_"
      } else{
        dataTemp[,1]<-newData[,1]; dataTemp[,2]<-"_"
      }
      ncvar_put(init,varid = vars[i],vals = dataTemp)
    } 
    # close the file. 
    nc_close(init)
  }
}


## create a .bat file to dump all inton text files so can fix the dimensions
batFile <- paste(basePath, "dumpChaosInput2text.bat", sep="")
cat("##\n", file=batFile, append=FALSE)
for(d in 1:nd){
  chaosDirection <- chaosDirections[d]
  for(s in 1:ns){
    chaosShift<-chaosShifts[s]
    thisNCfile <- paste(existingICfile, chaosDirection,gsub("\\.", "", chaosShift),sep="")
    
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
## create a .bat file to dump all inton text files so can fix the dimensions
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




