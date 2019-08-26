#the ROMS variables we have are for XX years
#this script will load the ROMS files (salinity, temperature and hydro), split into years, then create X number of new ROMS with the years bootstrapped (selected with replacement). For each new ROMS file set, a new RUN file is created with the commands stored in one file ready to copy over to Trubine. If doing 100 runs, output maybe 5-10 times per run (perhaps at years 5,10,20,30,50). Start the runs with initial conditions that have come from the outputs of a stable run.

basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

##A was the initial run, it was 20 12 year runs and it was unstable.
#B is a 5 year run, outputting every 5 days. The run is heading towards stable, but is only 5 years in and needs 150 years to reach stable, so it is not there yet
#B is intended to give insight into nutrients and primary producers. IC .nc file is CRAM_input_short_from5yr.nc

## TAKE 2: a 35 year burnin with no outputs, using the same one year of ROMS, followed by a 10 year period outputting once a year
version<-"A"

nlayer<-6; nboxes<-30

#read in the full ROMS files
ROMSpath<-paste(basePath,"inputs\\ROMS\\",sep="")

writePath<-paste(ROMSpath,"BootstrapROMS\\",version,"\\",sep="")

tempROMSFile<-paste(ROMSpath,"Chatham30_tempAll.nc",sep="")
hydroROMSFile<-paste(ROMSpath,"Chatham30_hydroAll.nc",sep="")
saltROMSFile<-paste(ROMSpath,"Chatham30_saltAll.nc",sep="")

fullHydro<-nc_open(hydroROMSFile)
fullTemp<-nc_open(tempROMSFile)
fullSalt<-nc_open(saltROMSFile)

vnames_salt<-names(fullSalt$var); vnames_temp<-names(fullTemp$var); vnames_hydro<-names(fullHydro$var)
saltData<-ncvar_get(fullSalt,"salinity"); tempData<-ncvar_get(fullTemp,"temperature")
nts<-dim(saltData)[3] #they are 12 hour time steps
hydroExchangeData<-ncvar_get(fullHydro,"exchange"); hydroDestBData<-ncvar_get(fullHydro,"dest_b"); hydroDestKData<-ncvar_get(fullHydro,"dest_k")

nyears<-nts/(365*2) #check it is a round number!
tsPerYear<-2*365
#split the data up by year
saltByYear<-NULL; tempByYear<-NULL; exchangeByYear<-NULL; destBByYear<-NULL;  destKByYear<-NULL
for(y in 1:nyears){
  startY<-(y-1)*tsPerYear+1; endY<-startY+tsPerYear-1
  # cat(paste(y, " has start year",startY,", endY ", endY," which in years are ",startY/(365*2)," ",endY/(365*2),"\n",sep=""))
  saltByYear[[y]]<-saltData[,,startY:endY]; tempByYear[[y]]<-tempData[,,startY:endY]; 
  exchangeByYear[[y]]<-hydroExchangeData[,,,startY:endY]
  destBByYear[[y]]<-hydroDestBData[,,,startY:endY]; destKByYear[[y]]<-hydroDestKData[,,,startY:endY]
}
#number of runs
nruns<-100; runLength<-20; burnin<-35 #runLength is number of years; burnin will repeat one ROMS - 2003 is closest to the mean based on temperature
baseYear<-2003; baseYearIndex<-grep(baseYear,seq(1996,2004))
yearsToSample<-seq(1,nyears)
storeRandomYears<-NULL
for(r in 1:nruns){
  storeRandomYears[[r]]<-sample(yearsToSample,replace=TRUE,size=runLength)
}

#################
##write out the years selected for each run
recordRunsFile<-paste(writePath,"RecordedYearsSelected.txt",sep="")
cat("Years selected for each run",file=recordRunsFile,append=FALSE)
for(r in 1:nruns){
  cat(paste("\nRun ",r,"\n",sep=""),file=recordRunsFile,append=TRUE)
  cat(storeRandomYears[[r]],file=recordRunsFile,append=TRUE)
}

#for each run, create the ROMS files
emptyArray<-array(NA,dim=c(nlayer,nboxes,(burnin+runLength)*365*2))
emptyArrayHydro<-array(NA,dim=c(21,5,30,(burnin+runLength)*365*2))
timeData<-seq(0,(burnin+runLength)*730-1,by=0.5)*(12*60*60)
#if more years than in existing ROMS, make multiple ones and add them to force.prm file
#how many files do we need?
nfiles<-ceiling((burnin+runLength)/nyears)
#how many years go into the last file?
lastFileYears<-(burnin+runLength)-trunc((burnin+runLength)/nyears)*nyears

for(r in 1:nruns){
  testFile<-paste(writePath,"Chatham_temp",r,"_1.nc",sep="")
  if(file.exists(testFile)){
    cat("ROMS files seem to exist - delete first to ensure we don't write over something that should be kept\n")
  }else{
    thisYearSelection<-c(rep(baseYearIndex,burnin),storeRandomYears[[r]])
    thisSaltData<-emptyArray; thisTempData<-emptyArray; thisExchangeData<-emptyArrayHydro; thisDestBData<-emptyArrayHydro; thisDestKData<-emptyArrayHydro
    thisDestB<-emptyArrayHydro; thisDestK<-emptyArrayHydro
    for(y in 1:runLength){
      thisSelectedYear<-thisYearSelection[y]
      startY<-(y-1)*tsPerYear+1; endY<-startY+tsPerYear-1
      thisSaltData[,,startY:endY]<-saltByYear[[thisSelectedYear]]
      thisTempData[,,startY:endY]<-tempByYear[[thisSelectedYear]]
      thisExchangeData[,,,startY:endY]<-exchangeByYear[[thisSelectedYear]]
      thisDestB[,,,startY:endY]<-destBByYear[[thisSelectedYear]]
      thisDestK[,,,startY:endY]<-destKByYear[[thisSelectedYear]]
    }
    for(f in 1:nfiles){
      ##TEMPERATURE FILES
      newTempNCfile<-paste(writePath,"Chatham_temp",r,"_",f,".nc",sep="")
      file.copy(from=tempROMSFile,to=newTempNCfile,overwrite=TRUE)
      newTempNC<-nc_open(newTempNCfile,write=TRUE)
      startY<-(f-1)*nts+1; endY<-min(startY+nts-1,dim(thisTempData)[3])
      thisValsData<-thisTempData[,,startY:endY]
      if(startY+nts-1>dim(thisValsData)[3]){
        valsData<-0*thisTempData[,,1:nts]
        valsData[,,1:dim(thisValsData)[3]]<-thisValsData
      }else{
        valsData<-thisValsData
      }
      ncvar_put(newTempNC,varid = "temperature",vals = valsData)
      nc_close(newTempNC)
      ##SALINITY FILES
      newSaltNCfile<-paste(writePath,"Chatham_salt",r,"_",f,".nc",sep="")
      file.copy(from=saltROMSFile,to=newSaltNCfile,overwrite=TRUE)
      newSaltNC<-nc_open(newSaltNCfile,write=TRUE)
      startY<-(f-1)*nts+1; endY<-min(startY+nts-1,dim(thisSaltData)[3])
      thisValsData<-thisSaltData[,,startY:endY]
      if(startY+nts-1>dim(thisValsData)[3]){
        valsData<-0*thisSaltData[,,1:nts]
        valsData[,,1:dim(thisValsData)[3]]<-thisValsData
      }else{
        valsData<-thisValsData
      }
      ncvar_put(newSaltNC,varid = "salinity",vals = valsData)
      nc_close(newSaltNC)
      ##HYDRO FILES vars "exchange" "dest_b"   "dest_k"  
      newHydroNCfile<-paste(writePath,"Chatham_hydro",r,"_",f,".nc",sep="")
      file.copy(from=hydroROMSFile,to=newHydroNCfile,overwrite=TRUE)
      newHydroNC<-nc_open(newHydroNCfile,write=TRUE)
      startY<-(f-1)*nts+1; endY<-min(startY+nts-1,dim(thisExchangeData)[4])
      #exchange data
      thisValsData<-thisExchangeData[,,,startY:endY]
      ###################
      ##SORT OUT DIMENSIONS FOR WHEN F = NFILES
      ###################
      if(startY+nts-1>dim(thisExchangeData)[4]){
        valsData<-0*thisExchangeData[,,,1:nts]
        valsData[,,,1:dim(thisValsData)[4]]<-thisValsData
      }else{
        valsData<-thisValsData
      }
      ncvar_put(newHydroNC,varid = "exchange",vals = valsData)
      # destb data
      thisValsData<-thisDestB[,,,startY:endY]
      if(startY+nts-1>dim(thisDestB)[4]){
        valsData<-0*thisDestB[,,,1:nts]
        valsData[,,,1:dim(thisValsData)[4]]<-thisValsData
      }else{
        valsData<-thisValsData
      }
      ncvar_put(newHydroNC,varid = "dest_b",vals = valsData)
      # destk data
      thisValsData<-thisDestK[,,,startY:endY]
      if(startY+nts-1>dim(thisDestK)[4]){
        valsData<-0*thisDestK[,,,1:nts]
        valsData[,,,1:dim(thisValsData)[4]]<-thisValsData
      }else{
        valsData<-thisValsData
      }
      ncvar_put(newHydroNC,varid = "dest_k",vals = valsData)
      #close the nc file
      nc_close(newHydroNC)
    }
  }
}

run_nruns<-50 # I ended up cutting this down to 50 runs as the nc files were taking up too much space! Can do 2 lots of 50 if needed but 50 may be enough for now

################################
#create run.prm file. start outputting from after burnin
tstart<-burnin*365; tstop<-(burnin+runLength)*365; tinc<-365 #output every year
newRunPrmFilename<-paste("CRAM_base_runROMS",version,".prm",sep="")
existingRunPrm<-paste(basePath,"CRAM_base_run_short.prm",sep=""); newRunPrm<-paste(basePath,newRunPrmFilename,sep="")
file.copy(from=existingRunPrm, to= newRunPrm, overwrite=TRUE)
runPrmLines<-readLines(newRunPrm)
##tstop
x<-grep("^tstop",runPrmLines)
newLine<-paste("tstop",tstop,"day",collapse="\t")
runPrmLines[x]<-newLine
## toutinc
x<-grep("^toutinc",runPrmLines)
newLine<-paste("toutinc",tinc,"day",collapse="\t")
runPrmLines[x]<-newLine
## toutstart
x<-grep("^toutstart",runPrmLines)
newLine<-paste("toutstart",tstart,"day",collapse="\t")
runPrmLines[x]<-newLine
##write the lines back out
writeLines(runPrmLines,newRunPrm)

###########################
##create force file for each ROMS bootstrap run
existingForcePrm<-paste(basePath,"inputs\\CRAM_force.prm",sep=""); 
for(r in 1:run_nruns){
  newForccePrm<-paste(basePath,"inputs\\ROMS\\BootstrapROMS\\",version,"\\CRAM_force_ROMS",version,r,".prm",sep="")
  file.copy(from=existingForcePrm,to=newForccePrm,overwrite = TRUE)
  forceLines<-readLines(newForccePrm)
  #hydro files
  x<-grep("^hd",forceLines)
  newHdLines<-paste("hd",seq(0,nfiles-1),".name inputs/ROMS/BootstrapROMS/",version,"/Chatham_hydro",r,"_",seq(1,nfiles),".nc\n",sep="")
  newForceLines<-c(forceLines[1:(x-1)],newHdLines,forceLines[(x+1):length(forceLines)])
  forceLines<-newForceLines
  
  x<-grep("nhdfiles",forceLines)
  forceLines[x]<-paste("nhdfiles",nfiles,collapse="\t")
  
  #temperature files
  x<-grep("^Temperature",forceLines)
  newHdLines<-paste("Temperature",seq(0,nfiles-1),".name inputs/ROMS/BootstrapROMS/",version,"/Chatham_temp",r,"_",seq(1,nfiles),".nc\n",sep="")
  newForceLines<-c(forceLines[1:(x-1)],newHdLines,forceLines[(x+1):length(forceLines)])
  forceLines<-newForceLines
  
  x<-grep("ntempfiles",forceLines)
  forceLines[x]<-paste("ntempfiles",nfiles,collapse="\t")
 
  #salinity files
  x<-grep("^Salinity",forceLines)
  newHdLines<-paste("Salinity",seq(0,nfiles-1),".name inputs/ROMS/BootstrapROMS/",version,"/Chatham_salt",r,"_",seq(1,nfiles),".nc\n",sep="")
  newForceLines<-c(forceLines[1:(x-1)],newHdLines,forceLines[(x+1):length(forceLines)])
  forceLines<-newForceLines
  
  x<-grep("nsaltfiles",forceLines)
  forceLines[x]<-paste("nsaltfiles",nfiles,collapse="\t")
  
  write(forceLines,newForccePrm)
}

#create RUN file for Turbine
thisICncFile<-"CRAM_input_short.nc" ##I need to update this with outs from a stable fun, just 5 years in for now
runFile<-paste(writePath,"RunROMSBootstrap_",version,sep="")
runText<-paste("#Run multiple runs with varying ROMS files, version ",version,", \n# ",sep="")
cat(runText,file=runFile,append=FALSE)

for(r in 1:nruns){
  runText<-paste("\nWD=\"$(pwd)\"\nRUN=\"../../bin/bin/atlantisMerged -i ",thisICncFile," 0 -o output.nc -r ",newRunPrmFilename," -f inputs/CRAM_force_ROMS",version,r,".prm -p inputs/CRAM_physics.prm -b CRAM_base_biol.prm -s CRAM_Groups.csv -q inputs/CRAM_Fisheries.csv -d base/outputROMSBootstrap",version,r,"\"\necho $RUN > RUN\nCMD=\"msub -l nodes=1 -l walltime=50:00:00 -l partition=slurm -l qos=standby -p -1000 -q large -o CRAMROMSBoostrap",version,i,".log.%j -e CRAMROMSbootstrap",version,r,".err.%j -S /bin/bash RUN\"\necho \"Running Atlantis ROMS bootstrap ",version,r," for CRAM on MOAB in directory:\" $WD\necho -n \"Job started at: \" ; date\necho $RUN\nCOMMAND=\"cd $WD ; $CMD\"\nssh turbine $COMMAND\nsleep 0.5",sep="")
  cat(runText,file=runFile,append=TRUE)
  
}

#create run instructions for turbine
filename<-paste(writePath,"CopyIntoTurbine.txt",sep="")
text<-paste("cp RunROMSBootstrap_",version, " RunROMSBootstrap_",version,".backup
dos2unix RunROMSBootstrap_",version,"
sed -i -r \"s/sleep 0.5/sleep 1/\" RunROMSBootstrap_",version,"
./RunROMSBootstrap_",version,sep="")
writeLines(text,filename)


#may need to run dos2unix on the force.prm files
filename<-paste(writePath,"CopyrIntoTurbine2fixForceFiles.txt",sep="")
writeLines("#\n",filename)
for(r in 1:nruns){
  thisForceFile<-paste("CRAM_force_ROMS",version,r,".prm",sep="")
  thistext<-paste("dos2unix ", thisForceFile,"\n",sep="")
  cat(thistext,file=filename,append=TRUE)
}

