##  read in and set up the full ROMS files

basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
version<-"A"

nlayer<-6; nboxes<-30
burnin<-35; runLength<-20
nruns<-20

#read in the full ROMS files
ROMSpath<-paste(basePath,"inputs\\ROMS\\",sep="")

inputROMSPath<-paste(ROMSpath,"BootstrapROMS\\",version,"\\",sep="")

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

##read in bootstrap year samples
startRun=1
bootstrapFile<-paste(inputROMSPath,"RecordedYearsSelected.txt",sep=""); bs_lines<-readLines(bootstrapFile)
storeRandomYears<-array(NA, dim=c(nruns,runLength))
indexRuns<-seq(startRun,(startRun+nruns-1))
for(r in 1:nruns){
  this_r<-indexRuns[r]
  x<-bs_lines==paste("Run ",this_r,sep=""); xi<-seq(1,length(bs_lines))[x]
  thisLine<-bs_lines[xi+1]; this_bs<-get_first_number(thisLine,n="all")
  storeRandomYears[r,]<-this_bs
}

baseYear<-2003; baseYearIndex<-grep(baseYear,seq(1996,2004))

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
  testFile<-paste(inputROMSPath,"Chatham_temp",r,"_1.nc",sep="")
  # if(file.exists(testFile)){
  #   cat("ROMS files seem to exist - delete first to ensure we don't write over something that should be kept\n")
  # }else{
    thisYearSelection<-c(rep(baseYearIndex,burnin),storeRandomYears[r,])
    thisSaltData<-emptyArray; thisTempData<-emptyArray; thisExchangeData<-emptyArrayHydro; thisDestBData<-emptyArrayHydro; thisDestKData<-emptyArrayHydro
    thisDestB<-emptyArrayHydro; thisDestK<-emptyArrayHydro
    for(y in 1:(runLength+burnin)){
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
      newTempNCfile<-paste(inputROMSPath,"Chatham_temp",r,"_",f,".nc",sep="")
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
      newSaltNCfile<-paste(inputROMSPath,"Chatham_salt",r,"_",f,".nc",sep="")
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
      newHydroNCfile<-paste(inputROMSPath,"Chatham_hydro",r,"_",f,".nc",sep="")
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
  # }
}


#create RUN file for Turbine
newRunPrmFilename<-"CRAM_base_runROMSA.prm"
thisICncFile<-"CRAM_input_short.nc" ##I need to update this with outs from a stable fun, just 5 years in for now
runFile<-paste(inputROMSPath,"RunROMSBootstrap_",version,"_",startRun,sep="")
runText<-paste("#Run multiple runs with varying ROMS files, version ",version,", \n# ",sep="")
cat(runText,file=runFile,append=FALSE)

for(r in 1:nruns){
  runText<-paste("\nWD=\"$(pwd)\"\nRUN=\"../../bin/bin/atlantisMerged -i ",thisICncFile," 0 -o output.nc -r ",newRunPrmFilename,
                 " -f inputs/CRAM_force_ROMS",version,r,
                 ".prm -p inputs/CRAM_physics.prm -b CRAM_base_biol.prm -s CRAM_Groups.csv -q inputs/CRAM_Fisheries.csv -d base/outputROMSBootstrap",
                 version,r,"\"\necho $RUN > RUN\nCMD=\"msub -l nodes=1 -l walltime=50:00:00 -l partition=slurm -l qos=standby -p -1000 -q large -o CRAMROMSBoostrap",
                 version,r,".log.%j -e CRAMROMSbootstrap",version,r,".err.%j -S /bin/bash RUN\"\necho \"Running Atlantis ROMS bootstrap ",version,r,
                 " for CRAM on MOAB in directory:\" $WD\necho -n \"Job started at: \" ; date\necho $RUN\nCOMMAND=\"cd $WD ; $CMD\"\nssh turbine $COMMAND\nsleep 0.5",sep="")
  cat(runText,file=runFile,append=TRUE)
  
}
# 
# ##read in a couple of temp ncs to check they are different
# f=6; r=48
# test1<-nc_open(paste(inputROMSPath,"test\\Chatham_temp",r,"_",f,".nc",sep="")); temp1<-ncvar_get(test1,"temperature")
# f=6; r=50
# test2<-nc_open(paste(inputROMSPath,"test\\Chatham_temp",r,"_",f,".nc",sep="")); temp2<-ncvar_get(test2,"temperature")
# 
# xx<-apply(temp1,3,mean,na.rm=TRUE)
# yy<-apply(temp2,3,mean,na.rm=TRUE)
# 

