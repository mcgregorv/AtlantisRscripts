##  read in and set up the full ROMS files
#bootstrap the ROMS years - and save them in an output file. Then create force.prm files for each run that determines the ROMS years to use and order
# for each model run, repeat 2003 ROMS 35 times for the burnin, then bootstrap ROMS years for a further 50 years. The run.prm file will be set up to start outputting from 35 years

basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
version<-"B"

nlayer<-6; nboxes<-30
burnin<-35; runLength<-50
nruns<-50; nROMSfiles<-burnin+runLength
#read in the full ROMS files
ROMSpath<-paste(basePath,"inputs\\ROMS\\",sep="")

inputROMSPath<-paste(ROMSpath,"BootstrapROMS\\",version,"\\",sep="")

ROMSyears<-seq(1996,2004); baseYear<-2003; ROMSyearIndex<-seq(1,length(ROMSyears))

baseForceFile<-paste(basePath,"inputs\\CRAM_forceBURNIN1865.prm",sep="")
forceLines<-readLines(baseForceFile)

storeSamplesROMSyears<-array(NA,dim=c(nruns,nROMSfiles))
for(r in 1:nruns){
  storeSamplesROMSyears[r,1:burnin]<-rep(baseYear,burnin)
  storeSamplesROMSyears[r,(burnin+1):nROMSfiles]<-sample(ROMSyears,size=runLength,replace=TRUE)
}
#write them out to store. then comment out so don't overwrite
# write.csv(storeSamplesROMSyears,paste(inputROMSPath,"ROMSsamples.csv",sep=""))

for(r in 1:nruns){
  thisForceFile<-paste(inputROMSPath,"CRAM_force_ROMS",version,r,".prm",sep="")
  thisForceLines<-forceLines
  thisROMSyears<-storeSamplesROMSyears[r,]
  thisROMSindex<-ROMSyearIndex[match(thisROMSyears,ROMSyears)]
  #edit nfiles for salinity, temp, and exchange
  thisVar<-"nhdfiles"; x<-grep(thisVar,thisForceLines)
  newLines<-paste(thisVar,nROMSfiles,sep=" "); thisForceLines[x]<-newLines
  ##
  thisVar<-"nsaltfiles"; x<-grep(thisVar,thisForceLines)
  newLines<-paste(thisVar,nROMSfiles,sep=" "); thisForceLines[x]<-newLines
  ##
  thisVar<-"ntempfiles"; x<-grep(thisVar,thisForceLines)
  newLines<-paste(thisVar,nROMSfiles,sep=" "); thisForceLines[x]<-newLines
  ## now add in the lines for each year to tell it which ROMS year to use
  ##hydro
  thisROMSfilenames<-paste(".name inputs/ROMS/Chatham_hydro_year",thisROMSyears, ".nc",sep="")
  thisROMSlines<-paste("hd",seq(0,(nROMSfiles-1)),thisROMSfilenames,sep="",collapse = "\n")
  x<-grep("hd0.name inputs",thisForceLines)
  y<-thisForceLines[x]; z<-x[grep("#",y,invert = TRUE)]
  thisForceLines[z]<-thisROMSlines
  ##temp
  thisROMSfilenames<-paste(".name inputs/ROMS/Chatham_temp_year",thisROMSyears, ".nc",sep="")
  thisROMSlines<-paste("Temperature",seq(0,(nROMSfiles-1)),thisROMSfilenames,sep="",collapse = "\n")
  x<-grep("Temperature0.name inputs",thisForceLines)
  y<-thisForceLines[x]; z<-x[grep("#",y,invert = TRUE)]
  thisForceLines[z]<-thisROMSlines
  ##salt
  thisROMSfilenames<-paste(".name inputs/ROMS/Chatham_salt_year",thisROMSyears, ".nc",sep="")
  thisROMSlines<-paste("Salinity",seq(0,(nROMSfiles-1)),thisROMSfilenames,sep="",collapse = "\n")
  x<-grep("Salinity0.name inputs",thisForceLines)
  y<-thisForceLines[x]; z<-x[grep("#",y,invert = TRUE)]
  thisForceLines[z]<-thisROMSlines
  
  
  writeLines(thisForceLines,thisForceFile)
}


#now create the run file

#create RUN file for Turbine
newRunPrmFilename<-paste("CRAM_base_runROMS",version,".prm",sep="")
thisICncFile<-"CRAM_input_short.nc"
runFile<-paste(inputROMSPath,"RunROMSBootstrap_",version,sep="")
runText<-paste("#Run multiple runs with varying ROMS files, version ",version,", \n# ",sep="")
cat(runText,file=runFile,append=FALSE)

for(r in 1:nruns){
  thisForceFile<-paste("CRAM_force_ROMS",version,r,".prm",sep="")
  runText<-paste("\nWD=\"$(pwd)\"\nRUN=\"../../bin/bin/atlantisMerged -i ",thisICncFile," 0 -o output.nc -r ",newRunPrmFilename,
                 " -f inputs/",thisForceFile,
                 " -p inputs/CRAM_physics.prm -b CRAM_BH_hybrid_biol.prm -s CRAM_Groups.csv -q inputs/CRAM_Fisheries.csv -d base/outputROMSBootstrap",
                 version,r,"\"\necho $RUN > RUN\nCMD=\"msub -l nodes=1 -l walltime=50:00:00 -l partition=slurm -l qos=standby -p -1000 -q large -o CRAMROMSBoostrap",
                 version,r,".log.%j -e CRAMROMSbootstrap",version,r,".err.%j -S /bin/bash RUN\"\necho \"Running Atlantis ROMS bootstrap ",version,r,
                 " for CRAM on MOAB in directory:\" $WD\necho -n \"Job started at: \" ; date\necho $RUN\nCOMMAND=\"cd $WD ; $CMD\"\nssh turbine $COMMAND\nsleep 0.5",sep="")
  cat(runText,file=runFile,append=TRUE)
  
}










