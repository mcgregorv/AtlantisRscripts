## catch histories were created in setUpScenario_ts_files_versionC1.R - they are full historic + 50 year future runs
## biol.prm files for eddy_scale are  created first

thisPath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

version<-"A"
version<-"B"
version<-"C"  # just like B, but has eddy strength change after the historic period, and has alt fishing scenarios

initFile<-"CRAM_input_short.nc";   thisFishPrm<-paste("CRAM_baseFish_run_200years.prm")

thisRunFile<-paste(thisPath,"RUN_EDDY_histAndFuture_scen",version,sep="")

## set up DF with runs specified
scenarios<-c("All0catch", "All50catch","All100catch", "All150catch", "Hoki50catch", "Hoki150catch"); 
scenarios<-c("All0catch") ## B
scenarios<-c("All0catch", "All50catch","All100catch", "All150catch") ## C
scenarios<-c("All0catch", "All50catch", "All80catch",  "All100catch",  "All120catch",  "All150catch") ## C


nscenarios<-length(scenarios)
eddysens<-c(0.01,0.05,0.1,0.5,1,5,10,50,100); 
eddysens<-c(0.025,0.05,0.075,0.1,0.25,0.5,0.75,1)
neddys<-length(eddysens); eddyIndex<-seq(1,neddys)
nruns<-nscenarios * neddys

lookup_df<-data.frame(array(NA, dim=c(nruns,3)))
colnames(lookup_df)<-c("Run", "Scenario", "Eddy")
lookup_df$Run<-seq(1,nruns)
lookup_df$Scenario<-rep( sort(rep(scenarios, neddys)))
lookup_df$Eddy<-rep(eddyIndex, (nscenarios))

baseBiol<-paste(thisPath,"CRAM_BH_hybrid_biol.prm", sep=""); biolLines<-readLines(baseBiol)

## create the biol.prm files
for(e in 1:neddys){
  thisEddy<-eddysens[e]
  thisFile<-paste(thisPath,"CRAM_BH_hybrid_biolEDDY",e,".prm", sep="")
  x<-grep("eddy_scale",biolLines)
  newLine<-paste("eddy_scale\t",thisEddy)
  newLines<-biolLines
  newLines[x]<-newLine
  writeLines(newLines, thisFile)
}

after<-"\necho $RUN > RUN
CMD=\"msub -l nodes=1 -l walltime=100:00:00 -l partition=slurm -l qos=standby -p -1000 -q large -o CRAMProjection.log.%j -e CRAMProjection.err.%j -S /bin/bash RUN\"
echo \"Running Atlantis Projection for CRAM on MOAB in directory:\" $WD
echo -n \"Job started at: \" ; date
echo $RUN
COMMAND=\"cd $WD ; $CMD\"
ssh turbine $COMMAND
sleep 0.5
WD=\"$(pwd)\"\n"
before<-"\nWD=\"$(pwd)\"\n"

cat("## hist and future eddy scalar runs\n", file=thisRunFile, append=FALSE)

for(r in 1:nruns){
  thisBiolPrm<-paste("CRAM_BH_hybrid_biolEDDY",lookup_df$Eddy[r], ".prm", sep="")
  thisForcePrm<-paste("CRAM_force_FS_C1",lookup_df$Scenario[r], ".prm", sep="")
  thisOUT<-paste("outputEDDY",version,"_",r, sep="")
  thisRunCommand<-paste("RUN=\"../../bin/bin/atlantisMerged -i ", initFile, " 0 -o output.nc -r ", thisFishPrm,
                        " -f inputs/", thisForcePrm, " -p inputs/CRAM_physics.prm -b ", thisBiolPrm, "  -h CRAM_harvest_200years.prm  -s CRAM_Groups.csv -q CRAM_Fisheries.csv -d EDDYsens/", thisOUT, "\"", sep="")
  
  cat(before, file=thisRunFile, append=TRUE)
  cat(thisRunCommand, file=thisRunFile, append=TRUE)
  cat(after, file=thisRunFile, append=TRUE)
}

