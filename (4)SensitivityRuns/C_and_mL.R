#what if we vary X_RS? Needs to change in biol.prm file (X_RS, KWSR and KWRR), and change the ratio RN/SN for all 'verts' in the initial conditions
#especially interested in this wrt recruitment
##so permuting with changes in C (attack rate) as well

basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
baseBiolPrm<-paste(basePath,"CRAM_base_biol_BH_mixedUnits075.prm",sep=""), steepnessText<-""
baseBiolPrm<-paste(basePath,"CRAM_base_biol_BH_mixedUnits.prm",sep=""); steepnessText<-"025"



biolLines<-readLines(baseBiolPrm)

baseICfile<-paste(basePath,"CRAM_input_short.txt",sep="")
IClines<-readLines(baseICfile)

groupsDF<-read.csv(paste(basePath,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

nboxes<-30; nlayers<-6

#mL will be multiplied by existing mL, result will need to be capped at something (set to 0.7 here)
max_mL<-0.7
test_mL<-c(0.01,0.1,0.5,1,1.5,5,10,100); nMruns<-length(test_mL)

testC<-c(1,10,100); nCruns<-length(testC)


for(m in 1:nMruns){
  this_mL_scalar<-test_mL[m]
  for(a in 1:nCruns){
    thisC<-testC[a]
    newBiolFile<-paste(basePath,"CRAM_biol_mL",m,"C",thisC,steepnessText,".prm",sep="")
    
    newBiolLines<-biolLines; 
    
    #scale mL in biolLines by this_mL_scalar
    x<-grep("_mL",newBiolLines); nx<-length(x)
    for(i in 1:nx){
      thisLine<-newBiolLines[x[i]+1]
      temp<-get_first_number(thisLine,"all")
      if(length(temp)==2){ #check if age structured
        temp2<-temp*this_mL_scalar; temp2[temp2>max_mL]<-max_mL
        newLine<-paste(temp2,collapse="\t")
        newBiolLines[x[i]+1]<-newLine
      }
    }
    #update C for all groups
    for(g in 1:ng){
      thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
      if(thisNumCohorts>1){
        thisVar<-paste("C_",thisCode,sep=""); x<-grep(thisVar,newBiolLines)
        newCline<-paste(rep(thisC,thisNumCohorts),collapse="\t")
        newBiolLines[x+1]<-newCline
      }
    }
    #write the new biolLines out
    writeLines(newBiolLines,newBiolFile)
  }
}

#create the run file


#creates the RUN file. biolprm files have already been created
runDF<-data.frame(matrix(NA,ncol=4,nrow=nCruns*nMruns)); colnames(runDF)<-c("Run","biolprm","ICfile","output")

for(r in 1:nMruns){
  for(a in 1:nCruns){
    thisR<-(r-1)*nCruns+a
    runDF$biolprm[thisR]<-paste("CRAM_biol_mL",r,"C",testC[a],steepnessText,".prm",sep="")
    runDF$ICfile[thisR]<-paste("CRAM_input_X_RSbase.nc",sep="")
    runDF$output[thisR]<-paste("Output_mL",r,"C",a,steepnessText,sep="")
    runDF$Run[thisR]<-thisR
  }
}

thisRunFile<-paste(basePath,"RunComparitive_mL_andC",steepnessText,".bat",sep="")
cat("#Runing 15 runs with varying 6 mL and 3 C\n",file=thisRunFile,append=FALSE)

for(r in 1:nMruns){
  for(a in 1:nCruns){
    thisR<-(r-1)*nCruns+a
    thisRunText<-paste("WD=\"$(pwd)\"
                       RUN=\"../../bin/bin/atlantisMerged -i ",runDF$ICfile[thisR]," 0 -o output.nc -r CRAM_base_run.prm -f inputs/CRAM_force.prm -p inputs/CRAM_physics.prm -b  ", runDF$biolprm[thisR]," -s CRAM_Groups.csv -q inputs/CRAM_Fisheries.csv -d ", runDF$output[thisR],"\"
                       echo $RUN > RUN
                       CMD=\"msub -l nodes=1 -l walltime=50:00:00 -l partition=slurm -l qos=standby -p -1000 -q large -o CRAM_X_RS",r,".log.%j -e CRAM_X_RS",r,".err.%j -S /bin/bash RUN\"
                       echo \"Running Atlantis Comparitive X_RS for CRAM on MOAB in directory:\" $WD
                       echo -n \"Job started at: \" ; date
                       echo $RUN
                       COMMAND=\"cd $WD ; $CMD\"
                       ssh turbine $COMMAND
                       sleep 1\n",sep="")
    cat(thisRunText,file=thisRunFile,append=TRUE)
  }
}
