#if we want to perturb one variable at a time, or subsets at a time, how might we want to structure this?
## what is the spread of variables by lifespan or trophic level? Depth range?

source(paste(DIR$`General functions`,"getVolDepth.R",sep=""))
source(paste(DIR$`General functions`,"read_boxes.R",sep=""))

getTracerFromNums<-function(x){
  thisCohort <- get_first_number(x)
  yy<-gsub(paste(thisCohort,"_Nums", sep=""), "", x)
  return(yy)
}

version<- ""
basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

existingICfile<-"CRAM_input_fromBase50yr"; 

init <- nc_open(paste(basePath,existingICfile,".nc",sep=""), write = F) 
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

vars<-unique(sort(var_names_init))
# skipVars<-c("nominal_dz", "volume", "dz", "numlayers", "topk")
skipVars <- c("canyon", "Chl_a", "CO2", "Denitrifiction",  "DiagNFlux", "DiagNGain", "DiagNLoss", "DON", "eddy", "eflux", "erosion_rate", "Filter_Other_Cover", 
              "flat", "hdsink", "hdsource", "Light", "Light_Adaptn_DF", "Light_Adaptn_MB", "Light_Adaptn_PL", "Light_Adaptn_PS", "Macroalgae_Cover", "MicroNut", 
              "MicroPB_Cover",   "Nitrification",  "porosity", "reef", "salt", "sedbiodens", "sedbiodepth", "seddetdepth", "sedirrigenh", "sedoxdepth", "sedturbenh", 
              "soft", "Stress", "Temp", "vflux", "water", "nominal_dz", "volume", "dz", "numlayers", "topk", "Oxygen", "NH3", "NO3", "Si", "Det_Si")
vars<-vars[!(vars %in% skipVars)]

## take out ResN and StructN
weightVars <- vars[grep("StructN|ResN", vars)]
vars <- vars[!(vars %in% weightVars)]
## 398 in total here

## take out any _N vars that have a matching 'Nums' var - in these the _N initial condition won't be used anyway
xx<-vars[grep("Nums", vars)]
numVars <- unlist(lapply(vars[grep("Nums", vars)], getTracerFromNums)); numNvars <- paste(unique(numVars),"_N", sep="")

## take the numNvars out of vars
vars<-vars[!(vars %in% numNvars)]
## 361 now

## how many of those left have _N (this will include _Nums)
allNvars <- vars[grep("_N", vars)] # 359 of these

test<-vars[!(vars %in% allNvars)]

## write vars out so can use them elsewhere
write.csv(vars, paste(DIR$'Table',"VARSforInitialConditions.csv", sep=""), row.names = FALSE)

# do 100 runs, and for each, vary each variable randomly
# read in rankings so can influence variance or initial conditions based on how well informed the group is
rankTable<-read.csv(paste(DIR$'Table',"interactionEffectsRANKINGs.csv", sep=""))
nageSdev<-0.05
## there are 4 informance rankings, so assign 4 standard deviations
ageDevs <- c(0.01,0.02,0.05,0.1)
resultingScalars <- rev(c(0.05, 0.1, 0.2, 0.5))
ageDevs <- (resultingScalars/1.96)
par(mfrow=c(2,2))
for(a in 1:length(ageDevs)){
  thisX=seq(-1,1,length.out=100); thisY<-dnorm(thisX, sd=ageDevs[a])
  plot(x=thisX, y=thisY, type="l")
  abline(v=c(-1,1)*resultingScalars[a], col=myOrange, lwd=2)
}

nvars<-length(vars)

infGrades <-sort(unique(rankTable$informPerformRank))

# sample the scalars for each species variable for each of (just 10 for now) runs
nruns <- 10; sampleVersion<-"KEY"
runSetupFile <- paste(basePath, "Chaos\\Setup\\", sep="")
thisSampleFile <- paste(runSetupFile,"chaosSampleSetup",sampleVersion,"n",nruns,".csv", sep="")

storeSampledScalars <- array(NA, dim=c(nruns, nvars))
thisSamples <- rep(NA, nvars)
doneVars<-c()
for(i in 1:nvars){
  thisVar<-vars[i]
  if(!(thisVar %in% doneVars)){
    # if 10th or higher for keystoneness, give SD of 0.25 (95% approx +/- 0.5); 0 otherwise
    asTest<-grep("_Nums", thisVar)
    if(length(asTest)==0){
      thisSD <- nageSdev
      allMyVarsIndex <- vars == 0
      
    }else{
      thisAge<-get_first_number(thisVar)
      thisName<-gsub(paste(thisAge,"_Nums", sep=""),"",thisVar)
      thisCode<-as.character(groupsDF$Code[str_trim(groupsDF$Name, side="both")==thisName])
      thisKey <- rankTable$KeyRank[rankTable$Code==thisCode]
      thisSD<-0
      if(thisKey<=10){ thisSD<-0.25}
      ## need to apply this to all the ageclasses of this group - so M is preserved
      thisNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisCode]
      allMyVars <- paste(thisName, seq(1, thisNumCohorts), "_Nums", sep="")
      allMyVarsIndex<-vars %in% allMyVars
      doneVars <- c(doneVars,allMyVars)
    }
    # sample from N(0,thisSD)
    thisScalar <- rnorm(nruns, 0, thisSD)
    storeSampledScalars[,allMyVarsIndex]<-thisScalar
  }
}
###### write the sampled scalars out
# write.csv(storeSampledScalars, file=thisSampleFile, row.names = FALSE)

## create initial conditions for each run
existingICfile<-"CRAM_input_fromBase50yr"; 
for(r in 1:nruns){
  newICFile<-paste(existingICfile,"sample",sampleVersion, r,sep="")
  
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
  
  ic_vars<-unique(sort(var_names_init))
  # skipVars<-c("nominal_dz", "volume", "dz", "numlayers", "topk")
  skipVars <- ic_vars[!(ic_vars %in% vars)]
  
  for(i in seq_along(vars)){ 
    thisVar <- vars[i]
    ## if thisVar has ResN or StructN, skip as don't want to change weights
    thisScale <- 1 + storeSampledScalars[r,i]  
    if(is.na(thisScale)){thisScale <- 1}
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



## create a .bat file to dump all inton text files so can fix the dimensions - run this file from command prompt once have created it
batFile <- paste(basePath,sampleVersion, "dumpChaosInput2text.bat", sep="")
cat("##\n", file=batFile, append=FALSE)
for(r in 1:nruns){
  thisNCfile<-paste(existingICfile,"sample",sampleVersion, r,sep="")
  
  thisLine <- paste("ncdump ", thisNCfile, ".nc > ", thisNCfile, ".txt \n", sep="")
  cat(thisLine, file=batFile, append=TRUE)
}

## fix the dimensions
for(r in 1:nruns){
  thisfile<-paste(basePath,existingICfile,"sample",sampleVersion, r,".txt",sep="")
  
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


## file to turn back into .nc files
## create a .bat file to dump all inton text files so can fix the dimensions - run this file from command prompt once have created it
batFile <- paste(basePath, sampleVersion, "dumpChaosText2input.bat", sep="")
cat("##\n", file=batFile, append=FALSE)
for(r in 1:nruns){
  thisNCfile<-paste(existingICfile,"sample",sampleVersion, r,sep="")
  
  thisLine <- paste("ncgen -o ", thisNCfile, ".nc  ", thisNCfile, ".txt \n", sep="")
  cat(thisLine, file=batFile, append=TRUE)
}



## set up the run file
baseICfile<- "CRAM_input"
baseRunCommand <- "../../bin/bin/atlantisMerged -i CRAM_input.nc 0 -o output.nc -r CRAM_base_run.prm -f inputs/CRAM_forceBURNIN1865.prm -p inputs/CRAM_physics.prm -b CRAM_BH_hybrid_biol.prm -s CRAM_Groups.csv -q CRAM_Fisheries.csv -d outputFolder"
fishRunCommand <- "../../bin/bin/atlantisMerged -i CRAM_input.nc 0 -o output.nc -r CRAM_baseFish_run.prm -f inputs/CRAM_forceBURNIN1865.prm -p inputs/CRAM_physics.prm -b CRAM_BH_hybrid_biol.prm  -h CRAM_harvest_short.prm -s CRAM_Groups.csv -q CRAM_Fisheries.csv -d outputFolder"
runFile<-paste(basePath, "RunChaosSample",sampleVersion, sep="")
cat("#Run base run and historical catches removed run with scaled ICs", file=runFile, append=FALSE)
for(r in 1:nruns){
  thisNCfile<-paste(existingICfile,"sample",sampleVersion, r,sep="")
  thisRunCommand <- gsub(baseICfile, thisNCfile, baseRunCommand)
  thisRunCommand <- gsub("outputFolder", paste("outputChaosSample",sampleVersion,r, sep="" ), thisRunCommand)
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
  thisRunCommand <- gsub("outputFolder", paste("outputChaosFISHSample",sampleVersion,r, sep="" ), thisRunCommand)
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

###################
## plot the resulting spreads, colored by keystoneness
colRamp<-c(myGrey, myOrange)

getColor<-function(x){
  thisCol<-myGrey
  if(!is.na(x)>0){
    if(x<=10){thisCol<-myOrange}
  }
  return(thisCol)
}
ksByVar <- rep(NA, nvars)
for(i in 1:nvars){
  thisVar <- vars[i]
  thisCohort <- get_first_number(thisVar); thisVar <-gsub(thisCohort, "", thisVar)
  thisName <- gsub("_N|_Nums", "", thisVar); thisName <- gsub("_", " ", thisName)
  thisCode <- rankTable$Code[grep(thisName, rankTable$Group)]
  thisKS<-rankTable$KeyRank[grep(thisName, rankTable$Group)]
  ksByVar[i]<-thisKS
}
colByPer<-unlist(lapply(ksByVar, getColor))
colByPer[is.na(colByPer)]<-myGrey

addedSampledScalars <- 1+storeSampledScalars

thisMax <- max(addedSampledScalars, na.rm=TRUE); thisMin <- min(addedSampledScalars, na.rm=TRUE)
pdf(paste(plotPath, ""))
plot(addedSampledScalars[1,], type="n", ylim=c(thisMin, thisMax), xlab="Variable", ylab="Scalar")
for(r in 1:nruns){
  points(addedSampledScalars[r,], pch=20, col=colByPer)
}

ksIndex<-ksByVar<=10 & !is.na(ksByVar)
storeHistBreaks <- NULL; storeHist <- NULL
ymax <- 0; xmin<-10; xmax=0
for(r in 1:nruns){
  thisHist <- hist(addedSampledScalars[r,ksIndex], plot=FALSE)
  storeHistBreaks[[r]]<-thisHist$breaks; storeHist[[r]]<-thisHist$counts
  if(max(thisHist$counts)>ymax){ymax<-max(thisHist$counts)}
  if(min(thisHist$breaks)<xmin){xmin <- min(thisHist$breaks)}
  if(max(thisHist$breaks > xmax)){xmax <- max(thisHist$breaks)}
}
xshift<-seq(-0.2,0.2, length.out=nruns)
colByRun <- colorRampPalette(colors=c(myGold,myOrange, "red", myPurple, myAqua))(nruns)
par(lend=1)
plot(1, type="n", ylim=c(0, ymax), xlim=c(xmin, xmax), ylab="Frequency", xlab="Scalar")
for(r in 1:nruns){
  points(x=storeHistBreaks[[r]][-1]+xshift[r], y=storeHist[[r]], type="h", lwd=5, col=colByRun[r])
}

## heat map it

## now to plot as heat map
getColor<-function(x,thisMax){
  thisCol<-"white"
  if(!is.na(x) & x>0){
    y<-round((x-thisMin)/(thisMax-thisMin),2)*100+1
    thisCol<-thisColRamp[y]
  }
  return(thisCol)
}
plotGrid<-function(x,y){
  thisX<-c(x-0.5,x-0.5,x+0.5,x+0.5); thisY<-c(y-0.5,y+0.5,y+0.5,y-0.5)
  thisCol<-plotColour[x,y]
  if(length(thisCol)>0){
    polygon(x=thisX,y=thisY,col=thisCol,border=NA)
  }
  return(NULL)
}
thisColRamp<-colorRampPalette(colors=c(myLightAqua,myAqua,"midnightblue"))(101)

plotData<-t(addedSampledScalars[,ksIndex])
thisMax<-max(plotData, na.rm=TRUE); thisMin<-min(plotData, na.rm=TRUE)

plotColour<-apply(plotData,c(1,2),getColor,thisMax)
tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotData)[1]),dim(plotData)[2]),"y"=sort(rep(seq(1,dim(plotData)[2]),dim(plotData)[1]))))
# pdf(paste(plotPath,"KemptonsQ_SSRsens.pdf",sep=""),height=4,width=10)
par(mar=c(6,4,1.5,1))
plot(x=seq(1,dim(plotData)[1]),y=rep(dim(plotData)[2],dim(plotData)[1]),type="n",xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,dim(plotData)[2]))
axis(at=seq(1,dim(plotData)[1]),labels = seq(1,dim(plotData)[1]),side=1,las=2)
axis(at=seq(1,dim(plotData)[2]),labels=1:nruns,side=2,las=1)
temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
box()
# dev

hist(addedSampledScalars[,ksIndex])
plot(as.double(storeSampledScalars))

