#read in the tracers (numbers), calculate spawning stock abundance, calculate rectruited abundance, then calculate steepness
this_run<-"base"
this_out<-paste("TEST2",sep="")

basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
# outPath<-paste(basePath,this_run,"\\","output",this_out,"_h",thisSteepness,"\\",sep="")
outPath<-paste(basePath,this_run,"\\","output",this_out,"\\",sep="")

groupsDF<-read.csv(paste(basePath,"CRAM_groups.csv", sep=""))
groupsAgeDF<-groupsDF[groupsDF$NumCohorts>1,]; nag<-dim(groupsAgeDF)[1]

biolLines<-readLines(paste(basePath,"CRAM_BH_hybrid_biol.prm",sep=""))


##read in biol.prm file and check out which recruit over the change of year
Time_Spawns<-rep(NA,nag); spawn_period<-rep(NA, nag)
for(g in 1:nag){
  thisCode<-groupsAgeDF$Code[g]
  thisVar<-paste(thisCode, "_Time_Spawn", sep=""); temp<-biolLines[grep(thisVar, biolLines)]
  Time_Spawns[g]<-get_first_number(temp)
  thisVar<-paste(thisCode, "_spawn_period", sep=""); temp<-biolLines[grep(thisVar, biolLines)]
  spawn_period[g]<-get_first_number(temp)
}
pdf(paste(plotPath,"RecruitTimings.pdf", sep=""), width=5)
par(las=1, lend=1)
plot(x=seq(1,365), y=rep(nag, 365), type="n", ylim=c(0,nag), yaxt="n", xlab="Day", ylab="")
axis(at=seq(1,nag), labels = as.character(groupsAgeDF$Code), side=2)
for(g in 1:nag){
  points(x=c(Time_Spawns[g], (Time_Spawns[g]+spawn_period[g])), y=rep(g,2), col=myBlue, lwd=2, type="l")
  if(Time_Spawns[g]+spawn_period[g] >365){ 
    secondPart<-(Time_Spawns[g]+spawn_period[g]) %% 365
    points(x=c(0, secondPart), y=rep(g,2), col=myPurple, lwd=2, type="l")
  }
}

logLines<-readLines(paste(outPath,"log.txt", sep=""))
logLines2<-readLines(paste(outPath,"log2.txt", sep=""))

testLines<-logLines[grep("External recruitment scaling of group BAL", logLines)]


getBox<-function(x){
  y<-unlist(str_split(x,"box"))[2]
  yy<-as.double(unlist(str_split(y,":"))[1])
  return(yy)
}

boxes<-unlist(lapply(testLines, getBox))
vdistribs<-unlist(lapply(testLines, get_first_number, n=4))

tapply(vdistribs, boxes, sum)


test_hdistrib<-c(0.1,0.15,0.2,0.05,0.5)
test_boxes<-c(1,2,2,3,3,3,4,5,5,5,5)
test_vdistrib<-c(1,0.8,0.2,0.8,0.1,0.1,1,0.4,0.3,0.2,0.1)
tapply(test_vdistrib, test_boxes, sum)

testRecs<-1000

recByCell<-array(NA,dim=c(4,5))
for(b in 1:5){
  for(l in 1:4){
    thisHprop<-test_hdistrib[b]
    temp<-test_vdistrib[test_boxes==b]
    thisVprop<-0
    if(length(temp)>=l){
      thisVprop<-temp[l]
    }
    thisProp<-thisVprop * thisHprop
    recByCell[l,b]<-thisProp*testRecs
  }
}

getBox<-function(x){
  xx<-unlist(str_split(x,"box"))[2]
  thisBox<-as.double(unlist(str_split(xx,"-"))[1])
  return(thisBox)
}
getHdistrib<-function(x){
  xx<-unlist(str_split(x,"hdistrib:"))[2]
  xx<-gsub(")","",xx)
  thisHdistrib<-as.double(xx)
  return(thisHdistrib)
}

thisCode<-"BAL"
for(g in 1:ng){
  ratio<-NA
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    temp<-logLines[grep("CHANGING NUMERATER", logLines)]
    testLines<-temp[grep(thisCode, temp)]
    
    temp<-logLines2[grep("CHANGING NUMERATER", logLines2)]
    testLines2<-temp[grep(thisCode, temp)]
    
    
    
    temp<-logLines2[grep("Bulkrecruits:", logLines2)]
    testLines<-temp[grep(thisCode, temp)]
    allBulkRecruits<-unique(unlist(lapply(testLines, get_first_number, n=5)))
    
    allNumRecs<-unlist(lapply(testLines, get_first_number, n=4)); allTimes<-unlist(lapply(testLines,get_first_number, n=1))
    
    # Time: 2.910000e+02, BAL box24-0 ngene: 0, qid: 1, num_rec: 2.986894e-02 (Bulkrecruits: 2.030639e+01, 
    # enviro_scalar: 1.000000e+00, vdistrib: 1.000000e-02, hdistrib: 1.470913e-01
    allBoxes<-unlist(lapply(testLines, getBox))
    numByBox<-tapply(allNumRecs, allBoxes, sum)
    
    vertDistrb<-unlist(lapply(testLines, get_first_number, n=7))
    horizDistrib<-unlist(lapply(testLines, getHdistrib))
    hDistribByBox<-tapply(horizDistrib, allBoxes, unique)
    
    testLines2<-logLines[grep(paste(thisCode," has larvae for qid",sep=""), logLines)]
    allLarvae<-unique(unlist(lapply(testLines2, get_first_number, n=3)))
   
    # Time: 2.920000e+02 BAL has larvae for qid: 1 set to 5.415037e+01 from totrecruit 5.41503703791823411962e+01
    
    ratio<-allLarvae[1]/allBulkRecruits[1]
    if(ratio != 1){
      cat(as.character(thisCode), ", ", ratio,"\n")
    }
  }
}

thisKWRR<-6288289; thisKWSR<-2445446; recWeight<-(thisKWRR+thisKWSR)*(mg_2_tonne*X_CN)

#check out larave starts
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    temp<-logLines[grep(paste(thisCode," get settlers", sep=""), logLines)]
    thisLarvae<-get_first_number(temp[1],n=6); thisBulkRec<-get_first_number(temp[1],7)
    if(thisLarvae>0 & thisBulkRec==0){cat(as.character(thisCode), ", ", thisLarvae,"\n")}
  }
}

