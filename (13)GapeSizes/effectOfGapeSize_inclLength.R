#read in tracers from the base model and store ready to check out gape sizes
this_run<-"base"
this_out<-"Base"
source(paste(DIR$'General functions',"nonZeroMean.R", sep=""))
source(paste(DIR$'General functions',"calc_length.R", sep=""))

mg_2_tonne<-2e-8; X_CN<-5.7

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")

## read them in  =c("storeTonnesTracers", "storeRNtracers", "storeSNtracers", "storeNumtracers", "dietSummary", "storeGapeScalars)
storeTracersFile<-paste(this_path,"GapeSizes\\baseModelTracers", sep="")
load(storeTracersFile)
# dimensions are: storeTonnesTracers (timesteps, groups); storeRNtracers, storeSNtracers, storeNumtracers (timesteps, group, age-class)
# dietSummary (predator, predator age-class, prey); colnames(storeGapeScalars)=(Code,KLP,KUP)

groupsDF<-read.csv(paste(this_path,"..\\CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

#length to weight pars
l2wpars<-read.csv(paste(this_path,"..\\inputs\\supporting\\length2weights.csv", sep=""))

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1900
fishingStartYear<-1900
modelStartYear<-1900

nts<-dim(storeTonnesTracers)[1]

predators<-as.character(groupsDF$Code[groupsDF$IsPredator==1])

storeGapeInteractions<-array(NA, dim=c(ng,10, ng,10))

for(g in 1:ng){
  thisCode<-as.character(groupsDF$Code[g]); predNumCohorts<-groupsDF$NumCohorts[g]
  if(predNumCohorts>1 & groupsDF$IsPredator[g]==1){
    #loop through the predator age-classes (aka cohorts) and check for gapesize-preysize overlap
    for(c in 1:predNumCohorts){
      tempDiet<-dietSummary[grep(thisCode,predators),c,]
      thisPrey<-rev(sort(tempDiet[tempDiet>0]))
      if(length(thisPrey)>0)
      # get prey that make up at least 95% of diet
      x<-cumsum(thisPrey); xx<-seq(1,length(thisPrey))[x == min(x[x>0.95])]
      thisTopPrey<-thisPrey[1:xx]
      #then just keep the prey that are age-structured
      thisPrey<-thisPrey[names(thisPrey) %in% as.character(groupsDF$Code[groupsDF$NumCohorts>1])]
      checkPreys<-thisPrey; npreys<-length(checkPreys)
      if(npreys>0){
        for(prey in 1:npreys){
          #get the size structure of the prey - just do the first one for now
          thisPrey<-names(checkPreys)[prey]
          prey_g<-grep(thisPrey, groupsDF$Code); preyNumCohorts<-groupsDF$NumCohorts[prey_g]
          if(preyNumCohorts>1 & groupsDF$IsPredator[prey_g]==1){
            preyWeights<-(storeRNtracers[,prey_g,] + storeSNtracers[,prey_g,]) * mg_2_tonne * X_CN
            prey_a<-l2wpars$a[l2wpars$Code==thisPrey]; prey_b<-l2wpars$b[l2wpars$Code==thisPrey]
            preyLengths<-apply(preyWeights, c(1,2), calc_length, this_a=prey_a, this_b=prey_b)
          
            thisKLP<-storeGapeScalars$KLP[g]; thisKUP<-storeGapeScalars$KUP[g]
            this_a<-l2wpars$a[l2wpars$Code==thisCode]; this_b<-l2wpars$b[l2wpars$Code==thisCode]
            thisSN<-storeSNtracers[,g,c]
            thisRN<-storeSNtracers[,g,c]
            thisWeight<-(thisRN + thisSN) * mg_2_tonne * X_CN
            thisLength<-unlist(lapply(thisWeight,  calc_length, this_a=this_a, this_b=this_b))
            
            gapeRange<-c(min(thisSN) * thisKLP, max(thisSN) * thisKUP)
            #for each prey age-class, which ones are within this range (and always of most of the time..?)
            for(prey_c in 1:preyNumCohorts){
              prey_SN<-storeSNtracers[,prey_g,prey_c]
              gapeIndex<-prey_SN >= gapeRange[1] & prey_SN<=gapeRange[2]
              gapeProp<-sum(gapeIndex)/length(gapeIndex)
              storeGapeInteractions[g,c, prey_g,prey_c]<-gapeProp
            }
          }
        }
      }
    }
  }
}

colByCohort<-colorRampPalette(colors=c(myGold, myGreen, myBlue, myRed,"red"))(10)
colByCohort_trans<-paste(colByCohort,"88", sep="")

# for hoki, say - which prey does it overlap with..?
g<-grep("LIN", groupsDF$Code); predNumCohorts<-groupsDF$NumCohorts[g]; thisCode<-as.character(groupsDF$Code[g])
test<-storeGapeInteractions[g,,,]
index<-apply(test,2,sum, na.rm=TRUE)>0
thisData<-test[,index,]
gapeInfo<-c(); 
myPrey<-groupsDF$Code[index]; nMyPreys<-length(myPrey)

plotPath<-paste(DIR$'Figures',sep="")
pdf(paste(plotPath,"gapeSizes_",thisCode,".pdf", sep=""))
par(mfrow=c(3,2), mar=c(4,4,2,1))
makeBlankPlot()
legend(legend=c("Within gape range","Outside gape range"), col=c(myBlue,"red"), pch=20, cex=1.5, x="center", bty="n")
makeBlankPlot()
legend(legend=seq(1,10),title="Predator\nage-class", col=colByCohort_trans, lwd=3, seg.len=3, x="left", bty="n")
legend(legend=c("Predator","Prey"), lty=c(1,2), col=c("black",myGrey),x="right",seg.len=3, bty="n")
for(prey in 1:nMyPreys){
  thisPrey<-as.character(myPrey[prey]); prey_g<-grep(thisPrey, groupsDF$Code); preyNumCohorts<-groupsDF$NumCohorts[prey_g]
  preyWeights<-(storeRNtracers[,prey_g,] + storeSNtracers[,prey_g,]) * mg_2_tonne * X_CN
  prey_a<-l2wpars$a[l2wpars$Code==thisPrey]; prey_b<-l2wpars$b[l2wpars$Code==thisPrey]
  preyLengths<-apply(preyWeights, c(1,2), calc_length, this_a=prey_a, this_b=prey_b)
  
  thisKLP<-storeGapeScalars$KLP[g]; thisKUP<-storeGapeScalars$KUP[g]
  this_a<-l2wpars$a[l2wpars$Code==thisCode]; this_b<-l2wpars$b[l2wpars$Code==thisCode]
  thisSN<-storeSNtracers[,g,]
  thisRN<-storeSNtracers[,g,]
  thisWeight<-(thisRN + thisSN) * mg_2_tonne * X_CN
  thisLength<-apply(thisWeight, c(1,2), calc_length, this_a=this_a, this_b=this_b)
  #grab the average length by cohort
  predLengthByAgeClass<-apply(thisLength,2,mean, na.rm=TRUE)
  preyLengthByAgeClass<-apply(preyLengths,2, mean, na.rm=TRUE)
  thisMin<-min(c(predLengthByAgeClass, preyLengthByAgeClass), na.rm=TRUE); thisMax<-max(c(predLengthByAgeClass, preyLengthByAgeClass), na.rm=TRUE)
  thisOverlap<-thisData[,prey,]
  ## in this case, cohorts 9 and 10 of hoki can eat cohorts 1 and 2 of squid - how does this look when we convert it to length?
  plot(x=preyLengthByAgeClass, y=preyLengthByAgeClass, type="n", ylim=c(thisMin, thisMax), ylab="Predator length (cm)", xlab="Prey length (cm)")
  for(i in 1:dim(thisOverlap)[1]){
    for(j in 1:dim(thisOverlap)[2]){
      thisCol<-"red"
      if(!is.na(thisOverlap[i,j])){
        if(thisOverlap[i,j]>0){
          thisCol<-myBlue
          #store the length prey/pred ratio
          gapeInfo<-c(gapeInfo,preyLengthByAgeClass[j]/predLengthByAgeClass[i])
        }
      }
      points(x=preyLengthByAgeClass[j], y=predLengthByAgeClass[i], pch=20, cex=1.5, col=thisCol)
    }
  }
  mtext(thisPrey,side=3,adj=0)
  plot(x=seq(1,predNumCohorts), y=predLengthByAgeClass, type="n", ylim=c(thisMin, thisMax), xlab="Age class", ylab="Length (cm)")
  #join the cohorts that are within range
  for(i in 1:dim(thisOverlap)[1]){
    for(j in 1:dim(thisOverlap)[2]){
      if(!is.na(thisOverlap[i,j])){
        if(thisOverlap[i,j]>0){
          # points(x=preyLengthByAgeClass[j], y=predLengthByAgeClass[i], type="l",col=myBlue_trans)
          segments(x0=j, y0=preyLengthByAgeClass[j], x1=i, predLengthByAgeClass[i], col=colByCohort_trans[i], lwd=3)
        } 
      }
    }
  }
  points(x=seq(1,preyNumCohorts), y=preyLengthByAgeClass[1:preyNumCohorts], type="l", lwd=2)
  points(x=seq(1,predNumCohorts), y=predLengthByAgeClass, type="l",lwd=2,lty=2)
}
  
  
dev.off()

## look at it as a heatmap
getColor<-function(x,thisMax){
  thisCol<-"white"
  if(!is.na(x) & x>0){
    y<-round((x)/(thisMax),2)*100+1
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

prey=3
thisPrey<-as.character(myPrey[prey]); prey_g<-grep(thisPrey, groupsDF$Code); preyNumCohorts<-groupsDF$NumCohorts[prey_g]
preyWeights<-(storeRNtracers[,prey_g,] + storeSNtracers[,prey_g,]) * mg_2_tonne * X_CN
prey_a<-l2wpars$a[l2wpars$Code==thisPrey]; prey_b<-l2wpars$b[l2wpars$Code==thisPrey]
preyLengths<-apply(preyWeights, c(1,2), calc_length, this_a=prey_a, this_b=prey_b)

#grab the average length by cohort
preyLengthByAgeClass<-apply(preyLengths,2, mean, na.rm=TRUE)
thisMin<-min(c(predLengthByAgeClass, preyLengthByAgeClass), na.rm=TRUE); thisMax<-max(c(predLengthByAgeClass, preyLengthByAgeClass), na.rm=TRUE)
thisOverlap<-thisData[,prey,]

plotData<-thisOverlap[1:predNumCohorts, 1:preyNumCohorts]
thisMax<-max(plotData, na.rm=TRUE); thisMin<-min(plotData, na.rm=TRUE)

plotColour<-apply(plotData,c(1,2),getColor,thisMax)
tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotData)[1]),dim(plotData)[2]),"y"=sort(rep(seq(1,dim(plotData)[2]),dim(plotData)[1]))))
pdf(paste(plotPath,"GapesizeByLength",thisCode,"_",thisPrey,".pdf",sep=""),height=4,width=10)
par(mar=c(6,4,1.5,1))
plot(x=seq(1,dim(plotData)[1]),y=rep(dim(plotData)[2],dim(plotData)[1]),type="n",xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0.5,dim(plotData)[2]+0.5))
axis(at=seq(1,dim(plotData)[1]),labels = seq(1,predNumCohorts),side=1,las=2)
axis(at=seq(1,dim(plotData)[2]),labels=seq(1,preyNumCohorts),side=2,las=1)
temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
box()
mtext(paste(thisCode," eating ",thisPrey, sep=""), side=3,adj=0)
mtext("Predator age class", side=1, adj=0.5,line=2.5)
mtext("Prey age class", side=2, adj=0.5,line=2.5)
dev.off()

load(paste(DIR$'Base',"data\\gapeSize\\fish_bio.RData",sep=""))
# brings in fish_bio
thisData<-fish_bio[fish_bio$Code=="LIN",]

