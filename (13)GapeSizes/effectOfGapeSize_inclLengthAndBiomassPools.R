#read in tracers from the base model and store ready to check out gape sizes
this_run<-"base"
this_out<-"Base"
source(paste(DIR$'General functions',"nonZeroMean.R", sep=""))
source(paste(DIR$'General functions',"calc_length.R", sep=""))

mg_2_tonne<-2e-8; X_CN<-5.7

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")

plotPath<-paste(DIR$'Figures',"GapeSizes\\",sep="")

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

storeGapeInteractions<-array(NA, dim=c(ng,10, ng,10)); StoreGapeInteractionEaten<-0*storeGapeInteractions

## store size (roughly) of biomass pool species - these are just approx, but sufficient for this purpose i think,. 
## Based on Fig by Peter for TBGB (Approx_size_ranges.jpg)
bpGroups<-as.character(groupsDF$Code[groupsDF$NumCohorts==1]); nbp<-length(bpGroups)
# bpSizes<-data.frame(matrix(NA, ncol=2, nrow=nbp)); colnames(bpSizes)<-c("Code","Size")
# bpSizes$Code<-bpGroups
# write.csv(bpSizes, paste(DIR$'Tables',"BiomassPoolIndividualSizes.csv", sep=""), row.names=FALSE)
bpSizes<-read.csv(paste(DIR$'Tables',"BiomassPoolIndividualSizes.csv", sep=""))

predatorMeanLengths<-array(NA, dim=c(ng, 10))

for(g in 1:ng){
  thisCode<-as.character(groupsDF$Code[g]); predNumCohorts<-groupsDF$NumCohorts[g]
  if(predNumCohorts>1 & groupsDF$IsPredator[g]==1){
    #loop through the predator age-classes (aka cohorts) and check for gapesize-preysize overlap
    for(c in 1:predNumCohorts){
      tempDiet<-dietSummary[grep(thisCode,predators),c,]
      thisPrey<-rev(sort(tempDiet[tempDiet>0]))
      # if(length(thisPrey)>0){
        thisKLP<-storeGapeScalars$KLP[g]; thisKUP<-storeGapeScalars$KUP[g]
        this_a<-l2wpars$a[l2wpars$Code==thisCode]; this_b<-l2wpars$b[l2wpars$Code==thisCode]
        thisSN<-storeSNtracers[,g,c]
        thisRN<-storeSNtracers[,g,c]
        thisWeight<-(thisRN + thisSN) * mg_2_tonne * X_CN
        thisLength<-unlist(lapply(thisWeight,  calc_length, this_a=this_a, this_b=this_b))
        predatorMeanLengths[g,c]<-nonZeroMean(thisLength)
        gapeRange<-c(min(thisSN) * thisKLP, max(thisSN) * thisKUP)
        
        # get prey that make up at least 95% of diet
        x<-cumsum(thisPrey); xx<-seq(1,length(thisPrey))[x == min(x[x>0.95])]
      thisTopPrey<-thisPrey[1:xx]
      #then just keep the prey that are age-structured
      # thisPrey<-thisPrey[names(thisPrey) %in% as.character(groupsDF$Code[groupsDF$NumCohorts>1])]
      checkPreys<-thisPrey; npreys<-length(checkPreys)
      if(npreys>0){
        for(prey in 1:npreys){
          #get the size structure of the prey - just do the first one for now
          thisPrey<-names(checkPreys)[prey]
          preyProp<-checkPreys[prey]
          prey_g<-grep(thisPrey, groupsDF$Code); preyNumCohorts<-groupsDF$NumCohorts[prey_g]
          if(preyNumCohorts==1){
            preyLength<-bpSizes$Size[bpSizes$Code==thisPrey]
            gapeProp<-1
            storeGapeInteractions[g,c, prey_g,prey_c]<-gapeProp
            StoreGapeInteractionEaten[g,c,prey_g,prey_c]<-gapeProp * preyProp * preyLength
            
          }
          if(preyNumCohorts>1){
            preyWeights<-(storeRNtracers[,prey_g,] + storeSNtracers[,prey_g,]) * mg_2_tonne * X_CN
            prey_a<-l2wpars$a[l2wpars$Code==thisPrey]; prey_b<-l2wpars$b[l2wpars$Code==thisPrey]
            if(length(prey_a)>0){
              preyLengths<-apply(preyWeights, c(1,2), calc_length, this_a=prey_a, this_b=prey_b)
              
               #for each prey age-class, which ones are within this range (and always of most of the time..?)
              for(prey_c in 1:preyNumCohorts){
                prey_SN<-storeSNtracers[,prey_g,prey_c]
                gapeIndex<-prey_SN >= gapeRange[1] & prey_SN<=gapeRange[2]
                gapeProp<-sum(gapeIndex)/length(gapeIndex)
                preyLength<-nonZeroMean(preyLengths[,prey_c])
                storeGapeInteractions[g,c, prey_g,prey_c]<-gapeProp
                StoreGapeInteractionEaten[g,c,prey_g,prey_c]<-(gapeProp * preyProp * preyLength)/preyNumCohorts
              }
            } 
            # else{
            #   cat("no length for ",thisPrey,"\n")
            # }
          }
        }
      }
    }
  }
}

meanPreyLengthByPred<-apply(StoreGapeInteractionEaten,c(1,2), sum, na.rm=TRUE)
propPreyPredLength<-meanPreyLengthByPred/predatorMeanLengths

for(g in 1:ng){
  thisCode<-as.character(groupsDF$Code[g]); thisNumCohorts<-groupsDF$NumCohorts[g]; 
  if(groupsDF$IsPredator[g]==1){
    thisProps<-propPreyPredLength[g,1:thisNumCohorts]
    if(sum(thisProps, na.rm=TRUE)>0){
      par(lend=1)
      plot(x=seq(1,thisNumCohorts), y=thisProps, type="h", lwd=5, col=myBlue, xlim=c(0,thisNumCohorts+1), xaxt="n", xlab="", ylab="Prey / predator length")
      axis(at=seq(1,thisNumCohorts), labels=seq(1,thisNumCohorts), side=1)
      mtext("Age class", side=1, adj=0.5, line=2)
      mtext(gsub("_"," ", groupsDF$Name[g]), side=3, adj=0)
    }
  }
}
## how do these compare to predator gape/length?
#read in gape size data
load(paste(DIR$'Base',"data\\gapeSize\\fish_bio.RData",sep=""))
# brings in fish_bio
thisBioData<-fish_bio[fish_bio$Code %in% groupsDF$Code,]
plot(1,xlim=c(0,0.5), ylim=c(0,0.5), type="n", xlab="Gape size / length", ylab="Prey length / predator length")
for(g in 1:ng){
  thisCode<-as.character(groupsDF$Code[g]); 
  thisNumCohorts<-groupsDF$NumCohorts[g]; 
  if(groupsDF$IsPredator[g]==1){
    predBioData<-thisBioData[thisBioData$Code==thisCode,]
    predBioProps<-as.double(predBioData$oral.gape)/as.double(predBioData$lgth)
    atLengthProps<-propPreyPredLength[g,]
    points(x=rep(mean(predBioProps), length(atLengthProps)), y=(atLengthProps), pch=20, col=myBlue_trans, cex=1.5)
  }
}
points(x=c(0,1), y=c(0,1), type="l", lwd=2, col="red", lty=4)

# now have length on x axis and compare prey/pred length with gape/length both on the y-axis
maxLength<-max(predatorMeanLengths, na.rm=TRUE)
par(mar=c(4,6,1,1))
plot(1, ylim=c(0,0.5), xlim=c(0,maxLength), type="n", xlab="Length (cm)", ylab="Gape / length\nPrey / predator length")
for(g in 1:ng){
  thisCode<-as.character(groupsDF$Code[g]); 
  thisNumCohorts<-groupsDF$NumCohorts[g]; 
  if(groupsDF$IsPredator[g]==1){
    predBioData<-thisBioData[thisBioData$Code==thisCode,]
    predBioLengths<-as.double(predBioData$lgth)
    predBioProps<-as.double(predBioData$oral.gape)/predBioLengths
    if(length(predBioProps)>0){
      points(x=predBioLengths, y=predBioProps, pch=20, col=myBlue_trans, cex=1.5)
    }
    atLengthProps<-propPreyPredLength[g,]
    atLengths<-predatorMeanLengths[g,]
    if(sum(atLengths,na.rm=TRUE)>0){
      points(x=atLengths, y=atLengthProps, pch=20, col=myOrange_trans, cex=1.5)
    }
  }
}

maxLength<-155
pdf(paste(plotPath,"GapePropByLengthALLSpecies.pdf", sep=""), width=6, height=4)
par(mar=c(4,6,1,1))
plot(1, ylim=c(0,0.5), xlim=c(0,maxLength), type="n", xlab="Length (cm)", ylab="Gape / length\nPrey / predator length")
for(g in 1:ng){
  thisCode<-as.character(groupsDF$Code[g]); 
  thisNumCohorts<-groupsDF$NumCohorts[g]; 
  if(groupsDF$IsPredator[g]==1){
    predBioData<-thisBioData[thisBioData$Code==thisCode,]
    predBioLengths<-as.double(predBioData$lgth)
    predBioProps<-as.double(predBioData$oral.gape)/predBioLengths
    atLengthProps<-propPreyPredLength[g,]
    atLengths<-predatorMeanLengths[g,]
    if(length(predBioProps)>0 & sum(atLengths,na.rm=TRUE)>0){
      points(x=predBioLengths, y=predBioProps, pch=20, col=myBlue_trans, cex=1.5)
      points(x=atLengths, y=atLengthProps, pch=20, col=myOrange_trans, cex=1.5)
    }
  }
}
legend(legend=c("Gape / length","Prey / predator length"), pch=20, pt.cex=2, cex=1, col=c(myBlue_trans, myOrange_trans), x="topright", bty="n")
dev.off()

pdf(paste(plotPath,"GapePropByLengthAndSpecies.pdf", sep=""), width=10)
par(mfrow=c(4,4), mar=c(4,6,1.5,0.5))
for(g in 1:ng){
  thisCode<-as.character(groupsDF$Code[g]); 
  thisNumCohorts<-groupsDF$NumCohorts[g]; 
  if(groupsDF$IsPredator[g]==1){
    predBioData<-thisBioData[thisBioData$Code==thisCode,]
    predBioLengths<-as.double(predBioData$lgth)
    predBioProps<-as.double(predBioData$oral.gape)/predBioLengths
    atLengthProps<-propPreyPredLength[g,]
    atLengths<-predatorMeanLengths[g,]
    if(length(predBioProps)>0 & sum(atLengths,na.rm=TRUE)>0){
      plot(1, ylim=c(0,0.5), xlim=c(0,maxLength), type="n", xlab="Length (cm)", ylab="Gape / length\nPrey / predator length")
      points(x=atLengths, y=atLengthProps, pch=20, col=myOrange_trans, cex=2)
      points(x=predBioLengths, y=predBioProps, pch=20, col=myBlue_trans, cex=2)
      mtext(gsub("_"," ", groupsDF$Name[g]), side=3, adj=0, cex=0.8)
    }
  }
}
makeBlankPlot()
legend(legend=c("Gape / length","Prey / predator length"), pch=20, pt.cex=2, cex=1, col=c(myBlue_trans, myOrange_trans), x="center", bty="n")
dev.off()



