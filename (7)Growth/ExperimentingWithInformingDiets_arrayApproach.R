#The idea: inform each predator/prey paring with a probability distribution which is informed by trophic level, spatial overlap and gape size
#Do multiple (as in 1000's?) of predator/prey matrices, where for each pairing you sample from its distribution (if there is no chance of a predator eating a prey, there is no distribution, it is just zero for that pairing. otherwise, there is a distribution which reflects how likely it is the predator eats the prey)
#calculate a 'likelihood' for each matrix that compares realised growth with desired, and same with natural mortality (ie the growth and mortality that would be required to retain initial conditions)
#can do outputs showing which predators or preys have the most trouble being fitted to; which ones it doesn't really matter
#can test which are driving it most - TL, spatial overlap, gapesize

#will need to handle canabilism - they are likely within the model to be in the same space due to coarse spatial scale

source(paste(DIR$'General functions',"get_interaction_trophicLevel.R",sep=""))
source(paste(DIR$'General functions',"get_first_number.R",sep=""))
source(paste(DIR$'General functions',"myRounding.R",sep=""))
source(paste(DIR$'General functions',"getSNRN.R",sep=""))
source(paste(DIR$'General functions',"nonZeroMean.R",sep=""))
source(paste(DIR$'General functions',"getICNumbers.R",sep=""))

source(paste(DIR$'General functions',"get_interaction_spatial_bySpace.R",sep=""))
source(paste(DIR$'General functions',"get_interaction_gapeSize_byCohort.R",sep=""))

this_run<-"Base"

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
plotPath<-paste(this_path,"Figures\\",this_run,"\\",sep="")

##we'll need the groups df 
groupsDF<-read.csv(paste(this_path,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
thisInitialConditionsFile<-paste(this_path,"CRAM_input_short.nc",sep="")
ThisIC.nc<-nc_open(thisInitialConditionsFile)
thisBiolFile<-paste(this_path,"CRAM_BH_hybrid_biol.prm",sep="")
biolLines<-readLines(thisBiolFile)
groupsTL<-read.csv(paste(this_path,"inputs\\biol_prm\\CRAM_trophicLevels_isotopes.csv",sep=""))

#get volume from ICs
volume<-ncvar_get(ThisIC.nc,"volume")[,,1]
nlayers<-dim(volume)[1]; nboxes<-dim(volume)[2]

#######################################
#first: the spatial overlap
#set up df
#for groups that are age-structured, need an adult and a juvenile column and row
index<-groupsDF$NumCohorts>1
unagedGroups<-as.character(groupsDF$Code[!index])
agedGroups<-paste(rep(groupsDF$Code[index],2),c("juv","ad"),sep="")
ppGroups<-sort(c(unagedGroups,agedGroups)); npg<-length(ppGroups)

#output ppGroups
write.csv(ppGroups,paste(DIR$'Tables',"SampledAvailabilities\\ppGroups.csv",sep=""),row.names = FALSE)

interaction_gapeSize_array<-array(NA,dim=c(ncol=ng, 10, nrow=ng, 10)); 

interaction_trophicLevel_array<-array(NA,dim=c(npg,npg)) #trophic level is same over space

interaction_spatial_array<-array(NA,dim=c(npg, npg, nlayers, nboxes))

## SPATIAL
for(i in 1:npg){
  thisPredVar<-ppGroups[i]
  cat(thisPredVar,"--")
  thisPredCode<-gsub('juv|ad',"",thisPredVar)
  thisPredAge<-gsub(thisPredCode,"",thisPredVar)
  for(j in 1:npg){
    thisPreyVar<-ppGroups[j]
    thisPreyCode<-gsub('juv|ad',"",thisPreyVar)
    thisPreyAge<-gsub(thisPreyCode,"",thisPreyVar)
    
    thisInt<-get_interaction_spatial_bySpace(predatorCode=thisPredCode,preyCode=thisPreyCode,predatorAge=thisPredAge,preyAge=thisPreyAge,groupsDF=groupsDF,ThisIC.nc=ThisIC.nc,plotPath=plotPath,biolLines=biolLines)
    
    interaction_spatial_array[i,j,,]<-thisInt
   }
}
#write it out
save(list=c("interaction_spatial_array"),file=paste(DIR$"Data","\\eating\\interaction_spatial_array",sep=""))
# load(paste(DIR$"Data","\\eating\\interaction_spatial_array",sep="")); ##to bring it back

##GAPE SIZE
for(i in 1:ng){
  thisPredCode<-groupsDF$Code[i]
  cat(as.character(thisPredCode),"--")
  for(j in 1:ng){
    thisPreyCode<-groupsDF$Code[j]
     ## test is dim prey cohorts, predator cohorts
    thisInt<-get_interaction_gapeSize_byCohort(predatorCode=thisPredCode,preyCode=thisPreyCode,groupsDF=groupsDF,ThisIC.nc=ThisIC.nc,plotPath=plotPath,biolLines=biolLines)
    interaction_gapeSize_array[i,,j,]<-thisInt
   }
}
#write it out
save(list=c("interaction_gapeSize_array"),file=paste(DIR$"Data","\\eating\\interaction_gapeSize_array",sep=""))
# load(paste(DIR$"Data","\\eating\\interaction_gapeSize_array",sep="")); ##to bring it back


##  TROPHIC LEVEL
for(i in 1:npg){
  thisPredVar<-ppGroups[i]
  cat(thisPredVar,"--")
  thisPredCode<-gsub('juv|ad',"",thisPredVar)
  thisPredAge<-gsub(thisPredCode,"",thisPredVar)
  for(j in 1:npg){
    thisPreyVar<-ppGroups[j]
    thisPreyCode<-gsub('juv|ad',"",thisPreyVar)
    thisPreyAge<-gsub(thisPreyCode,"",thisPreyVar)
    thisInt<-get_interaction_trophicLevel(predatorCode=thisPredCode,preyCode=thisPreyCode,predatorAge=thisPredAge,preyAge=thisPreyAge,groupsTL)
    interaction_trophicLevel_array[i,j]<-thisInt
  }
}
save(list=c("interaction_trophicLevel_array"),file=paste(DIR$"Data","\\eating\\interaction_trophicLevel_array",sep=""))
# load(paste(DIR$"Data","\\eating\\interaction_trophicLevel_array",sep="")); ##to bring it back

##test ZS 
zsIndex<-ppGroups=="ZS"
zsIntPred<-interaction_spatial_array[,zsIndex,,]; zsPredsInt<-apply(zsIntPred,1,sum)
zsIntPrey<-interaction_spatial_array[zsIndex,,,]; zsPreysInt<-apply(zsIntPrey,1,sum)


########################
#summarise gape size and spatial by npg
gapeSize_npg<-array(NA,dim=c(npg,npg)); 
for(i in 1:npg){
  thisPredVar<-ppGroups[i]
  cat(thisPredVar,"--")
  thisPredCode<-gsub('juv|ad',"",thisPredVar)
  thisPredAge<-gsub(thisPredCode,"",thisPredVar)
  thisPredIndex<-groupsDF$Code==thisPredCode
  thisPredNumCohorts<-groupsDF$NumCohorts[thisPredIndex]
  thisPredCohortIndex<-1
  if(thisPredNumCohorts>1){
    ##get age mature for predator, then index cohort for this predator age
    thisPredMat<-get_first_number(biolLines[grep(paste(thisPredCode,"_age_mat",sep=""),biolLines)]) 
    if(thisPredAge=="juv"){
      thisPredCohortIndex<-seq(1,thisPredMat)
    } else{
      thisPredCohortIndex<-seq((thisPredMat+1),thisPredNumCohorts)
    }
  }
  for(j in 1:npg){
    thisPreyVar<-ppGroups[j]
    thisPreyCode<-gsub('juv|ad',"",thisPreyVar)
    thisPreyAge<-gsub(thisPreyCode,"",thisPreyVar)
    thisPreyIndex<-groupsDF$Code==thisPreyCode
    thisPreyNumCohorts<-groupsDF$NumCohorts[thisPreyIndex]
    thisPreyCohortIndex<-1
    if(thisPreyNumCohorts>1){
      ##get age mature for predator, then index cohort for this predator age
      thisPreyMat<-get_first_number(biolLines[grep(paste(thisPreyCode,"_age_mat",sep=""),biolLines)]) 
      if(thisPreyAge=="juv"){
        thisPreyCohortIndex<-seq(1,thisPreyMat)
      } else{
        thisPreyCohortIndex<-seq((thisPreyMat+1),thisPreyNumCohorts)
      }
    }
    temp_gape<-interaction_gapeSize_array[thisPredIndex,thisPredCohortIndex,thisPreyIndex,thisPreyCohortIndex]
    potential_gape<-length(thisPredCohortIndex)*length(thisPreyCohortIndex)
    thisPropGape<-sum(temp_gape,na.rm=TRUE)/potential_gape
    gapeSize_npg[i,j]<-thisPropGape
    if(thisPreyNumCohorts==1 | thisPredNumCohorts==1){
      gapeSize_npg[i,j]<-1
    }
  }
}

potential_spatial<-nlayers*nboxes
spatial_npg<-apply(interaction_spatial_array,c(1,2),sum,na.rm=TRUE)/potential_spatial

# totalInteraction<-(spatial_npg+gapeSize_npg+interaction_trophicLevel_array)/3
totalInteraction<-(spatial_npg*gapeSize_npg*interaction_trophicLevel_array)

maxPP<-0.03
scaledInteraction<-(totalInteraction/max(totalInteraction))*maxPP

#read in original pPREY interactions (from Peter Horn) and add 0.05 to these
thisAddAvail<-0.05
origDiets<-read.csv(paste(this_path,"\\inputs\\biol_prm\\pPREY\\CRAM_biol_pPREY_original.csv",sep=""))
for(j in 2:ncol(origDiets)){
  thisVar<-colnames(origDiets)[j]
  predAge<-get_first_number(thisVar,n=2); preyAge<-get_first_number(thisVar)
  thisCode<-gsub(paste("pPREY|",predAge,"|",preyAge,sep=""),"",thisVar)
  thisPpred<-thisCode
  thisPredCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisCode]
  if(thisPredCohorts>1){
    if(predAge==1){
      thisPpred<-paste(thisCode,"juv",sep="")
    } else{
      thisPpred<-paste(thisCode,"ad",sep="")
    }
  }
  temp<-origDiets[,j]; index<-temp>0
  thisAvails<-temp[index]; thisPreys<-c(groupsDF$Code,c("DC","DL","DR"))[index]
  thisPpreys<-thisPreys
  if(length(thisPreys)>0){
    for(i in 1:length(thisPreys)){
      thisPrey<-thisPreys[i]
      thisPreyCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisPrey]
      if(thisPreyCohorts>1){
        if(preyAge==1){
          thisPpreys[i]<-paste(thisPrey,"juv",sep="")
        }else{
          thisPpreys[i]<-paste(thisPrey,"ad",sep="")
        }
      }
    }
  }
  colIndex<-ppGroups %in% thisPpreys; rowIndex<-ppGroups==thisPpred
  scaledInteraction[rowIndex,colIndex]<-scaledInteraction[rowIndex,colIndex]+thisAddAvail
  
}

# sampledAvails<-totalInteraction
sampledAvails<-scaledInteraction


#output the sampledAvails, then read it in (in another script) and use it + test wrt desired growth and mortality
write.csv(sampledAvails,paste(DIR$'Tables',"SampledAvailabilities\\ExpectedValues.csv",sep=""),row.names = FALSE)

#numbers are quite high - do a reduced version
# sampledAvails<-signif(totalInteraction*0.001,2)
# 
# write.csv(sampledAvails,paste(DIR$'Tables',"SampledAvailabilities\\ExpectedValues_scaledDown.csv",sep=""),row.names = FALSE)
# 
# sampledAvails<-signif(totalInteraction*10,2)
# 
# write.csv(sampledAvails,paste(DIR$'Tables',"SampledAvailabilities\\ExpectedValues_scaledUp.csv",sep=""),row.names = FALSE)
# 
