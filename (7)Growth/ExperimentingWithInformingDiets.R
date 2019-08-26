#The idea: inform each predator/prey paring with a probability distribution which is informed by trophic level, spatial overlap and gape size
#Do multiple (as in 1000's?) of predator/prey matrices, where for each pairing you sample from its distribution (if there is no chance of a predator eating a prey, there is no distribution, it is just zero for that pairing. otherwise, there is a distribution which reflects how likely it is the predator eats the prey)
#calculate a 'likelihood' for each matrix that compares realised growth with desired, and same with natural mortality (ie the growth and mortality that would be required to retain initial conditions)
#can do outputs showing which predators or preys have the most trouble being fitted to; which ones it doesn't really matter
#can test which are driving it most - TL, spatial overlap, gapesize

#will need to handle canabilism - they are likely within the model to be in the same space due to coarse spatial scale

source(paste(DIR$'General functions',"get_interaction_spatial.R",sep=""))
source(paste(DIR$'General functions',"get_interaction_gapeSize.R",sep=""))
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
thisBiolFile<-paste(this_path,"CRAM_base_biol.prm",sep="")
biolLines<-readLines(thisBiolFile)
groupsTL<-read.csv(paste(this_path,"inputs\\biol_prm\\CRAM_trophicLevels_isotopes.csv",sep=""))
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

interaction_gapeSize_df<-data.frame(matrix(NA,ncol=npg,nrow=npg)); colnames(interaction_gapeSize_df)<-(as.character(ppGroups)); rownames(interaction_gapeSize_df)<-as.character(ppGroups)

interaction_trophicLevel_df<-interaction_gapeSize_df

interaction_spatial_df<-interaction_gapeSize_df
# interaction_spatial_df<-data.frame(matrix(NA,ncol=npg,nrow=ng)); colnames(interaction_df)<-(as.character(groupsDF$Code)); rownames(interaction_df)<-as.character(groupsDF$Code)

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
    
    # test<-get_interaction_spatial_bySpace(predatorCode=thisPredCode,preyCode=thisPreyCode,predatorAge=thisPredAge,preyAge=thisPreyAge,groupsDF=groupsDF,ThisIC.nc=ThisIC.nc,plotPath=plotPath,biolLines=biolLines)
    
    thisInt<-get_interaction_spatial(predatorCode=thisPredCode,preyCode=thisPreyCode,predatorAge=thisPredAge,preyAge=thisPreyAge,groupsDF=groupsDF,ThisIC.nc=ThisIC.nc,plotPath=plotPath,biolLines=biolLines)
    
    interaction_spatial_df[i,j]<-thisInt
  }
}


 ##GAPE SIZE
for(i in 1:npg){
  thisPredVar<-ppGroups[i]
  cat(thisPredVar,"--")
  thisPredCode<-gsub('juv|ad',"",thisPredVar)
  thisPredAge<-gsub(thisPredCode,"",thisPredVar)
  for(j in 1:npg){
    thisPreyVar<-ppGroups[j]
    thisPreyCode<-gsub('juv|ad',"",thisPreyVar)
    thisPreyAge<-gsub(thisPreyCode,"",thisPreyVar)
    
    ## test is dim prey cohorts, predator cohorts
    # test<-get_interaction_gapeSize_byCohort(predatorCode=thisPredCode,preyCode=thisPreyCode,predatorAge=thisPredAge,preyAge=thisPreyAge,groupsDF=groupsDF,ThisIC.nc=ThisIC.nc,plotPath=plotPath,biolLines=biolLines)
    
    thisInt<-get_interaction_gapeSize(predatorCode=thisPredCode,preyCode=thisPreyCode,predatorAge=thisPredAge,preyAge=thisPreyAge,groupsDF=groupsDF,ThisIC.nc=ThisIC.nc,plotPath=plotPath,biolLines=biolLines)
    interaction_gapeSize_df[i,j]<-thisInt
  }
}

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
    interaction_trophicLevel_df[i,j]<-thisInt
  }
}

#get sampled availabilities
sampledAvails<-0*interaction_gapeSize_df
for(i in 1:dim(sampledAvails)[1]){
  thisPredVar<-ppGroups[i]
  thisPredCode<-gsub('juv|ad',"",thisPredVar)
  thisPredAge<-gsub(thisPredCode,"",thisPredVar)
  
  sampledAvails[i,]<-sampleFromPrey(predatorCode=thisPredCode,predatorAge=thisPredAge)
    
}

#output the sampledAvails, then read it in (in another script) and use it + test wrt desired growth and mortality
write.csv(sampledAvails,paste(DIR$'Tables',"SampledAvailabilities\\ExpectedValues.csv",sep=""))
