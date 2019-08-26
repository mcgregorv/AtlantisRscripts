#reads in sampled availabilities, then analyses it wrt intended growth and mortality
source(paste(DIR$'General functions',"getAssumedNaturalMortality.R",sep=""))
source(paste(DIR$'General functions',"getRequiredFoodPerDay.R",sep=""))
source(paste(DIR$'General functions',"getMortalityRequiredPerDay.R",sep=""))
source(paste(DIR$'General functions',"getICNumbers.R",sep=""))
source(paste(DIR$'General functions',"get_first_number.R",sep=""))
source(paste(DIR$'General functions',"getPreyBiomassAvail.R",sep=""))
source(paste(DIR$'General functions',"getSNRN.R",sep=""))
source(paste(DIR$'General functions',"nonZeroMean.R",sep=""))
source(paste(DIR$'General functions',"getAvailFood.R",sep=""))

getICBiomass<-function(groupName,ThisIC.nc,volume){
  thisVar<-paste(groupName,"_N",sep=""); ICbiomass<-ncvar_get(ThisIC.nc,thisVar)
  
  if(length(dim(ICbiomass))==3){ICbiomass<-ICbiomass[,,1]}
  if(length(dim(ICbiomass))==2){
    if(ICbiomass[1,2]>9e+36){
      ICbiomass<-ICbiomass[,1]
    }
  }
  if(length(dim(ICbiomass))==0){
    temp<-ICbiomass; ICbiomass<-0*volume
    ICbiomass[dim(volume)[1],]<-temp
  }
  return(ICbiomass)
}

ppGroups<-read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ppGroups.csv",sep="")); npg<-length(ppGroups[,1])
sampledAvails<-read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ExpectedValues.csv",sep=""))[,-1]
rownames(sampledAvails)<-ppGroups[,1]

this_run<-"Base"

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
plotPath<-paste(this_path,"Figures\\",this_run,"\\",sep="")

##we'll need the groups df 
groupsDF<-read.csv(paste(this_path,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
thisInitialConditionsFile<-paste(this_path,"CRAM_input_short.nc",sep=""); ThisIC.nc<-nc_open(thisInitialConditionsFile)
thisBiolFile<-paste(this_path,"CRAM_base_biol.prm",sep="")
biolLines<-readLines(thisBiolFile)
volume<-ncvar_get(ThisIC.nc,"volume")[,,1]; dz<-ncvar_get(ThisIC.nc,"dz")[,,1]; Seabird<-ncvar_get(ThisIC.nc,"Seabird_N")[,,1]; OrangeR<-ncvar_get(ThisIC.nc,"Orange_roughy_N")[,,1]; nlayers<-dim(volume)[1]; nboxes<-dim(volume)[2]

MbyGroup<-unlist(lapply(groupsDF$Code,getAssumedNaturalMortality,ThisIC.nc=ThisIC.nc,groupsDF=groupsDF,biolLines=biolLines))
reqFoodByGroupCohort<-array(NA,dim=c(ng,10,dim(volume)))
for(g in 1:ng){
  groupCode<-groupsDF$Code[g]
  thisRF<-getRequiredFoodPerDay(groupCode,groupsDF,biolLines,ThisIC.nc)
  reqFoodByGroupCohort[g,1:length(thisRF),,]<-thisRF
}

#to come after getAvailFood():
storeDeaths<-array(0,dim=c(ng,10,dim(volume)))
storeAssimilated<-array(0,dim=c(ng,10,dim(volume)))


#########################################
##NOTES NOTES NOTES #################
## getAvailFood() is REALLY slow! and clunky. Tidy it up.
##
## more notes further down


#get total available prey for each predator and store, so can use proportion eaten
totPreyAvail<-array(NA,dim=c(ng,10,dim(volume)))
for(P in 1:npg){
  thisPredVar<-ppGroups[P,1]
  thisPredCode<-gsub('juv|ad',"",thisPredVar)
  thisPredAge<-gsub(thisPredCode,"",thisPredVar)
  thisPredName<-str_trim(groupsDF$Name[groupsDF$Code==thisPredCode],side="both")
  thisNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisPredCode]
  thisRowName<-paste(thisPredCode,thisPredAge,sep="")
  thisPreyAvail<-sampledAvails[grep(thisRowName,rownames(sampledAvails)),] #rows are predators, columns are prey
  nonZeroPreys<-colnames(thisPreyAvail)[as.double(thisPreyAvail)>0]; npreys<-length(nonZeroPreys) #get names on columns with nonzero entries
  cat("PRED ",thisPredCode,"--")
  
  predBiomassIC<-getICBiomass(groupName = thisPredName,ThisIC.nc, volume)
  for(b in 1:nboxes){
    for(l in 1:nlayers){
      #check if we have predator biomass here
      if(predBiomassIC[l,b]>0){
        #getAvailFood returns value for each cohort withing predAge
        totPreyAvail[groupsDF$Code==thisPreyCode,firstCohort:lastCohort,l,b]<-getAvailFood(groupCode=thisPredCode,predAge=thisPredAge,groupsDF=groupsDF,sampledAvails=sampledAvails,biolLines=biolLines,volume=volume,thisLayer=l,thisBox=b)
      }
    }
  }
}
  
  

for(P in 1:npg){
  thisPredVar<-ppGroups[P,1]
  thisPredCode<-gsub('juv|ad',"",thisPredVar)
  thisPredAge<-gsub(thisPredCode,"",thisPredVar)
  thisPredName<-str_trim(groupsDF$Name[groupsDF$Code==thisPredCode],side="both")
  thisNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisPredCode]
  thisVar<-paste(thisPredName,"_N",sep=""); predBiomassIC<-ncvar_get(ThisIC.nc,thisVar)
  thisRowName<-paste(thisPredCode,thisPredAge,sep="")
  thisPreyAvail<-sampledAvails[grep(thisRowName,rownames(sampledAvails)),] #rows are predators, columns are prey
  nonZeroPreys<-colnames(thisPreyAvail)[as.double(thisPreyAvail)>0]; npreys<-length(nonZeroPreys) #get names on columns with nonzero entries
  
  thisE<-get_first_number(biolLines[grep(paste("E_",thisPredCode,sep=""),biolLines)])
  
  cat("PRED ",thisPredCode,"--")
  
  predBiomassIC<-getICBiomass(groupName = thisPredName,ThisIC.nc, volume)

  thisTotPrey<-array(NA,dim=c(10,nlayers,nboxes))
  firstCohort<-1; lastCohort<-1
  if(thisNumCohorts>1){
    #get age mature
    predAgeMature<-get_first_number(biolLines[grep(paste(thisPredCode,"_age_mat",sep=""),biolLines)]) 
    if(thisPredAge=="juv"){
      firstCohort<-1; lastCohort<-predAgeMature
    }else{
      firstCohort<-predAgeMature+1; lastCohort<-thisNumCohorts
    }
    #get KLP and KUP for predator
    thisVar<-paste("KLP_",thisPredCode,sep=""); thisLine<-biolLines[grep(thisVar,biolLines)]
    thisKLP<-get_first_number(thisLine)
    thisVar<-paste("KUP_",thisPredCode,sep=""); thisLine<-biolLines[grep(thisVar,biolLines)]
    thisKUP<-get_first_number(thisLine)
    predSN<-unlist(lapply(seq(1,thisNumCohorts),FUN=getSNRN,name=thisPredName,whichN="SN",ThisIC.nc=ThisIC.nc))
    lowerSNlimits<-thisKLP*predSN; upperSNlimits<-thisKUP*predSN
  }
    

    thisPropEat<-reqFoodByGroupCohort[groupsDF$Code==thisPredCode,firstCohort:lastCohort,,]/totPreyAvail[groupsDF$Code==thisPredCode,firstCohort:lastCohort,,]
    #########################################
    ##NOTES NOTES NOTES #################
    # pick up from here
    ##
    ## make sure track mortalities (in numbers for verts, biomass otherwise)
    ## and track assimilated, to compare with growth
    ## now do the eating
    for(b in 1:nboxes){
      for(l in 1:nlayers){
        if(predBiomassIC[l,b]>0){
          for(prey in 1:npreys){
            thisPreyVar<-nonZeroPreys[prey]
            thisPreyCode<-gsub('juv|ad',"",thisPreyVar)
            thisPreyAge<-gsub(thisPreyCode,"",thisPreyVar)
            thisPreyName<-str_trim(groupsDF$Name[groupsDF$Code==thisPreyCode],side="both")
            thisPreyBiomass<-getICBiomass(thisPreyName,ThisIC.nc,volume)
            if(thisPreyBiomass[l,b]>0){
              if(thisNumCohorts>1){
                  
                
               sumPreyBiomass<-rep(0,(lastCohort-firstCohort+1))
                for(c in firstCohort:lastCohort){
                  thisGapeL<-lowerSNlimits[c]; thisGapeU<-upperSNlimits[c]
                  for(prey in 1:npreys){
                    #loop through preys, for each calculate how much (mg N) within gape size in this cell, multiply by sampledAvail
                    thisPreyVar<-nonZeroPreys[prey]
                    # cat("prey ",thisPreyVar,",,")
                    thisPreyCode<-gsub('juv|ad',"",thisPreyVar)
                    thisPreyAge<-gsub(thisPreyCode,"",thisPreyVar)
                    thisPreyNumCohorts<-groupsDF$NumCohorts[groupsDF$Code==thisPreyCode]
                    thisPreyName<-str_trim(groupsDF$Name[groupsDF$Code==thisPreyCode],side="both")
                    thisPreyBiomass<-getPreyBiomassAvail(thisPreyCode,thisPreyName,thisPreyNumCohorts,thisPreyAge,ThisIC.nc,volume,thisLayer,thisBox,biolLines,thisGapeL,thisGapeU)
                    
                    sumPreyBiomass[c]<-sumPreyBiomass[c]+thisPreyBiomass
                  }
                }
                
            }
            
           
          }
        }
      }
    
    
  
  
}

#get proportion of avail food that will be eaten, based on desired growth
#kill prey, and store deaths. Store mg N assimilated
