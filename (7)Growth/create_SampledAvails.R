## take 2 - take 1 had everything avail to everythin to eat!

source(paste(DIR$'General functions',"get_first_number.R",sep=""))
source(paste(DIR$'General functions',"myRounding.R",sep=""))
source(paste(DIR$'General functions',"getSNRN.R",sep=""))
source(paste(DIR$'General functions',"nonZeroMean.R",sep=""))
source(paste(DIR$'General functions',"getICNumbers.R",sep=""))
source(paste(DIR$'General functions',"getRequiredFoodPerDay_bySpace.R",sep=""))

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

#get volume from ICs
volume<-ncvar_get(ThisIC.nc,"volume")[,,1]
nlayers<-dim(volume)[1]; nboxes<-dim(volume)[2]

ppGroups<-as.character(read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ppGroups.csv",sep=""))[,1]);
npg<-length(ppGroups)

load(paste(DIR$"Data","\\eating\\interaction_spatial_array",sep="")); ##to bring in interaction_spatial_array. dim = npg, npg, nlayers, nboxes
load(paste(DIR$"Data","\\eating\\interaction_gapeSize_array",sep="")); ## dim = ng, 10, ng, 10
load(paste(DIR$"Data","\\eating\\interaction_trophicLevel_array",sep="")); ##dim = npg, npg

#gape size only used at this stage for age-structured
ageCohortLinking<-read.csv(paste(DIR$'Tables',"SampledAvailabilities\\ageCohortLinking.csv",sep=""))

indexTo1<-function(x){
  y<-grep("ad",x,invert=TRUE)
  if(length(y)==0){y<-0}
  return(y)
}
temp<-array(eval(interaction_gapeSize_array),dim=c(ng,10,ng,10))
interaction_gapeSize_noNA<-apply(temp,c(1,2,3,4),NAtozero)
sampledAvails<-array(NA,dim=c(npg,npg))

juvCohortLink<-apply(ageCohortLinking[,-1],c(1,2),FUN=indexTo1)
adCohortLink<-apply(juvCohortLink,c(1,2),FUN=function(x){(x-1) %% 2})

for(g in 1:npg){
  thisVar<-ppGroups[g]
  thisCode<-gsub("ad|juv","",thisVar)
  thisAge<-gsub(thisCode,"",thisVar)
  temp<-ageCohortLinking[ageCohortLinking[,1]==thisCode,-1]
  thisCohorts<-seq(1,10)[temp==thisVar]
  thisCohorts<-thisCohorts[!is.na(thisCohorts)]
  # TROPHIC LEVEL
  thisTLint<-interaction_trophicLevel_array[g,]
  thisTL_ad<-thisTLint[grep("ad",ppGroups)]; thisTL_juv<-thisTLint[grep("ad",ppGroups,invert=TRUE)]
  ## GAPE SIZE
  thisGSint<-interaction_gapeSize_noNA[groupsDF$Code==thisCode,thisCohorts,,]
  if(length(dim(thisGSint))>2){
    thisGSint_mean<-apply(thisGSint,c(2,3),mean,na.rm=TRUE) # if have more than one cohort for thisVar, need to average over them
  } else{
    thisGSint_mean<-thisGSint 
  }
  thisGSint_juv<-apply(thisGSint_mean * juvCohortLink,1,mean,na.rm=TRUE)
  thisGSint_ad<-apply(thisGSint_mean * adCohortLink,1,mean,na.rm=TRUE)[groupsDF$NumCohorts>1]
  ## non age-structured need to have 1's in gape size
  if(thisAge==""){
    thisGSint_juv<-rep(1,length(thisGSint_juv))
  }
  ## SPATIAL INTERACTION
  thisSpatint<-interaction_spatial_array[g,,,]
  thisSpatialJuv<-apply(thisSpatint[grep("ad",ppGroups,invert=TRUE),,],1,mean,na.rm=TRUE)
  thisSpatialAd<-apply(thisSpatint[grep("ad",ppGroups),,],1,mean,na.rm=TRUE)
  ## combine them
  thisAdInt<-thisTL_ad*thisGSint_ad*thisSpatialAd
  # will have to fill the adult version in with zeros for non-age-structured prey
  # filledAdInt<-rep(0,ng); filledAdInt[groupsDF$NumCohorts>1]<-thisAdInt ## don;t actualy need this
  thisJuvInt<-thisTL_juv*thisGSint_juv*thisSpatialJuv
  # put them into sampledAvails
  sampledAvails[g,grep("ad",ppGroups,invert=TRUE)]<-thisJuvInt
  sampledAvails[g,grep("ad",ppGroups)]<-thisAdInt
}
write.csv(sampledAvails,paste(DIR$'Tables',"SampledAvailabilities\\ExpectedValues.csv",sep=""),row.names = FALSE)





