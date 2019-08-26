#read in tracers and check out temperature to see if correct ROMS being used
this_run<-"base"
this_out<-"FISH4"

modelYears<-1865:2015

basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
groupsDF<-read.csv(paste(basePath,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]


ThisNC.nc<-nc_open(paste(basePath,"base\\output",this_out,"\\output.nc", sep=""))

thisVol<-ncvar_get(ThisNC.nc, "volume")
nlayers<-dim(thisVol)[1]
baseVol<-ncvar_get(BaseNC.nc, "volume"); nts<-dim(baseVol)[3]

thisDiatom<-ncvar_get(ThisNC.nc, "Diatom_N")
baseDiatom<-ncvar_get(BaseNC.nc, "Diatom_N")

dynBoxes<-2:25

thisDiatomBiomass<-apply(thisDiatom[-nlayers,dynBoxes,] * thisVol[-nlayers, dynBoxes,], 3, sum) * mg_2_tonne * X_CN
baseDiatomBiomass<-apply(baseDiatom[-nlayers, dynBoxes, ] * baseVol[-nlayers, dynBoxes, ], 3, sum) * mg_2_tonne * X_CN

thisMeanDiatom<-apply(thisDiatom[-nlayers, dynBoxes,], 3, nonZeroMean)
baseMeanDiatom<-apply(baseDiatom[-nlayers, dynBoxes,], 3, nonZeroMean)

allTracers<-names(BaseNC.nc$var)

timeIndex<-1975-1865+1
basicInput<-read.csv(paste("C:\\Projects\\2018\\CREWE\\Model\\Base\\CREWEBasicInput_MON.csv", sep=""))

area_m2<-sum(baseVol[nlayers, dynBoxes,timeIndex])
area_km2<-area_m2 * 1e-6

storeBiomassPerHabKm2<-rep(NA, ng)
basicInput$CheckBiom<-NA
biomYears<-1865:2015

storeBiomassTS<-array(NA, dim=c(dim(basicInput)[1], nts))

g=grep("BC", groupsDF$Code)

for(g in 1:ng){
  thisCode<-groupsDF$Code[g]
  thisName<-str_trim(groupsDF$Name[g], side="both")
  thisRow<-basicInput$Code==thisCode
  thisHabProp<-basicInput[thisRow,"HabArea"]
  thisTracer<-paste(thisName,"_N", sep="")
  orh<-ncvar_get(BaseNC.nc,thisTracer)
  if(length(dim(orh))==3){
    orh_biomass<-apply(orh * baseVol, 3, sum) * mg_2_tonne * X_CN
  } else{
    orh_biomass <- apply(orh * baseVol[nlayers,,], 2, sum) * mg_2_tonne * X_CN
  }
  orh_b0_km2<-orh_biomass[timeIndex]/area_km2
  storeBiomassPerHabKm2[g]<-orh_b0_km2 * (1/thisHabProp)
  basicInput$CheckBiom[thisRow]<-orh_b0_km2 * (1/thisHabProp)
  storeBiomassTS[thisRow,]<-orh_biomass
}

write.csv(basicInput,file=paste("C:\\Projects\\2018\\CREWE\\Model\\Base\\BasicinputCheck.csv", sep=""), row.names = FALSE)

g=grep("ZM", groupsDF$Code)

## check catches
ALLcatchHistories<-read.csv(paste(DIR$'Tables',"ALLcatchHistories.csv",sep=""))

storeFs<-0*storeBiomassTS
for(g in 1:ng){
  thisCode<-as.character(groupsDF$Code[g]); 
  if(groupsDF$IsFished[g]==1){
    thisCatches<-ALLcatchHistories[,thisCode]*1e-3 # was kg, now tonnes
    thisBiomasses<-storeBiomassTS[g,]
    thisF<-thisCatches/thisBiomasses[modelYears %in% ALLcatchHistories[,1]]
    storeFs[g,modelYears %in% ALLcatchHistories[,1]]<-thisF
  }
}
save(list=c("storeFs"), file=paste(basePath, "\\base\\EWEbase\\CRAM_Fs", sep=""))

load(paste(basePath, "\\base\\EWEbase\\CRAM_Fs", sep=""))
# load(paste(basePath, "\\base\\EWEbase\\CRAM_biomasses", sep=""))

biolLines<-readLines(paste(basePath,"CRAM_BH_hybrid_biol.prm", sep=""))
groupsDF<-read.csv(paste(basePath,"CRAM_groups.csv", sep="")); ng<-dim(groupsDF)[1]
recruitWeights<-rep(NA, ng)
for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisVar<-paste("KWSR_",thisCode,sep=""); x<-grep(thisVar, biolLines)
  if(length(x)>0){
    thisKWSR<-get_first_number(biolLines[x])
  }
  thisVar<-paste("KWRR_",thisCode,sep=""); x<-grep(thisVar, biolLines)
  if(length(x)>0){
    thisKWRR<-get_first_number(biolLines[x])
  }
  recruitWeights[g]<-(thisKWRR + thisKWSR) * mg_2_tonne * X_CN
  
}
index<-groupsDF$NumCohorts>1
temp<-cbind(as.character(groupsDF$Code[index]), recruitWeights[index])

write.csv(temp, file=paste(DIR$'Tables',"WeightOfRecruits_tonnes.csv", sep=""), row.names = FALSE)



# matchInputNames<-gsub("\\(|)|'|-","",basicInput[,"Group.name"])
# groupsDF$LongName<-gsub("Invertibrate","Invertebrate",groupsDF$LongName)
# groupsDF$LongName[groupsDF$Code=="PIN"]<-"NZ fur seal"
# for(g in 1:ng){
#   thisCode<-groupsDF$Code[g]
#   if(groupsDF$NumCohorts[g]>1){
#     thisName<-str_trim(groupsDF$Name[g], side="both")
#     thisLongName<-str_trim(gsub("\\(|)|'|-|commercial","",groupsDF$LongName[g]),side="both")
#     thisRow<-grep(gsub("_"," ",thisName),matchInputNames)
#     if(length(thisRow)==0){
#       thisRow<-grep(thisLongName, matchInputNames)
#       if(length(thisRow)==0){
#         cat("No row for ",thisName,"\n")
#       }
#     }
#     if(length(thisRow)>0){
#       thisHabProp<-basicInput[thisRow,"Habitat.area..fraction."]
#       thisTracer<-paste(thisName,"_N", sep="")
#       orh<-ncvar_get(BaseNC.nc,thisTracer)
#       if(length(dim(orh))==3){
#         orh_biomass<-apply(orh * baseVol, 3, sum) * mg_2_tonne * X_CN
#       } else{
#         orh_biomass <- apply(orh * baseVol[nlayers,,], 2, sum) * mg_2_tonne * X_CN
#       }
#       
#       orh_b0_km2<-orh_biomass[timeIndex]/area_km2
#       storeBiomassPerHabKm2[g]<-orh_b0_km2 * (1/thisHabProp)
#       basicInput$CheckBiom[thisRow]<-orh_b0_km2 * (1/thisHabProp)
#     }
#   }
# }
