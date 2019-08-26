# read in tracers, and calculate proprotion of spatial overlap (by biomass)
#bring in the run, the observed data, the catch history and plot together
nboxes<-30

mg_2_tonne<-2e-8; X_CN<-5.7

catchYears<-seq(1900,2014); ny<-length(catchYears)

#################################### TRACERS
#plot all tracers for a given box and layer
this_run<-"Base"; this_out <- "Base"
get_age_mat <- function(x){
  #x is species group code
  thisVar<-paste(x,"_age_mat", sep="")
  y<-grep(thisVar,biolLines); z<-get_first_number(biolLines[y])
  return(z)
}
biolLines <- readLines(paste(this_path, "CRAM_BH_hybrid_biol.prm", sep=""))

burnin<-35 #number of years to skip in plot

this_path<-paste(DIR$'Base',"ATLANTISmodels\\ArchivedModels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")

ThisNC.nc <- nc_open(paste(outPath, "output.nc", sep=""))
thisVol <- ncvar_get(ThisNC.nc, "volume")
nlayers <- dim(thisVol)[1]; nts <- dim(thisVol)[3]; nboxes <- dim(thisVol)[2]

groupsDF<- read.csv(paste(this_path,"..\\CRAM_Groups.csv", sep="")); ng<-dim(groupsDF)[1]
asIndex <- groupsDF$NumCohorts>1 # for age-structured groups
asGroupsDF <- groupsDF[asIndex,]; nag <- dim(asGroupsDF)[1]

getBiomass <- function(Name, cohort){
  thisNumTracer <- paste(Name, cohort,"_Nums", sep=""); thisRNTracer <- paste(Name, cohort,"_ResN", sep="");
  thisSNtracer <- paste(Name, cohort,"_StructN", sep="");
  thisNums <- ncvar_get(ThisNC.nc, thisNumTracer); thisResN <- ncvar_get(ThisNC.nc, thisRNTracer)
  thisStructN <- ncvar_get(ThisNC.nc, thisSNtracer)
  thisBiomass <- (thisNums * (thisResN  + thisStructN)) * mg_2_tonne * X_CN
  return(thisBiomass)
}

storeBiomass <- array(NA, dim=c(nag, 2, nboxes))
store2DSpatialProportion <- 0*storeBiomass

for(g in 1:nag){
  thisName <- str_trim(asGroupsDF$Name[g], side="both"); thisCode <- as.character(asGroupsDF$Code[g])
  thisNumCohorts <- asGroupsDF$NumCohorts[g]
  this_age_mat <- get_age_mat(thisCode)
  juvACs <- 1:this_age_mat; adACs <- (this_age_mat + 1):thisNumCohorts # Age classes that are juvenile or adult
  juvBiomass <- array(NA, dim=c(length(juvACs), dim(thisVol))); adBiomass <- array(NA, dim=c(length(adACs), dim(thisVol)))
  for(a in 1:length(juvACs)){
    thisCohort <- juvACs[a]
    juvBiomass[a,,,] <- getBiomass(Name=thisName, cohort=thisCohort)
  }
  for(a in 1:(length(adACs))){
    thisCohort <- adACs[a]
    adBiomass[a,,,] <- getBiomass(Name=thisName, cohort=thisCohort)
  }
  # sum over depth layers and age classes, so have just by time and box, then get median over time, so have value by box
  medJuvBiomass <-apply( apply(juvBiomass, c(3,4), sum, na.rm=TRUE),1, median, na.rm=TRUE) 
  medAdBiomass <-apply( apply(adBiomass, c(3,4), sum, na.rm=TRUE),1, median, na.rm=TRUE) 
  storeBiomass[g,1,]<- medJuvBiomass; storeBiomass[g,2,] <- medAdBiomass
  thisJuvTotal <- sum(medJuvBiomass); store2DSpatialProportion[g,1,]<- medJuvBiomass/thisJuvTotal
  thisAdTotal <- sum(medAdBiomass); store2DSpatialProportion[g,2,]<- medAdBiomass/thisAdTotal
  # replaced with age-structure version
  # thisTracer <- paste(thisName, "_N", sep="")
  # thisData <- ncvar_get(ThisNC.nc, thisTracer)
  # if(length(dim(thisData))==3){
  #   thisBiom <- apply(thisData * thisVol, c(2,3), sum)*mg_2_tonne * X_CN
  # }else{
  #   thisBiom <- (thisData * thisVol[nlayers,,])*mg_2_tonne * X_CN
  # }
  # # take out the burnin, then take the median values
  # thisMedBiom <- apply(thisBiom[,burnin:nts],1,median, na.rm=TRUE)
  # storeBiomass[g,]<- thisMedBiom
  # thisTotal <- sum(thisMedBiom); store2DSpatialProportion[g,]<- thisMedBiom/thisTotal
}

save(list=c("store2DSpatialProportion", "storeBiomass", "asGroupsDF"), file=paste(this_path,"store2DSpatialProportion",sep=""))

