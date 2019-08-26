#reads in ICs and calculates BH alpha and beta for recruit_flag = 10 (SSB version)
source(paste(DIR$'General functions',"get_first_number.R",sep=""))
source(paste(DIR$'General functions',"myRounding.R",sep=""))
source(paste(DIR$'General functions',"getSNRN.R",sep=""))
source(paste(DIR$'General functions',"nonZeroMean.R",sep=""))
source(paste(DIR$'General functions',"getICNumbers.R",sep=""))

this_run<-"Base"

this_h<-0.7

calc_alpha<-function(R0,h){
  alpha<-(4*h*R0)/(5*h-1)
  return(alpha)
}
calc_beta<-function(B0,h){
  beta<-(B0*(1-h))/(5*h-1)
  return(beta)
}

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
plotPath<-paste(this_path,"Figures\\",this_run,"\\",sep="")

##we'll need the groups df 
groupsDF<-read.csv(paste(this_path,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
thisInitialConditionsFile<-paste(this_path,"CRAM_input_short.nc",sep="")
ThisIC.nc<-nc_open(thisInitialConditionsFile)
thisBiolFile<-paste(this_path,"CRAM_BH_hybrid_biol.prm",sep="")
biolLines<-readLines(thisBiolFile)

#get volume from ICs
volume<-ncvar_get(ThisIC.nc,"volume")[,,1]
nlayers<-dim(volume)[1]; nboxes<-dim(volume)[2]

for(g in 1:ng){
  thisCode<-groupsDF$Code[g]; thisName<-str_trim(groupsDF$Name[g],side="both"); thisNumCohorts<-groupsDF$NumCohorts[g]
  if(thisNumCohorts>1){
    thisVar<-paste("flagrecruit",thisCode,sep="")
    thisFlagRec<-get_first_number(biolLines[grep(thisVar,biolLines)])
    if(thisFlagRec==10){
      nums<-rep(0,thisNumCohorts); weights<-rep(0,thisNumCohorts)
      for(c in 1:thisNumCohorts){
        nums[c]<-getICNumbers(name=thisName,cohort=c,ThisIC.nc = ThisIC.nc)
        sn<-getSNRN(name=thisName,cohort=c, whichN = "SN",ThisIC.nc = ThisIC.nc)
        rn<-getSNRN(name=thisName,cohort=c, whichN = "RN", ThisIC.nc = ThisIC.nc)
        weights[c]<-sn+rn
      }
      thisSSB<-sum(weights*nums)
      ageClassSize<-get_first_number(biolLines[grep(paste(thisCode,"_AgeClassSize",sep=""),biolLines)])
      thisR0<-nums[1]/ageClassSize
      thisAlpha<-calc_alpha(R0=thisR0,h=this_h)
      thisBeta<-calc_beta(B0=thisSSB,h=this_h)
      if(thisAlpha>0 & thisBeta>0){
        thisVar<-paste("BHalpha_",thisCode,sep="")
        thisLine<-grep(thisVar,biolLines)
        newLine<-paste(thisVar,trunc(thisAlpha),collapse="\t")
        biolLines[thisLine]<-newLine
        # beta
        thisVar<-paste("BHbeta_",thisCode,sep="")
        thisLine<-grep(thisVar,biolLines)
        newLine<-paste(thisVar,trunc(thisBeta),collapse="\t")
        biolLines[thisLine]<-newLine
       
      }
    }
  }
}
writeLines(biolLines,thisBiolFile)
