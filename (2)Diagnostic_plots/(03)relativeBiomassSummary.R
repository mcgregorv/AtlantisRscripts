#set this run and path
this_run<-"base"
# this_out<-""

# this_out<-"Sensitivity\\KDENR_XXXE1"
this_out<-""
# # 

# outPath<-paste(DIR$'Base',"\\ATLANTIS\\Model_runs\\TBGB\\",this_run,"\\output",this_out,"\\",sep="")

# outPath<-paste(DIR$'Base',"\\ATLANTIS\\Model_runs\\TBGB\\",this_run,"\\Sensitivity\\KDENR\\outputKDENR_XXXD1\\",sep="")

this_path<-paste(DIR$'Base',"ATLANTISmodels\\",this_run,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")

plotPath<-paste(this_path,"..\\Figures\\",this_run,"\\",this_out,"",sep="")

daysTimeStep<-365
numStepsPerYear<-365/daysTimeStep
year0<-1880
fishingStartYear<-1900
modelStartYear<-1880

#read in groups
groups<-read.csv(paste(this_path,"..\\CRAM_Groups.csv",sep=""))

vertIndex<-groups$NumCohorts>=2
vertCodes<-as.character(groups$Code[vertIndex])

# ageInvertIndex<-groups$NumCohorts==2
# ageInvertCodes<-as.character(groups$Code[ageInvertIndex])

invertIndex<-groups$NumCohorts==1
invertCodes<-as.character(groups$Code[invertIndex])

# contants
nsecs <- 86400
ndays <- 365
burnin<-20 #years
startyear<-year0
g_per_ton <- 1e6
tons <- function(mgN) return(mgN * 5.7 * 20 / 1e9)

cat("### ------------ Reading in data                                         ------------ ###\n")
ncout<-"output"
nc_out <- nc_open(paste(outPath, ncout, ".nc", sep = ""))

biomass <- read.table(paste(outPath, ncout, "BiomIndx.txt", sep = ""), header = T)
tot_bio <- biomass[,c(1:(grep("Rel",colnames(biomass))[1]-1))]
colnames(tot_bio) <- c("Time", as.character(groups$Code), "DIN")
tot_bio$Time <- startyear-burnin + tot_bio$Time/365
tot_bio_l <- tot_bio %>% gather("Code", "value", -Time)

#get the vertebrates
pdf(paste(plotPath,"relativeBiomassVerts.pdf",sep=""),height=5)
plot(x=1,y=1,xlab="",ylab="Relative biomass",ylim=c(0,2),type="n",xlim=c(0,length(vertCodes)+1),xaxt="n")
abline(h=c(0.5,1.5),col="red",lty=2,lwd=2.5)
abline(h=c(1),col="red",lty=1,lwd=2.5)
par(las=2)
axis(at=seq(1,length(vertCodes)),labels=vertCodes,side=1)
for(v in 1:(length(vertCodes))){
  thisVert<-vertCodes[v]
  thisData<-tot_bio_l[tot_bio_l$Code==thisVert & tot_bio_l$Time>=startyear,]
  this1<-thisData$value[thisData$Time==startyear]
  this_rel_bio<-thisData$value/this1
  points(x=rep(v,length(this_rel_bio)),y=this_rel_bio,col=myBlue_trans,pch=16,cex=1.8)
  
}
dev.off()

#get the age-structured inverts
pdf(paste(plotPath,"relativeBiomassAgeInverts.pdf",sep=""),height=5)
plot(x=1,y=1,xlab="",ylab="Relative biomass",ylim=c(0,2),type="n",xlim=c(0,length(ageInvertCodes)+1),xaxt="n")
abline(h=c(0.5,1.5),col="red",lty=2,lwd=2.5)
abline(h=c(1),col="red",lty=1,lwd=2.5)
par(las=2)
axis(at=seq(1,length(ageInvertCodes)),labels=ageInvertCodes,side=1)
for(v in 1:(length(ageInvertCodes))){
  thisVert<-ageInvertCodes[v]
  thisData<-tot_bio_l[tot_bio_l$Code==thisVert & tot_bio_l$Time>=startyear,]
  this1<-thisData$value[thisData$Time==startyear]
  this_rel_bio<-thisData$value/this1
  points(x=rep(v,length(this_rel_bio)),y=this_rel_bio,col=myBlue_trans,pch=16,cex=1.8)
}
dev.off()

#get the inverts (not age-structured)
pdf(paste(plotPath,"relativeBiomassInverts.pdf",sep=""),height=5)
plot(x=1,y=1,xlab="",ylab="Relative biomass",ylim=c(0,5),type="n",xlim=c(0,length(invertCodes)+1),xaxt="n")
abline(h=c(0.5,1.5),col="red",lty=2,lwd=2.5)
abline(h=c(1),col="red",lty=1,lwd=2.5)
par(las=2)
axis(at=seq(1,length(invertCodes)),labels=invertCodes,side=1)
for(v in 1:(length(invertCodes))){
  thisVert<-invertCodes[v]
  thisData<-tot_bio_l[tot_bio_l$Code==thisVert & tot_bio_l$Time>=startyear,]
  this1<-thisData$value[thisData$Time==startyear]
  this_rel_bio<-thisData$value/this1
  points(x=rep(v,length(this_rel_bio)),y=this_rel_bio,col=myBlue_trans,pch=16,cex=1.8)
}
dev.off()

