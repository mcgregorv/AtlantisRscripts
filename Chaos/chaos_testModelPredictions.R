thisModelDesc<-"BP_glm"
thisModelDesc <-"ALL_glm"


load(file=paste(DIR$'Data', "Chaos\\", thisModelDesc,"Model", sep="")) # brings in this_model

plotPath <- paste(DIR$'Figures',"Chaos\\BP_glmModel", sep="")
groupsDF<- read.csv(paste(basePath,"..\\CRAM_groups.csv", sep="")); ng<-dim(groupsDF)[1]
