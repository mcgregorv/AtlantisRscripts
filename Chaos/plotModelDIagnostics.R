thisModelDesc<-"BP_glm"
thisModelDesc <-"ALL_glm"

load(file=paste(DIR$'Data', "Chaos\\", thisModelDesc,"Model", sep="")) # brings in this_model

plotPath <- paste(DIR$'Figures',"Chaos\\BP_glmModel", sep="")

pdf(paste(plotPath, "diagnostics.pdf"))
par(mfrow=c(2,2))
plot(this_model)
dev.off()
