thisFile<-paste(DIR$'Base',"reports\\(01)BaseReport\\targetAudience.csv",sep="")
ta<-read.csv(thisFile,sep=",")
sort(table(ta$'Journals.published.in'),decreasing=TRUE)[1:5]
# Fish and Fisheries                     ICES Journal of Marine Science 
# 7                                                  7 
# PLoS ONE Canadian Journal of Fisheries and Aquatic Sciences 
# 5                                                  4 
# Ecological Modelling 
# 3 