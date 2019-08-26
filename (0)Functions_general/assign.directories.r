assign.directories<-function(base="")
	  {
	  DIR<-list()
	  
	  DIR[["Base"]]<-base
	  
	  DIR[["R"]]<-make.filename("R",DIR[["Base"]],T)

	  DIR[["General functions"]]<-make.filename("(0)Functions_general",DIR[["R"]],T)
    
	  DIR[["Reports"]]<-make.filename("reports",DIR[["Base"]],T)
	
	  DIR[["Data"]]<-make.filename("data",DIR[["Base"]],T)

	  DIR[["Figures"]]<-make.filename("figures",DIR[["Base"]],T)

	  DIR[["Tables"]]<-make.filename("tables",DIR[["Base"]],T)
    
    DIR[["CASAL"]]<-make.filename("CASAL",DIR[["Base"]],T)

	  return(DIR)
}