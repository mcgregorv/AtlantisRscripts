getTotalNumbers<-function(Name,Cohort,startBox,stopBox){
  #assumes has availabilt init which is the biol initial conditions file lines
  thisVar<-paste("^",Name,Cohort,"_Nums",sep="")
  x<-grep(thisVar,init)
  thisLine<-init[(x+startBox):(x+stopBox)]
  xx<-as.double(get_first_number(thisLine,n="all"))
  totalNums<-sum(xx,na.rm=TRUE)
  return(totalNums)
}