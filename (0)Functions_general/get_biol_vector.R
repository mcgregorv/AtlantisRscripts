get_biol_vector<-function(thisLine,vector_length){
  #used to get the value of a variable in input file, where there is an example value (or any other value) in the same line which we need to ignore
  yy<-gsub("([^\\d*\\.?\\d*$])","#",thisLine,perl=TRUE)
  yyy<-unlist(str_split(yy,"#"))
  xPos<-grep("[^\\d]",yyy)[1:vector_length]
  thisVec<-as.numeric(yyy[xPos])
  return(thisVec)
}