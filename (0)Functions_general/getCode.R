getCode<-function(x){
  #gets species code from pPREY2XXX1 or pPREYXXX
  yy<-gsub("([\\d*\\.?\\d*$])","#",x,perl=TRUE)
  yyy<-unlist(str_split(yy,"#"))
  if(length(yyy)==1){
    code<-substr(x,start=6,stop=nchar(x))
  }else{
    code<-yyy[2]
  }
  return(code)
}