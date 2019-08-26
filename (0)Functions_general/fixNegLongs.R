fixNegLongs<-function(x){
  if(!is.na(x)){  
    if(x<0){
      diff<-abs(-180-x)
      newx<-180+diff
    } else{
      newx<-x
    }
  } else {
    newx<-x
  }
  return(newx)
}