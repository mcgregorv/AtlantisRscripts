myRounding<-function(x,fraction=0.25,direction="NULL"){ 
  #the whole bit 
  y<-trunc(x) 
  z<-x-y 
  if(direction=="down"){ 
    q<-trunc(z*(1/fraction))*fraction 
  } else if(direction=="up"){ 
    q<-ceiling(z*(1/fraction))*fraction 
  } else{ 
    q<-round(z*(1/fraction))*fraction 
  } 
  return(q+y) 
}
