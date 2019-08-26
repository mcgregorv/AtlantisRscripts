findMyPolygon<-function(x,y,boxes,nboxes){
  out<-NULL
  for(b in 1:nboxes){
    thisPolyY<-boxes[[b]]$y
    thisPolyX<-boxes[[b]]$x
    test<-point.in.polygon(point.x=x,point.y=y,pol.x=thisPolyX,pol.y=thisPolyY)
#     cat(test,"--")
    if(test>=1){
      out<-c(out,b)
    }    
  }
  return(out)
}