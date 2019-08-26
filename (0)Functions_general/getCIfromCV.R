getCIfromCV<-function(indices,cv){
  ##example inputs
  # indices=WCSI.OBS.COMB
  # cv=0.3
  sigma.log.space <- rep(sqrt( log(1 + cv^2)), length(indices))
  mu.log.space <-  log(indices)-0.5*(sigma.log.space^2)
  
  lower.CI<- exp(mu.log.space-1.96*sigma.log.space)    
  upper.CI  <-exp(mu.log.space + 1.96*sigma.log.space)
  
  out<-NULL
  out$'LowerCI'<-lower.CI
  out$'UpperCI'<-upper.CI
  return(out)
}
