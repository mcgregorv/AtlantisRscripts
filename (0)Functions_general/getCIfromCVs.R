getCIfromCVs<-function(indices,cvs){
  ##example inputs
  # indices=WCSI.OBS.COMB
  # cvs=rep(0.3, length(indices))
  sigma.log.space <- sqrt( log(1 + cvs^2))
  mu.log.space <-  log(indices)-0.5*(sigma.log.space^2)

  lower.CI<- exp(mu.log.space-1.96*sigma.log.space)    
  upper.CI  <-exp(mu.log.space + 1.96*sigma.log.space)
  
  out<-NULL
  out$'LowerCI'<-lower.CI
  out$'UpperCI'<-upper.CI
  return(out)
}