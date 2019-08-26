## 2 versions of this in Atlantis - the original one, adopted from PPBIM, and the one that more closely matches the standard one in the literature
# prey_avail is mg N per m^3 and assumes spatial, temporal, gape, habitat, etc all taken into account
# C is in m^3 per seccond - the input value is m^3 per day
# mum is mg N per second - input values are mg N per day
# ht is in seconds - input is in days (although actually takes in hta, htb)
numsec <- 24 * 60 * 60
plotPath <- paste(DIR$'Figures',"SETAs_", sep="")
mod_col<-myGrey; std_col<-"midnightblue"; cap_col<-"red"

predRes <- function(C, prey_avail, efficiencies, prey_types, ht, mum, predcase){
  # example parameters to show structure. Make sure prey_types is the same length as prey_avail, and that efficiencies is at least as long as max(prey_types)
  # C=0.1; mum=5e-4; ht=100; predcase=1
  # prey_avail<-runif(5); prey_types <- c(1,1,1,2,1); efficiencies <- c(0.5,0.2,0)
  total_prey <- sum(prey_avail)
  total_absorbed_prey <- 0
  nprey <- length(prey_avail)
  for(p in 1:nprey){
    thisEffIndex <- prey_types[p]
    thisPreyAvail <- prey_avail[p]
    thisE <- efficiencies[thisEffIndex]
    total_absorbed_prey <- total_absorbed_prey + (thisE * thisPreyAvail)
  }
  if(predcase==1){
    attackRate <- C/(1 + (C/mum)*total_absorbed_prey)
  } else{
    attackRate <- C/(1 + (C * ht * total_prey))
  }
}
X_CN = 5.7
############### do some plots to help understand the effects of the differences ##############
# first plot shows example with E=1 for all types, handling time = 1/mum, - should just see mum limit kick in at different place
# ht = 500 # seconds
test_prey <- seq(0,5,length.out = 1000)
Corig <- 10; C_secs <- Corig/numsec
orig_mum <- c(1, 1.5,2,5,10); test_mum <- orig_mum/numsec; test_ht <- 1/test_mum
prey_types <- rep(1, length(test_prey)); efficiencies <- 1
colByMum <- colorRampPalette(colors=c(myGold, myGreen, myAqua, myBlue))(length(test_mum))
E1 <- efficiencies[prey_types==1][1]
# m=length(test_mum)
# this_mum <- test_mum[m]; this_ht <- 1/this_mum
# test_predcase1 <- unlist(lapply(test_prey, FUN=function(x){predRes(C=C_secs, prey_avail = x, efficiencies = efficiencies, prey_types = prey_types, ht=this_ht, mum=this_mum, predcase=1)}))
# thisCap <- this_mum/E1
# test1_index <- test_predcase1>thisCap; test_predcase1[test1_index]<- thisCap
# test_predcase8 <- unlist(lapply(test_prey, FUN=function(x){predRes(C=C_secs, prey_avail = x, efficiencies = efficiencies, prey_types = prey_types, ht=this_ht, mum=this_mum, predcase=8)}))
# test_preyPrey8 <- test_prey * test_predcase8
# test8_index <- (test_preyPrey8)>thisCap; test_preyPrey8[test8_index]<- thisCap
# 
# plot(x=test_prey, y=test_prey * test_predcase1, type="n")
# for(m in 1:(length(test_mum))){
#   this_mum <- test_mum[m]; this_ht <- 1/this_mum
#   thisCap <- this_mum/E1
#   test_predcase1 <- unlist(lapply(test_prey, FUN=function(x){predRes(C=C_secs, prey_avail = x, efficiencies = efficiencies, prey_types = prey_types, ht=this_ht, mum=this_mum, predcase=1)}))
#   test1_index <- test_predcase1>thisCap; test_predcase1[test1_index]<- thisCap
#   points(x=test_prey, y=test_prey * test_predcase1, type="l", lwd=1.5, col=colByMum[m])
# }
# 
# plot(x=test_prey, y=test_preyPrey8, type="n")
# for(m in 1:(length(test_mum))){
#   this_mum <- test_mum[m]; this_ht <- 1/this_mum
#   thisCap <- this_mum/E1
#   test1_index <- test_predcase1>thisCap; test_predcase1[test1_index]<- thisCap
#   test_predcase8 <- unlist(lapply(test_prey, FUN=function(x){predRes(C=C_secs, prey_avail = x, efficiencies = efficiencies, prey_types = prey_types, ht=this_ht, mum=this_mum, predcase=8)}))
#   test_preyPrey8 <- test_prey * test_predcase8
#   test8_index <- (test_preyPrey8)>thisCap; test_preyPrey8[test8_index]<- thisCap
#   points(x=test_prey, y=test_preyPrey8, type="l", lwd=1.5, col=colByMum[m])
# }

m=3
this_mum <- test_mum[m]; this_ht <- 1/this_mum
thisCap <- this_mum/E1
test_predcase1 <- unlist(lapply(test_prey, FUN=function(x){predRes(C=C_secs, prey_avail = x, efficiencies = efficiencies, prey_types = prey_types, ht=this_ht, mum=this_mum, predcase=1)}))
test1_index <- test_predcase1>thisCap; test_predcase1[test1_index]<- thisCap
test_predcase8 <- unlist(lapply(test_prey, FUN=function(x){predRes(C=C_secs, prey_avail = x, efficiencies = efficiencies, prey_types = prey_types, ht=this_ht, mum=this_mum, predcase=8)}))

ht_scalars <- c(0.1,0.5, 0.8,1.2,1.5,2)
test_ht <- ht_scalars*this_ht
predcase8ByHt <- array(NA, dim=c(length(test_ht), length(test_prey))); predRate8ByHt <- predcase8ByHt
for(h in 1:(length(test_ht))){
  this_predcase <- unlist(lapply(test_prey, FUN=function(x){predRes(C=C_secs, prey_avail = x, efficiencies = efficiencies, prey_types = prey_types, ht=test_ht[h], mum=this_mum, predcase=8)}))
  test_preyPrey8 <- test_prey * this_predcase
  test8_index <- (test_preyPrey8)>thisCap; test_preyPrey8[test8_index]<- thisCap
  predcase8ByHt[h,]<-test_preyPrey8
  predRate8ByHt[h,] <- this_predcase
}
colByHt <- colorRampPalette(colors=c( myGreen, myBlue, myRed, myOrange))(length(test_ht))

png(paste(plotPath, "PreyAttackedRatesExample.png", sep=""),width=700, height=300)
par(mfrow=c(1,2), mar=c(4,4,1,1))
par(las=0, xpd=FALSE)
plot(x=test_prey, y=test_predcase1, col=mod_col, type="l", lwd=2.5, yaxt="n", xaxt="n", xlab="", ylab="", bty="n", ylim=c(0,max(test_predcase8)))
points(x=test_prey, y=test_predcase8, col=std_col, type="l", lwd=2, lty=2)
abline(h=thisCap, lty=3, col=cap_col, lwd=2)
mtext("Prey biomass attacked/available",side=2, line=1, cex=1.5)
mtext("Prey biomass available", side=1, line=1, cex=1.5)
par(xpd=NA)
arrows(x0=-0.2,y0=0, x1=-0.2, y1=max(test_predcase8), lwd=2)
arrows(x0=-0.2,y0=-0, x1=max(test_prey)*1.01, y1=-0, lwd=2)
legend(legend=c("Standard", "Modified"), col=c(std_col, mod_col), lwd=2, lty=c(2,1), x="right", bty="n", inset=0.1, seg.len=3, cex=1.2)

par(las=0, xpd=FALSE)
plot(x=test_prey, y=test_prey * test_predcase1, type="l", lwd=2, col=mod_col, ylim=c(0, thisCap), bty="n",
     yaxt="n", xaxt="n", ylab="", xlab="")
mtext("Prey biomass attacked",side=2, line=1, cex=1.5)
mtext("Prey biomass available", side=1, line=1, cex=1.5)
arrows(x0=0,y0=0, x1=0, y1=thisCap, lwd=2)
arrows(x0=0,y0=0, x1=max(test_prey)*1.01, y1=0, lwd=2)
test_preyPrey8 <- test_prey * test_predcase8
test8_index <- (test_preyPrey8)>thisCap; test_preyPrey8[test8_index]<- thisCap
points(x=test_prey, y=test_preyPrey8, type="l", lwd=2.5, col=std_col, lty=2)
# points(x=test_prey[test1_index], y=(test_prey * test_predcase1)[test1_index], lty=4, col="red", type="l")
abline(h=thisCap, col=cap_col, lty=3, lwd=2)
par(las=1)
axis(at=thisCap, label="Cap", side=2, col.axis=cap_col, cex=2, col.lab=myBlue, tick = FALSE, line=-0.5, cex.axis=1.2)
legend(legend=c("Standard", "Modified"), col=c(std_col, mod_col), lwd=2, lty=c(2,1), x="right", bty="n", inset=0.1, seg.len=3, cex=1.2)
dev.off()

png(paste(plotPath, "PreyAttackedRatesExample_incl_ht.png", sep=""),width=1000, height=300)
par(mfrow=c(1,3), mar=c(4,4,1,1), oma=c(0,0,1.2,0))
par(las=0, xpd=FALSE)
plot(x=test_prey, y=test_predcase1, col=mod_col, type="l", lwd=2.5, yaxt="n", xaxt="n", xlab="", ylab="", bty="n", ylim=c(0,max(test_predcase8)))
points(x=test_prey, y=test_predcase8, col=std_col, type="l", lwd=2, lty=2)
abline(h=thisCap, lty=3, col=cap_col, lwd=2)
mtext("A", side=3, adj=0, line=0.1, font=2, cex=1.2)
mtext("Prey biomass attacked/available",side=2, line=1, cex=1.5)
mtext("Prey biomass available", side=1, line=1, cex=1.5)
par(xpd=NA)
arrows(x0=-0.2,y0=0, x1=-0.2, y1=max(test_predcase8), lwd=2)
arrows(x0=-0.2,y0=-0, x1=max(test_prey)*1.01, y1=-0, lwd=2)
legend(legend=c("Standard", "Modified"), col=c(std_col, mod_col), lwd=2, lty=c(2,1), x="right", bty="n", inset=0.1, seg.len=3, cex=1.2)

par(las=0, xpd=FALSE)
plot(x=test_prey, y=test_prey * test_predcase1, type="l", lwd=2, col=mod_col, ylim=c(0, thisCap), bty="n",
     yaxt="n", xaxt="n", ylab="", xlab="")
mtext("B", side=3, adj=0, line=0.1, font=2, cex=1.2)
mtext("Prey biomass attacked",side=2, line=1, cex=1.5)
mtext("Prey biomass available", side=1, line=1, cex=1.5)
arrows(x0=0,y0=0, x1=0, y1=thisCap, lwd=2)
arrows(x0=0,y0=0, x1=max(test_prey)*1.01, y1=0, lwd=2)
test_preyPrey8 <- test_prey * test_predcase8
test8_index <- (test_preyPrey8)>thisCap; test_preyPrey8[test8_index]<- thisCap
points(x=test_prey, y=test_preyPrey8, type="l", lwd=2.5, col=std_col, lty=2)
abline(h=thisCap, col=cap_col, lty=3, lwd=2)
par(las=1)
axis(at=thisCap, label="Cap", side=2, col.axis=cap_col, cex=2, col.lab=myBlue, tick = FALSE, line=-0.5, cex.axis=1.2)
legend(legend=c("Standard", "Modified"), col=c(std_col, mod_col), lwd=2, lty=c(2,1), x="right", bty="n", inset=0.1, seg.len=3, cex=1.2)

par(las=0, xpd=FALSE)
test_preyPrey8 <- test_prey * test_predcase8
test8_index <- (test_preyPrey8)>thisCap; test_preyPrey8[test8_index]<- thisCap
par(las=0, xpd=FALSE)
plot(x=test_prey, y=test_preyPrey8, type="l", lwd=2, col=std_col, ylim=c(0, thisCap), bty="n", lty=2,
     yaxt="n", xaxt="n", ylab="", xlab="")
for(h in 1:length(test_ht)){
  points(x=test_prey, y=predcase8ByHt[h,], type="l", lwd=2, col=colByHt[h])
}
mtext("C", side=3, adj=0, line=0.1, font=2, cex=1.2)
mtext("Prey biomass attacked",side=2, line=1, cex=1.5)
mtext("Prey biomass available", side=1, line=1, cex=1.5)
arrows(x0=0,y0=0, x1=0, y1=thisCap, lwd=2)
arrows(x0=0,y0=0, x1=max(test_prey)*1.01, y1=0, lwd=2)
abline(h=thisCap, col=cap_col, lty=3, lwd=2)
par(las=1)
axis(at=thisCap, label="Cap", side=2, col.axis=cap_col, cex=2, col.lab=myBlue, tick = FALSE, line=-0.5, cex.axis=1.2)
legend(title="Handling time scalar", cex=1.2, legend=formatC(test_ht/this_ht,  zero.print = TRUE, digits=1, format="f"), col=colByHt, lwd=2, lty=1, x="bottomright", bty="n", inset=0.05, seg.len=3)

dev.off()

## chaning C..? change relative to mum (and hence ht for comparative ht=1/mum)
c_scalars <- c(0.2,0.5,1, 2, 5)
test_C <- c_scalars * this_mum*numsec; C_secs <- test_C/numsec; nCs <- length(test_C)
colByC <- colorRampPalette(colors=c(myGold, myAqua, "midnightblue"))(nCs)
predcase8ByC <- array(NA, dim=c(nCs, length(test_prey))); predRate8ByC <- predcase8ByHt
predcase1ByC <- predcase8ByC; predRate1ByC <- predRate8ByC
for(c in 1:nCs){
  this_predcase <- unlist(lapply(test_prey, FUN=function(x){predRes(C=C_secs[c], prey_avail = x, efficiencies = efficiencies, prey_types = prey_types, ht=this_ht, mum=this_mum, predcase=8)}))
  test_preyPrey8 <- test_prey * this_predcase
  test8_index <- (test_preyPrey8)>thisCap; test_preyPrey8[test8_index]<- thisCap
  predcase8ByC[c,]<-test_preyPrey8
  predRate8ByC[c,] <- this_predcase
  ## modified HII
  this_predcase1 <- unlist(lapply(test_prey, FUN=function(x){predRes(C=C_secs[c], prey_avail = x, efficiencies = efficiencies, prey_types = prey_types, ht=this_ht, mum=this_mum, predcase=1)}))
  test1_index <- this_predcase1>thisCap; this_predcase1[this_predcase1]<- thisCap # capped as rate
  predcase1ByC[c,]<-this_predcase1 * test_prey
  predRate1ByC[c,] <- this_predcase1
  
}

png(paste(plotPath, "PreyAttackedRatesExample_C_ht.png", sep=""),width=800, height=300)
par(mfrow=c(1,3), mar=c(4,4,1,1), oma=c(0,0,1.2,0), fig=c(0,0.4,0,1))
par(new=TRUE)
par(las=0, xpd=FALSE)
test_preyPrey8 <- test_prey * test_predcase8
test8_index <- (test_preyPrey8)>thisCap; test_preyPrey8[test8_index]<- thisCap
plot(x=test_prey, y=test_preyPrey8, type="l", lwd=2, col=std_col, ylim=c(0, thisCap), bty="n", lty=2,
     yaxt="n", xaxt="n", ylab="", xlab="")
for(c in 1:nCs){
  points(x=test_prey, y=predcase8ByC[c,], type="l", lwd=2, col=colByC[c])
}
mtext("A", side=3, adj=0, line=0.1, font=2, cex=1.2)
mtext("Prey biomass attacked",side=2, line=1, cex=1.5)
mtext("Prey biomass available", side=1, line=1, cex=1.5)
arrows(x0=0,y0=0, x1=0, y1=thisCap, lwd=2)
arrows(x0=0,y0=0, x1=max(test_prey)*1.01, y1=0, lwd=2)
abline(h=thisCap, col=cap_col, lty=3, lwd=2)
par(las=1)
axis(at=thisCap, label="Cap", side=2, col.axis=cap_col, cex=2, col.lab=myBlue, tick = FALSE, line=-0.5, cex.axis=1.2)
# legend(title="Clearance scalar", cex=1.2, legend=formatC(c_scalars,  zero.print = TRUE, digits=1, format="f"), col=colByC, lwd=2, lty=1, x="bottomright", bty="n", inset=0.05, seg.len=3)

# C:ht effects
c_index <- round(seq(1,nCs, length.out=3))
sample_cs <- c_scalars[c_index]
test_C <- sample_cs * this_mum*numsec; C_secs <- test_C/numsec;
sample_ht <- c(0.5,1,2); ht_lty<-c(2,1,4)
test_ht <- sample_ht*this_ht
predcase8ByC <- array(NA, dim=c(3,3, length(test_prey))); predRate8ByC <- predcase8ByC
for(c in 1:3){
  for(h in 1:3){
    this_predcase <- unlist(lapply(test_prey, FUN=function(x){predRes(C=C_secs[c], prey_avail = x, efficiencies = efficiencies, prey_types = prey_types, ht=test_ht[h], mum=this_mum, predcase=8)}))
    test_preyPrey8 <- test_prey * this_predcase
    test8_index <- (test_preyPrey8)>thisCap; test_preyPrey8[test8_index]<- thisCap
    predcase8ByC[c,h,]<-test_preyPrey8
    predRate8ByC[c,h,] <- this_predcase
  }
}

par(las=0, xpd=FALSE, mar=c(4,4,1,1), fig=c(0.4,0.8,0,1))
par(new=TRUE)
test_preyPrey8 <- test_prey * test_predcase8
test8_index <- (test_preyPrey8)>thisCap; test_preyPrey8[test8_index]<- thisCap
par(las=0, xpd=FALSE)
plot(x=test_prey, y=test_preyPrey8, type="n", lwd=2, col=std_col, ylim=c(0, thisCap), bty="n", lty=2,
     yaxt="n", xaxt="n", ylab="", xlab="")
for(c in 1:3){
  for(h in 1:3){
    points(x=test_prey, y=predcase8ByC[c,h,], type="l", lwd=2, col=colByC[c_index][c], lty=ht_lty[h])
  }
}
mtext("B", side=3, adj=0, line=0.1, font=2, cex=1.2)
mtext("Prey biomass attacked",side=2, line=1, cex=1.5)
mtext("Prey biomass available", side=1, line=1, cex=1.5)
arrows(x0=0,y0=0, x1=0, y1=thisCap, lwd=2)
arrows(x0=0,y0=0, x1=max(test_prey)*1.01, y1=0, lwd=2)
abline(h=thisCap, col=cap_col, lty=3, lwd=2)
par(las=1)
axis(at=thisCap, label="Cap", side=2, col.axis=cap_col, cex=2, col.lab=myBlue, tick = FALSE, line=-0.5, cex.axis=1.2)
# legend(title="Clearance scalar", cex=1.2, legend=formatC(c_scalars,  zero.print = TRUE, digits=1, format="f"), col=colByC, lwd=2, lty=1, x="bottomright", bty="n", inset=0.05, seg.len=3)
legendText <- c(paste(c_scalars[c_index][1],"  \n",sep=""), " ", " ", paste(c_scalars[c_index][2],"  \n",sep=""), " ", " ",paste(c_scalars[c_index][3],"  \n",sep=""), " ", " ")
legendCols <- sort(rep(colByC[c_index],3))
legendLtys <- rep(ht_lty,3)
par(fig=c(0.8,1,0,1))
par(new=TRUE)
par(mar=c(0,0,0,0))
makeBlankPlot()
legend(title="Clearance scalar", cex=1.8, legend=formatC(c_scalars,  zero.print = TRUE, digits=1, format="f"), col=colByC, lwd=2, lty=1, x="top", bty="n", inset=0.05, seg.len=3)
legend(title="Handling time scalar", cex=1.8, legend=formatC(test_ht/this_ht,  zero.print = TRUE, digits=1, format="f"), col=myGrey, lwd=2, lty=ht_lty, x="bottom", bty="n", inset=0.05, seg.len=3)

dev.off()


## do alt using 1/mum scalar with modified version






# changing handling time - prey attacked/available on y-axis, pray available on x-axis
# par(mfrow=c(1,2), mar=c(4,4,1,1))
# par(las=0, xpd=FALSE)
# plot(x=test_prey, y=test_predcase8, col=std_col, type="l", lwd=2.5, yaxt="n", xaxt="n", xlab="", ylab="", bty="n", ylim=c(0,max(test_predcase8)), lty=2)
# abline(h=thisCap, lty=3, col=cap_col, lwd=2)
# mtext("Prey biomass attacked/available",side=2, line=1, cex=1.5)
# mtext("Prey biomass available", side=1, line=1, cex=1.5)
# par(xpd=NA)
# arrows(x0=-0.2,y0=0, x1=-0.2, y1=max(test_predcase8), lwd=2)
# arrows(x0=-0.2,y0=-0, x1=max(test_prey)*1.01, y1=-0, lwd=2)
# for(h in 1:length(test_ht)){
#   points(x=test_prey, y=predcase8ByHt[h,]/test_prey, type="l", lwd=2, col=colByHt[h])
# }




###################################################
# HAK cohort 1s have SN = 1600, hta=28, htb=1.1
hake_hta <- 28; hake_htb <- 1.1; hake_SN <- 1600
hake_ht <- hake_hta * (hake_SN * X_CN) ^ (-hake_htb) * numsec

C=0.1; mum=5e-4; ht=100; predcase=1
prey_avail<-runif(5); prey_types <- c(1,1,1,2,1); efficiencies <- c(0.5,0.2,0)
# using predRes, mum won't affect predcase 8 at all, will only affect predcase 1
# the oposite for handling time, which doesn't affect pred-prey response when using predcase 1
# C affects both

orig_C <- seq(1,1000, length.out=50); test_C <- orig_C/numsec
test_predcase1 <- unlist(lapply(test_C, FUN=function(x){predRes(C=x, prey_avail = prey_avail, efficiencies = efficiencies, prey_types = prey_types, ht=ht, mum=mum, predcase=1)}))
test_predcase8 <- unlist(lapply(test_C, FUN=function(x){predRes(C=x, prey_avail = prey_avail, efficiencies = efficiencies, prey_types = prey_types, ht=ht, mum=mum, predcase=8)}))
axis1 <- pretty( seq(min(test_predcase1),max(test_predcase1), length.out = 5))

par(mar=c(4,4,1,4))
plot(x=orig_C, y=test_predcase8, pch=20, ylab="Attack rate", xlab="Clearance rate")
par(new=TRUE)
plot(x=orig_C, y=test_predcase1, pch=20, col=myBlue, yaxt="n", ylab="", xaxt="n", xlab="")
axis(at=axis1, labels=axis1, side=4, col.axis=myBlue, line=0)
legend(legend=c("Modified Holling Type II", "Holling Type II"), col=c(myBlue,"black"), pch=20, x="bottomright")
mtext("Attack rate", side=4, adj=0.5, line=3, col=myBlue)


# for a given clearance, vary total prey and plot these
orig_C <- c(1,10,100,1000); test_C <- orig_C/numsec; nC <- length(test_C)
test_prey <- seq(0,1,length.out=100)
c <- 3
test_predcase1 <- unlist(lapply(test_prey, FUN=function(x){predRes(C=test_C[c], prey_avail = x, efficiencies = efficiencies, prey_types = prey_types, ht=ht, mum=mum, predcase=1)}))
plot(x=test_prey, y=test_predcase1, type="l")
#amount of prey attacked
test_preyAttacked <- test_predcase1 * test_prey
plot(x=test_prey, y=test_preyAttacked, type="l")
par(mfrow=c(2,2))
store_testPred1 <- array(NA, dim=c(length(test_prey), nC)); store_testPred8 <- 0*store_testPred1
for(c in 1:nC){
  test_predcase1 <- unlist(lapply(test_prey, FUN=function(x){predRes(C=test_C[c], prey_avail = x, efficiencies = efficiencies, prey_types = prey_types, ht=ht, mum=mum, predcase=1)}))
  # plot(x=test_prey, y=test_predcase1, type="l")
  #amount of prey attacked
  test_preyAttacked <- test_predcase1 * test_prey
  store_testPred1[,c] <- test_preyAttacked
  # plot(x=test_prey, y=test_preyAttacked, type="l")
  # actual HII
  test_predcase8 <- unlist(lapply(test_prey, FUN=function(x){predRes(C=test_C[c], prey_avail = x, efficiencies = efficiencies, prey_types = prey_types, ht=ht, mum=mum, predcase=8)}))
  # plot(x=test_prey, y=test_predcase8, type="l")
  #amount of prey attacked
  test_preyAttacked <- test_predcase8 * test_prey
  store_testPred8[,c] <- test_preyAttacked
  # plot(x=test_prey, y=test_preyAttacked, type="l")
}
colByC <- colorRampPalette(colors=c(myGold, myGreen, myBlue))(nC)
thisYmax <- max(max(store_testPred1, na.rm=TRUE), max(store_testPred8, na.rm=TRUE))
plot(x=test_prey, y=store_testPred1[,c], ylim=c(0, thisYmax), type="l")
for(c in 1:nC){
  points(x=test_prey, y=store_testPred1[,c], type="l", col=colByC[c], lwd=1.5)
}
plot(x=test_prey, y=store_testPred8[,c], ylim=c(0, thisYmax), type="l")
for(c in 1:nC){
  points(x=test_prey, y=store_testPred8[,c], type="l", col=colByC[c], lwd=1.5)
}

## now bring in varying mum - which is applied differently for the 2 methods
E1 <- efficiencies[prey_types==1][1]
orig_mum <- c(1,5,10,50,100,1000); test_mum <- orig_mum/numsec
store_testPredMum1<-array(NA, dim=c(length(test_mum), length(test_prey), nC))
store_testPredMum8 <- 0*store_testPredMum1
for(m in 1:length(test_mum)){
  mum <- test_mum[m]
  for(c in 1:nC){
    thisCap <- mum/E1
    test_predcase1 <- unlist(lapply(test_prey, FUN=function(x){predRes(C=test_C[c], prey_avail = x, efficiencies = efficiencies, prey_types = prey_types, ht=ht, mum=mum, predcase=1)}))
    pred1Index <- test_predcase1>thisCap
    test_predcase1[pred1Index]<- thisCap
    preyPerPred1 <- test_predcase1 * test_prey
    test_predcase8 <- unlist(lapply(test_prey, FUN=function(x){predRes(C=test_C[c], prey_avail = x, efficiencies = efficiencies, prey_types = prey_types, ht=ht, mum=mum, predcase=8)}))
    preyPerPred8 <- test_predcase8 * test_prey
    pred8Index <- preyPerPred8>thisCap
    preyPerPred8[pred8Index]<- thisCap
    store_testPredMum1[m,,c] <- preyPerPred1
    store_testPredMum8[m,,c] <- preyPerPred8
  }
}
thisYmax <- max(max(store_testPredMum1, na.rm=TRUE), max(store_testPredMum8, na.rm=TRUE))

for(m in 1:length(test_mum)){
  plot(x=test_prey, y=store_testPredMum1[m,,1], type="l", ylim=c(0, thisYmax))
  for(c in 1:nC){
    points(x=test_prey, y=store_testPredMum1[m,,c], type="l", col=colByC[c], lwd=1.5)
  }
  plot(x=test_prey, y=store_testPredMum8[m,,1], type="l", ylim=c(0, thisYmax))
  for(c in 1:nC){
    points(x=test_prey, y=store_testPredMum8[m,,c], type="l", col=colByC[c], lwd=1.5)
  }
}

## lastly: handling time - only affects pred case 8
ht_test <- c(10,50,100,500); nH <- length(ht_test)
store_testPredMumHt8 <- array(NA, dim=c(nH, length(test_mum), length(test_prey), nC))
for(h in 1:nH){
  ht <- ht_test[h]
  for(m in 1:length(test_mum)){
    mum <- test_mum[m]
    for(c in 1:nC){
      thisCap <- mum/E1
      test_predcase8 <- unlist(lapply(test_prey, FUN=function(x){predRes(C=test_C[c], prey_avail = x, efficiencies = efficiencies, prey_types = prey_types, ht=ht, mum=mum, predcase=8)}))
      preyPerPred8 <- test_predcase8 * test_prey
      pred8Index <- preyPerPred8>thisCap
      preyPerPred8[pred8Index]<- thisCap
      store_testPredMumHt8[h,m,,c] <- preyPerPred8
    }
  }
}

thisYmax <- max(max(store_testPredMum1, na.rm=TRUE), max(store_testPredMumHt8, na.rm=TRUE))
maxByMum1 <- apply(store_testPredMum1, 1, max); maxByMum8 <- apply(store_testPredMumHt8,2, max)
maxByMum <- apply(cbind(maxByMum1, maxByMum8), 1, max)
for(m in 1:length(test_mum)){
  plot(x=test_prey, y=store_testPredMum1[m,,1], type="l", ylim=c(0, maxByMum[m]))
  for(c in 1:nC){
    points(x=test_prey, y=store_testPredMum1[m,,c], type="l", col=colByC[c], lwd=1.5)
  }
  plot(x=test_prey, y=store_testPredMum8[m,,1], type="l", ylim=c(0, maxByMum[m]))
  for(c in 1:nC){
    for(h in 1:nH){
      points(x=test_prey, y=store_testPredMumHt8[h,m,,c], type="l", col=colByC[c], lwd=1.5, lty=h)
    }
  }
}




axis1 <- pretty( seq(min(test_predcase1),max(test_predcase1), length.out = 5))




## an aside, to check growth pars
# fprintf(llogfp,
#         "EAT parameters: cohort = %d, sp = %.20e, C_sp = %.20e, mum_sp = %.20e,  ht_sp = %.20e, E1_sp = %.20e, E2_sp = %.20e\n",
#         cohort, sp, C_sp, mum_sp,  ht_sp, E1_sp, E2_sp);
logLines <- readLines(paste(DIR$'Base',"AtlantisModels\\Base\\outputTestGrowthPars\\log.txt", sep=""))
x <- grep("EAT parameters:", logLines)
growthLogLines <- logLines[x]

# test also initial conditions to get SN for handling time
ThisNC.nc <- nc_open(paste(DIR$'Base',"AtlantisModels\\Base\\outputTestGrowthPars\\output.nc", sep=""))
thisData <- ncvar_get(ThisNC.nc, "Hake1_StructN")
# # HAK cohort 1s have SN = 1600, hta=28, htb=1.1
# X_CN = 5.7
# hake_hta <- 28; hake_htb <- 1.1; hake_SN <- 1600
# hake_ht <- hake_hta * (hake_SN * X_CN) ^ (-hake_htb) * numsec
# 
# handling time is in seconds - ht = hta * (SN * X_CN ) ^ (-htb) *numsec
# C is in seconds
# mum is in seconds

# what is a typical mg N per m^3..?
allVars <- sort(unique(names(ThisNC.nc$var)))
x <- grep('_N', allVars)
temp<-allVars[x]
xx <- grep("_Nums", temp, invert = TRUE)
Nvars <- temp[xx]; nNV <- length(Nvars)
storeNs <- c()
for(v in 1:nNV){
  thisVar <- Nvars[v]
  thisData <- ncvar_get(ThisNC.nc, thisVar)
  if(length(dim(thisData))==3){
    thisVector <- as.vector(thisData[,,1])
  } else{
    thisVector <- as.vector(thisData[,1])
  }
  thisMax <- max(thisVector, na.rm=TRUE)
  if(thisMax>1000){cat("Check on ", thisVar,"--")}
  thisVector <- thisVector[thisVector>1e-8]
  storeNs <- c(storeNs,thisVector)
}
hist(storeNs)






