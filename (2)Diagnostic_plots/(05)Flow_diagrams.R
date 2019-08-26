create_diamond<-function(center,height=1,length=2){
  if(length(center)==1){
    xCenter<-center; yCenter=center;
  }else{
    xCenter<-center[1]
    yCenter<-center[2]
  }
  x1<-xCenter-0.5*length; x2<-xCenter; x3<-xCenter+0.5*length; x4<-xCenter
  y1<-yCenter;  y2<-yCenter+0.5*height; y3<-yCenter;  y4<-yCenter-0.5*height
  
  out<-cbind("x"=c(x1,x2,x3,x4),"y"=c(y1,y2,y3,y4))
  return(out)
}

create_arrow<-function(center,height=2,length=3,direction="right"){
  xCenter<-center[1]
  yCenter<-center[2]
  x1<-xCenter-0.5*length; x2<-xCenter+0.15*length; x3<-x2; x4<-xCenter+0.5*length;
  x5<-x3; x6<-x3; x7<-x1; x8<-xCenter-0.15*length
  y1<-yCenter+0.3*height; y2<-y1; y3<-yCenter+0.5*height; y4<-yCenter
  y5<-yCenter-0.5*height; y6<-yCenter-0.3*height; y7<-y6
  
  xRight<-c(x1,x2,x3,x4,x5,x6,x7)
  yRight<-c(y1,y2,y3,y4,y5,y6,y7)
  
  xLeft<-c(x1,x8,x8,x4,x4,x8,x8)
  yLeft<-c(y4,y5,y6,y6,y1,y1,y3)
  
  if(direction=="right"){
    out<-cbind("x"=xRight,"y"=yRight)
  } else if(direction=="left"){
    out<-cbind("x"=xLeft,"y"=yLeft)
  } else if(direction=="up"){
    out<-cbind("x"=yRight,"y"=xRight)
  } else{
    out<-cbind("x"=yLeft,"y"=xLeft)
  }

  return(out)
}

create_rectangle<-function(xCenter,yCenter,length=1,height=1){
  x1<-xCenter-0.5*length;   x2<-x1
  x3<-xCenter+0.5*length;   x4<-x3
  y1<-yCenter-0.5*height
  y2<-yCenter+0.5*height
  y3<-y2;y4<-y1
  out<-cbind("x"=c(x1,x2,x3,x4),"y"=c(y1,y2,y3,y4))
  return(out)
}


thisTextCex<-0.9

# pdf(paste(DIR$'Figures',"FlowChartBH_InitialPars.pdf",sep=""),height=6)
jpeg(paste(DIR$'Figures',"FlowChartBH_InitialPars.jpg",sep=""),height=500,quality=400)

par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(x=1,y=1,ylim=c(0,5),xlim=c(-0.5,5),xaxt="n",yaxt="n",xlab="",ylab="",type="n",bty="n")
a1<-create_arrow(center=c(0.75,4),direction="right",height=1,length=1.9)
polygon(a1,col=myBlue_trans,border=myBlue)
text("IN\nCatch histories\nBiological parameters",x=0.5,y=4.02,adj=0.5,cex=thisTextCex)

d1<-create_diamond(c(2.8,4),height=1.5,length=1.5*(10/7))
polygon(d1,col=myOrange_trans,border=myOrange)
text("Bayesian\nage-structured\nmodel",x=2.8,y=4,cex=thisTextCex)

a2<-create_arrow(center=c(2.6,2.8),height=1,length=1.0,direction="down")
polygon(a2,col=myOrange_trans,border=myOrange)
text(expression("B"[0]),x=2.8,y=2.7,cex=thisTextCex)
text(expression("R"[0]),x=2.8,y=2.5,cex=thisTextCex)


d2<-create_diamond(center=c(1.7,1.8),height=1.5,length=1.5*(10/7))
polygon(d2,col=myOrange_trans,border=myOrange)
text("Beverton-Holt\nequations",x=1.7,y=1.8,cex=thisTextCex)

d3<-create_diamond(center=c(3.9,1.8),height=1.5,length=1.5*(10/7))
polygon(d3,col=myOrange_trans,border=myOrange)
text("B0 distributed over\nage-classes by\nnumbers and mgN\nper individual",x=3.9,y=1.8,cex=thisTextCex)

a3<-create_arrow(center=c(0.0,1.8),length=1.2,height=1)
polygon(a3,col=myBlue_trans,border=myBlue)
text("IN\nSteepness (h)",x=-0.1,y=1.8,cex=thisTextCex)

r1<-create_rectangle(xCenter=1.7,yCenter=0.6,length=1.8,height=0.8)
polygon(r1,col=myGreen_trans,border=myGreen)
text("ATLANTIS PARAMETERS\nBeverton-Holt recruitment",x=1.7,y=0.6,cex=thisTextCex,family="serif")
text("a, b",x=1.7,y=0.35,cex=thisTextCex*1.2,family="serif")


r2<-create_rectangle(xCenter=3.9,yCenter=0.6,length=1.8,height=0.8)
polygon(r2,col=myGreen_trans,border=myGreen)
text("ATLANTIS PARAMETERS\nInitial Conditions\n",x=3.9,y=0.6,cex=thisTextCex)

text(expression("N"[0[a]]),x=3.7,y=0.4,cex=thisTextCex*1.2,family="serif")
text(expression(mu[w[a]]),x=4.1,y=0.4,cex=thisTextCex*1.2,family="serif")
##
dev.off()


#DIET STUDIES TO PREDATOR PREY
pdf(paste(DIR$'Figures',"FlowChartDIET.pdf",sep=""),height=6)
par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(x=1,y=1,ylim=c(0,5),xlim=c(-0.5,5),xaxt="n",yaxt="n",xlab="",ylab="",type="n",bty="n")
a1<-create_arrow(center=c(0.75,4),direction="right",height=1,length=1.9)
polygon(a1,col=myBlue_trans,border=myBlue)
text("IN\nDiet studies\nanalysed data",x=0.5,y=4.02,adj=0.5,cex=thisTextCex)

d1<-create_diamond(c(2.8,4),height=1.5,length=1.5*(10/7))
polygon(d1,col=myOrange_trans,border=myOrange)
text("Summarise w.r.t\nmodel groups",x=2.8,y=4,cex=thisTextCex)

a2<-create_arrow(center=c(2.6,2.8),height=1,length=1.0,direction="down")
polygon(a2,col=myOrange_trans,border=myOrange)
text("")
text(expression("B"[0]),x=2.8,y=2.7,cex=thisTextCex)
text(expression("R"[0]),x=2.8,y=2.5,cex=thisTextCex)




# thisTextCex<-0.8
# 
# # pdf(paste(DIR$'Figures',"\\Flowchart_B0R0.pdf",sep=""),height=3)
# par(mar=c(0,0,0,0),oma=c(0,0,0,0))
# plot(x=1,y=1,ylim=c(3,6),xlim=c(0,5.5),xaxt="n",yaxt="n",xlab="",ylab="",type="n",bty="n")
# polygon(poly1)
# text("IN\nCatch histories\nBiological parameters",x=0.75,y=4.5,adj=0.5,cex=thisTextCex)
# polygon(poly2)
# text("Bayesian age-based model",x=1.8,y=4.5,adj=0,cex=thisTextCex)
# polygon(poly3)
# text("OUT\nUnfished biomass\nRecruitment at\nunfished biomass",x=4.5,y=4.6,cex=thisTextCex)
# dev.off()
