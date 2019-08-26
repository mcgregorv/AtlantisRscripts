library(igraph)
library(RColorBrewer)
library(plotrix)
library(shape)


this_out<-"TestZG"; runFolder="TBGB_JP2"
thisDesc <- paste(runFolder, this_out,sep="")

this_path = paste(DIR$'Base', "TBGB\\",runFolder,"\\",sep="")
outPath<-paste(this_path,"output",this_out,"\\",sep="")

plotPath<-paste(DIR$'Base',"TBGB\\Figures\\Testing\\Diet\\",thisDesc,sep="")

# read in diet snapshot and create dat
snapshotDiet <- read.csv(paste(outPath, "realiseDietSnapshot1880_1960.csv", sep=""))
dat <- melt(snapshotDiet); colnames(dat)<-c("Pred", "Prey", "value")
col.b <- brewer.pal(9,'GnBu')
tl  <-  read.csv(paste(outPath,'TBGBGroupsTL.csv', sep=""))[,-1]
## I will remove bacterias and detritus
## I dont want them in the plot
tl  <- tl[- which(tl[, 1] %in% c('BB','DR', 'DL', 'PB')), ]
dat <- dat[- which(dat[, 2] %in% c('BB','DR', 'DL', 'PB')), ]
dat<-dat[!is.na(dat$value),]
## MAtcing Trophic lvl with the
dat$TL   <- NA
dat$TLprey <- NA
for(fg in 1 : nrow(tl)){
  pred.fg              <- which(dat$Pred %in% tl$Code[fg])
  prey.fg              <- which(dat$Prey %in% tl$Code[fg])
  dat$TL[pred.fg]      <- tl$TL[fg]
  dat$TLprey[prey.fg]  <- tl$TL[fg]
  if(fg == nrow(tl)){
    dat2 <- data.frame(Pred = c('MA', 'PL', 'PS'),
                       value = c(1, 1, 1),
                       Prey = c('MA', 'PL', 'PS'),
                       TL   = c(2, 2, 2),
                       TLprey = c(2, 2, 2))
    dat <- rbind(dat,dat2)
  }
}

## sorting the talbe based on the trophic level
tl.p <- tl$TL[order(tl$TL, decreasing = TRUE)]

pdf(paste(plotPath,'foodweb.pdf', sep=""), width = 10,height = 10)
par(oma = c(1, 1,1, 1), mar = c(0,0,0,0))
plotDat <- dat[with(dat, order(TL, decreasing = TRUE)), c(1, 2, 4)]
col.p <- ifelse(plotDat$TL < 2, 1, ifelse(plotDat$TL < 3, 3,
                                      ifelse(plotDat$TL < 4, 5, ifelse(plotDat$TL < 4.5, 7, 9))))
col.p2 <- ifelse(tl.p < 2, 1,
                 ifelse(tl.p < 3, 2,
                        ifelse(tl.p < 4, 4,
                               ifelse(tl.p < 4.5, 6,
                                      ifelse(tl.p < 6, 7,9)))))
dat2 <- graph_from_data_frame(plotDat)
dat2 <- simplify(dat2, remove.multiple = TRUE, remove.loops = TRUE,
                 edge.attr.comb = igraph_opt("edge.attr.comb"))
for( i in 1 : 100){ ## permuting to get the edges properly
  dat2 <- permute(dat2, order(tl$TL[match(names(degree(dat2)), tl$Code)], decreasing = TRUE))
  dat2order <- names(degree(dat2))
  test<- dat2order == tl$Code[order(tl$TL, decreasing = TRUE)]
  # if(length(dat2order[!test])<10){
  #   cat(" almost there--", length(dat2order[!test]),"-, ")
  # }
}
plot(dat2, edge.arrow.size = .3, curve = 5, vertex.color=col.b[col.p2], frame.color = col.b[9],
     label.color = col.b[1], edge.color = col.b[9], layout = layout_in_circle)
legend('topleft', legend = c('1 - 1.9', '2 - 2.9', '3 - 3.9', '4 - 4.49', '4.5 - 4.99'), title = 'Trophic Level', cex = 1.3,
       pch = 20, bty = 'n', pt.cex = 3, col = col.b[c(1, 2, 4, 6, 7)])
dev.off()

# do with focus group colored differently
plotDat <- dat; plotDat$PredCode <- dat$Pred; plotDat$PreyCode <- dat$Prey
mygraph<- graph_from_data_frame(plotDat, directed=FALSE)
mySimpleGraph <- simplify(mygraph,  edge.attr.comb="sum")
thisEdgeCol <- rep(myGrey_trans, length(E(mygraph)))
thisEdgeCol[E(mygraph)$PredCode==focusGroup]<-col.b[6]
thisEdgeCol[E(mygraph)$PreyCode==focusGroup]<-col.b[9]
plot(mygraph, edge.arrow.size = .3, curve = 5, vertex.color=col.b[col.p2], frame.color = col.b[col.g],
     label.color = col.b[1], edge.color = thisEdgeCol, layout = layout_in_circle)

E(mySimpleGraph)
E(mygraph)



pdf(paste(plotPath,focusGroup,'foodweb.pdf', sep=""), width = 10,height = 10)
par(oma = c(1,1,1,1), mar = c(0,0,0,0))
col.p <- ifelse(datSca$TL < 2, 1, ifelse(datSca$TL < 3, 3,
                                      ifelse(datSca$TL < 4, 5, ifelse(datSca$TL < 4.5, 7, 9))))
col.p2 <- ifelse(tl.p < 2, 1,
                 ifelse(tl.p < 3, 2,
                        ifelse(tl.p < 4, 4,
                               ifelse(tl.p < 4.5, 6,
                                      ifelse(tl.p < 6, 7,9)))))
# col.g <- ifelse(dat$Pred ==focusGroup | dat$Prey==focusGroup, 2, 1) # 2 colours - connected to SCA or not
col.g <- ifelse(datSca$Pred ==focusGroup,3, ifelse(datSca$Prey==focusGroup, 2, 1)) # 3 colours - prey=sca, pred=sca, neither

keepIndex <- as.character(datSca$Pred)==as.character(datSca$Prey) | datSca$Pred==focusGroup | datSca$Prey==focusGroup
dat2 <- graph_from_data_frame(dat)
dat2 <- simplify(dat2, remove.multiple = TRUE, remove.loops = TRUE,
                 edge.attr.comb = igraph_opt("edge.attr.comb"))
for( i in 1 : 10){ ## permuting to get the edges properly
  dat2 <- permute(dat2, order(tl$TL[match(names(degree(dat2)), tl$Code)], decreasing = TRUE))
}
plot(dat2, edge.arrow.size = .3, curve = 5, vertex.color=col.b[col.p2], frame.color = col.b[col.g],
     label.color = col.b[1], edge.color = colLines[col.g], layout = layout_in_circle)
legend('topleft', legend = c('1 - 1.9', '2 - 2.9', '3 - 3.9', '4 - 4.49', '4.5 - 4.99'), title = 'Trophic Level', cex = 1.3,
       pch = 20, bty = 'n', pt.cex = 3, col = col.b[c(1, 2, 4, 6, 7)])
dev.off()










############################################
scaDat <- melt(snapshotDiet); colnames(scaDat)<-c("Pred", "Prey", "value")
scaDat <- scaDat[!is.na(scaDat$value) & (scaDat$Pred=="SCA" | scaDat$Prey=="SCA"),]

scaDat$TL   <- NA
scaDat$TLprey <- NA
for(fg in 1 : nrow(tl)){
  pred.fg              <- which(scaDat$Pred %in% tl$Code[fg])
  prey.fg              <- which(scaDat$Prey %in% tl$Code[fg])
  scaDat$TL[pred.fg]      <- tl$TL[fg]
  scaDat$TLprey[prey.fg]  <- tl$TL[fg]
  if(fg == nrow(tl)){
    dat2 <- data.frame(Pred = c('MA', 'PL', 'PS'),
                       value = c(1, 1, 1),
                       Prey = c('MA', 'PL', 'PS'),
                       TL   = c(2, 2, 2),
                       TLprey = c(2, 2, 2))
    scaDat <- rbind(scaDat,dat2)
  }
}

## sorting the talbe based on the trophic level
tl.p <- tl$TL[order(tl$TL, decreasing = TRUE)]

# png(paste(plotPath,'SCAfoodweb.png', sep=""), width = 800,height = 800)
par(oma = c(1, 5, 5, 1), mar = c(1, 1, 1, 1))
dat <- scaDat[with(scaDat, order(TL, decreasing = TRUE)), c(1, 2, 4)]
col.p <- ifelse(dat$TL < 2, 1, ifelse(dat$TL < 3, 3,
                                      ifelse(dat$TL < 4, 5, ifelse(dat$TL < 4.5, 7, 9))))
col.p2 <- ifelse(tl.p < 2, 1,
                 ifelse(tl.p < 3, 2,
                        ifelse(tl.p < 4, 4,
                               ifelse(tl.p < 4.5, 6,
                                      ifelse(tl.p < 6, 7,9)))))
dat2 <- graph_from_data_frame(dat)
dat2 <- simplify(dat2, remove.multiple = TRUE, remove.loops = TRUE,
                 edge.attr.comb = igraph_opt("edge.attr.comb"))
for( i in 1 : 10){ ## permuting to get the edges properly
  dat2 <- permute(dat2, order(tl$TL[match(names(degree(dat2)), tl$Code)], decreasing = TRUE))
}
plot(dat2, edge.arrow.size = .3, curve = 5, vertex.color=col.b[col.p2], frame.color = col.b[9],
     label.color = col.b[1], edge.color = col.b[9], layout = layout_in_circle)
legend('topleft', legend = c('1 - 1.9', '2 - 2.9', '3 - 3.9', '4 - 4.49', '4.5 - 4.99'), title = 'Trophic Level', cex = 1.3,
       pch = 20, bty = 'n', pt.cex = 3, col = col.b[c(1, 2, 4, 6, 7)])
# dev.off()



