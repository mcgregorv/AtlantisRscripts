library(igraph)
library(RColorBrewer)
library(plotrix)
library(shape)
col.b <- brewer.pal(9,'GnBu')
dat <- read.csv('fw.csv')
tl  <-  read.csv('TL.csv')
## I will remove bacterias and detritus
## I dont want them in the plot
tl  <- tl[- which(tl[, 1] %in% c('BB','DR', 'DL', 'PB')), ]
dat <- dat[- which(dat[, 2] %in% c('BB','DR', 'DL', 'PB')), ]

## MAtcing Trophic lvl with the
dat$TL   <- NA
dat$TLprey <- NA
for(fg in 1 : nrow(tl)){
    pred.fg              <- which(dat$Pred %in% tl$FG[fg])
    prey.fg              <- which(dat$Prey %in% tl$FG[fg])
    dat$TL[pred.fg]      <- tl$Tlevel[fg]
    dat$TLprey[prey.fg]  <- tl$Tlevel[fg]
    if(fg == nrow(tl)){
        dat2 <- data.frame(Pred = c('MA', 'LPH', 'SPH'),
                           value = c(1, 1, 1),
                           Prey = c('MA', 'LPH', 'SPH'),
                           TL   = c(2, 2, 2),
                           TLprey = c(2, 2, 2))
        dat <- rbind(dat,dat2)
    }
}

## sorting the talbe based on the trophic level
tl.p <- tl$Tlevel[order(tl$Tlevel, decreasing = TRUE)]



png('try.png', width = 800,height = 800)
par(oma = c(1, 5, 5, 1), mar = c(1, 1, 1, 1))
dat <- dat[with(dat, order(TL, decreasing = TRUE)), c(1, 2, 4)]
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
    dat2 <- permute(dat2, order(tl$Tlevel[match(names(degree(dat2)), tl$FG)], decreasing = TRUE))
}
plot(dat2, edge.arrow.size = .3, curve = 5, vertex.color=col.b[col.p2], frame.color = col.b[9],
     label.color = col.b[1], edge.color = col.b[9], layout = layout_in_circle)
legend('topleft', legend = c('1 - 1.9', '2 - 2.9', '3 - 3.9', '4 - 4.49', '4.5 - 5'), title = 'Trophic Level', cex = 1.3,
       pch = 20, bty = 'n', pt.cex = 2, col = col.b[c(1, 2, 4, 6, 7)])
dev.off()
