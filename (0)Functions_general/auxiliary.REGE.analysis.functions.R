
# function to get color labels
colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
  }
  n
}


make.igraph.object = function(App){
  library(igraph)
  # for TrophInd function
  library(NetIndices)
  pp <- graph.adjacency(App, mode="directed", 
         add.rownames=TRUE, diag=FALSE)
  pp$layout <- layout.kamada.kawai(pp)
  pp$edge.arrow.size=0.4
  pp <- set.graph.attribute(pp, "edge.color", "black")
  pp$edge.color="black"
  pp$edge.width=0.4
  pp
}


define.link = function(sub.block, method="any"){
  # Make a link if there are any non-zero elements in the block
  # Returned link values are either 0 or 1
  if(method=="any") link = as.numeric( any( c(as.logical(sub.block)) ) )
  link
}



subcluster.labels = function(hc, num.clusters){
  #
  # Bunch together those groups in the same cluster
  # If there is only one member in the cluster retain 
  # the original G1, G2, etc, short label. Otherwise make up
  # a new label NX where X is the cluster number.
  #
  # input:  hc is result from hierarchial clustering (hclust) on 
  #         similarity matrix    
  # 
  # output is:
  #
  # > sub.clusters[[1]]
  # $names
  # [1] "Birds (G1)"         "Baleen_whales (G3)"
  # 
  # $label
  # [1] "N1"
  # 
  # Usage:    form.subcluster(hc, num.clusters)
  #
  #  
  clus = sort( cutree(hc, k=num.clusters) )
  sub.clusters = list()
  for (cluster in 1:num.clusters){
    # names: e.g. Hoki (G5), want just G5 part for short.label
    this.cluster = names(clus[clus==cluster])
    if(length(this.cluster)==1){
      pos.left.bracket  = grepRaw("(", this.cluster, fixed=T)
      pos.right.bracket = grepRaw(")", this.cluster, fixed=T)
      this.short.label = substr(this.cluster,
        start=pos.left.bracket + 1, stop=pos.right.bracket - 1) 
      this.long.label  = this.cluster
    } else {
      this.short.label= paste0("N",cluster)
      this.long.label = this.short.label
    }
    sub.clusters[[cluster]] = list(names=this.cluster, 
      short.label=this.short.label, long.label=this.long.label)
  }
  
  return(sub.clusters)
}


sub.adjacency.matrix = function(hc, num.clusters, method="any"){
  #
  # clus (from cutree) gives the cluster number that a group is in:
  #
  # Birds (G1)            Baleen_whales (G3) 
  #         1                             1 
  #  
  clus = sort( cutree(hc, k=num.clusters) )
  sub.adj = array(data=NA, dim=c(num.clusters,num.clusters))
  for(row.num in 1:num.clusters){
    for(col.num in 1:num.clusters){
      row.label = names(clus[clus==row.num])
      long.col.label = names(clus[clus==col.num])
      pos.left.bracket  = sapply(long.col.label, 
                                 function(x) grepRaw("(", x, fixed=T))
      pos.right.bracket = sapply(long.col.label, 
                                 function(x) grepRaw(")", x, fixed=T))
      col.label = substr(long.col.label, start=pos.left.bracket + 1, 
                    stop=pos.right.bracket - 1)
      
      sub.block = App[row.label, col.label]
      sub.adj[row.num, col.num]  = define.link(sub.block,method=method)
    }
  }
  
 labels = subcluster.labels(hc, num.clusters)
 dimnames(sub.adj) = list( sapply(labels, function(x) x$long.label), 
     sapply(labels, function(x) x$short.label) ) 
 
 list(matrix=sub.adj, labels=labels)
}



plot.trophic.levels = function(App){
  
  num.groups = nrow(App)
  
  troph.App<-TrophInd(App)
  troph.level <- round(troph.App$TL -1,2)
  troph.level <- troph.level^1.8
  troph.level.int <- round(troph.level)
  
  pp = make.igraph.object(App)
  vertex.degree <- degree(pp, mode="total")
  
  pull.out.names = function(x) unlist(strsplit(x, split=" "))[1]
  long.names = sapply(row.names(App), pull.out.names)
  
  lay = matrix(nrow=num.groups, ncol=2)
  # lay[,1] = (vertex.degree/max(vertex.degree))^1.0
  lay[,1] = pp$layout[,1]
  lay[,2] = troph.level
  
  old.par = par(mar=c(.1,.1,.1,.1))
  plot.igraph(pp,layout=lay,
              # xlim=c(-1,1),
              # ylim=c(0,5),
              rescale=TRUE,
              # vertex.shape="rectangle",
              vertex.shape="none",
              vertex.color=labelColors[clusMember],
              vertex.size=20,
              vertex.size2 = 8,
              vertex.label=long.names,
              vertex.label.cex=0.8,
              edge.color="black",
              edge.arrow.size=0.5,
              edge.width=0.4,
              edge.curved=TRUE,
  )
  par(old.par)

}


make.qual.matrix = function(Adj){
  #
  # Adj:  adjacency matrix. Filled with 0 and 1 values, and as used
  #       for REGE analysis input
  #
  
  qual.matrix = -Adj + t(Adj)
  diag(qual.matrix) = -diag(Adj)

  qual.matrix
}


make.maple.input = function(qual.matrix, filename="maple.input.txt"){

#
#  Input is a qualitative community matrix consisting of
#  0 and 1 values. It is assumed the columns have some 
#  sort of short label (e.g G1, G2, etc)
#
#  qual.matrix:  community qualitative matrix
#
# 
#  The output is the Maple form for the matrix (ready for 
#  input to the online loop analysis site)
#
#  http://ipmnet.org/loop/loopanalysis.aspx
#
#   n:=5:A:=array(1..n,1..n,[[-1,-1,0,0,0],[1,0,-1,0,0],
#   [0,1,0,-1,0],[0,0,1,0,-1],[0,0,0,1,0]]);[A,B,C,D,E]
#
#

  file.create(filename)
  first.line = paste("n:=",nrow(qual.matrix),":A:=array(1..n,1..n,[", sep="")
  cat(first.line, file=filename, append=FALSE)
  num.rows = nrow(qual.matrix)
  for(row.index in 1:nrow(qual.matrix)){
    this.row = qual.matrix[row.index,]
    this.line = paste0("[",paste(this.row,collapse=","),"]")
    if(row.index<num.rows)  this.line = paste0(this.line,",")
    if(row.index==num.rows) this.line = paste0(this.line,"]);")
    cat(this.line, file=filename, append=TRUE)
  }
  last.line = paste(colnames(qual.matrix),collapse=",")
  last.line = paste0("[",last.line,"]")
  cat(last.line, file=filename, append=TRUE)
  cat("\nMaple input matrix written to file: ",filename,"\n\n")
}




