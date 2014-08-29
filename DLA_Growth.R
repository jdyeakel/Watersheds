rm(list=c(ls()))
library(igraph)
library(Rcpp)

L=32
sourceCpp('src/DLA_func.cpp')
output <- DLA_func(L*L,30)
link.tree <- output[[1]]
dla.tree <- output[[2]]

g <- graph.edgelist(link.tree)

# plot(g, layout=dla.tree,
#      vertex.color="lightblue",edge.color="darkgrey",edge.arrow.mode=0,edge.arrow.size=0.3,
#      vertex.label=NA,vertex.size=1)


# Calculate area
# First we need to orient all directed edges to flow towards the confluence
# Identify Confluence and tributaries
size <- vcount(g)
tribs <- which(degree(g)==1)
conf <- tribs[1]
tribs <- tribs[2:length(tribs)]
lt <- length(tribs)
#Build (accurate) directed net
dadj.m <- matrix(0,size,size)
for (i in 1:lt) {
  trib.i <- tribs[i]
  #find the path from trib to confluence (w.o.r. to flow)
  path <- get.shortest.paths(g,trib.i,conf,mode="all")[[1]]
  lp <- length(path)
  #Build the directed adjacency matrix
  for (j in 1:(lp-1)) {
    to.node <- path[j]
    from.node <- path[j+1]
    dadj.m[from.node,to.node] <- 1
  }  
}
dg <- graph.adjacency(dadj.m)

#Assign River Area (second try)
gamma <- 0.5
r.area <- numeric(size)
paths <- shortest.paths(dg,conf,mode="out")
longest.path <- which(paths==max(paths))
lpath <- paths[longest.path[1]]
for (i in seq(lpath,0,-1)) {
  nodes.dist <- which(paths == i)
  num.nodes <- length(nodes.dist)
  for (j in nodes.dist) {
    incoming.nodes <- neighbors(dg,j,mode="out")
    if (length(incoming.nodes) > 0) {
      area.neighbors <- r.area[incoming.nodes]
      r.area[j] <- sum(area.neighbors) + 1
    } else {
      r.area[j] <- 1
    }    
  }
}

edge.weight <- numeric(length(link.tree[,1]))
for (i in 1:length(link.tree[,1])) {
  edge1 <- link.tree[i,1]
  edge2 <- link.tree[i,2]
  edge.weight[i] <- min(c(r.area[edge1],r.area[edge2]))
}
plot(g,layout=dla.tree,vertex.size=0.01,
     vertex.color="transparent",vertex.frame.color=NA,edge.color="dodgerblue3",edge.arrow.mode=0,
     edge.width=log(edge.weight+1)/2,vertex.label=NA)

