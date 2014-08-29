rm(list=c(ls()))
library(igraph)
library(Rcpp)

#This is to analyze the properties of multiple Tree-like nets...

## Load Source Files
##

sourceCpp('src/DLA_func.cpp')
L = 32
dla <- DLA_func(size=(L*L),boundary=20)


sourceCpp('src/DLA_func2.cpp')
L = 32
dla <- DLA_func2(size=(L*L),boundary=80)

sourceCpp('src/DLA_func_radial.cpp')
L = 32
dla <- DLA_func_radial(size=(L*L),boundary=20)

dla.obj <- dla[[1]]
g <- graph.edgelist(dla.obj)
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
dla.area <- r.area
edge.weight <- numeric(length(dla.obj[,1]))
for (i in 1:length(dla.obj[,1])) {
  edge1 <- dla.obj[i,1]
  edge2 <- dla.obj[i,2]
  edge.weight[i] <- min(c(r.area[edge1],r.area[edge2]))
}
plot(g,layout=dla[[2]],vertex.size=0.01,
     vertex.color="transparent",vertex.frame.color=NA,edge.color="dodgerblue3",edge.arrow.mode=0,
     edge.width=log(edge.weight+1)/2,vertex.label=NA)



sourceCpp('src/dyn_func.cpp')
sourceCpp('src/DLA_func2.cpp')
source("RivOCN_func.R")

early.riv.net <- list()
late.riv.net <- list()
dla.net <- list()
L <- 32

for (i in 1:50) {
  
  print(paste("Repetition = ",i,"/50"))
  riv.net <- RivOCN_func(L=L,tic_max=10000)
  dla <- DLA_func2(size=(L*L),boundary=20)
  
  early.riv.net[[i]] <- riv.net$early.net
  late.riv.net[[i]] <- riv.net$late.net
  dla.net[[i]] <- graph.edgelist(dla[[1]])
  
}
save.image("tree_net_properties.RData")


