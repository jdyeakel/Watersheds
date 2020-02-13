library(RColorBrewer)
library(igraph)
library(intergraph)
library(Rcpp)

sourceCpp('/Users/justinyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/src/OCN_func.cpp')
L <- 10
its <- 10
cout <- OCN_func(L,its)

edgelist = cout[[3]] + 1
coords = cout[[2]] + 1

g <- graph_from_edgelist(cout[[3]]+1)

#LEVELS
#CM
#SC
#AG
#RN
#FD

library(igraph)
library(OCNet)
# OCN <- aggregate_OCN(landscape_OCN(OCN_250_T), thrA = 4)
OCN <- aggregate_OCN(landscape_OCN(create_OCN(50, 50)), thrA = 4)
# draw_simple_OCN(OCN)
g <- OCN_to_igraph(OCN, level = "FD")

link.tree = get.edgelist(g)
edge.weight <- numeric(length(link.tree[,1]))
r.area <- log(OCN$FD$A+1)
for (i in 1:length(link.tree[,1])) {
  edge1 <- link.tree[i,1]
  edge2 <- link.tree[i,2]
  edge.weight[i] <- min(c(r.area[edge1],r.area[edge2]))
}
#River plot
plot(g,layout=matrix(c(OCN$FD$X,OCN$FD$Y), ncol = 2, nrow = OCN$FD$nNodes),vertex.size=0.01,
     vertex.color="transparent",vertex.frame.color=NA,edge.color="dodgerblue3",edge.arrow.mode=0,edge.width=edge.weight,vertex.label=NA)
#River NETWORK plot
plot(g, layout = matrix(c(OCN$FD$X,OCN$FD$Y), ncol = 2, nrow = OCN$FD$nNodes),vertex.size=log(OCN$FD$A+1),vertex.color="dodgerblue3",edge.color="darkgrey",edge.arrow.size=0.3,vertex.label=NA)

#Area values
OCN$FD$A

plot(g,layout=coords)




