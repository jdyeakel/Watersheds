rm(list=c(ls()))

library(igraph)
#setwd("/Users/jyeakel/Dropbox/PostDoc/2014_Empirical_Watersheds")

L <- 10


#Build Nearest Neighbor Information (Periodic Boundary conditions for LR)
riv.nn <- array(list(),c(L,L))
for (i in 1:L) {
  for (j in 1:L) {
    # In the center of the matrix
    if ((1 < i) && (i < L) && (1 < j) && (j < L)) { 
      nn.ij <- matrix(0,4,2)
      nn.ij[1,] <- c(i+0,j+1)
      nn.ij[2,] <- c(i+0,j-1)
      nn.ij[3,] <- c(i-1,j+0)
      nn.ij[4,] <- c(i+1,j+0)
    }
    # On the right edge of matrix
    if ((1 < i) && (i < L) && (1 < j) && (j == L)) { 
      nn.ij <- matrix(0,3,2)
      nn.ij[1,] <- c(i+0,j-1)
      nn.ij[2,] <- c(i-1,j+0)
      nn.ij[3,] <- c(i+1,j+0)
    }
    # On the Left edge of matrix
    if ((1 < i) && (i < L) && (1 == j) && (j < L)) { 
      nn.ij <- matrix(0,3,2)
      nn.ij[1,] <- c(i+0,j+1)
      nn.ij[2,] <- c(i-1,j+0)
      nn.ij[3,] <- c(i+1,j+0)
    }
    # On the Top edge of matrix
    if ((1 == i) && (i < L) && (1 < j) && (j < L)) { 
      nn.ij <- matrix(0,3,2)
      nn.ij[1,] <- c(i+0,j+1)
      nn.ij[2,] <- c(i+0,j-1)
      nn.ij[3,] <- c(i+1,j+0)
    }
    # On the Bottom edge of matrix
    if ((1 < i) && (i == L) && (1 < j) && (j < L)) { 
      nn.ij <- matrix(0,3,2)
      nn.ij[1,] <- c(i+0,j+1)
      nn.ij[2,] <- c(i+0,j-1)
      nn.ij[3,] <- c(i-1,j+0)
    }
    # On the Top-Left Corner of matrix
    if ((1 == i) && (i < L) && (1 == j) && (j < L)) { 
      nn.ij <- matrix(0,2,2)
      nn.ij[1,] <- c(i+0,j+1)
      nn.ij[2,] <- c(i+1,j+0)
    }
    # On the Top-Right Corner of matrix
    if ((1 == i) && (i < L) && (1 < j) && (j == L)) { 
      nn.ij <- matrix(0,2,2)
      nn.ij[1,] <- c(i+0,j-1)
      nn.ij[2,] <- c(i+1,j+0)
    }
    # On the Bottom-Left Corner of matrix
    if ((1 < i) && (i == L) && (1 == j) && (j < L)) { 
      nn.ij <- matrix(0,2,2)
      nn.ij[1,] <- c(i+0,j+1)
      nn.ij[2,] <- c(i-1,j+0)
    }
    # On the Bottom-Right Corner of matrix
    if ((1 < i) && (i == L) && (1 < j) && (j == L)) { 
      nn.ij <- matrix(0,2,2)
      nn.ij[1,] <- c(i+0,j-1)
      nn.ij[2,] <- c(i-1,j+0)
    }
    riv.nn[[i,j]] <- nn.ij
  }	
}


tic <- 0
coords <- matrix(0,2,L*L)
for (i in 1:L) {
  for (j in 1:L) {
    tic <- tic + 1
    coords[,tic] <- c(i,j)
  }	
}

##################################
# Construct RANDOM SPANNING TREE #
# Use a random walk algorithm to explore a grid and
# create a fully connected spanning tree.
# Serves as a starting condition for OCN sim annealing alg
riv.mat <- matrix(0,L,L)
adj.m <- matrix(0,L*L,L*L)
g <- graph.adjacency(adj.m)
node <- 1
riv.mat[node] <- 1
node.coord <- which(riv.mat == 1, arr.ind=T)
#span.tree <- node
endtree.tic <- 0
stop.tic <- 0
empty.nodes <- rep(1,10)
while ((stop.tic == 0) && (length(empty.nodes) > 0)) {
  #Reset remove local tree order
  remove.local.tree <- FALSE
  advance.endtree.tic <- FALSE
  #Initiate local.tree (if it does not exist)
  if (exists("local.tree")==FALSE) {
    local.tree <- node
  }
  
  #If on first tree, define spanning tree as local tree
  if (endtree.tic == 0) {span.tree <- local.tree}
  
  #Find nearest neighbors to current node
  node.nn <- riv.nn[[node.coord[1],node.coord[2]]]
  #eliminate backtrack option (delete prev.node option)
  if (exists("prev.node")) {
    del.nn <- which(apply(node.nn, 1, function(x) all(x == prev.node.coord)) == TRUE)
    if (length(del.nn) > 0) {node.nn <- node.nn[-del.nn,]}
  }
  # choose random nn to attach to & define: connect.node, connect.node.coords
  if (is.matrix(node.nn) == TRUE) {
    connect.node.coords <- node.nn[sample(seq(1,length(node.nn[,1])),1),]
  } else {connect.node.coords <- node.nn}
  check.mat <- matrix(0,L,L)
  check.mat[connect.node.coords[1],connect.node.coords[2]] <- 1
  connect.node <- which(check.mat == 1)
  
  #Build 'Proposed graph' and check if connection is acceptable
  prop.adj.m <- adj.m
  prop.adj.m[node,connect.node] <- 1
  prop.adj.m[connect.node,node] <- 1
  proposed.g <- graph.adjacency(prop.adj.m,mode="undirected")
  #connect.degree <- degree(proposed.g)[connect.node]
  
  #If on FIRST local tree and proposed node hits it, then:
  if ((connect.node %in% local.tree) && (endtree.tic == 0)) {
    #first tic is the end of the first spanning tree RW
    advance.endtree.tic <- TRUE
    #Reset local tree (later)
    remove.local.tree <- TRUE
    #Find new starting node but do not accept connect.node
    empty.nodes <- which(degree(g) == 0)
    if (length(empty.nodes) > 0) {
      node <- sample(empty.nodes,1)
      check.mat <- matrix(0,L,L)
      check.mat[node] <- 1
      node.coord <- which(check.mat == 1, arr.ind=T)
      #Now eliminate previous node (because we are restarting a new search)
      rm(prev.node)
      rm(prev.node.coord)
    }
  }
  
  #If on >1 local tree, and proposed node hits the SPANNING TREE
  if ((connect.node %in% span.tree) && (endtree.tic > 0)) {
    #Find new starting node and ACCEPT connect.node
    #This records how many local trees are ended.
    #first tic is the end of the first spanning tree RW
    advance.end.tree.tic <- TRUE
    #Update Graph and then find NEW STARTING POINT
    #Accept the connection
    adj.m <- prop.adj.m
    g <- graph.adjacency(adj.m,mode="undirected")
    #collapse local tree and spanning tree
    span.tree <- unique(c(span.tree,local.tree))
    #Reset local tree (later)
    remove.local.tree <- TRUE
    #Find new starting node somewhere along the spanning tree
    empty.nodes <- which(degree(g) == 0)
    if (length(empty.nodes) > 0) {
      node <- sample(empty.nodes,1)
      check.mat <- matrix(0,L,L)
      check.mat[node] <- 1
      node.coord <- which(check.mat == 1, arr.ind=T)
      #Now eliminate previous node (because we are restarting a new search)
      rm(prev.node)
      rm(prev.node.coord)
    }
  } 
  
  #If the proposed node DOES NOT HIT the spanning tree, but DOES HIT the local tree
#   if ((!prop.node %in% span.tree) && (prop.node %in% local.tree) && (endtree.tic > 0)) {
#     #Redraw and do over ((this may not require a separate statement))
#   }
  
  #If the proposed node is NOT in the spanning tree, or the local tree (open case)
  if ((!connect.node %in% span.tree) && (!connect.node %in% local.tree)) {
    #Continue growing and ACCEPT connect.node
    adj.m <- prop.adj.m
    g <- graph.adjacency(adj.m,mode="undirected")
    #span.tree <- cbind(span.tree,connect.node)
    local.tree <- c(local.tree,connect.node)
    #Identify previous node so that it is not 'reselected'
    prev.node <- node
    check.mat <- matrix(0,L,L)
    check.mat[prev.node] <- 1
    prev.node.coord <- which(check.mat == 1, arr.ind=T)
    #The new node is now the 'connected node'
    node <- connect.node
    check.mat <- matrix(0,L,L)
    check.mat[node] <- 1
    node.coord <- which(check.mat == 1, arr.ind=T)
    no.components <- no.clusters(g)
    print(no.components)
  }
  
  # If there are NO NEAREST NEIGHBORS to choose from (DEAD END), select a point on local tree, and branch off.
  #Find matrix locations for each nearest neighbor
  if (is.matrix(node.nn)) {
    nn.fill <- all(apply(node.nn, 1, function(x) {
          check.mat <- matrix(0,L,L)
          check.mat[x[1],x[2]] <-1
          nn.check <- which(check.mat==1)
      #nn.check <- x[2]*L - (L-x[1]) #finds the matrix positon!
      ((nn.check %in% local.tree) || (nn.check %in% span.tree))
    }))
  } else {
    check.mat <- matrix(0,L,L)
    check.mat[node.nn[1],node.nn[2]] <-1
    nn.check <- which(check.mat==1)
    nn.fill <- ((nn.check %in% local.tree) || (nn.check %in% span.tree))
  }
  if (nn.fill) {
    #Find new starting node somewhere along the spanning tree
    #empty.nodes <- which(degree(g) == 0)
    #if (length(empty.nodes) > 0) {
    node <- sample(local.tree,1)
    check.mat <- matrix(0,L,L)
    check.mat[node] <- 1
    node.coord <- which(check.mat == 1, arr.ind=T)
    #Now eliminate previous node (because we are restarting a new search)
    rm(prev.node)
    rm(prev.node.coord)
    #}
  }
#   else {
#     node <- node.nn
#     check.mat <- matrix(0,L,L)
#     check.mat[node] <- 1
#     node.coord <- which(check.mat == 1, arr.ind=T)
#     rm(prev.node)
#     rm(prev.node.coord)
#   }

  #Remove local tree?
  if (remove.local.tree) {rm(local.tree)}
  if (advance.endtree.tic) {endtree.tic <- endtree.tic + 1}
  #Modify stop.tic?
  
  if (no.components == 1) {
    stop.tic <- 1
  } else {
    stop.tic <- 0
  }
  
  
}


plot(g,layout=t(coords),vertex.size=1,vertex.label=NA,edge.arrow.mode=0)





