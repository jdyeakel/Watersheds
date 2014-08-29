EdenGrowth_func3 <- function(size,boundary) {
  L <- size*2
  #Set Seed
  f.eden.tree <- matrix(0,size,2)
  f.eden.tree[1,] <- c(1,1)
  f.eden.tree[2,] <- c(1,2)
  link.tree <- matrix(0,(size-1),2)
  link.tree[1,] <- c(1,2)
  num.node <- 2
  node.matrix <- matrix(0,L*L,L*L); node.matrix[1,1] <- 1; node.matrix[1,2] <- 2
  nn.matrix <- matrix(c(0,1,1,1,0,-1,-1,-1,1,1,0,-1,-1,-1,0,1),8,2)
  
  #Initiate Growth Process
  while (length(which(f.eden.tree[,1]>0)) < size) {
    eden.tree <- f.eden.tree[which(f.eden.tree[,1]>0),]
    set.seed <- matrix(0,100,2)
    max.bound <- max(eden.tree) + boundary
    set.seed.draw <- t(apply(set.seed,1,function(x){sample(seq(1,max.bound),2)}))
    #Select Coords not in eden.tree
    set.seed.test <- apply(set.seed.draw,1,function(x){
      all((x[1] %in% eden.tree[,1]) && (x[2] %in% eden.tree[,2]))
      })
    seed <- set.seed.draw[sample(which(set.seed.test == FALSE),1),]
    
    #Begin Random Walk until it hits eden.tree or leaves the Zone
    check <- TRUE
    while (check) {
      #Advance RW
      seed.new <- seed + sample(seq(-1,1,1),2,replace=TRUE)
      euclid.dist <- as.numeric(apply(eden.tree,1,function(x){dist(rbind(x,seed.new))}))
      min.euclid.dist <- min(euclid.dist)
      # Are you within bounds?
      if ((min.euclid.dist < (boundary+10)) && (all(seed.new >= 1))) {
        
        #If the proposed seed is not part of the tree
        in.tree.test <- any(apply(eden.tree,1,function(x){all(seed.new[1]==x[1],seed.new[2]==x[2])}))
        
        if (in.tree.test == FALSE) {
          #hit.test <- apply(eden.tree,1,function(x){dist(rbind(x,seed.new))})
          hit.test <- euclid.dist
          #Find the node with  0 < dist < sqrt(2)
          #hit.node <- which(apply(cbind(hit.test <= sqrt(2),hit.test > 0),1,all))
          # <= sqrt(2) includes diagonal attachments
          # < sqrt(2) does not allow diagonal attachments
          
          #Evaluate the nearest neighbors of the seed.new
          nn.seed <- t(apply(nn.matrix,1,function(x){c(seed.new[1]+x[1],seed.new[2]+x[2])}))
          #Find nodes in the nearest neighbors
          hit.node <- which(apply(nn.seed,1,function(x){node.matrix[x[1],x[2]]})>0)
          
          #If you hit the tree
          if (length(hit.node) > 0) {
            #Name new node
            num.node <- num.node + 1
            check <- FALSE
            #Add seed.pre to eden.tree
            f.eden.tree[num.node,] <- seed.new
            #Note the link
            link.tree[(num.node-1),] <- c(node.matrix[nn.seed[hit.node[1],1],nn.seed[hit.node[1],2]],num.node)
            #link.tree <- rbind(link.tree,c(min(hit.node),num.node))
            #Record the node on the matrix
            node.matrix[seed.new[1],seed.new[2]] <- num.node
            
            print(length(which(eden.tree[,1]>0)))
          } else {
            #If you do not hit the tree
            seed <- seed.new
          }
        } else {check <- FALSE} #Pick New Seed and start over
      } else {check <- FALSE} #Pick New Seed and start over
    }
  }
  g <- graph.edgelist(link.tree)
  # Calculate area
  # First we need to orient all directed edges to flow towards the confluence
  # Identify Confluence and tributaries
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
  return(list('eden.tree'=eden.tree,'edge.list'=link.tree,'area'=r.area))  
}

# g <- graph.edgelist(link.tree)
# plot(f.eden.tree,pch=15,cex=1,xlim=c(0,10),ylim=c(0,10))
# plot(g, layout=eden.tree,
#      vertex.color="lightblue",edge.color="darkgrey",edge.arrow.mode=0,edge.arrow.size=0.3,
#      vertex.label=NA,vertex.size=1)

