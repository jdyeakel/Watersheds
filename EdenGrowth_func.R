EdenGrowth_func <- function(size) {
  L <- size*2
  #Set Seed
  eden.tree <- rbind(c(1,1),c(1,2))
  link.tree <- c(1,2)
  num.node <- 2
  
  #Initiate Growth Process
  while (length(eden.tree[,1]) < size) {
    set.seed <- matrix(0,100,2)
    max.bound <- max(eden.tree) + 10
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
      euclid.dist <- min(as.numeric(apply(eden.tree,1,function(x){dist(rbind(x,seed.new))})))
      # Are you within bounds?
      if ((euclid.dist < 15) && (all(seed.new > 1))) {
        
        #       hit.test <- apply(eden.tree,1,function(x){
        #         all((seed.new[1] %in% x[1]) && (seed.new[2] %in% x[2]))
        #       })
        #If the proposed seed is not part of the tree
        in.tree.test <- any(apply(eden.tree,1,function(x){all(seed.new[1]==x[1],seed.new[2]==x[2])}))
        if (in.tree.test == FALSE) {
          hit.test <- apply(eden.tree,1,function(x){dist(rbind(x,seed.new))})
          #Find the node with  0 < dist < sqrt(2)
          hit.node <- which(apply(cbind(hit.test < sqrt(2),hit.test > 0),1,function(x){all(x)}))
          # <= sqrt(2) includes diagonal attachments
          # < sqrt(2) does not allow diagonal attachments
          
          #If you hit the tree
          if (length(hit.node) > 0) {
            #Name new node
            num.node <- num.node + 1
            check <- FALSE
            #Add seed.pre to eden.tree
            eden.tree <- rbind(eden.tree,seed.new)
            #Note the link
            link.tree <- rbind(link.tree,c(min(hit.node),num.node))
            
            #           #Remove duplicates
            #           remove <- which(duplicated(eden.tree))
            #           if (length(remove)>0){
            #             eden.tree <- eden.tree[-remove,]
            #             link.tree <- link.tree[-remove,]
            #           }
            
            print(length(link.tree[,1]))
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
# plot(eden.tree,pch=15,cex=1,xlim=c(0,50),ylim=c(0,50))
# plot(g, layout=eden.tree,
#      vertex.color="lightblue",edge.color="darkgrey",edge.arrow.mode=0,edge.arrow.size=0.3,
#      vertex.label=NA,vertex.size=log(r.area +1))

