RivOCN_func <- function(L,tic_max) {
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
  
  riv.mat <- matrix(0,L,L)
  adj.m <- matrix(0,L*L,L*L)
  #g <- graph.adjacency(adj.m)
  node <- 1
  riv.mat[node] <- 1
  node.coord <- which(riv.mat == 1, arr.ind=T)
  #empty.nodes <- which(degree(g)==0)
  tic <- 0
  span.tree <- 0
  while (length(span.tree) < (L*L)) { #empty.nodes > 0
    span.tree <- unique(span.tree)
    tic <- tic + 1
    if (tic == 1) {span.tree <- node} else {span.tree <- c(span.tree,node)}
    #Find nearest neighbors to node
    node.nn.coord <- riv.nn[[node.coord[1],node.coord[2]]]
    node.nn <- apply(node.nn.coord, 1, function(x) {
      check.mat <- matrix(0,L,L)
      check.mat[x[1],x[2]] <-1
      which(check.mat==1)
    })
    #Find neighboring nodes that are OPEN
    open.nodes <- node.nn[which(!node.nn %in% span.tree == TRUE)]
    #IF THERE ARE OPEN NODES
    #Randomly select an open node
    if (length(open.nodes) > 0) {
      random.draw <- round(runif(1,min=1,max=length(open.nodes)),0)
      connect.node <- open.nodes[random.draw]
      
      #Now establish link between current node and connect.node
      #Update the spanning tree
      adj.m[node,connect.node] <- 1
      adj.m[connect.node,node] <- 1
      #g <- graph.adjacency(adj.m,mode="undirected")
      
      #Now label the connect.node the current node
      node <- connect.node
      check.mat <- matrix(0,L,L)
      check.mat[node] <- 1
      node.coord <- which(check.mat == 1, arr.ind=T)
    } else {
      
      #IF THERE ARE NO OPEN NODES
      #Randomly select a node from the span.tree (with open neighbors) and try again.
      #Find nodes on span.tree with OPEN NEIGHBORS
      #Make a list of neighbors for each node on span.tree
      nn.span <- list(length(span.tree))
      for (i in 1:length(span.tree)) {
        span.node <- span.tree[i]
        nn.span[[i]] <- riv.nn[[span.node]]
      }
      
      #Determine which neighbors are open:
      nn.open.test <- lapply(nn.span,function(x){
        apply(x,1,function(y){
          check.mat <- matrix(0,L,L)
          check.mat[y[1],y[2]] <- 1
          nn.check <- which(check.mat==1)
          nn.check %in% span.tree
        })
      })
      span.tree.open <- span.tree[which(unlist(lapply(
        nn.open.test,function(x){any(x==FALSE)}))==TRUE)]
      #Min=2 does not allow drawing the confluence node as a starting point
      random.draw <- round(runif(1,min=2,max=length(span.tree.open)),0)
      
      #Define new starting position (do not draw links on this step)
      node <- span.tree.open[random.draw]
      check.mat <- matrix(0,L,L)
      check.mat[node] <- 1
      node.coord <- which(check.mat == 1, arr.ind=T)
    }
    
    #Current alg. stops when span.tree fills nodes.
    #This way, we don't have to use igraph (which takes a long time)
    #empty.nodes <- which(degree(g)==0)
  }
  g <- graph.adjacency(adj.m,mode="undirected")
  ###########################################
  # BEGIN OPTIMAL CHANNEL NETWORK ALGORITHM #
  ###########################################
  
  # First we need to orient all directed edges to flow towards the confluence
  # Identify Confluence and tributaries
  tribs <- which(degree(g)==1)
  conf <- tribs[1]
  tribs <- tribs[2:length(tribs)]
  lt <- length(tribs)
  #Build (accurate) directed net
  dadj.m <- matrix(0,L*L,L*L)
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
  r.area <- numeric(L*L)
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
  early.r.area <- r.area
  H <- sum(r.area^gamma)		
  early.dg <- dg
  
  
#     plot(dg,layout=t(coords),vertex.size=log(early.r.area+1)/2,vertex.label=NA,
#        edge.arrow.mode=1,edge.arrow.size=0.5)
#   
  
  
  # Build Late River Network with Simulated Annealing Algorithm
  
  tic <- 0
  tic_max <- tic_max
  r.H <- numeric()
  #minutes <- time
  #pt.max <- minutes*60 # how long to run simulation in seconds
  #pt1 <- proc.time()
  #pt.elapsed <- 0
  #while (pt.elapsed < pt.max) {
  while (tic < tic_max) {
    tic <- tic + 1
    new.adj.m <- dadj.m
    
    edge.list <- get.edgelist(dg)
    #Delete links connected to conluence (to retain confluence)
    to.sample <- intersect(which(edge.list[,1] != conf),which(edge.list[,2] != conf))
    r <- sample(to.sample,1)
    change.link <- edge.list[r,]
    #Find coordinates of ego node
    ego <- change.link[1]
    change.link.m <- matrix(0,L,L); change.link.m[ego] <- 1
    node.i <- as.numeric(which(change.link.m == 1, arr.ind = TRUE))
    
    #Find nearest nieghbors
    nn.i <- riv.nn[[node.i[1],node.i[2]]]
    # delete nearest neighbors that have links to ego
    lnn <- length(nn.i[,1])
    tick <- 0; nn.del <- numeric()
    for (k in 1:lnn) {
      conn <- riv.nn[[node.i[1],node.i[2]]][k,]
      link.m <- matrix(0,L,L); link.m[conn[1],conn[2]] <- 1
      link <- which(link.m == 1)
      
      if ((dadj.m[ego,link] == 1) || (dadj.m[link,ego] == 1)) {
        tick <- tick + 1
        nn.del[tick] <- k
      }
    }
    if (length(nn.del) > 0) {
      nn.i <- nn.i[-nn.del,]
    }
    
    # If there are nearest neighbors, perform next step. If not, continue while loop
    if (length(nn.i) > 0) {
      # choose random nearest neighbor
      if (is.matrix(nn.i) == TRUE) {
        lnn <- length(nn.i[,1])
        seq.nn <- sample(seq(1,lnn),lnn,replace=F)
        test.no.components <- numeric(lnn)
        new.link <- numeric(lnn)
        ticn <- 0
        for (n in seq.nn) {
          ticn <- ticn + 1
          new.link.m <- matrix(0,L,L); new.link.m[nn.i[n,1],nn.i[n,2]] <- 1
          new.link[ticn] <- which(new.link.m == 1)
          # Make sure the change does not increase components or make loops
          ### Find minimum components!
          # Test change
          test.adj.m <- dadj.m
          test.adj.m[change.link[1],change.link[2]] <- 0
          test.adj.m[ego,new.link[ticn]] <- 1
          test.g <- graph.adjacency(test.adj.m)
          #Convert test.g to undirected graph
          #test.g <- as.undirected(test.g)
          test.no.components[ticn] <- no.clusters(test.g)
        }	
      } else {
        new.link.m <- matrix(0,L,L); new.link.m[nn.i[1],nn.i[2]] <- 1
        new.link <- which(new.link.m == 1)
        # Make sure the change does not increase components or make loops
        ### Find minimum components!
        # Test change
        test.adj.m <- dadj.m
        test.adj.m[change.link[1],change.link[2]] <- 0
        test.adj.m[ego,new.link] <- 1
        test.g <- graph.adjacency(test.adj.m)
        #Convert test.g to undirected graph
        #test.g <- as.undirected(test.g)
        test.no.components <- no.clusters(test.g)
      }		
      #which(test.no.components == 1)
      potential.moves <- which(test.no.components == 1)
      
      # If the change does not increase components, then test RIVER AREA
      if (length(potential.moves) > 0) {
        
        #Choose applicable link
        rl <- sample(potential.moves,1)
        #Eliminate old link
        new.adj.m[change.link[1],change.link[2]] <- 0
        # Build new link in proposed adjacency matrix
        new.adj.m[ego,new.link[rl]] <- 1
        # Build new proposed graph
        new.g <- graph.adjacency(new.adj.m)
        
        # Establish New River Flow
        # Identify Confluence and tributaries
        tribs <- which(degree(new.g)==1)
        tribs <- tribs[-1]
        lt <- length(tribs)
        #Build (accurate) directed net
        new.dadj.m <- matrix(0,L*L,L*L)
        for (n in 1:lt) {
          trib.i <- tribs[n]
          #find the path from trib to confluence (w.o.r. to flow)
          path <- get.shortest.paths(new.g,trib.i,conf,mode="all")[[1]]
          lp <- length(path)
          #Build the directed adjacency matrix
          for (m in 1:(lp-1)) {
            to.node <- path[m]
            from.node <- path[m+1]
            new.dadj.m[from.node,to.node] <- 1
          }	
        }
        new.dg <- graph.adjacency(new.dadj.m)
        
        if (no.clusters(new.dg) == 1) {
          
          # #Assign River Area
          # new.r.area <- numeric(L*L) + 1
          # for (n in 1:(L*L)) {
          # if (all(n != tribs)) {
          # #Find all paths flowing into node i
          # paths.upriver <- get.shortest.paths(new.dg,n,mode="out")
          # paths.upriver <- do.call(c,paths.upriver)
          # #What are the unique paths?
          # unique.upriver <- unique(paths.upriver)
          # # Do not need '+1' because node i is included in this list
          # new.r.area[n] <- length(unique.upriver)
          # }
          # }
          
          # Assign River Area (NEW)
          new.r.area <- numeric(L*L) + 1
          paths <- shortest.paths(new.dg,conf,mode="out")
          longest.path <- which(paths==max(paths))
          lpath <- paths[longest.path[1]]
          for (n in seq(lpath,0,-1)) {
            nodes.dist <- which(paths == n)
            num.nodes <- length(nodes.dist)
            for (m in nodes.dist) {
              incoming.nodes <- neighbors(new.dg,m,mode="out")
              if (length(incoming.nodes) > 0) {
                area.neighbors <- new.r.area[incoming.nodes]
                new.r.area[m] <- sum(area.neighbors) + 1
              } else {
                new.r.area[m] <- 1
              }		
            }
          }
          
          # Calculate H(s')
          new.H <- sum(new.r.area^gamma)
          
          #######################
          # SIMULATED ANNEALING #
          #######################
          
          # Cooling Schedule
          T0 <- tic_max
          alpha <- 0.99
          temp <- T0*alpha^tic
          draw <- runif(1,c(0,1))
          prob <- min(exp(-(new.H - H)/temp),1) #constraints max value to 1
          
          # If draw is less than prob, accept regardless of potential
          if (draw < prob) {
            late.r.area <- new.r.area
            H <- new.H
            dadj.m <- new.dadj.m
            dg <- new.dg
            late.dg <- dg
            #print(c(H,no.clusters(dg)))
          } else {
            # Accept new river config if potential is lower
            if (new.H < H) {
              late.r.area <- new.r.area
              H <- new.H
              dadj.m <- new.dadj.m
              dg <- new.dg
              late.dg <- dg
              #print(c(H,no.clusters(dg)))
            }
          }
        } # End if no.clusters == 1
        
      } # End if potential moves > 0
      
    } # End if there are nearest neighbors
    #Record current H
    r.H[tic] <- H
    print.seq <- round(seq(0,tic_max,(tic_max/50)),0)
    if (tic %in% print.seq) {
      print(paste("iteration =",tic))
    }
    #Record elapsed time
    #pt2 <- proc.time()
    #pt <- pt2-pt1
    #pt.elapsed <- pt[3]
  } # End While loop
  riv.net.prop <- list("late.net"=late.dg,
                       "late.area"=late.r.area,
                       "early.net"=early.dg,
                       "early.area"=early.r.area,
                       "coords"=coords,
                       "r.H"=r.H,
                       "neighbors"=riv.nn)
return(riv.net.prop)
} #End function




