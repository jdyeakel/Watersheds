
library(RColorBrewer)
library(igraph)
library(intergraph)
library(Rcpp)
sourceCpp('/Users/justinyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/src/ESL_meta2.cpp')
sourceCpp('/Users/justinyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/src/ESL_meta.cpp')

#Build nondirected (p=0) and directed graphs (p=1) to start
#Length of simulation
t.term <- 5000
#Number of nodes
n <- 100
#Number of links
l <- 200
repetitions <- 10
p.single <- 1

graphs <- list()
dir.graphs <- list()

# #Build initial nets first
# #Full
# full.adj.m <- matrix(1,n,n); diag(full.adj.m) <- 0
# full.net <- graph.adjacency(full.adj.m)
# graphs[[1]] <- full.net

#Lattice
lattice.net <- graph.lattice(dimvector=c(sqrt(n),sqrt(n)),directed=TRUE,mutual=TRUE)
graphs[[1]] <- lattice.net
edgel <- get.edgelist(lattice.net)
edgel.new <- edgel
skip <- 0
tic <- 0
del.link <- numeric()
for (k in 1:length(edgel[,1])) {
  if ((k %in% skip) == FALSE) {
    link <- edgel[k,]
    pos1 <- which(apply(edgel,1,function(x){(x[1]==link[1] && x[2]==link[2])})==TRUE)
    pos2 <- which(apply(edgel,1,function(x){(x[1]==link[2] && x[2]==link[1])})==TRUE)
    skip <- c(skip,pos1,pos2)
    draw.single <- runif(1)
    if (draw.single < p.single) {
      tic <- tic + 1
      del.link[tic] <- sample(c(pos1,pos2),1)
    }
  }
}
if (length(del.link > 0)) {
  edgel.new <- edgel.new[-del.link,]
}
lattice.net.dir <- graph.edgelist(edgel.new)
adj.m <- get.adjacency(lattice.net.dir)
#isSymmetric(as.matrix(adj.m))
dir.graphs[[1]] <- lattice.net.dir

#Random net
random.net <- erdos.renyi.game(n,l,type="gnm",directed=FALSE,loops=FALSE)
nc <- no.clusters(random.net)
while (nc > 1) {
  random.net <- erdos.renyi.game(n,l,type="gnm",directed=FALSE,loops=FALSE)
  nc <- no.clusters(random.net)
}
graphs[[2]] <- random.net
edgelist_in <- get.edgelist(random.net)
edgelist_out <- matrix(0,length(edgelist_in[,1]),2)
edgelist_out[,1] <- edgelist_in[,2]; edgelist_out[,2] <- edgelist_in[,1]
edgelist_dir <- rbind(edgelist_in,edgelist_out)
random.net.dir <- graph.edgelist(edgelist_dir,directed=TRUE) #FULLY BIDIRECTIONAL GRAPH
edgel <- get.edgelist(random.net.dir)
edgel.new <- edgel
skip <- 0
tic <- 0
del.link <- numeric()
for (k in 1:length(edgel[,1])) {
  if ((k %in% skip) == FALSE) {
    link <- edgel[k,]
    pos1 <- which(apply(edgel,1,function(x){(x[1]==link[1] && x[2]==link[2])})==TRUE)
    pos2 <- which(apply(edgel,1,function(x){(x[1]==link[2] && x[2]==link[1])})==TRUE)
    skip <- c(skip,pos1,pos2)
    draw.single <- runif(1)
    if (draw.single < p.single) {
      tic <- tic + 1
      del.link[tic] <- sample(c(pos1,pos2),1)
    }
  }
}
if (length(del.link > 0)) {
  edgel.new <- edgel.new[-del.link,]
}
random.net.dir <- graph.edgelist(edgel.new)
adj.m <- get.adjacency(random.net.dir)
#isSymmetric(as.matrix(adj.m))
dir.graphs[[2]] <- random.net.dir

#Scale Free
#Scale-Free ::
degs <- sample(1:n, n, replace=TRUE, prob=(1:n)^-1.8)
edgec <- 1
while (edgec != 2*l){
  degs <- sample(1:n, n, replace=TRUE, prob=(1:n)^-1.8)
  edgec <- sum(degs)
}
if (sum(degs) %% 2 != 0) { degs[1] <- degs[1] + 1 }
sf.net <- degree.sequence.game(out.deg=degs, method="vl")
graphs[[3]] <- sf.net
edgelist_in <- get.edgelist(sf.net)
edgelist_out <- matrix(0,length(edgelist_in[,1]),2)
edgelist_out[,1] <- edgelist_in[,2]; edgelist_out[,2] <- edgelist_in[,1]
edgelist_dir <- rbind(edgelist_in,edgelist_out)
sf.net.dir <- graph.edgelist(edgelist_dir,directed=TRUE) #FULLY BIDIRECTIONAL GRAPH
#all(degree(sf.net) == degs)
#adj.m <- get.adjacency(sf.net)
#isSymmetric(adj.m)
edgel <- get.edgelist(sf.net.dir)
edgel.new <- edgel
skip <- 0
tic <- 0
del.link <- numeric()
for (k in 1:length(edgel[,1])) {
  if ((k %in% skip) == FALSE) {
    link <- edgel[k,]
    pos1 <- which(apply(edgel,1,function(x){(x[1]==link[1] && x[2]==link[2])})==TRUE)
    pos2 <- which(apply(edgel,1,function(x){(x[1]==link[2] && x[2]==link[1])})==TRUE)
    skip <- c(skip,pos1,pos2)
    draw.single <- runif(1)
    if (draw.single < p.single) {
      tic <- tic + 1
      del.link[tic] <- sample(c(pos1,pos2),1)
    }
  }
}
if (length(del.link > 0)) {
  edgel.new <- edgel.new[-del.link,]
}
sf.net.dir <- graph.edgelist(edgel.new)
adj.m <- get.adjacency(sf.net.dir)
dir.graphs[[3]] <- sf.net.dir

MRG <- list()
MRG.dir <- list()

#Loop through Graphs
for (i in 1:3) {
  
  MR <- list()
  MR.dir <- list()
  
  #Loop through Rescue Effects
  m.seq <- seq(0,1,0.1)
  for (resc in 1:length(m.seq)) {
    
    m <- m.seq[resc]
    
    print(paste("Graphs = ",i,"/3 :: Rescue Effect = ",m,sep=""))
    
    #Establish adjacency matrix
    if (i==1) {
      adj.m <- get.adjacency(lattice.net)
      adj.m.dir <- get.adjacency(lattice.net.dir)
    }
    if (i==2) {
      adj.m <- get.adjacency(random.net)
      adj.m.dir <- get.adjacency(random.net.dir)
    }
    if (i==3) {
      adj.m <- get.adjacency(sf.net)
      adj.m.dir <- get.adjacency(sf.net.dir)
    }
    
    #Make list of nearest neighbors for each node
    nn <- list()
    for (j in 1:n) {
      arow <- adj.m[j,]
      nn[[j]] <- which(arow==1)
    }
    nn_cpp <- lapply(nn,function(x){x-1})
    
    nn.dir <- list()
    for (j in 1:n) {
      arow <- adj.m.dir[j,]
      nn.dir[[j]] <- which(arow==1)
    }
    nn_cpp_dir <- lapply(nn.dir,function(x){x-1})
    
    #Extinction sequence (prob(ext) for Large patch)
    ext.seq <- seq(0,0.8,0.04) #Change to seq(0,0.4,0.02) for future runs... or seq(0,0.4,0.01) for the final product
    
    m.e <- matrix(0,length(ext.seq),repetitions)
    m.s <- matrix(0,length(ext.seq),repetitions)
    m.l <- matrix(0,length(ext.seq),repetitions)
    
    m.e.dir <- matrix(0,length(ext.seq),repetitions)
    m.s.dir <- matrix(0,length(ext.seq),repetitions)
    m.l.dir <- matrix(0,length(ext.seq),repetitions)
    
    #Iterate over different values of large-patch extinction probabilities
    #NOTE :: These iterations are now all done in C++
    
    #for (j in 1:length(ext.seq)) {
    #for (rep in 1:replicates) {
    
    #Define intial conditions
    X <- matrix(0,n,t.term)
    Meta <- matrix(0,3,t.term)
    #Initiate starting vector
    initial <- c(n - 2*round(n/3,0),round(n/3,0),round(n/3,0))
    X[,1] <- c(rep(0,initial[1]),rep(1,initial[2]),rep(2,initial[3]))
    Meta[,1] <- initial
    
    #Values of constants
    #Probabilities
    c <- 0.2
    r <- 0.2
    #el <- ext.seq[j]
    #es <- 2*el
    
    #Perform Simulation
    ESL.out <- ESL_meta2(n, t.term, repetitions, X, Meta, nn_cpp, c, r, m, ext.seq)
    ESL.out.dir <- ESL_meta2(n, t.term, repetitions, X, Meta, nn_cpp_dir, c, r, m, ext.seq)
    
    #Save summary statistics
    #These are matrices of mean proportions with ext.seq=rows and repetitions=columns...
    for (j in 1:length(ext.seq)) {
      m.e[j,] <- ESL.out[[1]][[j]][1,]/n
      m.s[j,] <- ESL.out[[1]][[j]][2,]/n
      m.l[j,] <- ESL.out[[1]][[j]][3,]/n
    } 
    for (j in 1:length(ext.seq)) {
      m.e.dir[j,] <- ESL.out.dir[[1]][[j]][1,]/n
      m.s.dir[j,] <- ESL.out.dir[[1]][[j]][2,]/n
      m.l.dir[j,] <- ESL.out.dir[[1]][[j]][3,]/n
    } 
    
    M <- list()
    M.dir <- list()
    
    M[[1]] <- m.e
    M[[2]] <- m.s
    M[[3]] <- m.l
    
    M.dir[[1]] <- m.e.dir
    M.dir[[2]] <- m.s.dir
    M.dir[[3]] <- m.l.dir
    
    MR[[resc]] <- M
    
    MR.dir[[resc]] <- M.dir
    
  } # end resc
  
  MRG[[i]] <- MR
  MRG.dir[[i]] <- MR.dir
  
} # end i
image_filename <- paste("/Users/justinyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/Results/ESL_DirectedDifference2.RData",sep="")
save.image(image_filename)
