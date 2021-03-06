rm(list=c(ls()))

library(RColorBrewer)
library(igraph)
library(intergraph)
library(Rcpp)
sourceCpp('/Users/justinyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/src/ESL_meta2.cpp')
sourceCpp('/Users/justinyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/src/ESL_meta.cpp')

t.term <- 1000
n <- 100
l <- 200
X <- matrix(0,n,t.term)
Meta <- matrix(0,3,t.term)
#Initiate starting vector
initial <- c(n - 2*round(n/3,0),round(n/3,0),round(n/3,0))
X[,1] <- c(rep(0,initial[1]),rep(1,initial[2]),rep(2,initial[3]))
Meta[,1] <- initial

#Full Network :: every node is connected to every other node
adj.m <- matrix(1,n,n); diag(adj.m) <- 0




#Lattice Network
lattice.net <- graph.lattice(dimvector=c(sqrt(n),sqrt(n)))
adj.m <- get.adjacency(lattice.net)
isSymmetric(as.matrix(adj.m))

#Random Network :: nodes are connected randomly




#Exponential Network ::
degs <- sample(1:n, n, replace=TRUE, prob=exp(-0.5*(1:n)))
edgec <- 1
while (edgec != 2*l){
  degs <- sample(1:n, n, replace=TRUE, prob=(1:n)^-1.8)
  edgec <- sum(degs)
}
if (sum(degs) %% 2 != 0) { degs[1] <- degs[1] + 1 }
exp.net <- degree.sequence.game(degs, method="vl")
all(degree(exp.net) == degs)
adj.m <- get.adjacency(exp.net)
isSymmetric(adj.m)

#Watts-Strogratz ::
ws.net <- watts.strogatz.game(dim=1,size=n,nei=2,p=0.1)
adj.m <- get.adjacency(ws.net)
isSymmetric(as.matrix(adj.m))

#Scale-Free ::
degs <- sample(1:n, n, replace=TRUE, prob=(1:n)^-1.8)
edgec <- 1
while (edgec != 2*l){
  degs <- sample(1:n, n, replace=TRUE, prob=(1:n)^-1.8)
  edgec <- sum(degs)
}
if (sum(degs) %% 2 != 0) { degs[1] <- degs[1] + 1 }
sf.net <- degree.sequence.game(degs, method="vl")
all(degree(sf.net) == degs)
adj.m <- get.adjacency(sf.net)
isSymmetric(adj.m)






colors <- brewer.pal(3,"Set1"); trans <- ""
plot(graph.adjacency(adj.m,mode="undirected"),vertex.size=2,vertex.label=NA,vertex.color=colors[2])
###########################
###########################
###########################
###########################
###########################
sourceCpp('/Users/jdyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/src/ESL_meta.cpp')

t.term <- 10000
n <- 100
l <- 200
X <- matrix(0,n,t.term)
Meta <- matrix(0,3,t.term)
#Initiate starting vector
initial <- c(n - 2*round(n/3,0),round(n/3,0),round(n/3,0))
X[,1] <- c(rep(0,initial[1]),rep(1,initial[2]),rep(2,initial[3]))
Meta[,1] <- initial

#THIS CHECKS OUT 2/11/20
#Make list of nearest neighbors for each i
nn <- list()
for (i in 1:n) {
  arow <- adj.m[i,]
  nn[[i]] <- which(arow==1)
}

#Probabilities
c <- 0.2
r <- 0.2
m <- 0.00
el <- 0.05
es <- 2*el

#### C++ Version
nn_cpp = lapply(nn,function(x){x-1})
ESL.out <- ESL_meta(n, t.term, X, Meta, nn_cpp, c, r, m, es, el)

colors <- brewer.pal(3,"Set1"); trans <- ""
plot(ESL.out[[2]][1,]/n,type="l",col=paste(colors[1],trans,sep=""),ylim=c(0,1),xlim=c(0,t.term),lwd=2)
lines(ESL.out[[2]][2,]/n,col=paste(colors[2],trans,sep=""),lwd=2)
lines(ESL.out[[2]][3,]/n,col=paste(colors[3],trans,sep=""),lwd=2)


# Analysis of C++ simulation
# Across different types of networks
#NOTE 2/11/20 - 
sourceCpp('src/ESL_meta2.cpp')
#Length of simulation
t.term <- 5000
#Number of nodes
n <- 100
#Number of links
l <- 200
repetitions <- 10

MG <- list()
graphs <- list()

#Full Network
full.adj.m <- matrix(1,n,n); diag(full.adj.m) <- 0
full.net <- graph.adjacency(full.adj.m)
graphs[[1]] <- full.net

#Lattice Network
lattice.net <- graph.lattice(dimvector=c(sqrt(n),sqrt(n)))
#isSymmetric(as.matrix(adj.m))
graphs[[2]] <- lattice.net

#Random Network
random.net <- erdos.renyi.game(n,l,type="gnm",directed=FALSE,loops=FALSE)
nc <- no.clusters(random.net)
while (nc > 1) {
  random.net <- erdos.renyi.game(n,l,type="gnm",directed=FALSE,loops=FALSE)
  nc <- no.clusters(random.net)
}
graphs[[3]] <- random.net

#Scale Free network
#Scale-Free ::
degs <- sample(1:n, n, replace=TRUE, prob=(1:n)^-1.8)
edgec <- 1
while (edgec != 2*l){
  degs <- sample(1:n, n, replace=TRUE, prob=(1:n)^-1.8)
  edgec <- sum(degs)
}
if (sum(degs) %% 2 != 0) { degs[1] <- degs[1] + 1 }
sf.net <- degree.sequence.game(degs, method="vl")
#all(degree(sf.net) == degs)
#isSymmetric(adj.m)
graphs[[4]] <- sf.net

for (i in 1:4) {
  print(paste("Graph ",i,"/4"),sep="")
  
  net <- graphs[[i]]
  adj.m <- get.adjacency(net)
  
  #Make list of nearest neighbors for each node
  nn <- list()
  for (j in 1:n) {
    arow <- adj.m[j,]
    nn[[j]] <- which(arow==1)
  }
  nn_cpp <- lapply(nn,function(x){x-1})
  
  #Extinction sequence (prob(ext) for Large patch)
  ext.seq <- seq(0,0.4,0.02)
  
  m.e <- matrix(0,length(ext.seq),repetitions)
  m.s <- matrix(0,length(ext.seq),repetitions)
  m.l <- matrix(0,length(ext.seq),repetitions)
  
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
  m <- 0.1
  #el <- ext.seq[j]
  #es <- 2*el
  
  #Perform Simulation
  ESL.out <- ESL_meta2(n, t.term, repetitions, X, Meta, nn_cpp, c, r, m, ext.seq)
  
  #Save summary statistics
  for (j in 1:length(ext.seq)) {
    m.e[j,] <- ESL.out[[1]][[j]][1,]/n
    m.s[j,] <- ESL.out[[1]][[j]][2,]/n
    m.l[j,] <- ESL.out[[1]][[j]][3,]/n
  }  
  
  #} #end r
  #} #end j
  
  M <- list()
  
  M[[1]] <- m.e
  M[[2]] <- m.s
  M[[3]] <- m.l
  
  MG[[i]] <- M
  
}# end i
image_filename <- paste("/Users/justinyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/Results/ESL_m",m,".RData",sep="")
save.image(image_filename)

load("/Users/justinyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/Results/ESL_m0.RData")

pdf_filename <- paste("/Users/justinyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/Results2/ESL_m",m,".pdf",sep="")

colors <- brewer.pal(3,"Set1"); trans <- ""
dev.off()
pdf(file=pdf_filename
    ,width=4,height=8)
op <- par(mfrow = c(4,1),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,3,1,1) + 0.1,
          mgp = c(2, 1, 0))
for (i in 1:4) {
  
  m.e <- MG[[i]][[1]]
  m.s <- MG[[i]][[2]]
  m.l <- MG[[i]][[3]]
  if (i == 4) {
    plot(ext.seq/c,apply(m.e,1,mean),type="l",col=colors[1],xlim=c(0,2),ylim=c(0,1),xlab="e/c",ylab="Proportion")
  } else {plot(ext.seq/c,apply(m.e,1,mean),type="l",col=colors[1],xlim=c(0,2),ylim=c(0,1), xaxt='n', ann=FALSE)}
  polygon(x=c(ext.seq/c,rev(ext.seq/c)),
          y=c(apply(m.e,1,mean)+apply(m.e,1,sd),rev(apply(m.e,1,mean)-apply(m.e,1,sd))),col=paste(colors[1],"50",sep=""),border=NA)
  lines(ext.seq/c,apply(m.s,1,mean),col=colors[2])
  polygon(x=c(ext.seq/c,rev(ext.seq/c)),
          y=c(apply(m.s,1,mean)+apply(m.s,1,sd),rev(apply(m.s,1,mean)-apply(m.s,1,sd))),col=paste(colors[2],"50",sep=""),border=NA)
  lines(ext.seq/c,apply(m.l,1,mean),col=colors[3])
  polygon(x=c(ext.seq/c,rev(ext.seq/c)),
          y=c(apply(m.l,1,mean)+apply(m.l,1,sd),rev(apply(m.l,1,mean)-apply(m.l,1,sd))),col=paste(colors[3],"50",sep=""),border=NA)
  lines(ext.seq/c,apply(m.s,1,mean)+apply(m.l,1,mean),col="black",lty=2)
}
par(op)
dev.off()



#Difference in overall persistence between Lattice and SF networks
load("/Users/justinyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/Results/ESL_m0.RData")
Diffm0 = apply(MG[[4]][[2]] + MG[[4]][[3]],1,mean) - apply(MG[[3]][[2]] + MG[[3]][[3]],1,mean)
load("/Users/justinyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/Results/ESL_m0.1.RData")
Diffm1 = apply(MG[[4]][[2]] + MG[[4]][[3]],1,mean) - apply(MG[[3]][[2]] + MG[[3]][[3]],1,mean)
pdf(file="/Users/justinyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/Results/ESL_Diff.pdf",width=5.5,height=5.0)
plot(ext.seq/c,Diffm0,type="l",lwd=2,col=colors[2],xlim=c(0,2),ylim=c(-0.2,0.7),xlab="e/c",ylab="Difference in Cum. Proportion (SF-Rand)")
lines(ext.seq/c,Diffm1,lwd=2,col=colors[3])
lines(seq(0,2,0.1),rep(0,21),lty=3)
dev.off()




#Lattice Network
n <- 16 #Must be a 'square-rootable number'
#Probability that a single direction edge exists between 2 vertices
p.single <- 1

lattice.net <- graph.lattice(dimvector=c(sqrt(n),sqrt(n)),directed=TRUE,mutual=TRUE)
adj.m <- get.adjacency(lattice.net)
edgel <- get.edgelist(lattice.net)
edgel.new <- edgel
skip <- 0
tic <- 0
del.link <- numeric()
for (i in 1:length(edgel[,1])) {
  if ((i %in% skip) == FALSE) {
    link <- edgel[i,]
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


#####################################
#With distributions of Directed Links
#####################################
rm(list=c(ls()))

#LOOP over values for p.single
sourceCpp('src/ESL_meta2.cpp')



#Length of simulation
t.term <- 5000
#Number of nodes
n <- 100
#Number of links
l <- 200
repetitions <- 10
p.single.vec <- seq(0,1,0.1)

graphs <- list()

# #Build initial nets first
# #Full
# full.adj.m <- matrix(1,n,n); diag(full.adj.m) <- 0
# full.net <- graph.adjacency(full.adj.m)
# graphs[[1]] <- full.net

#Lattice
lattice.net <- graph.lattice(dimvector=c(sqrt(n),sqrt(n)),directed=TRUE,mutual=TRUE)
graphs[[1]] <- lattice.net

#Random net
random.net <- erdos.renyi.game(n,l,type="gnm",directed=FALSE,loops=FALSE)
nc <- no.clusters(random.net)
while (nc > 1) {
  random.net <- erdos.renyi.game(n,l,type="gnm",directed=FALSE,loops=FALSE)
  nc <- no.clusters(random.net)
}
graphs[[2]] <- random.net

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

MGP <- list()

for (ps in 1:length(p.single.vec)) {
  
  p.single <- p.single.vec[ps]
  
  MG <- list()
  dir.graphs <- list()
  for (i in 1:3) {
    
    print(paste("Pr(single link)=",p.single,"/1.0 :: ","Graph ",i,"/4"),sep="")
    
    #Select Network
#     if (i == 1) {
#       #Full Network :: every node is connected to every other node
#       #adj.m <- get.adjacency(full.net)
#       edgel <- get.edgelist(full.net)
#       edgel.new <- edgel
#       skip <- 0
#       tic <- 0
#       del.link <- numeric()
#       for (k in 1:length(edgel[,1])) {
#         if ((k %in% skip) == FALSE) {
#           link <- edgel[k,]
#           pos1 <- which(apply(edgel,1,function(x){(x[1]==link[1] && x[2]==link[2])})==TRUE)
#           pos2 <- which(apply(edgel,1,function(x){(x[1]==link[2] && x[2]==link[1])})==TRUE)
#           skip <- c(skip,pos1,pos2)
#           draw.single <- runif(1)
#           if (draw.single < p.single) {
#             tic <- tic + 1
#             del.link[tic] <- sample(c(pos1,pos2),1)
#           }
#         }
#       }
#       if (length(del.link > 0)) {
#         edgel.new <- edgel.new[-del.link,]
#       }
#       full.net.dir <- graph.edgelist(edgel.new)
#       adj.m <- get.adjacency(full.net.dir)
#       #graphs[[i]] <- full.net.dir
#     }
    if (i == 1) {
      edgel <- get.edgelist(lattice.net)
      edgel.new <- edgel
      skip <- 0
      tic <- 0
      del.link <- numeric()
      for (k in 1:length(edgel[,1])) {
        if ((k %in% skip) == FALSE) {
          link <- edgel[k,]
          #Find positions in edgelist where node A -> B and B -> A
          pos1 <- which(apply(edgel,1,function(x){(x[1]==link[1] && x[2]==link[2])})==TRUE)
          pos2 <- which(apply(edgel,1,function(x){(x[1]==link[2] && x[2]==link[1])})==TRUE)
          #Add these to skip list
          skip <- c(skip,pos1,pos2)
          draw.single <- runif(1)
          if (draw.single < p.single) {
            tic <- tic + 1
            #Chose one of those links to delete (at random)
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
      dir.graphs[[i]] <- lattice.net.dir
    
      # plot(lattice.net,vertex.size=0.25)
      # plot(lattice.net.dir,vertex.size=0.25)
    }
    if (i == 2) {
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
          #Find positions in edgelist where node A -> B and B -> A
          pos1 <- which(apply(edgel,1,function(x){(x[1]==link[1] && x[2]==link[2])})==TRUE)
          pos2 <- which(apply(edgel,1,function(x){(x[1]==link[2] && x[2]==link[1])})==TRUE)
          #Add these to skip list
          skip <- c(skip,pos1,pos2)
          draw.single <- runif(1)
          if (draw.single < p.single) {
            tic <- tic + 1
            #Chose one of those links to delete (at random)
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
      dir.graphs[[i]] <- random.net.dir
      # plot(random.net,vertex.size=0.25)
      # plot(random.net.dir,vertex.size=0.25)
    }
    if (i == 3) {
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
          #Find positions in edgelist where node A -> B and B -> A
          pos1 <- which(apply(edgel,1,function(x){(x[1]==link[1] && x[2]==link[2])})==TRUE)
          pos2 <- which(apply(edgel,1,function(x){(x[1]==link[2] && x[2]==link[1])})==TRUE)
          #Add these to skip list
          skip <- c(skip,pos1,pos2)
          draw.single <- runif(1)
          if (draw.single < p.single) {
            tic <- tic + 1
            #Chose one of those links to delete (at random)
            del.link[tic] <- sample(c(pos1,pos2),1)
          }
        }
      }
      if (length(del.link > 0)) {
        edgel.new <- edgel.new[-del.link,]
      }
      sf.net.dir <- graph.edgelist(edgel.new)
      adj.m <- get.adjacency(sf.net.dir)
      dir.graphs[[i]] <- sf.net.dir
      # plot(sf.net,vertex.size=0.25)
      # plot(sf.net.dir,vertex.size=0.25)
    }
    
    #Make list of nearest neighbors for each node
    nn <- list()
    for (j in 1:n) {
      arow <- adj.m[j,]
      nn[[j]] <- which(arow==1)
    }
    nn_cpp <- lapply(nn,function(x){x-1})
    
    #Extinction sequence (prob(ext) for Large patch)
    ext.seq <- seq(0,0.4,0.02) #Change to seq(0,0.4,0.02) for future runs... or seq(0,0.4,0.01) for the final product
    
    m.e <- matrix(0,length(ext.seq),repetitions)
    m.s <- matrix(0,length(ext.seq),repetitions)
    m.l <- matrix(0,length(ext.seq),repetitions)
    
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
    m <- 0.0
    #el <- ext.seq[j]
    #es <- 2*el
    
    #Perform Simulation
    ESL.out <- ESL_meta2(n, t.term, repetitions, X, Meta, nn_cpp, c, r, m, ext.seq)
    
    #Save summary statistics
    for (j in 1:length(ext.seq)) {
      m.e[j,] <- ESL.out[[1]][[j]][1,]/n
      m.s[j,] <- ESL.out[[1]][[j]][2,]/n
      m.l[j,] <- ESL.out[[1]][[j]][3,]/n
    }  
    
    #} #end r
    #} #end j
    
    M <- list()
    
    M[[1]] <- m.e
    M[[2]] <- m.s
    M[[3]] <- m.l
    
    MG[[i]] <- M
    
  }# end i
  
  
  MGP[[ps]] <- MG
  
  
}# end ps

image_filename <- paste("/Users/justinyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/Results/ESLdir_m",m,".RData",sep="")
save.image(image_filename)

load("/Users/justinyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/Results/ESLdir_m0.RData")

#Exporing Results for 3D plots in Mathematica
name_var <- c("E","S","L")
name_net <- c("Full","Lattice","SF")
for (i in 1:3) { #nets
  for (j in 1:3) { #ESL
    #Assemble table
    tab <- matrix(0,length(ext.seq)*length(p.single.vec),3)
    tic <- 0
    for (k in 1:length(p.single.vec)) {
      temp_tab <- apply(MGP[[k]][[i]][[j]],1,mean)
      tab[seq(length(temp_tab)*tic+1,length(temp_tab)*(tic+1),1),3] <- temp_tab
      tab[seq(length(temp_tab)*tic+1,length(temp_tab)*(tic+1),1),2] <- p.single.vec[k]
      tic <- tic + 1
    }
    tab[,1] <- rep(ext.seq/c,length(p.single.vec))
    filename <- paste("/Users/justinyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/Results/dir_text/dir_m",m,"_",name_net[i],"_",name_var[j],".csv",sep="")
    write.table(tab,file=filename,col.names=FALSE,row.names=FALSE,sep=",")
  }
}

# Contour plot of the effect of directedness on persistence
pdf_filename <- paste("/Users/justinyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/Results2/ESLdir_m",m,".pdf",sep="")
dev.off()
pdf(file=pdf_filename,width=3.5,height=6)
op <- par(mfrow = c(3,1),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,3,1,1) + 0.1,
          mgp = c(2, 1, 0))
colors <- brewer.pal(11,"Spectral")
for (i in 1:3) { # Loop over different networks
  for (ps in 1:length(p.single.vec)) {
    tot.pers <- apply(MGP[[ps]][[i]][[2]] + MGP[[ps]][[i]][[3]],1,mean)
    if (ps == 1) {
      if (i == 3) {
        plot(ext.seq/c,tot.pers,type="l",lwd = 2,col=colors[ps],xlim=c(0,2),ylim=c(0,1))
      } else {
        plot(ext.seq/c,tot.pers,type="l",lwd = 2,col=colors[ps],xlim=c(0,2),ylim=c(0,1), xaxt='n', ann=FALSE)
      }
      if (i == 1){
        legend(1.7,1.1,legend = p.single.vec,col=colors,bty='n',pch=16)
      }
    } else {
      lines(ext.seq/c,tot.pers,type="l",lwd=2,col=colors[ps])
    }
  }
}
par(op)
dev.off()


#Plotting the Graphs
colors <- brewer.pal(3,"Set1")
# op <- par(mfrow = c(1,3),
#           oma = c(0,0,0,0) + 0.1,
#           mar = c(0,3,1,1) + 0.1)
layout(matrix(c(1,2,3), 1, 3, byrow = TRUE), 
   widths=c(1,1,1), heights=c(0.5,0.5,0.5))
par(oma = c(3, 4, 1, 1), mar = c(1, 1, 0, 1)) #,mai=c(0.6,0.6,0,0.1)
plot(graphs[[1]],vertex.size=degree(graphs[[1]])+2,vertex.label=NA,vertex.color=colors[3],edge.color="lightblue",edge.width=1,edge.arrow.size=0.6)
plot(graphs[[2]],vertex.size=degree(graphs[[2]]),vertex.label=NA,vertex.color=colors[1],edge.color="lightblue",edge.arrow.size=0.6)
plot(graphs[[3]],vertex.size=degree(graphs[[3]])+2,vertex.label=NA,vertex.color=colors[2],edge.color="lightblue",edge.width=1,edge.arrow.size=0.6)
#plot(dir.graphs[[1]],vertex.size=degree(graphs[[1]]),vertex.label=NA,vertex.color=colors[1],edge.color="lightblue",edge.arrow.size=0.6)
#plot(dir.graphs[[2]],vertex.size=degree(graphs[[2]])+2,vertex.label=NA,vertex.color=colors[2],edge.color="lightblue",edge.width=2,edge.arrow.size=0.6)
#plot(dir.graphs[[3]],vertex.size=degree(graphs[[3]])+2,vertex.label=NA,vertex.color=colors[3],edge.color="lightblue",edge.width=2,edge.arrow.size=0.6)
# par(op)

load("/Users/justinyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/Results/ESLdir_m0.1.RData")
pdf_filename <- paste("/Users/justinyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/Results/ESLdir_m",m,".pdf",sep="")






############################################################################################
#Determine the DIFFERENCE between p.single=1 and p.single=0 across ranges of e/c and rescue
############################################################################################




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
    
    #} #end r
    #} #end j
    
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

load(image_filename)

colors <- brewer.pal(11,"Spectral")


for (i in 1:3) {
  pdf_filename <- paste("/Users/justinyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/Results2/ESL_DirectedDifference_Graph",i,".pdf",sep="")
  pdf(file=pdf_filename,width=12,height=8)
  op <- par(mfrow = c(3,4),
            oma = c(5,4,0,0) + 0.1,
            mar = c(0,3,3,1) + 0.1,
            mgp = c(2, 1, 0))
  cum.diff.matrix <- matrix(0,length(ext.seq),length(m.seq))
  for (resc in 1:length(m.seq)) {
    cum.prop.undir <- apply(MRG[[i]][[resc]][[2]],1,mean) + apply(MRG[[i]][[resc]][[3]],1,mean)
    cum.prop.dir <- apply(MRG.dir[[i]][[resc]][[2]],1,mean) + apply(MRG.dir[[i]][[resc]][[3]],1,mean)
    cum.diff <- cum.prop.undir - cum.prop.dir
    #cum.diff.nonzero <- cum.diff[which(cum.diff != 0)]
    cum.diff.matrix[,resc] <- cum.diff
    plot(ext.seq/c,cum.prop.undir,type="l",ylim=c(0,1),col=colors[resc],lwd=2)
    lines(ext.seq/c,cum.prop.dir,lty=2,col=colors[resc],lwd=2)
  }
  par(op)
  dev.off()
}


#BoxPlot
colors <- brewer.pal(11,"RdYlBu")
pdf_filename <- paste("/Users/justinyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/Results/ESL_DirectedDifference_Boxplots.pdf",sep="")
pdf(file=pdf_filename,width=8,height=12)
op <- par(mfrow = c(3,1),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,3,1,1) + 0.1,
          mgp = c(2, 1, 0))
for (i in 1:3) {
  cum.diff.matrix <- matrix(0,length(ext.seq),length(m.seq))
  for (resc in 1:length(m.seq)) {
    cum.prop.undir <- apply(MRG[[i]][[resc]][[2]],1,mean) + apply(MRG[[i]][[resc]][[3]],1,mean)
    cum.prop.dir <- apply(MRG.dir[[i]][[resc]][[2]],1,mean) + apply(MRG.dir[[i]][[resc]][[3]],1,mean)
    cum.diff <- cum.prop.undir - cum.prop.dir
    #cum.diff.nonzero <- cum.diff[which(cum.diff != 0)]
    cum.diff.matrix[,resc] <- cum.diff
  }
  boxplot(cum.diff.matrix,boxwex=0.5,col=colors,names=m.seq,ylim=c(0,1),ylab=c("Cum. proportion differences (undirected-directed)"))
}
par(op)
dev.off()












#### R Version



pr.Grow <- matrix(0,n,t.term)
pr.NoGrowRescue <- matrix(0,n,t.term)
pr.NoGrowNoRescueStay <- matrix(0,n,t.term)
pr.NoGrowNoRescueExtinct <- matrix(0,n,t.term)
pr.Colonize <- matrix(0,n,t.term)

nn.state <- list()

for (t in 1:(t.term-1)) {
  
  XE <- which(X[,t]==0)
  num.E <- length(XE)
  
  XS <- which(X[,t]==1)
  num.S <- length(XS)
  
  XL <- which(X[,t]==2)
  num.L <- length(XL)
  
  #Determine state for each nearest neighbor of i
  for (i in 1:n) {
    nn.state[[i]] <- X[nn[[i]],t]
  }
  
  #Empty Dynamics
  for (i in XE) {
    #How many Large nodes are nearest neighbors to node i?
    num.L <- length(which(nn.state[[i]] == 2))
    
    pr.ES <- 1-(1-c)^num.L
    
    draw.ES <- runif(1)
    
    if (draw.ES < pr.ES) {
      X[i,t+1] <- 1
    } else {X[i,t+1] <- 0}
    
    pr.Colonize[i,t] <- pr.ES
  }
  
  #Small Dynamics
  for (i in XS) {
    #How many Large nodes are nearest neighbors to node i?
    num.L <- length(which(nn.state[[i]] == 2))
    
    pr.SL.grow <- r
    draw.SL.grow <- runif(1)
    if (draw.SL.grow < pr.SL.grow) {
      X[i,t+1] <- 2
    } else {
      pr.SL.rescue <- (1-(1-m)^num.L)
      draw.SL.rescue <- runif(1)
      if (draw.SL.rescue < pr.SL.rescue) {
        X[i,t+1] <- 2
      } else {
        pr.S.stay <- (1 - es)
        draw.S.stay <- runif(1)
        if (draw.S.stay < pr.S.stay) {
          X[i,t+1] <- 1
        } else {X[i,t+1] <- 0}
      }
    }
    #What can happen to you as an S?
    #The sum better Sam-Hell equal 1!
    pr.Grow[i,t] <- pr.SL.grow
    pr.NoGrowRescue[i,t] <- (1-pr.SL.grow)*pr.SL.rescue
    pr.NoGrowNoRescueStay[i,t] <- (1-pr.SL.grow)*(1-pr.SL.rescue)*(1-es)
    pr.NoGrowNoRescueExtinct[i,t] <- (1-pr.SL.grow)*(1-pr.SL.rescue)*es
  }
 
  
  
  
  #Large Dynamics
  for (i in XL) {
    pr.LS.shrink <- el
    draw.LS.shrink <- runif(1)
    if (draw.LS.shrink < pr.LS.shrink) {
      X[i,t+1] <- 1
    } else {
      X[i,t+1] <- 2
    }
  }

  XE.next <- which(X[,t+1]==0)
  num.E.next <- length(XE.next)
  
  XS.next <- which(X[,t+1]==1)
  num.S.next <- length(XS.next)
  
  XL.next <- which(X[,t+1]==2)
  num.L.next <- length(XL.next)
  
  Meta[,t+1] <- c(num.E.next,num.S.next,num.L.next)

}

#par(mfrow=c(2,1))
#E = Red
#S = Blue
#L = Green
colors <- brewer.pal(3,"Set1"); trans <- ""
plot(Meta[1,]/n,type="l",col=paste(colors[1],trans,sep=""),ylim=c(0,1),xlim=c(0,t.term),lwd=3)
lines(Meta[2,]/n,col=paste(colors[2],trans,sep=""),lwd=3)
lines(Meta[3,]/n,col=paste(colors[3],trans,sep=""),lwd=3)
#Rescue = green
#Stay = orange
#Extinct = blue
#Colonize = pink
colors <- brewer.pal(4,"Set2"); trans <- ""
plot(apply(pr.NoGrowRescue,2,mean),type=2,col=paste(colors[1],trans,sep=""),ylim=c(0,0.2),xlim=c(0,t.term),lwd=3)
lines(apply(pr.NoGrowNoRescueStay,2,mean),col=paste(colors[2],trans,sep=""),lwd=3)
lines(apply(pr.NoGrowNoRescueExtinct,2,mean),col=paste(colors[3],trans,sep=""),lwd=3)
lines(apply(pr.Colonize,2,mean),col=paste(colors[4],trans,sep=""),lwd=3)

colors <- brewer.pal(4,"Set2"); trans <- ""
plot(apply(pr.NoGrowRescue,2,sd),type=2,col=paste(colors[1],trans,sep=""),ylim=c(0,0.5),xlim=c(0,t.term),lwd=3)
lines(apply(pr.NoGrowNoRescueStay,2,sd),col=paste(colors[2],trans,sep=""),lwd=3)
lines(apply(pr.NoGrowNoRescueExtinct,2,sd),col=paste(colors[3],trans,sep=""),lwd=3)
lines(apply(pr.Colonize,2,sd),col=paste(colors[4],trans,sep=""),lwd=3)



par(mfrow=c(1,1))
colors <- brewer.pal(3,"Set1"); trans <- ""
plot((Meta[2,]+Meta[3,])/n,type=2,col=paste(colors[2],trans,sep=""),ylim=c(0,1),xlim=c(0,t.term),lwd=3)



