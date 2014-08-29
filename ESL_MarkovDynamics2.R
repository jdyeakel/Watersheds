rm(list=c(ls()))

library(RColorBrewer)
library(igraph)
library(intergraph)
library(Rcpp)

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
isSymmetric(adj.m)

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
sourceCpp('src/ESL_meta.cpp')

t.term <- 5000
n <- 100
l <- 200
X <- matrix(0,n,t.term)
Meta <- matrix(0,3,t.term)
#Initiate starting vector
initial <- c(n - 2*round(n/3,0),round(n/3,0),round(n/3,0))
X[,1] <- c(rep(0,initial[1]),rep(1,initial[2]),rep(2,initial[3]))
Meta[,1] <- initial

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
el <- 0.2
es <- 2*el

#### C++ Version
nn_cpp = lapply(nn,function(x){x-1})
ESL.out <- ESL_meta(n, t.term, X, Meta, nn_cpp, c, r, m, es, el)

colors <- brewer.pal(3,"Set1"); trans <- ""
plot(ESL.out[[2]][1,]/n,type="l",col=paste(colors[1],trans,sep=""),ylim=c(0,1),xlim=c(0,t.term),lwd=3)
lines(ESL.out[[2]][2,]/n,col=paste(colors[2],trans,sep=""),lwd=3)
lines(ESL.out[[2]][3,]/n,col=paste(colors[3],trans,sep=""),lwd=3)


# Analysis of C++ simulation
# Across different types of networks

#Length of simulation
t.term <- 5000
#Number of nodes
n <- 100
#Number of links
l <- 200

MG <- list()
graphs <- list()
for (i in 1:3) {
  #Select Network
  if (i == 1) {
    #Full Network :: every node is connected to every other node
    adj.m <- matrix(1,n,n); diag(adj.m) <- 0
    graphs[[i]] <- graph.adjacency(adj.m)
  }
  if (i == 2) {
    lattice.net <- graph.lattice(dimvector=c(sqrt(n),sqrt(n)))
    adj.m <- get.adjacency(lattice.net)
    #isSymmetric(as.matrix(adj.m))
    graphs[[i]] <- lattice.net
  }
  if (i == 3) {
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
    adj.m <- get.adjacency(sf.net)
    #isSymmetric(adj.m)
    graphs[[i]] <- sf.net
  }
  
  #Make list of nearest neighbors for each node
  nn <- list()
  for (j in 1:n) {
    arow <- adj.m[j,]
    nn[[j]] <- which(arow==1)
  }
  nn_cpp <- lapply(nn,function(x){x-1})
  
  #Extinction sequence (prob(ext) for Large patch)
  ext.seq <- seq(0,1,0.02)
  
  replicates <- 10
  m.e <- matrix(0,length(ext.seq),replicates)
  m.s <- matrix(0,length(ext.seq),replicates)
  m.l <- matrix(0,length(ext.seq),replicates)
  
  #Iterate over different values of large-patch extinction probabilities
  for (j in 1:length(ext.seq)) {
    for (rep in 1:replicates) {
      
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
      el <- ext.seq[j]
      es <- 2*el
      
      #Perform Simulation
      ESL.out <- ESL_meta(n, t.term, X, Meta, nn_cpp, c, r, m, es, el)
      
      #Save summary statistics
      m.e[j,rep] <- mean(ESL.out[[2]][1,round(t.term/2,0):t.term]/n)
      m.s[j,rep] <- mean(ESL.out[[2]][2,round(t.term/2,0):t.term]/n)
      m.l[j,rep] <- mean(ESL.out[[2]][3,round(t.term/2,0):t.term]/n)
      
    } #end r
  } #end j
  
  M <- list()
  
  M[[1]] <- m.e
  M[[2]] <- m.s
  M[[3]] <- m.l
  
  MG[[i]] <- M
  
}# end i
image_filename <- paste("/Users/justinyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/Results/ESL_m",m,".RData",sep="")
save.image(image_filename)

load("/Users/justinyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/Results/ESL_m0.RData")
pdf_filename <- paste("/Users/justinyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/Results/ESL_m",m,".pdf",sep="")

dev.off()
pdf(file=pdf_filename
    ,width=3.5,height=6)
op <- par(mfrow = c(3,1),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,3,1,1) + 0.1,
          mgp = c(2, 1, 0))
for (i in 1:3) {
  
  m.e <- MG[[i]][[1]]
  m.s <- MG[[i]][[2]]
  m.l <- MG[[i]][[3]]
  if (i == 3) {
    plot(ext.seq/c,apply(m.e,1,mean),type="l",col=colors[1],xlim=c(0,2),ylim=c(0,1),xlab="e/c",ylab="Proportion")
  } else {plot(ext.seq/c,apply(m.e,1,mean),type="l",col=colors[1],xlim=c(0,1.5),ylim=c(0,1), xaxt='n', ann=FALSE)}
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

#Plotting the Graphs
plot(lattice.net,vertex.size=degree(lattice.net),vertex.label=NA,vertex.color=colors[2],edge.color="lightblue")
plot(sf.net,vertex.size=degree(sf.net)+2,vertex.label=NA,vertex.color=colors[3],edge.color="lightblue",edge.width=2)



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



