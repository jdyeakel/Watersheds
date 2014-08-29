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
X[,1] <- c(rep("E",initial[1]),rep("S",initial[2]),rep("L",initial[3]))
Meta[,1] <- initial

#Full Network :: every node is connected to every other node
adj.m <- matrix(1,n,n); diag(adj.m) <- 0

#Lattice Network
lattice.net <- graph.lattice(dimvector=c(sqrt(n),sqrt(n)))
adj.m <- get.adjacency(lattice.net)
isSymmetric(adj.m)

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





#River Network

colors <- brewer.pal(3,"Set1"); trans <- ""
plot(graph.adjacency(adj.m,mode="undirected"),vertex.size=2,vertex.label=NA,vertex.color=colors[2])
###########################
###########################
###########################
###########################
###########################


#Probabilities
c <- 0.2
r <- 0.2
m <- 0.1
el <- 0.2
es <- 2*el

pr.Grow <- matrix(0,n,t.term)
pr.NoGrowRescue <- matrix(0,n,t.term)
pr.NoGrowNoRescueStay <- matrix(0,n,t.term)
pr.NoGrowNoRescueExtinct <- matrix(0,n,t.term)
pr.Colonize <- matrix(0,n,t.term)

#Make list of nearest neighbors for each i
nn <- list()
for (i in 1:n) {
  arow <- adj.m[i,]
  nn[[i]] <- which(arow==1)
}
nn.state <- list()

for (t in 1:(t.term-1)) {
  
  XE <- which(X[,t]=="E")
  num.E <- length(XE)
  
  XS <- which(X[,t]=="S")
  num.S <- length(XS)
  
  XL <- which(X[,t]=="L")
  num.L <- length(XL)
  
  #Determine state for each nearest neighbor of i
  for (i in 1:n) {
    nn.state[[i]] <- X[nn[[i]],t]
  }
  
  #Empty Dynamics
  for (i in XE) {
    #How many Large nodes are nearest neighbors to node i?
    num.L <- length(which(nn.state[[i]] == "L"))
    
    pr.ES <- 1-(1-c)^num.L
    
    draw.ES <- runif(1)
    
    if (draw.ES < pr.ES) {
      X[i,t+1] <- "S"
    } else {X[i,t+1] <- "E"}
    
    pr.Colonize[i,t] <- pr.ES
  }
  
  #Small Dynamics
  for (i in XS) {
    #How many Large nodes are nearest neighbors to node i?
    num.L <- length(which(nn.state[[i]] == "L"))
    
    pr.SL.grow <- r
    draw.SL.grow <- runif(1)
    if (draw.SL.grow < pr.SL.grow) {
      X[i,t+1] <- "L"
    } else {
      pr.SL.rescue <- (1-(1-m)^num.L)
      draw.SL.rescue <- runif(1)
      if (draw.SL.rescue < pr.SL.rescue) {
        X[i,t+1] <- "L"
      } else {
        pr.S.stay <- (1 - es)
        draw.S.stay <- runif(1)
        if (draw.S.stay < pr.S.stay) {
          X[i,t+1] <- "S"
        } else {X[i,t+1] <- "E"}
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
      X[i,t+1] <- "S"
    } else {
      X[i,t+1] <- "L"
    }
  }

  XE.next <- which(X[,t+1]=="E")
  num.E.next <- length(XE.next)
  
  XS.next <- which(X[,t+1]=="S")
  num.S.next <- length(XS.next)
  
  XL.next <- which(X[,t+1]=="L")
  num.L.next <- length(XL.next)
  
  Meta[,t+1] <- c(num.E.next,num.S.next,num.L.next)

}

par(mfrow=c(2,1))
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
plot(apply(pr.NoGrowRescue,2,mean),type="l",col=paste(colors[1],trans,sep=""),ylim=c(0,0.2),xlim=c(0,t.term),lwd=3)
lines(apply(pr.NoGrowNoRescueStay,2,mean),col=paste(colors[2],trans,sep=""),lwd=3)
lines(apply(pr.NoGrowNoRescueExtinct,2,mean),col=paste(colors[3],trans,sep=""),lwd=3)
lines(apply(pr.Colonize,2,mean),col=paste(colors[4],trans,sep=""),lwd=3)

colors <- brewer.pal(4,"Set2"); trans <- ""
plot(apply(pr.NoGrowRescue,2,sd),type="l",col=paste(colors[1],trans,sep=""),ylim=c(0,0.5),xlim=c(0,t.term),lwd=3)
lines(apply(pr.NoGrowNoRescueStay,2,sd),col=paste(colors[2],trans,sep=""),lwd=3)
lines(apply(pr.NoGrowNoRescueExtinct,2,sd),col=paste(colors[3],trans,sep=""),lwd=3)
lines(apply(pr.Colonize,2,sd),col=paste(colors[4],trans,sep=""),lwd=3)



par(mfrow=c(1,1))
colors <- brewer.pal(3,"Set1"); trans <- ""
plot((Meta[2,]+Meta[3,])/n,type="l",col=paste(colors[2],trans,sep=""),ylim=c(0,1),xlim=c(0,t.term),lwd=3)



