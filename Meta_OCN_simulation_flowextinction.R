rm(list=c(ls()))
library(igraph)
source("RivOCN_func.R")
#source("EdenGrowth_func.R")
source("EdenGrowth_func2.R")
L <- 32
riv.net <- RivOCN_func(L=L,tic_max=10000)
save(riv.net,file="riv_net400.RData")
load("riv_net400.RData")
plot(riv.net$late.net,layout=t(riv.net$coords),vertex.size=log(riv.net$late.area+1),
     vertex.color="lightblue",edge.color="darkgrey",edge.arrow.mode=1,edge.arrow.size=0.3,
     vertex.label=NA)

#Build network using Eden growth model
#eden <- EdenGrowth_func(size=(300))
eden <- EdenGrowth_func3(size=(L*L),boundary=20)
save(eden,file="DLA_net400_B30.RData")
load("DLA_net400_B30.RData")
#write.table(eden$edge.list,"DLA_net400_B30.txt",sep="\t",col.names=F,row.names=F)


g <- graph.edgelist(eden$edge.list)
edge.weight <- numeric(length(eden$edge.list[,1]))
for (i in 1:length(eden$edge.list[,1])) {
  edge1 <- eden$edge.list[i,1]
  edge2 <- eden$edge.list[i,2]
  edge.weight[i] <- min(c(eden$area[edge1],eden$area[edge2]))
}
plot(g,layout=eden$eden.tree,vertex.size=0.01,
     vertex.color="lightblue",edge.color="dodgerblue3",edge.arrow.mode=0,
     edge.width=log(edge.weight+1)+2,vertex.label=NA)

plot(g, layout=eden$eden.tree,
     vertex.color="lightblue",edge.color="darkgrey",edge.arrow.mode=0,edge.arrow.size=0.3,
     vertex.label=NA,vertex.size=log(eden$area+1))
plot(eden$eden.tree,pch=15,cex=1,xlim=c(0,60),ylim=c(0,60))



# River network with edge weights equal to river area
late.edges <- get.edgelist(riv.net$late.net)
late.edge.weight <- numeric(length(late.edges[,1]))
early.edges <- get.edgelist(riv.net$early.net)
early.edge.weight <- numeric(length(early.edges[,1]))
for (i in 1:length(late.edges[,1])) {
  edge1 <- late.edges[i,1]
  edge2 <- late.edges[i,2]
  late.edge.weight[i] <- min(c(riv.net$late.area[edge1],riv.net$late.area[edge2]))
  edge1 <- early.edges[i,1]
  edge2 <- early.edges[i,2]
  early.edge.weight[i] <- min(c(riv.net$early.area[edge1],riv.net$early.area[edge2]))
}
plot(riv.net$late.net,layout=t(riv.net$coords),vertex.size=0.01,
     vertex.color="lightblue",edge.color="dodgerblue3",edge.arrow.mode=0,
     edge.width=log(late.edge.weight+1)+2,vertex.label=NA)

plot(riv.net$early.net,layout=t(riv.net$coords),vertex.size=0.01,
     vertex.color="lightblue",edge.color="dodgerblue3",edge.arrow.mode=0,
     edge.width=log(early.edge.weight+1)+2,vertex.label=NA)


#Other Interesting layouts
#layout.kamada.kawai
plot(riv.net$late.net,layout=layout.kamada.kawai(riv.net$late.net),vertex.size=0.01,
     vertex.color="lightblue",edge.color="dodgerblue3",edge.arrow.mode=0,
     edge.width=log(late.edge.weight+1)+2,vertex.label=NA)

########################
# Random trees vs. OCN

#Todo: Feb 28 2014
#1 Build in repetitions
for (Z in 1:3) {
  print(paste("Z=",Z))
  if (Z == 1) {
    net <- riv.net$late.net #River networks
    area <- riv.net$late.area #River area
    s.area <- area/max(area) #Standardized river area
  }
  if (Z == 2) {
    net <- riv.net$early.net #River networks
    area <- riv.net$early.area #River area
    s.area <- area/max(area) #Standardized river area
  }
  if (Z == 3) {
    net <- graph.edgelist(eden$edge.list)
    area <- eden$area
    s.area <- area/max(area) #Standardized river area
  }
  
  #Simple Levins metapopulation model
  #Colonization and extinction probabilities
  node.trajectories <- list()
  mean.node.occupancy <- list()
  reps <- 20
  ext.vector <- seq(0,0.5,0.01)
  #ext.vector <- 0.3
  mean.occupancy <- numeric(length(ext.vector))
  sd.occupancy <- numeric(length(ext.vector))
  for (n in 1:length(ext.vector)) {
    mean.occupancy.rep <- numeric(reps)
    mean.node.occupancy.rep <- matrix(0,L*L,reps)
    for (R in 1:reps) {
      pr.ext <- numeric(L*L) + ext.vector[n]
      pr.col <- 0.6
      pr.spill <- 0
      start.cond <- sample(c(0,1),L*L,replace=TRUE)
      start.cond <- numeric(L*L)+1
      t.term <- 1000
      t.mid <- round(t.term/1.5,0)
      node.cond <- matrix(0,L*L,t.term)
      node.cond[,1] <- start.cond
      for (t in 1:(t.term-1)) {
        #Which nodes are occupied
        node.occupied <- which(node.cond[,t] == 1)
        #If there are occupied nodes in the system...
        if (length(node.occupied) > 0) {
          #EXTINCTIONS
          extinct <- 1-as.numeric(runif(length(node.occupied)) < pr.ext)
          #COLONIZATIONS
          #Find neighboring nodes for each occupied node
          node.nn <- lapply(node.occupied,function(x){neighbors(net,x,mode="total")})
          #Draw colonization probabilities for each neighbor
          colonize <- lapply(node.nn,function(x){as.numeric(runif(length(x)) < pr.col)})
          nodes.colonize <- unlist(node.nn)*unlist(colonize)
          nodes.colonize <- nodes.colonize[nodes.colonize != 0]
          #Take old nodes and update
          #Step 1: apply colonizations
          new.cond <- node.cond[,t]
          new.cond[nodes.colonize] <- 1
          #Step 2: apply extinctions
          new.cond <- new.cond * extinct
          #Spill Extinctions
          spill <- runif(1) < pr.spill
          if (spill) {
            #Choose random starting node
            node.rand <- sample(seq(1,L*L),1)
            #Find path from random node to confluence
            spill.nodes <- unlist(get.shortest.paths(net,node.rand,1,mode="all"))
            #eliminate those populations
            new.cond[spill.nodes] <- 0
          }
          #Update
          node.cond[,t+1] <- new.cond
        } else {node.cond[,t+1] <- node.cond[,t]}
      } # End time loop
      occupancy <- apply(node.cond,2,sum)/(L*L)
      mean.occupancy.rep[R] <- mean(occupancy[t.mid:t.term])
      mean.node.occupancy.rep[,R] <- apply(node.cond,1,function(x){sum(x)/t.term})
    } # End repetitions loop
    mean.occupancy[n] <- mean(mean.occupancy.rep)
    sd.occupancy[n] <- sd(mean.occupancy.rep)
    mean.node.occupancy[[n]] <- apply(mean.node.occupancy.rep,1,mean)
    print(n)
  }
  if (Z==1) {
    occupancy.river <- mean.occupancy
    occupancy.river.sd <- sd.occupancy
  }
  if (Z==2) {
    occupancy.rand <- mean.occupancy
    occupancy.rand.sd <- sd.occupancy
  }
  if (Z==3) {
    occupancy.eden <- mean.occupancy
    occupancy.eden.sd <- sd.occupancy
  }
}
save.image("ec_sim.RData")
plot(ext.vector/pr.col,occupancy.river,type="l",lwd=3,col="dodgerblue3",
     ylab="Relative abundance",xlab="e/c",xlim=c(0,1))
polygon(x=c(ext.vector/pr.col,rev(ext.vector/pr.col)),
        y=c(occupancy.river+occupancy.river.sd,rev(occupancy.river-occupancy.river.sd)),
        col="#1874CD50",border=NA)
polygon(x=c(ext.vector/pr.col,rev(ext.vector/pr.col)),
        y=c(occupancy.rand+occupancy.rand.sd,rev(occupancy.rand-occupancy.rand.sd)),
        col="#CD262650",border=NA)
points(ext.vector/pr.col,occupancy.rand,type="l",lwd=3,col="firebrick3")
polygon(x=c(ext.vector/pr.col,rev(ext.vector/pr.col)),
        y=c(occupancy.eden+occupancy.eden.sd,rev(occupancy.eden-occupancy.eden.sd)),
        col="#228B2250",border=NA)
points(ext.vector/pr.col,occupancy.eden,type="l",lwd=3,col="forestgreen")


###########################
# Treat probability of extinction as a function of river area

net <- riv.net$late.net #River networks
area <- riv.net$late.area #River area
s.area <- area/max(area) #Standardized river area

net <- riv.net$early.net #River networks
area <- riv.net$early.area #River area
s.area <- area/max(area) #Standardized river area

net <- graph.edgelist(eden$edge.list)
area <- eden$area
s.area <- area/max(area) #Standardized river area

#Define curvature of extinction risk/Area relationship
b <- 1

#Simple Levins metapopulation model
#Colonization and extinction probabilities
node.trajectories <- list()
mean.node.occupancy <- list()
max.extinction <- seq(0,0.6,0.05)
max.extinction <- 0.3
ratio.vector <- numeric(length(max.extinction))
mean.pr.ext <- numeric(length(ratio.vector))
mean.occupancy <- numeric(length(ratio.vector))
for (n in 1:length(ratio.vector)) {
  #Define maximum extinction risk
  a <- max.extinction[n]
  pr.ext <- a*exp(-b*s.area)
  mean.pr.ext[n] <- median(pr.ext)
  pr.col <- 0.6
  #Calculate e/c ratio
  ratio.vector[n] <- mean.pr.ext[n]/pr.col
  start.cond <- sample(c(0,1),L*L,replace=TRUE)
  start.cond <- numeric(L*L)+1
  t.term <- 1000
  t.mid <- round(t.term/2,0)
  node.cond <- matrix(0,L*L,t.term)
  node.cond[,1] <- start.cond
  for (t in 1:(t.term-1)) {
    #Which nodes are occupied
    node.occupied <- which(node.cond[,t] == 1)
    #EXTINCTIONS
    extinct <- 1-as.numeric(runif(length(node.occupied)) < pr.ext)
    #COLONIZATIONS
    #Find neighboring nodes for each occupied node
    node.nn <- lapply(node.occupied,function(x){neighbors(net,x,mode="total")})
    #Draw colonization probabilities for each neighbor
    colonize <- lapply(node.nn,function(x){as.numeric(runif(length(x)) < pr.col)})
    nodes.colonize <- unlist(node.nn)*unlist(colonize)
    nodes.colonize <- nodes.colonize[nodes.colonize != 0]
    #Take old nodes and update
    #Step 1: apply colonizations
    new.cond <- node.cond[,t]
    new.cond[nodes.colonize] <- 1
    #Step 2: apply extinctions
    new.cond <- new.cond * extinct
    #Update
    node.cond[,t+1] <- new.cond
  }
  node.trajectories[[n]] <- node.cond
  occupancy <- apply(node.cond,2,sum)/(L*L)
  mean.occupancy[n] <- mean(occupancy[500:1000])
  mean.node.occupancy[[n]] <- apply(node.cond,1,function(x){sum(x)/t.term})
  print(n)
}
occupancy.river <- mean.occupancy; ratio.vector.river <- ratio.vector
occupancy.rand <- mean.occupancy; ratio.vector.rand <- ratio.vector


###################### 
# PLOTS
# Metapopulation trajectory plot
plot(apply(node.cond,2,sum)/(L*L),type="l",ylim=c(0,1))
palette(gray(seq(0,1,0.01)))
plot(net,layout=t(riv.net$coords),vertex.size=log(area+1),
     vertex.color=node.cond[,t.term],edge.color="darkgrey",edge.arrow.mode=1,
     edge.arrow.size=0.3,vertex.label=NA)
palette(rainbow(100))
palette(gray(seq(0,1,0.01)))
ec <- 1
plot(net,layout=t(riv.net$coords),vertex.size=log(area+1),
     vertex.color=(100-mean.node.occupancy[[ec]]*100),edge.color="darkgrey",edge.arrow.mode=0,
     edge.arrow.size=0.3,vertex.label=NA)
hist(area*node.cond[,t.term])
hist(area)

ec <- 1
plot(net,layout=t(riv.net$coords),vertex.size=log(area+1),
     vertex.color=(100-mean.node.occupancy[[ec]]*100),edge.color="dodgerblue3",edge.arrow.mode=0,
     edge.width=log(late.edge.weight+1)+2,vertex.label=NA)

hist(degree(riv.net$early.net))

plot(ratio.vector.river,occupancy.river,type="o",col="dodgerblue3",
     ylab="Relative abundance",xlab="e/c",xlim=c(0,1))
points(ratio.vector.rand,occupancy.rand,type="o",col="firebrick3")

#Degree vs. node occupancy
plot(degree(net),mean.node.occupancy[[1]],pch=16)

#Distance from confluence vs. node occupancy
plot(unlist(lapply(seq(1,L*L),function(x){shortest.paths(net,x,1,mode="all")})),
     mean.node.occupancy[[1]],pch=16)
