rm(list=c(ls()))
library(igraph)
source("RivOCN_func.R")
L <- 32
riv.net <- RivOCN_func(L=L,time=5)

plot(riv.net$late.net,layout=t(riv.net$coords),vertex.size=log(riv.net$late.area+1),
     vertex.color="lightblue",edge.color="darkgrey",edge.arrow.mode=1,edge.arrow.size=0.3,
     vertex.label=NA)


########################
# Random trees vs. OCN

net <- riv.net$early.net #River networks
area <- riv.net$early.area #River area

net <- riv.net$late.net #River networks
area <- riv.net$late.area #River area
s.area <- area/max(area) #Standardized river area

#Simple Levins metapopulation model
#Colonization and extinction probabilities
node.trajectories <- list()
mean.node.occupancy <- list()
ext.vector <- seq(0,0.6,0.02)
ext.vector <- 0.3
mean.occupancy <- numeric(length(ext.vector))
for (n in 1:length(ext.vector)) {
  pr.ext <- numeric(L*L) + ext.vector[n]
  pr.col <- 0.6
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
occupancy.river <- mean.occupancy
occupancy.random <- mean.occupancy

###########################
# Treat probability of extinction as a function of river area

net <- riv.net$early.net #River networks
area <- riv.net$early.area #River area
s.area <- area/max(area) #Standardized river area

net <- riv.net$late.net #River networks
area <- riv.net$late.area #River area
s.area <- area/max(area) #Standardized river area

#Define curvature of extinction risk/Area relationship
b <- 1

#Simple Levins metapopulation model
#Colonization and extinction probabilities
node.trajectories <- list()
mean.node.occupancy <- list()
max.extinction <- seq(0,0.6,0.05)
#max.extinction <- 0.3
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

plot(ratio.vector.river,occupancy.river,type="o",col="dodgerblue3",
     ylab="Relative abundance",xlab="e/c",xlim=c(0,1))
points(ratio.vector.rand,occupancy.rand,type="o",col="firebrick3")

#Degree vs. node occupancy
plot(degree(net),mean.node.occupancy[[1]],pch=16)

#Distance from confluence vs. node occupancy
plot(unlist(lapply(seq(1,L*L),function(x){shortest.paths(net,x,1,mode="all")})),
     mean.node.occupancy[[1]],pch=16)
