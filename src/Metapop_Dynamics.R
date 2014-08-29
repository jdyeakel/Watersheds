

###########################
# METAPOPULATION DYNAMICS #
###########################

rm(list=c(ls()))
library(igraph)
library(Rcpp)


## Load Source Files
##
sourceCpp('src/dyn_func.cpp')
sourceCpp('src/DLA_func.cpp')
source("RivOCN_func.R")

# Set network size ~ size = L^2
L <- 32
riv.net <- RivOCN_func(L=L,tic_max=10000)
save(riv.net,file="riv_net1024_2.RData")
load("riv_net1024.RData")
#load("riv_net400.RData")

dla <- DLA_func(size=(L*L),boundary=20)
save(dla,file="DLA_net1024_B20_2.RData")
load("DLA_net1024_B30.RData")
#load("DLA_net400_B30.RData")

t.term <- 5000 #10000
reps <- 10 #20
#ext.vector <- 0
pr.col <- 0.5
ext.vector <- seq(0,0.6,0.005) #0.005

fme <- 0.2
pr.fme.vector <- c(0,fme,0,fme)
#Probability of transecting colonizations
tc <- 0.2
pr.tc.vector <- c(0,0,tc,tc)

occupancy.river <- list()
occupancy.river.sd <- list()
per.node.river <- list()
per.node.river.sd <- list()
occupancy.random <- list()
occupancy.random.sd <- list()
occupancy.dla <- list()
occupancy.dla.sd <- list()
#Factorial Design
#1) Dynamic sim ~ normal
#2) Dynamic sim ~ with flow mediated extinctions (pr.fmc)
#3) Dynamic sim ~ with probabilistic transecting colonizations (pr.tc)
#4) Dynamic sim ~ with flow mediated extinctions + probabilistic transecting colonizations
F <- 3
for (f in 1:F) {
  print(paste("Iteration:",f))
  #Set attributes
  pr.fme <- pr.fme.vector[f]
  pr.tc <- pr.tc.vector[f]
  
  #Late Net
  net <- riv.net$late.net
  #Build list of nearest neighbors in the network
  nn.list <- lapply(seq(1,vcount(net),1),function(x){neighbors(net,x,mode="all")-1}) #-1 to convert to C++ inidces
  downstream.list <- lapply(get.shortest.paths(net,1,mode="out"),function(x){x-1}) #-1 to convert to C++ inidces
  output <- dyn_func(t.term, reps, nn.list, downstream.list, ext.vector, pr.col, pr.fme, pr.tc)
  occupancy.river[[f]] <- output[[1]]
  occupancy.river.sd[[f]] <- output[[2]]
  per.node.river[[f]] <- output[[3]]
  per.node.river.sd[[f]] <- output[[4]]
  
  #Early Net
  net <- riv.net$early.net
  #Build list of nearest neighbors in the network
  nn.list <- lapply(seq(1,vcount(net),1),function(x){neighbors(net,x,mode="all")-1}) #-1 to convert to C++ inidces
  downstream.list <- lapply(get.shortest.paths(net,1,mode="out"),function(x){x-1}) #-1 to convert to C++ inidces
  output <- dyn_func(t.term, reps, nn.list, downstream.list, ext.vector, pr.col, pr.fme, pr.tc)
  occupancy.random[[f]] <- output[[1]]
  occupancy.random.sd[[f]] <- output[[2]]
  
  #DLA Net
  net <- graph.edgelist(dla[[1]])
  #Build list of nearest neighbors in the network
  nn.list <- lapply(seq(1,vcount(net),1),function(x){neighbors(net,x,mode="all")-1}) #-1 to convert to C++ inidces
  downstream.list <- lapply(get.shortest.paths(net,1,mode="out"),function(x){x-1}) #-1 to convert to C++ inidces
  output <- dyn_func(t.term, reps, nn.list, downstream.list, ext.vector, pr.col, pr.fme, pr.tc)
  occupancy.dla[[f]] <- output[[1]]
  occupancy.dla.sd[[f]] <- output[[2]]
  
}

# add inset



save.image("ec_sim_Net1024_T5K_R10_E005_F3.RData")
load()

dev.off()
pdf(file="/Users/jyeakel/Dropbox/PostDoc/2014_Empirical_Watersheds/Manuscript/Figure_ec_fme_tc.pdf"
    ,width=6,height=5)
plot(ext.vector/pr.col,occupancy.river[[1]],type="l",lwd=2,col="dodgerblue3",
     ylab="Relative abundance",xlab="e/c",xlim=c(0,1))
polygon(x=c(ext.vector/pr.col,rev(ext.vector/pr.col)),
        y=c(occupancy.river[[1]]+occupancy.river.sd[[1]],rev(
          unlist(lapply(occupancy.river[[1]]-occupancy.river.sd[[1]],function(x){max(x,0)})))),
        col="#1874CD50",border=NA)
for (i in 2:F) {
  polygon(x=c(ext.vector/pr.col,rev(ext.vector/pr.col)),
          y=c(occupancy.river[[i]]+occupancy.river.sd[[i]],rev(
            unlist(lapply(occupancy.river[[i]]-occupancy.river.sd[[i]],function(x){max(x,0)})))),
          col="#1874CD50",border=NA)
  lines(ext.vector/pr.col,occupancy.river[[i]],type="l",lwd=2,lty=i,col="dodgerblue3")
}

for (i in 1:F) {
polygon(x=c(ext.vector/pr.col,rev(ext.vector/pr.col)),
        y=c(occupancy.random[[i]]+occupancy.random.sd[[i]],rev(
          unlist(lapply(occupancy.random[[i]]-occupancy.random.sd[[i]],function(x){max(x,0)})))),
        col="#CD262650",border=NA)
lines(ext.vector/pr.col,occupancy.random[[i]],type="l",lwd=2,lty=i,col="firebrick3")
}

for (i in 1:F) {
polygon(x=c(ext.vector/pr.col,rev(ext.vector/pr.col)),
        y=c(occupancy.dla[[i]]+occupancy.dla.sd[[i]],rev(
          unlist(lapply(occupancy.dla[[i]]-occupancy.dla.sd[[i]],function(x){max(x,0)})))),
        col="#228B2250",border=NA)
lines(ext.vector/pr.col,occupancy.dla[[i]],type="l",lwd=2,lty=i,col="forestgreen")
}

par(fig = c(.65, 0.9, .6, 0.8), mar=c(0,0,0,0), new=TRUE)
#plot(rnorm(5,1), col=2) # inset bottomright
plot(ext.vector/pr.col,occupancy.river[[1]],type="l",lwd=2,col="dodgerblue3",
     ylab="Relative abundance",xlab="e/c",xlim=c(0.25,0.45),ylim=c(0.6,0.8),cex.axis=0.8)
polygon(x=c(ext.vector/pr.col,rev(ext.vector/pr.col)),
        y=c(occupancy.river[[1]]+occupancy.river.sd[[1]],rev(
          unlist(lapply(occupancy.river[[1]]-occupancy.river.sd[[1]],function(x){max(x,0)})))),
        col="#1874CD50",border=NA)
for (i in 2:F) {
  polygon(x=c(ext.vector/pr.col,rev(ext.vector/pr.col)),
          y=c(occupancy.river[[i]]+occupancy.river.sd[[i]],rev(
            unlist(lapply(occupancy.river[[i]]-occupancy.river.sd[[i]],function(x){max(x,0)})))),
          col="#1874CD50",border=NA)
  lines(ext.vector/pr.col,occupancy.river[[i]],type="l",lwd=2,lty=i,col="dodgerblue3")
}

for (i in 1:F) {
  polygon(x=c(ext.vector/pr.col,rev(ext.vector/pr.col)),
          y=c(occupancy.random[[i]]+occupancy.random.sd[[i]],rev(
            unlist(lapply(occupancy.random[[i]]-occupancy.random.sd[[i]],function(x){max(x,0)})))),
          col="#CD262650",border=NA)
  lines(ext.vector/pr.col,occupancy.random[[i]],type="l",lwd=2,lty=i,col="firebrick3")
}

for (i in 1:F) {
  polygon(x=c(ext.vector/pr.col,rev(ext.vector/pr.col)),
          y=c(occupancy.dla[[i]]+occupancy.dla.sd[[i]],rev(
            unlist(lapply(occupancy.dla[[i]]-occupancy.dla.sd[[i]],function(x){max(x,0)})))),
          col="#228B2250",border=NA)
  lines(ext.vector/pr.col,occupancy.dla[[i]],type="l",lwd=2,lty=i,col="forestgreen")
}
dev.off()

i = 6
df <- cbind(degree(net),per.node.river[[1]][,i],per.node.river[[2]][,i])
df <- as.data.frame(df)
colnames(df) <- c("degree","pn1","pn2")
plot(unlist(lapply(split(df$pn1,df$degree),mean)),type="b",ylim=c(0.2,0.8))
points(unlist(lapply(split(df$pn2,df$degree),mean)),type="b")


#Per-node occupancy vs. degree
par(mfrow=c(1,1))
plot(degree(net),per.node.river[[1]][,6],xlim=c(0,5),ylim=c(0,0.7))
#Per-node occupancy vs. degree with fme
plot(degree(net),per.node.river[[2]][,6],xlim=c(0,5),ylim=c(0,0.7))

plot(per.node.river[[1]],per.node.river[[2]],pch=".",xlab="without fme",ylab="with fme")

#Plotting mean per-node occupancy on network
col <- RColorBrewer::brewer.pal(9, "YlOrRd")
pal <- colorRampPalette(col)
num.cols <- 100
gradient <- pal(num.cols)
values <- per.node.river[[1]][,3]
plot(net,layout=t(riv.net$coords),vertex.size=log(riv.net$late.area+1),
     vertex.color=col[values],edge.color="darkgrey",edge.arrow.mode=1,
     edge.arrow.size=0,vertex.label=NA)



#################################################
# METAPOPULATION DYNAMICS AS A FUNCTION OF AREA #
#################################################

#Includes area-dependent extinctions
#Flow-mediated colonization

rm(list=c(ls()))
library(igraph)
library(Rcpp)


## Load Source Files
##
sourceCpp('src/dyn_area_func.cpp')
sourceCpp('src/DLA_func.cpp')
source("RivOCN_func.R")

# Set network size ~ size = L^2
L <- 32
riv.net <- RivOCN_func(L=L,tic_max=10000)
save(riv.net,file="riv_net1024_2.RData")
load("riv_net1024.RData")
#load("riv_net400.RData")

dla <- DLA_func(size=(L*L),boundary=200)
save(dla,file="DLA_net1024_B30.RData")
load("DLA_net1024_B30.RData")
# Calculate area
# First we need to orient all directed edges to flow towards the confluence
# Identify Confluence and tributaries
g <- graph.edgelist(dla[[1]])
size <- vcount(g)
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
dla.area <- r.area
edge.weight <- numeric(length(dla[[1]][,1]))
for (i in 1:length(dla[[1]][,1])) {
  edge1 <- dla[[1]][i,1]
  edge2 <- dla[[1]][i,2]
  edge.weight[i] <- min(c(r.area[edge1],r.area[edge2]))
}
plot(g,layout=dla[[2]],vertex.size=0.01,
     vertex.color="transparent",vertex.frame.color=NA,edge.color="dodgerblue3",edge.arrow.mode=0,
     edge.width=log(edge.weight+1)/2,vertex.label=NA)

#load("DLA_net400_B30.RData")

t.term <- 5000 #10000
reps <- 20 #20
#ext.vector <- 0
pr.col <- 0.5
max.ext.vector <- seq(0,0.8,0.005) #0.005


fme <- 0.2
pr.fme.vector <- c(0,fme,0,fme)
#Probability of transecting colonizations
tc <- 0.2
pr.tc.vector <- c(0,0,tc,tc)

occupancy.river <- list()
occupancy.river.sd <- list()
per.node.river <- list()
per.node.river.sd <- list()
occupancy.random <- list()
occupancy.random.sd <- list()
per.node.random <- list()
per.node.random.sd <- list()
occupancy.dla <- list()
occupancy.dla.sd <- list()
per.node.dla <- list()
per.node.dla.sd <- list()

F <- 2
for (f in 3:3) {
  print(paste("Iteration:",f))
  #Set attributes
  pr.fme <- pr.fme.vector[f]
  pr.tc <- pr.tc.vector[f]
  
  #Late Net
  net <- riv.net$late.net
  area <- riv.net$late.area
  #Build list of nearest neighbors in the network
  nn.list <- lapply(seq(1,vcount(net),1),function(x){neighbors(net,x,mode="all")-1}) #-1 to convert to C++ inidces
  upstream.nn.list <- lapply(seq(1,vcount(net),1),function(x){neighbors(net,x,mode="out")-1}) #-1 to convert to C++ inidces
  downstream.nn.list <- lapply(seq(1,vcount(net),1),function(x){neighbors(net,x,mode="in")-1}) #-1 to convert to C++ inidces
  downstream.list <- lapply(get.shortest.paths(net,1,mode="out"),function(x){x-1}) #-1 to convert to C++ inidces
  output <- dyn_area_func(t.term, reps, nn.list, downstream.nn.list, upstream.nn.list, downstream.list, area, max.ext.vector, pr.col, pr.fme, pr.tc)
  occupancy.river[[f]] <- output[[1]]
  occupancy.river.sd[[f]] <- output[[2]]
  per.node.river[[f]] <- output[[3]]
  per.node.river.sd[[f]] <- output[[4]]
  ratio.river <- output[[5]]
  
  
  #Early Net
  net <- riv.net$early.net
  area <- riv.net$early.area
  #Build list of nearest neighbors in the network
  nn.list <- lapply(seq(1,vcount(net),1),function(x){neighbors(net,x,mode="all")-1}) #-1 to convert to C++ inidces
  upstream.nn.list <- lapply(seq(1,vcount(net),1),function(x){neighbors(net,x,mode="out")-1}) #-1 to convert to C++ inidces
  downstream.nn.list <- lapply(seq(1,vcount(net),1),function(x){neighbors(net,x,mode="in")-1}) #-1 to convert to C++ inidces
  downstream.list <- lapply(get.shortest.paths(net,1,mode="out"),function(x){x-1}) #-1 to convert to C++ inidces
  output <- dyn_area_func(t.term, reps, nn.list, downstream.nn.list, upstream.nn.list, downstream.list, area, max.ext.vector, pr.col, pr.fme, pr.tc)
  occupancy.random[[f]] <- output[[1]]
  occupancy.random.sd[[f]] <- output[[2]]
  per.node.random[[f]] <- output[[3]]
  per.node.random.sd[[f]] <- output[[4]]
  ratio.random <- output[[5]]
  
  #DLA Net
  net <- graph.edgelist(dla[[1]])
  area <- dla.area
  #Build list of nearest neighbors in the network
  nn.list <- lapply(seq(1,vcount(net),1),function(x){neighbors(net,x,mode="all")-1}) #-1 to convert to C++ inidces
  upstream.nn.list <- lapply(seq(1,vcount(net),1),function(x){neighbors(net,x,mode="out")-1}) #-1 to convert to C++ inidces
  downstream.nn.list <- lapply(seq(1,vcount(net),1),function(x){neighbors(net,x,mode="in")-1}) #-1 to convert to C++ inidces
  downstream.list <- lapply(get.shortest.paths(net,1,mode="out"),function(x){x-1}) #-1 to convert to C++ inidces
  output <- dyn_area_func(t.term, reps, nn.list, downstream.nn.list, upstream.nn.list, downstream.list, area, max.ext.vector, pr.col, pr.fme, pr.tc)
  occupancy.dla[[f]] <- output[[1]]
  occupancy.dla.sd[[f]] <- output[[2]]
  per.node.dla[[f]] <- output[[3]]
  per.node.dla.sd[[f]] <- output[[4]]
  ratio.dla <- output[[5]]
  
}

save.image("ec_area_sim_Net1024_T5K_R20_E005.RData")

pdf(file="/Users/jyeakel/Dropbox/PostDoc/2014_Empirical_Watersheds/Manuscript/Figure_ec_area_fme.pdf"
    ,width=6,height=5)

ratio <- apply(ratio.river,2,mean)
plot(ratio,occupancy.river[[1]],type="l",lwd=2,col="dodgerblue3",
     ylab="Relative abundance",xlab="e/c",xlim=c(0,1),xpd=F)
polygon(x=c(ratio,rev(ratio)),
        y=c(occupancy.river[[1]]+occupancy.river.sd[[1]],rev(
          unlist(lapply(occupancy.river[[1]]-occupancy.river.sd[[1]],function(x){max(x,0)})))),
        col="#1874CD50",border=NA)
if (F>1) {
  for (i in 2:F) {
    polygon(x=c(ratio,rev(ratio)),
            y=c(occupancy.river[[i]]+occupancy.river.sd[[i]],rev(
              unlist(lapply(occupancy.river[[i]]-occupancy.river.sd[[i]],function(x){max(x,0)})))),
            col="#1874CD50",border=NA)
    lines(ratio,occupancy.river[[i]],type="l",lwd=2,lty=i,col="dodgerblue3")
  }
}

ratio <- apply(ratio.random,2,mean)
for (i in 1:F) {
  polygon(x=c(ratio,rev(ratio)),
          y=c(occupancy.random[[i]]+occupancy.random.sd[[i]],rev(
            unlist(lapply(occupancy.random[[i]]-occupancy.random.sd[[i]],function(x){max(x,0)})))),
          col="#CD262650",border=NA)
  lines(ratio,occupancy.random[[i]],type="l",lwd=2,lty=i,col="firebrick3")
}

ratio <- apply(ratio.dla,2,mean)
for (i in 1:F) {
  polygon(x=c(ratio,rev(ratio)),
          y=c(occupancy.dla[[i]]+occupancy.dla.sd[[i]],rev(
            unlist(lapply(occupancy.dla[[i]]-occupancy.dla.sd[[i]],function(x){max(x,0)})))),
          col="#228B2250",border=NA)
  lines(ratio,occupancy.dla[[i]],type="l",lwd=2,lty=i,col="forestgreen")
}

dev.off()




pdf(file="/Users/jyeakel/Dropbox/PostDoc/2014_Empirical_Watersheds/Manuscript/Figure_ec_area_fme2.pdf"
    ,width=6,height=5)

ratio <- max.ext.vector/(pr.col)
plot(ratio,occupancy.river[[1]],type="l",lwd=2,col="dodgerblue3",
     ylab="Relative abundance",xlab="max_e/max_c",xlim=c(0,max(ratio)))
polygon(x=c(ratio,rev(ratio)),
        y=c(occupancy.river[[1]]+occupancy.river.sd[[1]],rev(
          unlist(lapply(occupancy.river[[1]]-occupancy.river.sd[[1]],function(x){max(x,0)})))),
        col="#1874CD50",border=NA)
if (F>1) {
  for (i in 2:F) {
    polygon(x=c(ratio,rev(ratio)),
            y=c(occupancy.river[[i]]+occupancy.river.sd[[i]],rev(
              unlist(lapply(occupancy.river[[i]]-occupancy.river.sd[[i]],function(x){max(x,0)})))),
            col="#1874CD50",border=NA)
    lines(ratio,occupancy.river[[i]],type="l",lwd=2,lty=i,col="dodgerblue3")
  }
}

for (i in 1:F) {
  polygon(x=c(ratio,rev(ratio)),
          y=c(occupancy.random[[i]]+occupancy.random.sd[[i]],rev(
            unlist(lapply(occupancy.random[[i]]-occupancy.random.sd[[i]],function(x){max(x,0)})))),
          col="#CD262650",border=NA)
  lines(ratio,occupancy.random[[i]],type="l",lwd=2,lty=i,col="firebrick3")
}

for (i in 1:F) {
  polygon(x=c(ratio,rev(ratio)),
          y=c(occupancy.dla[[i]]+occupancy.dla.sd[[i]],rev(
            unlist(lapply(occupancy.dla[[i]]-occupancy.dla.sd[[i]],function(x){max(x,0)})))),
          col="#228B2250",border=NA)
  lines(ratio,occupancy.dla[[i]],type="l",lwd=2,lty=i,col="forestgreen")
}

dev.off()


boxplot(ratio.river,names=max.ext.vector,boxwex=0.25,at=seq(1,length(max.ext.vector),1),col="dodgerblue3",outline=FALSE)
boxplot(ratio.random,names=max.ext.vector,boxwex=0.25,at=seq(1,length(max.ext.vector),1)-0.25,outline=FALSE,add=T,col="firebrick3")
boxplot(ratio.dla,names=max.ext.vector,boxwex=0.25,at=seq(1,length(max.ext.vector),1)+0.25,outline=FALSE,add=T,col="forestgreen")


plot(max.ext.vector,apply(ratio.river,2,sd),col="dodgerblue3",type='l',ylim=c(0,2),ylab="Mean e_i/c_i")
points(max.ext.vector,apply(ratio.random,2,sd),col="firebrick3",type='l')
points(max.ext.vector,apply(ratio.dla,2,sd),col="forestgreen",type='l')


library(vioplot)
cdat <- as.list(as.data.frame(ratio.river))
names(cdat)[1] <- "x"
do.call(vioplot,cdat)
vioplot(ratio.river,names=max.ext.vector,boxwex=0.25,at=seq(1,length(max.ext.vector),1),col="dodgerblue3",outline=F)


ratio.river.temp <- ratio.river[,-1]
library(reshape2)
mdat <- melt(ratio.river.temp)
library(ggplot2)
ggplot(mdat,aes(x=Var2,y=value,group=Var2))+geom_violin()


plot(ratio.river[,5],pch=16,col="dodgerblue3",ylim=c(0,0.1))
points(ratio.random[,5],pch=16,col="firebrick3")
points(ratio.dla[,5],pch=16,col="forestgreen")




#Fractal/Geometric Properties of networks
#####
#log P(A>a) vs. Log a
area <- riv.net$late.area
area <- riv.net$early.area
area <- dla.area

min_area <- min(area)
max_area <- max(area)
a_seq <- seq(1,max(area),1)
PrAa <- numeric(length(a_seq))
for (i in 1:length(a_seq)) {
  a <- a_seq[i]
  PrAa[i] <- length(which(area > a)) / length(area)
}
max_measure <- 35
lmodel <- lm(log(PrAa)[2:max_measure] ~ log(a_seq)[2:max_measure])
slope = round(lmodel[[1]][[2]],4)
plot(log(a_seq),log(PrAa),pch=16,ylab="Log Pr(A>a)",xlab="Log a",cex=0.5,main=slope)
abline(lmodel)

#Correlation dimension of DLA
pairwise.dis <- apply(dla[[1]],1,function(x){apply(dla[[1]],1,function(y){dist(rbind(x,y))})})
pairwise.dis <- c(pairwise.dis.m)
seq_dis <-seq(max(round(pairwise.dis,0)),1,-1)
count <- numeric(length(seq_dis))
D2 <- numeric(length(seq_dis))
tic <- 0
for (i in seq_dis) {
  tic <- tic + 1
  r <- i
  count[i] <- length(which(pairwise.dis <= r))
  D2[i] <- log(count[i])/log(r)
}
plot(seq_dis,D2,pch='.',ylim=c(0,3))


