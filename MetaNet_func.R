MetaNet_func <- (net,scale) {
  riv.net <- net[[1]]
  riv.area <- net[[2]]
  
  #Define features of watershed net
  tribs <- which(degree(riv.net) == 1)
  conf <- which(riv.area == max(riv.area))
  tribs <- tribs[-conf]
  l.tribs <- length(tribs)
  paths <- shortest.paths(riv.net,conf,mode="out")
  longest.path <- which(paths==max(paths))
  lpath <- paths[longest.path[1]]
  
  rem.nodes <- seq(1,vcount(riv.net),1)
  
  for (i in seq(lpath,1,-1)) {
    #find nodes at this order
    lnodes <- which(paths == i)
    for (j in 1:length(lnodes)) {
      node <- lnodes[j]
      #Determine upriver nodes
      upriver.nodes <- neighbors(riv.net,node,mode="out")
      if (length(upriver.nodes) > 0) {
        #Which of these nodes have not been included in the meta-network?
        include.nodes <- upriver.nodes[which(upriver.nodes %in% rem.nodes)]
        #Calculate cumulative area
        cum.area <- sum(c(riv.area[include.nodes]),riv.area[node])
        if (cum.area > scale) {
          #Delete these nodes from the remaining nodes list
          rem.nodes <- rem.nodes[-c(include.nodes,node)]
          #Add single node to the metapopulation network
          
          
          
        }
      }
    }
    
  }
  
}