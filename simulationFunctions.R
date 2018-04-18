# this calls functions to simulate a graph from fastRG:

source("https://raw.githubusercontent.com/karlrohe/fastRG/master/fastRDPG.R")
library(tibble)

getGraph = function(N=100000, K=5){
  # simulates a degree corrected stochastic blockmodel with N nodes and K blocks.
  #  returns an igraph
  pi = rexp(K) +1
  pi = pi/sum(pi) 
  pi = -sort(-pi)
  B = matrix(rexp(K^2)+1, nrow=K)
  diag(B) = diag(B)+ mean(B)*K
  
  # these are the edges in the social graph
  graphData = dcsbm(rgamma(N,shape = 2,scale = .4), pi,B,avgDeg = 10,returnParameters = T, returnEdgeList = T)  
  block = apply(graphData$X,1, which.max)  # these are the block labels for the graph
  G = graph_from_edgelist(graphData$el,directed = F)
  positiveDegreeNodes = 1:N %in% (as.vector(graphData$el))
  G = delete.vertices(G, V(G)[degree(G) ==0])
  V(G)$block = as.character(block[positiveDegreeNodes])
  V(G)$HIV = V(G)$block == "5"
  return(G)
}





# lambda = 2
# n = 1000
# here are functions to simulate an RDS:

simpleRDS=function(G, lambda, n, inflation = 3){
  # Given an igraph G, a poisson rate parameter lambda, and a sample size n, this function simulates an RDS process.
  #   inflation needs to be big enough such that we can get enough samples.  In sparse graphs G, it needs to be bigger.
  #   when inflation is larger, the first step takes more time:
  Tr = treeSamp(lambda, n*inflation)  # this function samples a tree as an igraph.
  thisSeed = sample(V(G), 1)[1]  # this samples a seed node from G. 
  
  # this function samples without replacement, filling the nodes of Tr with node id's from G.  
  # Tr has 3*n nodes (instead of n) because sometimes the sampling cannot be peformed...
  #    a node must refer 3 friends in Tr, but that node doesn't have that many friends (who have not yet participated)
  # X will only have 3*n elements, many of which are NA. 
  X = RDSwithTree(G,Tr, seed = thisSeed,repl=F, n)  
  TrSub = induced_subgraph(Tr, which(complete.cases(X))[1:n])  # this removes the nodes of Tr that were not used in X.  
  Xsub = X[complete.cases(X)][1:n]  # discard the NA's.  
  el = get.edgelist(TrSub)  # this records which elements of X referred which other elements. 
  
  return(  # this constructs the table "x" that is used in atb.  
    tibble(
      id = c(as.numeric(thisSeed), Xsub[el[,2]]),
      recruiter.id  = c(-1, Xsub[el[,1]]),
      network.size = degree(graph = G,v = Xsub),
      block = V(G)$block[Xsub],
      HIV = V(G)$HIV[Xsub]
    ) 
  )
}


treeSamp = function(lambda, sampSizei, seeds = NULL){
  Gi = graph(c(1,1))
  active = 1
  n=1
  # if the number of seeds is specified,
  if(length(seeds)>0){
    newActive = c()
    addnum = ceiling(seeds)
    if(addnum > 0) {
      newids = n + 1:addnum
      n = max(newids)
      newActive = c(newActive,newids)
      Gi = add.vertices(graph = Gi,nv = addnum)
      newedges = t(cbind(rep(1, addnum), newids))
      Gi = add.edges(graph = Gi,edges = newedges)
    }
    active = newActive
    newActive = c()
  }
  
  while(n<sampSizei){
    if(length(active) == 0){  # this conditions on survival.  
      return(treeSamp(lambda, sampSizei, seeds))
    }
    newActive = c()
    for(j in active){
      addnum = rpois(n = 1, lambda= lambda)
      if(addnum > 0) {
        newids = n + 1:addnum
        n = max(newids)
        newActive = c(newActive,newids)
        Gi = add.vertices(graph = Gi,nv = addnum)
        newedges = t(cbind(rep(j, addnum), newids))
        Gi = add.edges(graph = Gi,edges = newedges)
      }
    }
    active = newActive
    newActive = c()
  }
  Gi = delete.edges(graph=Gi,edges = c(1,1))
  if(n>sampSizei) Gi = delete.vertices(graph=Gi, v = (sampSizei+1):n)
  return(Gi)
}



RDSwithTree = function(Gi,tree, seed = c(),repl=F, sampSizei){
  # Gi is social graph
  # tree is referral topology
  # seed is ego to start
  # set (repl = T) if sampling is with replacement.
  if(is.vector(tree)){
    offDist = tree
    tree = treeSamp(roff,sampSizei)
  }
  if(is.function(tree)) tree = treeSamp(tree,sampSizei)
  
  n = length(V(tree))
  X = rep(NA,n)  # this records the sample.  the ith element will contain the node id (in G) for node i in the tree.
  if(length(seed) ==0) seed = sample(1:length(V(Gi)),size = 1,prob = degree(Gi))
  
  X[1] = seed
  for(i in 1:length(V(tree))){
    # some nodes in tree do not have observations in X.  
    #ensure that X[i] is filled.
    if(!is.na(X[i])){                                           # if i is observed in X
      children = neighbors(graph = tree,i,mode = "out")         # find the children in the tree
      if(length(children)>0){                                   # if this is non-empty,
        offsp = sampChildren(ego = X[i],X,Gi, 
                             numRefs = length(children), repl)  # sample the neighborhood.
        if(length(offsp) >0){                                   # if the sample is nonempty
          children = children[1:length(offsp)]                  # fill the appropriate number of space in X
          X[children] = offsp
        }
      }
    }
    
    if(sum(complete.cases(X))>sampSizei) return(X)
  }
  return(X)
}



sampChildren = function(ego, X, Gi, numRefs, repl){
  # this function is used inside RDSwithTree
  neigh = neighbors(graph = Gi,v = ego)              # find the neighbors
  if(repl) return(sample(neigh,numRefs,replace = T))# if it is with replacement, just do it.
  if(!repl){                                        # otherwise,
    neigh = setdiff(neigh, X[complete.cases(X)])    # remove previously sampled nodes,
    numRefs = min(numRefs, length(neigh)) 
    return(sample(neigh, numRefs, replace = F))
  }
}





