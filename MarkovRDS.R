library(dplyr)
library(lazyeval)



atb = function(x, outcome = "HIV", blockVariable = NULL, Bcount, rawSamples = F, verbose = T, pretty = T){
  
  # x is a tibble or data.frame with the following columns:
  #  id, recruiter.id, network.size (or network.size.variable), "outcome" and "blockVariable"
  # outcome is the column name of the measurement for which we want to estimate the population mean. 
  # blockVariable is the column name of a feature (as.character) which pseudo-stratifies the sampling.  
  #    By "pseudo-stratifies", these variables should be highly auto-correlated in referrals.  
  #    Future code will help to select these features. 
  #    If it is not specified, then it is set to outcome.
  if(length(blockVariable)==0) blockVariable = outcome
  
  # Bcount is the number of bootstrap samples. 
  # rawSamples = T returns an (n x Bcount) matrix identifying each of the bootstrap samples.  
  #    The elements of the returned matrix are row indices corresponding to rows of x.  
  # verbose = F should be used in simulation experiments. 
  # pretty = F returns the Bcount-many bootstrapped point estimates.
  
  #  TODO:  automatically fix the "chronology issue" described next.
  #  TODO:  create a procedure to help pick the blockVariable
  #  TODO:  allow atb to use other estimators (not just Volz-Heck) by taking a function as an argument
  #  TODO:  should we allow outcome to be a vector? 
  
  # This function presumes that the order of the nodes in x is chronological.
  #   More specifically, for all pairs where i refers j, the row for i must appear in x before the row for j.  
  # here is a check:
  
  parentRowID = match(x$recruiter.id, x$id)
  parentRowID[is.na(parentRowID)] = 0
  badChildren = which(parentRowID > 1:length(x$id))
  if( length(badChildren) > 0){
    cat("This function only works if all recruiter.id's appear first as id's, \n")
    cat("  this will happen if (for example) the data is not arranged chronologically. \n")
    cat("  Unfortunately, there appear to be some problems with the ordering. \n")
    cat("The currently identified problem rows are: \n \n      ")
    cat(c(badChildren, parentRowID[badChildren]))
    cat("\n \n ")
    cat("Here are those rows in the original data: \n \n      ")
    print(x[c(badChildren, parentRowID[badChildren]),])
    return(NA)
  }
  
  
  
  
  # check to see if any id's are repeating:
  tab = table(x$id)
  if(max(tab)>1){
    cat("\n")
    cat("**These id's are repeating:** \n \n")
    tmp  = min(c(5, ncol(x)))
    bad = which(x$id %in% x$id[anyDuplicated(x$id)])
    print(x[bad,])
    cat("
This is a problem because we don't know who did the recruiting.
If any repeating id's did not do any recruiting, then assign them arbitrary new id's and re-run.
        ")
    return(NULL)
  }
  
  
  
  #  A key component of atb is constructing a matrix Q that counts the number of transitions 
  #   from one block (i.e. strata) to the next.  The function Q2rds computes Q and assigns it to attr(x,'Q')
  if(length(attr(x,'blockVariable')) ==0) x = Q2rds(x, blockVariable)
  blockVariable = attr(x,'blockVariable')
  
  Qsbm = attr(x,'Q')
  Psbm = diag(1/rowSums(Qsbm))%*%Qsbm 
  rownames(Psbm) = rownames(Qsbm)
  
  z = as.matrix(x[,blockVariable])  
  
  # ensure the outcome variable takes multiple values:
  if(length(unique(as.matrix(x[,outcome]))) < 2){
    cat("Outcome variable only takes one value. Addressing such situations needs an entirely different approach.")
    return(NA)
  }
  
  #  sometimes the "degree" of each participant is in "network.size.variable".  othertimes, it is in "network.size".  
  if(suppressWarnings(length(x$network.size.variable))>0) degs = x[,"network.size.variable"]
  if(suppressWarnings(length(x$network.size.variable))==0) degs = x$network.size
  
  # refs[i] counts the number of times that row i refers participants.
  tab=table(x$recruiter.id)
  refs = tab[match(x$id, as.numeric(names(tab)))]
  refs[is.na(refs)]=0

  n = nrow(x)
  
  # offSp[[i]] gives the rows of x that i refers.
  offSp = lapply(x$id, function(s)  return(which(s ==x$recruiter.id)))
  
  
  # There are two steps to the sampling.  First, sample the block label (z) for each node in the bootstrap sample.
  #   b.z will hold those values for all n x Bcount observations.  
  b.z = matrix(NA, nrow = n, ncol = Bcount)
  
  # Second, given a block label (z) from b.z, select an individual from the original sample.
  #   Select that person uniformly from the set of participants in the group z. 
  #   This yields a matrix bsIndices:
  
  bsIndices = matrix(NA, nrow = n, ncol = Bcount)
  
  
  # to perform the first step (sampling the z labels), there are three sub-steps.  
  # A) sample the "mother node" (explained below). 
  # B) sample the seed nodes
  # C) sample the rest of the nodes. 
  
  # In order to sample the labels for the "n observations", first sample the node type of the "mother node" uniformly from all types. 
  
  # #### A) MOTHER NODE ######
  # If there is only one seed node, then this is irrelevant (but we do it anyways).
  # If there is more than one seed node, this is an important step and it is a bit strange.  
  # To ensure that seed nodes are not treated as iid samples from the stationary distribution (something which is typically not possible for RDS),
  #   this code presumes that all seeds are referred by a "new node" which is best considered as "the researcher team"
  #   This ensures that seed nodes are dependent samples.... 
  
  # Perhaps it is 'more natural' to re-sample the seeds from the original seeds.
  #  However, such a procedure would likely have a smaller variance.  
  #  this smaller variance would come from the fact that the initial distribution is heavily biased (based upon the social network position of the research team).  
  #  Imagine if the initial distribution was entirely on one class... 
  #  this would result in a smaller variance because it does not consider the possibility of all seed nodes coming from a different class. 
  
  
  # the mother node's type will be ignored in the computation of the estimator.  It is only here to make the seeds dependent. 
  
  motherType = sample(z,Bcount, replace = T)  
  
  
  
  # #### B) SEED NODES ######
  #  sample each seed as if it is one transition away from the Mother node. 
  
  # here, seeds are defined to be the rows of x for which their recruiter.id was not a participant.  
  #    (There are certainly ways for this step to fail if the data is messy.)
  seedRows = which(!(x$recruiter.id %in%  x$id))  
  
  # because this step could fail, make a print out for 
  if(verbose){
    cat("\n")
    cat("The seed rows are identified to be indices: \n")
    cat(seedRows)
    cat("\n \n \n")
    tmp = min(5, ncol(x))
    cat("They look like this: \n")
    print(x[seedRows, 1:tmp])
    cat("\n If you believe that this is incorrect, please ensure that which(!(x$recruiter.id %in%  x$id)) identifies your seeds. \n \n")
    
    cat("The matrix of transitions between blocks is:\n")
    print(round(Qsbm))
    cat("\n")
    
    cat("HIV prevalence in each block is:\n")
    x$blockLabel = as.matrix(x[,blockVariable]) %>% as.vector
    x$outcome.local =  as.matrix(x[,outcome]) %>% as.vector
    print(x %>% group_by(blockLabel) %>% summarize(blockSize = n(), AverageOutcome = mean(outcome.local)))
    
  }
  
  un = sort(unique(z))
  
  # here is where we sample the seed nodes. 
  # b.z has Bcount columns.   we loop over the unique values of z and fill in the appropriate columns:  b.z[seedRows,motherType ==un[b]] 
  #   For each class b, select the columns of b.z where the mother type is that class.
  #   Then, resample the seedRows to have the appropriate distribution (specified by Psbm[motherType,]):
  
  for(b in 1:length(un)) b.z[seedRows,motherType ==un[b]] = sample(
    x = colnames(Psbm), 
    size = length(seedRows), 
    replace = T,
    prob = Psbm[which(colnames(Psbm) == un[b]),])
  
  
  # C) sample the rest of the nodes. 
  
  # the rest of the sampling of b.z is very similar to sampling the seeds above,
  #   but we add an outer loop through the rows of x for which refs[i]>0 (i.e. they referred people).
  #   For each row i that refers participants, we fill in the values of  b.z[offSp[[i]],b.z[i,] ==un[b]].
  
  # The sequence of this loop relies on the fact that if i refers j, then j>i.   TODO: fix this.  
  
  i = 1
  for(i in which(refs>0)){
    offSpring = offSp[[i]]  # these are the rows of b.z that we are going to sample. These are offspring of row i.
    
    #  for each block/class/stratum b, select the columns (i.e. bootstrap samples) such that 
    #   the class of row i is b.  In code:  b.z[i,] ==un[b]
    #   Then, we are sampling rows "offSpring" from the distribution Psbm[b,].
    for(b in 1:length(un)) b.z[offSpring,b.z[i,] ==un[b]] = sample(
      x = colnames(Psbm), 
      size = length(offSpring)*sum(b.z[i,] ==un[b]), 
      replace = T,
      prob = Psbm[which(colnames(Psbm) == un[b]),])
  }
  
  
  # Now that we have sampled all the block/class/stratum z in the matrix b.z, we can now sample individual rows of x, 
  #    conditional on that block assignment:
  
  for(u in rownames(Psbm)){
    # for each of the node type in b.z, sample an observation from x with that type.
    these = which(b.z == u)
    bsIndices[these] = sample(which(x[,blockVariable] == u), size = length(these),replace = T)
    
  }
  
  if(rawSamples) return(bsIndices)
  
  # compute the Volz-Heckathorn point estimate for each bootstrap sample, i.e. using the samples from each column of bsIndices:
  
  b.vh= apply(bsIndices, 2, function(samp) return(vh(z = as.matrix(x[samp,outcome]), degs[samp])))
  if(!pretty)return(b.vh)
  
  if(pretty){
    cat("The Volz-Heckathorn estimate is:\n")
    cat(round(vh(as.matrix(x[,outcome]), degs),2))
    cat("\n \n")
    cat("The 90% confidence interval is:\n")
    print(round(quantile(b.vh, probs = c(.05,.95)),2))
    
    
    # print("A qq-plot examines whether the estimator appears approximately normally distributed.")
    # print("Do not be scared by this plot and feel free to ignore it.")
    # print("The qq-plot can be used to examine why the confidence interval is not symmetric.")
    # print("If it is not symmetric, fear not; the bootstrap confidence interval above is still valid.")
    # qqnorm(scale(b.vh))
    # abline(0,1)
  }
  
}






vh = function(z,degs){
  # computes the vh estimator.
  degH = 1/mean(1/degs)
  return(mean(z/degs)*degH)
}


Q2rds = function(x, blockVariable = "z"){
  
  #   if(!is.rds.data.frame(x)){
  #     print("Error: Please pass treeboot an rds.data.frame.")
  #     return(NA)
  #   }
  #   
  id = as.numeric(as.matrix(x[,"id"]))
  
  crosswalk = function(bigNum){
    return(match(bigNum, id))
  }
  
  pid = suppressWarnings(as.numeric(as.matrix(x[,"recruiter.id"])))
  z = as.character(as.matrix(x[,blockVariable]))
  if(length(unique(z)) == 1) return(NA)
  
  
  seedRows = which(!(x$recruiter.id %in%  x$id))
  
  pz = as.character(z[crosswalk(pid)])
  pz[seedRows] = "seed"
  
  
  n = length(id)
  
  # the follwoing makes the plug in estimate of Qsbm, 
  # Qsbm = table(as.data.frame(cbind(pz,z)[-1,]))   # add .01 to each element.
  
  cases = unique(c(pz,z))
  cases = cases[complete.cases(cases)]
  cases = setdiff(cases,"seed")
  
  nk =length(cases)
  columnIndex = matrix(1:nk, ncol = nk, nrow = nk, byrow = T)
  rowIndex  = matrix(1:nk, ncol = nk, nrow = nk)
  nullEdges = cbind(cases[columnIndex], cases[rowIndex])
  refEdges = cbind(pz,z)[-seedRows,]
  Qsbm = table(as.data.frame(rbind(nullEdges,refEdges))) - matrix(1, nrow = nk, ncol = nk)  # this -matrix(1,...) removes the nullEdges
  Qsbm = Qsbm + matrix(.05, nrow = nk, ncol = nk)  # this adds .05 to each element of Qsbm
  attr(x,'Q') = Qsbm
  attr(x,'blockVariable') = blockVariable
  return(x)
}








atb_power = function(n,m, A, seeds = .9, pi = NULL,ybar = NULL, Bcount= 1000, for.simulation =F, rawOutput=F){
  
  # n                 # sample size
  # m                 # average number of coupons returned by participants (will be simulated as Poisson distributed)
  # A                 # this is a K x K matrix where the row normalized version creates a Markov transition matrix on 1 ... K
  #                        alternatively, this can be a single number which is the probability of a referral within the same class,
  #                        with this parameterization, all other probabilities are uniform. 
  # seeds             # if(seeds >= 1)
  #                        then this computes (1 - the probability of chain death), 
  #                        i.e. a proxy for reaching the target sample size before chain death. 
  #                        This computation presumes a Galton-Watson referral process with Poisson(m) offspring distribution.
  #                   # if(seeds<1)
  #                        then this computes the number of initial referrals needed so that (1 - the probability of chain death) = seeds. 
  #                   # then, the referral tree is initialized with this many seeds (i.e. wave 0, but connected by a mother node... see discussion in atb above.)
  # pi                # this is a K vector. Element j is the proportion of individuals in group j.
  #                     if not specified, it is set to the stationary distribution.
  # ybar              # this is a K vector. Element j is the proportion of individuals in block j that are HIV+. 
  #                     if not specified it is set to c(0,1).
  
  if(length(ybar) == 0){
    ybar  = c(0,1)
    K = 2
    if(nrow(as.matrix(A))>2) cat("If A is larger than 2x2, you must specify ybar (the block average of y)")
  }
  if(length(K) ==0)  K = max(nrow(as.matrix(A)), length(pi), length(ybar))
  if(K ==1){
    cat("somehow, K = 1.  This is an error. \n")
    return(NA)
  }
  
  
  if(length(A) == 1){
    pIN = A
    pOUT = (1 - pIN)/(K-1)
    A = diag(rep(pIN - pOUT, K))
    A = A + pOUT
  }
  
  P = diag(1/rowSums(A))%*%A
  ei=eigen(t(P),symmetric = F)
  
  if(length(pi) == 0) pi = Re(ei$vec[,1])
  
  if(rawOutput) for.simulation=T
  
  # if the eigenvalues are complex (P not reversible),
  #   then this is a conjecture that will (for the time being) be conservative:
  criticalThreshold = 1/abs(ei$val[2])^2  
  if(!for.simulation) if(m > criticalThreshold){
    cat("
        Based upon the specification for A, the critical threshold for m is at: ")
    cat(round(criticalThreshold,1))
    cat(".  
        You have specified a value of m which exceeds this value; this is bad.
        This code will still compute a width of uncertainty.
        However, it will be very wide (relative to the specified sample size).
        Moreover, increasing n is unlikely to help much.  In order to decrease the 'standard error' by 50%,
        you will require ~")
    gamma = log(abs(ei$val[2]), m)   
    cat(round(.5^(1/gamma)))
    cat("x as many samples; because the calculation of this number presumes 
        an infinite population, this is potentially an absurd number. \n  
        One way to avoid this problem is to reduce m,
        or for small populations, increase m and aim for a census.
        Alternatively, use an estimator which is less sensitive to this threshold (forthcoming).
        Alternatively, you could consider a novel design technique,
        such as the one described in https://arxiv.org/abs/1606.00387
        If you are sure of your specification of A and m, I (karl) would love to talk more; 
        feel free to send me an electronic mail message:  my gmail username is karlrohe.\n \n ")
  }
  
  # extinction probability for poisson(lambda) offspring satifies, exp(-lambda(1-pExtinct))=pExtinct
  #   so, given pExtinct, lambda = log(pExtinct)(pExtinct-1)
  extinctionProb = seq(0,1, by = .001)[-c(1,1001)]  # this is extinction probability *for one seed*
  mExtinct = log(extinctionProb)/(extinctionProb-1) 
  # fixing m, find the extinction Prob.  
  thisIndex = which(mExtinct<m)[1]
  q= extinctionProb[thisIndex]
  

  
  if(seeds > 1){
    if(!for.simulation) cat("For the specified value of m =", m, "and the number of seeds =", seeds,"
      the probability of chain-death is:")
    prob = round(q^seeds,2)
    if(!for.simulation)  if(prob<.05) cat(" less than 5%")
    if(!for.simulation)  if(prob>.05) cat(round(q^seeds,2))
  }
  

  if(seeds < 1){
    if(!for.simulation)  cat("\nFor the specified value of m =", m, "and the desired probability of chain-death <",1-seeds,",
you should select at least this many seeds:", ceiling(log(1-seeds, q)))
    seeds = ceiling(log(1-seeds, q))
  }
  
  if(!for.simulation)  cat("
      This probability calculation is based upon 
      independent Poisson referrals (i.e. wrong) 
      and should only be used as a rough guide.
      ")
  
  Tr = treeSamp(lambda=m, n, seeds)
  
  b.z = matrix(NA, n,Bcount)
  b.z[1,] = sample(K, size = Bcount, replace = T, prob = pi)
  
  tmp = as_edgelist(Tr)
  id = c(1,tmp[,2])
  recruiter.id = c(NA, tmp[,1])
  
  
  # offSp[[i]] gives the rows of x that i refers.
  offSp = lapply(id, function(s)  return(which(s ==recruiter.id)))
  refs = do.call(c,lapply(offSp, length))  # how many referrals from each node
  
  i = 1
  for(i in which(refs>0)){
    offSpring = offSp[[i]]  # these are the rows of b.z that we are going to sample. These are offspring of row i.
    
    #  for each block/class/stratum b, select the columns (i.e. bootstrap samples) such that 
    #   the class of row i is b.  In code:  b.z[i,] ==un[b]
    #   Then, we are sampling rows "offSpring" from the distribution Psbm[b,].
    for(b in 1:K) b.z[offSpring,b.z[i,] == b] = sample(
      x = 1:K, 
      size = length(offSpring)*sum(b.z[i,] == b), 
      replace = T,
      prob = P[b,])
  }
  
  #if pi is not the stationary distribution of P
  #   Then, we need to adjust with inverse-probability weighting. 
  #   When the process is reversible and the degrees are homogeneous within blocks, this is equivalent to Volz-Heckathorn.
  stationary = Re(ei$vec[,1])  # by Perron-Frobeneius, this should be real valued ... any Im is machine error.
  stationary = stationary/sum(stationary)
  pi = pi / sum(pi)
  
  # there are pi individuals in the population,
  #   but they are selected with probability stationary.  
  # So, we need to adjust by pi/stationary.
  ypi = ybar * pi/stationary
  
  b.mu = apply(b.z,2, function(index) return(mean(ypi[index])))
  if(rawOutput) return(b.mu)
  
  interval = quantile(b.mu, c(.05,.9))
  # dat = c(diff(interval), interval)
  dat = diff(interval)
  names(dat)[1]  = "Width of 90% confidence interval"
  
  if(!for.simulation)  cat("\n",names(dat), "\n" , round(dat,2), "\n \n")
  if(for.simulation) return(dat)
}





