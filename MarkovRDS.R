library(dplyr)
library(lazyeval)
library(tibble)
library(Matrix)

# This file contains code to perform 
# 1) the assisted treebootstrap, as described in https://arxiv.org/abs/1505.05461
# 2) feasible GLS, as described in https://arxiv.org/abs/1708.04999



#######################################################################
########  functions to perform theassisted treebootstrap for RDS ######
#######################################################################


atb = function(x, outcome = "HIV", blockVariable = NULL, Bcount, rawSamples = F, verbose = T, pretty = T, glsBoot = F, symmetric = F){
  
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
  # symmetric = T symmetrizes Q (the estimate of B) before doing the resampling.  
  #               this is useful for testing reversibility (samples under null = reversible)
  
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
  if(symmetric) Qsbm = Qsbm + t(Qsbm)
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
  
  
  if(!glsBoot){# compute the Volz-Heckathorn point estimate for each bootstrap sample, i.e. using the samples from each column of bsIndices:
    
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
  if(glsBoot){
    b.gls = apply(bsIndices,2, bootGLSfunction, x = dat, outcome = "HIV", blockVariable = "block", deg = T)
    
    
    
    # b.vh= apply(bsIndices, 2, function(samp) return(vh(z = as.matrix(x[samp,outcome]), degs[samp])))
    if(!pretty)return(b.gls)
    
    if(pretty){
      cat("The sbmGLS estimate is:\n")
      cat(
        round(
          sbmGLS(x, outcome, blockVariable,deg= T), 2))
      
      cat("\n \n")
      cat("The 90% confidence interval is:\n")
      print(round(quantile(b.gls, probs = c(.05,.95)),2))
      
      
      # print("A qq-plot examines whether the estimator appears approximately normally distributed.")
      # print("Do not be scared by this plot and feel free to ignore it.")
      # print("The qq-plot can be used to examine why the confidence interval is not symmetric.")
      # print("If it is not symmetric, fear not; the bootstrap confidence interval above is still valid.")
      # qqnorm(scale(b.vh))
      # abline(0,1)
    }
  }
  
}


vh = function(z,degs){
  # computes the vh estimator.
  degH = 1/mean(1/degs)
  return(mean(z/degs)*degH)
}

bootGLSfunction = function(bootID, x, outcome, blockVariable, deg){
  # calls sbmGLS, for use inside of atb.
  x[,-(1:2)]= x[bootID,-(1:2)]
  return(sbmGLS(x=x, outcome=outcome, blockVariable=blockVariable, deg = deg))
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







#######################################################################
########  functions to perform feasible GLS estimation for RDS ########
#######################################################################

# There are three types of estimators below, 
# Delta, auto, and SBM. 
# Delta is not recommended.  It is included for completeness. 
# Each of these estimators is described in 
# "Generalized Least Squares can overcome the critical threshold 
#    in Respondent-Driven Sampling" by Roch and Rohe 2017. 






######
### The basic structure of fGLS is to:
#  1)  (somehow) estimate \lambda_2, <y,f_2> and any other (lam,coef) pairs. (autocorrelation, delta-first diff, sbm)
#  2)  use these to construct the auto-covariance function 
#          - this is called acHat below
#  3)  with the referral tree and the auto-covariance function, construct the nxn covariance matrix 
#          - this is called makeCov below
#  4)  with the covariance matrix, compute the GLS estimator by first solving system of equations:  Sigma x = 1.  
#          - this is called gls below



acHat = function(coefs, lams, diam){
  #  given the (estimated) eigen-stuff (i.e. coefs and lams) and 
  #      the maximum distance between two nodes in the referral tree (diam)
  #  this function return the estimated auto-covariance function.
  #  The output est[t+1] gives the covariance between Y_sigma and Y_tau 
  #   for sigma and tau at distance t appart in the referral tree.
  
  # given the eigenvalues (lams) and the <y,f_\ell>^2 coefs,
  #  compute the covariance as a function of t.  
  #  returns a vector for which this is est[t+1]
  #   the length of the vector should be (diam),
  #   the diameter of the tree + 1.
  
  if(abs(lams[1]-1)<10^(-14)) lams = lams[-1]
  if(length(lams)!=length(coefs)){
    print("ERROR:  In function acHat, do not pass the first eigenvalue. Then, ensure that lams and coefs are vectors of the same length.")
    return(NA)
  }
  K = length(lams)
  mat = matrix(NA, nrow = diam+1, ncol = K)
  for(j in 1:K) mat[,j] = lams[j]^(0:diam)
  est = mat %*% coefs
  return(est)
}



makeCov = function(Tr, ac, ridge = 0){
  # ac is the (estimated) auto-covariance function; e.g.  use acHat
  # Tr is the tree structure.
  # 
  # should be used as follows:
  #  lams, coefs estimated via first difference, auto correlation, or SBM.
  #  makeCov(Tr, acHat(coefs, lams, diameter(Tr, directed = F)))
  #  then, this function plays an essential role in compute the GLS estimator.  
  
  dmat = shortest.paths(Tr,mode = "all") + 1 # this is an igraph function that finds pairwise distances in the tree.  for indexing... +1.
  n = nrow(dmat)
  Cov = diag(0, n)
  
  for(t in 1:length(ac)) Cov[which(dmat==t)] = ac[t]
  diag(Cov) = (1+ridge) * diag(Cov)
  return(Cov)
}


gls = function(Tr = c(), ac = c(), Y, Cov = c()){
  # after estimating lams, coefs via first difference, auto correlation, or SBM... 
  # this can estimate the mean via:
  #  gls(Y= Yvh, Cov= makeCov(Tr, acHat(coefs, lams, diameter(Tr, directed = F))))
  # could 
  if(sum(is.na(Y))==0) if(length(ac)>0){
    if(sd(ac[-length(ac)]/ac[-1])<10^-8){  # if the Rank-2 model.  
      lam = ac[2]/ac[1]
      x = (1- lam*(degree(Tr,mode = "all") - 1))/(1 + lam)
      return(1/(sum(x))*as.vector(x%*%Y))
    }
  }
  if(length(Cov) ==0) Cov = makeCov(Tr, ac)  
  x = solve(a = Cov[complete.cases(Y),complete.cases(Y)],b = rep(1, length(Y[complete.cases(Y)])))
  return(1/(sum(x))*as.vector(x%*%Y[complete.cases(Y)]))
}







################
#### DELTA ####
################
# This estimator is not that good.  


diffGLS = function(Tr, Y){
  # this adds laplace smoothing in this first line:
  n = sum(complete.cases(Y))
  lam = (d2(Tr,Y) - d1(Tr,Y))/(d1(Tr,Y)+1/sqrt(n))  # the value of 1/sqrt(n) is appropriate for y \in {0,1}
  ac = acHat(coefs = 1,lam,diameter(Tr, directed = F))  
  return(c(gls(Tr = Tr, ac = ac, Y = Y), lam))
}


d1 = function(Tr, Y){
  # computes the first difference, squared, over the tree.
  eg = ends(Tr, 1:length(E(Tr)))
  # this fixes the indexing, so that eg doesn't contain node id's, but rather row # in x.
  eg = cbind(match(as.numeric(eg[,1]), as.vector(V(Tr)$name)),match(as.numeric(eg[,2]), as.vector(V(Tr)$name)))
  
  return(mean((Y[eg[,1]] - Y[eg[,2]])^2, na.rm=T))
}


d2 = function(Tr,Y){
  # computes the second difference, squared, over the tree.
  at = get.adjacency(as.undirected(Tr)) 
  a2 = at%*%at
  a2 = a2 - diag(diag(a2))
  eg = which(a2!=0, arr.ind = T)
  return(mean((Y[eg[,1]] - Y[eg[,2]])^2, na.rm = T))
}



################
#### auto ####
################




# estimate the autocovariance of adjacent entries in the tree.
autoInnerRDS = function(Tr, Y, mu = c()){
  # given a tree Tr and values Y, compute auto correlation between Y_sigma  and  Y_sigma', 
  #  where sigma' is the parent of sigma.  
  # if given a previous estimate of mu, pass this.  otherwise, uses sample average.  
  if(length(mu)==0) mu = mean(Y, na.rm=T)
  Yc = as.vector(scale(Y - mu, center = F))
  eg = ends(Tr, 1:length(E(Tr)))
  # fix the indexing (as done elsewhere in this document)
  eg = cbind(match(as.numeric(eg[,1]), as.vector(V(Tr)$name)),match(as.numeric(eg[,2]), as.vector(V(Tr)$name)))
  return(mean(Yc[eg[,1]] * Yc[eg[,2]], na.rm=T))
}



autoGLS = function(Tr,Y, lamBounds = c(.1,.95)){
  # this code computes the results for the AUTO estimator of \lambda. 
  #  it allows for "non-convergent" iterative-sequences of fGLS estimators (estimate mu, then \lambda, then mu, then \lambda).  
  
  # returns muHat, lamHat, if(converged), 
  if(sd(Y, na.rm= T)==0) return(mean(Y))
  
  lamIter = autoInnerRDS(Tr,Y)
  delt = 100
  while(delt > .0001 & lamIter < .99 & lamIter > 0){
    ac = acHat(coefs = 1,lamIter,diameter(Tr, directed = F))  
    
    # the one-dimensional GLS estimator has a closed form solution that is much faster to compute.  
    muHatGLS= gls(Tr = Tr, ac = ac, Y = Y)  
    lamOld = lamIter
    lamIter = autoInnerRDS(Tr,Y, mu = muHatGLS)
    
    delt = abs(lamOld  -lamIter)
  }
  if(delt < .01 & lamIter < .99 & lamIter > 0) return(c(muHatGLS, lamIter, T))
  
  # otherwise, do the iterative search.
  
  
  lamSeq = seq(lamBounds[1],lamBounds[2], by = .01)  # the resolution could be increased... by = .001 would slow the calculation.
  
  len = length(lamSeq)
  lamGradient = cbind(lamSeq, NA)
  muHatGLS = rep(NA, len)
  
  for(i in 1:len){
    
    # this uses coefs =1 because in this "one dimensional" model muHatGLS does not change with coef.
    ac = acHat(coefs = 1,lamSeq[i],diameter(Tr, directed = F))  
    
    # the one-dimensional GLS estimator has a closed form solution that is much faster to compute.  
    muHatGLS[i]= gls(Tr = Tr, ac = ac, Y = Y)  
    lamGradient[i,2] = autoInnerRDS(Tr, scale(Y- muHatGLS[i], center =F))  # this is the AUTO estimator for \lambda
    
    # sig2 = mean((Y-muHatGLS[i])^2)
    # lamGradient[i,3]= 1 - d1(Tr,Y)/(2*sig2)   # this is the DIFF estimator for \lambda
  }
  
  # this finds the "fixed points", even when they don't exist...
  gradient = lamGradient[,-1]-lamGradient[,1]
  # fixed points have gradient = 0; this just finds the minimum in abs.
  fixedPoint = which.min(abs(gradient))
  
  
  results = rep(F, 8)
  #   names(results) =c(
  #     "muHatAuto", "lamHatAuto", "AutoConverge")
  
  results[1] = muHatGLS[fixedPoint]
  results[2] =  lamGradient[fixedPoint,1]
  
  # it is possible that the code returns an estimate of \lambda at the boundary of the region investigated.
  #  if this is the case, more care is required.  
  # if it is on the lower boundary, broaden the range and repeat. (not implemented automatically)
  # if it is on the upper boundary, plot the gradient, plot muHatGLS vector above.  Investigate.  Use other estimators?
  
  return(c(results, F))
}



##################
###### SBM  ######
##################


# Qest = function(Tr, memb){
#   eg = ends(Tr, 1:length(E(Tr)))
#   Q = table(as.data.frame(cbind(memb[eg[,1]], memb[eg[,2]])))
#   if(nrow(Q) != ncol(Q)){
#     print("ERROR:  Q is rectangular")
#     return(NA)
#   }
#   
#     if(!diagnostics) return(c(mean(Y, na.rm = T),NA))
#     if(diagnostics){
#       return(list(muhatGLS = mean(Y, na.rm = T), eigenSummaries = data.frame(eigenValues = rep(NA, length(unique(memb))-1 ), yfInnerprod = rep(NA, length(unique(memb))-1 ))))
#     }
#   }
#   Q = (Q + t(Q))/2
#   
#   
# }
# 



sbmGLS = function(x, outcome, blockVariable, diagnostics = F, Q = NULL, deg= T){
  # x is a tibble or data.frame with the following columns:
  #  id, recruiter.id, network.size (or network.size.variable)
  # outcome is the name of a column in x that contains the measurement for which we want to estimate the population mean. 
  # blockVariable is the column name of a feature (as.character) which pseudo-stratifies the sampling.  
  #    By "pseudo-stratifies", these variables should be highly auto-correlated in referrals.  
  #    Future code will help to select these features. 
  #    If it is not specified, then it is set to outcome.
  # diagnostics = T outputs additional information
  # Q (i.e. the matrix which counts of transitions between blocks) can be specified.  If Q = NULL, then the code below computes Q.
  # If deg = T and network.size.variable are specified, then this uses adjusted Volz-Heckathorn weights... 
  #   the adjustment to VH is that the normalizing constant which estimates E(1/deg(X)) is computed via GLS, rather than via a sample average (as in the original VH).
  
  # using 
  #    the tree, 
  #    Y (which could be re-weighted via VH), 
  #    and a vector of memberships (e.g. demographic variables such as gender or race),
  # this function returns the fGLS estimate. 
  # Without an apriori good choice for memb, see the diagnostic functions.  
  #   Alternatively, use Y (first ensure it is discrete, e.g. "high" vs "low")
  # diagnostics = T returns more information about the spectral properties used in RSE and other diagnostics.
  # Q can be supplied.  If not, it is computed with memb.  
  # if deg is specified, this function will compute Volz-Heckathorn weights and use the GLS correction for the normalizing constant
  #   standard weights are computed as:
  #   d =degree(g)[s]
  #   H = 1/mean(1/d)
  #   w = H/d
  #   Yvh = w*Y
  #   $correction=H.gls = mean(1/d)/sbmGLS(Tr,Y = 1/deg,memb)
  
  # construct the tree:
  Tr = graph.edgelist(cbind(as.character(x[,"recruiter.id"][[1]]),as.character(x[,"id"][[1]])))
  
  makeNewSeed = F
  if(length(which(degree(Tr, mode = "in") == 0))>1){
    tmpRow = x[1,]
    tmpRow[] = NA
    tmpRow[1] = names(which(degree(Tr, mode = "in") == 0))
    x = rbind(tmpRow, x) %>% as_tibble
    makeNewSeed = T
  }
  if(!makeNewSeed) Tr = delete_vertices(Tr,v = 1)
  
  Y = x[,outcome][[1]]
  memb = x[,blockVariable][[1]]
  
  # Set the deg variable. 
  # if there are no degrees specified, set them to all be one. 
  if(length(deg)==0) deg = rep(1, length(Y))
  if(length(deg)>1 && length(deg) < length(memb)){
    print("deg is incorrectly specified")
    return(NA)
  }
  if(length(deg)==1){  
    if(!deg){
      deg = rep(1, length(Y))
    }else{
      if("network.size" %in% names(x)) deg = x[,"network.size"][[1]]
      if("network.size.variable" %in% names(x)) deg = x[,"network.size.variable"][[1]]
    }
  } 
  
  # if there is no variation in the observed Y's, output that value.
  if(sd(Y, na.rm=T)==0){
    if(!diagnostics) return(c(mean(Y, na.rm = T), NA))
    if(diagnostics) return(list(
      muhatGLS = mean(Y), 
      eigenSummaries = data.frame(
        eigenValues = NA, 
        yfInnerprod = NA),
      rse = 1,
      correction = NA
    ))
  }
  
  # these are the referral edges.
  eg = ends(Tr, 1:length(E(Tr)))
  # this fixes the indexing, so that eg doesn't contain node id's, but rather row # in x.
  eg = cbind(match(as.numeric(eg[,1]), as.vector(V(Tr)$name)),match(as.numeric(eg[,2]), as.vector(V(Tr)$name)))
  
  
  # function can specify the Q matrix.  If not, it is computed.
  if(length(Q) ==0){
    transitions = cbind(memb[eg[,1]], memb[eg[,2]])  # edge list for referral tree
    Q = table(as.data.frame(rbind(transitions,transitions[,2:1])))  # the funny rbind is to make Q symmetric.
    if(nrow(Q)==1){  # if there is no variation in Q, then we cannot make any GLS adjustments.
      if(!diagnostics) return(c(mean(Y, na.rm = T), NA))
      if(diagnostics) return(list(
        muhatGLS = mean(Y), # if degrees were specified, this is VH estimator.
        eigenSummaries = data.frame(
          eigenValues = NA, 
          yfInnerprod = NA),
        rse = 1,
        correction = NA
      ))
    }
    #     if(nrow(Q) != ncol(Q)){
    #       if(!diagnostics) return(c(mean(Y, na.rm = T),NA))
    #       if(diagnostics){
    #         return(list(muhatGLS = mean(Y, na.rm = T), eigenSummaries = data.frame(eigenValues = rep(NA, length(unique(memb))-1 ), yfInnerprod = rep(NA, length(unique(memb))-1 ))))
    #       }
    #     }
    # Q = (Q + t(Q))/2
  }
  # Added "regularization" or laplace smoothing here:
  L = diag(1/sqrt(rowSums(Q)+1))%*%Q%*%diag(1/sqrt(colSums(Q)+1))  
  eiL = eigen(L)
  
  # if L is approximately rank 1, then no need for GLS... just leads to numerical problems.  #  TODO:  add a diagnostic check to see if GLS is necessary.
  if(abs(eiL$values[2]) < 10^(-5)){   
    if(!diagnostics) return(c(mean(Y, na.rm  =T),0))
    if(diagnostics) return(list(
      muhatGLS = mean(Y), 
      eigenSummaries = data.frame(
        eigenValues = eiL$values[-1], 
        yfInnerprod = sqrt(coefs[-1]/var(Y, na.rm=T))),
      rse = 1,
      correction = NA
    ))
    #     if(diagnostics){
    #       return(list(muhatGLS = mean(Y, na.rm=T), eigenSummaries = data.frame(eigenValues = rep(0, length(unique(memb))-1 ), yfInnerprod = rep(NA, length(unique(memb))-1 ))))
    #     } 
    #     
  }
  
  
  # this is the function that is described in the paper.
  m = sum(Q)
  Z = model.matrix(~as.factor(memb) -1)
  if(makeNewSeed) Z = rbind(rep(0, ncol(Z)),Z)
  fhat = sqrt(m) * diag(1/sqrt(rowSums(Q))) %*% eiL$vec
  # Y = Y[-1]  # remove the "false seed"
  Yna = Y; Yna[is.na(Yna)] = 0
  Ybar = t(Y[complete.cases(Y)])%*%Z[complete.cases(Y),]/sum(complete.cases(Y))
  coefs = (Ybar%*%fhat)^2
  ac = acHat(coefs[-1], eiL$values[-1], diam = diameter(Tr,directed = F)+1)
  # if(var(Y)>ac[1]) ac[1] = var(Y)
  ac[1] = var(Y, na.rm=T) + ac[1] 
  
  # deg = deg[-1]
  Di = 1/(deg); Di[is.na(Di)] = 0 
  DiBar = t(Di[complete.cases(Di)])%*%Z[complete.cases(Di),]/sum(complete.cases(Di))
  DiCoefs = (DiBar%*%fhat)^2
  acD = acHat(DiCoefs[-1], eiL$values[-1], diam = diameter(Tr,directed = F)+1)
  # if(var(Y)>ac[1]) ac[1] = var(Y)
  acD[1] = var(Di, na.rm=T) + acD[1] 
  
  correction = mean(Di)/gls(Tr, acD, Di)
  
  if(!diagnostics) return(gls(Tr, ac, Y)*correction) 
  if(diagnostics){
    return(list(
      muhatGLS = gls(Tr, ac, Y) * correction, 
      eigenSummaries = data.frame(
        eigenValues = eiL$values[-1], 
        yfInnerprod = coefs[-1]),
      rse = rse(Tr=Tr, Y, lams = eiL$values[-1], beta2s = coefs[-1]),
      correction = correction
    ))}
}





rse = function(Tr, Yvh, lams, beta2s){
  # compute the ratio of standard errors (rse) for GLS / sampleAverage.  
  ac = acHat(beta2s, lams, diam = diameter(Tr,directed = F)+1)
  if(length(beta2s)>1) if(var(Yvh) > 5* ac[1] ) ac[1] = 2* ac[1]
  Sig = makeCov(Tr, ac)
  return(sqrt(sum(solve(Sig,rep(1, nrow(Sig))))^(-1) / sum(Gz(Tr, z = lams)[,2]*beta2s)))
}



#### Diagnostic plotting function #####










summary.gls4rds = function(x, outcome, blockVariable, legendTF = T, leftPlot = T, plotTitle = " ", bottomPlot = T, textSize=1){
  # summary.gls4rds = function(Yvh, Tr, memb = c(), legendTF = T, leftPlot = T, plotTitle = " ", bottomPlot = T, textSize=1){
  
  
  #   we want code that will plot the ratio of log_2 SE(VH) / SE(GLS).  
  #         that is log(sqrt(G(z)*n)/sqrt((1+z)/(1-z)),2)
  #   over the range [-.95,.95], unless we have an estimate of \lambda_2 outside of this...
  #   In addition to that line, we want marks on the x-axis corresponding to 
  #   \lambda_auto, \lambda_\Delta, \lambda_{SBM-Y}
  #   Finally, if there is a partition z, for each \lambda_j, j>1 of Q_L, 
  #   make a stick proportional to <y,f_j>^2 / var(y) at \lambda_j.  
  #        The sum of the stick heights should be the maximal value of the line.  
  #        Also, add a line at zero to reflect how much variation in y is not explained by the spectrum of Q_L.
  
  
  ####  1)  compute all eigenvalues and <y,f_j>^2 for SBM-z
  ####  2)  compute G(z) and log(sqrt(G(z)*n)/sqrt((1+z)/(1-z)),2) on the appropriate range
  ####  3)  plot the function, plot the marks, plot the sticks. 
  
  
  Tr = graph.edgelist(cbind(as.character(x[,"recruiter.id"][[1]]),as.character(x[,"id"][[1]])))
  
  makeNewSeed = F
  if(length(which(degree(Tr, mode = "in") == 0))>1){
    tmpRow = x[1,]
    tmpRow[] = NA
    tmpRow[1] = names(which(degree(Tr, mode = "in") == 0))
    x = rbind(tmpRow, x) %>% as_tibble
    makeNewSeed = T
  }
  Y = x[,outcome][[1]]
  memb = x[,blockVariable][[1]]
  
  if(!makeNewSeed) Tr = delete_vertices(Tr,v = 1)
  
  
  ####  1)  compute all eigenvalues and <y,f_j>^2 for SBM-z
  n = length(Y)
  muAUTO = autoGLS(Tr,Y)
  muDIFF = diffGLS(Tr,Y)
  muSBMy = sbmGLS(x, outcome, outcome,diagnostics = T, deg= F)
  if(length(memb) > 0) muSBMfull = sbmGLS(x, outcome, blockVariable,diagnostics = T, deg= T)
  
  eigVals = c(muAUTO[2], muDIFF[2], muSBMy$eigenSummaries$eigenValues)
  if(length(memb) > 0) eigVals = c(muAUTO[2], muDIFF[2], muSBMy$eigenSummaries$eigenValues,muSBMfull$eigenSummaries[,1])
  
  rng = range(eigVals)
  if(rng[1] > 0) rng[1] =  0
  if(rng[1] < 0) rng[1] =  - sqrt(abs(rng[1]))
  if(rng[2] < .95) rng[2] = .95
  if(rng[2] > .95) rng[2] = .999
  #   
  rng =  c(0,.97)
  ####  2)  compute G(z) and log(sqrt(G(z)*n)/sqrt((1+z)/(1-z)),2) on the appropriate range
  
  gz = Gz(Tr,z = c(seq(rng[1], rng[2], len= 100), eigVals))
  zseq = gz[,1]
  gz = gz[,2]
  Ratio = (sqrt(gz*n)/sqrt((1+zseq)/(1-zseq)))^(-1)
  
  
  
  
  ####  3)  plot the function, plot the marks, plot the sticks. 
  par(las = 1)
  ylab = " "
  if(leftPlot){
    ylab = "ratio of SE's"
  }
  xlab = " " 
  if(bottomPlot) xlab = "estimated eigenvalue"
  
  # plot(zseq[1:100], Ratio[1:100], type = 'l', col = "grey", xlab = xlab, ylab = ylab, log = "y", main  = plotTitle)
  if(leftPlot) plot(zseq[1:100], Ratio[1:100], ylim = c(min(Ratio[1:100])-.05,1 + .05),type = 'l', col = "grey", xlab = xlab, ylab = ylab, main  = plotTitle)
  if(!leftPlot) plot(zseq[1:100], Ratio[1:100], ylim = c(min(Ratio[1:100])-.05,1 + .05),type = 'l', col = "grey", xlab = xlab, ylab = ylab, main  = plotTitle, yaxt ="n")
  points(eigVals[1],rse(Tr, Y, eigVals[1], var(Y)), pch = "a", cex = textSize)
  points(eigVals[2],rse(Tr, Y, eigVals[2], var(Y)), pch = 2, cex = textSize)
  points(eigVals[3],rse(Tr, Y, eigVals[3], var(Y)), pch = "y", cex = textSize)
  
  if(length(memb) > 0){
    nostics = muSBMfull$eigenSummaries
    yaxis = muSBMfull$rse
    yaxis = rse(Tr,Y, nostics[,1], nostics[,2])
    for(k in 1:nrow(nostics)){
      points(nostics[k,1], yaxis, pch ="z", cex = textSize)
    }
  }
  lines(c(-1,1), c(1,1), col = "grey", lty=3)
  # legend("topleft", pch = list(97,2,122,121), legend = c(expression(hat(mu)[auto]), expression(hat(mu)[Delta]), expression(hat(mu)[SBM-y]), expression(hat(mu)[SBM-z])))
  if(legendTF){
    # if(length(memb) > 0) 
    legend("bottomleft", pch = list(97,2,121,122), legend = c("auto", expression(Delta), "SBM-y", "SBM-z"))
    # if(length(memb) == 0) legend("bottomleft", pch = list(97,2,121), legend = c("auto", expression(Delta), "SBM-y"))
    #   
    #   if((length(memb) >0) & (mean(Y==0) >.05) & !quiet){
    #     #   Still testing this feature.
    #     tmp =fisher.test(Y, memb)
    #     names(tmp)
    #     
    #     print("For the null hypothesis that the labels Y are independent from the memberships,")
    #     if(tmp$p.value < .05){
    #       paste("Fisher's exact test rejects at level .05, with p-value = ", tmp$p.value,".", sep ="")
    #       print("This suggests that the supplied memberships are relevant and provides evidence")
    #       print("in support of using SBM-z")
    #     }
    #     if(tmp$p.value > .05){
    #       paste("Fisher's exact test does not reject at level .05, with p-value = ", round(tmp$p.value,2),".", sep ="")
    #       print("This suggests that the supplied memberships are not relevant for this outcome.")
    #       print("SBM-z is unlikely to diminish the variability of the standard estimator.")
    #   }
    #   
    
    print(paste("muHat_auto =", round(muAUTO[1],2)))
    print(paste("muHat_Delta =",round(muDIFF[1],2)))
    print(paste("muHat_SBM-y =",round(muSBMy$muhatGLS,2)))
    if(length(memb) > 0) print(paste("muHat_SBM-z =",round(muSBMfull$muhatGLS,2)))
    
  }
  
  
}







Gz = function(treei = c(), dsi = c(), z = seq(0,.95,len = 100), plott = T){
  # given a referral tree, compute the function G(z), where G is the probability generating function for the following random variable D:
  #  select two nodes I and J uniformly from the tree and define D to be the distance between I and J.  
  #       the smallest value is D=0, the largest value is the diameter of the tree. 
  #  We are most interested in this function for values between -1 and 1.  
  
  
  good = F
  if(length(dsi)>0){
    ni = (1+sqrt(1+8*sum(dsi)))/2  # solve quadratic formula to get number of nodes from sum of pairs: sum(dsi) 
    good = T
  }
  if(length(treei)>0){
    dsi = distance_table(treei,directed = F)$res
    ni = length(V(treei))
    good = T
  }
  if(good ==F){
    print("you need to pass the referral tree or a table giving the 'histogram' of pairwise distances in the tree.")
    return(NA)
  }
  gz = rep(NA, length(z))
  pathLengths=c(ni,2*dsi)
  for(i in 1:length(z)) gz[i] = z[i]^(0:length(dsi))%*% pathLengths/ni^2
  
  dat = cbind(z, as.vector(gz))
  #   which(diff(dat[,2])/diff(dat[,1]) > 1)
  #   
  #   colnames(dat) = c("z", "G(z)")
  #   if(plott) plot(dat, type = "l", las = 1, log = "y")
  #   plot(dat[-1,1], diff(log(dat[,2]))/diff(dat[,1]))
  #   plot(dat[-(1:2),1], diff(diff(log(dat[,2])))/diff(dat[-1,1]))
  #   
  #   
  #   plot(dat[,2])
  #   
  return(dat)
}












