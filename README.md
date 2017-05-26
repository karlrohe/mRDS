MarkovRDS: A collection of techniques for Respondent-Driven Sampling (RDS) based upon the Markovian model
==============================

There are currently two key procedures in this package, `atb` and `atb_power`.  `atb` performs the assisted treebootstrap procedure described in the paper "A critical threshold for design effects in network sampling". `atb_power` is similar, but it should be used for sample size and power calculations; it is parameterized by a priori guesses of the key quantities in atb.  

This package is underdevelopment. More functions, for different sorts of operations, will be added over time.  While moderate efforts will be made to make the revised packages "backwards compatible," there are absolutely no guarantees.  

Installation in R
------------

```R
source("https://raw.githubusercontent.com/karlrohe/mRDS/master/MarkovRDS.R")
```



Functions 
------------
```R
atb(x, outcome = "HIV", blockVariable = NULL, Bcount, rawSamples = F, verbose = T, pretty = T)
atb_power(n,m, A, seeds = .9, pi = NULL, ybar = NULL, Bcount= 1000, for.simulation =F, rawOutput=F)
```


Arguments 
------------
```R
x                 # same format as given by the package RDS. 
                  #   (requires columns id, recruiter.id, network.size, outcome, blockVariable)
outcome           # name of a variable in x.  Population mean of this quantity is the target.
blockVariable     # name of a variable in x.  Used to construct the groups. If not specified, 
                  #   then set to outcome.  
Bcount            # number of bootstrap samples

rawSamples        # if TRUE, returns an (n x Bcount) matrix identifying each of the bootstrap samples.  
                  #    The elements of the returned matrix are row indices corresponding to rows of x.  
verbose           # verbose = F should be used in simulation experiments. 
pretty            # pretty = F returns the Bcount-many bootstrapped point estimates.


n                 # sample size
m                 # average number of coupons returned by participants 
                  #   (will be simulated as Poisson distributed)
A                 # this is a K x K matrix where the row normalized version creates a 
                  #   Markov transition matrix on 1 ... K.
                  #   Alternatively, this can be a single number which is the probability of a 
                  #   referral within the same class, with this parameterization, 
                  #   all other probabilities are uniform. 
seeds             # if 1 or more, this will be the number of seeds used to initialize the process and
                  #    the results are conditional on the process reaching the target population size n.
                  # if seeds < 1, then it is a probability.  the number of seeds will be choosen so that 
                  #    the probability of reaching n samples is at least that probability.
pi                # this is a K vector. Element j is the proportion of individuals in group j.
                  #   if not specified, it is set to uniform on 1 ... K.
ybar              # set to a K vector. 
                  #   Element j is the proportion of individuals in block j that are HIV+. 
                  #   if not specified, it is set to c(0,1).
Bcount            # this is the number of simulations that are performed. 
for.simulation    # if True, then there will be less text printed to the terminal window and 
                  #   the function returns the width of the confidence interval
rawOutput         # if True, then there will be less text printed to the terminal window and 
                  #   the function returns Bcount many simulations of the Volz-Heckathorn estimator
```

Details
------------
`atb` performs the resampling procedure described in ``A critical threshold for design effects in network sampling''.  This requires all of the data that is typically collected as part of an RDS and a specification of which variable will serve as the group labels. 

`atb_power` it is acceptable to only specify `n,m' and `A` as a single probability in (0,1).  This is equivalent to setting `pi = c(1/2, 1/2), ybar = c(0,1)`.  This function runs Bcount simulations from the model parameterized by the arguments.  For each of those simulations, it compute the Volz-Heckathorn estimator.  Then, it computes the 90% interval between the 5th and 95th percentiles of those simulated values.  Increasing n, decreasing m, or making A closer to .5 will all make the confidence interval more narrow.  In practice, we have little control over A (for an exception, see https://arxiv.org/abs/1606.00387).

Values
------------
`atb` returns a 90% confidence interval for all specified outcomes.

`atb` returns an estimate of the 90% confidence interval that would be obtained in an experiment with the specified settings.  

All confidence intervals are constructed as the 5th and 95th percentiles of the bootstrap distribution.  While 95% intervals are common in many fields, due to the massive model uncertainty that we have for RDS (most importantly, referrals are not random, but are instead modeled as random), this package will only report 90% intervals.

For `atb_power` see the example usage below.  I

Example Usage
-------------

In order to perform a bootstrap on an example, first simulate a social network with fastRG, then simulate the RDS process on the social network.  With the resulting data, compute a 90% confidence interval for the population prevalence of HIV.  



```R
source("https://raw.githubusercontent.com/karlrohe/mRDS/master/simulationFunctions.R")
# this imports functions (1) getGraph and (2) simpleRDS that 
# (1) simulate a graph and (2) an RDS on that graph.

# simulate a 100k node graph from the degree corrected stochastic block model:
G = getGraph()
summary(G)
# the target value is:
V(G)$HIV %>% mean
# the blocks are 
V(G)$block %>% table


# To estimate the target value, take an RDS of G, where the number of referrals per participant is 
#   Poisson(\lambda=4). 
# Simulation is done without replacement.
# Sample size is 1000
dat = simpleRDS(G,4, 1000)
str(dat)

# Under the default settings, atb does not return point estimates.
#   Instead, it prints the summaries to the screen.
atb(dat, blockVariable = "block", outcome = "HIV", Bcount = 1000)

# for simulation...
# this returns VH estimator for each of B runs:
bootEstimates = atb(dat, blockVariable = "block", outcome = "HIV", Bcount = 1000, pretty = F, verbose = F)

# this returns the indices of dat that are sampled 
bootSamples = atb(dat, blockVariable = "block", outcome = "HIV", Bcount = 1001, pretty = F, verbose = F, rawSamples = T)
dim(bootSamples)
# for example, here is the first "bootstrapped" data set:
bootData = dat[bootSamples[,1],] 
# note that "recruiter.id" is the recruiter from the original sample, not the bootstrapped sample. 

```


Or, perform a power calculation
```R
# If you set seeds as a probability...

atb_power(n=1000, m=2, A=.8, seeds=.9, Bcount= 1000)

# If you set seeds as an integer... 
#   atb_power estimates the probability of chain death 
#   (based upon the same Galton-Watson model).

atb_power(n=1000, m=2, A=.8, seeds=10, Bcount=1000)

# If the settings of m and A exceed the critical threshold, then there is a warning displayed:

atb_power(n=1000, m=3, A=.8, seeds=10, Bcount=1000)

# if for.simulation=T, then the function is quiet and returns the width
width = atb_power(n=1000, m=3, A=.8, seeds=10, Bcount=1000, for.simulation=T)
width

# you can provide more detailed specifications of A, the proportion of individuals in each strata (pi), and the mean outcome within each strata (ybar)
K = 5
A = matrix(runif(5^2), nrow = 5); diag(A) = 3*rowSums(A)
pi = runif(5)
ybar = rbinom(5,1,.4); while(sum(ybar) ==0) ybar = rbinom(5,1,.4)
ybar 
atb_power(n = 1000, m=2, A=A, seeds=.9, pi=pi, ybar=ybar, Bcount=1000)

```


Upcoming papers/procedures/software:
1)  selecting blockVariables
2)  new estimators for population mean that are less sensitive to critical threshold. 
3)  Diagnostic plots
