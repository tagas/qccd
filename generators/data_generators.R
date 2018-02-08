library(MASS)

# Parts of these generatos are inspired by the supplementary material for 
# J. Peters, J. Mooij, D. Janzing, B. Sch\"olkopf: 
# "Causal Discovery with Continuous Additive Noise Models", JMLR 2014
# 
# 
# Functions are modified to adapt for our experiments, so we include their copyright notice.
# Copyright (c) 2013  Jonas Peters [peters@stat.math.ethz.ch]
# All rights reserved.
# 


sample_AN <- function(n, mu = 0, sd = 1){
  ran = rnorm(n, mean = mu, sd)
  
  # noiseExp > 1 not Gaussian nose
  noiseExp = 1
  
  noiseVar = runif(n, 1, 2)
  noiseVar_ch = runif(n, 1, 2)
  
  #1. generate X_pa
  X_pa = (sqrt(noiseVar) * abs(ran)) ^ (noiseExp) * sign(ran) 
  
  #2. generate a random function of X_pa using Gaussian Process RBF(Gaussian) Kernel with bandwidth=1
  kernPa <- computeGaussKernel(X_pa, 1, 1)
  fpa <- mvrnorm(1, rep(0, n), kernPa)
  
  #Generate random noise for the child
  ran <- rnorm(n, mean = mu, sd)
  noisetmp <- (0.2 * sqrt(noiseVar_ch) * abs(ran)) ^ (noiseExp) * sign(ran)
  X_child = fpa + noisetmp
  X = cbind(X_pa, X_child)
  X
}


sample_ANs = function(n){
  p = 2
  G = matrix(0, 2,2)
  G[1,2] = 1
  
  parsFuncType = list(kap=1,sigmax=1,sigmay=1,output=FALSE)
  noiseType = "normalRandomVariances"
  parsNoise = list(noiseExp=1,varMin=1,varMax=2)
  
  pair = sampleDataFromMonotoneSigmoid(n,G,parsFuncType, noiseType,parsNoise )  
  pair
}


sample_MN_g <- function(n, mu = 0, sd = 1){
  ran = rnorm(n, mean = mu, sd)
  
  # noiseExp>1 not Gaussian nose
  noiseExp = 1
  
  noiseVar = runif(n, 1, 2)
  noiseVar_ch = runif(n, 1, 2)
  
  #1. generate X_pa
  X_pa = (sqrt(noiseVar) * abs(ran)) ^ (noiseExp) * sign(ran) 
  
  #2. generate a random function of X_pa using Gaussian Process RBF(Gaussian) Kernel with bandwidth=1
  kernPa <- computeGaussKernel(X_pa, 1, 1)
  fpa <- mvrnorm(1, rep(0, n), kernPa)
  
  #Generate random noise for the child
  ran <- rnorm(n, mean = mu, sd)
  noisetmp <- ran
  X_child = fpa * noisetmp
  X = cbind(X_pa, X_child)
  X
}


sample_MN_u = function(n){
  p = 2
  G = matrix(0, 2,2)
  G[1,2] = 1
  
  parsFuncType = list(kap=1,sigmax=1,sigmay=1,output=FALSE)
  noiseType = "normalRandomVariances"
  parsNoise = list(noiseExp=1,varMin=1,varMax=2)
  
  pair = sampleDataFromMonotoneSigmoid_multiNoise(n,G,parsFuncType, noiseType,parsNoise )  
  pair
}


sample_HN <- function(n, mu = 0, sd = 1){
  ran = rnorm(n, mean = mu, sd)
  
  noiseExp = 1
  
  noiseVar = runif(n, 1, 2)
  noiseVar_ch = runif(n, 1, 2)
  
  X_pa = (sqrt(noiseVar) * abs(ran)) ^ (noiseExp) * sign(ran) 
  kernPa <- computeGaussKernel(X_pa, 1, 1)
  fpa <- mvrnorm(1, rep(0, n), kernPa)
  
  ran <- rnorm(n,0, 1)
  noisetmp <- (0.2 * sqrt(noiseVar_ch) * abs(ran)) ^ (noiseExp) * sign(ran)
  a = runif(1, 0.1,0.9)
  b = runif(1, 0.1,0.9)
  X_child = a*fpa + noisetmp*(1 + b*X_pa)
  X = cbind(X_pa, X_child)
  X
}


sample_HNs = function(n_size){
  p = 2
  G = matrix(0, 2,2)
  G[1,2] = 1
  
  parsFuncType = list(kap=1,sigmax=1,sigmay=1,output=FALSE)
  noiseType = "normalRandomVariances"
  parsNoise = list(noiseExp=1,varMin=1,varMax=2)
  
  pair = sampleDataFromMonotoneSigmoid_HN(n_size,G,parsFuncType, noiseType,parsNoise )  
  pair
}


computeGaussKernel <- function(x, sigmay, sigmax)
{
  if (is.matrix(x) == FALSE){
    x<-as.matrix(x) }
  
  n <- dim(x)[1]
  xnorm<-as.matrix(dist(x,method = "euclidean", diag=TRUE, upper=TRUE))
  xnorm<-xnorm^2
  
  KX <- sigmay * exp(-xnorm/(2*sigmax^2))
  
  return(KX)
}


computeCausOrder <- function(G)
 {
  p <- dim(G)[2]
  remaining <- 1:p
  causOrder <- rep(NA,p)
  for(i in 1:(p-1))
  {
    root <- min(which(colSums(G) == 0))
    causOrder[i] <- remaining[root]
    remaining <- remaining[-root]
    G <- G[-root,-root]
  }
  causOrder[p] <- remaining[1]    
  return(causOrder)
}


sampleDataFromMonotoneSigmoid <- function(n, G, parsFuncType, noiseType, parsNoise)
  # INPUTS:   n:  number of samples
  #           G:  adjacency matrix of Graph to simulate from
  #           
  #           
  # OUTPUTS:  X:      sampled data
  #
  # This function samples from modified (MONOTONE) sigmoid function 
  # c*b*(x+a)/(1+abs(b*(x+a))) where the choice of a,b,c is random.
  
{
  p <- dim(G)[2]
  X <- matrix(NA,n,p)
  # determine the causal Order which is needed for sampling
  causOrder <- computeCausOrder(G)
  
  if(parsFuncType$output)
  {
    show(causOrder)
  }
  
  # sample noise variances
  noiseVar <- runif(p,parsNoise$varMin,parsNoise$varMax)
  
  # loop through each node according to the causal order
  for(node in causOrder)
  {
    if(parsFuncType$output)
    {
      cat("generating GP for node ", node, "\r")
    }
    paOfNode <- which(G[,node] == 1)
    # simulation of noise at source nodes
    if(length(paOfNode) ==0)
    {
      if(noiseType == "normalRandomVariances" || noiseType == "normalRandomVariancesFixedExp")
      {
        ran <- rnorm(n)
        noisetmp <- (sqrt(noiseVar[node]) * abs(ran))^(parsNoise$noiseExp) * sign(ran)
      } else
      {
        error("This noiseType is not implemented yet.")
      }
      X[,node] <- noisetmp
    } else
    {
      nuPa <- length(paOfNode)
      X[,node] <- rep(0,n)
      
      # If kap>0 there is an additive model component
      if(parsFuncType$kap>0)
      {
        for(pa in paOfNode)
        {
          a.sig <- runif(n=1, min=-2, max=2)
          bern <- rbinom(1,1,0.5)
          b.sig <- bern*runif(n=1, min=0.5, max=2) + (1-bern)*runif(n=1, min=-2, max=-0.5)
          c.sig <- rexp(n=1,rate=4)+1
          X[,node] <- X[,node] + c.sig*(b.sig*(X[,pa]+a.sig))/(1+abs(b.sig*(X[,pa]+a.sig)))
          if(parsFuncType$output)
          {
            plot(X[,pa],c.sig*(b.sig*(X[,pa]+a.sig))/(1+abs(b.sig*(X[,pa]+a.sig))))
          }
        }
      }
      
      # Additive noise
      if(noiseType == "normalRandomVariances" || noiseType == "normalRandomVariancesFixedExp")
      {
        ran <- rnorm(n)
        noisetmp <- (0.2*sqrt(noiseVar[node]) * abs(ran))^(parsNoise$noiseExp) * sign(ran)
      } else
      {
        error("This noiseType is not implemented yet.")
      }
      X[,node] <- X[,node] + noisetmp       
    }
  }
  
  return(X)
}


sampleDataFromMonotoneSigmoid_multiNoise <- function(n, G, parsFuncType, noiseType, parsNoise)
  # INPUTS:   n:  number of samples
  #           G:  adjacency matrix of Graph to simulate from
  #           
  #           
  # OUTPUTS:  X:      sampled data
  #
  # This function samples from modified (MONOTONE) sigmoid function 
  # c*b*(x+a)/(1+abs(b*(x+a))) where the choice of a,b,c is random.
  # Copyright (c) 2013  Jonas Peters  [peters@stat.math.ethz.ch]
  # All rights reserved.  See the file COPYING for license terms. 
  
{
  p <- dim(G)[2]
  X <- matrix(NA,n,p)
  # determine the causal Order which is needed for sampling
  causOrder <- computeCausOrder(G)
  
  if(parsFuncType$output)
  {
    show(causOrder)
  }
  
  # sample noise variances
  noiseVar <- runif(p,parsNoise$varMin,parsNoise$varMax)
  
  # loop through each node according to the causal order
  for(node in causOrder)
  {
    if(parsFuncType$output)
    {
      cat("generating GP for node ", node, "\r")
    }
    paOfNode <- which(G[,node] == 1)
    # simulation of noise at source nodes
    if(length(paOfNode) ==0)
    {
      if(noiseType == "normalRandomVariances" || noiseType == "normalRandomVariancesFixedExp")
      {
        ran <- rnorm(n)
        noisetmp <- (sqrt(noiseVar[node]) * abs(ran))^(parsNoise$noiseExp) * sign(ran)
      } else
      {
        error("This noiseType is not implemented yet.")
      }
      X[,node] <- noisetmp
    } else
    {
      nuPa <- length(paOfNode)
      X[,node] <- rep(0,n)
      
      # If kap>0 there is an additive model component
      if(parsFuncType$kap>0)
      {
        for(pa in paOfNode)
        {
          a.sig <- runif(n=1, min=-2, max=2)
          bern <- rbinom(1,1,0.5)
          b.sig <- bern*runif(n=1, min=0.5, max=2) + (1-bern)*runif(n=1, min=-2, max=-0.5)
          c.sig <- rexp(n=1,rate=4)+1
          X[,node] <- X[,node] + c.sig*(b.sig*(X[,pa]+a.sig))/(1+abs(b.sig*(X[,pa]+a.sig)))
          if(parsFuncType$output)
          {
            plot(X[,pa],c.sig*(b.sig*(X[,pa]+a.sig))/(1+abs(b.sig*(X[,pa]+a.sig))))
          }
        }
      }
      
      # Additive noise
      if(noiseType == "normalRandomVariances" || noiseType == "normalRandomVariancesFixedExp")
      {
         noisetmp <- runif(n,-1,1) 
        
      } else
      {
        error("This noiseType is not implemented yet.")
      }
      X[,node] <- X[,node] * noisetmp       
    }
  }
  
  return(X)
}








sampleDataFromMonotoneSigmoid_HN <- function(n, G, parsFuncType, noiseType, parsNoise)
  # INPUTS:   n:  number of samples
  #           G:  adjacency matrix of Graph to simulate from
  #           
  #           
  # OUTPUTS:  X:      sampled data
  #
  # This function samples from modified (MONOTONE) sigmoid function 
  # c*b*(x+a)/(1+abs(b*(x+a))) where the choice of a,b,c is random.
  
{
  p <- dim(G)[2]
  X <- matrix(NA,n,p)
  # determine the causal Order which is needed for sampling
  causOrder <- computeCausOrder(G)
  
  if(parsFuncType$output)
  {
    show(causOrder)
  }
  
  # sample noise variances
  noiseVar <- runif(p,parsNoise$varMin,parsNoise$varMax)
  
  # loop through each node according to the causal order
  for(node in causOrder)
  {
    if(parsFuncType$output)
    {
      cat("generating GP for node ", node, "\r")
    }
    paOfNode <- which(G[,node] == 1)
    # simulation of noise at source nodes
    if(length(paOfNode) ==0)
    {
      if(noiseType == "normalRandomVariances" || noiseType == "normalRandomVariancesFixedExp")
      {
        ran <- rnorm(n)
        noisetmp <- (sqrt(noiseVar[node]) * abs(ran))^(parsNoise$noiseExp) * sign(ran)
      } else
      {
        error("This noiseType is not implemented yet.")
      }
      X[,node] <- noisetmp
    } else
    {
      nuPa <- length(paOfNode)
      X[,node] <- rep(0,n)
      
      # If kap>0 there is an additive model component
      if(parsFuncType$kap>0)
      {
        for(pa in paOfNode)
        {
          a.sig <- runif(n=1, min=-2, max=2)
          bern <- rbinom(1,1,0.5)
          b.sig <- bern*runif(n=1, min=0.5, max=2) + (1-bern)*runif(n=1, min=-2, max=-0.5)
          c.sig <- rexp(n=1,rate=4)+1
          X[,node] <- X[,node] + c.sig*(b.sig*(X[,pa]+a.sig))/(1+abs(b.sig*(X[,pa]+a.sig)))
          
          if(parsFuncType$output)
          {
            plot(X[,pa],c.sig*(b.sig*(X[,pa]+a.sig))/(1+abs(b.sig*(X[,pa]+a.sig))))
          }
        }
        noise_w = runif(1, 0.1,0.9)
        var_from_pa = (noise_w*X[,pa]+1)
      }
      
      # Additive noise
      if(noiseType == "normalRandomVariances" || noiseType == "normalRandomVariancesFixedExp")
      {
        ran <- rnorm(n)
        noisetmp <- (0.2 * sqrt(noiseVar[node]) * abs(ran)) ^ (parsNoise$noiseExp) * sign(ran)
        
      } else
      {
        error("This noiseType is not implemented yet.")
      }
      
      a = runif(1, 0.1,0.9)
      X[,node] <- a*X[,node] + noisetmp*(var_from_pa)      
    }
  }
  
  return(X)
}



