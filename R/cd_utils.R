# causal inference method
# 
# implementation details:
# if 'X -> Y'  output cd  = 1
# if 'Y -> X'  output cd = 0
# epsilon: confidence (or score)
library(rvinecopulib)
library(statmod)
## Proper scoring rule for predicted quantiles
quantileScoring <- function(actual, pred, prob = 0.95) {
  mean((as.numeric(actual <= pred) - prob) * (pred - actual))
}


QCCD <- function(pair){
  # changing the seed might slightly affect the results due to the random jittering
  # in the rank function
  set.seed(0)
  n <- length(pair[,1])
  X = pair[,1]
  Y = pair[,2]
  x = X
  y = Y
  
  # Recover pseudo-observations and estimate the copula non-parametrically
  u1 <- rank(x, ties.method = "random")/(n + 1)
  u2 <- rank(y, ties.method = "random")/(n + 1)
  a <- acepack::ace(x, y)
  cop <- bicop(data = cbind(u1,u2),
               family_set = "tll",
               nonpar_method = "constant",
               mult = n^(1/6 - 1/5) * abs(cor(x, y)) / abs(cor(a$tx, a$ty)))
  
  uw <- gauss.quad.prob(3)
  h <- sapply(uw$nodes, function(uu) {
    
    u1p <- predict(object = cop, newdata = cbind(uu, u2), what = "hinv2")
    u2p <- predict(object = cop, newdata = cbind(u1, uu), what = "hinv1")
    xp <- quantile(x, u1p)
    yp <- quantile(y, u2p)
    
    h1 <- quantileScoring(x, xp, uu)
    h2 <- quantileScoring(y, yp, uu)
    
    rel_sc = h1/(h1 + h2)
    rel_sc
  })
  r1 <- sum(uw$weights[!is.na(h)]*h[!is.na(h)])/sum(uw$weights[!is.na(h)])   
  
  cd = ifelse(r1 > 0.5, 1, 0)
  return(list(cd = cd, epsilon = r1))
  
}


QCCDwrap <- function(X, Y){
  res = QCCD(cbind(X,Y))
  if(!is.na(res$cd)) {
    cd = ifelse(res$cd == 1, "->", "<-")
  } else{
    cd = "--"
  }
  list(cd = cd, eps = res$eps)
}


# code for reading in data pairs is borrowed from 
# Marx, A. and Vreeken, J. Telling Cause from Effect using MDL-based Local and Global Regression.
# In ICDM, 2017
uv = c(1:51,56:70,72:104,106)
ref = read.table("../data/tuebingen_benchmark/README_polished.tab", sep="\t", header=F, stringsAsFactors = F)
meta = read.table("../data/tuebingen_benchmark/pairmeta.txt")
ref.uv = ref[uv, ]

readI = function(i, r=ref){
  f = paste(c("../data/tuebingen_benchmark/", r$V1[i], ".txt"), collapse="")
  t = read.table(f, sep=" ", header=F, stringsAsFactors = F)
  return(t)
}

