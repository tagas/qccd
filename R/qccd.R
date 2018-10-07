# causal inference method
# 
# implementation details:
# if 'X -> Y'  output cd = 1
# if 'Y -> X'  output cd = 0
# epsilon: confidence (or score)
# 

library(rvinecopulib)
library(statmod)


## Proper scoring rule for predicted quantiles
quantileScoring <- function(actual, pred, prob = 0.95) {
  mean((as.numeric(actual <= pred) - prob) * (pred - actual))
}

## Causal discovery method
QCCD <- function(pair) {
  
  # to get reproducible jittering results for discreet data
  set.seed(0)
  
  n <- nrow(pair)
  # Recover pseudo-observations and estimate the copula non-parametrically
  u <- apply(pair, 2, function(x) rank(x, ties.method = "random")/(n + 1))
  cop <- bicop(data = u,
               family_set = "tll",
               nonpar_method = "constant")
  
  # deal with discrete data
  pair_scaled <- qnorm(u)

    # integrate over quantiles
    if(n < 200){
         uw <- gauss.quad.prob(1)
    } else {
      uw <- gauss.quad.prob(3)
    }

  cls <- sapply(uw$nodes, function(uu) {
    u_pred <- cbind(predict(object = cop,
                            newdata = cbind(uu, u[, 2]),
                            what = "hinv2"),
                    predict(object = cop,
                            newdata = cbind(u[, 1], uu),
                            what = "hinv1"))
    
    # marginal and conditional quantiles
    marg_q <- sapply(1:2, function(i) quantile(pair_scaled[,i], uu))
    cond_q <- sapply(1:2, function(i) quantile(pair_scaled[,i], u_pred[, i]))
    
    # code lengths
    cl_marginal <- sapply(1:2, function(i)
      quantileScoring(pair_scaled[, i], marg_q[i], uu))
    cl_conditional <- sapply(1:2, function(i)
      quantileScoring(pair_scaled[, i], cond_q[, i], uu))
    
    c(cl_marginal, cl_conditional)
  })
 
  sel <- !apply(is.na(cls), 2, any)
  uw$weights <- uw$weights[sel] / sum(uw$weights[sel])
  cls <- apply(na.omit(t(cls)) * uw$weights, 2, sum)
  
  dx_to_y <- (cls[1] + cls[4])/sum(cls[1:2])
  dy_to_x <- (cls[2] + cls[3])/sum(cls[1:2])
  
  cd <- ifelse(dy_to_x > dx_to_y, 1, 0)
  epsilon <-  (-dx_to_y + dy_to_x )
  
  return(list(cd = cd, epsilon = epsilon))
}


# wrapper used for the real data pairs
QCCDwrap <- function(X, Y){
  res = QCCD(cbind(X,Y))
  if(!is.na(res$cd)) {
    cd = ifelse(res$cd == 1, "->", "<-")
  } else{
    cd = "--"
  }
  list(cd = cd, eps = res$eps)
}
