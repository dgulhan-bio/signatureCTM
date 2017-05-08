LL <- function(x, mu, sigma){
  R = dnorm(x, mean = mu, sd = sigma)
  return(-sum(log(R)))
}

calculate_bic <- function(A, x, y){
  library('stats4')
  
  ndim = length(x) 
  ndata = dim(A)[[1]]

  diff = A %*% x - y
  
  sig = sd(diff)
  loglikelihood = sum(log(dnorm(diff, mean = 0, sd = sig)))
  
  coefs <- coef(mle(minuslogl = LL, start = list( sigma = sig, mu = 0), fixed = list( x = diff)))
  maxll <- sum(log(dnorm(diff, mean = coefs["mu"], sd= coefs["sigma"])))
 # print(sprintf('loglik = %.3f  max loglik = %.3f',loglikelihood, maxll))
 # loglikelihood = loglikelihood/maxll
  
  bic = -2*loglikelihood + ndim*(log(ndata) - log(2*pi))
  return(bic)

}
