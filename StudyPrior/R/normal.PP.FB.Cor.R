#' Full Bayes Power Prior with Correlated Weight Parameters using Rao-Blackwell
#'
#' @param x vector of means of historical studies
#' @param verbose Print messages
#' @param mixture.size Number of components in the mixture. If d.prior.cor=0, will be rounded to nearest value a^length(x) to make the grid evenly spaced.
#' @param mc.cores Number of cores for parallel
#' @param d.prior.cor Correlation parameter for transformed multivariate normal prior on weights
#' @param sd vector of standard deviations of historical studies
#'
#' @return A density function
#'
#'
normal.PP.FB.COR <- function(x, sd,  verbose=FALSE, mixture.size=1000, d.prior.cor=0, mc.cores=1){
  n.hist <- length(x)
  
  sigma <- matrix(d.prior.cor, ncol=n.hist, nrow = n.hist)
  diag(sigma) <- rep(1, n.hist)
  
  D <- t(pnorm(rmvnorm(mixture.size, mean=rep(0,n.hist), sigma=sigma) ))
  
  v <- sd^2/D  #for convenience
  all.v <- apply(v, 2, function(V) 1/sum(1/V))
  
  all.mean <-   colSums(x/v) * all.v
  
  pars <- matrix(c(all.mean, sqrt(all.v)), ncol=2)
  
  mix <- create.mixture.prior("normal", pars, weights=rep(1/mixture.size,mixture.size))

  f <- function(p,X) eval.mixture.prior(p, mix)
  
  return(f)
}
#   
#   
#   if(d.prior.cor==1){#pooled case
#     sumx <- sum(x)
#     sumnx <- sum(n)-sumx
#     D <- seq(0,1,len=mixture.size)
#     
#     
#     pars <- outer(D, c(sumx, sumnx))+1
#     
#     mix <- create.mixture.prior("beta",
#                                 pars,
#                                 weights=rep(1/mixture.size,mixture.size))
#     
#     
#   } else if(d.prior.cor==0){#independent case
#     leng <- round(mixture.size^(1/n.hist))
#     mixture.size <- leng^n.hist
#     
#     D <- as.matrix(expand.grid(rep(list(seq(0,1,len=leng)),n.hist)))
#     pars <- D %*% matrix(c(x, n-x), ncol=2) +1
#       
#     mix <- create.mixture.prior("beta", pars, weights=rep(1/mixture.size,mixture.size))
#     
#     
#   } else{ #anything else
#     sigma <- matrix(d.prior.cor, ncol=n.hist, nrow = n.hist)
#     diag(sigma) <- rep(1, n.hist)
#     
#     D <- pnorm(rmvnorm(mixture.size, mean=rep(0,n.hist), sigma=sigma) )
#     
#     pars <- D %*% matrix(c(x, n-x), ncol=2) +1
#     
#     mix <- create.mixture.prior("beta", pars, weights=rep(1/mixture.size,mixture.size))
#   }
#   
#   
#   f <- function(p,X) eval.mixture.prior(p, mix)
#   
#   return(f)
#   
# }
