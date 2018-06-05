#' Full Bayes Power Prior with Correlated Weight Parameters using Rao-Blackwell
#'
#' @param x number of historical successes
#' @param n number historical patients
#' @param verbose Print messages
#' @param mixture.size Number of components in the mixture.
#' @param mc.cores Number of cores for parallel
#' @param p.prior.a shape1 parameter for initial beta prior on probability
#' @param p.prior.b shape2 parameter for initial beta prior on probability
#' @param d.prior.cor Correlation parameter for transformed multivariate normal prior on weights
#' @param mix return the mixture object
#'
#' @return A density function
#' @examples \donttest{
#' xh <- c(30,40,50)
#' nh <- c(90,95,110)
#' fb <- binom.PP.FB.COR(x=xh, n=nh, mixture.size=500, mix=TRUE) #fit a full Bayes power prior with 500 samples
#' mean.mixture.prior(fb) #calculate the mean
#' }
#'
#'
binom.PP.FB.COR <- function(x, n, verbose=FALSE, mixture.size=1000, d.prior.cor=0, p.prior.a=1, p.prior.b=1, mc.cores=1, mix=FALSE){
  n.hist <- length(x)
  
  
  if(d.prior.cor==1){#pooled case
    sumx <- sum(x)
    sumnx <- sum(n)-sumx
    D <- seq(0,1,len=mixture.size)
    
    
    pars <- outer(D, c(sumx, sumnx))+rep(c(p.prior.a,p.prior.b),each=mixture.size)
    
    mixture <- create.mixture.prior("beta",
                                pars,
                                weights=rep(1/mixture.size,mixture.size))
    
    
  # } else if(d.prior.cor==0){#independent case
  #   leng <- round(mixture.size^(1/n.hist))
  #   mixture.size <- leng^n.hist
  #   
  #   D <- as.matrix(expand.grid(rep(list(seq(0,1,len=leng)),n.hist)))
  #   pars <- D %*% matrix(c(x, n-x), ncol=2) +1
  #     
  #   mix <- create.mixture.prior("beta", pars, weights=rep(1/mixture.size,mixture.size))
  #   
    
  } else{ #anything else
    sigma <- matrix(d.prior.cor, ncol=n.hist, nrow = n.hist)
    diag(sigma) <- rep(1, n.hist)
    
    D <- pnorm(rmvnorm(mixture.size, mean=rep(0,n.hist), sigma=sigma) )
    
    pars <- D %*% matrix(c(x, n-x), ncol=2) +rep(c(p.prior.a,p.prior.b),each=mixture.size)
    
    mixture <- create.mixture.prior("beta", pars, weights=rep(1/mixture.size,mixture.size))
  }
  
  if(mix) return(mixture) else return(function(p,X) eval.mixture.prior(p, mixture))
}
