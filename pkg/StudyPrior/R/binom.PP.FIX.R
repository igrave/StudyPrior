#' Fixed Weight Power Prior for Binomial Data
#'
#' @param x historical events
#' @param n historical trials
#' @param d fixed weight
#' @param p.prior.a shape1 parameter for beta prior on p
#' @param p.prior.b shape2 parameter for beta prior on p
#' @param mix Create a mixture
#'
#' @return Probability density function for paramater p
#'
#' @examples \donttest{
#' xh <- c(30,40,50)
#' nh <- c(90,95,110)
#' fix <- binom.PP.FIX(x=xh, n=nh, d=c(.2,.5,.7), mix=TRUE) # set different weights
#' fix.5 <- binom.PP.FIX(x=xh, n=nh, d=0.5, mix=TRUE) #use 0.5 for all studies
#' }
binom.PP.FIX <- function(x, n, d,  p.prior.a=1, p.prior.b=1, mix=FALSE){

  if(length(n) != length(x)) stop("Length of x and n do not match") 
  if(length(d)==1) d <- rep(d, length(x))
    
  if(mix){
    
    pars <- matrix(c(p.prior.a + sum(x*d), p.prior.b + sum(d*(n-x))),ncol=2 )
    
    create.mixture.prior("beta", pars, weights=1)
    
  } else{
    function(p,...) dbeta(p,p.prior.a + sum(x*d), p.prior.b + sum(d*(n-x)))
  }
  
}
