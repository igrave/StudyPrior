#' Calculate mean squared error based on the mode of the posterior
#'
#' @param prior Prior to calculate posterior. Specify posterior instead if available.
#' @param prob.range Range of values to calculate MSE over
#' @param length Number of values to calculate MSE for
#' @param n.binom Number of patients in new trial
#' @param mc.cores Number of cores for parallel
#' @param posterior Posterior density 
#'
#' @return A vector of error values
#'
#' @examples \donttest{
#' xh <- c(30,40,50)
#' nh <- c(90,95,110)
#' fix <- binom.PP.FIX(x=xh, n=nh, d=c(.2,.5,.7), mix=TRUE) # set different weights
#' post <- posterior.mixture.prior(34,75, mixture.prior = fix)  
#' calc.MSE.mode(posterior = fix, prob.range=c(0,1), n.binom=100)
#' }
#' 
calc.MSE.mode <- function(prior, prob.range=c(.5,1), length=20, n.binom=30, mc.cores=1, posterior){

  P <- seq(prob.range[1],prob.range[2],len=length)



  if(missing(posterior)){
  MSE.for.x <- parallel::mclapply(0:n.binom, function(Xs){


      if(inherits(prior, "function")){
        post <- function(p,g=1) prior(p,Xs)*dbinom(x=Xs, size=n.binom, prob=p)/g
        #mode
         p <- optimise(post, interval=c(0,1), maximum=TRUE)
        #   square error
        sq.err <- (p$maximum-P)^2
        return(sq.err)
      } else if(inherits(prior, "mixture.prior")){
        post.list <- posterior.mixture.prior(Xs, n.binom,  mixture.prior=prior)
        return((mean.mixture.prior(post.list)-P)^2 + var.mixture.prior(post.list))
      } else if(inherits(prior, "list")){
        post.list <- posterior.mixture.prior(Xs, n.binom,  mixture.prior=prior[[Xs+1]])
        return((mean.mixture.prior(post.list)-P)^2 + var.mixture.prior(post.list))
      }
    }, mc.cores = mc.cores)
  } else if(!missing(posterior)){
    MSE.for.x <- parallel::mclapply(posterior, function(post){
      p <- optimise(post, interval=c(0,1), maximum=TRUE)
      sq.err <- (p$maximum-P)^2
      return(sq.err)
      }, mc.cores = mc.cores)
  }


  MSE.for.x <- matrix(unlist(MSE.for.x), nrow=length)


  sapply(seq_along(P), function(i){
    sum(MSE.for.x[i,] * dbinom(0:n.binom, size=n.binom, prob=P[i]))})
}

