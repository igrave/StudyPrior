#' Test if object is a mixture.prior
#'
#' @param x Object 
#'
#' @return TRUE if \code{x} is a \code{mixture.prior}, else FALSE
#'  
#'
is.mixture.prior <- function(x) {
  inherits(x, "mixture.prior")
}


#' Create mixture model
#'
#' @param type Type of mixture, "normal" or "beta"
#' @param pars Matrix of parameter values
#' @param weights Mixture weights
#'
#' @return A \code{mixture.prior} object
#'  
#'@examples \donttest{
#'pars <- matrix(c(2.4,1,3.7,1),ncol=2)
#'
#' #A mixture of 0.5 Be(2.4,3.7) + 0.5 Be(1,1) 
#'mymix <- create.mixture.prior(type="beta", pars = pars, weights=c(0.5,0.5))
#' }

create.mixture.prior <- function(type, pars, weights=rep(1/nrow(pars),nrow(pars))){
  

  mix <- list(fun = switch(type, "normal" = dnorm,"beta" = dbeta))

  attr(mix,"type") <- type
  attr(mix,"par.type") <- switch(type,
                                "normal" = c("mean","sd"),
                                "beta" = c("shape1","shape2") )

  
  attr(mix,"pars") <- pars
  attr(mix,"weights") <- weights/sum(weights)
  
  
  if(length( attr(mix,"weights")) == nrow( attr(mix,"pars"))){
    attr(mix,"degree") <- length(attr(mix,"weights"))
  } else(stop("Degree mismatch between weights and parameters"))
  
  class(mix) <- c("mixture.prior", type)
  return(mix)
}


#' Update parameters or weights in mixture model
#'
#' @param ... Unused
#' @param object mixture.prior to update
#' @param pars New paramers 
#' @param weights New weights
#'
#' @return A \code{mixture.prior} object with updated weights and parameters
#'  
#' @method update mixture.prior
update.mixture.prior <- function(object, ..., pars, weights ){
  
  
  if(!missing(pars)){
    attr(object,"pars") <- pars
  }
  if(!missing(weights)){
    attr(object,"weights") <- weights/sum(weights)
  }
  
  if(length( attr(object,"weights")) == nrow( attr(object,"pars"))){
    attr(object,"degree") <- length(attr(object,"weights"))
  } else(stop("Degree mismatch between weights and parameters"))
  
 
  return(object)
}

#' Add a vague componenet to mixture model to make robust
#'
#' @param object mixture.prior to update 
#' @param weight weight of vague component in updated model
#'
#' @return A \code{mixture.prior} object with updated weights and parameters
#'  
make.robust <- function(object, weight ){
  
  if(inherits(object, "beta")) {
    
    w <- attr(object,"weights")
    attr(object,"weights") <- c(w*(1-weight), weight)
    
    attr(object,"pars") <- rbind(attr(object, "pars") , c(1,1))
    
    if(length( attr(object,"weights")) == nrow( attr(object,"pars"))){
      attr(object,"degree") <- length(attr(object,"weights"))
    } else(stop("Degree mismatch between weights and parameters"))
    
    
    return(object)
  } else(stop("Only implemented for Beta mixture models."))
}



#' Get density of mixture
#'
#' @param x Vector of values to evaluate density at
#' @param mixture.prior Mixture prior
#' @param subset Vector of values specifying which mixture components should be included, eg c(1,2)
#'
#' @return Densities at of prior at x 
#'   
#' @examples \donttest{
#' xh <- c(30,40,50)
#' nh <- c(90,95,110)
#' 
#' #fit a full Bayes power prior with 500 samples
#' fb <- binom.PP.FB(x=xh, n=nh, mixture.size=500, mix=TRUE) 
#' eval.mixture.prior(seq(0,1,0.1), mixture.prior = fb) #Evalute the prior
#' }

# eval.mixture.prior <- function(x, mixture.prior, subset){
#   
#   if(missing(subset)){
#     sapply(x, function(X){
#       attr(mixture.prior,"weights") %*%
#         mixture.prior$fun(X,
#                           attr(mixture.prior,"pars")[,1],
#                           attr(mixture.prior,"pars")[,2])
#     })
#   } else {
#     sapply(x, function(X){
#       attr(mixture.prior,"weights")[subset] %*%
#         mixture.prior$fun(X,
#                           attr(mixture.prior,"pars")[subset,1],
#                           attr(mixture.prior,"pars")[subset,2])
#     })
#   }
#   
# }


eval.mixture.prior <- function(x, mixture.prior, subset){
  if(missing(subset)){
    
    X <- rep(x, each=attr(mixture.prior,"degree"))
    FUN <- mixture.prior$fun
    Z <- FUN(X,
        attr(mixture.prior,"pars")[,1],
        attr(mixture.prior,"pars")[,2])
    Z2 <- matrix(Z,
      nrow = attr(mixture.prior,"degree"), ncol=length(x) )
    return(as.numeric(attr(mixture.prior,"weights") %*% Z2))
    
    
    sapply(x, function(X){
      attr(mixture.prior,"weights") %*%
        mixture.prior$fun(X,
                          attr(mixture.prior,"pars")[,1],
                          attr(mixture.prior,"pars")[,2])
    })
  } else {
    sapply(x, function(X){
      attr(mixture.prior,"weights")[subset] %*%
        mixture.prior$fun(X,
                          attr(mixture.prior,"pars")[subset,1],
                          attr(mixture.prior,"pars")[subset,2])
    })
  }
  
}

#' Plot mixture model
#'
#' @param x a mixture prior object to plot
#' @param at Vector of values to evaluate density at
#' @param stack Plot the mixture components as stacked regions
#' @param add lines to an existing plot
#' @param ... Additional arguments to \code{\link{lines}}
#'
#'  
#'@examples \donttest{
#' xh <- c(30,40,50)
#' nh <- c(90,95,110)
#' 
#' #fit a full Bayes power prior with 500 samples
#' fb <- binom.PP.FB(x=xh, n=nh, mixture.size=500, mix=TRUE) 
#' 
#' #fit a full Bayes power prior with 500 samples and correlation 0.9
#' fb2 <- binom.PP.FB(x=xh, n=nh, mixture.size=500, mix=TRUE,d.prior.cor=0.9) 
#' plot.mixture.prior(fb2, xlab="p",ylab="Density")
#' plot.mixture.prior(fb, lines.only=TRUE, col=2)
#' }


plot.mixture.prior <- function(x,..., at=seq(0,1,0.005), stack=FALSE, add=FALSE){
  
  
  if(add) lines.only <- TRUE else lines.only<-FALSE
  
  mixture.prior <- x
  
  parameter <- at
  
  density <- eval.mixture.prior(parameter, mixture.prior)
  if(!lines.only) plot(parameter,density, type='n', ...)


  if(stack==TRUE){
    
    z <- order(attr(mixture.prior, "weights"),
               decreasing = TRUE)

    for(i in 1:(attr(mixture.prior,"degree")-1)){

      lines(parameter,eval.mixture.prior(parameter, mixture.prior,subset=z[1:i]), type='l', ...)

    }
  }
  lines(parameter, density, type='l', ...)
}

#' Calculate the equivalent sample size of mixture model
#' 
#' @param mixture.prior Mixture prior object
#'
#' @return Sample size
#'  
#'@examples \donttest{
#' xh <- c(30,40,50)
#' nh <- c(90,95,110)
#' fb <- binom.PP.FB(x=xh, n=nh, mixture.size=500, mix=TRUE)
#' fb2 <- binom.PP.FB(x=xh, n=nh, mixture.size=500, mix=TRUE,d.prior.cor=0.9) 
#' ess.mixture.prior(fb2)
#' ess.mixture.prior(fb)
#' }
ess.mixture.prior <- function(mixture.prior){
  if(inherits(mixture.prior,"beta")) sum(attr(mixture.prior,"weights")%*% attr(mixture.prior,"pars")) else NA
}


#' Calculate the mean of mixture model
#'
#' @param x Mixture prior object
#' @param ... Ignored
#'
#' @return Mean of prior
#'  
#' @examples \donttest{
#' xh <- c(30,40,50)
#' nh <- c(90,95,110)
#' 
#'  #fit a full Bayes power prior with 500 samples
#' fb <- binom.PP.FB(x=xh, n=nh, mixture.size=500, mix=TRUE)
#' mean.mixture.prior(fb) #calculate the mean
#' }

mean.mixture.prior <- function(x, ...){
  if(inherits(x,"beta")) {
    as.numeric((attr(x,"pars")[,1]/rowSums(attr(x,"pars"))) %*% attr(x,"weights"))
  } else if(inherits(x,"normal")) {
    as.numeric(attr(x,"pars")[,1] %*% attr(x,"weights"))
  }
}

#' Calculate the variance of mixture model
#'
#' @param mixture.prior Mixture prior object
#'
#' @return Variance of prior
#'  
#' @examples \donttest{
#' xh <- c(30,40,50)
#' nh <- c(90,95,110)
#' 
#'  #fit a full Bayes power prior with 500 samples
#' fb <- binom.PP.FB(x=xh, n=nh, mixture.size=500, mix=TRUE)
#' 
#'  #fit a full Bayes power prior with 500 samples with correlation
#' fb2 <- binom.PP.FB(x=xh, n=nh, mixture.size=500, mix=TRUE,d.prior.cor=0.9)
#' 
#' # calculate the variances of the distributions
#' var.mixture.prior(fb2)
#' var.mixture.prior(fb)
#' }
#'
var.mixture.prior <- function(mixture.prior){

  m <- apply(attr(mixture.prior,"pars"), 1, function(x) x[1]/sum(x))
  v <- apply(attr(mixture.prior,"pars"), 1, function(x) x[1]*x[2]/( sum(x)^2 * (sum(x)+1) ))
  w <- attr(mixture.prior,"weights")

  if(length(v)==1) return(v)

  mw <- m*w
  #by law of total variance for partitioned space
  svw <- sum(v*w)
  sm2w <- sum(m^2*(1-w)*w)
  
  # bigsum <- 2* sum( unlist(sapply(2:length(m), function(i) sapply(1:(i-1), function(j) mw[i]*mw[j]))))
  mat <- outer(mw,mw)
  diag(mat) <- 0
  bigsum <- sum(mat)
     svw + sm2w -  bigsum
     

}
  
  
  
  #' Calculate posterior mixture model
  #'
#' @param xs Number of successes in new trial
#' @param ns Number of patients in new trial
#' @param mixture.prior Mixture prior object
#' @param sd Standard deviation for normal model
#'
#' @return A \code{mixture.prior} object
#'   
#'@examples \donttest{
#' xh <- c(30,40,50)
#' nh <- c(90,95,110)
#' 
#' #fit a full Bayes power prior with 500 samples
#' fb <- binom.PP.FB(x=xh, n=nh, mixture.size=500, mix=TRUE) 
#' 
#' #calulate the posterior
#' post.fb <- posterior.mixture.prior(xs=51, ns=100, mixture.prior=fb)
#' 
#' # Compare with plot
#' plot(fb, xlab="p",ylab="Density", ylim=c(0,12))
#' #Use lines.only to add to existing plot
#' plot(post.fb, lines.only=TRUE, col=3)
#' }
#'
posterior.mixture.prior <- function(xs, ns, sd, mixture.prior){
  pars <- attr(mixture.prior,'pars')
  # new.pars <- pars 
  old.weights <- attr(mixture.prior,'weights')
  # new.weights <- numeric(length = length(old.weights))
  
  if(inherits(mixture.prior,"beta")){
    new.pars <- pars + matrix(c(xs,ns-xs),
                              nrow=attr(mixture.prior,"degree"),
                              ncol=2,
                              byrow=TRUE)
    
    if(length(xs)>1 | length(ns)>1) stop("posterior.mixture.prior is not vectorized for xs and ns.")
    
    lik <- dbetabinom.ab(x=xs, size=ns, shape1=pars[,1], shape2=pars[,2])
    ws <- old.weights*lik
    new.weights <- ws/sum(ws)
  
    } else if(inherits(mixture.prior,"normal")){
      
      all.sd <-  sqrt(1/( 1/pars[,2]^2 + 1/sd^2 ))
      
      all.mean <- (pars[,1]/pars[,2]^2  +  xs/sd^2 ) * all.sd^2
      
      new.pars <- matrix(c(all.mean,all.sd),
                         nrow=attr(mixture.prior,"degree"),
                         ncol=2,
                         byrow=FALSE)
      
      ws  <- mapply(dnorm, x = xs, mean = pars[,1], sd = sqrt(pars[,2]^2+sd^2)) * old.weights
      new.weights <- ws / sum(ws)
       
    }
      
  
  
  flp <- update.mixture.prior(mixture.prior,
                              pars=new.pars,
                              weights= new.weights)
  return(flp)
}


#' Print mixture model
#'
#' @param x Mixture prior object
#' @param ... Ignored
#'
#'  
#'
print.mixture.prior <- function(x, ...){
  
  w <- attr(x,"weights")
  p <- attr(x,"pars")
  par.type <- attr(x,"par.type")
  df <- as.data.frame(p)
  po <- rowSums(df)*w
  df <- cbind(df,w, po)
  colnames(df) <- c(par.type, "weights", "ESS")
  rownames(df) <- NULL
  cat("\n",attr(x,"type"),"mixture\n")
  print(df)
  cat("\nTotal sample size: ",ess.mixture.prior(x),"\n")
}


#' Sample from the mixture distribution
#'
#' @param n Number of samples
#' @param mixture.prior Mixture prior object
#'
#' @return A vector of n samples from \code{mixture.prior}
#'  
#'
sample.mixture.prior <- function(n, mixture.prior){
  w <- attr(mixture.prior, "weights")
  p <- attr(mixture.prior, "pars")
  degree <-  attr(mixture.prior, "degree")
  each.dist <- sample(1:degree, size=n, replace=TRUE, prob=w)
  
  if(inherits(mixture.prior,"beta")) {
    unlist(
      sapply(1:degree, function(d){
        rbeta(sum(each.dist==d), shape1=p[d,1], shape2=p[d,2])
      }))
  }
  
  if(inherits(mixture.prior,"normal")) {
    unlist(
      sapply(1:degree, function(d){
        rnorm(sum(each.dist==d), mean=p[d,1], sd=p[d,2])
      }))}
}
  

#' Probabilty function
#'
#' @param q Quantile
#' @param mixture.prior Mixture prior object
#'
#' @return Probability for the given quantile.
#'  
#'
p.mixture.prior <- function(q, mixture.prior){
  w <- attr(mixture.prior,"weights")
  L <- attr(mixture.prior,"degree")
  p <- attr(mixture.prior,"pars")
  
  PFUN <- if(inherits(mixture.prior,"beta")) pbeta else if (inherits(mixture.prior,"normal")) pnorm
  
  sapply(q, function(q){
    w%*% PFUN(q, p[,1],p[,2])
    })
}



#' Quantile function
#'
#' @param p Probability
#' @param mixture.prior Mixture prior object
#'
#' @return Quantile for the given probability
#'  
#'
q.mixture.prior <- function(p, mixture.prior){
  if(any(p > 1 | p < 0)) stop("p must be between 0 and 1")
  

  
  sapply(p, function(p){
    #special cases
    if(inherits(mixture.prior,"beta")){
      if(p == 0) return(0)
      if(p == 1) return(1)
    } else if (inherits(mixture.prior,"normal")){
      if(p == 0) return(-Inf)
      if(p == 1) return(Inf)
    } 
    
    # optimise( function(q){(p.mixture.prior(q,mixture.prior)-p)^2},
              # interval=c(0,1))$minimum
    uniroot(function(q) (p.mixture.prior(q,mixture.prior)-p),
            interval=c(0,1), extendInt = "yes")$root
  })
}




