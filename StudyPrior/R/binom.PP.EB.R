#' Empirical Bayes Power Prior for Binomial Data
#'
#' @param x number of historical successes
#' @param n number historical patients
#' @param X vector of new successes to computer prior for
#' @param N number of new patients
#' @param verbose TRUE
#' @param mc.cores number of cores for parallel
#' @param p.prior.a shape1 parameter of initial beta prior for successes
#' @param p.prior.b shape2 parameter of initial beta prior for successes 
#' @param max.dn Set a limit on the sum of the weights*n (FALSE or a numerical value)
#' @param mix if TRUE return a mixture.prior object otherwise the function
#'
#' @return Either a list of mixture.prior objects which can be evaluated with eval.mixture.prior or a density function of the probability parmater p given X.
#'
#' @examples \donttest{
#' xh <- c(30,40,50)
#' nh <- c(90,95,110)
#' eb <- binom.PP.EB(xh, nh, X=34,N=75, mix=TRUE) # EB power prior for new data 34/75.
#' }
binom.PP.EB <- function(x, n, X, N, verbose=FALSE, mc.cores=1, p.prior.a=1, p.prior.b=1, max.dn = FALSE, mix=FALSE){
  #if X isn't specified we calculate it for all, set flag too
  if(missing(X)) {
    X <- 0:N
    X.only <- FALSE
  } else X.only=TRUE

 ds <-
   mclapply(mc.cores=mc.cores,
            X,
            function(X){


    lik.d <- function(d) VGAM::dbetabinom.ab(X, N, p.prior.a+sum(d*x), p.prior.b+sum(d*(n-x)))

    # opd <- optimr::optimr(par = rep(.005, n.hist),
    #                       fn = lik.d,
    #                       lower=rep(0, n.hist),
    #                       upper=rep(1, n.hist),
    #                       method = "L-BFGS-B",
    #                       control=list(maximize=TRUE,
    #                                    fnscale=1.0e-20))
    opd <- BB::spg(par = rep(.005, length(x)),
               fn = lik.d,
               lower=rep(0, length(x)),
               upper=rep(1, length(x)),
               control=list(maximize=TRUE,
                            trace=FALSE))


    if(opd$convergence!=0) print(opd)

    
    
    d <- opd$par

    if(max.dn & sum(d*n) > max.dn ) d <- d * max.dn/sum(d*n)

    return(d)
  })
 
 if(mix){ #return mixture object
   ret.obj <- lapply(ds, function(d){
     obj <- create.mixture.prior("beta", pars = matrix(c(p.prior.a + sum(x*d),p.prior.b + sum(d*(n-x))),nrow=1))
     attr(obj, "powers") <- d
     obj
   })
    if(length(ret.obj)==1) ret.obj <- ret.obj[[1]]
return(ret.obj)
 
  } else{ #return simple function
   f <- function(p,X) {
     d <- ds[[X+1]]
     dbeta(p,p.prior.a + sum(x*d), p.prior.b + sum(d*(n-x)))
   }
   return(f) 
 }
}
