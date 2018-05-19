#' Empirical Bayes Power Prior for Binomial Data with the power parameters esimated separately
#'
#' @param x number of historical successes
#' @param n number historical patients
#' @param X number of new successes
#' @param N number of new patients
#' @param verbose Print messages
#' @param p.prior.a shape1 parameter of initial beta prior for successes
#' @param p.prior.b shape2 parameter of initial beta prior for successes 
#' @param mix if TRUE return a mixture.prior object otherwise the function
#'
#' @return Either a list of mixture.prior objects which can be evaluated with eval.mixture.prior or a density function of the probability parmater p given X.
#' @export
#'
#'
binom.PP.EB.Sep <- function(x, n, X, N, verbose=FALSE, p.prior.a=1, p.prior.b=1, mix=FALSE){

  if(missing(X)) {
    X <- 0:N
    X.only <- FALSE
  } else X.only=TRUE
  
  
  n.hist <- length(x)

  ddbinom <- function(x, size, prob, delta) dbinom(x,size,prob)^delta
  ds <-
    lapply(X, function(X){
      lik.d <- function(d, xh, nh) VGAM::dbetabinom.ab(X, N, 1+(d*xh), 1+(d*(nh-xh)))
      opd <-sapply(1:n.hist, function(i){
        optimize(f = lik.d,
                 lower=0,
                 upper=1,
                 maximum = TRUE,
                 xh=x[i],
                 nh=n[i])$maximum
      })
      opd
})
  
  if(mix){ #return mixture object
    return(
      lapply(ds, function(d)
        create.mixture.prior("beta", pars = matrix(c(p.prior.a + sum(x*d), p.prior.b + sum(d*(n-x))),nrow=1)
        ))
    ) 
  } else{ #return simple function
    f <- function(p,X) {
      d <- ds[[X+1]]
      dbeta(p, p.prior.a + sum(x*d), p.prior.b + sum(d*(n-x)))
    }
    return(f) 
  }
}
