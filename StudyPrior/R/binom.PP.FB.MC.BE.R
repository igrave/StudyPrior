#' Full Bayes Power Prior with Beta Prior on Weights
#'
#' @param x number of historical successes
#' @param n number historical patients
#' @param length  Number of points to evaluate density at
#' @param mc.cores Number of cores for parallel
#' @param samples Number of Monte Carlo samples
#' @param focus List of triples specifying regions to focus on, eg peaks. Specified as list(c(lower, upper, length))
#' @param verbose Print messages
#' @param d.prior.a shape1 parameter for beta prior on weights
#' @param d.prior.b shape2 parameter for beta prior on weights
#' @param p.prior.a shape1 parameter for beta prior on probability
#' @param p.prior.b shape2 parameter for beta prior on probability
#'
#' @return A density function
#' 
#' @examples \donttest{
#' x <- 21
#' n <- 37
#' 
#' #Full Bayes with delta~Be(1,1)
#' pp.fullbayes.11 <- binom.PP.FB.MC.BE(x = x, n = n, d.prior.a = 1, d.prior.b = 1)
#' 
#' #Full Bayes with delta~Be(.5,.5)
#' pp.fullbayes.55 <- binom.PP.FB.MC.BE(x = x, n = n, d.prior.a = .5, d.prior.b = .5)
#' 
#' }

binom.PP.FB.MC.BE <- function(x, n, verbose=FALSE, length=30,
                              d.prior.a=1, d.prior.b=1, p.prior.a=1, p.prior.b=1,
                              mc.cores=1, samples=10000, focus ){
  n.hist <- length(x)

  #Where to calculate density at
  p <- seq(0.0001, .9999, len=length)
  if(!missing(focus)){
    p2 <- unlist(sapply(focus, function(f) seq(from=f[1], to=f[2], length.out = f[3])))
    p <- unique(sort(c(p,p2)))
  }



  dens <- mclapply(p, function(PROB){
    # print(PROB)

    eval.f <- function(d) dbeta(PROB, p.prior.a+sum(d*x), p.prior.b+sum(d*(n-x)))

    mean(apply(
      matrix(rbeta(samples*n.hist, d.prior.a, d.prior.b),ncol=n.hist),
      1,eval.f))
  },mc.cores=mc.cores)

  f <- splinefun(smooth.spline(x=p, y=dens))

  k <- integrate(f, 0,1)$value

  f <- splinefun(x=p, y=unlist(dens)/k)

  g <- function(p,X) ifelse(0<=p&p<=1, f(p),0)

  return(g)

}
