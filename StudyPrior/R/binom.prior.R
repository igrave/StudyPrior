

#' Prior distribution for binomial data
#'
#' @param x Vector of historical events
#' @param n Vector of historical sample sizes
#' @param type Which prior to provide
#' @param ... Parameters to provide to prior
#' @param verbose Print messages
#'
#' @return A density function
#'
#' @examples \donttest{
#' xh <- c(30,40,50)
#' nh <- c(90,95,110)
#' eb <- binom.prior("PP.EB", xh, nh, X=34,N=75)
#' fb <- binom.prior("PP.FB", xh, nh)
#'
#' p <- seq(0,1,len=100)
#' plot(p,fb(p), ty='l', ylim=c(0,15))
#' lines(p, eb(p,X=0), col=2)
#' }
#'

binom.prior <- function(type = c("MAP.FB", "PP.FB","PP.EB","Beta.EB","PP.EB.Sep", "PP.Fix"),
                        x, n,
                        verbose=FALSE,
                        ...){

  type <- match.arg(type)

  f.to.call <- switch(type,
                      "MAP.FB" = binom.MAP.FB,
                      # "MAP.EB" = binom.MAP.EB,
                      "PP.FB" = binom.PP.FB.MC.BE,
                      "PP.EB" = binom.PP.EB,
                      # "PP.EB.Beta" = binom.PP.EB.Beta,
                      "PP.EB.Sep" = binom.PP.EB.Sep,
                      "PP.Fix" = binom.PP.FIX,
                      # "PP.Cor" = binom.PP.FB,
                      "Beta.EB" = binom.Beta.EB
                      )
v <- verbose

  do.call(f.to.call, args=alist(x=x,n=n, verbose=v, ...))

}
