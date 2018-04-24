#' Empirical Bayes Power Prior for Binomial Data Esimated Separately
#'
#' @param x number of historical successes
#' @param n number historical patients
#' @param X number of new successes
#' @param N number of new patients
#' @param verbose Print messages
#'
#' @return A function of the probability parmater p
#' @export
#'
#'
binom.PP.EB.Sep.Mix <- function(x, n, X, N, verbose=FALSE, p.prior.a=1, p.prior.b=1){

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
  
  
  lapply(ds, function(d) {
    create.mixture.prior("beta", 
                         pars = matrix(c(p.prior.a + sum(x*d),p.prior.b + sum(d*(n-x))),nrow=1))
  })
  
}
