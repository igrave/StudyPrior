#' Full Bayes Bias Prior for Binomial Data
#'
#' @param x number of historical successes
#' @param n number historical patients
#' @param verbose Print messages
#' @param mc.cores number of cores for parallel
#'
#' @return A function of the probability parmater p
#'
#'
binom.Bias.FB <- function(x, n, verbose=FALSE,  mc.cores=1){
  n.hist <- length(x)

  dat <- data.frame(x, n, z=1:n.hist)
  
  check.inla() #check for inla functions
  
  # if(missing(tau.prior)){
    prior <-  list(prior= "logtnormal", param=c(0,1))
  # } else {
  #   prior <- tau.prior
  # }

    f<- INLA::f
    


  result <- INLA::inla(x ~ 1 + f(z, model="iid", hyper = list(prec = prior)),
                       data = dat,
                       family = "binomial",
                       control.fixed = list(mean.intercept = 0, prec.intercept = 1/1000),
                       Ntrials=n,
                       control.predictor = list(compute=TRUE, link=1))

k <- integrate(Vectorize(function(p) prod(sapply(result$marginals.fitted.values, INLA::inla.dmarginal, x=p))),
               lower=0,
               upper=1)$value


  f <- function(p) sapply(p, function(p) prod(sapply(result$marginals.fitted.values, INLA::inla.dmarginal, x=p))/k)



}
