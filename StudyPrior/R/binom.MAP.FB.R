#' Full Bayesian Meta-Analytic Prior for Binomial Data
#'
#' @param x number of historical successes
#' @param n number historical patients
#' @param tau.prior optional prior on heterogeneity parameter
#' @param verbose Print messages
#'
#' @return A function of the probability parmater p
#'
#' @examples \dontrun{
#' xh <- c(30,40,50)
#' nh <- c(90,95,110)
#' 
#' #This model uses INLA
#' library(INLA)
#' # Fit a MAP prior
#' map <- binom.MAP.FB(xh,nh, tau.prior = list(prior= "logtnormal", param=c(0,0.5)))
#' 
#' # Approximate with beta densities and add a 20% robust component
#' mix.map <- conj.approx(map, "beta", robust=0.2) 
#'  }
binom.MAP.FB <- function(x, n, tau.prior,verbose=FALSE){
  n.hist <- length(x)
  dat <- data.frame(x=c(x,NA), n=c(n,NA), z=1:(n.hist+1))

  
  check.inla() #check for inla functions
  
  
  if(missing(tau.prior)){
    prior <-  list(prior= "logtnormal", param=c(0,1))
  } else {
    prior <- tau.prior
  }

  formula <- 

  result <- INLA::inla(x ~ 1 + f(z, model="iid", hyper = list(theta = prior)),
                       data = dat,
                       family = "binomial",
                       control.fixed = list(mean.intercept = 0, prec.intercept = 1/sqrt(pi^2 /3)),
                       Ntrials=n,
                       control.predictor = list(compute=TRUE, link=1))


  # marg.string <- paste0(capture.output(dput(result$marginals.fitted.values[n.hist+1])),
                        # collapse = "")


#
#   fun.string <- paste0(
#     "function(p) {
#     marg <-",marg.string,"
#     inla.dmarginal(p, marg)
#   }", collapse="")
#
# f <- eval(parse(text=fun.string))

X <- result$marginals.fitted.values[[n.hist+1]][,1]
Y <- result$marginals.fitted.values[[n.hist+1]][,2]

 f <-  splinefun(c(0,min(X)/3, 2*min(X)/3, X, .3+.7*max(X),0.7+0.3*max(X),1),
                 c(0,0,0,Y,0,0,0),
                 method = "monoH")

 f <-  splinefun(c(0, X,1),
                 c(0,Y,0),
                 method = "monoH")
 
 
  rm(dat, formula, n, n.hist, prior, result, tau.prior, verbose, x)


   function(p,X) {
    dens <- rep(0,length(p))
    i <- which(0<p&p<1)
    dens[i] <- f(p[i])
    dens
  }

}
