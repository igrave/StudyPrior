
#' Calculate power
#'
#' @param prior Prior to calculate posterior
#' @param prob.range Range of values to calculate over
#' @param length Number of values to calculate for
#' @param n.binom.control Number of patients in new trial's contral arm
#' @param n.binom.treatment Number of patients in new trial's treatment arm
#' @param treatment.difference Predefined treatment difference
#' @param level Significance level for calculation of sig.mat
#' @param LR Calculate the likelihood ratio of power/type I error
#' @param sig.mat Precalculated significance matrix (See sig.matrix)
#'
#' @return Vector or power values
#' @examples \donttest{
#' xh <- c(30,40,50)
#' nh <- c(90,95,110)
#' fix <- binom.PP.FIX(x=xh, n=nh, d=c(.2,.5,.7), mix=TRUE) # set different weights
#' mat <- sig.matrix(n.control=50, n.treatment=75, prior = fix)
#' 
#' # Calculate the power to detect difference of 0.15
#' pow <- calc.power(sig.mat=mat, prob.range=c(.5,.85),
#'   n.binom.control = 50, n.binom.treatment=75, treatment.difference=0.15)
#' 
#' # Calculate the type 1 error by setting treatment.difference = 0
#' pow <- calc.power(sig.mat=mat, prob.range=c(.5,.85),
#'   n.binom.control = 50, n.binom.treatment=75, treatment.difference=0)
#' 
#' }
calc.power <- function(prior, prob.range=c(.5,.9), length=20, n.binom.control=30,
                       n.binom.treatment=n.binom.control, treatment.difference=0.1, level=0.95, sig.mat, LR=FALSE){
# the probability to detect a given true difference
# detection by
  if(treatment.difference+prob.range[2] > 1) stop("Unable to calculate power for difference at upper end of probability range (prob.range+treatment.difference > 1).")

  if(missing(sig.mat)) sig.mat <- sig.matrix(n.binom.control,
                                             n.binom.treatment, level, prior,
                                             treat.beta.prior.par=c(1,1))


  
  
  probs <- sapply(seq(prob.range[1],prob.range[2], len=length),
                  function(P){
                    Mc <- matrix(dbinom(0:n.binom.control, n.binom.control, prob=P),nrow=1)
                    Mt <- if(LR){
                      matrix(dbinom(0:n.binom.treatment, n.binom.treatment, prob=P+treatment.difference)/
                               dbinom(0:n.binom.treatment, n.binom.treatment, prob=P+0),ncol=1)
                    } else matrix(dbinom(0:n.binom.treatment, n.binom.treatment, prob=P+treatment.difference),ncol=1)
                    
     Mc %*% sig.mat %*% Mt
      
  })

  names(probs) <- seq(prob.range[1],prob.range[2], len=length)

  return(probs)
}
