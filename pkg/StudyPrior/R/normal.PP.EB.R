
#' Empirical Bayes power prior for normal likelihoods 
#'
#' @param x vector of historical means
#' @param sd vector of historical standard deviations
#' @param X new study mean
#' @param SD new study standard deviation
#'
#' @return A density function of the mean parameter
#'
normal.PP.EB <- function(x, sd, X, SD){

  r <- 1 / sd^2
  R<- 1 / SD^2
  
  marg <- function(d){
    dr <- d*r
    sumdr <- sum(dr)
    sumdrx <- sum(d*r*x)
    sumdrX <- sum(d*r*X)
    
    exp( - (R*(sumdrx - sumdrX)^2) /
             (2*(sumdr)*(R+sumdr))) *
      sqrt( (sumdr) / (R+sumdr) )
  }
  
  n.hist <- length(x)
  
  opt <- 
  optimx::optimx( rep(0.8, n.hist),
                  marg,
                  control =list(maximize=TRUE),
                  method = "L-BFGS-B",
                  lower = rep(0, n.hist),
                  upper = rep(1, n.hist))

  d <- as.numeric(opt[1,seq(n.hist)])
  
# See http://www.tina-vision.net/docs/memos/2003-003.pdf for product of gaussian densities

  # eq (4) and (5)
all.sd2 <- 1/sum(r*d)
all.mean <- sum(x*r*d) * all.sd2
  

# The density:
  fun <- function(mu) dnorm(mu,all.mean, sqrt(all.sd2))  
  
  attr(fun, "parameters") <- list(x = x, sd = sd,
                                  X = X, SD = SD,
                                  d = d,
                                  all.mean = all.mean,
                                  all.sd2 = all.sd2,
                                  optim = opt)
  
  fun
} 

