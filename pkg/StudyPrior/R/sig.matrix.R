
#' Calculate a matrix of significance test results for all possible values
#'
#' @param n.control Number of patients in the control group
#' @param n.treatment Number of patients in the treatment group
#' @param level Significance level to test (default 0.975)
#' @param prior Prior on binomial probability parameter of the control group
#' @param treat.beta.prior.par Paramters for beta prior on treatment group probability parameter
#' @param mc.cores Number of cores to use for mclapply
#' @param posterior List of posterior functions, one for each outcome in the control arm. (Default length=\code{n.control+1})
#' @param check.xt Optional argument to specify which outcomes in the treament arm should be considered. Otherwise \code{0:n.treatment}
#' @param check.xs Optional argument to specify which outcomes in the control arm should be considered. Otherwise \code{0:n.control}
#' @param approx Use \code{approxfun} to approximate the posterior if a mixture.prior is supplied. 
#'
#' @return A matrix of TRUE/FALSE values of dimensions (n.control+1, n.treatment+1) representing significant tests
#' @export
#'
#' Calculations are done using fast functions when \code{treat.beta.prior.par == c(1,1)} and the posterior is supplied as an \code{mixture.prior} object.

sig.matrix <- function(n.control, n.treatment, level=0.975, prior, posterior, treat.beta.prior.par=c(1,1), mc.cores=1, check.xt, check.xs, approx=TRUE) {

  if(missing(check.xt)) check.xt <- 0:n.treatment
  if(missing(check.xs)) check.xs <- 0:n.control

  use.posterior = !missing(posterior)

  
 if(treat.beta.prior.par == c(1,1) && use.posterior && inherits(posterior[[1]],"mixture.prior")){
   return(sigmat2(n.control, n.treatment, level, posterior, check.xs, check.xt))}
  
  ZZ.list <-   mclapply(check.xs,
    function(Xs){
      # print(Xs)
      post <-
        if(!use.posterior){
          if(inherits(prior, "function")){
            post <- function(p,g=1) prior(p,Xs)*dbinom(x=Xs, size=n.control, prob=p)/g
            f <- splinefun(smooth.spline(seq(0.001,.999,len=1000), pmax(0,post(seq(0.001,.999,len=1000)))))

            K <- adaptIntegrate(f, lowerLimit = 0, upperLimit = 1, maxEval = 2e5)$integral
            function(p,g=K) f(p)/g
            # formals(post) <- alist(p = , g = K)
            # K <- K *integrate(post, lower=0, upper=1)$value

          } else if(inherits(prior, "mixture.prior")){
            post.list <- posterior.mixture.prior(Xs, n.control,  mixture.prior=prior)
            function(p) eval.mixture.prior(p, post.list)

          } else if(inherits(prior, "list")){
            post.list <- posterior.mixture.prior(Xs, n.control,  mixture.prior=prior[[Xs+1]])
            function(p) eval.mixture.prior(p, post.list)
        }
        } else if(use.posterior) {
            if(inherits(posterior[[1]],"mixture.prior")) {
               if(approx){
                 approxfun(seq(0,1,len=300), eval.mixture.prior(seq(0,1,len=300), posterior[[Xs+1]]))
               } else  function(p) eval.mixture.prior(p, posterior[[Xs+1]])
            } 
          else posterior[[Xs+1]]
          }

    ZZ <-unlist(
      lapply(check.xt, function(xT){
       # print(paste0('.',xT))
      res <-  try({
        # if(debug ) browser()

        unsure <- TRUE
        #start quick and dirty
        tol <- 0.06

        this.int <- NA

        while(unsure){
          # this.int.res <-  adaptIntegrate(
          #   function(p0) {
          #     pmax(0,
          #         dbeta(p0,
          #               treat.beta.prior.par[1]+xT,
          #               treat.beta.prior.par[2] + n.treatment-xT) -
          #           post(p0))
          #     },
          #   lower=0,upper=1, rel.tol = tol, subdivisions = 200)

          
          this.int.res <-  adaptIntegrate(
            function(p0) {
              pbeta(p0,
                    treat.beta.prior.par[1]+xT,
                    treat.beta.prior.par[2] + n.treatment-xT,
                    lower.tail = FALSE) *
                post(p0)
            },
            lowerLimit = 0,upperLimit = 1, tol = tol, maxEval = 2e5)
          this.int.res$value <- this.int.res$integral
          this.int.res$abs.error <- this.int.res$integral*this.int.res$error

          this.int <- this.int.res$value

          if( (this.int - this.int.res$abs.error) < level & (this.int + this.int.res$abs.error)  > level){
            # print("We were very close! Trying again -------------------")
            # print(this.int)
            tol <- tol/3
            # print(paste0("Tol: ",tol))
          } else {
            unsure <- FALSE
          }

        }

        this.int > level
      })
      if(inherits(res,"try-error")) browser()#save(xT, post, file=paste0('ERROR-2-',Xs,'-',xT,'.Rda'))
      return(res)
    }
    # , mc.cores=mc.cores, mc.preschedule = FALSE
    ) )

    # print(ZZ)
# browser()
    ZZ
  }
  , mc.cores=mc.cores, mc.preschedule = FALSE
  )

  matrix(unlist(ZZ.list),
         byrow=TRUE,
         nrow=n.control+1,
         ncol=n.treatment+1,
         dimnames = list(Control = 0:n.control,
                         Treatment = 0:n.treatment)
  )

}




# A faster method
sigmat2 <- function(n.control, n.treatment, level, posterior, check.xs, check.xt){
  lapply(check.xs, function(xs){
    # print(paste("XS",xs))
    w <- attr(posterior[[xs+1]],"weights")
    par <- attr(posterior[[xs+1]],"pars")
    OK <- length(w) > 20
    
    sapply(check.xt, function(xt) {
      prXY(xt+1,n.treatment-xt+1, par[,1],par[,2], approx.allowed=OK, force=TRUE)  %*% w
    })-> prob 
    # browser()
    #these are too close
    # print(1-prob)
    redo <- (which(abs(1-prob-level)/level < 0.03))
    #recalculate without approximation

    prob[redo] <- sapply(check.xt[redo], function(xt) {
      prXY(xt+1,n.treatment-xt+1, par[,1],par[,2], approx.allowed=FALSE, force=FALSE)  %*% w
    })
    
    prob
    #end loop over mixture
  }) -> all.probs #end xs
  # browser()
  matrix(1-unlist(all.probs)>level, nrow=length(check.xs),ncol = length(check.xt),byrow = TRUE)}



#Analytical (and fast approximation for) calculation  of probability beta RV is larger than another
prXY <- function(a,b,c,d, approx.allowed=TRUE, force=FALSE){
  if(force | (a>5 &&b>10 && approx.allowed)){ 
    # cat("A")
    #for large values do a normal approximation based on moments
    pnorm((c/(c+d)-a/(a+b))/
          sqrt(a*b/((a+b)^2 * (a+b+1)) + c*d/((c+d)^2 * (c+d+1))))
  } else {
    #finite sum of beta functions
    sapply(1:length(c),function(i) {
      j <- a:(a+b-1)
         1 / beta(c[i],d[i]) * sum(choose(a+b-1,j)*beta(a+b+d[i]-j-1, c[i]+j))}
         )
    }
  }

