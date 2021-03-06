#' #' Empirical Bayes Meta-Analytic Prior for Binomial Data
#' #'
#' #' @param x number of historical successes
#' #' @param n number historical patients
#' #' @param verbose Print messages
#' #' @param mc.cores number of cores for parallel
#' #' @param X vector of new successes to computer prior for
#' #' @param N number of new patients
#' #' @param upper Upper limit of variance parameter
#' #'
#' #' @return A function of the probability parmater p
#' #' @examples \donttest{
#' #' xh <- c(30,40,50)
#' #' nh <- c(90,95,110)
#' #'  map <- binom.MAP.EB(xh,nh, X=34,N=75)
#' #'  }
#' #'
#' binom.MAP.EB <- function(x, n, X, N, verbose=FALSE, upper = 4,  mc.cores=1){
#'   n.hist <- length(x)
#'   n.new <- 0
#' 
#'   
#'   check.inla() #check for inla functions
#'   
#' 
#'  f <- mcmapply(mc.cores=mc.cores,
#'                X=X, N=N,
#'     FUN=function(X,N, upper=upper){
#' 
#'       if(!(missing(X)|missing(N))) {
#'         x <- c(x,X)
#'         n <- c(n,N)
#'         n.new <- 1
#'       }
#' 
#'       dat <- data.frame(x=x, n=n, z=1:(n.hist+n.new))
#' 
#' 
#' 
#'       if(missing(upper)) upper <- 4
#' 
#'       logiflatAB <- paste0("expression:
#'   upper = ",upper,";
#'   sd = exp(-x/2);
#'   A = (sd > upper ) ? 0.0000001 : 1/upper;
#'   B = abs(exp(-x/2)/2);
#'   logdens = log(A) + log(B);
#'   return(logdens)", collapse="")
#' 
#'       prior <-  list(prior= logiflatAB, initial=-3)
#'       # prior <-  list(prior= "logtnormal", param=c(0,1))
#'       # prior <-  list(prior= "flat", param=numeric(0))
#' 
#' 
#'    
#'       ##########################################
#' 
#'       result <- INLA::inla(x ~ 1 + f(z, model="iid", hyper = list(theta = prior)),
#'                            data = dat,
#'                            family = "binomial",
#'                            control.fixed = list(mean.intercept = 0, prec.intercept = 1/1000),
#'                            Ntrials=n,
#'                            verbose = verbose,
#'                            control.inla = list(int.strategy = "eb"))
#'       
#'       result <- INLA::inla.hyperpar(result)
#'       
#'       mode_tau <-  INLA::inla.mmarginal(INLA::inla.tmarginal(function(x) 1/x^.5,
#'                                                        result$marginals.hyperpar[[1]],
#'                                                        n=2000)[50:1950,])
#'       
#'       
#'       print(paste(X,mode_tau))
#'       
#'       
#'       VXN <- mode_tau^2
#'       
#'     
#' 
#'       dat <- data.frame(x=c(x[1:n.hist],NA), n=c(n[1:n.hist],NA), z=1:(n.hist+1))
#' 
#'       resultEB = INLA::inla(x ~ 1 + f(z, model="iid",
#'                                             hyper = list(theta = list(fixed=TRUE,
#'                                                                       initial=log(1/VXN)#1/mode_tau$maximum
#'                                             ))),
#'                             data = dat,
#'                             family = "binomial",
#'                             control.fixed = list(mean.intercept = 0, prec.intercept = 1/1000),
#'                             Ntrials=n,
#'                             control.predictor = list(compute=TRUE, link=1))
#' 
#'       #   plot( resultEB$marginals.fixed$`(Intercept)`)
#'       #   points( resultEB$marginals.fitted.values$fitted.predictor.4)
#'       #   plot( resultEB$marginals.fitted.values$fitted.predictor.4)
#'       #
#'       #
#'       #   # for some reason the link is wrong
#'       #   plot(inla.tmarginal(function(x) 1/(1+exp(-x)), resultEB$marginals.fixed$`(Intercept)`),xlim=c(0,1))
#'       # points(inla.tmarginal(function(x) 1/(1+exp(-x)), resultEB$marginals.fitted.values$fitted.predictor.4))
#'       ind <- round(resultEB$marginals.fitted.values[[n.hist+1]][,1],7)%in% c(0,1)
#' 
#'       A <- resultEB$marginals.fitted.values[[n.hist+1]][!ind,1]
#'       B <- resultEB$marginals.fitted.values[[n.hist+1]][!ind,2]
#' 
#' 
#'       return(list(X=A,Y=B))
#' 
#'       # f <- INLA:::inla.sfmarginal(inla.smarginal(marginal=list(x=A,y=B)))
#'       # 
#'       # 
#'       # function(p)
#'       # {
#'       #   n = length(p)
#'       #   d = numeric(n)
#'       #   for (i in 1:n) {
#'       #     if (p[i] >= f$range[1] && p[i] <= f$range[2]) {
#'       #       d[i] = exp(f$fun(p[i]))
#'       #     }
#'       #     else {
#'       #       d[i] = 0
#'       #     }
#'       #   }
#'       #   return(d)
#'       # }
#'       
#'       
#' })
#' 
#' 
#'  rm(x,n,X,N, n.hist, n.new, upper, verbose, mc.cores)
#' 
#'  function(p,X) {
#'    dens <- rep(0,length(p))
#'    i <- which(0<p&p<1)
#'    dens[i] <- splinefun(f[1,X+1][[1]], f[2,X+1][[1]])(p[i])
#'    dens
#'  }
#' }
