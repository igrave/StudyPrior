#'  Prior effective sample size
#'
#'  Estimate the prior effecive sample size for Beta-mixture
#'  priors based on Morita, Thall, Mueller (2008, Biometrics)
#'  evaluated at the prior mode.
#'
#' @param mixture A mixture.prior object
#' @param c  factor such that Beta( mu/c , (1-mu)/c ) is non0informative, where mu is mode of prior
#' @param LEN Number of points to evaluate at. Fewer is faster
#'
#' @return Estimate of the effective sample size.
#'

ess.morita.mixture.prior <-  function(mixture, c = 100, LEN=1000) {
  ## prior effective sample size ESS for Beta-mixture priors
  ## based on
  ## Morita, Thall, Mueller (MTM) 2008 Biometrics
  ## only difference: evaluated at mode of prior rather than at mean

  ## K: number of mixture components
  ## w: mixture weights w[k] k=1,...,K
  ## alpha, beta: components Beta(alpha[k],beta[k]) k=1,...,K
  ## c: factor such that Beta( mu/c , (1-mu)/c ) is noninformative
  ##    mu is mode of prior
  attrs <- attributes(mixture)
  K <- attrs$degree
  w <- attrs$weights
  alpha <- attrs$pars[,1]
  beta <- attrs$pars[,2]


  # initial estimate of prior mode
  p.eps=0.00001
  p = seq(0+p.eps,1-p.eps,p.eps)
  p.n = length(p)
  dens = rep(0,p.n)
  for (k in 1:K) {
    dens.k = w[k]*dbeta(p,alpha[k],beta[k])
    dens = dens + dens.k
  }
  mode = p[(1:p.n)[dens==max(dens)]]
  mode = mean(mode)
  if (mode < 0.001)     mode=0.001
  if (mode > (1-0.001)) mode=1-0.001

  # print("PART1")

  # 2nd derivative of log( prior ) at mode
  # numerical: local quadratic fit
  # refined estimate of mode
  delta=0.0001
  x = seq(mode-delta,mode+delta,length.out=LEN)
  xn = length(x)
  y = rep(0,xn)
  for (k in 1:K) {
    y.k = w[k]*dbeta(x,alpha[k],beta[k])
    y = y + y.k
  }
  y = log(y)
  mode = x[(1:xn)[y==max(y)]]
  mode = mean(mode)
  if (mode < 0.001)     mode=0.001
  if (mode > (1-0.001)) mode=1-0.001
  x = x - mode
  x2 = x^2
  fit = lm( y ~ x + x2 )
  deriv2.prior = fit$coefficients[3]

  # print("PART2")

  # 2nd derivative of log( posterior ) at mode
  # numerical: local quadratic fit
  # for all sample size up to ESSmax
  # for all possible number of responders y
  #
  ED2P = NULL
  ESSmax = round(max(alpha+beta))+10
  SS = 1:ESSmax


  alpha0 = mode / c
  beta0 = (1-mode)/c
  delta=0.001
  x = seq(mode-delta,mode+delta,length.out=LEN)
  xn = length(x)

  xmat <- matrix(1, ncol=3, nrow=xn)
  xmat[,2] <- x - mode
  xmat[,3] <- xmat[,2]^2


  for (m in SS) {
    y = 0:m
    deriv2.post = y*0

    # deriv2.post <-
    #   sapply(0:m, function(yn){
    #
    #     alpha0.post = alpha0 + yn
    #     beta0.post = beta0 + m - yn
    #     y = dbeta(x,alpha0.post,beta0.post, log=TRUE)
    #      # if(any(!is.finite(y))) browser()
    #     fit = .lm.fit( y=y, x=xmat)
    #     fit$coefficients[3]
    #   })
    #

    alpha0.post = alpha0 + y
    beta0.post = beta0 + m - y

    deriv2.post <-
      mapply(alpha0.post=alpha0.post,
             beta0.post=beta0.post,
             FUN=function(alpha0.post, beta0.post){


               y = dbeta(x,alpha0.post,beta0.post, log=TRUE)
               # if(any(!is.finite(y))) browser()
               fit = .lm.fit( y=y, x=xmat)
               fit$coefficients[3]
             })

    # print(paste(m,'/',ESSmax))

    #
    # prior predictive distr
    ypred = rep(0,m+1)
    for (yn in 0:m) {
      # for (k in 1:K) {
        # ypred.k = w[k]*dbetabinom.ab(yn,m,alpha[k],beta[k])
      ypred[yn+1]= sum(w*dbetabinom.ab(yn,m,alpha,beta))
        # ypred[yn+1] =  ypred.k
        # ypred[yn+1] = ypred[yn+1] + ypred.k
      # }
    }
    #
    # expected 2nd derivative of log( posterior ) at mode
    expect.deriv2.post = sum(ypred*deriv2.post)
    ED2P =c(ED2P,expect.deriv2.post)
  }
  #
  if (min(ED2P) >= deriv2.prior) ESS.MTM = ESSmax
  if (min(ED2P) <  deriv2.prior) ESS.MTM = min(SS[ED2P<deriv2.prior])
  #
  ESS.MTM
}
