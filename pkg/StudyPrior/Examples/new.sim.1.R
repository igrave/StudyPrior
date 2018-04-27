# New simulation


library(parallel)
library(StudyPrior)
inla.setOption(num.threads=2)




NH <- 100
NN <- 75
NT <- 200

gen.models <- function(i){
  
  print(i)
  
  set.seed(2019)
  rnorm(500*i)
  
  N <- 5
  n <- rep(NH,N)
  
  z <- rnorm(N, 0.65, 0.05)
  
  x <- mapply(rbinom, size=n, n=1, prob=z)
  
  #EB Combined
  F.PEP <- binom.prior("PP.EB", x = x, n=n, X=0:NH, N=NH, verbose=FALSE, mc.cores=1)
  #EB pool
  F.SGL <- binom.prior("PP.EB", x = sum(x), n=sum(n), X=0:NH, N=NH, verbose=FALSE, mc.cores=1)
  #EP separate
  F.PSP <- binom.prior("PP.EB.Sep", x = x, n=n, X=0:NH, N=NH)
  
  #PP FB
  F.PFP <- binom.PP.FB.COR("PP.FB", x = x, n=n, mixture.size = 1000, d.prior.cor = 0)
  ## With feedback
  F.COR <- binom.PP.FB.COR("PP.FB", x = x, n=n, mixture.size = 1000, d.prior.cor = 0.5)
  
  #PP FB pool
  F.SUM <- binom.PP.FB.COR("PP.FB", x = sum(x), n=sum(n), mixture.size = 1000, d.prior.cor = 0)
  
  #MAP
  F.MFP <- binom.prior("MAP.FB", x = x, n=n)
  
  F.ROB <- conj.approx(distr = binom.prior("MAP.FB", x = x, n=n),
                       type = "beta",
                       robust = 0.5)
  
  # No historical
  F.FX0 <- binom.prior("PP.Fix", x = x, n=n, d=0)
  
  
  
  
  
  save(F.PEP, F.SGL, F.PSP, F.PFP, F.COR, F.SUM, F.MFP, F.ROB, x,n,
       file=paste0("new.sim.1/models_",i,".rda"))
  
}


mclapply(1:9, gen.models, mc.cores=3)