# RandomModelSim

# New simulation


library(parallel)
library(StudyPrior)
library(INLA)
inla.setOption(num.threads=2)



NH <- 100
NN <- 75
NT <- 200


gen.models <- function(I){
  
  print(I)
  
  set.seed(2019)
  rnorm(500*I)
  
  N <- 5
  nh <- rep(NH,N)
  
  z <- rnorm(N, 0.65, 0.05)
  
  xh <- mapply(rbinom, size=nh, n=1, prob=z)
  
  binom.PP.EB(xh, nh, N=NN, mix=TRUE) -> a
  binom.PP.EB.Sep(xh, nh, N=NN, mix=TRUE) -> b
  binom.PP.EB(sum(xh), sum(nh), N=NN, mix=TRUE) -> c
  
  
  lapply(a, function(A) create.mixture.prior("beta",
                                             pars = rbind(attr(A,"pars"),c(1,1)),
                                             weights = c(0.5*attr(A,"weights"),0.5))) -> d
  
  binom.PP.EB(xh, nh, N=NN, mix=TRUE, max.dn = 200-75) -> e
  
  binom.PP.EB(xh, nh, N=NN, mix=TRUE, max.dn = 50) -> f
  
  # binom.PP.FB.COR(xh, nh, d.prior.cor = 0, mix=TRUE) -> d
  # binom.PP.FB.COR(xh, nh, d.prior.cor = .5, mix=TRUE) -> e
  # binom.PP.FB.COR(xh, nh, d.prior.cor = 1, mix=TRUE) -> f
  
  binom.MAP.FB(xh, nh) -> map.1
  conj.approx(distr = map.1,type = 'beta', max.degree = 10) -> g
  
  create.mixture.prior("beta",
                       pars = rbind(attr(g,"pars"),c(1,1)),
                       weights = c(0.5*attr(g,"weights"),0.5)) -> h
  
  create.mixture.prior("beta",pars=matrix(c(1,1),ncol=2))->k
  
  list.prior <- list(a,b,c,d,e,f,g,h,k)
  
  lapply(list.prior, function(priors){
    if(length(priors)>1) {
      lapply(0:NN, function(x) posterior.mixture.prior(xs=x, n=NN, mixture.prior = priors[[x+1]]))  
    }else{
      lapply(0:NN, function(x) posterior.mixture.prior(xs=x, n=NN, mixture.prior = priors))  
    }
  }) -> list.posterior
 
  
  save(file=paste0('models_',formatC(I, width=4, flag="0"),'.rda'),
       list.prior, list.posterior,xh,nh,  NH, NN, NT )
  
}





calc.oc <- function(I){
  try( {load(file=paste0('models_',formatC(I, width=4, flag="0"),'.rda'))
  
    
   ess <- lapply(list.posterior, function(p) sapply(p, ess.mixture.prior))
    
    mse <- lapply(list.posterior,
                  function(pr) calc.MSE.mean(posterior=pr, prob.range=c(0,1), length = 100,  n.binom=NN))
    
  
    bias <- lapply(list.posterior,
                   function(pr) calc.bias(posterior = pr, prob.range=c(0,1), length = 100, n.binom=NN))
        
    SIGMAT <- lapply(list.posterior, function(pr) sig.matrix(posterior = pr, n.control=NN, n.treatment = NT,
                                                   check.xt=00:NT, check.xs=0:NN,
                                                   level=0.975))
    
    
    pow <- lapply(SIGMAT,
                    function(S) calc.power(sig.mat=S, n.binom.control = NN, n.binom.treatment = NT,
                                           prob.range = c(0,0.85), length =200, treatment.difference=0.12))
    
    lr <- lapply(SIGMAT,
                    function(S) calc.power(sig.mat=S, n.binom.control = NN, n.binom.treatment = NT,
                                           prob.range = c(0,0.85), length =200, treatment.difference=0.12, LR=TRUE))
    
    t1e <- lapply(SIGMAT,
                    function(S) calc.power(sig.mat=S, n.binom.control = NN, n.binom.treatment = NT, 
                                           prob.range = c(0,0.85), length = 200, treatment.difference = 0))
 
    
    save(file=paste0('oc_',formatC(I, width=4, flag="0"),'.rda'),
         ess,mse, bias, SIGMAT, pow, t1e,lr, xh,nh,  NH, NN, NT )
  })
}


mclapply(1:200, gen.models, mc.cores=20, mc.preschedule = FALSE)
# mclapply(1:1000, calc.oc, mc.cores=30)
  

TF2 <- rep(FALSE,200)
for(I in 1:200) TF2[I] <-  !file.exists(paste0('models_',formatC(I, width=4, flag="0"),'.rda'))
recalc2 <- which(TF2)
recalc2
mclapply(recalc2, gen.models, mc.cores=20)
#   


TF3 <- rep(FALSE,200)
for(I in 1:200) TF3[I] <-  !file.exists(paste0('oc_',formatC(I, width=4, flag="0"),'.rda'))
recalc3 <- which(TF3)
mclapply(recalc3, calc.oc, mc.cores=30, mc.preschedule = FALSE)


load(paste0('oc_0001.rda'))
pow.all <- pow
t1e.all <- t1e
for(I in 2:200) {
  load(paste0('oc_',formatC(I, width=4, flag="0"),'.rda'))
  for(i in 1:9){
    pow.all[[i]] <- pow[[i]]+pow.all[[i]]
    t1e.all[[i]] <- t1e[[i]]+t1e.all[[i]]
  }
}
