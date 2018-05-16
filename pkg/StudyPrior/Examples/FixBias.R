


fix.models <- function(I){
  try( {load(file=paste0('models_',formatC(I, width=4, flag="0"),'.rda'))
    
    print(I)
    
    binom.PP.EB.Mix(xh, nh, N=NN) -> a
    binom.PP.EB.Sep.Mix(xh, nh, N=NN) -> b
    binom.PP.EB.Mix(sum(xh), sum(nh), N=NN) -> c
    
    list.prior[[1]] <- a
    list.prior[[2]] <- b
    list.prior[[3]] <- c
    
    
    lapply(list.prior[1:3], function(priors){
      if(length(priors)>1) {
        lapply(0:NN, function(x) posterior.mixture.prior(xs=x, n=NN, mixture.prior = priors[[x+1]]))  
      }else{
        lapply(0:NN, function(x) posterior.mixture.prior(xs=x, n=NN, mixture.prior = priors))  
      }
    }) -> list.posterior13
    
    list.posterior[1] <- list.posterior13[1]
    list.posterior[2] <- list.posterior13[2]
    list.posterior[3] <- list.posterior13[3]

    
    save(file=paste0('models_fix_',formatC(I, width=4, flag="0"),'.rda'),
         list.prior, list.posterior,xh,nh,  NH, NN, NT )
    
  
        
  })}
    





fix.oc <- function(I){
  try( {load(file=paste0('models_fix_',formatC(I, width=4, flag="0"),'.rda'))
    load(file=paste0('oc_',formatC(I, width=4, flag="0"),'.rda'))
    
    ess <- lapply(list.posterior, function(p) sapply(p, ess.mixture.prior))
    
    mse[1:3] <- lapply(list.posterior[1:3],
                  function(pr) calc.MSE.mean(posterior=pr, prob.range=c(0,1), length = 100,  n.binom=NN))
    
    
    bias[1:3] <- lapply(list.posterior[1:3],
                   function(pr) calc.bias(posterior = pr, prob.range=c(0,1), length = 100, n.binom=NN))
    
    SIGMAT[1:3] <- lapply(list.posterior[1:3], function(pr) sig.matrix(posterior = pr, n.control=NN, n.treatment = NT,
                                                             check.xt=00:NT, check.xs=0:NN,
                                                             level=0.975))
    
    
    pow[1:3] <- lapply(SIGMAT[1:3],
                  function(S) calc.power(sig.mat=S, n.binom.control = NN, n.binom.treatment = NT,
                                         prob.range = c(0,0.85), length =200, treatment.difference=0.12))
    
    lr[1:3] <- lapply(SIGMAT[1:3],
                 function(S) calc.power(sig.mat=S, n.binom.control = NN, n.binom.treatment = NT,
                                        prob.range = c(0,0.85), length =200, treatment.difference=0.12, LR=TRUE))
    
    t1e[1:3] <- lapply(SIGMAT[1:3],
                  function(S) calc.power(sig.mat=S, n.binom.control = NN, n.binom.treatment = NT, 
                                         prob.range = c(0,0.9), length = 200, treatment.difference = 0))
    
    save(file=paste0('oc_',formatC(I, width=4, flag="0"),'.rda'),
         ess,mse, bias, SIGMAT, pow, t1e,lr, xh,nh,  NH, NN, NT )  
  })
}



TF3 <- rep(FALSE,1000)
for(I in 1:1000) TF3[I] <-  !file.exists(paste0('models_fix_',formatC(I, width=4, flag="0"),'.rda'))
recalc3 <- which(TF3)
mclapply(recalc3, fix.models, mc.cores=5)



TF3 <- rep(FALSE,1000)
for(I in 1:1000) TF3[I] <-  !file.exists(paste0('oc_',formatC(I, width=4, flag="0"),'.rda'))
recalc3 <- which(TF3)
mclapply(recalc3, fix.oc, mc.cores=37)








  nmodel <- n.models<- 9
  ratio.all <- lr.all  <- t1e.all <- pow.all <- matrix(0, ncol=nmodel,nrow=200)
  bias.all <- mse.all <- matrix(0, ncol=nmodel,nrow=100)
  
  ess.all <- matrix(0, ncol=nmodel,nrow=75+1)
  cover.all <- matrix(0, ncol=nmodel,nrow=3999)
  n.used <- 0
  for(I in 1:700){
    if(file.exists(paste0('oc_',formatC(I, width=4, flag="0"),'.rda'))){
      load(paste0('oc_',formatC(I, width=4, flag="0"),'.rda'))
      i <-I
      print(i)
      
      ratio.list <- list()
      for(i in 1:n.models){
        x=seq(0,0.85,len=200)
        t1efun <- splinefun(x=seq(0,0.9,len=200), y=t1e[[i]])
        powfun <- splinefun(x=seq(0,0.85,len=200), y=pow[[i]])
        ratio.list[[i]] <- powfun(x)/t1efun(x)
      }
      
      ratio.all  <- matrix(unlist(ratio.list),ncol=nmodel)+ratio.all
      lr.all <- matrix(unlist(lr),ncol=nmodel)+lr.all
      pow.all <- matrix(unlist(pow),ncol=nmodel)+pow.all
      t1e.all <- matrix(unlist(t1e),ncol=nmodel)+t1e.all
      bias.all <- matrix(unlist(bias),ncol=nmodel)+bias.all
      mse.all <- matrix(unlist(mse),ncol=nmodel)+mse.all
      ess.all <- matrix(unlist(ess),ncol=nmodel)+ess.all
      
      # print(paste(ess[[1]][30], ess.all[30,1]/i))  
      n.used <- n.used + 1
    }
    
  }
  
  ess.list <- lapply(1:nmodel, function(i) ess.all[,i]/n.used)
  mse.list <- lapply(1:nmodel, function(i) mse.all[,i]/n.used)
  t1e.list <- lapply(1:nmodel, function(i) t1e.all[,i]/n.used)
  pow.list <- lapply(1:nmodel, function(i) pow.all[,i]/n.used)
  cover.list <- lapply(1:nmodel, function(i) cover.all[,i]/n.used)
  bias.list <- lapply(1:nmodel, function(i) bias.all[,i]/n.used)
  lr.list <- lapply(1:nmodel, function(i) lr.all[,i]/n.used)
  ratio.list <- lapply(1:nmodel, function(i) ratio.all[,i]/n.used)
  save(ess.list, mse.list, t1e.list, pow.list, cover.list, bias.list, lr.list, ratio.list,
       file=paste0("Study_new.rda"))
  
  




















