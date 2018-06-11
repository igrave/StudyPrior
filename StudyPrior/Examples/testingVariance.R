mean.variance <- matrix(0, nrow=101,ncol=9)
Nrep<-9
library(StudyPrior)
library(parallel)
sav.var <- 
mclapply(mc.cores=30,1:1000, function(i){
  set.seed(i)
  rnorm(100)
  print(paste("starting",i))
  het <- rnorm(1,0.65,0.05)
  nh <- c(100,100,100,100,100)
  xh <-rbinom(5, nh,p=het)
  
  
  xh/nh
  
  binom.PP.EB.Mix(xh, nh, N=100) -> a
  binom.PP.EB.Sep.Mix(xh, nh, N=100) -> b
  binom.PP.EB.Mix(sum(xh), sum(nh), N=100) -> c
  binom.PP.FB.COR(xh, nh, d.prior.cor = 0, mix=TRUE) -> d
  binom.PP.FB.COR(xh, nh, d.prior.cor = .5, mix=TRUE) -> e
  binom.PP.FB.COR(xh, nh, d.prior.cor = 1, mix=TRUE) -> f
  
  binom.MAP.FB(xh, nh) -> map.1
  conj.approx(distr = map.1,type = 'beta', max.degree = 10) -> g
  
  create.mixture.prior("beta",
                       pars = rbind(attr(g,"pars"),c(1,1)),
                       weights = c(0.5*attr(g,"weights"),0.5)) -> h

  create.mixture.prior("beta",pars=matrix(c(1,1),ncol=2))->i
  
  list.prior <- list(a,b,c,d,e,f,g,h,i)
  
  
  lapply(list.prior, function(priors){
    if(length(priors)>1) {
      lapply(0:100, function(x) posterior.mixture.prior(xs=x, n=100, mixture.prior = priors[[x+1]]))  
    }else{
      lapply(0:100, function(x) posterior.mixture.prior(xs=x, n=100, mixture.prior = priors))  
    }
  }) -> list.posterior
  
  
  
  sapply(seq_along(list.posterior), function(i) {
    print(i)
    sapply(list.posterior[[i]], var.mixture.prior )
  }) -> variances
  variances 
})

save(file="testingVariance.rda", sav.var)
for(i in 1:1000) mean.variance <- mean.variance+ sav.var[[i]]/1000

pdf()
png()
plot(sav.var[[1]][,7], ty='n', ylim=c(0,0.01))
for(i in 1:1000){
  lines(sav.var[[i]][,7], col=rgb(0,0,0,0.1))
  lines(sav.var[[i]][,8], col=rgb(0,0,1,0.1))
}
dev.off()

lapply(sav.var, function(x) any(x[,7]>0.006)) -> largevar

which(unlist(largevar)) -> strange.i

matplot(mean.variance, ty='l', col=1:9, lty=1)
legend("topleft", col=1:9, legend = letters[1:9], lty=1)

var.mixture.prior <- compiler::cmpfun(var.mixture.prior)

str(variances)
  matplot(variances, ty='l', col=1:8, lty=1, ylim = c(0,0.01))
  legend("topleft", col=1:8, legend = letters[1:8], lty=1)
  
system.time(sapply(1:101, function(i) var.mixture.prior(list.posterior[[4]][[i]]))) 
system.time(sapply(1:101, function(i) vp(list.posterior[[4]][[i]]))) 
vp <- compiler::cmpfun(var.mixture.prior)


#50.87s

var.mixture.prior(list.posterior[[4]][[30]])


bigsum*2
m <- outer(mw,mw)
diag(m) <- 0
sum(m)
