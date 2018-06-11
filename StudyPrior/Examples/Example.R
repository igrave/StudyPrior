library(parallel)
library(StudyPrior)
inla.setOption(num.threads=2)


NN <- 24

if(TRUE){

  i<- 1
  N <- 3
  n <- c(111,122,88)
  
  x <- c(70,71,50)
  
  
  x/n
  Ns <- NN
  
  cat(i,".\n")
  F.MFP <- binom.prior("MAP.FB", x = x, n=n)
  cat(i,"..\n")
  F.PFP <- binom.PP.FB.COR("PP.FB", x = x, n=n, mixture.size = 2000, d.prior.cor = 0)
  ## With feedback
  cat(i,"....\n")
  F.PEP <- binom.prior("PP.EB", x = x, n=n, X=0:Ns, N=Ns, verbose=FALSE, mc.cores=1)
  cat(i,".....\n")
  F.FX0 <- binom.prior("PP.Fix", x = x, n=n, d=0)
  cat(i,"......\n")
  
  F.COR <- binom.PP.FB.COR("PP.FB", x = x, n=n, mixture.size = 2000, d.prior.cor = 0.5)
  cat(i,".......\n")
  F.PSP <- binom.prior("PP.EB.Sep", x = x, n=n, X=0:Ns, N=Ns)
  cat(i,"........\n")
  F.SGL <- binom.prior("PP.EB", x = sum(x), n=sum(n), X=0:Ns, N=Ns, verbose=FALSE, mc.cores=1)
  cat(i,".........\n")
  
  F.CR9 <- binom.PP.FB.COR("PP.FB", x = x, n=n, mixture.size = 2000, d.prior.cor = 0.9)
  print(".")
  
  F.SUM <- binom.PP.FB.COR("PP.FB", x = sum(x), n=sum(n), mixture.size = 2000, d.prior.cor = 0)
  print("..")
  F.ROB <- conj.approx(distr = binom.prior("MAP.FB", x = x, n=n),
                       type = "beta",
                       robust = 0.1)
  
  
  save(F.MFP, F.PFP,F.PEP,F.FX0, F.COR, F.PSP, F.SGL,F.CR9, F.SUM, F.ROB, n,x,
       file="models_ex_.rda")

}





calc <- function(i){
  print(i)
  
  
  clusterEvalQ(cl, {
    library(StudyPrior)
    library(parallel)
    mc.cores<- 1
    CORES <- 1
    load(file=paste0("models_ex_.rda"))
  
  })
  
  
  Nt <- i
  
  Ns <- NN
  CORES=1
  
  Calc.posterior.all <- function(prior, N){
    lapply(0:N, function(X) calc.posterior(prior, X, N))
  }
  
  posteriors <- parLapply(cl,list(F.MFP, F.PFP,F.PEP,F.FX0, F.COR, F.PSP, F.SGL,F.CR9,F.SUM),
                      Calc.posterior.all, N=Ns)
  
  for.clust <- new.env()
  for.clust$posteriors <- posteriors
  for.clust$Nt <- i
  for.clust$Ns <- NN
  
  clusterExport(cl, varlist=c("Nt","Ns","posteriors"), envir = for.clust)
    
  print("mse")
    
    mse <- parLapply(cl,posteriors,
                  function(pr) calc.MSE.mean(posterior=pr, prob.range=c(0,1), length = 100, mc.cores=CORES, n.binom=Ns))
    
    mse <- c(mse, list(
      calc.MSE.mean(prior=F.ROB, prob.range=c(0,1), length = 100, mc.cores=CORES, n.binom=Ns)
    ))

    print("bias")
    
    tc <- system.time(
      bias <- c(parLapply(cl,posteriors,
                     function(pr) {
                       print("Starting!")
                       calc.bias(posterior = pr, prob.range=c(0,1), length = 100, n.binom=Ns, mc.cores=CORES)
                     }
      ), list(calc.bias(prior = F.ROB, prob.range=c(0,1), length = 100, n.binom=Ns, mc.cores=CORES))
      )
      )

    do.sigmat <- function(pr) sig.matrix(posterior = pr, n.control=Ns, n.treatment = Nt,
                                         check.xt=00:Nt, check.xs=0:Ns,
                                         level=0.975, mc.cores=CORES)

    print("sigmat")
    
    SIGMAT <- c(parLapply(cl,posteriors, function(p){
      print("SIGMATING!")
      do.sigmat(p)
    } ),
    list(sig.matrix(prior = F.ROB, n.control=Ns, n.treatment = Nt,
                    check.xt=0:Nt, check.xs=0:Ns,
                    level=0.975, mc.cores=CORES))
    )


print("power")
    
    pow <- parLapply(cl,SIGMAT,
                    function(S) calc.power(sig.mat=S, n.binom.control = Ns, n.binom.treatment = Nt,
                                           prob.range = c(0,0.85), length =200, treatment.difference=0.12))

    print("t1e")
    
    t1e <- parLapply(cl,SIGMAT,
                    function(S) calc.power(sig.mat=S, n.binom.control = Ns, n.binom.treatment = Nt,
                                           prob.range = c(0,0.9), length = 200, treatment.difference = 0))

    print("cover")
    
    cover <-  c(parLapply(cl,posteriors,
                       function(pr) calc.coverage(posterior=pr, level = 0.95, n.control = Ns, smooth = 0.03)),
                list(calc.coverage(prior = F.ROB, level = 0.95, n.control = Ns, smooth = 0.03)))

    n.fix <- n
    x.fix <- x
    save(n.fix, x.fix,  bias,  cover, t1e, pow, mse, SIGMAT,  file=paste0("study_ex_",i,".rda"))

  }


stopCluster(cl)

cl <- makeCluster(4)

system.time(calc(24))

# system.time(calc(50))


par(mfrow=c(2,2))
load(file="study_ex_24.rda")
COLS <- c( rainbow(7),gray(seq(0,.75,len=3)))


p1 <- seq(0,0.85, len=200)
plot(p1,t1e[[1]],ty='n', ylim=c(0,.2), xlab='p', ylab="Type 1 Error")
for(i in seq_along(t1e)) lines(p1, t1e[[i]], col=COLS[i], lwd=2)
legend("topleft", lty=1, col=COLS,
       legend = c("F.MFP", "F.PFP","F.PEP","F.FX0", "F.COR", "F.PSP", "F.SGL","F.CR9","F.SUM","F.ROB"), lwd=2)
abline(h=0.025)

p2 <- seq(0,0.9, len=200)
plot(p2,pow[[1]],ty='n', ylim=c(0,1), xlab='p', ylab="Power")
for(i in seq_along(pow)) lines(p2,pow[[i]], col=COLS[i], lwd=2)
legend("topleft", lty=1, col=COLS,
       legend = c("F.MFP", "F.PFP","F.PEP","F.FX0", "F.COR", "F.PSP", "F.SGL","F.CR9","F.SUM","F.ROB"), lwd=2)
abline(h=0.025)


############################################################
# load(file="study_ex_50.rda")
# # COLS <- rainbow(10)
# 
# # par(mfrow=c(1,2))
# p1 <- seq(0,0.85, len=200)
# plot(p1,t1e[[1]],ty='n', ylim=c(0,.2), xlab='p', ylab="Type 1 Error")
# for(i in seq_along(t1e)) lines(p1, t1e[[i]], col=COLS[i], lwd=2)
# legend("topleft", lty=1, col=COLS,
#        legend = c("F.MFP", "F.PFP","F.PEP","F.FX0", "F.COR", "F.PSP", "F.SGL","F.CR9","F.SUM","F.ROB"), lwd=2)
# abline(h=0.025)
# 
# p2 <- seq(0,0.9, len=200)
# plot(p2,pow[[1]],ty='n', ylim=c(0,1), xlab='p', ylab="Power")
# for(i in seq_along(pow)) lines(p2,pow[[i]], col=COLS[i], lwd=2)
# legend("topleft", lty=1, col=COLS,
#        legend = c("F.MFP", "F.PFP","F.PEP","F.FX0", "F.COR", "F.PSP", "F.SGL","F.CR9","F.SUM","F.ROB"), lwd=2)
# abline(h=0.025)
