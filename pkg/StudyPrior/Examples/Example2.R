library(parallel)
library(foreach)
library(StudyPrior)
inla.setOption(num.threads=2)


clinical <- data.frame(Trial = c("Schorr", "Freire", "Kollef", "Ramirez"),
                       Year =  c(    2005,      2010,    2012,      2013),
                       CC =    c(      70,        71,      50,        18),
                       nCC =   c(     111,       122,      88,        24),
                       AM =    c(      19,        15,      13,         7),
                       nAM =   c(     111,       122,      88,        34))


nh <- clinical$nCC[-4]
xh <-clinical$CC[-4]
ns <- clinical$nCC[4]
i<-1
  set.seed(300*1)
  N <- 3
  n <- nh
  
  x <- xh
  
  
  Ns <- ns
  F.MFP <- binom.prior("MAP.FB", x = x, n=n)
  # F.MEP <- binom.prior("MAP.EB", x = x, n=n, X=0:Ns, N=Ns, mc.cores=1, verbose=FALSE)
  
  F.PFP <- binom.PP.FB.COR("PP.FB", x = x, n=n, mixture.size = 1000, d.prior.cor = 0)
  ## With feedback
  cat(i,"....\n")
  F.PEP <- binom.prior("PP.EB", x = x, n=n, X=0:Ns, N=Ns, verbose=FALSE, mc.cores=1)
  cat(i,".....\n")
  F.FX0 <- binom.prior("PP.Fix", x = x, n=n, d=0)
  cat(i,"......\n")
  F.COR <- binom.PP.FB.COR("PP.FB", x = x, n=n, mixture.size = 1000, d.prior.cor = 0.5)
  cat(i,".......\n")
  F.PSP <- binom.prior("PP.EB.Sep", x = x, n=n, X=0:Ns, N=Ns)
  cat(i,"........\n")
  F.SGL <- binom.prior("PP.EB", x = sum(x), n=sum(n), X=0:Ns, N=Ns, verbose=FALSE, mc.cores=1)
  cat(i,".........\n")

  F.CR9 <- binom.PP.FB.COR("PP.FB", x = x, n=n, mixture.size = 1000, d.prior.cor = 0.9)

  F.SUM <- binom.PP.FB.COR("PP.FB", x = sum(x), n=sum(n), mixture.size = 1000, d.prior.cor = 0)

    F.ROB <- conj.approx(distr = binom.prior("MAP.FB", x = x, n=n),
                       type = "beta",
                       robust = 0.5)
print("saving")
save(F.MFP, F.PFP,F.PEP,F.FX0, F.COR, F.PSP, F.SGL, n,x,F.CR9,F.SUM,F.ROB,
     file=paste0("models_example.rda"))



Calc.posterior.all <- function(prior, N, mc.cores){
  mclapply(0:N, function(X) calc.posterior(prior, X, N), mc.cores = mc.cores)
}



# library(foreach)
# library(doParallel)

####################################################
CORES <- 1


calc1 <- function(i){
  print(i)
  try({
    load(file=paste0("models_example.rda"))

    Ns 

    posteriors <-
    lapply(list(F.MFP, F.PFP,F.PEP,F.FX0, F.COR, F.PSP, F.SGL,F.CR9,F.SUM),
      function(p) Calc.posterior.all(prior=p, N=Ns, mc.cores=CORES)
    )

    mse <- lapply(posteriors,
                  function(pr) calc.MSE.mean(posterior=pr, prob.range=c(0,1), length = 100, mc.cores=CORES, n.binom=Ns))
    
    mse <- c(mse, list(
      calc.MSE.mean(prior=F.ROB, prob.range=c(0,1), length = 100, mc.cores=CORES, n.binom=Ns)
    ))
    
    tc <- system.time(
      bias <- c(lapply(posteriors,
                       function(pr) {
                         print("Starting!")
                         calc.bias(posterior = pr, prob.range=c(0,1), length = 100, n.binom=Ns, mc.cores=CORES)
                       }
      ), list(calc.bias(prior = F.ROB, prob.range=c(0,1), length = 100, n.binom=Ns, mc.cores=CORES))
      )
    )
    
    do.sigmat <- function(pr) sig.matrix(posterior = pr, n.control=Ns, n.treatment = Ns,
                                         check.xt=00:Ns, check.xs=0:Ns,
                                         level=0.975, mc.cores=CORES)
    
    
    SIGMAT <- c(lapply(posteriors, function(p){
      print("SIGMATING!")
      do.sigmat(p)
    } ),
    list(sig.matrix(prior = F.ROB, n.control=Ns, n.treatment = Ns,
                    check.xt=00:Ns, check.xs=0:Ns,
                    level=0.975, mc.cores=CORES))
    )
    
    
NT <- Ns

    pow <- mclapply(SIGMAT,
                    function(S) calc.power(sig.mat=S, n.binom.control = Ns,  n.binom.treatment = NT,
                                           prob.range = c(0,0.85), length =200, treatment.difference=0.12),
                    mc.cores=CORES)

    t1e <- mclapply(SIGMAT,
                    function(S) calc.power(sig.mat=S, n.binom.control = Ns,  n.binom.treatment = NT,
                                           prob.range = c(0,0.9), length = 200, treatment.difference = 0),
                    mc.cores=CORES)


    cover <-  mclapply(posteriors,
                       function(pr) calc.coverage(posterior=pr, level = 0.95, n.control = Ns, smooth = 0.03),
                       mc.cores=CORES)

    n.fix <- n
    x.fix <- x
    save(n.fix, x.fix,  bias,  cover, t1e, pow, mse, SIGMAT,  file=paste0("study_example.rda"))

  })
}

