library(parallel)
library(foreach)
library(StudyPrior)
inla.setOption(num.threads=2)
recalc <- 1:500

NH <- 100
NN <- 75
NT <- 200

if(FALSE){
mclapply(mc.cores=20, recalc, function(i){
  # lapply( recalc, function(i){
  # for(i in recalc){
  print(i)

  set.seed(300*i)
  N <- 5
  n <- rep(NH,N)
  
  z <- rnorm(N, 0.65, 0.05)
  
  x <- mapply(rbinom, size=n, n=1, prob=z)
  
  
  x/n
  Xs <- rbinom(1,100,0.6)
  Ns <- NN
  cat(i,".\n")
  F.MFP <- binom.prior("MAP.FB", x = x, n=n)
  cat(i,"..\n")
  # F.MEP <- binom.prior("MAP.EB", x = x, n=n, X=0:Ns, N=Ns, mc.cores=1, verbose=FALSE)
  cat(i,"...\n")
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
  save(F.MFP, F.PFP,F.PEP,F.FX0, F.COR, F.PSP, F.SGL, n,x,
       file=paste0("Rand-75/models_f_",i,".rda"))
}
)

}

if(FALSE){
recalc <- 1:500
mclapply(mc.cores=20, recalc, function(i){
  # lapply( recalc, function(i){
  # for(i in recalc){
  print(i)

  set.seed(300*i)
  N <- 5
  n <- rep(NH,N)
  
  z <- rnorm(N, 0.65, 0.05)
  
  x <- mapply(rbinom, size=n, n=1, prob=z)
  Xs <- rbinom(1,100,0.6)
  
  Ns <-NN
  F.CR9 <- binom.PP.FB.COR("PP.FB", x = x, n=n, mixture.size = 1000, d.prior.cor = 0.9)
print(".")
  # F.C95 <- binom.prior("PP.Cor", x = x, n=n, d.prior.cor=0.95, samples=1000, length=512)
  F.SUM <- binom.PP.FB.COR("PP.FB", x = sum(x), n=sum(n), mixture.size = 1000, d.prior.cor = 0)
print("..")
  F.ROB <- conj.approx(distr = binom.prior("MAP.FB", x = x, n=n),
                       type = "beta",
                       robust = 0.1)
print("saving")
  save(F.CR9,F.SUM,F.ROB, n, x,
       file=paste0("Rand-75/models_g_",i,".rda"))
})
}
#
# mclapply(mc.cores=25, recalc, function(i){
#   load(file=paste0("Rand-5/models_f_",i,".rda"))
#   print(i)
#
#   set.seed(300*i)
#   N <- 5
#   n <- rep(50,N)
#   z <- 0.65
#   x <- mapply(rbinom, size=n, n=1, prob=z)
#   x/n
#   Xs <- rbinom(1,100,0.6)
#   Ns <- 200
#   cat(i,".\n")
#   F.MFP <- binom.prior("MAP.FB", x = x, n=n)
#   cat(i,"..\n")
#   F.MEP <- binom.prior("MAP.EB", x = x, n=n, X=0:Ns, N=Ns, mc.cores=1, verbose=FALSE)
#   F.PFP <- binom.prior("PP.FB", x = x, n=n, samples=5000, length=512, mc.cores=1, verbose=FALSE)
#   save(F.MFP, F.MEP,F.PFP,F.PEP,F.FX0, F.COR, F.PSP, F.SGL, n,x,
#        file=paste0("Rand-5/models_f_",i,".rda"))
# }
# )


if(TRUE){

Calc.posterior <- function(prior, X, N){
  pr <- function(p) prior(p,X)
  k <- adaptIntegrate(function(p) pr(p)*dbinom(X,N,p),
                      lowerLimit = 0,
                      upperLimit = 1)$integral
  function(p) pr(p)*dbinom(X,N,p)/k
}

Calc.posterior.all <- function(prior, N, mc.cores){
  mclapply(0:N, function(X) Calc.posterior(prior, X, N), mc.cores = mc.cores)
}



# library(foreach)
# library(doParallel)

####################################################
CORES <- 1
recalc <- 1:500

calc1 <- function(i){
  print(i)
  try({
    load(file=paste0("Rand-75/models_f_",i,".rda"))

    Ns <-NN

    posteriors <-
    lapply(list(
      F.MFP, F.PFP,F.PEP,F.FX0, F.COR, F.PSP, F.SGL),
      Calc.posterior.all, N=Ns, mc.cores=CORES
    )

    # posteriors2 <-
    #   lapply(list(
    #     F.MFP),
    #     Calc.posterior.all, N=Ns, mc.cores=CORES
    #   )

    mse <- lapply(posteriors,
      function(pr) calc.MSE.mean(posterior=pr, prob.range=c(0,1), length = 100, mc.cores=CORES, n.binom=Ns))


    tc <- system.time(
      bias <- lapply(posteriors,
        function(pr) {
         print("Starting!")
           calc.bias(posterior = pr, prob.range=c(0,1), length = 100, n.binom=Ns, mc.cores=CORES)
        }
    ))

    do.sigmat <- function(pr) sig.matrix(posterior = pr, n.control=Ns, n.treatment = NT,
                                         check.xt=00:NT, check.xs=0:Ns,
                                         level=0.975, mc.cores=CORES )


    SIGMAT <- lapply(posteriors, function(p){
      print("SIGMATING!")
        do.sigmat(p)
    } )



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
    save(n.fix, x.fix,  bias,  cover, t1e, pow, mse, SIGMAT,  file=paste0("Rand-75/study_Rand_",i,".rda"))

  })
}




TF <- rep(FALSE,500)
for(i in 1:500) TF[i] <-  !file.exists(paste0("Rand-75/study_Rand_",i,".rda"))
recalc <- which(TF)

if(TRUE)  mclapply(recalc,calc1, mc.cores=30)


# For new priors  -------------------------------------------------
CORES <- 1
recalc <- 1:500
# recalc <- 501:750
# recalc <- 751:1000

calc <- function(i){
  print(i)
  try({
    load(file=paste0("Rand-75/models_g_",i,".rda"))

    Ns <- NN

    posteriors <-lapply(list(F.CR9,F.SUM),
                          Calc.posterior.all, N=Ns, mc.cores=CORES)

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



    pow <- mclapply(SIGMAT,
                    function(S) calc.power(sig.mat=S, n.binom.control = Ns,
                                           prob.range = c(0,0.85), length =200, treatment.difference=0.12),
                    mc.cores=CORES)

    t1e <- mclapply(SIGMAT,
                    function(S) calc.power(sig.mat=S, n.binom.control = Ns,
                                           prob.range = c(0,0.9), length = 200, treatment.difference = 0),
                    mc.cores=CORES)


    cover <-  c(mclapply(posteriors,
                       function(pr) calc.coverage(posterior=pr, level = 0.95, n.control = Ns, smooth = 0.03),
                       mc.cores=CORES),
                list(calc.coverage(prior = F.ROB, level = 0.95, n.control = Ns, smooth = 0.03)))

    n.fix <- n
    x.fix <- x
    save(n.fix, x.fix,  bias,  cover, t1e, pow, mse, SIGMAT,  file=paste0("Rand-75/study_grand_",i,".rda"))

  })
}



TF2 <- rep(FALSE,500)
for(i in 1:500) TF2[i] <-  !file.exists(paste0("Rand-75/study_grand_",i,".rda"))
recalc2 <- which(TF2)


mclapply(recalc2, calc, mc.cores=30)

}


