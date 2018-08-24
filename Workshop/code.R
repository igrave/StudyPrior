## ----setup, include=FALSE------------------------------------------------
## knitr::opts_chunk$set(echo = TRUE, message=FALSE, eval=FALSE)

## ----Practical_1.data----------------------------------------------------
## 
## library(StudyPrior)
## library(INLA)
## 
## 
## x.hist <- c(6,8,17,28,26,8,22,8,6,16,53)
## n.hist <- c(33,45,74,103,140,49,83,59,22,109,213)
## p.hist <- x.hist/n.hist
## 
## knitr::kable(data.frame(Study = 1:11,
##                         Failures = paste0(x.hist,'/',n.hist),
##                         Percent = paste0(round(p.hist*100),'%')),
##              col.names = c("Study","Failures/Controls","Failure %"),
##              align='c')
## 

## ----Practical_1.setup---------------------------------------------------
## x.c <- 12
## n.c <- 80
## 
## x.t <- 8
## n.t <- 80

## ----Practical_1.q1------------------------------------------------------
## my.choice <- 3    # Enter a study number between 1 and 11
## 
## x.h <- x.hist[[my.choice]]
## n.h <- n.hist[[my.choice]]

## ----Practical_1.q1a-----------------------------------------------------
## #approximate with a Riemann sum
## p <- seq(0,1,len=2000)
## mean(dbeta(p, 1+x.c, 1+n.c-x.c) *
##         pbeta(p, 1+x.t, 1+n.t-x.t, lower.tail = TRUE))
## 
## # numerical integration
## integrate(function(p) dbeta(p, 1+x.c, 1+n.c-x.c) *
##             pbeta(p, 1+x.t, 1+n.t-x.t, lower.tail = TRUE),
##           lower=0, upper=1)
## 
## 
## #Using the package StudyPrior
## q1.nohist <- create.mixture.prior("beta",pars=matrix(c(1,1),ncol=2))
## q1.nohist.post <- posterior.mixture.prior(x.c, n.c, mix=q1.nohist)
## 
## 1-prob.control.smaller(xt=x.t, nt=n.t, posterior=q1.nohist.post)
## 
## #Plot posteriors for control and active rate
## plot(q1.nohist.post)
## q1.active <- create.mixture.prior("beta",pars=matrix(c(1,1),ncol=2))
## q1.active.post <- posterior.mixture.prior(x.t, n.t, mix=q1.active)
## plot(q1.active.post, col=2, add=TRUE)

## ----Practical_1.q1b-----------------------------------------------------
## q1.pool <- binom.PP.FIX(x=x.h, n=n.h, d=1, mix = TRUE)
## q1.pool.post <- posterior.mixture.prior(x.c, n.c, mix=q1.pool)
## 1-prob.control.smaller(xt=x.t,nt=n.t,posterior=q1.pool.post)
## 
## plot(q1.pool.post)

## ----Practical_1.q1c-----------------------------------------------------
## q1.pp08 <- binom.PP.FIX(x=x.h, n=n.h, d=0.8, mix = TRUE)
## q1.pp08.post <- posterior.mixture.prior(x.c, n.c, mix=q1.pp08)
## 
## 1-prob.control.smaller(xt=x.t,nt=n.t,posterior=q1.pp08.post)
## 
## q1.pp01 <- binom.PP.FIX(x=x.h, n=n.h, d=0.1, mix = TRUE)
## q1.pp01.post <- posterior.mixture.prior(x.c, n.c, mix=q1.pp01)
## 
## 1-prob.control.smaller(xt=x.t,nt=n.t,posterior=q1.pp01.post)
## 
## plot(q1.pp01.post)

## ----Practical_1.q1d-----------------------------------------------------
## f.a0 <- function(a0,Nh,Rh,N,R,a,b) {
##   logf <- {lgamma(a0*Nh+a+b) + lgamma(a0*Rh+R+a) + lgamma(a0*(Nh-Rh)+N-R+b)}-
##   {lgamma(a0*Rh+a) + lgamma(a0*(Nh-Rh)+b) + lgamma(a0*Nh+N+a+b)}
##   exp(logf)
## }
## 
## k <- integrate(function(a0) f.a0(a0, n.h, x.h, n.c, x.c, a=1, b=1), lower=0,upper=1)
## 
## p <- seq(0,1,len=2000)
## plot(p, f.a0(p, n.h, x.h, n.c, x.c, a=1, b=1)/k$value, type='l',
##      xlab="a0",ylab="posterior density")
## 
## mean.a0 <- integrate(function(a0) a0 * f.a0(a0, n.h, x.h, n.c, x.c, a=1, b=1)/k$value,
##                      lower=0,upper=1)
## mean.a0$value
## 
## q1.pp.mean <- binom.PP.FIX(x.h, n.h, d=mean.a0$value, mix=TRUE)
## 
## q1.pp.mean.post <- posterior.mixture.prior(x.c,n.c, mixture.prior = q1.pp.mean)
## plot(q1.pp.mean.post)
## 1-prob.control.smaller(xt=x.t,nt=n.t,posterior=q1.pp.mean.post)

## ----Practical_1.q1e-----------------------------------------------------
## q1.joint <- binom.PP.FB(x.h, n.h,  d.prior.a=1, d.prior.b=1, mix = TRUE)
## q1.joint.01 <- binom.PP.FB(x.h, n.h, d.prior.a=0.02, d.prior.b=0.02, mix = TRUE)
## 
## q1.joint.post <- posterior.mixture.prior(x.c,n.c, mixture.prior = q1.joint)
## q1.joint.01.post <- posterior.mixture.prior(x.c,n.c, mixture.prior = q1.joint.01)
## 
## plot(q1.joint.post)
## plot(q1.joint.01.post, add=TRUE, col=2)
## 
## 1-prob.control.smaller(xt=x.t,nt=n.t,posterior=q1.joint.post)
## 1-prob.control.smaller(xt=x.t,nt=n.t,posterior=q1.joint.01.post)

## ----Practical_1.q2------------------------------------------------------
## plot(q1.nohist.post, ylim=c(0,14), xlim=c(0,0.4))
## plot(q1.pool.post, add=TRUE, col=2)
## plot(q1.pp08.post, add=TRUE, col=3)
## plot(q1.pp01.post, add=TRUE, col=4)
## plot(q1.pp.mean.post, add=TRUE, col=5)
## plot(q1.joint.post, add=TRUE, col=6)
## plot(q1.joint.01.post, add=TRUE, col=7)
## legend("topright", legend=c("No hist", "Pooled", "PP FIX 0.8", "PP FIX 0.1", "PP MEAN", "PP FULL 1", "PP FULL 0.2"), lty=rep(1,7), col=1:7)
## 

## ----Practical_1.q3------------------------------------------------------
## my.choice <- 7
## 
## #Repeat as before

## ----Practical_1.q1f-----------------------------------------------------
## q1.pool.robust <- make.robust(q1.pool, 0.5)
## plot(q1.pool)
## plot(q1.pool.robust, add=TRUE, col=2)
## 
## q1.pool.robust.post <- posterior.mixture.prior(x.c,n.t, mix=q1.pool.robust)
## 
## plot(q1.pool.post, lty=2, add=TRUE)
## plot(q1.pool.robust.post, add=TRUE, col=2, lty=2)
## 
## 1-prob.control.smaller(xt=x.t,nt=n.t,posterior=q1.pool.post)
## 1-prob.control.smaller(xt=x.t,nt=n.t,posterior=q1.pool.robust.post)

## ----Practical_1.q1g-----------------------------------------------------
## q1.EB <- binom.PP.EB(x.h, n.h,X = x.c, n.c, mix = TRUE)
## attr(q1.EB,"powers")
## 
## q1.EB.post <- posterior.mixture.prior(x.c, n.c, mix=q1.EB)
## 1-prob.control.smaller(xt=x.t,nt=n.t,posterior=q1.EB.post)

## ----Practical_2.q1------------------------------------------------------
## 
## q2.MAP <- binom.MAP.FB(x.hist, n.hist)

## ----Practical_2.q2a-----------------------------------------------------
## MAP4 <- conj.approx(q2.MAP, type="beta", max.degree = 4)
## MAP3 <- conj.approx(q2.MAP, type="beta", max.degree = 3)
## MAP2 <- conj.approx(q2.MAP, type="beta", max.degree = 2)
## MAP1 <- conj.approx(q2.MAP, type="beta", max.degree = 1)
## 
## p <- seq(0,1,len=500)
## plot(p, q2.MAP(p), type='l')
## 
## plot(MAP1, col=2, add=TRUE)
## plot(MAP2, col=3, add=TRUE)
## plot(MAP3, col=4, add=TRUE)
## plot(MAP4, col=5, add=TRUE)
## 
## 

## ----Practical_2.q2b-----------------------------------------------------
## 
## KLD <- function(prior){
##   integrate(function(p) q2.MAP(p) *(log(q2.MAP(p)) - log(eval.mixture.prior(p,mixture=prior))),
##             lower=0.1, upper=.9)
## }
## 
## KLD(MAP1)
## KLD(MAP2)
## KLD(MAP3)
## KLD(MAP4)
## 

## ----Practical_2.q3a-----------------------------------------------------
## q2.MAP4.post <- posterior.mixture.prior(x.c,n.c, mix=MAP4)

## ----Practical_2.q3b-----------------------------------------------------
## q2.MAP.robust.01 <- make.robust(MAP4, 0.1)
## q2.MAP.robust.01.post <- posterior.mixture.prior(x.c,n.c, mix=q2.MAP.robust.01)
## 
## q2.MAP.robust.05 <- make.robust(MAP4, 0.5)
## q2.MAP.robust.05.post <- posterior.mixture.prior(x.c,n.c, mix=q2.MAP.robust.05)

## ----Practical_2.q3c-----------------------------------------------------
## q2.fix <- binom.PP.FIX(x.hist, n.hist, d=0.4 , mix = TRUE)
## q2.fix.diff <- binom.PP.FIX(x.hist, n.hist, d=rep(c(0.8,0.4), times=c(5,6)), mix = TRUE)
## 
## q2.fix.post <- posterior.mixture.prior(x.c, n.c, mix=q2.fix)
## q2.fix.diff.post <- posterior.mixture.prior(x.c, n.c, mix=q2.fix.diff)

## ----Practical_2.q3d-----------------------------------------------------
## q2.FB <- binom.PP.FB(x.hist, n.hist,
##                      d.prior.a = 1, d.prior.b=1,
##                      p.prior.a = 1, p.prior.b = 1,
##                      mix = TRUE)
## 
## q2.FB.post <- posterior.mixture.prior(x.c, n.c, mix=q2.FB)

## ----Practical_2.q3e-----------------------------------------------------
## q2.EB <- binom.PP.EB(x.hist, n.hist, x.c, n.c, mix = TRUE)
## attr(q2.EB,"powers")
## 
## q2.EB.post <- posterior.mixture.prior(x.c, n.c, mix=q2.EB)

## ----Practical_2.q3prob--------------------------------------------------
## 1-prob.control.smaller(x.t,n.t, posterior=q2.MAP4.post)
## 1-prob.control.smaller(x.t,n.t, posterior=q2.MAP.robust.01.post)
## 1-prob.control.smaller(x.t,n.t, posterior=q2.MAP.robust.05.post)
## 
## 1-prob.control.smaller(x.t,n.t, posterior=q2.fix.post)
## 1-prob.control.smaller(x.t,n.t, posterior=q2.fix.diff.post)
## 
## 1-prob.control.smaller(x.t,n.t, posterior=q2.FB.post)
## 1-prob.control.smaller(x.t,n.t, posterior=q2.EB.post)

## ----Practical_2.q3plot--------------------------------------------------
## plot(q2.MAP4.post, ylim=c(0,25), xlim=c(0,0.5), xlab=expression(theta))
## plot(q2.MAP.robust.01.post, col=2, add=TRUE)
## plot(q2.MAP.robust.05.post, col=3, add=TRUE)
## 
## plot(q2.fix.post, add=TRUE, col=4)
## plot(q2.fix.diff.post, add=TRUE, col=5)
## 
## plot(q2.FB.post, add=TRUE, col=6)
## plot(q2.EB.post, add=TRUE, col=7)
## 
## legend("topright", legend=c("MAP", "Robust MAP 0.1", "Robust MAP 0.5", "Fixed PP (1)",
##                             "Fixed PP (2)", "Joint PP", "Empirical PP",
##                             "Historical control rates", "Current control rate"), lty=c(rep(1,7), NA, NA), col=c(1:7, 1,2), pch=c(rep(NA, 7), 16,16))
## points(y=rep(0,length(x.hist)), x=x.hist/n.hist,  pch=16)
## points(y=0, x=x.c/n.c, pch=16, col=2)
## 

## ------------------------------------------------------------------------
## # x.c <- 36
## # n.c <- 320
## # x.t <- 32
## # n.t <- 320
## #Continue as before

## ----Practical_3.q1------------------------------------------------------
## # q2.MAP
## q3.map <- MAP4

## ----Practical_3.q1a-----------------------------------------------------
## power <- 0.8
## type1e <- 0.05
## p1 <- 0.22 #control failure
## p2 <-  0.15 # experimental treatment failure
## 
## (qnorm(1-power)+qnorm(type1e/2))^2 * ((p1*(1-p1))+p2*(1-p2))/(p1-p2)^2
## 
## power.prop.test(n=NULL, p1,p2, type1e, power)
## 

## ----Practical_3.q1b-----------------------------------------------------
## q3.map.ess <- ess.morita.mixture.prior(q3.map,LEN = 500) #takes about 2 minutes
## 
## q3.map.ess

## ----Practical_3.q1c-----------------------------------------------------
## N.new <- 482
## #Calculate the posterior since we need it for the next questions
## q3.map.post <- lapply(0:N.new, posterior.mixture.prior, ns=N.new, mixture.prior=q3.map)
## 
## # sig.matrix calculates the test [Pr(Control < Treat)>level] for all values of Xt, Xc
## # so we need  1-0.975=level for Treat<Control and then the TRUEs <-> FALSEs
## SIG <- !sig.matrix(n.control=N.new, n.treatment=N.new, level=0.025, posterior = q3.map.post)
## 
## # Calculate power to detect a treatment difference of -0.05
## q3c.power <- calc.power(sig.mat = SIG, treatment.difference = 0.15-0.22, prob.range = c(0.1,0.5), length=41, n.binom.control = N.new)

## ----Practical_3.q1c_plot------------------------------------------------
## # Plot power as a function of assumed control failure rate for new trial
## plot(y=q3c.power, x=seq(0.1,0.5,len=41), xlab="Control Failure Rate (new trial)", ylab="Power", ty='b')
## abline(v=0.22, lty=2)  # mark value of control rate assumed for new trial
## x=seq(0.1,0.5,len=41)
## q3c.power[x==0.22]   # Power of study using MAP prior plus 482 current control subjects
## 
## # Calculate power for standard design with no historical data
## q3.nohist.post <- lapply(0:N.new, posterior.mixture.prior, ns=N.new, mixture.prior=q1.nohist)
## 
## SIG.nohist <- !sig.matrix(n.control=N.new, n.treatment=N.new, level=0.025, posterior = q3.nohist.post)
## 
## q3c.power.nohist <- calc.power(sig.mat = SIG.nohist, treatment.difference = 0.15-0.22, prob.range = c(0.1,0.5), length=41, n.binom.control = N.new)
## 
## # Add power curve to plot
## lines(y=q3c.power.nohist, x=seq(0.1,0.5,len=41), col=2)

## ----Practical_3.q1d-----------------------------------------------------
## # We can reuse the result of sig.matrix() since it doesn't depend on the treatment parameter or null/alternative hypothesis
## 
## #we use calc.power with treatment.difference=0 to find the type 1 error
## q3c.t1e <- calc.power(sig.mat = SIG, treatment.difference = 0,
##                         prob.range = c(0.1,0.5),
##                         length=41, n.binom.control = N.new)

## ----Practical_3.q1d_plot------------------------------------------------
## plot(y=q3c.t1e, x=seq(0.1,0.5,len=41), xlab="Control Failure Rate (new trial)", ylab="Type I Error", ty='b', ylim=c(0, 0.1))
## abline(v=0.22, lty=2)  # mark value of control rate assumed for new trial
## x=seq(0.1,0.5,len=41)
## q3c.t1e[x==0.22]   # Type 1 error of study using MAP prior plus 482 current control subjects
## 
## #Type 1 error for standard design without historical data
## q3c.t1e.nohist <- calc.power(sig.mat = SIG.nohist, treatment.difference = 0,
##                         prob.range = c(0.1,0.5),
##                         length=41, n.binom.control = N.new)
## 
## #Add to plot
## lines(y=q3c.t1e.nohist, x=seq(0.1,0.5,len=41), col=2)
## 

## ----Practical_3.q2------------------------------------------------------
## q3.rob.05 <- make.robust(q3.map, 0.5)
## q3.rob.01 <- make.robust(q3.map, 0.1)
## q3.EB <- lapply(0:N.new, function(X) binom.PP.EB(x.hist, n.hist, X, N.new, mix = TRUE))
## q3.FB <- q2.FB
## q3.fix <- q2.fix
## 
## q3.rob.05.post <- lapply(0:N.new, posterior.mixture.prior, ns=N.new, mixture.prior=q3.rob.05)
## q3.rob.01.post <- lapply(0:N.new, posterior.mixture.prior, ns=N.new, mixture.prior=q3.rob.01)
## q3.fix.post <- lapply(0:N.new, posterior.mixture.prior, ns=N.new, mixture.prior=q3.fix)
## q3.FB.post <- lapply(0:N.new, posterior.mixture.prior, ns=N.new, mixture.prior=q3.FB)
## q3.EB.post <- lapply(0:N.new, function(i) {posterior.mixture.prior( xs=i, ns=N.new, mixture.prior=q3.EB[[i+1]])})
## 
## SIG.05 <- !sig.matrix(n.control=N.new, n.treatment=N.new, level=0.025, posterior = q3.rob.05.post)
## 
## q3c.power.05 <- calc.power(sig.mat = SIG.05, treatment.difference = 0.15-0.22, prob.range = c(0.1,0.5), length=41, n.binom.control = N.new)
## q3c.t1e.05 <- calc.power(sig.mat = SIG.05, treatment.difference = 0, prob.range = c(0.1,0.5), length=41, n.binom.control = N.new)
## 
## 
## SIG.01 <- !sig.matrix(n.control=N.new, n.treatment=N.new, level=0.025, posterior = q3.rob.01.post)
## 
## q3c.power.01 <- calc.power(sig.mat = SIG.01, treatment.difference = 0.15-0.22, prob.range = c(0.1,0.5), length=41, n.binom.control = N.new)
## q3c.t1e.01 <- calc.power(sig.mat = SIG.01, treatment.difference = 0, prob.range = c(0.1,0.5), length=41, n.binom.control = N.new)
## 
## SIG.fix <- !sig.matrix(n.control=N.new, n.treatment=N.new, level=0.025, posterior = q3.fix.post)
## 
## q3c.power.fix <- calc.power(sig.mat = SIG.fix, treatment.difference = 0.15-0.22, prob.range = c(0.1,0.5), length=41, n.binom.control = N.new)
## q3c.t1e.fix <- calc.power(sig.mat = SIG.fix, treatment.difference = 0, prob.range = c(0.1,0.5), length=41, n.binom.control = N.new)
## 
## SIG.FB <- !sig.matrix(n.control=N.new, n.treatment=N.new, level=0.025, posterior = q3.FB.post)
## 
## q3c.power.FB <- calc.power(sig.mat = SIG.FB, treatment.difference = 0.15-0.22, prob.range = c(0.1,0.5), length=41, n.binom.control = N.new)
## q3c.t1e.FB <- calc.power(sig.mat = SIG.FB, treatment.difference = 0, prob.range = c(0.1,0.5), length=41, n.binom.control = N.new)
## 
## 
## SIG.EB <- !sig.matrix(n.control=N.new, n.treatment=N.new, level=0.025, posterior = q3.EB.post)
## 
## q3c.power.EB <- calc.power(sig.mat = SIG.EB, treatment.difference = 0.15-0.22, prob.range = c(0.1,0.5), length=41, n.binom.control = N.new)
## q3c.t1e.EB <- calc.power(sig.mat = SIG.EB, treatment.difference = 0, prob.range = c(0.1,0.5), length=41, n.binom.control = N.new)
## 
## 
## # Plot power as a function of assumed control failure rate for new trial, and overlay power curves for each prior
## plot(y=q3c.power, x=seq(0.1,0.5,len=41), xlab="Control Failure Rate (new trial)", ylab="Power", ty='l', ylim=c(0,1))
## lines(y=q3c.power.05, x=seq(0.1,0.5,len=41), col=2)
## lines(y=q3c.power.01, x=seq(0.1,0.5,len=41), col=3)
## lines(y=q3c.power.fix, x=seq(0.1,0.5,len=41), col=4)
## lines(y=q3c.power.FB, x=seq(0.1,0.5,len=41), col=5)
## lines(y=q3c.power.EB, x=seq(0.1,0.5,len=41), col=6)
## lines(y=q3c.power.nohist, x=seq(0.1,0.5,len=41), col=7)
## abline(v=0.22, lty=2)  # mark value of control rate assumed for new trial
## legend("bottomleft", legend=c("MAP", "Robust MAP 0.5", "Robust MAP 0.1", "Fixed PP", "FB PP", "EB PP", "No historical"), lty=rep(1,7), col=1:7, cex=0.6)
## 
## # Plot type 1 error as a function of assumed control failure rate for new trial, and overlay type 1 error curves for each prior
## plot(y=q3c.t1e, x=seq(0.1,0.5,len=41), xlab="Control Failure Rate (new trial)", ylab="Type 1 Error", ty='l', ylim=c(0,0.1))
## lines(y=q3c.t1e.05, x=seq(0.1,0.5,len=41), col=2)
## lines(y=q3c.t1e.01, x=seq(0.1,0.5,len=41), col=3)
## lines(y=q3c.t1e.fix, x=seq(0.1,0.5,len=41), col=4)
## lines(y=q3c.t1e.FB, x=seq(0.1,0.5,len=41), col=5)
## lines(y=q3c.t1e.EB, x=seq(0.1,0.5,len=41), col=6)
## lines(y=q3c.t1e.nohist, x=seq(0.1,0.5,len=41), col=7)
## abline(v=0.22, lty=2)  # mark value of control rate assumed for new trial
## legend("topright", legend=c("MAP", "Robust MAP 0.5", "Robust MAP 0.1", "Fixed PP", "FB PP", "EB PP", "No historical"), lty=rep(1,7), col=1:7, cex=0.6)
## 

