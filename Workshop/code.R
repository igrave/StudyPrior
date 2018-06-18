## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message=FALSE, eval=FALSE)

## ----Practical_1.data----------------------------------------------------

library(StudyPrior)
library(INLA)


x.hist <- c(6,8,17,28,26,8,22,8,6,16,53)
n.hist <- c(33,45,74,103,140,49,83,59,22,109,213)
p.hist <- x.hist/n.hist

knitr::kable(data.frame(Study = 1:11,
                        Failures = paste0(x.hist,'/',n.hist),
                        Percent = paste0(round(p.hist*100),'%')),
             col.names = c("Study","Failures/Controls","Failure %"),
             align='c')


## ----Practical_1.setup---------------------------------------------------
x.c <- 12
n.c <- 80

x.t <- 8
n.t <- 80

## ----Practical_1.q1------------------------------------------------------
my.choice <- 3

x.h <- x.hist[[my.choice]]
n.h <- n.hist[[my.choice]]

## ----Practical_1.q1a-----------------------------------------------------
#approximate with a Riemann sum
p <- seq(0,1,len=2000)
mean(dbeta(p, 1+x.c, 1+n.c-x.c) *
        pbeta(p, 1+x.t, 1+n.t-x.t, lower.tail = FALSE))

# numerical integration
integrate(function(p) dbeta(p, 1+x.c, 1+n.c-x.c) *
            pbeta(p, 1+x.t, 1+n.t-x.t, lower.tail = FALSE),
          lower=0, upper=1)


#Using the package StudyPrior
q1.nohist <- create.mixture.prior("beta",pars=matrix(c(1,1),ncol=2))
q1.nohist.post <- posterior.mixture.prior(x.c, n.c, mix=q1.nohist)

prob.control.smaller(xt=x.t, nt=n.t, posterior=q1.nohist.post)

## ----Practical_1.q1b-----------------------------------------------------
q1.pool <- binom.PP.FIX(x=x.h, n=n.h, d=1, mix = TRUE)
q1.pool.post <- posterior.mixture.prior(x.c, n.c, mix=q1.pool)
prob.control.smaller(xt=x.t,nt=n.t,posterior=q1.pool.post)

## ----Practical_1.q1c-----------------------------------------------------
q1.pp08 <- binom.PP.FIX(x=x.h, n=n.h, d=0.8, mix = TRUE)
q1.pp08.post <- posterior.mixture.prior(x.c, n.c, mix=q1.pp08)

prob.control.smaller(xt=x.t,nt=n.t,posterior=q1.pp08.post)

q1.pp01 <- binom.PP.FIX(x=x.h, n=n.h, d=0.1, mix = TRUE)
q1.pp01.post <- posterior.mixture.prior(x.c, n.c, mix=q1.pp01)

prob.control.smaller(xt=x.t,nt=n.t,posterior=q1.pp01.post)

## ----Practical_1.q1d-----------------------------------------------------
f.a0 <- function(a0,Nh,Rh,N,R,a,b) {
  logf <- {lgamma(a0*Nh+a+b) + lgamma(a0*Rh+R+a) + lgamma(a0*(Nh-Rh)+N-R+b)}-
  {lgamma(a0*Rh+a) + lgamma(a0*(Nh-Rh)+b) + lgamma(a0*Nh+N+a+b)}
  exp(logf)
}

k <- integrate(function(a0) f.a0(a0, n.h, x.h, n.c, x.c, a=1, b=1), lower=0,upper=1)

plot(p, f.a0(p, n.h, x.h, n.c, x.c, a=1, b=1)/k$value, type='l',
     xlab="a0",ylab="posterior density")

mean.a0 <- integrate(function(a0) a0 * f.a0(a0, n.h, x.h, n.c, x.c, a=1, b=1)/k$value,
                     lower=0,upper=1) 
mean.a0$value

q1.pp.mean <- binom.PP.FIX(x.h, n.h, d=mean.a0$value, mix=TRUE)

q1.pp.mean.post <- posterior.mixture.prior(x.c,n.c, mixture.prior = q1.pp.mean)
plot(q1.pp.mean.post)
prob.control.smaller(xt=x.t,nt=n.t,posterior=q1.pp.mean.post)

## ----Practical_1.q1e-----------------------------------------------------
q1.joint <- binom.PP.FB(x.h, n.h,  d.prior.a=1, d.prior.b=1, mix = TRUE)
q1.joint.01 <- binom.PP.FB(x.h, n.h, d.prior.a=0.02, d.prior.b=0.02, mix = TRUE)

q1.joint.post <- posterior.mixture.prior(x.c,n.c, mixture.prior = q1.joint)
q1.joint.01.post <- posterior.mixture.prior(x.c,n.c, mixture.prior = q1.joint.01)

plot(q1.joint.post)
plot(q1.joint.01.post, add=TRUE, col=2)

prob.control.smaller(xt=x.t,nt=n.t,posterior=q1.joint.post)
prob.control.smaller(xt=x.t,nt=n.t,posterior=q1.joint.01.post)

## ----Practical_1.q1g-----------------------------------------------------
q1.pool.robust <- make.robust(q1.pool, 0.5)
plot(q1.pool)
plot(q1.pool.robust, add=TRUE, col=2)

q1.pool.robust.post <- posterior.mixture.prior(x.c,n.t, mix=q1.pool.robust)

plot(q1.pool.post, lty=2, add=TRUE)
plot(q1.pool.robust.post, add=TRUE, col=2, lty=2)

prob.control.smaller(xt=x.t,nt=n.t,posterior=q1.pool.post)
prob.control.smaller(xt=x.t,nt=n.t,posterior=q1.pool.robust.post)

## ----Practical_1.q2------------------------------------------------------
my.choice <- 7

#Repeat as before

## ----Practical_2.q1------------------------------------------------------

q2.MAP <- binom.MAP.FB(x.hist, n.hist)

## ----Practical_2.q2a-----------------------------------------------------
MAP4 <- conj.approx(q2.MAP, type="beta", max.degree = 4)
MAP3 <- conj.approx(q2.MAP, type="beta", max.degree = 3)
MAP2 <- conj.approx(q2.MAP, type="beta", max.degree = 2)
MAP1 <- conj.approx(q2.MAP, type="beta", max.degree = 1)

p <- seq(0,1,len=500)
plot(p, q2.MAP(p), type='l')

plot(MAP1, col=2, add=TRUE)
plot(MAP2, col=3, add=TRUE)
plot(MAP3, col=4, add=TRUE)
plot(MAP4, col=5, add=TRUE)



## ----Practical_2.q2b-----------------------------------------------------

KLD <- function(prior){
  integrate(function(p) q2.MAP(p) *(log(q2.MAP(p)) - log(eval.mixture.prior(p,mixture=prior))),
            lower=0.1, upper=.9)
}

KLD(MAP1)
KLD(MAP2)
KLD(MAP3)
KLD(MAP4)


## ----Practical_2.q3a-----------------------------------------------------
q2.MAP4.post <- posterior.mixture.prior(x.c,n.c, mix=MAP4)

## ----Practical_2.q3b-----------------------------------------------------
q2.MAP.robust <- make.robust(MAP4, 0.1)
q2.MAP.robust.post <- posterior.mixture.prior(x.c,n.c, mix=q2.MAP.robust)

## ----Practical_2.q3c-----------------------------------------------------
q2.fix <- binom.PP.FIX(x.hist, n.hist, d=0.4 , mix = TRUE)
q2.fix.diff <- binom.PP.FIX(x.hist, n.hist, d=rep(c(0.8,0.4), times=c(5,6)), mix = TRUE)

q2.fix.post <- posterior.mixture.prior(x.c, n.c, mix=q2.fix)
q2.fix.diff.post <- posterior.mixture.prior(x.c, n.c, mix=q2.fix.diff)

## ----Practical_2.q3d-----------------------------------------------------
q2.FB <- binom.PP.FB(x.hist, n.hist,
                     d.prior.a = 1, d.prior.b=1,
                     p.prior.a = 1, p.prior.b = 1,
                     mix = TRUE)

q2.FB.post <- posterior.mixture.prior(x.c, n.c, mix=q2.FB)

## ----Practical_2.q3e-----------------------------------------------------
q2.EB <- binom.PP.EB(x.hist, n.hist,X = x.c, n.c, mix = TRUE)
attr(q2.EB,"powers")

q2.EB.post <- posterior.mixture.prior(x.c, n.c, mix=q2.EB)

## ----Practical_2.q3e2----------------------------------------------------
prob.control.smaller(x.t,n.t, posterior=q2.MAP4.post)
prob.control.smaller(x.t,n.t, posterior=q2.MAP.robust.post)

prob.control.smaller(x.t,n.t, posterior=q2.fix.post)
prob.control.smaller(x.t,n.t, posterior=q2.fix.diff.post)

prob.control.smaller(x.t,n.t, posterior=q2.FB.post)
prob.control.smaller(x.t,n.t, posterior=q2.EB.post)


plot(MAP4)
plot(q2.MAP.robust, col=2, add=TRUE)

plot(q2.MAP4.post, add=TRUE, lty=2)
plot(q2.MAP.robust.post, col=2, lty=2, add=TRUE)

plot(q2.FB.post, add=TRUE, col=3, lty=3)
plot(q2.EB.post, add=TRUE, col=4, lty=3)

points(y=rep(0,length(x.hist)), x=x.hist/n.hist,  col="grey")
points(y=0, x=x.c/n.c)


## ------------------------------------------------------------------------
# x.c <- 36
# n.c <- 320
# x.t <- 32
# n.t <- 320
#Continue as before

## ----Practical_3.q1------------------------------------------------------
# q2.MAP
q3.map <- MAP4

## ----Practical_3.q1a-----------------------------------------------------
power <- 0.8
type1e <- 0.05
p1 <- 0.22 #control failure
p2 <-  0.15 # experimental treatment failure

(qnorm(1-power)+qnorm(type1e/2))^2 * ((p1*(1-p1))+p2*(1-p2))/(p1-p2)^2

power.prop.test(n=NULL, p1,p2, type1e, power)


## ------------------------------------------------------------------------


## ----Practical_3.q1c-----------------------------------------------------
N.new <- 482
#Calculate the posterior since we need it for the next questions
q3.map.post <- lapply(0:N.new, posterior.mixture.prior, ns=N.new, mixture.prior=q3.map)

# sig.matrix calculates the test [Pr(Control < Treat)>level] for all values of Xt, Xc
# so we need  1-0.975=level for Treat<Control and then the TRUEs <-> FALSEs
SIG <- !sig.matrix(n.control=N.new, n.treatment=N.new, level=0.025, posterior = q3.map.post)

q3c.power <- calc.power(sig.mat = SIG, treatment.difference = 0.15-0.22, prob.range = c(0.1,0.5),
           length=41, n.binom.control = N.new)

plot(y=q3c.power, x=seq(0.1,0.5,len=41), xlab="Control Probability", ylab="Power", ty='b')
abline(v=0.22, lty=2)

## ----Practical_3.q1d-----------------------------------------------------
# We can reuse the result of sig.matrix() since it doesn't depend on the treatment parameter or null/alternative hypothesis

#we use calc.power with treatment.difference=0 to find the type 1 error
q3c.t1e <- calc.power(sig.mat = SIG, treatment.difference = 0,
                        prob.range = c(0.1,0.5),
                        length=40, n.binom.control = N.new)

plot(y=q3c.t1e, x=seq(0.1,0.5,len=40), xlab="Control Probability", ylab="Type I Error", ty='b')
abline(v=0.22, lty=2)

## ----Practical_3.q2------------------------------------------------------
q3.rob.05 <- make.robust(q3.map, 0.5)
q3.rob.01 <- make.robust(q3.map, 0.1)

q3.rob.05.post <- lapply(0:N.new, posterior.mixture.prior, ns=N.new, mixture.prior=q3.rob.05)
q3.rob.01.post <- lapply(0:N.new, posterior.mixture.prior, ns=N.new, mixture.prior=q3.rob.01)

SIG.05 <- !sig.matrix(n.control=N.new, n.treatment=N.new, level=0.025, posterior = q3.rob.01.post) 


q3c.power.05 <- calc.power(sig.mat = SIG.05, treatment.difference = 0.15-0.22, prob.range = c(0.1,0.5),
           length=41, n.binom.control = N.new)

# and so on

## ----Practical_3B.q1a----------------------------------------------------
n.ci <- n.ti <- 241
x.ci <- rbinom(1,241,p1)
x.ti <- rbinom(1,241,p2)

## ----Practical_3B.q1b----------------------------------------------------
f.a0 <- function(a0,Nh,Rh,N,R,a,b) {
  logf <- {lgamma(a0*Nh+a+b) + lgamma(a0*Rh+R+a) + lgamma(a0*(Nh-Rh)+N-R+b)}-
  {lgamma(a0*Rh+a) + lgamma(a0*(Nh-Rh)+b) + lgamma(a0*Nh+N+a+b)}
  exp(logf)
}

k <- integrate(function(a0) f.a0(a0, n.h, x.h, n.ci, x.ci, a=1, b=1), lower=0,upper=1)

plot(p, f.a0(p, n.h, x.h, n.ci, x.ci, a=1, b=1)/k$value, type='l',
     xlab="a0",ylab="posterior density")

mean.a0 <- integrate(function(a0) a0 * f.a0(a0, n.h, x.h, n.ci, x.ci, a=1, b=1)/k$value,
                     lower=0,upper=1) 

mean.a0$value

q3b.pp.interim <- posterior.mixture.prior(x.ci, n.ci,
                                          mixture=binom.PP.FIX(x.h, n.h, d = mean.a0$value, mix = TRUE))



## ----Practical_3B.q1c----------------------------------------------------

n.c2 <- 481 - 282
x.c2 <- 0:n.c2
x.t <- 0:482
n.t <- 0482

q3c.post <- lapply(x.c2, posterior.mixture.prior, ns=n.c2, mixture.prior=q3b.pp.interim)

#now we calculate the significance, but only for the post-interim values for control, since the interim values are included in the prior
SIG.q3c <- !sig.matrix(posterior = q3c.post, 
                      check.xt = x.t, check.xs = x.c2,
                      n.control=n.c2, n.treatment=n.t,
                      level = 0.0275)


q3c.power <- calc.power(sig.mat = SIG.q3c, treatment.difference = 0.15-0.22,
                        prob.range = c(0.1,0.5), length=41, 
                        n.binom.control = n.c2, n.binom.treatment = 482)

## ----Practical_3B.q1d----------------------------------------------------

n.c2 <- 481 - 282
x.c2 <- 0:n.c2
x.t <- 0:482
n.t <- 0482


#Do the same but set the power to 0 to ignroe the historical data
q3d.pp.interim <- posterior.mixture.prior(x.ci, n.ci,
                                          mixture=binom.PP.FIX(x.h, n.h, d = 0, mix = TRUE))
q3d.post <- lapply(x.c2, posterior.mixture.prior, ns=n.c2, mixture.prior=q3d.pp.interim)

#now we calculate the significance, but only for the post-interim values for control, since the interim values are included in the prior
SIG.q3d <- !sig.matrix(posterior = q3d.post, 
                      check.xt = x.t, check.xs = x.c2,
                      n.control=n.c2, n.treatment=n.t,
                      level = 0.0275)


q3d.power <- calc.power(sig.mat = SIG.q3d, treatment.difference = 0.15-0.22,
                        prob.range = c(0.1,0.5), length=41, 
                        n.binom.control = n.c2, n.binom.treatment = 482)

