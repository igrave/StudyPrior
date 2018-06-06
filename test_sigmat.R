



fb <- binom.PP.FB.COR(x=c(64,82,42),n=c(100,100,100),mixture.size = 1000, mix=TRUE)

fp <- eval.mixture.prior(p, fb)

plow <- p.mixture.prior(p, fb)


fp <- dbeta(p, 1+0, 1+10)
phigh <- pbeta(p, 1+10, 1+0, lower.tail = FALSE)

sum(fp*phigh/10000)


StudyPrior:::prXY(a = 1,11,11,1, approx.allowed = FALSE)
StudyPrior:::prXY(a = 60,51,31,51, force = TRUE)



POST<- lapply(0:Nc, posterior.mixture.prior, ns=Nc, mixture.prior=fb)


v1.time <- system.time(v1 <- sig.matrix(Nc,Nt, posterior = POST, approx = TRUE))
v2.time <- system.time(v2 <- sigmat3(Nc,Nt, posterior = POST))
v3.time <- system.time(v3 <- sigmat3.c(Nc,Nt, posterior = POST))
v3.time <- system.time(v3 <- sigmat3.c(Nc,Nt, posterior = POST,lenp = 10000))



Nc <- 50
Nt <- 200
system.time(a <- outer(0:Nc,0:Nt, Vectorize(function(x,y) StudyPrior:::prXY(1+x,1+Nc-x,1+y,1+Nt-y, approx=FALSE))))
system.time(a2 <- outer(0:Nc,0:Nt, Vectorize(function(x,y) StudyPrior:::prXY(1+x,1+Nc-x,1+y,1+Nt-y))))

system.time(b <- outer(0:Nc,0:Nt, Vectorize(function(x,y) sum(dbeta(p, 1+x, 1+Nc-x) * pbeta(p, 1+y, 1+Nt-y, lower=FALSE)/1000))))


sigmat3 <- function(Nc,Nt, posterior, lenp=1000){
  
  p <-  seq(0,1,len=lenp)
  MAT1 <- vapply(POST, FUN.VALUE = numeric(lenp), function(m) eval.mixture.prior(p,m))
  MAT2 <- vapply(0:Nt, FUN.VALUE = numeric(lenp), function(x) pbeta(p, 1+x, 1+Nt-x, lower=FALSE))
  
  return( t(MAT1) %*% (MAT2) > .975*lenp  )
}


sigmat3.c <- compiler::cmpfun(sigmat3)
