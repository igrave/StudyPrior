
normal.PP.EB <- function(x, sd, X, SD){

  r <- 1 / sd^2
  R<- 1 / SD^2
  
  marg <- function(d){
    dr <- d*r
    sumdr <- sum(dr)
    sumdrx <- sum(d*r*x)
    sumdrX <- sum(d*r*X)
    
    exp( - (R*(sumdrx - sumdrX)^2) /
             (2*(sumdr)*(R+sumdr))) *
      sqrt( (sumdr) / (R+sumdr) )
  }
  
  n.hist <- length(x)
  
  opt <- 
  optimx::optimx( rep(0.8, n.hist),
                  marg,
                  control =list(maximize=TRUE),
                  method = "L-BFGS-B",
                  lower = rep(0, n.hist),
                  upper = rep(1, n.hist))

  d <- opt[1,seq(n.hist)]
  
# See http://www.tina-vision.net/docs/memos/2003-003.pdf for product of gaussians

  # eq (4)
all.sd2 <- 1/sum(r,R)
all.mean <- sum()
  
  function(mu){

    dnorm(mu,
          )
    
  1/sqrt(2*pi) *
      dr <- d*r
      sumdr <- sum(dr)
      sumdrx <- sum(d*r*x)
      sumdrX <- sum(d*r*X)

      exp( - (R*(sumdrx - sumdrX)^2) /
             (2*(sumdr)*(R+sumdr))) *
        sqrt( (sumdr) / (R+sumdr) )
    }

} 


x <-c(1,2,50,2.5)
s <- c(1,1,1.4,0.8)
a <- normal.PP.EB(x,s , 2.5, 1.2)

sum(a * x)


