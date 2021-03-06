% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun.list.rw.R
\name{posterior.mixture.prior}
\alias{posterior.mixture.prior}
\title{Calculate posterior mixture model}
\usage{
posterior.mixture.prior(xs, ns, sd, mixture.prior)
}
\arguments{
\item{xs}{Number of successes in new trial}

\item{ns}{Number of patients in new trial}

\item{sd}{Standard deviation for normal model}

\item{mixture.prior}{Mixture prior object}
}
\value{
A \code{mixture.prior} object
}
\description{
Calculate posterior mixture model
}
\examples{
\donttest{
xh <- c(30,40,50)
nh <- c(90,95,110)

#fit a full Bayes power prior with 500 samples
fb <- binom.PP.FB(x=xh, n=nh, mixture.size=500, mix=TRUE) 

#calulate the posterior
post.fb <- posterior.mixture.prior(xs=51, ns=100, mixture.prior=fb)

# Compare with plot
plot(fb, xlab="p",ylab="Density", ylim=c(0,12))
#Use lines.only to add to existing plot
plot(post.fb, lines.only=TRUE, col=3)
}

}
