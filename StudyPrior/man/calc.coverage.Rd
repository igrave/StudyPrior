% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Coverage.R
\name{calc.coverage}
\alias{calc.coverage}
\title{Coverage of Posterior using local smoothing with beta kernels}
\usage{
calc.coverage(prior, level, n.control, smooth, posterior = NULL, ...)
}
\arguments{
\item{prior}{Prior to calculate posterior. Specify posterior instead if available.}

\item{level}{Significance level}

\item{n.control}{Number of patients in new trial}

\item{smooth}{Local smoothing parameter for beta kernels}

\item{posterior}{Posterior density}

\item{...}{Additional arguments for plot}
}
\value{
A vector of coverage values
}
\description{
Coverage of Posterior using local smoothing with beta kernels
}
\examples{
\donttest{
xh <- c(30,40,50)
nh <- c(90,95,110)
fix <- binom.PP.FIX(x=xh, n=nh, d=c(.2,.5,.7), mix=TRUE) # set different weights
calc.coverage(prior = fix, level=0.95, n.control=75, smooth=0.05)
}

}
