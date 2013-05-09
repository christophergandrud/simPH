\name{SpinBounds}
\alias{SpinBounds}
\title{Internal function for finding the lower and upper bound of the shortest probability interval.}
\usage{
  SpinBounds(x, conf = 0.95, LowUp = NULL, bw = 0, lb = 0,
    ub = Inf, l = NA, u = NA)
}
\arguments{
  \item{x}{the name of the simulated quantities of interest
  variable}

  \item{conf}{numeric confidence level}

  \item{LowUp}{numeric specifying if you want to find the
  lower or upper bound of the interval. Possible options
  include \code{1} (Low) or \code{2} (High).}

  \item{bw}{scalar, the bandwidth of the weighting kernel
  in terms of sample points. If not specified, sqrt(n) will
  be used, where n is the sample size.}

  \item{lb,ub}{scalars, the lower and upper bounds of the
  distribution. If specified, a pseudo-sample point equal
  to the corresponding bound will be added. By default,
  \code{lb = 0} as this is the relevant lower bound for
  quantities of interest from survival models.
  \code{\link{coxsimInteract}} modifies \code{lb} to
  \code{-Inf} if \code{qi = "Marginal Effect"}.}

  \item{l,u}{scalars, weighting centers (if provided).}
}
\description{
  \code{SpinBounds} is an internal function used by
  \code{simPH}'s simulation commands to find the shortest
  probability interval. It is largely drawn from Liu's
  (2013) \code{SPIn} command (version 1.1), with two
  modifications. First it returns just the lower or upper
  bound of the interval, rather than the whole SPIn object.
  Second, if there is no variation in the interval (as the
  hazard ratio for a fitted value of 0, e.g. all simulation
  values are 1) it takes the "central interval" of the
  simulations. Effectively it reduces the number of
  simulations by \code{nsim} * \code{ci}.
}
\details{
  \code{SpinBounds} internal function for finding the lower
  bound from SPIn. This function and documentation is
  largely based directly on the original \code{SPIn}
  command by Ying Liu.
}
\references{
  Liu, Ying, Andrew Gelman, and Tian Zheng. 2013.
  ''Simulation-Efficient Shortest Probablility Intervals.''
  Arvix. http://arxiv.org/pdf/1302.2142v1.pdf.
}

