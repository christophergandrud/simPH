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
}
\description{
  \code{SpinBounds} internal function for finding the lower
  bound from SPIn
}
\references{
  Liu, Ying, Andrew Gelman, and Tian Zheng. 2013.
  “Simulation-Efficient Shortest Probablility Intervals.”
  Arvix. http://arxiv.org/pdf/1302.2142v1.pdf.
}

