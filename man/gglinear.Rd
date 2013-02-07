\name{gglinear}
\alias{gglinear}
\title{Plot simulated linear time-constant hazards.}
\usage{
  gglinear(obj, qi = "Relative Hazard", from = NULL,
    to = NULL, xlab = NULL, ylab = NULL, title = NULL,
    smoother = "auto", colour = "#A6CEE3",
    spalette = "Set1", leg.name = "", lsize = 2, psize = 1,
    palpha = 0.1, ...)
}
\arguments{
  \item{obj}{a simlinear object}

  \item{qi}{character string indicationg what quantity of
  interest you would like to calculate. Can be
  \code{'Relative Hazard'}, \code{'First Difference'},
  \code{'Hazard Ratio'}, or \code{'Hazard Rate'}. Default
  is \code{qi = 'Relative Hazard'}.}

  \item{xlab}{a label for the plot's x-axis.}

  \item{ylab}{a label of the plot's y-axis. The default
  uses the value of \code{qi}.}
}
\description{
  \code{gglinear} uses ggplot2 to plot the quantities of
  interest from \code{simlinear} objects.
}
\references{
  Licht, Amanda A. 2011. “Change Comes with Time:
  Substantive Interpretation of Nonproportional Hazards in
  Event History Analysis.” Political Analysis 19: 227–43.
}
\seealso{
  \code{\link{coxsimLinear}} and \code{\link{ggplot2}}
}

