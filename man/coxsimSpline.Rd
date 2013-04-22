\name{coxsimSpline}
\alias{coxsimSpline}
\title{Simulated quantities of interest for penalised splines from \code{coxph} models.}
\usage{
  coxsimSpline(obj, bspline, bdata, qi = "Relative Hazard",
    Xj = 1, Xl = 0, nsim = 1000, ci = 0.95, spin = FALSE)
}
\arguments{
  \item{obj}{a \code{coxph} fitted model object with a
  penalised spline.}

  \item{bspline}{a character string of the full
  \code{\link{psline}} call used in \code{obj}.}

  \item{bdata}{a numeric vector of splined variable's
  values.}

  \item{qi}{quantity of interest to simulate. Values can be
  \code{"Relative Hazard"}, \code{"First Difference"},
  \code{"Hazard Ratio"}, and \code{"Hazard Rate"}. The
  default is \code{qi = "Relative Hazard"}. Think carefully
  before using \code{qi = "Hazard Rate"}. You may be
  creating very many simulated values which can be very
  computationally intensive to do. Adjust the number of
  simulations per fitted value with \code{nsim}.}

  \item{Xj}{numeric vector of values of X to simulate for.}

  \item{Xl}{numeric vector of values to compare \code{Xj}
  to. Note if \code{qi = "Relative Hazard"} or \code{code =
  "Hazard"} only \code{Xj} is relevant.}

  \item{nsim}{the number of simulations to run per value of
  X. Default is \code{nsim = 1000}.}

  \item{ci}{the proportion of middle simulations to keep.
  The default is \code{ci = 0.95}, i.e. keep the middle 95
  percent. If \code{spin = TRUE} then \code{ci} is the
  convidence level of the shortest probability interval.
  Any value from 0 through 1 may be used.}

  \item{spin}{logical, whether or not to keep only the
  shortest proability interval rather than the middle
  simulations.}
}
\value{
  a simspline object
}
\description{
  Simulates relative hazards, first differences, hazard
  ratios, and hazard rates for penalised splines from Cox
  Proportional Hazards models. These can be plotted with
  \code{\link{ggspline}}.

  Currently does not support simulating hazard rates form
  multiple stratified models.
}
\details{
  \code{coxsimSpline} simulates quantities of interest from
  penalised splines using multivariate normal
  distributions.
}
\examples{
# Load Carpenter (2002) data
data("CarpenterFdaData")

# Load survival package
library(survival)

# Run basic model
# From Keele (2010) replication data
M1 <- coxph(Surv(acttime, censor) ~  prevgenx + lethal + deathrt1 + acutediz + hosp01  + pspline(hospdisc, df = 4) + pspline(hhosleng, df = 4) + mandiz01 + femdiz01 + peddiz01 + orphdum + natreg + vandavg3 + wpnoavg3 + pspline(condavg3, df = 4) + pspline(orderent, df = 4) + pspline(stafcder, df = 4), data = CarpenterFdaData)

# Simulate Relative Hazards for orderent
Sim1 <- coxsimSpline(M1, bspline = "pspline(orderent, df = 4)",
                       bdata = CarpenterFdaData$orderent,
                       qi = "Relative Hazard",
                       Xj = seq(2, 30, by = 3))

# Simulate Hazard Rates for orderent
Sim2 <- coxsimSpline(M1, bspline = "pspline(orderent, df = 4)",
                       bdata = CarpenterFdaData$orderent,
                       qi = "Hazard Rate",
                       Xj = seq(2, 53, by = 3),
                       nsim = 100)
}
\references{
  Luke Keele, "Replication data for: Proportionally
  Difficult: Testing for Nonproportional Hazards In Cox
  Models", 2010, http://hdl.handle.net/1902.1/17068 V1
  [Version]

  King, Gary, Michael Tomz, and Jason Wittenberg. 2000.
  “Making the Most of Statistical Analyses: Improving
  Interpretation and Presentation.” American Journal of
  Political Science 44(2): 347–61.

  Liu, Ying, Andrew Gelman, and Tian Zheng. 2013.
  “Simulation-Efficient Shortest Probablility Intervals.”
  Arvix. http://arxiv.org/pdf/1302.2142v1.pdf.
}
\seealso{
  \code{\link{gg}}, \code{\link{pspline}},
  \code{\link{survival}}, \code{\link{strata}}, and
  \code{\link{coxph}}
}

