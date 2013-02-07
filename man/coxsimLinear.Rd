\name{coxsimLinear}
\alias{coxsimLinear}
\title{Simulate hazards for linear time-constant covariates from Cox Proportional Hazards models.}
\usage{
  coxsimLinear(obj, b, Xj, Xl, nsim = 1000, ci = "95")
}
\arguments{
  \item{b}{character string name of the coefficient you
  would like to simulate.}

  \item{X}{numeric vector of values of X to simulate
  relative hazards for.}

  \item{nsim}{the number of simulations to run per value of
  X. Default is \code{nsim = 1000}.}

  \item{ci}{the proportion of middle simulations to keep.
  The default is \code{ci = "95"}, i.e. keep the middle 95
  percent. Other possibilities include: \code{"90"},
  \code{"99"}, \code{"all"}.}
}
\description{
  Simulates relative hazards, first differences, hazard
  ratios, and hazard rates for linear time-constant
  covariates from Cox Proportional Hazard models.
}
\details{
  \code{coxsimLinear} simulates relative hazards, first
  differences, and hazard ratios for time-constant
  covariates from models estimated with
  \code{\link{coxph}}.

  \code{obj} obj a coxph fitted model object.
}

