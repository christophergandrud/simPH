\name{coxsimInteract}
\alias{coxsimInteract}
\title{Simulate quantities of interest for linear multiplicative interactions.}
\usage{
  coxsimInteract(obj, b1, b2, qi = "Relative Hazard",
    X1 = NULL, X2 = NULL, nsim = 1000, ci = "95")
}
\arguments{
  \item{obj}{a coxph fitted model object with a linear
  multiplicative interaction.}

  \item{b1}{character string of the first constitutive
  variable's name.}

  \item{b2}{character string of the second constitutive
  variable's name.}

  \item{qi}{quantities of interest to simulate. Values can
  be \code{"Relative Hazard"}, \code{"Marginal Effect"},
  \code{"First Difference"}, \code{"Hazard Ratio"}, and
  \code{"Hazard Rate"}. The default is \code{qi = "Relative
  Hazard"}. If \code{qi = "Hazard Rate"} and the
  \code{coxph} model has strata, then hazard rates for each
  strata will also be calculated.}

  \item{X1}{numeric vector of fitted values of \code{b1} to
  simulate for. If \code{qi = "Marginal Effect"} then only
  \code{X2} can be set.}

  \item{X2}{numeric vector of fitted values of \code{b2} to
  simulate for. If \code{qi = "Marginal Effect"} then only
  \code{X2} can be set.}

  \item{nsim}{the number of simulations to run per value of
  X. Default is \code{nsim = 1000}.}

  \item{ci}{the proportion of middle simulations to keep.
  The default is \code{ci = "95"}, i.e. keep the middle 95
  percent. Other possibilities include: \code{"90"},
  \code{"99"}, \code{"all"}.}
}
\value{
  a siminteract class object
}
\description{
  \code{coxsimInteract} simulates quantities of interest
  for linear multiplicative interactions.
}
\references{
  Brambor, Thomas, William Roberts Clark, and Matt Golder.
  2006. “Understanding Interaction Models: Improving
  Empirical Analyses.” Political Analysis 14(1): 63–82.
}

