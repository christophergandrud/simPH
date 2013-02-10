\name{coxsimLinear}
\alias{coxsimLinear}
\title{Simulate hazards for linear time-constant covariates from Cox Proportional Hazards models.}
\usage{
  coxsimLinear(obj, b, qi = "Relative Hazard", Xj = 1,
    Xl = 0, nsim = 1000, ci = "95")
}
\arguments{
  \item{obj}{a coxph fitted model object.}

  \item{b}{character string name of the coefficient you
  would like to simulate.}

  \item{qi}{quantity of interest to simulate. Values can be
  \code{"Relative Hazard"}, \code{"First Difference"},
  \code{"Hazard Ratio"}, and \code{"Hazard Rate"}. The
  default is \code{qi = "Relative Hazard"}. If \code{qi =
  "Hazard Rate"} and the \code{coxph} model has strata,
  then hazard rates for each strata will also be
  calculated.}

  \item{Xj}{numeric vector of values of X to simulate for.}

  \item{Xl}{numeric vector of values to compare \code{Xj}
  to. Note if \code{qi = "Relative Hazard"} only \code{Xj}
  is relevant.}

  \item{nsim}{the number of simulations to run per value of
  X. Default is \code{nsim = 1000}.}

  \item{ci}{the proportion of middle simulations to keep.
  The default is \code{ci = "95"}, i.e. keep the middle 95
  percent. Other possibilities include: \code{"90"},
  \code{"99"}, \code{"all"}.}
}
\value{
  a simlinear object
}
\description{
  Simulates relative hazards, first differences, hazard
  ratios, and hazard rates for linear time-constant
  covariates from Cox Proportional Hazard models.
}
\details{
  \code{coxsimLinear} simulates relative hazards, first
  differences, and hazard ratios for time-constant
  covariates from models estimated with \code{\link{coxph}}
  using the multivariate normal distribution.
}
\examples{
# Load Carpenter (2002) data
# data("CarpenterFdaData")

# Load survival package
# library(survival)

# Run basic model
# M1 <- coxph(Surv(acttime, censor) ~ prevgenx + lethal + deathrt1 + acutediz + hosp01  + hhosleng + mandiz01 + femdiz01 + peddiz01 + orphdum +natreg + I(natreg^2) + I(natreg^3) + vandavg3 + wpnoavg3 + condavg3 + orderent + stafcder, data = CarpenterFdaData)

# Simulate values
# Sim1 <- coxsimLinear(M1, b = "stafcder", qi = "Hazard Ratio", Xj = c(1237, 1600), Xl = c(1000, 1000), ci = "99")
}
\references{
  Licht, Amanda A. 2011. “Change Comes with Time:
  Substantive Interpretation of Nonproportional Hazards in
  Event History Analysis.” Political Analysis 19: 227–43.

  King, Gary, Michael Tomz, and Jason Wittenberg. 2000.
  “Making the Most of Statistical Analyses: Improving
  Interpretation and Presentation.” American Journal of
  Political Science 44(2): 347–61.
}
\seealso{
  \code{\link{gglinear}}, \code{\link{survival}},
  \code{\link{strata}}, and \code{\link{coxph}}
}

