\name{coxsimPoly}
\alias{coxsimPoly}
\title{Simulate quantities of interest for a range of values for a polynomial nonlinear effect from Cox Proportional Hazards models.}
\usage{
  coxsimPoly(obj, b, qi = "Hazard Ratio", pow = 2,
    Xj = NULL, Xl = NULL, nsim = 1000, ci = 0.95,
    spin = FALSE)
}
\arguments{
  \item{obj}{a coxph fitted model object with a polynomial
  coefficient.}

  \item{b}{character string name of the coefficient you
  would like to simulate.}

  \item{qi}{quantity of interest to simulate. Values can be
  \code{"Hazard Ratio"}, \code{"First Difference"}, and
  \code{"Hazard Rate"}. The default is \code{qi = "Hazard
  Ratio"}. If \code{qi = "Hazard Rate"} and the
  \code{coxph} model has strata, then hazard rates for each
  strata will also be calculated.}

  \item{pow}{numeric polynomial used in \code{coxph}.}

  \item{Xj}{numeric vector of values of X to simulate for.}

  \item{Xl}{numeric vector of values to compare \code{Xj}
  to. If \code{NULL}, then it is authomatically set to 0.}

  \item{nsim}{the number of simulations to run per value of
  \code{Xj}. Default is \code{nsim = 1000}.}

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
  a simpoly class object.
}
\description{
  Simulates quantities of interest for polynomial covariate
  effects.

  Note, you must use \code{\link{I}} to create the
  polynomials.
}
\details{
  \code{coxsimPoly} simulates quantities of interest for
  polynomial covariate effects estimated from Cox
  Proportional Hazards models.
}
\examples{
# Load Carpenter (2002) data
data("CarpenterFdaData")

# Load survival package
library(survival)

# Run basic model
M1 <- coxph(Surv(acttime, censor) ~ prevgenx + lethal + deathrt1 + acutediz +
				hosp01  + hhosleng + mandiz01 + femdiz01 + peddiz01 + orphdum +
				natreg + I(natreg^2) + I(natreg^3) + vandavg3 + wpnoavg3 +
				condavg3 + orderent + stafcder, data = CarpenterFdaData)

# Simulate simpoly First Difference
Sim1 <- coxsimPoly(M1, b = "natreg", qi = "First Difference",
						pow = 3, Xj = seq(1, 150, by = 5))

# Simulate simpoly Hazard Ratio with spin probibility interval
Sim2 <- coxsimPoly(M1, b = "natreg", qi = "Hazard Ratio"
						pow = 3, Xj = seq(1, 150, by = 5), spin = TRUE)
}
\references{
  Keele, Luke. 2010. ''Proportionally Difficult: Testing
  for Nonproportional Hazards in Cox Models.'' Political
  Analysis 18(2): 189–205.

  Carpenter, Daniel P. 2002. ''Groups, the Media, Agency
  Waiting Costs, and FDA Drug Approval.'' American Journal
  of Political Science 46(3): 490–505.
}
\seealso{
  \code{\link{simGG}}, \code{\link{survival}},
  \code{\link{strata}}, and \code{\link{coxph}}
}

