\name{coxsimPoly}
\alias{coxsimPoly}
\title{Simulate non-linear hazards for a range of values}
\usage{
  coxsimPoly(obj, b, pow = 2, X, nsim = 1000, ci = "95")
}
\arguments{
  \item{obj}{a coxph fitted model object with a polynomial
  coefficient.}

  \item{b}{character string name of the coefficient you
  would like to simulate.}

  \item{pow}{numeric polynomial used in \code{coxph}.}

  \item{X}{numeric vector of values of X to simulate
  relative hazards for.}

  \item{nsim}{the number of simulations to run per value of
  X. Default is \code{nsim = 1000}.}

  \item{ci}{the proportion of middle simulations to keep.
  The default is \code{ci = "95"}, i.e. keep the middle 95
  percent. Other possibilities include: \code{"90"},
  \code{"99"}, \code{"all"}.}
}
\value{
  a simpoly class object.
}
\description{
  Simulates relative hazards for polynomial covariate
  effects.

  Note, you must use \code{\link{I}} to create the
  polynomials.
}
\details{
  \code{coxsimPoly} simulates hazards for polynomial
  covariate effects.
}
\examples{
# Load Carpenter (2002) data
# data("CarpenterFdaData")

# Load survival package
# library(survival)

# Run basic model
# M1 <- coxph(Surv(acttime, censor) ~ prevgenx + lethal + deathrt1 + acutediz + hosp01  + hhosleng + mandiz01 + femdiz01 + peddiz01 + orphdum +natreg + I(natreg^2) + I(natreg^3) + vandavg3 + wpnoavg3 + condavg3 + orderent + stafcder, data = CarpenterFdaData)

# Simulate simpoly class object
# simM1 <- coxsimPoly(M1, b = "natreg", pow = 3, X = seq(1, 150, by = 5))
}
\references{
  Keele, Luke. 2010. “Proportionally Difficult: Testing for
  Nonproportional Hazards in Cox Models.” Political
  Analysis 18(2): 189–205.

  Carpenter, Daniel P. 2002. “Groups, the Media, Agency
  Waiting Costs, and FDA Drug Approval.” American Journal
  of Political Science 46(3): 490–505.
}
\seealso{
  \code{\link{ggpoly}}, \code{\link{rmultinorm}},
  \code{\link{survival}}, \code{\link{strata}}, and
  \code{\link{coxph}}
}

