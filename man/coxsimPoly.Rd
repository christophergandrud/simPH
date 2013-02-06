\name{coxsimPoly}
\alias{coxsimPoly}
\title{Simulate non-linear hazards for a range of values}
\usage{
  coxsimPoly(obj, b, pow, from, to, by)
}
\arguments{
  \item{obj}{a coxph fitted model object with a polynomial
  coefficient.}
}
\description{
  Simulates hazards for polynomial covariate effects.
}
\details{
  \code{coxsimPoly} simulates hazards for polynomial
  covariate effects.
}
\examples{
# Load Carpenter (2002) data
data("CarpenterFdaData")

# Load survival package
library(survival)
}
\references{
  Keele, Luke. 2010. “Proportionally Difficult: Testing for
  Nonproportional Hazards in Cox Models.” Political
  Analysis 18(2): 189–205.

  Carpenter, Daniel P. 2002. “Groups, the Media, Agency
  Waiting Costs, and FDA Drug Approval.” American Journal
  of Political Science 46(3): 490–505.
}

