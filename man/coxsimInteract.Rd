\name{coxsimInteract}
\alias{coxsimInteract}
\title{Simulate quantities of interest for linear multiplicative interactions from \code{coxph} models.}
\usage{
  coxsimInteract(obj, b1, b2, qi = "Marginal Effect",
    X1 = NULL, X2 = NULL, means = FALSE, nsim = 1000,
    ci = 0.95, spin = FALSE)
}
\arguments{
  \item{obj}{a \code{coxph} fitted model object with a
  linear multiplicative interaction.}

  \item{b1}{character string of the first constitutive
  variable's name. Note \code{b1} and \code{b2} must be
  entered in the order in which they are entered into the
  \code{coxph} model.}

  \item{b2}{character string of the second constitutive
  variable's name.}

  \item{qi}{quantities of interest to simulate. Values can
  be \code{"Marginal Effect"}, \code{"First Difference"},
  \code{"Relative Hazard"}, and \code{"Hazard Rate"}. The
  default is \code{qi = "Relative Hazard"}. If \code{qi =
  "Hazard Rate"} and the \code{coxph} model has strata,
  then hazard rates for each strata will also be
  calculated.}

  \item{X1}{numeric vector of fitted values of \code{b1} to
  simulate for. If \code{qi = "Marginal Effect"} then only
  \code{X2} can be set. If you want to plot the results,
  \code{X1} should have more than one value.}

  \item{X2}{numeric vector of fitted values of \code{b2} to
  simulate for.}

  \item{means}{logical, whether or not to use the mean
  values to fit the hazard rate for covaraiates other than
  \code{b1}.}

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
  a siminteract class object
}
\description{
  \code{coxsimInteract} simulates quantities of interest
  for linear multiplicative interactions using multivariate
  normal distributions.
}
\details{
  Simulates marginal effects, first differences, hazard
  ratios, and hazard rates for linear multiplicative
  interactions. Marginal effects are calculated as in
  Brambor et al. (2006) with the addition that we take the
  exponent, so that it resembles a hazard ratio. For an
  interaction between variables \eqn{X} and \eqn{Z} then
  the marginal effect for \eqn{X} is: \deqn{ME_{X} =
  exp(\beta_{X} + \beta_{XZ}Z).}

  Note that for First Differences the comparison is not
  between two values of the same variable but two values of
  the constitute variable and 0.
}
\examples{
# Load Carpenter (2002) data
data("CarpenterFdaData")

# Load survival package
library(survival)

# Run basic model
M1 <- coxph(Surv(acttime, censor) ~ lethal*prevgenx, data = CarpenterFdaData)

# Simulate Marginal Effect of lethal for multiple values of prevgenx
Sim1 <- coxsimInteract(M1, b1 = "lethal", b2 = "prevgenx",
						  X2 = seq(2, 115, by = 2), spin = TRUE)

# Change the order of the covariates to make a more easily
# interpretable relative hazard graph.
M2 <- coxph(Surv(acttime, censor) ~ prevgenx*lethal, data = CarpenterFdaData)

# Simulate Relative Hazard of lethal for multiple values of prevgenx
Sim2 <- coxsimInteract(M2, b1 = "prevgenx", b2 = "lethal",
                    X1 = seq(2, 115, by = 2),
                    X2 = c(0, 1),
                    qi = "Relative Hazard", ci = 0.9)

# Simulate First Difference
Sim3 <- coxsimInteract(M2, b1 = "prevgenx", b2 = "lethal",
                       X1 = seq(2, 115, by = 2),
                       X2 = c(0, 1),
                       qi = "First Difference", spin = TRUE)
}
\references{
  Brambor, Thomas, William Roberts Clark, and Matt Golder.
  2006. ''Understanding Interaction Models: Improving
  Empirical Analyses.'' Political Analysis 14(1): 63–82.

  King, Gary, Michael Tomz, and Jason Wittenberg. 2000.
  ''Making the Most of Statistical Analyses: Improving
  Interpretation and Presentation.'' American Journal of
  Political Science 44(2): 347–61.

  Liu, Ying, Andrew Gelman, and Tian Zheng. 2013.
  ''Simulation-Efficient Shortest Probablility Intervals.''
  Arvix. http://arxiv.org/pdf/1302.2142v1.pdf.
}
\seealso{
  \code{\link{simGG}}, \code{\link{survival}},
  \code{\link{strata}}, and \code{\link{coxph}},
}

