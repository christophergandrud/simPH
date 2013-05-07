\name{simGG.simpoly}
\alias{simGG.simpoly}
\title{Plot simulated polynomial hazards.}
\usage{
  \method{simGG}{simpoly} (obj, xlab = NULL, ylab = NULL,
    title = NULL, smoother = "auto", lcolour = "#2B8CBE",
    pcolour = "#A6CEE3", lsize = 2, psize = 1,
    palpha = 0.1, ...)
}
\arguments{
  \item{obj}{a simpoly object}

  \item{xlab}{a label for the plot's x-axis.}

  \item{ylab}{a label of the plot's y-axis. The default is
  \code{ylab = "Relative Hazard \n"}.}

  \item{title}{the plot's main title}

  \item{smoother}{what type of smoothing line to use to
  summarize the plotted coefficient.}

  \item{lcolour}{character string colour of the smoothing
  line. The default is hexadecimal colour \code{lcolour =
  '#2B8CBE'}. Works if \code{strata = FALSE}.}

  \item{pcolour}{character string colour of the simulated
  points for relative hazards. Default is hexadecimal
  colour \code{pcolour = '#A6CEE3'}.}

  \item{lsize}{size of the smoothing line. Default is 2.
  See \code{\link{ggplot2}}.}

  \item{psize}{size of the plotted simulation points.
  Default is \code{psize = 1}. See \code{\link{ggplot2}}.}

  \item{palpha}{point alpha (e.g. transparency). Default is
  \code{palpha = 0.05}. See \code{\link{ggplot2}}.}

  \item{...}{other arguments passed to specific methods}
}
\value{
  a ggplot2 object.
}
\description{
  \code{simGG.simpoly} uses ggplot2 to plot simulated
  relative hazards from a \code{simpoly} class object.
}
\examples{
# Load Carpenter (2002) data
# data("CarpenterFdaData")

# Load survival package
library(survival)

# Run basic model
# M1 <- coxph(Surv(acttime, censor) ~ prevgenx + lethal + deathrt1 +
  #           acutediz + hosp01  + hhosleng + mandiz01 + femdiz01 +
  #           peddiz01 + orphdum + natreg + I(natreg^2) + vandavg3 +
  #           wpnoavg3 + condavg3 + orderent + stafcder,
  #          data = CarpenterFdaData)

# Simulate simpoly class object
# Sim1 <- coxsimPoly(M1, b = "natreg", pow = 3, X = seq(1, 150, by = 5))

# Plot simulations
# simGG(Sim1)
}
\seealso{
  \code{\link{coxsimPoly}} and \code{\link{ggplot2}}
}

